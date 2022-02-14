#include "aligner_chainer.hpp"

#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"

namespace mtg {
namespace graph {
namespace align {

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

typedef DeBruijnGraph::node_index node_index;

constexpr size_t nid = std::numeric_limits<uint32_t>::max();

typedef std::tuple<Alignment::Column /* label */,
                   ssize_t /* coordinate */,
                   size_t /* seed clipping */,
                   Alignment::node_index /* first node of seed */,
                   ssize_t /* seed length */,
                   score_t /* chain score */,
                   uint32_t /* previous seed index */,
                   uint32_t /* current seed index */> TableElem;
typedef std::vector<TableElem> ChainDPTable;

std::tuple<ChainDPTable, size_t, size_t>
chain_seeds(const DBGAlignerConfig &config,
            std::string_view query,
            const std::vector<Alignment> &seeds);

std::pair<size_t, size_t>
call_seed_chains_both_strands(std::string_view forward,
                              std::string_view reverse,
                              const DeBruijnGraph &graph,
                              const DBGAlignerConfig &config,
                              std::vector<Alignment>&& fwd_seeds,
                              std::vector<Alignment>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback,
                              const std::function<bool()> &terminate) {
    if (terminate())
        return {};

    // filter out empty seeds
    std::vector<Alignment> both_seeds[2];
    both_seeds[0].reserve(fwd_seeds.size());
    both_seeds[1].reserve(bwd_seeds.size());
    for (auto&& fwd_seed : fwd_seeds) {
        if (!fwd_seed.empty())
            both_seeds[0].emplace_back(std::move(fwd_seed));
    }
    for (auto&& bwd_seed : bwd_seeds) {
        if (!bwd_seed.empty())
            both_seeds[1].emplace_back(std::move(bwd_seed));
    }

    fwd_seeds = std::vector<Alignment>();
    bwd_seeds = std::vector<Alignment>();

    // perform chaining on the forward, and the reverse-complement seeds
    ChainDPTable dp_tables[2];
    size_t num_seeds;
    size_t num_nodes;
    std::tie(dp_tables[0], num_seeds, num_nodes)
        = chain_seeds(config, forward, both_seeds[0]);

    size_t num_seeds_bwd;
    size_t num_nodes_bwd;
    std::tie(dp_tables[1], num_seeds_bwd, num_nodes_bwd)
        = chain_seeds(config, reverse, both_seeds[1]);

    num_seeds += num_seeds_bwd;
    num_nodes += num_nodes_bwd;

    // construct chains by backtracking
    std::vector<std::tuple<score_t, size_t, ssize_t>> starts;
    sdsl::bit_vector both_used[2] {
        sdsl::bit_vector(dp_tables[0].size(), false),
        sdsl::bit_vector(dp_tables[1].size(), false)
    };

    for (size_t j : { 0, 1 }) {
        for (size_t i = 0; i < dp_tables[j].size(); ++i) {
            starts.emplace_back(
                std::get<5>(dp_tables[j][i]),
                j,
                -static_cast<ssize_t>(i)
            );
        }
    }

    // use heap sort to keep extracting chains until terminate() returns true
    std::make_heap(starts.begin(), starts.end());

    Chain last_chain;
    score_t last_chain_score = std::numeric_limits<score_t>::min();

    for (auto it = starts.rbegin(); it != starts.rend() && !terminate(); ++it) {
        std::pop_heap(starts.begin(), it.base());
        score_t chain_score = std::get<0>(*it);
        const auto &dp_table = dp_tables[std::get<1>(*it)];
        const auto &seeds = both_seeds[std::get<1>(*it)];
        auto &used = both_used[std::get<1>(*it)];
        size_t i = -std::get<2>(*it);
        if (used[i])
            continue;

        // iterate through the DP table, adding seeds to the chain
        Chain chain;

        while (i != nid) {
            const auto &[label, coord, clipping, node, length, score, pred, seed_i] = dp_table[i];
            used[i] = true;
            chain.push_back(seeds[seed_i]);
            chain.back().label_columns.clear();
            chain.back().label_coordinates.clear();
            chain.back().label_columns.assign(1, label);
            chain.back().label_coordinates.resize(1);
            chain.back().label_coordinates[0].assign(1, coord);
            i = pred;
        }

        if (chain.empty())
            continue;

        // clean chain by merging overlapping seeds
        for (size_t i = chain.size() - 1; i > 0; --i) {
            const char *prev_end = chain[i - 1].get_query().data() + chain[i - 1].get_query().size();
            const char *cur_begin = chain[i].get_query().data();
            if (prev_end + chain[i].get_offset() >= cur_begin + graph.get_k() - 1) {
                Alignment trim = chain[i];
                Alignment merge = chain[i - 1];
                merge.trim_end_clipping();
                trim.trim_query_prefix(prev_end - cur_begin, graph, config);
                trim.trim_clipping();
                assert(trim.get_offset() == graph.get_k() - 1);
                merge.append(std::move(trim));
                if (merge.size()) {
                    assert(merge.is_valid(graph, &config));
                    std::swap(merge, chain[i - 1]);
                    chain[i] = Alignment();
                }
            }
        }

        auto chain_it = std::remove_if(chain.begin(), chain.end(),
                                       [](const auto &a) { return a.empty(); });

        chain.erase(chain_it, chain.end());

        if (last_chain.empty()) {
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
            continue;
        }

        // if this chain has the same seeds as the last one, merge their coordinate sets
        if (last_chain == chain) {
            for (size_t i = 0; i < chain.size(); ++i) {
                Alignment::LabelSet label_union;
                Alignment::CoordinateSet coord_union;
                auto a_begin = last_chain[i].label_columns.begin();
                auto a_end = last_chain[i].label_columns.end();
                auto a_c_begin = last_chain[i].label_coordinates.begin();
                auto b_begin = chain[i].label_columns.begin();
                auto b_end = chain[i].label_columns.end();
                auto b_c_begin = chain[i].label_coordinates.begin();
                while (a_begin != a_end && b_begin != b_end) {
                    if (*a_begin < *b_begin) {
                        label_union.push_back(*a_begin);
                        coord_union.push_back(*a_c_begin);
                        ++a_begin;
                        ++a_c_begin;
                    } else if (*a_begin > *b_begin) {
                        label_union.push_back(*b_begin);
                        coord_union.push_back(*b_c_begin);
                        ++b_begin;
                        ++b_c_begin;
                    } else {
                        label_union.push_back(*a_begin);
                        coord_union.emplace_back();
                        std::set_union(a_c_begin->begin(), a_c_begin->end(),
                                       b_c_begin->begin(), b_c_begin->end(),
                                       std::back_inserter(coord_union.back()));
                        ++a_begin;
                        ++a_c_begin;
                        ++b_begin;
                        ++b_c_begin;
                    }
                }
                label_union.insert(label_union.end(), a_begin, a_end);
                label_union.insert(label_union.end(), b_begin, b_end);
                coord_union.insert(coord_union.end(), a_c_begin, last_chain[i].label_coordinates.end());
                coord_union.insert(coord_union.end(), b_c_begin, chain[i].label_coordinates.end());
                std::swap(last_chain[i].label_columns, label_union);
                std::swap(last_chain[i].label_coordinates, coord_union);
            }
        } else {
            callback(std::move(last_chain), last_chain_score);
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
        }
    }

    if (!terminate() && last_chain.size())
        callback(std::move(last_chain), last_chain_score);

    return std::make_pair(num_seeds, num_nodes);
}

std::tuple<ChainDPTable, size_t, size_t>
chain_seeds(const DBGAlignerConfig &config,
            std::string_view query,
            const std::vector<Alignment> &seeds) {
    if (seeds.empty())
        return {};

    size_t num_nodes = 0;
    size_t num_seeds = 0;

    ssize_t query_size = query.size();

    ChainDPTable dp_table;
    dp_table.reserve(num_seeds);

    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty() || seeds[i].label_columns.empty()
                || seeds[i].label_coordinates.empty())
            continue;

        size_t seed_count = 0;
        for (const auto &tuple : seeds[i].label_coordinates) {
            seed_count += tuple.size();
        }

        if (seed_count > config.max_num_seeds_per_locus)
            continue;

        num_seeds += seed_count;

        for (size_t j = 0; j < seeds[i].label_coordinates.size(); ++j) {
            Alignment::Column c = seeds[i].label_columns[j];
            for (auto coord : seeds[i].label_coordinates[j]) {
                dp_table.emplace_back(c,
                                      coord,
                                      seeds[i].get_clipping(),
                                      seeds[i].get_nodes()[0],
                                      seeds[i].get_query().size(),
                                      seeds[i].get_query().size(),
                                      nid,
                                      i);
            }
        }
    }

    // sort seeds by label, then by decreasing reference coordinate
    std::sort(dp_table.begin(), dp_table.end(), std::greater<TableElem>());

    // scoring function derived from minimap2
    // https://academic.oup.com/bioinformatics/article/34/18/3094/4994778
    auto gap_penalty = [sl=config.min_seed_length](size_t len) -> score_t {
        return !len ? 0 : (len * sl + 99) / 100 + (LOG2(len) / 2);
    };

    size_t bandwidth = 200;

    for (size_t i = 0; i < dp_table.size() - 1; ++i) {
        const auto &[prev_label, prev_coord, prev_clipping, prev_node, prev_length,
                     prev_score, prev_pred, prev_seed_i] = dp_table[i];
        if (!prev_clipping)
            continue;

        size_t end = std::min(bandwidth, dp_table.size() - i) + i;

        for (size_t j = i + 1; j < end; ++j) {
            auto &[label, coord, clipping, node, length, score, pred, seed_i] = dp_table[j];
            if (label != prev_label)
                break;

            if (clipping >= prev_clipping)
                continue;

            ssize_t dist = prev_clipping - clipping;
            if (dist > query_size || prev_clipping + prev_length == clipping + length)
                continue;

            ssize_t coord_dist = prev_coord - coord;

            if (coord_dist > query_size)
                break;

            score_t cur_score = prev_score + std::min(std::min(length, dist), coord_dist)
                - gap_penalty(std::abs(coord_dist - dist));

            if (cur_score >= score) {
                score = cur_score;
                pred = i;
            }
        }
    }

    return std::make_tuple(std::move(dp_table), num_seeds, num_nodes);
}

template <class AlignmentCompare>
void construct_alignment_chain(const DeBruijnGraph &graph,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> &best_score,
                               const std::function<void(Alignment&&)> &callback);

template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query,
                                        const DBGAlignerConfig &config,
                                        const DeBruijnGraph &graph) {
    if (alignments.size() < 2 || !config.post_chain_alignments)
        return std::move(alignments);

    if (std::any_of(alignments.begin(), alignments.end(), [](const auto &a) {
        return a.label_coordinates.size();
    })) {
        throw std::runtime_error("Post-chaining alignments with coordinates not supported");
    }

    DBGAlignerConfig no_chain_config { config };
    no_chain_config.post_chain_alignments = false;
    AlignmentAggregator<AlignmentCompare> aggregator(no_chain_config);

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_orientation(),
                               a.get_clipping() + a.get_query().size(),
                               a.get_clipping(),
                               b.get_score(),
                               a.get_sequence().size())
            < std::make_tuple(b.get_orientation(),
                              b.get_clipping() + b.get_query().size(),
                              b.get_clipping(),
                              a.get_score(),
                              b.get_sequence().size());
    });

    auto run = [&](bool rev_compl, auto begin, auto end) {
        std::string_view this_query = rev_compl ? rc_query : query;
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = begin; it != end; ++it) {
            size_t end_pos = it->get_query().data() + it->get_query().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain<AlignmentCompare>(
                    graph, config, this_query, Alignment(*it), it + 1, end, best_score,
                    [&](Alignment&& chain) {
                      aggregator.add_alignment(std::move(chain));
                    }
                );
            }
        }
    };

    // recursively construct chains
    auto split_it = std::find_if(alignments.begin(), alignments.end(),
                                 [](const auto &a) { return a.get_orientation(); });
    run(false, alignments.begin(), split_it);
    run(true, split_it, alignments.end());

    return aggregator.get_alignments();
}

template <class AlignmentCompare>
void construct_alignment_chain(const DeBruijnGraph &graph,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> &best_score,
                               const std::function<void(Alignment&&)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query().data();
    const char *chain_end = chain_begin + chain.get_query().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    size_t k = graph.get_k();
    score_t score = chain.get_score();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        if (it->get_offset())
            continue;

        const char *next_begin = it->get_query().data();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        const char *next_end = next_begin + it->get_query().size();

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        if (chain.label_columns.size()
                && !utils::share_element(it->label_columns.begin(),
                                         it->label_columns.end(),
                                         chain.label_columns.begin(),
                                         chain.label_columns.end())) {
            continue;
        }

        Alignment aln(*it);

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, graph, config);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().data().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query().data(), graph, config)
            );

            if (aln.empty() || aln.get_sequence().size() < k
                    || (aln.get_cigar().data().begin()
                            + static_cast<bool>(aln.get_clipping()))->first != Cigar::MATCH) {
                continue;
            }

            assert(aln.get_query().data() == chain.get_query().data() + chain.get_query().size());

            if (overlap < k - 1)
                aln.insert_gap_prefix(-overlap, graph, config);
        }

        if (aln.empty())
            continue;

        score_t next_score = score + aln.get_score();
        if (next_score > best_score[next_end - query.data()]) {
            best_score[next_end - query.data()] = next_score;

            Alignment next_chain(chain);
            next_chain.trim_end_clipping();
            bool modified = next_chain.append(std::move(aln));
            assert(next_chain.get_score() == next_score);
            assert(next_chain.is_valid(graph, &config));
            if (next_chain.size()) {
                called |= modified;
                construct_alignment_chain<AlignmentCompare>(
                    graph, config, query, std::move(next_chain),
                    it + 1, end, best_score, callback
                );
            }
        }
    }

    if (!called)
        callback(std::move(chain));
}

template std::vector<Alignment> chain_alignments<LocalAlignmentLess>(std::vector<Alignment>&&,
                                                                     std::string_view,
                                                                     std::string_view,
                                                                     const DBGAlignerConfig&,
                                                                     const DeBruijnGraph&);

} // namespace align
} // namespace graph
} // namespace mtg
