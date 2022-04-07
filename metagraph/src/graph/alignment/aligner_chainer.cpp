#include "aligner_chainer.hpp"

#include "aligner_seeder_methods.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"

#include "common/aligned_vector.hpp"

namespace mtg {
namespace graph {
namespace align {

using common::logger;

typedef DeBruijnGraph::node_index node_index;

constexpr uint32_t nid = std::numeric_limits<uint32_t>::max();

typedef std::tuple<Alignment::Column /* label */,
                   ssize_t /* coordinate */,
                   int32_t /* seed clipping */,
                   int32_t /* seed end */,
                   score_t /* chain score */,
                   uint32_t /* current seed index */> TableElem;
static_assert(sizeof(TableElem) == 32);
typedef AlignedVector<TableElem> ChainDPTable;

std::tuple<ChainDPTable, std::vector<uint32_t>, size_t, size_t>
chain_seeds(const IDBGAligner &aligner,
            const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds);

std::pair<size_t, size_t>
call_seed_chains_both_strands(const IDBGAligner &aligner,
                              std::string_view forward,
                              std::string_view reverse,
                              const DBGAlignerConfig &config,
                              std::vector<Seed>&& fwd_seeds,
                              std::vector<Seed>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback,
                              const std::function<bool(Alignment::Column)> &skip_column) {
    fwd_seeds.erase(std::remove_if(fwd_seeds.begin(), fwd_seeds.end(),
                                   [](const auto &a) {
                                       return a.empty() || a.label_columns.empty();
                                   }),
                    fwd_seeds.end());
    bwd_seeds.erase(std::remove_if(bwd_seeds.begin(), bwd_seeds.end(),
                                   [](const auto &a) {
                                       return a.empty() || a.label_columns.empty();
                                   }),
                    bwd_seeds.end());

    if (fwd_seeds.empty() && bwd_seeds.empty())
        return { 0, 0 };


    bool has_labels = dynamic_cast<const ILabeledAligner*>(&aligner);

    // filter out empty seeds
    std::vector<Seed> both_seeds[2];
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

    fwd_seeds = std::vector<Seed>();
    bwd_seeds = std::vector<Seed>();

    // perform chaining on the forward, and the reverse-complement seeds
    ChainDPTable dp_tables[2];
    std::vector<uint32_t> seed_backtraces[2];

    logger->trace("Chaining forward seeds");
    size_t num_seeds;
    size_t num_nodes;
    std::tie(dp_tables[0], seed_backtraces[0], num_seeds, num_nodes)
        = chain_seeds(aligner, config, forward, both_seeds[0]);

    logger->trace("Chaining reverse complement seeds");
    size_t num_seeds_bwd;
    size_t num_nodes_bwd;
    std::tie(dp_tables[1], seed_backtraces[1], num_seeds_bwd, num_nodes_bwd)
        = chain_seeds(aligner, config, reverse, both_seeds[1]);

    num_seeds += num_seeds_bwd;
    num_nodes += num_nodes_bwd;

    // construct chains by backtracking
    std::vector<std::tuple<score_t, uint32_t, ssize_t>> starts;
    starts.reserve(dp_tables[0].size() + dp_tables[1].size());
    sdsl::bit_vector both_used[2] {
        sdsl::bit_vector(dp_tables[0].size(), false),
        sdsl::bit_vector(dp_tables[1].size(), false)
    };

    for (size_t i = 0; i < dp_tables[0].size(); ++i) {
        starts.emplace_back(dp_tables[0][i].chain_score, 0, -static_cast<ssize_t>(i));
    }
    for (size_t i = 0; i < dp_tables[1].size(); ++i) {
        starts.emplace_back(dp_tables[1][i].chain_score, 1, -static_cast<ssize_t>(i));
    }

    // use heap sort to keep extracting chains
    std::make_heap(starts.begin(), starts.end());

    Chain last_chain;
    score_t last_chain_score = std::numeric_limits<score_t>::min();

    for (auto it = starts.rbegin(); it != starts.rend(); ++it) {
        std::pop_heap(starts.begin(), it.base());
        auto [chain_score, j, neg_i] = *it;
        const auto &dp_table = dp_tables[j];
        const auto &seeds = both_seeds[j];
        const auto &seed_backtrace = seed_backtraces[j];
        auto &used = both_used[j];
        uint32_t i = -neg_i;
        if (used[i])
            continue;

        // iterate through the DP table, adding seeds to the chain
        std::vector<std::pair<Seed, int64_t>> chain_seeds;

        while (i != nid) {
            const auto &[label, coord, clipping, end, score, seed_i] = dp_table[i];
            if (skip_column(label))
                break;

            used[i] = true;
            chain_seeds.emplace_back(seeds[seed_i], coord);
            chain_seeds.back().first.label_columns.clear();
            chain_seeds.back().first.label_coordinates.clear();
            if (has_labels) {
                chain_seeds.back().first.label_columns.assign(1, label);
                if (aligner.has_coordinates()) {
                    chain_seeds.back().first.label_coordinates.resize(1);
                    chain_seeds.back().first.label_coordinates[0].assign(1, coord);
                }
            }
            i = seed_backtrace[i];
        }

        if (chain_seeds.empty())
            continue;

        // clean chain by merging overlapping seeds
        if (aligner.has_coordinates()) {
            for (size_t i = chain_seeds.size() - 1; i > 0; --i) {
                auto &cur_seed = chain_seeds[i].first;
                auto &prev_seed = chain_seeds[i - 1].first;

                assert(cur_seed.size());
                assert(prev_seed.size());
                assert(prev_seed.get_clipping() < cur_seed.get_clipping());
                assert(prev_seed.get_end_clipping() > cur_seed.get_end_clipping());

                size_t prev_end = prev_seed.get_clipping()
                                    + prev_seed.get_query_view().size();
                if (prev_end > cur_seed.get_clipping()) {
                    // they overlap
                    size_t coord_dist = cur_seed.label_coordinates[0][0]
                                            + cur_seed.get_query_view().size()
                                            - prev_seed.label_coordinates[0][0]
                                            - prev_seed.get_query_view().size();
                    size_t dist = cur_seed.get_clipping()
                                    + cur_seed.get_query_view().size() - prev_end;

                    if (dist == coord_dist && cur_seed.get_nodes().size() >= dist) {
                        prev_seed.expand({ cur_seed.get_nodes().end() - dist,
                                           cur_seed.get_nodes().end() });
                        cur_seed = Seed();
                    }
                }
            }
        }

        chain_seeds.erase(std::remove_if(chain_seeds.begin(), chain_seeds.end(),
                                         [](const auto &a) { return a.first.empty(); }),
                          chain_seeds.end());
        assert(chain_seeds.size());

        for (size_t i = chain_seeds.size() - 1; i > 0; --i) {
            assert(chain_seeds[i].first.get_clipping()
                        > chain_seeds[i - 1].first.get_clipping());
            assert(chain_seeds[i].first.get_end_clipping()
                        < chain_seeds[i - 1].first.get_end_clipping());
            chain_seeds[i].second -= chain_seeds[i - 1].second;
            assert(chain_seeds[i].second > 0);
        }

        chain_seeds[0].second = 0;

        Chain chain;
        chain.reserve(chain_seeds.size());
        std::transform(chain_seeds.begin(), chain_seeds.end(), std::back_inserter(chain),
                       [&](const auto &c) {
                           return std::make_pair(Alignment(c.first, config), c.second);
                       });

        if (last_chain.empty()) {
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
            continue;
        }

        if (last_chain != chain) {
            callback(std::move(last_chain), last_chain_score);
            std::swap(last_chain, chain);
            last_chain_score = chain_score;
            continue;
        }

        if (chain[0].first.label_columns.empty())
            continue;

        // if this chain has the same seeds as the last one, merge their coordinate sets
        for (size_t i = 0; i < chain.size(); ++i) {
            Alignment::Columns columns;
            if (chain[i].first.label_coordinates.size()) {
                assert(last_chain[i].first.label_columns.size()
                        == last_chain[i].first.label_coordinates.size());
                assert(chain[i].first.label_columns.size()
                        == chain[i].first.label_coordinates.size());
                Alignment::CoordinateSet coord_union;
                auto add_col_coords = [&](auto col, const auto &coords) {
                    columns.push_back(col);
                    coord_union.push_back(coords);
                };
                utils::match_indexed_values(
                    last_chain[i].first.label_columns.begin(),
                    last_chain[i].first.label_columns.end(),
                    last_chain[i].first.label_coordinates.begin(),
                    chain[i].first.label_columns.begin(),
                    chain[i].first.label_columns.end(),
                    chain[i].first.label_coordinates.begin(),
                    [&](auto col, const auto &coords, const auto &other_coords) {
                        columns.push_back(col);
                        coord_union.emplace_back();
                        std::set_union(coords.begin(), coords.end(),
                                       other_coords.begin(), other_coords.end(),
                                       std::back_inserter(coord_union.back()));
                    },
                    add_col_coords, add_col_coords
                );
                std::swap(last_chain[i].first.label_coordinates, coord_union);
            } else {
                assert(chain[i].first.label_columns.size());
                std::set_union(last_chain[i].first.label_columns.begin(),
                               last_chain[i].first.label_columns.end(),
                               chain[i].first.label_columns.begin(),
                               chain[i].first.label_columns.end(),
                               std::back_inserter(columns));
            }
            std::swap(last_chain[i].first.label_columns, columns);
        }
    }

    if (last_chain.size())
        callback(std::move(last_chain), last_chain_score);

    return std::make_pair(num_seeds, num_nodes);
}

std::tuple<ChainDPTable, std::vector<uint32_t>, size_t, size_t>
chain_seeds(const IDBGAligner &aligner,
            const DBGAlignerConfig &config,
            std::string_view query,
            std::vector<Seed> &seeds) {
    const auto *labeled_aligner = dynamic_cast<const ILabeledAligner*>(&aligner);
    if (seeds.empty() || !labeled_aligner || !aligner.has_coordinates())
        return {};

    AnnotationBuffer &buffer = labeled_aligner->get_annotation_buffer();

    size_t num_nodes = 0;

    ssize_t query_size = query.size();

    ChainDPTable dp_table;
    dp_table.reserve(seeds.size());

    logger->trace("Sorting seeds");
    std::sort(seeds.begin(), seeds.end(), [](const auto &a, const auto &b) {
        return a.get_clipping() > b.get_clipping();
    });

    tsl::hopscotch_map<Alignment::Column, size_t> label_sizes;

    for (size_t i = 0; i < seeds.size(); ++i) {
        if (!(i % config.row_batch_size)) {
            size_t end = std::min(i + config.row_batch_size, seeds.size());
            std::vector<node_index> nodes;
            nodes.reserve(end - i);
            for (size_t j = i; j < end; ++j) {
                nodes.emplace_back(seeds[j].get_nodes()[0]);
            }
            buffer.prefetch_coords(nodes);
        }

        auto fetch_coords = buffer.get_coords(seeds[i].get_nodes()[0], false,
                                              &seeds[i].label_columns);
        assert(fetch_coords);

        for (size_t j = 0; j < fetch_coords->size(); ++j) {
            Alignment::Column c = seeds[i].label_columns[j];
            size_t num_seeds = std::min((*fetch_coords)[j].size(),
                                        config.max_num_seeds_per_locus);
            auto rbegin = (*fetch_coords)[j].rbegin();
            auto rend = rbegin + num_seeds;
            std::for_each(rbegin, rend, [&](ssize_t coord) {
                ++label_sizes[c];
                dp_table.emplace_back(c, coord + seeds[i].get_offset(), seeds[i].get_clipping(),
                                      seeds[i].get_clipping() + seeds[i].get_query_view().size(),
                                      seeds[i].get_query_view().size(), i);
            });
        }
    }

    size_t num_seeds = dp_table.size();
    std::vector<uint32_t> backtrace(dp_table.size(), nid);
    if (dp_table.empty())
        return std::make_tuple(std::move(dp_table), std::move(backtrace), num_seeds, num_nodes);

    logger->trace("Sorting {} anchors", dp_table.size());
    // sort seeds by label, then by decreasing reference coordinate
    std::sort(dp_table.begin(), dp_table.end(), std::greater<TableElem>());
    logger->trace("Chaining anchors");

    // scoring function derived from minimap2
    // https://academic.oup.com/bioinformatics/article/34/18/3094/4994778
#define LOG2(X) ((unsigned) (8 * sizeof(unsigned long long) - __builtin_clzll((X)) - 1))
    auto gap_penalty = [sl=config.min_seed_length](size_t len) -> score_t {
        return !len ? 0 : (len * sl + 127) / 128 + (LOG2(len) / 2);
    };

    size_t bandwidth = 201;

    size_t cur_label_end = 0;
    size_t i = 0;
    while (cur_label_end < dp_table.size()) {
        cur_label_end += label_sizes[std::get<0>(dp_table[i])];
        for ( ; i < cur_label_end; ++i) {
            const auto &[prev_label, prev_coord, prev_clipping, prev_end,
                         prev_score, prev_seed_i] = dp_table[i];

            if (!prev_clipping)
                continue;

            size_t end = std::min(bandwidth, cur_label_end - i) + i;
            ssize_t coord_cutoff = prev_coord - query_size;

            for (size_t j = i + 1; j < end; ++j) {
                auto &[label, coord, clipping, end, score, seed_i] = dp_table[j];
                assert(label == prev_label);

                if (coord_cutoff > coord)
                    break;

                if (clipping >= prev_clipping || end >= prev_end)
                    continue;

                int32_t dist = prev_clipping - clipping;
                if (dist > query_size)
                    continue;

                int32_t coord_dist = prev_coord - coord;
                score_t cur_score = prev_score + std::min({ end - clipping, dist, coord_dist })
                    - gap_penalty(std::abs(coord_dist - dist));

                if (cur_score >= score) {
                    score = cur_score;
                    backtrace[j] = i;
                }
            }
        }
    }

    return std::make_tuple(std::move(dp_table), std::move(backtrace), num_seeds, num_nodes);
}

template <class AlignmentCompare>
void construct_alignment_chain(size_t node_overlap,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback);

template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query,
                                        const DBGAlignerConfig &config,
                                        size_t node_overlap) {
    if (alignments.size() < 2 || !config.post_chain_alignments)
        return std::move(alignments);

    for (const auto &a : alignments) {
        if (a.label_coordinates.size())
            throw std::runtime_error("Post-chaining alignments with coordinates not supported");
    }

    DBGAlignerConfig no_chain_config { config };
    no_chain_config.post_chain_alignments = false;
    AlignmentAggregator<AlignmentCompare> aggregator(no_chain_config);

    alignments.erase(std::remove_if(alignments.begin(), alignments.end(), [&](Alignment &a) {
        if (!a.get_clipping() && !a.get_end_clipping()) {
            aggregator.add_alignment(std::move(a));
            return true;
        }

        return false;
    }), alignments.end());

    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        return std::make_tuple(a.get_orientation(),
                               a.get_clipping() + a.get_query_view().size(),
                               a.get_clipping(),
                               b.get_score(),
                               a.get_sequence().size())
            < std::make_tuple(b.get_orientation(),
                              b.get_clipping() + b.get_query_view().size(),
                              b.get_clipping(),
                              a.get_score(),
                              b.get_sequence().size());
    });

    DEBUG_LOG("Chaining alignments:\n{}", fmt::join(alignments, "\t\n"));


    auto run = [&](std::string_view this_query, auto begin, auto end) {
        std::vector<score_t> best_score(this_query.size() + 1, 0);
        for (auto it = begin; it != end; ++it) {
            size_t end_pos = it->get_query_view().data() + it->get_query_view().size()
                                - this_query.data();
            if (it->get_score() > best_score[end_pos]) {
                best_score[end_pos] = it->get_score();
                construct_alignment_chain<AlignmentCompare>(
                    node_overlap, config, this_query, Alignment(*it), it + 1, end, &best_score,
                    [&](Alignment&& chain) { aggregator.add_alignment(std::move(chain)); }
                );
            }
        }
    };

    // recursively construct chains
    auto split_it = std::find_if(alignments.begin(), alignments.end(),
                                 [](const auto &a) { return a.get_orientation(); });
    run(query, alignments.begin(), split_it);
    run(rc_query, split_it, alignments.end());

    return aggregator.get_alignments();
}

// TODO: rewrite this to not use recursion
template <class AlignmentCompare>
void construct_alignment_chain(size_t node_overlap,
                               const DBGAlignerConfig &config,
                               std::string_view query,
                               Alignment&& chain,
                               typename std::vector<Alignment>::iterator begin,
                               typename std::vector<Alignment>::iterator end,
                               std::vector<score_t> *best_score,
                               const std::function<void(Alignment&&)> &callback) {
    assert(begin <= end);
    assert(chain.size());

    const char *chain_begin = chain.get_query_view().data();
    const char *chain_end = chain.get_query_view().data() + chain.get_query_view().size();
    if (begin == end || chain_end == query.data() + query.size()) {
        callback(std::move(chain));
        return;
    }

    score_t score = chain.get_score();

    bool called = false;
    for (auto it = begin; it != end; ++it) {
        // TODO: handle this case later
        if (it->get_offset())
            continue;

        const char *next_begin = it->get_query_view().data();
        const char *next_end = it->get_query_view().data() + it->get_query_view().size();

        assert(chain_begin - chain.get_clipping() == next_begin - it->get_clipping());
        assert(it->get_orientation() == chain.get_orientation());

        if (next_begin <= chain_begin || next_end == chain_end)
            continue;

        if (chain.label_columns.size()
                && !utils::share_element(it->label_columns.begin(),
                                         it->label_columns.end(),
                                         chain.label_columns.begin(),
                                         chain.label_columns.end())) {
            continue;
        }

        Alignment aln = *it;

        if (next_begin >= chain_end) {
            // no overlap
            aln.insert_gap_prefix(next_begin - chain_end, node_overlap, config);

        } else {
            // trim, then fill in dummy nodes
            assert(chain.get_end_clipping());

            // first trim front of the incoming alignment
            size_t overlap = std::min(
                static_cast<size_t>((chain.get_cigar().data().end() - 2)->second),
                aln.trim_query_prefix(chain_end - it->get_query_view().data(),
                                      node_overlap, config)
            );

            if (aln.empty() || aln.get_sequence().size() <= node_overlap
                    || (aln.get_cigar().data().begin()
                            + static_cast<bool>(aln.get_clipping()))->first != Cigar::MATCH) {
                continue;
            }

            assert(aln.get_query_view().data()
                    == chain.get_query_view().data() + chain.get_query_view().size());

            if (overlap < node_overlap) {
                aln.insert_gap_prefix(-overlap, node_overlap, config);
            } else {
                aln.trim_clipping();
            }
        }

        assert(!aln.empty());

        score_t next_score = score + aln.get_score();
        if (next_score <= (*best_score)[next_end - query.data()])
            continue;

        (*best_score)[next_end - query.data()] = next_score;
        // use append instead of splice because any clipping in aln represents
        // internally clipped characters
        Alignment next_chain = chain;
        next_chain.trim_end_clipping();
        bool changed = next_chain.append(std::move(aln));
        if (next_chain.size()) {
            assert(next_chain.get_score() == next_score);
            construct_alignment_chain<AlignmentCompare>(
                    node_overlap, config, query, std::move(next_chain),
                    it + 1, end, best_score, callback);
            called |= changed;
        }
    }

    if (!called)
        callback(std::move(chain));
}

template
std::vector<Alignment> chain_alignments<LocalAlignmentLess>(std::vector<Alignment>&&,
                                                            std::string_view,
                                                            std::string_view,
                                                            const DBGAlignerConfig&,
                                                            size_t);

} // namespace align
} // namespace graph
} // namespace mtg
