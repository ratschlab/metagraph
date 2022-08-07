#include "dbg_aligner.hpp"

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "aligner_labeled.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/graph_extensions/unitigs.hpp"

namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

const DBGSuccinct* get_dbg_succ(const DeBruijnGraph &graph) {
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    return dynamic_cast<const DBGSuccinct*>(&(canonical ? canonical->get_graph() : graph));
}

AlignmentResults IDBGAligner::align(std::string_view query) const {
    AlignmentResults result;
    align_batch({ Query{ std::string{}, query } },
        [&](const std::string&, AlignmentResults&& alignment) {
            std::swap(result, alignment);
        }
    );

    return result;
}

template <class Seeder, class Extender, class AlignmentCompare>
DBGAligner<Seeder, Extender, AlignmentCompare>
::DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
      : graph_(graph), config_(config) {
    if (!config_.min_seed_length)
        config_.min_seed_length = graph_.get_k();

    if (config_.min_seed_length < graph.get_k() && !get_dbg_succ(graph))
        config_.min_seed_length = graph.get_k();

    if (!config_.max_seed_length)
        config_.max_seed_length = graph_.get_k();

    std::tie(config_.min_seed_length, config_.max_seed_length)
        = std::make_pair(std::min(config_.min_seed_length, config_.max_seed_length),
                         std::max(config_.min_seed_length, config_.max_seed_length));

    assert(config_.max_seed_length >= config_.min_seed_length);
    assert(config_.num_alternative_paths);
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
        && "primary graphs must be wrapped into canonical");

    if (!config_.check_config_scores())
        throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");

    // extensions should not trim characters from chain seed prefixes
    if (config_.chain_alignments)
        config_.allow_left_trim = false;
}

/**
 * Partition the alignment at the last k-mer. Return a pair containing the
 * alignment of all but the last k-mers, and the alignment of the last k-mer.
 */
std::pair<Alignment, Alignment> split_seed(const DeBruijnGraph &graph,
                                           const DBGAlignerConfig &config,
                                           const Alignment &alignment) {
    assert(alignment.is_valid(graph, &config));
    if (alignment.get_sequence().size() < graph.get_k() * 2)
        return std::make_pair(Alignment(), alignment);

    auto ret_val = std::make_pair(alignment, alignment);
    assert(std::find(alignment.get_nodes().begin(), alignment.get_nodes().end(), DeBruijnGraph::npos)
            == alignment.get_nodes().end());
    ret_val.first.trim_reference_suffix(graph.get_k(), config, false);

    // ensure that there's no DELETION at the splice point
    size_t trim_nodes = graph.get_k();
    if (ret_val.first.size()) {
        auto it = ret_val.first.get_cigar().data().rbegin();
        if (it->first == Cigar::CLIPPED)
            ++it;

        if (it->first == Cigar::DELETION) {
            trim_nodes += it->second;
            ret_val.first.trim_reference_suffix(it->second, config, false);
        }

        assert(ret_val.first.size());
    }

    if (ret_val.first.empty())
        return std::make_pair(Alignment(), alignment);

    ret_val.second.trim_reference_prefix(alignment.get_sequence().size() - trim_nodes,
                                         graph.get_k() - 1, config, true);

    assert(ret_val.first.is_valid(graph, &config));
    assert(ret_val.second.is_valid(graph, &config));

    return ret_val;
}

void filter_seed(const Alignment &prev, Alignment &a) {
    if (!prev.has_annotation()) {
        DEBUG_LOG("Skipping seed {}", a);
        a = Alignment();
    } else if (prev.label_coordinates.empty()) {
        Vector<Alignment::Column> diff;
        std::set_difference(a.get_columns().begin(),
                            a.get_columns().end(),
                            prev.get_columns().begin(),
                            prev.get_columns().end(),
                            std::back_inserter(diff));
        if (diff.empty()) {
            DEBUG_LOG("Skipping seed {}", a);
            a = Alignment();
        } else {
            a.set_columns(std::move(diff));
        }
    } else {
        Vector<Alignment::Column> diff;
        Vector<Alignment::Tuple> diff_coords;
        utils::match_indexed_values(
            a.get_columns().begin(), a.get_columns().end(), a.label_coordinates.begin(),
            prev.get_columns().begin(), prev.get_columns().end(), prev.label_coordinates.begin(),
            [&](auto col, const auto &coords, const auto &other_coords) {
                Alignment::Tuple set;
                // filter_seed: clear the seed a if it has no unexplored labels or coordinates
                // relative to the seed prev
                std::set_difference(coords.begin(), coords.end(),
                                    other_coords.begin(), other_coords.end(),
                                    std::back_inserter(set));
                if (set.size()) {
                    diff.push_back(col);
                    diff_coords.push_back(std::move(set));
                }
            }
        );
        if (diff.empty()) {
            DEBUG_LOG("Skipping seed {}", a);
            a = Alignment();
        } else {
            a.set_columns(std::move(diff));
            std::swap(a.label_coordinates, diff_coords);
        }
    }

}

// Extend the alignment first until it reaches the end of the alignment second.
// Return whether the connection was successful. If not, then replace first with second
// and push first to the back of partial_alignments.
template <class Extender>
bool align_connect(const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config,
                   Alignment &first,
                   Alignment &second,
                   int64_t coord_dist,
                   Extender &extender,
                   std::vector<Alignment> &partial_alignments) {
    auto [left, next] = split_seed(graph, config, first);
    coord_dist += second.get_sequence().size() + next.get_sequence().size()
            - first.get_sequence().size();
    assert(coord_dist > 0);

    auto extensions = extender.get_extensions(next,
        config.ninf,                         // min_path_score
        true,                                // force_fixed_seed
        coord_dist,                          // target_length
        second.get_nodes().back(),           // target_node
        false,                               // trim_offset_after_extend
        second.get_end_clipping(),           // trim_query_suffix
        first.get_score() - next.get_score() // added_xdrop
    );

    if (extensions.size() && extensions[0].get_end_clipping() < first.get_end_clipping()) {
        assert(extensions[0].get_nodes().front() == next.get_nodes().front());
        left.splice(std::move(extensions[0]));
        std::swap(left, first);
        assert(first.is_valid(graph, &config));

        if (first.get_score() >= second.get_score())
            return true;
    }

    logger->trace("Extension score too low, restarting from seed {}", second);
    partial_alignments.emplace_back(first);
    std::swap(first, second);
    return false;
}

typedef boss::BOSS::edge_index edge_index;
typedef std::pair<edge_index, edge_index> BOSSRange;

std::pair<std::unique_ptr<boss::BOSS>, std::vector<BOSSRange>>
extract_subgraph(const DeBruijnGraph &graph,
                 const std::vector<IDBGAligner::Query> &seq_batch,
                 const DBGAlignerConfig &config,
                 size_t batch_size) {
    if (config.min_seed_length >= graph.get_k())
        return {};

    const auto *base_dbg_succ = get_dbg_succ(graph);
    if (!base_dbg_succ)
        return {};

    size_t query_total_length = 0;
    bool add_rev_comp = graph.get_mode() == DeBruijnGraph::CANONICAL
                                    || config.forward_and_reverse_complement;
    boss::BOSSConstructor constructor(config.min_seed_length - 1, add_rev_comp, 0,
                                      "", 1, batch_size, kmer::ContainerType::VECTOR, "");
    for (const auto &[header, query] : seq_batch) {
        constructor.add_sequence(query);
        query_total_length += query.size() + (query.size() * add_rev_comp);
    }
    auto boss = std::make_unique<boss::BOSS>(&constructor);

    logger->trace("Getting matching ranges");
    typedef boss::BOSS::TAlphabet TAlphabet;
    const auto &base_boss = base_dbg_succ->get_boss();
    std::vector<std::tuple<edge_index, edge_index, SmallVector<TAlphabet>>> initial_edge_stack;
    for (TAlphabet s = 1; s < boss->alph_size; ++s) {
        edge_index rl = boss->get_F(s) + 1 < boss->get_W().size()
             ? boss->get_F(s) + 1
             : boss->get_W().size(); // lower bound
        edge_index ru = s + 1 < boss->alph_size
             ? boss->get_F(s + 1)
             : boss->get_W().size() - 1; // upper bound

        if (rl <= ru)
            initial_edge_stack.emplace_back(rl, ru, SmallVector<TAlphabet>{ s });
    }

    std::vector<std::tuple<edge_index, edge_index, size_t>> edge_stack;
    std::vector<std::pair<BOSSRange, size_t>> traversal_stack;
    while (initial_edge_stack.size()) {
        auto [first, last, spelling] = initial_edge_stack.back();
        initial_edge_stack.pop_back();
        if (spelling.size() >= base_boss.get_indexed_suffix_length()) {
            auto [base_first, base_last, match_end]
                = base_boss.index_range(spelling.begin(), spelling.end());
            if (match_end == spelling.end()) {
                edge_stack.emplace_back(first, last, spelling.size());
                traversal_stack.emplace_back(std::make_pair(base_first, base_last), spelling.size());
            }
            continue;
        }

        spelling.push_back(0);
        for (TAlphabet s = 1; s < boss->alph_size; ++s) {
            edge_index next_first = first;
            if (spelling.size() <= boss->get_k()) {
                edge_index next_last = last;
                if (boss->tighten_range(&next_first, &next_last, s)) {
                    spelling.back() = s;
                    initial_edge_stack.emplace_back(next_first, next_last, spelling);
                }
            } else if (edge_index next_last = boss->pick_edge(last, s)) {
                spelling.back() = s;
                initial_edge_stack.emplace_back(next_last, next_last, spelling);
            }
        }
    }

    std::vector<BOSSRange> matching_ranges(boss->num_edges() + 1);
    size_t total_matches = 0;
    size_t total_ranges = 0;
    while (traversal_stack.size()) {
        auto [base_range, l] = traversal_stack.back();
        traversal_stack.pop_back();

        assert(edge_stack.size());
        auto [first, last, length] = edge_stack.back();
        assert(l == length);
        edge_stack.pop_back();

        if (l == config.min_seed_length) {
            matching_ranges[last] = base_range;
            total_matches += last - first + 1;
            ++total_ranges;
            continue;
        }

        if (++length == config.min_seed_length) {
            if (base_range.first != base_range.second) {
                boss->call_outgoing(last, [&](edge_index next_last) {
                    if (TAlphabet s = boss->get_W(next_last)) {
                        auto [next_base_first, next_base_last] = base_range;
                        if (base_boss.tighten_range(&next_base_first, &next_base_last, s % boss->alph_size)) {
                            edge_stack.emplace_back(next_last, next_last, length);
                            traversal_stack.emplace_back(
                                std::make_pair(next_base_first, next_base_last), length);
                        }
                    }
                });
            } else if (TAlphabet s = base_boss.get_W(base_range.second)) {
                s %= base_boss.alph_size;
                if (auto next_last = boss->pick_edge(last, s)) {
                    auto next_base_last = base_boss.fwd(base_range.second, s);
                    auto next_base_first = base_boss.pred_last(next_base_last - 1) + 1;
                    edge_stack.emplace_back(next_last, next_last, length);
                    traversal_stack.emplace_back(
                        std::make_pair(next_base_first, next_base_last), length);
                }
            }
            continue;
        }

        if (base_range.first == base_range.second) {
            if (TAlphabet s = base_boss.get_W(base_range.second)) {
                s %= base_boss.alph_size;
                if (boss->tighten_range(&first, &last, s)) {
                    auto next_base_last = base_boss.fwd(base_range.second, s);
                    auto next_base_first = base_boss.pred_last(next_base_last - 1) + 1;
                    edge_stack.emplace_back(first, last, length);
                    traversal_stack.emplace_back(
                        std::make_pair(next_base_first, next_base_last), length);
                }
            }
        } else if (first == last) {
            if (TAlphabet s = boss->get_W(last)) {
                s %= boss->alph_size;
                if (base_boss.tighten_range(&base_range.first, &base_range.second, s)) {
                    auto next_last = boss->fwd(last, s);
                    auto next_first = boss->pred_last(next_last - 1) + 1;
                    edge_stack.emplace_back(next_first, next_last, length);
                    traversal_stack.emplace_back(base_range, length);
                }
            }
        } else {
            for (TAlphabet s = 1; s < boss->alph_size; ++s) {
                edge_index next_first = first;
                edge_index next_last = last;

                if (boss->tighten_range(&next_first, &next_last, s)) {
                    auto [next_base_first, next_base_last] = base_range;
                    if (base_boss.tighten_range(&next_base_first, &next_base_last, s)) {
                        edge_stack.emplace_back(next_first, next_last, length);
                        traversal_stack.emplace_back(
                            std::make_pair(next_base_first, next_base_last), length);
                    }
                }
            }
        }
    }

    logger->trace("Found {} overlapping ranges, totalling {} k-mers",
                  total_ranges, total_matches);
    return std::make_pair(std::move(boss), std::move(matching_ranges));
}

template <class Seeder, class Extender, class AlignmentCompare>
auto DBGAligner<Seeder, Extender, AlignmentCompare>
::build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                const std::vector<AlignmentResults> &wrapped_seqs) const -> BatchSeeders {
    assert(seq_batch.size() == wrapped_seqs.size());

    BatchSeeders result;
    result.reserve(seq_batch.size());
    size_t batch_size = 0;
    for (const auto &[header, query] : seq_batch) {
        batch_size += query.size();
    }

    auto [sample_boss, matching_ranges]
        = extract_subgraph(graph_, seq_batch, config_, batch_size * 16);

    ProgressBar progress_bar(seq_batch.size(), "Seeding sequences",
                             std::cerr, !common::get_verbose());
    for (size_t i = 0; i < seq_batch.size(); ++i, ++progress_bar) {
        const auto &[header, query] = seq_batch[i];
        std::string_view this_query = wrapped_seqs[i].get_query(false);
        assert(this_query == query);

        std::vector<node_index> nodes;
        if (config_.max_seed_length >= graph_.get_k()) {
            nodes = map_to_nodes_sequentially(graph_, query);
        } else if (this_query.size() >= graph_.get_k()) {
            nodes.resize(this_query.size() - graph_.get_k() + 1);
        }

        std::shared_ptr<ISeeder> seeder;
        if constexpr(std::is_base_of_v<ISuffixSeeder, Seeder>) {
            seeder = std::make_shared<Seeder>(graph_, this_query, false,
                                              std::vector<node_index>(nodes), config_,
                                              sample_boss, matching_ranges);
        } else {
            seeder = std::make_shared<Seeder>(graph_, this_query, false,
                                              std::vector<node_index>(nodes), config_);
        }

        if (this_query.size() * config_.min_exact_match > seeder->get_num_matches())
            seeder = std::make_shared<ManualMatchingSeeder>(std::vector<Seed>{}, 0, config_);

        std::shared_ptr<ISeeder> seeder_rc;
        std::vector<node_index> nodes_rc;

#if ! _PROTEIN_GRAPH
        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            nodes_rc = nodes;
            std::string dummy(query);
            if (config_.max_seed_length >= graph_.get_k()) {
                reverse_complement_seq_path(graph_, dummy, nodes_rc);
                assert(dummy == wrapped_seqs[i].get_query(true));
            }

            assert(nodes_rc.size() == nodes.size());

            std::string_view reverse = wrapped_seqs[i].get_query(true);

            if constexpr(std::is_base_of_v<ISuffixSeeder, Seeder>) {
                seeder_rc = std::make_shared<Seeder>(graph_, reverse, true,
                                                     std::move(nodes_rc), config_,
                                                     sample_boss, matching_ranges);
            } else {
                seeder_rc = std::make_shared<Seeder>(graph_, reverse, true,
                                                     std::move(nodes_rc), config_);
            }
            if (reverse.size() * config_.min_exact_match > seeder_rc->get_num_matches())
                seeder_rc = std::make_shared<ManualMatchingSeeder>(std::vector<Seed>{}, 0, config_);
        }
#endif
        result.emplace_back(std::move(seeder), std::move(seeder_rc));
    }

    return result;
}

template <class Seeder, class Extender, class AlignmentCompare>
void DBGAligner<Seeder, Extender, AlignmentCompare>
::align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
              const AlignmentCallback &callback,
              size_t first_seq_offset) const {
    std::vector<AlignmentResults> paths;
    paths.reserve(seq_batch.size());
    for (const auto &[header, query] : seq_batch) {
        paths.emplace_back(query);
    }

    auto seeders = build_seeders(seq_batch, paths);
    if (!config_.chain_alignments) {
        size_t old_seed_count = 0;
        for (const auto &[seeder, seeder_rc] : seeders) {
            if (seeder)
                old_seed_count += seeder->get_seeds().size();

            if (seeder_rc)
                old_seed_count += seeder_rc->get_seeds().size();
        }
        size_t new_seed_count = cluster_seeds(*this, seeders, old_seed_count);
        logger->trace("Seed count:\tbefore clustering: {}\tafter clustering: {}",
                      old_seed_count, new_seed_count);
    }

    assert(seeders.size() == seq_batch.size());

    for (size_t i = 0; i < seq_batch.size(); ++i) {
        const auto &[header, query] = seq_batch[i];
        auto &[seeder, seeder_rc] = seeders[i];
        AlignmentAggregator<AlignmentCompare> aggregator(config_);

        size_t num_seeds = 0;
        size_t num_explored_nodes = 0;
        size_t num_extensions = 0;

        auto add_alignment = [&](Alignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            aggregator.add_alignment(std::move(alignment));
        };

        auto get_min_path_score = [&](const Alignment &seed) {
            return std::max(config_.min_path_score,
                            seed.has_annotation()
                                ? aggregator.get_score_cutoff(seed.get_columns())
                                : aggregator.get_global_cutoff());
        };

        std::string_view this_query = paths[i].get_query(false);
        assert(this_query == query);

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            std::string_view reverse = paths[i].get_query(true);
            std::tie(num_seeds, num_extensions, num_explored_nodes) =
                align_both_directions(this_query, reverse, *seeder, *seeder_rc,
                                      add_alignment, get_min_path_score);

        } else {
            std::tie(num_seeds, num_extensions, num_explored_nodes) =
                align_core(*seeder, Extender(*this, this_query), add_alignment,
                           get_min_path_score, false);
        }
#else
        if (config_.chain_alignments) {
            std::string_view reverse = paths[i].get_query(true);
            std::tie(num_seeds, num_extensions, num_explored_nodes) =
                align_both_directions(this_query, reverse, *seeder, *seeder_rc,
                                      add_alignment, get_min_path_score);
        } else {
            std::tie(num_seeds, num_extensions, num_explored_nodes) =
                align_core(*seeder, Extender(*this, this_query), add_alignment,
                           get_min_path_score, false);
        }
#endif

        size_t aligned_labels = aggregator.num_aligned_labels();
        score_t best_score = std::numeric_limits<score_t>::min();
        size_t query_coverage = 0;

        auto alignments = aggregator.get_alignments();
        if (config_.post_chain_alignments) {
            bool &post_chain = const_cast<bool&>(config_.post_chain_alignments);
            bool &allow_left_trim = const_cast<bool&>(config_.allow_left_trim);
            bool old_allow_left_trim = allow_left_trim;
            allow_left_trim = false;
            post_chain = false;
            size_t n_ext;
            size_t n_exp;
            std::tie(alignments, n_ext, n_exp)
                = chain_alignments<AlignmentCompare>(*this, std::move(alignments),
                                                     paths[i].get_query(false),
                                                     paths[i].get_query(true));
            post_chain = true;
            allow_left_trim = old_allow_left_trim;
            num_extensions += n_ext;
            num_explored_nodes += n_exp;
        }

        for (auto&& alignment : alignments) {
            assert(alignment.is_valid(graph_, &config_));
            if (alignment.get_score() > best_score) {
                best_score = alignment.get_score();
                query_coverage = alignment.get_cigar().get_coverage();
            }
            paths[i].emplace_back(std::move(alignment));
        }

        double explored_nodes_d = num_explored_nodes;
        double explored_nodes_per_kmer =
            explored_nodes_d / (query.size() - graph_.get_k() + 1);
        if (common::get_verbose()) {
            std::string label_trace;
            if (dynamic_cast<const ILabeledAligner*>(this)) {
                label_trace = fmt::format("\tlabels: {}\texplored nodes/k-mer/label: {:.2f}",
                                          aligned_labels,
                                          aligned_labels ? explored_nodes_per_kmer / aligned_labels : 0);
            }

            logger->trace("{}\t{}\tlength: {}\tcovered: {}\tbest score: {}\tseeds: {}\t"
                    "extensions: {}\texplored nodes: {}\texplored nodes/extension: {:.2f}\t"
                    "explored nodes/k-mer: {:.2f}{}",
                    first_seq_offset + i, header, query.size(), query_coverage,
                    best_score, num_seeds, num_extensions, num_explored_nodes,
                    num_explored_nodes ? explored_nodes_d / num_extensions : 0,
                    explored_nodes_per_kmer, label_trace);
        }

        callback(header, std::move(paths[i]));
    };
}

// Generates seeds and extends them. If force_fixed_seed is true, then
// all alignments must have the seed as a prefix. Otherwise, only the first
// node of the seed is used as an alignment starting node.
template <class Seeder, class Extender>
std::tuple<size_t, size_t, size_t>
align_core(const Seeder &seeder,
           Extender&& extender,
           const std::function<void(Alignment&&)> &callback,
           const std::function<score_t(const Alignment&)> &get_min_path_score,
           bool force_fixed_seed) {
    auto seeds = seeder.get_alignments();

    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty())
            continue;

        score_t min_path_score = get_min_path_score(seeds[i]);

        for (auto&& extension : extender.get_extensions(seeds[i], min_path_score,
                                                        force_fixed_seed)) {
            callback(std::move(extension));
        }

        for (size_t j = i + 1; j < seeds.size(); ++j) {
            if (seeds[j].size() && !extender.check_seed(seeds[j]))
                filter_seed(seeds[i], seeds[j]);
        }
    }

    return { seeds.size(), extender.num_extensions(), extender.num_explored_nodes() };
}

// Construct a full alignment from a chain by aligning the query agaisnt
// the graph in the regions of the query in between the chain seeds.
template <class Seeder, class Extender, class AlignmentCompare>
void DBGAligner<Seeder, Extender, AlignmentCompare>
::extend_chain(std::string_view query,
               std::string_view query_rc,
               Chain&& chain,
               size_t &num_extensions,
               size_t &num_explored_nodes,
               const std::function<void(Alignment&&)> &callback,
               bool extend_ends) const {
    Extender extender(*this, query);
    assert(chain.size());

    Alignment cur = std::move(chain[0].first);
    assert(cur.is_valid(graph_, &config_));
    Alignment best = cur;

    std::vector<Alignment> partial_alignments;
    std::vector<int64_t> coord_offsets { 0 };

    int64_t coord_offset = 0;
    for (size_t i = 1; i < chain.size(); ++i) {
        assert(chain[i].first.is_valid(graph_, &config_));
        coord_offset += chain[i].second;
        assert(coord_offset > 0);
        if (!align_connect(graph_, config_, cur, chain[i].first, coord_offset,
                           extender, partial_alignments)) {
            if (AlignmentCompare()(best, partial_alignments.back()))
                best = partial_alignments.back();

            coord_offsets.push_back(coord_offset);
            coord_offset = 0;
        }

        if (AlignmentCompare()(best, cur))
            best = cur;
    }
    num_extensions += chain.size() - 1;

    if (partial_alignments.size()) {
        partial_alignments.emplace_back(cur);
        assert(partial_alignments.size() == coord_offsets.size());
        Alignment *first = &partial_alignments[0];
        DEBUG_LOG("Partial alignments:\n\t{}", fmt::join(partial_alignments, "\n\t"));
        for (size_t i = 1; i < partial_alignments.size(); ++i) {
            Alignment &next = partial_alignments[i];
            if (next.get_sequence().size() < graph_.get_k()) {
                first = &next;
                continue;
            }

            int64_t num_unknown = coord_offsets[i] - first->get_sequence().size();
            if (num_unknown > 0 && next.get_clipping() > num_unknown) {
                auto merged = *first;
                auto next_fixed = next;
                next_fixed.trim_offset();
                merged.splice_with_unknown(
                    std::move(next_fixed),
                    num_unknown,
                    graph_.get_k() - 1,
                    config_
                );
                assert(merged.is_valid(graph_, &config_));
                if (merged.size()) {
                    std::swap(*first, merged);
                } else {
                    first = &next;
                    continue;
                }
            } else {
                first = &next;
                continue;
            }

            if (AlignmentCompare()(best, *first))
                best = *first;
        }

        if (AlignmentCompare()(best, *first))
            best = *first;
    }

    assert(!best.empty());

    if (best.get_end_clipping() && extend_ends) {
        DEBUG_LOG("Extending back");
        auto extensions = extender.get_extensions(best, config_.ninf, true);
        if (extensions.size()
                && extensions[0].get_end_clipping() < best.get_end_clipping()
                && extensions[0].get_score() > best.get_score()) {
            std::swap(best, extensions[0]);
        }
        DEBUG_LOG("done");
    }

    assert(!best.empty());
    best.trim_offset();
    assert(best.is_valid(graph_, &config_));

    // for now, backwards alignment not supported for amino acids
#if ! _PROTEIN_GRAPH

    if (best.get_clipping() && extend_ends) {
        DEBUG_LOG("Extending front");
        RCDBG rc_dbg(std::shared_ptr<const DeBruijnGraph>(
                        std::shared_ptr<const DeBruijnGraph>(), &graph_));
        const DeBruijnGraph &rc_graph = graph_.get_mode() != DeBruijnGraph::CANONICAL
            ? rc_dbg : graph_;

        auto rev = best;
        rev.reverse_complement(rc_graph, query_rc);
        if (rev.size() && rev.get_nodes().back()) {
            assert(rev.get_end_clipping());
            Extender extender_rc(*this, query_rc);
            extender_rc.set_graph(rc_graph);
            auto extensions = extender_rc.get_extensions(rev, config_.ninf, true);
            ++num_extensions;
            num_explored_nodes += extender_rc.num_explored_nodes() - rev.size();
            if (extensions.size()
                    && extensions[0].get_end_clipping() < rev.get_end_clipping()) {
                extensions[0].reverse_complement(rc_graph, query);
                if (extensions[0].size())
                    std::swap(best, extensions[0]);
            }
        }
        DEBUG_LOG("done");
    }

#else

    std::ignore = query_rc;

#endif

    num_extensions += extender.num_extensions();
    num_explored_nodes += extender.num_explored_nodes();

    assert(!best.empty());
    callback(std::move(best));

    for (auto&& aln : partial_alignments) {
        if (!aln.empty())
            callback(std::move(aln));
    }
}

template <class Seeder, class Extender, class AlignmentCompare>
std::tuple<size_t, size_t, size_t>
DBGAligner<Seeder, Extender, AlignmentCompare>
::align_both_directions(std::string_view forward,
                        std::string_view reverse,
                        const ISeeder &forward_seeder,
                        const ISeeder &reverse_seeder,
                        const std::function<void(Alignment&&)> &callback,
                        const std::function<score_t(const Alignment&)> &get_min_path_score) const {
    size_t num_seeds = 0;
    size_t num_extensions = 0;
    size_t num_explored_nodes = 0;

    if (config_.chain_alignments) {
        if (!has_coordinates()) {
            logger->error("Chaining only supported for seeds with coordinates. Skipping seed chaining.");
            exit(1);
        }

        auto fwd_seeds = forward_seeder.get_seeds();

#if ! _PROTEIN_GRAPH
        auto bwd_seeds = reverse_seeder.get_seeds();
#else
        std::vector<Seed> bwd_seeds;
        std::ignore = reverse_seeder;
#endif

        if (fwd_seeds.empty() && bwd_seeds.empty())
            return std::make_tuple(num_seeds, num_extensions, num_explored_nodes);

        AlignmentAggregator<AlignmentCompare> aggregator(config_);
        tsl::hopscotch_set<Alignment::Column> all_columns;
        for (const auto &seed : fwd_seeds) {
            const auto &columns = seed.get_columns();
            all_columns.insert(columns.begin(), columns.end());
        }
        for (const auto &seed : bwd_seeds) {
            const auto &columns = seed.get_columns();
            all_columns.insert(columns.begin(), columns.end());
        }

        try {
            tsl::hopscotch_set<Alignment::Column> finished_columns;
            size_t this_num_explored;
            std::tie(num_seeds, this_num_explored) = call_seed_chains_both_strands(
                forward, reverse, config_, std::move(fwd_seeds), std::move(bwd_seeds),
                [&](Chain&& chain, score_t score) {
                    if (config_.num_alternative_paths <= 1
                            && finished_columns.size() == all_columns.size()) {
                        throw std::bad_function_call();
                    }

                    double exact_match_fraction
                        = static_cast<double>(get_num_char_matches_in_seeds(chain.begin(), chain.end()))
                            / forward.size();

                    logger->trace("Chain: score: {}; exact match fraction: {}",
                                  score, exact_match_fraction);

#ifndef NDEBUG
                    for (const auto &[chain, dist] : chain) {
                        DEBUG_LOG("\t{}\tdist: {}", chain, dist);
                    }
#endif

                    if (exact_match_fraction < config_.min_exact_match)
                        throw std::bad_function_call();

                    extend_chain(chain[0].first.get_orientation() ? reverse : forward,
                                 chain[0].first.get_orientation() ? forward : reverse,
                                 std::move(chain), num_extensions, num_explored_nodes,
                                 [&](Alignment&& aln) {
                        const auto &cur_columns = aln.get_columns();
                        if (!aggregator.add_alignment(std::move(aln))) {
                            finished_columns.insert(cur_columns.begin(), cur_columns.end());
                        }
                    });
                },
                [&](Alignment::Column column) { return finished_columns.count(column); }
            );
            num_explored_nodes += this_num_explored;
        } catch (const std::bad_function_call&) {}

        for (Alignment &alignment : aggregator.get_alignments()) {
            if (alignment.get_score() < get_min_path_score(alignment))
                continue;

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL && alignment.get_orientation()) {
                Alignment rev(alignment);
                rev.reverse_complement(graph_, forward);
                if (rev.size())
                    std::swap(rev, alignment);
            }

            callback(std::move(alignment));
        }

        return std::make_tuple(num_seeds, num_extensions, num_explored_nodes);
    }

    Extender forward_extender(*this, forward);
    Extender reverse_extender(*this, reverse);

#if ! _PROTEIN_GRAPH

    auto fwd_seeds = forward_seeder.get_alignments();
    auto bwd_seeds = reverse_seeder.get_alignments();

    RCDBG rc_dbg(std::shared_ptr<const DeBruijnGraph>(
                    std::shared_ptr<const DeBruijnGraph>(), &graph_));
    bool use_rcdbg = graph_.get_mode() != DeBruijnGraph::CANONICAL
                        && config_.forward_and_reverse_complement;

    const DeBruijnGraph &rc_graph = use_rcdbg ? rc_dbg : graph_;

    auto is_reversible = [this](const Alignment &alignment) {
        return graph_.get_mode() == DeBruijnGraph::CANONICAL
            && alignment.get_orientation()
            && !alignment.get_offset();
    };

    auto aln_both = [&](std::string_view query,
                        std::string_view query_rc,
                        std::vector<Alignment>&& seeds,
                        Extender &fwd_extender,
                        Extender &bwd_extender,
                        const std::function<void(Alignment&&)> &callback) {
        fwd_extender.set_graph(graph_);
        bwd_extender.set_graph(rc_graph);
        num_seeds += seeds.size();

        if (seeds.empty())
            return;

        for (size_t i = 0; i < seeds.size(); ++i) {
            if (seeds[i].empty())
                continue;

            score_t min_path_score = config_.min_cell_score;
            auto extensions = fwd_extender.get_extensions(seeds[i], min_path_score, false);

            std::vector<Alignment> rc_of_alignments;
            for (Alignment &path : extensions) {
                if (path.get_score() >= get_min_path_score(path)) {
                    if (is_reversible(path)) {
                        Alignment out_path = path;
                        out_path.reverse_complement(graph_, query_rc);
                        assert(out_path.size());
                        callback(std::move(out_path));
                    } else {
                        callback(Alignment(path));
                    }
                }

                if (!path.get_clipping() || path.get_offset())
                    continue;

                path.reverse_complement(rc_graph, query_rc);

                if (path.empty()) {
                    DEBUG_LOG("This local alignment cannot be reversed, skipping");
                    continue;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(path.get_end_clipping());
                assert(path.is_valid(rc_graph, &config_));

                rc_of_alignments.emplace_back(std::move(path));
            }

            align_core(ManualSeeder(std::move(rc_of_alignments)), bwd_extender,
                [&](Alignment&& path) {
                    if (use_rcdbg || is_reversible(path)) {
                        path.reverse_complement(rc_graph, query);
                        if (path.empty())
                            return;

                        if (auto *filter = dynamic_cast<SeedFilteringExtender*>(&fwd_extender)) {
                            for (node_index node : path.get_nodes()) {
                                filter->filter_nodes(node, path.get_clipping(),
                                                     query.size() - path.get_end_clipping());
                            }
                        }
                    }

                    assert(path.is_valid(graph_, &config_));

                    callback(std::move(path));
                },
                get_min_path_score,
                true /* alignments must have the seed as a prefix */
            );

            for (size_t j = i + 1; j < seeds.size(); ++j) {
                if (seeds[j].size() && !fwd_extender.check_seed(seeds[j]))
                    filter_seed(seeds[i], seeds[j]);
            }
        }
    };

    size_t fwd_num_matches = forward_seeder.get_num_matches();
    size_t bwd_num_matches = reverse_seeder.get_num_matches();

    if (fwd_num_matches >= bwd_num_matches) {
        aln_both(forward, reverse, std::move(fwd_seeds),
                 forward_extender, reverse_extender, callback);
        if (bwd_num_matches >= fwd_num_matches * config_.rel_score_cutoff) {
            aln_both(reverse, forward, std::move(bwd_seeds),
                     reverse_extender, forward_extender, callback);
        }
    } else {
        aln_both(reverse, forward, std::move(bwd_seeds),
                 reverse_extender, forward_extender, callback);
        if (fwd_num_matches >= bwd_num_matches * config_.rel_score_cutoff) {
            aln_both(forward, reverse, std::move(fwd_seeds),
                     forward_extender, reverse_extender, callback);
        }
    }

    return std::make_tuple(
        num_seeds,
        num_extensions + forward_extender.num_extensions()
                        + reverse_extender.num_extensions(),
        num_explored_nodes + forward_extender.num_explored_nodes()
                            + reverse_extender.num_explored_nodes()
    );

#else

    throw std::runtime_error("Reverse complements not defined for amino acids");

#endif
}

template class DBGAligner<>;
template class DBGAligner<SuffixSeeder<ExactSeeder>, LabeledExtender>;

} // namespace align
} // namespace graph
} // namespace mtg
