#include "dbg_aligner.hpp"

#include <progress_bar.hpp>

#include "aligner_aggregator.hpp"
#include "aligner_labeled.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/graph_extensions/path_index.hpp"

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
    }
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

        std::shared_ptr<ISeeder> seeder
            = std::make_shared<Seeder>(graph_, this_query, false,
                                       std::vector<node_index>(nodes), config_);

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
            seeder_rc = std::make_shared<Seeder>(graph_, reverse, true,
                                                 std::move(nodes_rc), config_);
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
    assert(seeders.size() == seq_batch.size());

    for (size_t i = 0; i < seq_batch.size(); ++i) {
        const auto &[header, query] = seq_batch[i];
        auto &[seeder, seeder_rc] = seeders[i];
        AlignmentAggregator<AlignmentCompare> aggregator(config_);

        size_t num_seeds = 0;
        size_t num_explored_nodes = 0;
        size_t num_extensions = 0;

        auto get_min_path_score = [&](const Alignment &seed) {
            return std::max(config_.min_path_score,
                            seed.has_annotation()
                                ? aggregator.get_score_cutoff(seed.get_columns())
                                : aggregator.get_global_cutoff());
        };

        auto add_alignment = [&](Alignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            if (config_.allow_jump || alignment.get_score() >= get_min_path_score(alignment)) {
                aggregator.add_alignment(std::move(alignment));
            }
        };

        std::string_view this_query = paths[i].get_query(false);
        std::string_view reverse = paths[i].get_query(true);
        assert(this_query == query);

        if (!config_.chain_alignments
                && graph_.get_extension_threadsafe<IPathIndex>()) {
            if (graph_.get_mode() != DeBruijnGraph::BASIC)
                seeder_rc = std::make_shared<ManualSeeder>();

            auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                chain_and_filter_seeds(*this, seeder, Extender(*this, this_query),
                    Extender(*this, reverse));
            num_seeds += num_seeds_c;
            num_extensions += num_extensions_c;
            num_explored_nodes += num_explored_nodes_c;
            auto alignments = seeder->get_alignments();
            std::for_each(std::make_move_iterator(alignments.begin()),
                          std::make_move_iterator(alignments.end()),
                          add_alignment);
            seeder = std::make_shared<ManualSeeder>();
            if (seeder_rc) {
                auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                    chain_and_filter_seeds(*this, seeder_rc, Extender(*this, reverse),
                        Extender(*this, this_query));
                num_seeds += num_seeds_c;
                num_extensions += num_extensions_c;
                num_explored_nodes += num_explored_nodes_c;
                auto alignments = seeder_rc->get_alignments();
                std::for_each(std::make_move_iterator(alignments.begin()),
                              std::make_move_iterator(alignments.end()),
                              add_alignment);
                seeder_rc = std::make_shared<ManualSeeder>();
            }
        }

#if ! _PROTEIN_GRAPH
        if (seeder_rc) {
            auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                align_both_directions(this_query, reverse, *seeder, *seeder_rc,
                                      add_alignment, get_min_path_score);
            num_seeds += num_seeds_c;
            num_extensions += num_extensions_c;
            num_explored_nodes += num_explored_nodes_c;

        } else {
            auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                align_core(*seeder, Extender(*this, this_query), add_alignment,
                           get_min_path_score, false);
            num_seeds += num_seeds_c;
            num_extensions += num_extensions_c;
            num_explored_nodes += num_explored_nodes_c;
        }
#else
        if (config_.chain_alignments) {
            std::string_view reverse = paths[i].get_query(true);
            auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                align_both_directions(this_query, reverse, *seeder, *seeder_rc,
                                      add_alignment, get_min_path_score);
            num_seeds += num_seeds_c;
            num_extensions += num_extensions_c;
            num_explored_nodes += num_explored_nodes_c;
        } else {
            auto [num_seeds_c, num_extensions_c, num_explored_nodes_c] =
                align_core(*seeder, Extender(*this, this_query), add_alignment,
                           get_min_path_score, false);
            num_seeds += num_seeds_c;
            num_extensions += num_extensions_c;
            num_explored_nodes += num_explored_nodes_c;
        }
#endif

        size_t aligned_labels = aggregator.num_aligned_labels();
        score_t best_score = std::numeric_limits<score_t>::min();
        size_t query_coverage = 0;

        auto alignments = aggregator.get_alignments();
        if (alignments.size() && (config_.allow_jump || config_.allow_label_change)) {
            std::vector<Alignment> chained_alignments;
            auto add_aln = [&](auto&& alignment) {
                if (chained_alignments.empty()) {
                    chained_alignments.emplace_back(std::move(alignment));
                    return;
                }

                if (alignment.get_score() > chained_alignments[0].get_score())
                    chained_alignments.clear();

                chained_alignments.emplace_back(std::move(alignment));
            };

            std::vector<Alignment> alns;
            alns.reserve(alignments.size());
            for (const Alignment &aln : alignments) {
                if (!aln.get_orientation())
                    alns.emplace_back(aln);
            }
            for (const Alignment &aln : alignments) {
                if (aln.get_orientation())
                    alns.emplace_back(aln);
            }
            assert(alns.size() == alignments.size());
            alignments.resize(1);

            chain_alignments(*this, alns, add_aln);

            if (chained_alignments.size()
                    && chained_alignments[0].get_score() > alignments[0].get_score()) {
                std::swap(chained_alignments, alignments);
            }
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

        extender.extend_seed_end(seeds[i], callback, force_fixed_seed,
                                 get_min_path_score(seeds[i]));

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
::extend_chain(Chain&& chain,
               SeedFilteringExtender &extender,
               const std::function<void(Alignment&&)> &callback,
               bool try_connect) const {
    assert(chain.size());

    Alignment cur = std::move(chain[0].first);
    assert(cur.size());
    assert(cur.is_valid(graph_, &config_));
    Alignment best = cur;

    std::vector<Alignment> partial_alignments;
    std::vector<bool> is_disconnected { false };
    std::vector<int64_t> coord_offsets { 0 };

    auto is_coord_jump = [&](const Alignment &a, const Alignment &b) {
        auto b_first_end = b.get_query_view().begin() + graph_.get_k() - b.get_offset();
        return a.get_query_view().end() == b_first_end;
    };

    int64_t coord_offset = 0;
    for (size_t i = 1; i < chain.size(); ++i) {
        assert(chain[i].first.is_valid(graph_, &config_));
        std::vector<Alignment> alignments;

        if (chain[i].second >= std::numeric_limits<uint32_t>::max()) {
            // connected disjoint regions
            coord_offset += chain[i].second - std::numeric_limits<uint32_t>::max();
        } else {
            coord_offset += chain[i].second;

            if (!is_coord_jump(cur, chain[i].first)) {
                assert(coord_offset > 0);
                // config_.allow_label_change = (chain[i].first.label_columns != chain[i - 1].first.label_columns);
                alignments = extender.connect_seeds(cur, chain[i].first, coord_offset);
            }
        }

        if (alignments.size()) {
            // TODO: what if there are multiple alignments?
            assert(alignments[0].size());
            std::swap(alignments[0], cur);
            if (AlignmentCompare()(best, cur))
                best = cur;
        } else {
            DEBUG_LOG("Extension not found, restarting from seed {}", chain[i].first);
            partial_alignments.emplace_back(cur);
            is_disconnected.emplace_back(chain[i].second >= std::numeric_limits<uint32_t>::max());
            std::swap(cur, chain[i].first);
            if (AlignmentCompare()(best, partial_alignments.back()))
                best = partial_alignments.back();

            coord_offsets.push_back(coord_offset);
            coord_offset = 0;
        }
        assert(cur.size());
        assert(cur.is_valid(graph_, &config_));
    }

    if (try_connect && partial_alignments.size()) {
        partial_alignments.emplace_back(cur);
        is_disconnected.emplace_back(chain.back().second >= std::numeric_limits<uint32_t>::max());
        assert(partial_alignments.size() == coord_offsets.size());
        Alignment *first = &partial_alignments[0];
        DEBUG_LOG("Partial alignments:\n\t{}", fmt::join(partial_alignments, "\n\t"));
        for (size_t i = 1; i < partial_alignments.size(); ++i) {
            Alignment &next = partial_alignments[i];
            if (next.get_sequence().size() < graph_.get_k()) {
                first = &next;
                continue;
            }

            // int64_t num_unknown = coord_offsets[i] - first->get_sequence().size();
            ssize_t overlap = first->get_query_view().end() - next.get_query_view().begin();
            if (is_coord_jump(*first, next) || is_disconnected[i]) {
                auto merged = *first;
                auto next_fixed = next;
                if (overlap > 0)
                    next_fixed.trim_query_prefix(overlap, graph_.get_k() - 1, config_, false);

                next_fixed.insert_gap_prefix(-overlap, graph_.get_k() - 1, config_);
                merged.splice(std::move(next_fixed));
                assert(merged.is_valid(graph_, &config_));
                if (merged.size()) {
                    std::swap(*first, merged);
                } else {
                    first = &next;
                    continue;
                }
            // } else if (num_unknown > 0 && next.get_clipping() > num_unknown && overlap <= 0) {
            //     auto merged = *first;
            //     auto next_fixed = next;
            //     next_fixed.trim_offset();
            //     assert(merged.is_valid(graph_, &config_));
            //     assert(next_fixed.is_valid(graph_, &config_));
            //     merged.splice_with_unknown(std::move(next_fixed), num_unknown,
            //                                graph_.get_k() - 1, config_);
            //     assert(merged.is_valid(graph_, &config_));
            //     if (merged.size()) {
            //         std::swap(*first, merged);
            //     } else {
            //         first = &next;
            //         continue;
            //     }
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
    callback(std::move(best));

    if (!try_connect) {
        for (auto&& aln : partial_alignments) {
            if (!aln.empty())
                callback(std::move(aln));
        }
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

    std::vector<Alignment> fwd_alignments;
    std::vector<Alignment> bwd_alignments;
    size_t fwd_num_matches = 0;
    size_t bwd_num_matches = 0;

    Extender forward_extender(*this, forward);
    Extender reverse_extender(*this, reverse);

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

        tsl::hopscotch_set<Alignment::Column> finished_columns;
        size_t num_nodes;
        std::tie(num_seeds, num_nodes) = call_seed_chains_both_strands(
            *this, forward, reverse, config_, std::move(fwd_seeds), std::move(bwd_seeds),
            [&](Chain&& chain, score_t score) {
                logger->trace("Chain: score: {}", score);

#ifndef NDEBUG
                for (const auto &[chain, dist] : chain) {
                    DEBUG_LOG("\t{}\tdist: {}", chain, dist);
                }
#endif

                extend_chain(std::move(chain),
                             chain[0].first.get_orientation()
                                 ? reverse_extender : forward_extender,
                             [&](Alignment&& aln) {
                    const auto &cur_columns = aln.get_columns();
                    if (!aggregator.add_alignment(std::move(aln))) {
                        finished_columns.insert(cur_columns.begin(), cur_columns.end());
                    }
                });
            },
            [&](Alignment::Column column) { return finished_columns.count(column); },
            [&]() { return config_.num_alternative_paths <= 1
                            && finished_columns.size() == all_columns.size(); }
        );
        num_explored_nodes += num_nodes;

        for (Alignment &alignment : aggregator.get_alignments()) {
            if (alignment.get_score() < get_min_path_score(alignment))
                continue;

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL && alignment.get_orientation()) {
                Alignment rev(alignment);
                rev.reverse_complement(graph_, forward);
                if (rev.size())
                    std::swap(rev, alignment);
            }

#if ! _PROTEIN_GRAPH
            auto &seeds = alignment.get_orientation() ? bwd_alignments : fwd_alignments;
            seeds.emplace_back(std::move(alignment));
#else
            callback(std::move(alignment));
#endif
        }

#if _PROTEIN_GRAPH
        return std::make_tuple(num_seeds,
                               num_extensions + forward_extender.num_extensions()
                                               + reverse_extender.num_extensions(),
                               num_explored_nodes + forward_extender.num_explored_nodes()
                                                   + reverse_extender.num_explored_nodes());
#endif

    } else {
        fwd_alignments = forward_seeder.get_alignments();
        bwd_alignments = reverse_seeder.get_alignments();
        fwd_num_matches = forward_seeder.get_num_matches();
        bwd_num_matches = reverse_seeder.get_num_matches();
    }

#if ! _PROTEIN_GRAPH

    auto aln_both = [&](std::string_view query,
                        std::vector<Alignment>&& seeds,
                        Extender &fwd_extender,
                        Extender &bwd_extender,
                        const std::function<void(Alignment&&)> &callback) {
        num_seeds += seeds.size();

        if (seeds.empty())
            return;

        for (size_t i = 0; i < seeds.size(); ++i) {
            if (seeds[i].empty())
                continue;

            score_t min_path_score = config_.min_cell_score;
            auto extensions = fwd_extender.get_extensions(seeds[i], min_path_score, true);

            std::vector<Alignment> alignments_to_front_extend;
            for (Alignment &path : extensions) {
                if (path.get_score() >= get_min_path_score(path))
                    callback(Alignment(path));

                if (path.get_clipping() && !path.get_offset())
                    alignments_to_front_extend.emplace_back(std::move(path));
            }

            for (const Alignment &seed : alignments_to_front_extend) {
                bwd_extender.rc_extend_rc(seed, [&](Alignment&& aln) {
                    assert(aln.is_valid(graph_, &config_));
                    for (node_index node : aln.get_nodes()) {
                        fwd_extender.filter_nodes(node, aln.get_clipping(),
                                                  query.size() - aln.get_end_clipping());
                    }

                    callback(std::move(aln));
                }, true, get_min_path_score(seed));
            }

            for (size_t j = i + 1; j < seeds.size(); ++j) {
                if (seeds[j].size() && !fwd_extender.check_seed(seeds[j]))
                    filter_seed(seeds[i], seeds[j]);
            }
        }
    };

    if (fwd_num_matches >= bwd_num_matches) {
        aln_both(forward, std::move(fwd_alignments), forward_extender,
                 reverse_extender, callback);
        if (bwd_num_matches >= fwd_num_matches * config_.rel_score_cutoff) {
            aln_both(reverse, std::move(bwd_alignments),
                     reverse_extender, forward_extender, callback);
        }
    } else {
        aln_both(reverse, std::move(bwd_alignments),
                 reverse_extender, forward_extender, callback);
        if (fwd_num_matches >= bwd_num_matches * config_.rel_score_cutoff) {
            aln_both(forward, std::move(fwd_alignments),
                     forward_extender, reverse_extender, callback);
        }
    }

    return std::make_tuple(num_seeds,
                           num_extensions + forward_extender.num_extensions()
                                           + reverse_extender.num_extensions(),
                           num_explored_nodes + forward_extender.num_explored_nodes()
                                               + reverse_extender.num_explored_nodes());

#else

    throw std::runtime_error("Reverse complements not defined for amino acids");

#endif
}

template class DBGAligner<>;
template class DBGAligner<SuffixSeeder<ExactSeeder>, LabeledExtender>;

} // namespace align
} // namespace graph
} // namespace mtg
