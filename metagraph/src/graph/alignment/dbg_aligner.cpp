#include "dbg_aligner.hpp"

#include "aligner_aggregator.hpp"

namespace mtg {
namespace graph {
namespace align {

IDBGAligner::DBGQueryAlignment IDBGAligner::align(std::string_view query,
                                                  bool is_reverse_complement) const {
    DBGQueryAlignment result(query);
    std::string empty_header;
    align_batch(
        [&](const QueryCallback &callback) {
            callback(empty_header, query, is_reverse_complement);
        },
        [&](std::string_view, DBGQueryAlignment&& alignment) {
            result = std::move(alignment);
        }
    );

    return result;
}

void IDBGAligner
::align_batch(const std::vector<std::pair<std::string, std::string>> &seq_batch,
              const AlignmentCallback &callback) const {
    align_batch([&](const QueryCallback &query_callback) {
        for (const auto &[header, seq] : seq_batch) {
            query_callback(header, seq, false /* orientation of seq */);
        }
    }, callback);
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_core(std::string_view query,
             const ISeeder<node_index> &seeder,
             IExtender<node_index> &extender,
             const LocalAlignmentCallback &callback,
             const MinScoreComputer &get_min_path_score) const {
    bool filter_seeds = dynamic_cast<const ExactSeeder<node_index>*>(&seeder);

    std::vector<DBGAlignment> seeds = seeder.get_seeds();
    std::sort(seeds.begin(), seeds.end(), LocalAlignmentLess());

    for (DBGAlignment &seed : seeds) {
        score_t min_path_score = get_min_path_score(seed);

        // check if this seed has been explored before in an alignment and discard
        // it if so
        if (filter_seeds) {
            size_t found_count = 0;
            std::pair<size_t, size_t> idx_range {
                seed.get_clipping(),
                seed.get_clipping() + graph_.get_k() - seed.get_offset()
            };
            for (node_index node : seed) {
                auto emplace = visited_nodes_.emplace(node, idx_range);
                auto &range = emplace.first.value();
                if (emplace.second) {
                } else if (range.first > idx_range.first || range.second < idx_range.second) {
                    DEBUG_LOG("Node: {}; Prev_range: [{},{})", node, range.first, range.second);
                    range.first = std::min(range.first, idx_range.first);
                    range.second = std::max(range.second, idx_range.second);
                    DEBUG_LOG("Node: {}; cur_range: [{},{})", node, range.first, range.second);
                } else {
                    ++found_count;
                }

                if (idx_range.second - idx_range.first == graph_.get_k())
                    ++idx_range.first;

                ++idx_range.second;
            }

            if (found_count == seed.size()) {
                DEBUG_LOG("Skipping seed: {}", seed);
                continue;
            }
        }

        if (seed.get_query().data() + seed.get_query().size()
                == query.data() + query.size()) {
            if (seed.get_score() >= min_path_score) {
                seed.trim_offset();
                assert(seed.is_valid(graph_, &config_));
                DEBUG_LOG("Alignment: {}", seed);
                callback(std::move(seed));
            }

            continue;
        }

        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

        extender.initialize(seed);
        auto extensions = extender.get_extensions(min_path_score);

        // if the ManualSeeder is not used, then add nodes to the visited_nodes_
        // table to allow for seed filtration
        if (filter_seeds) {
            extender.call_visited_nodes([&](node_index node, size_t begin, size_t end) {
                auto emplace = visited_nodes_.emplace(node, std::make_pair(begin, end));
                auto &range = emplace.first.value();
                if (!emplace.second) {
                    range.first = std::min(range.first, begin);
                    range.second = std::max(range.second, end);
                }
            });
        }

        if (extensions.empty() && seed.get_score() >= min_path_score) {
            seed.extend_query_end(query.data() + query.size());
            seed.trim_offset();
            assert(seed.is_valid(graph_, &config_));
            DEBUG_LOG("Alignment (seed): {}", seed);
            callback(std::move(seed));
        }

        for (auto&& extension : extensions) {
            assert(extension.is_valid(graph_, &config_));
            callback(std::move(extension));
        }
    }
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_one_direction(DBGQueryAlignment &paths,
                      bool orientation_to_align,
                      const ISeeder<node_index> &seeder,
                      IExtender<node_index>&& extender) const {
    std::string_view query = paths.get_query(orientation_to_align);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(query, seeder, extender, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_best_direction(DBGQueryAlignment &paths,
                       const ISeeder<node_index> &seeder,
                       const ISeeder<node_index> &seeder_rc,
                       IExtender<node_index>&& extender,
                       IExtender<node_index>&& extender_rc) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query(true);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(forward, seeder, extender, alignment_callback, get_min_path_score);
        align_core(reverse, seeder_rc, extender_rc, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_both_directions(DBGQueryAlignment &paths,
                        const ISeeder<node_index> &forward_seeder,
                        const ISeeder<node_index> &reverse_seeder,
                        IExtender<node_index>&& forward_extender,
                        IExtender<node_index>&& reverse_extender,
                        const AlignCoreGenerator &rev_comp_core_generator) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query(true);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        auto get_forward_alignments = [&](std::string_view query,
                                          std::string_view query_rc,
                                          const ISeeder<node_index> &seeder,
                                          IExtender<node_index> &extender) {
            std::vector<DBGAlignment> rc_of_alignments;

            DEBUG_LOG("Extending in forwards direction");
            align_core(query, seeder, extender, [&](DBGAlignment&& path) {
                score_t min_path_score = get_min_path_score(path);

                if (path.get_score() >= min_path_score)
                    alignment_callback(DBGAlignment(path));

                if (!path.get_clipping())
                    return;

                auto rev = path;
                rev.reverse_complement(graph_, query_rc);
                if (rev.empty()) {
                    DEBUG_LOG("Alignment cannot be reversed, returning");
                    if (path.get_score() >= min_path_score)
                        alignment_callback(std::move(path));

                    return;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(rev.get_end_clipping());
                rev.trim_end_clipping();
                assert(rev.is_valid(graph_, &config_));

                // Pass the reverse complement of the forward alignment
                // as a seed for extension
                rc_of_alignments.emplace_back(std::move(rev));
            }, [&](const auto&) { return config_.min_cell_score; });

            return rc_of_alignments;
        };

        std::vector<DBGAlignment> rc_of_reverse = get_forward_alignments(
            reverse, forward, reverse_seeder, reverse_extender
        );
        std::vector<DBGAlignment> rc_of_forward = get_forward_alignments(
            forward, reverse, forward_seeder, forward_extender
        );

        auto extend_reverse = [&](std::string_view query_rc,
                                  const ISeeder<node_index> &seeder,
                                  std::vector<DBGAlignment>&& rc_of_alignments) {
            DEBUG_LOG("Extending in reverse direction");
            rev_comp_core_generator(query_rc, seeder, std::move(rc_of_alignments),
                                    [&](const auto &seeder_rc, auto&& extender_rc) {
                align_core(query_rc, seeder_rc, extender_rc,
                           alignment_callback, get_min_path_score);
            });
        };

        extend_reverse(forward, reverse_seeder, std::move(rc_of_reverse));
        extend_reverse(reverse, forward_seeder, std::move(rc_of_forward));
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_aggregate(DBGQueryAlignment &paths,
                  const AlignmentGenerator &alignment_generator) const {
    AlignmentAggregator<node_index, AlignmentCompare> path_queue(
        paths.get_query() /* forward */,
        paths.get_query(true) /* reverse complement */,
        config_
    );

    alignment_generator(
        [&](DBGAlignment&& alignment) { path_queue.add_alignment(std::move(alignment)); },
        [&](const DBGAlignment &seed) { return path_queue.get_min_path_score(seed); }
    );

    path_queue.call_alignments([&](auto&& alignment) {
        assert(alignment.is_valid(graph_, &config_));
        paths.emplace_back(std::move(alignment));
    });
}

template class SeedAndExtendAlignerCore<>;

} // namespace align
} // namespace graph
} // namespace mtg
