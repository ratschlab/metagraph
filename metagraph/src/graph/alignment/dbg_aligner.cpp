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
             IExtender<node_index>&& extender,
             const LocalAlignmentCallback &callback,
             const MinScoreComputer &get_min_path_score) const {
    std::vector<DBGAlignment> seeds;
    seeder.call_seeds([&](DBGAlignment&& seed) {
        assert(seed.is_valid(graph_, &config_));
        seeds.emplace_back(std::move(seed));
    });

    for (auto &seed : seeds) {
        bool inserted = false;
        std::pair<size_t, size_t> idx_range {
            seed.get_clipping(),
            seed.get_clipping() + seed.get_query().size()
        };
        for (node_index node : seed) {
            auto emplace = visited_nodes_.emplace(node, idx_range);
            auto &range = emplace.first.value();
            if (emplace.second) {
                inserted = true;
            } else if (range.first > idx_range.first || range.second < idx_range.second) {
                range.first = std::min(range.first, idx_range.first);
                range.second = std::max(range.second, idx_range.second);
                inserted = true;
            }
        }

        if (!inserted) {
#ifndef NDEBUG
            mtg::common::logger->trace("Skipping seed: {}", seed);
#endif
            continue;
        }

#ifndef NDEBUG
        mtg::common::logger->trace("Seed: {}", seed);
#endif
        score_t min_path_score = get_min_path_score(seed);

        if (seed.get_query_end() == query.data() + query.size()) {
            if (seed.get_score() >= min_path_score) {
                seed.trim_offset();
                assert(seed.is_valid(graph_, &config_));
#ifndef NDEBUG
                mtg::common::logger->trace("Alignment: {}", seed);
#endif
                callback(std::move(seed));
            }

            continue;
        }

        bool extended = false;
        extender.initialize(seed);
        extender([&](DBGAlignment&& extension, auto start_node) {
            if (!start_node && !extended) {
                // no good extension found
                if (seed.get_score() >= min_path_score) {
                    seed.extend_query_end(query.data() + query.size());
                    seed.trim_offset();
                    assert(seed.is_valid(graph_, &config_));
#ifndef NDEBUG
                    mtg::common::logger->trace("Alignment: {}", seed);
#endif
                    callback(std::move(seed));
                }
                extended = true;
                return;
            }

            assert(extension.is_valid(graph_, &config_));
            extension.extend_query_end(query.data() + query.size());

            if (extension.get_clipping() || start_node != seed.back()) {
                // if the extension starts at a different position
                // from the seed end, then it's a new alignment
                extension.extend_query_begin(query.data());
                extension.trim_offset();
                assert(extension.is_valid(graph_, &config_));
#ifndef NDEBUG
                mtg::common::logger->trace("Alignment: {}", extension);
#endif
                callback(std::move(extension));
                return;
            }

            assert(extension.get_offset() == graph_.get_k() - 1);
            auto next_path = seed;
            next_path.append(std::move(extension));
            next_path.trim_offset();
            assert(next_path.is_valid(graph_, &config_));

#ifndef NDEBUG
            mtg::common::logger->trace("Alignment: {}", next_path);
#endif
            callback(std::move(next_path));
            extended = true;
        }, min_path_score);

        // if !extended, then the seed was not extended because of early cutoff

        extender.call_explored_nodes([&](node_index node, size_t begin, size_t end) {
            auto emplace = visited_nodes_.emplace(node, std::make_pair(begin, end));
            auto &range = emplace.first.value();
            if (!emplace.second) {
                range.first = std::min(range.first, begin);
                range.second = std::max(range.second, end);
            }
        });
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
        align_core(query, seeder, std::move(extender),
                   alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_best_direction(DBGQueryAlignment &paths,
                       const ISeeder<node_index> &forward_seeder,
                       const ISeeder<node_index> &reverse_seeder,
                       IExtender<node_index>&& forward_extender,
                       IExtender<node_index>&& reverse_extender) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query(true);

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
        align_core(forward, forward_seeder, std::move(forward_extender),
                   alignment_callback, get_min_path_score);

        align_core(reverse, reverse_seeder, std::move(reverse_extender),
                   alignment_callback, get_min_path_score);

    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_both_directions(DBGQueryAlignment &paths,
                        const ISeeder<node_index> &forward_seeder,
                        IExtender<node_index>&& forward_extender,
                        IExtender<node_index>&& reverse_extender,
                        const SeederGenerator &seeder_generator) const {
    std::string_view forward = paths.get_query();
    std::string_view reverse = paths.get_query(true);

    std::vector<DBGAlignment> reverse_seeds;

    align_aggregate(paths, [&](const auto &alignment_callback,
                               const auto &get_min_path_score) {
#ifndef NDEBUG
        mtg::common::logger->trace("Aligning forwards");
#endif

        // First get forward alignments
        align_core(forward, forward_seeder, std::move(forward_extender),
            [&](DBGAlignment&& path) {
                score_t min_path_score = get_min_path_score(path);

                // If the alignment starts from the beginning of the query,
                // there's no sequence left for aligning backwards.
                if (!path.get_clipping()) {
                    if (path.get_score() >= min_path_score)
                        alignment_callback(std::move(path));

                    return;
                }

                auto rev = path;
                rev.reverse_complement(graph_, reverse);
                if (rev.empty()) {
#ifndef NDEBUG
                    mtg::common::logger->trace("Alignment cannot be reversed, returning");
#endif
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
                reverse_seeds.emplace_back(std::move(rev));
            },
            [&](const auto&) {
                // ignore the min path score for the forward alignment,
                // since it may have a score that is too low before it is
                // extended backwards
                return config_.min_cell_score;
            }
        );

#ifndef NDEBUG
        mtg::common::logger->trace("Aligning backwards");
#endif

        // Then use the reverse complements of the forward alignments as seeds
        seeder_generator(forward_seeder, std::move(reverse_seeds),
                         [&](const auto &reverse_seeder) {
            align_core(reverse, reverse_seeder, std::move(reverse_extender),
                [&](DBGAlignment&& path) {
                    // If the path originated from a backwards alignment (forward
                    // alignment of a reverse complement) and did not skip the first
                    // characters (so it is unable to be reversed), change back
                    // to the forward orientation
                    if (path.get_orientation()) {
                        auto forward_path = path;
                        forward_path.reverse_complement(graph_, forward);
                        if (!forward_path.empty()) {
                            path = std::move(forward_path);
#ifndef NDEBUG
                        } else {
                            mtg::common::logger->trace(
                                "Backwards alignment cannot be reversed, returning");
#endif
                        }

                        // for (node_index node : forward_path) {
                        //     visited_nodes_.emplace(node);
                        // }
                    }

                    assert(path.is_valid(graph_, &config_));
                    alignment_callback(std::move(path));
                },
                get_min_path_score
            );
        });
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
