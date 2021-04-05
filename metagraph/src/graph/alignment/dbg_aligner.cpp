#include "dbg_aligner.hpp"

#include "common/algorithms.hpp"

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
             const MinScoreComputer &get_min_path_score) {
    bool filter_seeds = dynamic_cast<const ExactSeeder<node_index>*>(&seeder);

    std::vector<DBGAlignment> seeds = seeder.get_seeds();

    for (const DBGAlignment &seed : seeds) {
        if (seed.get_score() >= get_min_path_score(seed)) {
            DBGAlignment out_seed(seed);
            out_seed.extend_query_end(query.data() + query.size());
            out_seed.trim_offset();
            assert(out_seed.is_valid(graph_, &config_));
            DEBUG_LOG("Alignment (seed): {}", out_seed);
            callback(std::move(out_seed));
        }
    }

    if (filter_seeds
            && std::any_of(seeds.begin(), seeds.end(),
                           [&](const auto &a) { return a.target_columns.size(); })) {
        constexpr auto get_end_ptr = [](const DBGAlignment &a) {
            return a.get_query().data() + a.get_query().size();
        };

        std::sort(seeds.begin(), seeds.end(), [&](const auto &a, const auto &b) {
            return get_end_ptr(a) < get_end_ptr(b);
        });

        for (size_t i = 0; i + 1 < seeds.size(); ++i) {
            DBGAlignment &seed = seeds[i];
            if (get_end_ptr(seeds[i + 1]) - get_end_ptr(seed) == 1
                    && seed.target_columns.size()
                    && seeds[i + 1].target_columns.size()) {
                uint64_t cnt = utils::count_intersection(
                    seed.target_columns.begin(), seed.target_columns.end(),
                    seeds[i + 1].target_columns.begin(), seeds[i + 1].target_columns.end()
                );

                if (!cnt || (cnt == seed.target_columns.size()
                                && graph_.get_mode() == DeBruijnGraph::CANONICAL)) {
                    DEBUG_LOG("Skipping seed: {}", seed);
                    seed = DBGAlignment();
                }
            }
        }
    }

    std::sort(seeds.begin(), seeds.end(), LocalAlignmentGreater());

    std::vector<sdsl::bit_vector> seed_nodes_covered;
    std::vector<Vector<uint64_t>> seed_target_updates(seeds.size());
    for (size_t i = 0; i < seeds.size(); ++i) {
        seed_nodes_covered.emplace_back(seeds[i].size(), false);
    }

    for (size_t i = 0; i < seeds.size(); ++i) {
        DBGAlignment &seed = seeds[i];
        if (seed.empty())
            continue;

        score_t min_path_score = get_min_path_score(seed);
        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

        extender.initialize(seed);

        auto extensions = extender.get_extensions(min_path_score);

        // if the ManualSeeder is not used, then add nodes to the visited_nodes_
        // table to allow for seed filtration
        if (filter_seeds) {
            extender.call_visited_nodes([&](node_index node, size_t begin, size_t end) {
                for (size_t j = i + 1; j < seeds.size(); ++j) {
                    if (seeds[j].empty())
                        continue;

                    assert(!seeds[j].get_offset() || seeds[j].size() == 1);
                    size_t s_begin = seeds[j].get_clipping();
                    size_t s_end = s_begin + graph_.get_k() - seeds[j].get_offset();

                    assert(sdsl::util::cnt_one_bits(seed_nodes_covered[j]) < seeds[j].size());
                    bool found = true;
                    for (size_t l = 0; l < seeds[j].size(); ++l) {
                        assert(s_end <= seeds[j].get_clipping() + seeds[j].get_query().size());
                        if (seeds[j][l] == node && begin <= s_begin && end >= s_end) {
                            found = true;
                            seed_nodes_covered[j][l] = true;
                        }

                        assert(s_end - s_begin <= graph_.get_k());
                        if (s_end - s_begin == graph_.get_k())
                            ++s_begin;

                        ++s_end;
                    }

                    Vector<uint64_t> diff;
                    if (found) {
                        std::set_difference(seeds[j].target_columns.begin(),
                                            seeds[j].target_columns.end(),
                                            seed.target_columns.begin(),
                                            seed.target_columns.end(),
                                            std::back_inserter(diff));
                        Vector<uint64_t> target_union;
                        std::set_union(seed_target_updates[j].begin(), seed_target_updates[j].end(),
                                       diff.begin(), diff.end(), std::back_inserter(target_union));
                        std::swap(seed_target_updates[j], target_union);
                    }

                    if (sdsl::util::cnt_one_bits(seed_nodes_covered[j]) == seeds[j].size()) {
                        assert(found);
                        if (diff.empty()) {
                            seeds[j] = DBGAlignment();
                            seed_nodes_covered[j] = sdsl::bit_vector{};
                        } else {
                            std::swap(seed_target_updates[j], seeds[j].target_columns);
                            sdsl::util::set_to_value(seed_nodes_covered[j], false);
                        }
                        seed_target_updates[j] = Vector<uint64_t>{};
                    }
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
::align_one_direction(bool orientation_to_align,
                      const ISeeder<node_index> &seeder,
                      IExtender<node_index> &extender) {
    std::string_view query = paths_.get_query(orientation_to_align);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        align_core(query, seeder, extender, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_best_direction(const ISeeder<node_index> &seeder,
                       const ISeeder<node_index> &seeder_rc,
                       IExtender<node_index> &extender,
                       IExtender<node_index> &extender_rc) {
    std::string_view forward = paths_.get_query();
    std::string_view reverse = paths_.get_query(true);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
        align_core(forward, seeder, extender, alignment_callback, get_min_path_score);
        align_core(reverse, seeder_rc, extender_rc, alignment_callback, get_min_path_score);
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_both_directions(const ISeeder<node_index> &forward_seeder,
                        const ISeeder<node_index> &reverse_seeder,
                        IExtender<node_index> &forward_extender,
                        IExtender<node_index> &reverse_extender,
                        const SeederGenerator &rev_comp_core_generator) {
    std::string_view forward = paths_.get_query();
    std::string_view reverse = paths_.get_query(true);

    align_aggregate([&](const auto &alignment_callback, const auto &get_min_path_score) {
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

                if (!path.get_clipping() || path.get_offset())
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

        DEBUG_LOG("Extending in reverse direction");
        rev_comp_core_generator(std::move(rc_of_reverse), [&](const auto &seeder) {
            align_core(forward, seeder, forward_extender,
                       alignment_callback, get_min_path_score);
        });

        rev_comp_core_generator(std::move(rc_of_forward), [&](const auto &seeder) {
            align_core(reverse, seeder, reverse_extender,
                       alignment_callback, get_min_path_score);
        });
    });
}

template <class AlignmentCompare>
void SeedAndExtendAlignerCore<AlignmentCompare>
::align_aggregate(const AlignmentGenerator &alignment_generator) {
    alignment_generator(
        [&](DBGAlignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            aggregator_.add_alignment(std::move(alignment));
        },
        [&](const DBGAlignment &seed) { return aggregator_.get_min_path_score(seed); }
    );
}

template class SeedAndExtendAlignerCore<>;

} // namespace align
} // namespace graph
} // namespace mtg
