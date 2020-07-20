#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include <priority_deque.hpp>

#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {

class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~IDBGAligner() {}

    virtual DBGQueryAlignment align(const std::string_view query) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};


template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>>
class DBGAligner : public IDBGAligner {
  public:
    DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual ~DBGAligner() {}

    DBGQueryAlignment align(const std::string_view query) const {
        if (graph_.is_canonical_mode()) {
            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            return align_both_directions(query);
        } else {
            return config_.forward_and_reverse_complement
                ? align_forward_and_reverse_complement(query)
                : align_one_direction(query, false);
        }
    }

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  protected:
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&)> SeedGenerator;
    typedef const std::function<void(const std::function<void(DBGAlignment&&)>&,
                                     const std::function<score_t(const DBGAlignment&)>&)> AlignmentGenerator;

    // Generate seeds, then extend them
    void align(const std::string_view query,
               const SeedGenerator &seed_generator,
               const std::function<void(DBGAlignment&&)> &callback,
               const std::function<score_t(const DBGAlignment&)> &get_min_path_score) const {
        assert(config_.check_config_scores());

        auto extend = build_extender();
        extend.initialize_query(query);

        seed_generator([&](DBGAlignment&& seed) {
            assert(seed.is_valid(graph_, &config_));

            score_t min_path_score = get_min_path_score(seed);

            if (seed.get_score() >= min_path_score) {
                DBGAlignment path(seed);
                path.extend_query_end(query.data() + query.size());
                callback(std::move(path));
            }

            if (seed.get_query_end() == query.data() + query.size())
                return;

            extend.initialize(seed);
            extend([&](DBGAlignment&& extension, auto start_node) {
                assert(extension.is_valid(graph_, &config_));
                extension.extend_query_end(query.data() + query.size());

                if (extension.get_clipping() || start_node != seed.back()) {
                    // if the extension starts at a different position
                    // from the seed end, then it's a new alignment
                    extension.extend_query_begin(query.data());
                    callback(std::move(extension));
                    return;
                }

                auto next_path = seed;
                next_path.append(std::move(extension));
                assert(next_path.is_valid(graph_, &config_));

                callback(std::move(next_path));
            }, min_path_score);
        });
    }

  private:
    virtual Seeder build_seeder() const { return Seeder(graph_, config_); }
    virtual Extender build_extender() const { return Extender(graph_, config_); }

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    DBGQueryAlignment align_one_direction(const std::string_view query,
                                          bool orientation) const {
        auto seeder = build_seeder();
        DBGQueryAlignment paths(query);

        if (orientation) {
            std::swap(const_cast<std::string&>(paths.get_query()),
                      const_cast<std::string&>(paths.get_query_reverse_complement()));
        }

        const auto &query_alignment = orientation ? paths.get_query_reverse_complement()
                                                  : paths.get_query();
        assert(query_alignment == query);

        seeder.initialize(query_alignment, orientation);

        align_aggregate([&](const auto &alignment_callback,
                            const auto &get_min_path_score) {
            align(query_alignment, [&](const auto &callback) {
                seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
            }, alignment_callback, get_min_path_score);

        }, [&](DBGAlignment&& path) { paths.emplace_back(std::move(path)); });

        return paths;
    }

    // Align both forwards and backwards from a given seed. Procedure
    // 1. Given each seed, extend forward to produce an alignment A
    // 2. Reverse complement the alignment to get A', treated like a new seed
    // 3. Extend A' forwards
    // 4. Reverse complement A' to get the final alignment A''
    DBGQueryAlignment align_both_directions(const std::string_view query) const {
        auto seeder = build_seeder();
        DBGQueryAlignment paths(query);

        seeder.initialize(paths.get_query(), false);

        align_aggregate([&](const auto &alignment_callback,
                            const auto &get_min_path_score) {
            // The outer loop generates uses the reverse complements of the forward alignments
            // as seeds for extension
            align(paths.get_query_reverse_complement(),
                [&](const auto &reverse_seed_callback) {

                    // Inner loop aligns the forward strand of the query
                    align(paths.get_query(),
                        [&](const auto &forward_seed_callback) {
                            seeder.call_seeds([&](DBGAlignment&& seed) {
                                forward_seed_callback(std::move(seed));
                            });
                        },
                        [&](DBGAlignment&& path) {
                            if (!path.get_clipping()) {
                                // If the alignment starts from the beginning of the query,
                                // there's no sequence left for aligning backwards
                                if (path.get_score() >= get_min_path_score(path))
                                    alignment_callback(std::move(path));

                            } else {
                                // Add seed to the list of alignments
                                if (path.get_score() >= get_min_path_score(path))
                                    alignment_callback(DBGAlignment(path));

                                // If the alignment skipped the first characters in the first
                                // node of its path, then don't align backwards
                                if (path.get_offset())
                                    return;

                                path.reverse_complement(
                                    seeder.get_graph(),
                                    paths.get_query_reverse_complement()
                                );

                                // Remove any character skipping from the end so that the
                                // alignment can proceed
                                assert(path.get_end_clipping());
                                path.trim_end_clipping();
                                assert(path.is_valid(graph_, &config_));

                                // Pass the reverse complement of the forward alignment
                                // as a seed for extension
                                reverse_seed_callback(std::move(path));
                            }
                        },
                        [&](const auto&) {
                            // ignore the min path score for the forward alignment,
                            // since it may have a score that is too low before it is
                            // extended backwards
                            return config_.min_cell_score;
                        }
                    );

                }, alignment_callback, get_min_path_score);
            },
            [&](DBGAlignment&& path) {
                // If the path originated from a backwards alignment (forward alignment
                // of a reverse complement) and did not skip the first characters
                // (so it is unable to be reversed), change back to the forward orientation
                if (!path.get_offset() && path.get_orientation())
                    path.reverse_complement(seeder.get_graph(), paths.get_query());

                assert(path.is_valid(graph_, &config_));
                paths.emplace_back(std::move(path));
            }
        );

        return paths;
    }

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    DBGQueryAlignment align_forward_and_reverse_complement(const std::string_view query) const {
        auto seeder = build_seeder();
        DBGQueryAlignment paths(query);
        const auto &forward = paths.get_query();
        const auto &reverse = paths.get_query_reverse_complement();

        align_aggregate([&](const auto &alignment_callback,
                            const auto &get_min_path_score) {
            seeder.initialize(paths.get_query(), false);
            align(forward, [&](const auto &callback) {
                seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
            }, alignment_callback, get_min_path_score);

            seeder.initialize(paths.get_query_reverse_complement(), true);
            align(reverse, [&](const auto &callback) {
                seeder.call_seeds([&](DBGAlignment&& seed) { callback(std::move(seed)); });
            }, alignment_callback, get_min_path_score);

        }, [&](DBGAlignment&& path) { paths.emplace_back(std::move(path)); } );

        return paths;
    }

    // Given alignments generated by a generator, add them to a priority queue
    // and output the top ones.
    virtual void align_aggregate(const AlignmentGenerator &alignment_generator,
                                 const std::function<void(DBGAlignment&&)> &callback) const {
        boost::container::priority_deque<DBGAlignment,
                                         std::vector<DBGAlignment>,
                                         AlignmentCompare> path_queue;

        alignment_generator(
            [&](DBGAlignment&& alignment) {
                if (path_queue.size() < config_.num_alternative_paths) {
                    // add to the queue
                    path_queue.emplace(std::move(alignment));
                } else if (alignment.get_score() > path_queue.minimum().get_score()) {
                    // if the queue is full and the current alignment is better,
                    // replace the bottom element and bubble it up the heap
                    path_queue.update(path_queue.begin(), std::move(alignment));
                }
            },
            [&](const DBGAlignment &) {
                return path_queue.size() ? path_queue.minimum().get_score()
                                         : config_.min_path_score;
            }
        );

        while (path_queue.size()) {
            assert(path_queue.maximum().is_valid(graph_, &config_));
            callback(DBGAlignment(path_queue.maximum()));
            path_queue.pop_maximum();
        }
    }

    const DeBruijnGraph& graph_;
    DBGAlignerConfig config_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
