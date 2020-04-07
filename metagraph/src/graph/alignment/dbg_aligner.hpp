#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "common/bounded_priority_queue.hpp"
#include "graph/representation/base/sequence_graph.hpp"


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
        return config_.forward_and_reverse_complement
            ? align_forward_and_reverse_complement(query)
            : align_one_direction(query, false);
    }

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  protected:
    typedef BoundedPriorityQueue<DBGAlignment,
                                 std::vector<DBGAlignment>,
                                 AlignmentCompare> AlignmentQueue;

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
        AlignmentQueue path_queue(config_.num_alternative_paths);

        alignment_generator(
            [&](DBGAlignment&& alignment) { path_queue.emplace(std::move(alignment)); },
            [&](const DBGAlignment &) {
                return path_queue.size() ? path_queue.bottom().get_score()
                                         : config_.min_path_score;
            }
        );

        while (path_queue.size()) {
            auto path = path_queue.pop_top();
            assert(path.is_valid(graph_, &config_));
            callback(std::move(path));
        }
    }

    const DeBruijnGraph& graph_;
    DBGAlignerConfig config_;
};


#endif // __DBG_ALIGNER_HPP__
