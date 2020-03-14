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

    virtual DBGQueryAlignment align(const std::string &query) const = 0;

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

    DBGQueryAlignment align(const std::string &query) const {
        Seeder seeder(graph_, config_);
        return config_.forward_and_reverse_complement
            ? align_forward_and_reverse_complement(query, seeder)
            : align_one_direction(query, false, seeder);
    }

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  private:
    typedef BoundedPriorityQueue<DBGAlignment,
                                 std::vector<DBGAlignment>,
                                 AlignmentCompare> AlignmentQueue;

    DBGQueryAlignment align_one_direction(const std::string &query,
                                          bool orientation,
                                          Seeder &seeder) const {
        DBGQueryAlignment paths(query);

        if (orientation) {
            std::swap(const_cast<std::string&>(paths.get_query()),
                      const_cast<std::string&>(paths.get_query_reverse_complement()));
        }

        const auto& query_alignment = orientation ? paths.get_query_reverse_complement()
                                                  : paths.get_query();
        assert(query_alignment == query);

        seeder.initialize(query_alignment, orientation);

        align(query_alignment,
              [&](DBGAlignment&& path) {
                  paths.emplace_back(std::move(path));
              },
              orientation,
              config_.min_path_score,
              seeder);

        return paths;
    }

    DBGQueryAlignment align_forward_and_reverse_complement(const std::string &query,
                                                           Seeder &seeder) const {
        DBGQueryAlignment final_paths(query);
        const auto &forward = final_paths.get_query();
        const auto &rev_comp = final_paths.get_query_reverse_complement();

        std::vector<DBGAlignment> all_paths;

        seeder.initialize(forward, false);

        align(forward,
              [&](DBGAlignment&& alignment) {
                  all_paths.emplace_back(std::move(alignment));
              },
              false,
              config_.min_path_score,
              seeder);

        auto size = all_paths.size();

        seeder.initialize(rev_comp, true);

        align(rev_comp,
              [&](DBGAlignment&& alignment) {
                  all_paths.emplace_back(std::move(alignment));
              },
              true,
              size >= config_.num_alternative_paths
                  ? all_paths[config_.num_alternative_paths - 1].get_score() - 1
                  : config_.min_path_score,
              seeder);

        size_t max_left = std::min(size, config_.num_alternative_paths);
        size_t max_right = std::min(all_paths.size() - size, config_.num_alternative_paths);
        std::merge(std::make_move_iterator(all_paths.begin()),
                   std::make_move_iterator(all_paths.begin() + max_left),
                   std::make_move_iterator(all_paths.begin() + size),
                   std::make_move_iterator(all_paths.begin() + size + max_right),
                   std::back_inserter(final_paths),
                   std::not_fn(AlignmentCompare()));

        final_paths.erase(final_paths.begin()
                              + std::min(final_paths.size(), config_.num_alternative_paths),
                          final_paths.end());

        return final_paths;
    }

    // Align a sequence to the graph
    void align(std::string_view query,
               const std::function<void(DBGAlignment&&)> &callback,
               bool orientation,
               score_t min_path_score,
               const Seeder &seeder) const {
        assert(config_.check_config_scores());
        min_path_score = std::max(min_path_score, config_.min_cell_score);

        AlignmentQueue path_queue(config_.num_alternative_paths);

        Extender extend(graph_, config_);
        extend.initialize_query(query);

        seeder.call_seeds([&](DBGAlignment&& seed) {
            assert(seed.get_query().data() >= query.data());
            assert(seed.get_query_end() <= query.data() + query.size());
            assert(seed.get_query_end() > seed.get_query().data());
            assert(seed.get_clipping() == seed.get_query().data() - query.data());
            assert(seed.is_valid(graph_, &config_));

            if (seed.get_query_end() == query.data() + query.size()
                    && seed.get_score() >= config_.min_path_score) {
                path_queue.emplace(std::move(seed));
                return;
            }

            extend.initialize(seed);
            bool extended = false;
            extend(seed, query, [&](DBGAlignment&& extension, auto start_node) {
                assert(extension.is_valid(graph_, &config_));
                extension.extend_query_end(query.data() + query.size());

                if (extension.get_clipping() || start_node != seed.back()) {
                    // if the extension starts at a different position
                    // from the seed end, then it's a new alignment
                    extension.extend_query_begin(query.data());
                    path_queue.emplace(std::move(extension));
                    return;
                }

                extended = true;

                auto next_path = seed;
                next_path.append(std::move(extension));
                assert(next_path.is_valid(graph_, &config_));

                path_queue.emplace(std::move(next_path));
            }, orientation, min_path_score);

            if (!extended && seed.get_score() >= config_.min_path_score) {
                seed.extend_query_end(query.data() + query.size());
                assert(seed.is_valid(graph_, &config_));
                path_queue.emplace(std::move(seed));
            }
        });

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
