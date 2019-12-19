#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "common/bounded_priority_queue.hpp"
#include "graph/base/sequence_graph.hpp"


class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~IDBGAligner() {}

    virtual DBGQueryAlignment
    align(const std::string &query) const = 0;

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

        seeder.initialize(query, orientation);

        align(query_alignment.c_str(), query_alignment.c_str() + query_alignment.size(),
              [&](DBGAlignment&& path) {
                  paths.emplace_back(std::move(path));
              },
              orientation,
              config_.min_path_score,
              seeder
        );

        return paths;
    }

    DBGQueryAlignment align_forward_and_reverse_complement(const std::string &query,
                                                           Seeder &seeder) const {
        DBGQueryAlignment final_paths(query);
        const auto &forward = final_paths.get_query();
        const auto &rev_comp = final_paths.get_query_reverse_complement();

        std::vector<DBGAlignment> all_paths;

        seeder.initialize(forward, false);

        align(forward.c_str(), forward.c_str() + forward.size(),
              [&](DBGAlignment&& alignment) {
                  all_paths.emplace_back(std::move(alignment));
              },
              false,
              config_.min_path_score,
              seeder);

        auto size = all_paths.size();

        seeder.initialize(rev_comp, true);

        align(rev_comp.c_str(), rev_comp.c_str() + rev_comp.size(),
              [&](DBGAlignment&& alignment) {
                  all_paths.emplace_back(std::move(alignment));
              },
              true,
              size >= config_.num_alternative_paths
                  ? all_paths[config_.num_alternative_paths - 1].get_score() - 1
                  : config_.min_path_score,
              seeder);

        std::merge(std::make_move_iterator(all_paths.begin()),
                   std::make_move_iterator(all_paths.begin() + size),
                   std::make_move_iterator(all_paths.begin() + size),
                   std::make_move_iterator(all_paths.end()),
                   std::back_inserter(final_paths),
                   std::not_fn(AlignmentCompare()));

        final_paths.erase(final_paths.begin()
                              + std::min(final_paths.size(), config_.num_alternative_paths),
                          final_paths.end());

        return final_paths;
    }

    // Align a sequence to the graph
    void align(const char *query_begin,
               const char *query_end,
               const std::function<void(DBGAlignment&&)> &callback,
               bool orientation,
               score_t min_path_score,
               const Seeder &seeder) const {
        assert(config_.check_config_scores());
        min_path_score = std::max(min_path_score, config_.min_cell_score);

        BoundedPriorityQueue<DBGAlignment> path_queue(config_.num_alternative_paths);

        // a different comparator may be used when picking the next seed to extend
        BoundedPriorityQueue<DBGAlignment,
                             std::vector<DBGAlignment>,
                             AlignmentCompare> partial_paths(config_.queue_size);

        Extender extend(graph_, config_);

        // compute perfect match scores for all suffixes
        // used for branch and bound checks below
        std::vector<score_t> partial_sum(query_end - query_begin);
        std::transform(query_begin, query_end,
                       partial_sum.begin(),
                       [&](char c) { return config_.get_row(c)[c]; });

        std::partial_sum(partial_sum.rbegin(), partial_sum.rend(), partial_sum.rbegin());
        assert(config_.match_score(query_begin, query_end) == partial_sum.front());
        assert(config_.get_row(*(query_end - 1))[*(query_end - 1)] == partial_sum.back());

        std::vector<DBGAlignment> next_paths;
        for (auto it = query_begin; it + graph_.get_k() <= query_end; ++it) {
            if (partial_sum[it - query_begin] < min_path_score)
                break;

            bool full_seed = false;

            auto seeds = seeder(it, query_end, it - query_begin, orientation);
            assert(seeds.size() <= config_.max_num_seeds_per_locus);

            for (auto&& seed : seeds) {
                assert(seed.size());
                assert(seed.get_score() > config_.min_cell_score);
                assert(seed.is_valid(graph_, &config_));

                if (seed.get_num_matches() >= graph_.get_k())
                    full_seed = true;

                if (seed.get_query_end() == query_end) {
                    path_queue.emplace(std::move(seed));
                } else {
                    partial_paths.emplace(std::move(seed));
                }
            }

            if (partial_paths.empty() && full_seed)
                break;

            // a seed has been found
            while (partial_paths.size()) {
                next_paths.clear();
                auto partial_path = partial_paths.top();

                assert(partial_path.get_score() > config_.min_cell_score);
                assert(partial_path.get_query_end() >= partial_path.get_query_begin());
                assert(partial_path.size());

                if (partial_path.get_query_end() == query_end) {
                    // if at the end, go for the next path
                    if (partial_path.get_score() >= min_path_score)
                        path_queue.emplace(std::move(partial_path));

                    partial_paths.pop();

                    continue;
                }

                extend.initialize(partial_path);

                // continue extending until the path is depleted
                while (partial_path.get_query_end() < query_end) {
                    auto cur_paths = extend(
                        partial_path,
                        query_end,
                        partial_sum.data()
                            + (partial_path.get_query_end() - 1 - query_begin),
                        orientation,
                        min_path_score
                    );

                    // The graph stored in extend may differ from the one stored
                    // in the DBGAligner and Seeder
                    assert(std::all_of(
                        cur_paths.begin(), cur_paths.end(),
                        [&](const auto &c) {
                            return c.is_valid(extend.get_graph(), &extend.get_config());
                        }
                    ));

                    assert(std::all_of(
                        cur_paths.begin(), cur_paths.end(),
                        [&](const auto &c) {
                            return partial_path.get_score() + c.get_score()
                                >= min_path_score;
                        }
                    ));

                    assert(std::all_of(
                        cur_paths.begin(), cur_paths.end(),
                        [&](const auto &c) {
                            return partial_path.get_query_end() + c.get_clipping()
                                    == c.get_query_begin();
                        }
                    ));

                    // only one extension. take it if good enough
                    if (cur_paths.size() == 1) {
                        partial_path.append(std::move(cur_paths.front()));
                    } else {
                        next_paths.assign(std::make_move_iterator(cur_paths.begin()),
                                          std::make_move_iterator(cur_paths.end()));
                        break;
                    }
                }

                assert(next_paths.size() != 1);

                if (next_paths.empty() && partial_path.get_score() >= min_path_score) {
                    if (full_seed)
                        ++it;

                    partial_path.extend_query_end(query_end);
                    path_queue.emplace(std::move(partial_path));

                    if (path_queue.size() == config_.num_alternative_paths) {
                        min_path_score = std::max(min_path_score,
                                                  path_queue.bottom().get_score());
                    }
                }

                bool picked = false;

                for (auto&& path : next_paths) {
                    assert(partial_path.get_score() + path.get_score()
                                > config_.min_cell_score);
                    assert(partial_path.get_query_end() == path.get_query_begin());

                    auto path_extend = partial_path;
                    path_extend.append(std::move(path));

                    if (path_extend.get_query_end() == query_end) {
                        if (full_seed)
                            ++it;

                        path_queue.emplace(std::move(path_extend));

                        if (path_queue.size() == config_.num_alternative_paths)
                            min_path_score = std::max(min_path_score,
                                                      path_queue.bottom().get_score());
                    } else {
                        if (!picked) {
                            partial_paths.update(partial_paths.top_it(),
                                                 std::move(path_extend));
                        } else {
                            partial_paths.emplace(std::move(path_extend));
                        }

                        picked = true;
                    }
                }

                if (!picked) {
                    assert(partial_paths.size());
                    partial_paths.pop();
                }
            }
        }

        while (path_queue.size()) {
            callback(path_queue.pop_top());
        }
    }

    const DeBruijnGraph& graph_;
    DBGAlignerConfig config_;
};


#endif // __DBG_ALIGNER_HPP__
