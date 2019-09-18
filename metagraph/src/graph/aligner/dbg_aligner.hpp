#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <functional>

#include "sequence_graph.hpp"
#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
#include "bounded_priority_queue.hpp"
#include "reverse_complement.hpp"


template <class AlignmentCompare = std::less<Alignment<DeBruijnGraph::node_index>>,
          class ColumnCompare = std::less<typename Alignment<DeBruijnGraph::node_index>::Column>>
class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    explicit DBGAligner(const DeBruijnGraph &graph,
                        const DBGAlignerConfig &config,
                        const Seeder<node_index> &seed = suffix_seeder<node_index>,
                        const Extender<node_index> &extend = default_extender<node_index, ColumnCompare>)
          : graph_(graph),
            config_(config),
            seed_(seed),
            extend_(extend) {}

    explicit DBGAligner(const DeBruijnGraph &graph,
                        const Config &config,
                        const Seeder<node_index> &seed = suffix_seeder<node_index>,
                        const Extender<node_index> &extend = default_extender<node_index, ColumnCompare>)
          : DBGAligner(graph,
                       DBGAlignerConfig(config, graph),
                       seed,
                       extend) {}

    DBGQueryAlignment align(const std::string &query,
                            score_t min_path_score = std::numeric_limits<score_t>::min()) const {
        return config_.forward_and_reverse_complement
            ? align_forward_and_reverse_complement(query, min_path_score)
            : align_one_direction(query, false, min_path_score);
    }

    DBGQueryAlignment
    extend_mapping_forward_and_reverse_complement(const std::string &query,
                                                  score_t min_path_score
                                                      = std::numeric_limits<score_t>::min(),
                                                  const MapExtendSeederBuilder<node_index> &seeder_builder
                                                      = build_unimem_seeder<node_index>) const {
        std::vector<node_index> nodes;
        nodes.reserve(query.size() - graph_.get_k() + 1);
        graph_.map_to_nodes_sequentially(query.begin(),
                                         query.end(),
                                         [&](auto node) { nodes.emplace_back(node); });
        assert(nodes.size() == query.size() - graph_.get_k() + 1);

        auto seeder = seeder_builder(nodes, graph_);

        auto paths = DBGAligner(graph_, config_, seeder, extend_).align_one_direction(
            query, false, min_path_score
        );

        auto size = paths.size();

        std::vector<node_index> rc_nodes;
        rc_nodes.reserve(nodes.size());
        graph_.map_to_nodes_sequentially(paths.get_query_reverse_complement().begin(),
                                         paths.get_query_reverse_complement().end(),
                                         [&](auto node) { rc_nodes.emplace_back(node); });
        assert(rc_nodes.size() == nodes.size());

        auto rc_seeder = seeder_builder(rc_nodes, graph_);

        DBGAligner(graph_, config_, rc_seeder, extend_).align(
            paths.get_query_reverse_complement().begin(),
            paths.get_query_reverse_complement().end(),
            [&](DBGAlignment&& alignment) { paths.emplace_back(std::move(alignment)); },
            true,
            size >= config_.num_alternative_paths
                ? paths[config_.num_alternative_paths - 1].get_score() - 1
                : min_path_score
        );

        std::inplace_merge(paths.begin(),
                           paths.begin() + size,
                           paths.end(),
                           std::greater<DBGAlignment>());

        paths.erase(paths.begin() + std::min(paths.size(), config_.num_alternative_paths),
                    paths.end());

        return paths;
    }

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }
    const Seeder<node_index>& get_seeder() const { return seed_; }
    const Extender<node_index>& get_extender() const { return extend_; }

  private:
    // A convenience function
    DBGQueryAlignment
    align_one_direction(const std::string &query,
                        bool orientation = false,
                        score_t min_path_score = std::numeric_limits<score_t>::min()) const {
        DBGQueryAlignment paths(query);

        if (orientation)
            std::swap(const_cast<std::string&>(paths.get_query()),
                      const_cast<std::string&>(paths.get_query_reverse_complement()));

        const auto& query_alignment = orientation ? paths.get_query_reverse_complement()
                                                  : paths.get_query();

        assert(query_alignment == query);

        align(query_alignment.begin(),
              query_alignment.end(),
              [&](auto&& path) { paths.emplace_back(std::move(path)); },
              orientation,
              min_path_score);

        return paths;
    }

    DBGQueryAlignment
    align_forward_and_reverse_complement(const std::string &query,
                                         score_t min_path_score
                                             = std::numeric_limits<score_t>::min()) const {
        auto paths = align_one_direction(query, false, min_path_score);
        auto size = paths.size();

        align(paths.get_query_reverse_complement().begin(),
              paths.get_query_reverse_complement().end(),
              [&](DBGAlignment&& alignment) { paths.emplace_back(std::move(alignment)); },
              true,
              size >= config_.num_alternative_paths
                  ? paths[config_.num_alternative_paths - 1].get_score() - 1
                  : min_path_score);

        std::inplace_merge(paths.begin(),
                           paths.begin() + size,
                           paths.end(),
                           std::greater<DBGAlignment>());

        paths.erase(paths.begin() + std::min(paths.size(), config_.num_alternative_paths),
                    paths.end());

        return paths;
    }

    // Align a sequence to the graph
    template <class StringIt>
    void align(StringIt query_begin_it,
               StringIt query_end_it,
               const std::function<void(DBGAlignment&&)> &callback,
               bool orientation = false,
               score_t min_path_score = std::numeric_limits<score_t>::min()) const {
        const char *query_begin = &*query_begin_it;
        const char *query_end = &*query_end_it;

        min_path_score = std::max(min_path_score, config_.min_cell_score);

        BoundedPriorityQueue<DBGAlignment> path_queue(config_.num_alternative_paths);

        BoundedPriorityQueue<DBGAlignment,
                             std::vector<DBGAlignment>,
                             AlignmentCompare> partial_paths(config_.queue_size);

        std::vector<DBGAlignment> next_paths;
        for (auto it = query_begin; it + graph_.get_k() <= query_end; ++it) {
            bool full_seed = false;

            auto seeds = seed_(graph_,
                               config_,
                               it,
                               query_end,
                               it - query_begin,
                               orientation);
            assert(seeds.size() <= config_.max_num_seeds_per_locus);

            for (auto&& seed : seeds) {
                assert(seed.size());
                assert(seed.get_score() > config_.min_cell_score);
                assert(seed.is_valid(graph_));

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

                // continue extending until the path is depleted
                do {
                    next_paths.clear();
                    extend_(graph_,
                            partial_path,
                            &next_paths,
                            query_end,
                            config_,
                            orientation,
                            min_path_score);

                    assert(std::all_of(next_paths.begin(),
                                       next_paths.end(),
                                       [&](const auto &c) {
                                           return partial_path.get_score() + c.get_score()
                                               >= min_path_score;
                                       }));

                    assert(std::all_of(next_paths.begin(),
                                       next_paths.end(),
                                       [&](const auto &c) {
                                           return partial_path.get_query_end()
                                                   == c.get_query_begin();
                                       }));

                    assert(std::all_of(next_paths.begin(),
                                       next_paths.end(),
                                       [&](const auto &c) {
                                           return (c.get_nodes().size() || c.get_query_end()
                                                       > partial_path.get_query_end());
                                       }));

                    // only one extension. take it if good enough
                    if (next_paths.size() == 1)
                        partial_path.append(std::move(next_paths.front()));

                } while (partial_path.get_query_end() < query_end && next_paths.size() == 1);

                if (next_paths.size() <= 1) {
                    assert(next_paths.empty() || partial_path.get_query_end() == query_end);
                    if (partial_path.get_score() >= min_path_score) {
                        if (full_seed)
                            it = partial_path.get_query_end() - 1;

                        path_queue.emplace(std::move(partial_path));

                        if (path_queue.size() == config_.num_alternative_paths)
                            min_path_score = std::max(min_path_score,
                                                      path_queue.bottom().get_score());
                    }

                    next_paths.clear();
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
                            it = path_extend.get_query_end() - 1;

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

                next_paths.clear();
            }
        }

        while (path_queue.size()) {
            //assert(path_queue.top().is_valid(graph_));
            callback(path_queue.pop_top());
        }
    }

    const DeBruijnGraph& graph_;
    DBGAlignerConfig config_;
    const Seeder<node_index> seed_;
    const Extender<node_index> extend_;
};


#endif // __DBG_ALIGNER_HPP__
