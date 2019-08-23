#include "dbg_aligner.hpp"

#include "reverse_complement.hpp"
#include "bounded_priority_queue.hpp"


DBGAligner::DBGAligner(const DeBruijnGraph &graph,
                       const DBGAlignerConfig &config,
                       const Seeder &seed,
                       const Extender &extend,
                       const PriorityFunction &priority_function)
      : graph_(graph),
        config_(config),
        seed_(seed),
        extend_(extend),
        priority_function_(priority_function) { }

DBGAligner::DBGAligner(const DeBruijnGraph &graph,
                       const Config &config,
                       const Seeder &seed,
                       const Extender &extend,
                       const PriorityFunction &priority_function)
      : DBGAligner(graph,
                   DBGAlignerConfig(config, graph),
                   seed,
                   extend,
                   priority_function) { }


std::vector<DBGAligner::DBGAlignment>
DBGAligner::align(const std::string &query,
                  bool orientation,
                  score_t min_path_score) const {
    std::vector<DBGAlignment> paths;
    align(query.begin(),
          query.end(),
          [&](auto&& path) { paths.emplace_back(std::move(path)); },
          orientation,
          min_path_score);

    return paths;
}

template <class StringIt>
void DBGAligner::align(StringIt query_begin_it,
                       StringIt query_end_it,
                       const std::function<void(DBGAlignment&&)> &callback,
                       bool orientation,
                       score_t min_path_score) const {
    const char *query_begin = &*query_begin_it;
    const char *query_end = &*query_end_it;

    min_path_score = std::max(min_path_score, config_.min_cell_score);

    BoundedPriorityQueue<DBGAlignment> path_queue(config_.num_alternative_paths);

    // TODO: priority function here
    BoundedPriorityQueue<DBGAlignment> partial_paths(config_.queue_size);

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
        assert(path_queue.top().is_valid(graph_));
        callback(path_queue.pop_top());
    }
}

template void DBGAligner
::align<const char*>(const char*,
                     const char*,
                     const std::function<void(DBGAlignment&&)> &,
                     bool,
                     score_t) const;

template void DBGAligner
::align<char*>(char*,
               char*,
               const std::function<void(DBGAlignment&&)> &,
               bool,
               score_t) const;

template void DBGAligner
::align<std::string::const_iterator>(std::string::const_iterator,
                                     std::string::const_iterator,
                                     const std::function<void(DBGAlignment&&)> &,
                                     bool,
                                     score_t) const;

template void DBGAligner
::align<std::string::iterator>(std::string::iterator,
                               std::string::iterator,
                               const std::function<void(DBGAlignment&&)> &,
                               bool,
                               score_t) const;


std::vector<DBGAligner::DBGAlignment> DBGAligner
::align_forward_and_reverse_complement(const std::string &query,
                                       const std::string &reverse_complement_query,
                                       score_t min_path_score) const {
    assert(query.size() == reverse_complement_query.size());

    auto paths = align(query, false, min_path_score);
    auto size = paths.size();

    align(reverse_complement_query.begin(),
          reverse_complement_query.end(),
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

std::vector<DBGAligner::DBGAlignment> DBGAligner
::extend_mapping_forward_and_reverse_complement(const std::string &query,
                                                const std::string &reverse_complement_query,
                                                score_t min_path_score,
                                                const MapExtendSeederBuilder &seeder_builder) const {
    assert(query.size() == reverse_complement_query.size());

    std::vector<node_index> nodes;
    nodes.reserve(query.size() - graph_.get_k() + 1);
    graph_.map_to_nodes_sequentially(query.begin(),
                                     query.end(),
                                     [&](auto node) { nodes.emplace_back(node); });
    assert(nodes.size() == query.size() - graph_.get_k() + 1);

    auto seeder = seeder_builder(nodes, graph_);

    std::vector<node_index> rc_nodes;
    rc_nodes.reserve(nodes.size());
    graph_.map_to_nodes_sequentially(reverse_complement_query.begin(),
                                     reverse_complement_query.end(),
                                     [&](auto node) { rc_nodes.emplace_back(node); });
    assert(rc_nodes.size() == nodes.size());

    auto rc_seeder = seeder_builder(rc_nodes, graph_);

    auto paths = DBGAligner(graph_, config_, seeder, extend_, priority_function_).align(
        query, false, min_path_score
    );

    auto size = paths.size();

    DBGAligner(graph_, config_, rc_seeder, extend_, priority_function_).align(
        reverse_complement_query.begin(),
        reverse_complement_query.end(),
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
