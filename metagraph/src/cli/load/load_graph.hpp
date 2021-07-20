#ifndef __LOAD_GRAPH_HPP__
#define __LOAD_GRAPH_HPP__

#include <string>

#include "common/logger.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

Config::GraphType parse_graph_type(const std::string &filename);

template <class Graph>
std::shared_ptr<Graph> load_critical_graph_from_file(const std::string &filename) {
    auto graph = std::make_shared<Graph>(2);
    if (!graph->load(filename)) {
        mtg::common::logger->error("Cannot load graph from file '{}'", filename);
        exit(1);
    }
    return graph;
}

std::shared_ptr<graph::DeBruijnGraph> load_critical_dbg(const std::string &filename);

std::shared_ptr<graph::DeBruijnGraph>
primary_to_canonical(std::shared_ptr<graph::DeBruijnGraph> graph,
                     size_t cache_size = 1024);

std::shared_ptr<graph::DeBruijnGraph>
make_cached_graph(std::shared_ptr<graph::DeBruijnGraph> graph,
                  const Config &config,
                  size_t cache_size = 1024);

} // namespace cli
} // namespace mtg

#endif // __LOAD_GRAPH_HPP__
