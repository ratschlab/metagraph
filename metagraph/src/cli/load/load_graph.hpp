#ifndef __LOAD_GRAPH_HPP__
#define __LOAD_GRAPH_HPP__

#include <string>

#include "common/logger.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "cli/config/config.hpp"


Config::GraphType parse_graph_type(const std::string &filename);

template <class Graph>
std::shared_ptr<Graph> load_critical_graph_from_file(const std::string &filename) {
    auto graph = std::make_shared<Graph>(2);
    if (!graph->load(filename)) {
        mg::common::logger->error("Cannot load graph from file '{}'", filename);
        exit(1);
    }
    return graph;
}

std::shared_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename);

#endif // __LOAD_GRAPH_HPP__
