#ifndef __LOAD_GRAPH_HPP__
#define __LOAD_GRAPH_HPP__

#include <string>

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace cli {

/** Loads |filename| with the graph type's extension (|.dbg|, |.orhashdbg|, etc.). On failure,
 *  logs the resolved path and likely causes (permissions, etc.); implementations may also log. */
template <class Graph>
std::shared_ptr<Graph> load_critical_graph_from_file(const std::string &filename) {
    auto graph = std::make_shared<Graph>(2);
    if (!graph->load(filename)) {
        const std::string on_disk = utils::make_suffix(filename, Graph::kExtension);
        common::logger->error("Cannot load graph from '{}': {}", on_disk,
                              utils::file_read_failure_detail(on_disk));
        exit(1);
    }
    return graph;
}

std::shared_ptr<graph::DeBruijnGraph> load_critical_dbg(const std::string &filename);

} // namespace cli
} // namespace mtg

#endif // __LOAD_GRAPH_HPP__
