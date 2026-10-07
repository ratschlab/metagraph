#ifndef __LOAD_GRAPH_HPP__
#define __LOAD_GRAPH_HPP__

#include <string>

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {

namespace graph {
class DBGSuccinct;
} // namespace graph

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

// What mask_dummy_edges masked (the figures of `metagraph stats --count-dummy`) and its time
struct DummyMaskCounts {
    uint64_t edges = 0;         // BOSS edges, the graph's max_index
    uint64_t kmers = 0;         // edges left valid: the graph's k-mers
    uint64_t sink_dummy = 0;    // masked edges with W = $ (the main dummy edge 1 excluded)
    uint64_t source_dummy = 0;  // the other masked edges: the dummy tree, its root 1 included
    double seconds = 0;
};

/**
 * The dummy-edge mask of a succinct graph (DESIGN-pattern-search.md §4), built in |graph|
 * exactly as `metagraph build --mask-dummy` builds it: DBGSuccinct::mask_dummy_kmers without
 * pruning, which marks the source dummies by traversing their tree and the sink dummies
 * (W = $) and removes no edge, so node ids, and the annotation rows that follow them, stay as
 * they are. Replaces a mask |graph| had. Used by `transform --mask-dummy`, which writes the
 * mask beside the graph, and by --pattern-build-mask, which builds it at load; one function,
 * so that the two cannot drift from each other or from the build's.
 */
DummyMaskCounts mask_dummy_edges(graph::DBGSuccinct *graph, size_t num_threads);

} // namespace cli
} // namespace mtg

#endif // __LOAD_GRAPH_HPP__
