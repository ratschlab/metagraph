#ifndef __LOAD_GRAPH_HPP__
#define __LOAD_GRAPH_HPP__

#include <stdexcept>
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

/**
 * What a loader does with a file of the index it cannot load:
 *   EXIT   logs the failure and exits with status 1: a command's or a server's start-up load,
 *          with nothing to serve without the index (and the server's start-up loads run inside
 *          OpenMP regions, which an exception must not leave)
 *   THROW  throws it as a std::runtime_error whose message names the file and the step: a
 *          server's load of an index for one request (`in_ram`), which that request is
 *          answered with while the server runs on with its resident indexes
 */
enum class LoadFailure { EXIT, THROW };

/** Loads |filename| with the graph type's extension (|.dbg|, |.orhashdbg|, etc.). On failure,
 *  logs the resolved path and likely causes (permissions, etc.) and exits, or throws that
 *  message (|on_failure|); implementations may also log. */
template <class Graph>
std::shared_ptr<Graph> load_critical_graph_from_file(const std::string &filename,
                                                     LoadFailure on_failure = LoadFailure::EXIT) {
    auto graph = std::make_shared<Graph>(2);
    if (!graph->load(filename)) {
        const std::string on_disk = utils::make_suffix(filename, Graph::kExtension);
        const std::string message = fmt::format("Cannot load graph from '{}': {}", on_disk,
                                                utils::file_read_failure_detail(on_disk));
        if (on_failure == LoadFailure::THROW) {
            throw std::runtime_error(message);
        }
        common::logger->error("{}", message);
        exit(1);
    }
    return graph;
}

std::shared_ptr<graph::DeBruijnGraph>
load_critical_dbg(const std::string &filename, LoadFailure on_failure = LoadFailure::EXIT);

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

/**
 * After |graph| was written to |filename| (DBGSuccinct::serialize writes a .edgemask only for a
 * graph that has a mask, and leaves an older one in place): when |graph| is a succinct graph
 * without a mask, removes a .edgemask beside it, which another graph left there and which would
 * be loaded with this one (a stale mask of the right size passes the load's check), warning
 * that it did, or that it could not (the user is told to delete it). Called by the commands
 * that write a new graph over a prefix: extend, build --in-ram and merge; not by transform,
 * which rewrites the same graph (the mask beside it is the graph's).
 */
void remove_stale_dummy_mask(const graph::DeBruijnGraph &graph, const std::string &filename);

} // namespace cli
} // namespace mtg

#endif // __LOAD_GRAPH_HPP__
