#include "load_graph.hpp"

#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using mtg::common::logger;

using GraphLoader = std::shared_ptr<DeBruijnGraph> (*)(const std::string&, LoadFailure);

template <class Graph>
std::shared_ptr<DeBruijnGraph> load_as(const std::string &filename, LoadFailure on_failure) {
    return load_critical_graph_from_file<Graph>(filename, on_failure);
}

#define GRAPH_LOADER_ENTRY(GRAPH) { GRAPH::kExtension, &load_as<GRAPH> }
static const std::vector<std::pair<std::string, GraphLoader>> kGraphLoaders = {
    GRAPH_LOADER_ENTRY(DBGSuccinct),
    GRAPH_LOADER_ENTRY(DBGHashOrdered),
    GRAPH_LOADER_ENTRY(DBGHashString),
    GRAPH_LOADER_ENTRY(DBGHashFast),
    GRAPH_LOADER_ENTRY(DBGBitmap),
    GRAPH_LOADER_ENTRY(DBGSSHash),
};

void remove_stale_dummy_mask(const DeBruijnGraph &graph, const std::string &filename) {
    const auto *succinct = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!succinct || succinct->get_mask())
        return;
    const std::string stale = utils::remove_suffix(filename, DBGSuccinct::kExtension)
                                + DBGSuccinct::kDummyMaskExtension;
    std::error_code error;
    if (!std::filesystem::exists(stale, error))
        return;
    if (std::filesystem::remove(stale, error)) {
        logger->warn("Removed the dummy-edge mask {} left beside the output: it was not this "
                     "graph's (the graph written has no mask)", stale);
    } else {
        logger->warn("The dummy-edge mask {} beside the output is not this graph's and could "
                     "not be removed ({}): delete it, or rebuild it with `metagraph transform "
                     "--mask-dummy --force`", stale, error.message());
    }
}

DummyMaskCounts mask_dummy_edges(DBGSuccinct *graph, size_t num_threads) {
    assert(graph);
    Timer timer;
    graph->mask_dummy_kmers(num_threads, false);

    DummyMaskCounts counts;
    counts.seconds = timer.elapsed();
    counts.edges = graph->max_index();
    counts.kmers = graph->num_nodes();
    // every sink dummy is an edge with W = $ and every such edge is masked; mark_sink_dummy_edges
    // counts them by one rank (the main dummy edge 1, the root of the source tree, excluded)
    counts.sink_dummy = graph->get_boss().mark_sink_dummy_edges();
    const uint64_t masked = counts.edges - counts.kmers;
    if (masked < counts.sink_dummy) {
        throw std::logic_error("the dummy-edge mask leaves " + std::to_string(masked)
                               + " edges masked, fewer than the " + std::to_string(counts.sink_dummy)
                               + " sink dummies");
    }
    counts.source_dummy = masked - counts.sink_dummy;
    return counts;
}

std::shared_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename,
                                                 LoadFailure on_failure) {
    for (const auto &[extension, load] : kGraphLoaders) {
        if (utils::ends_with(filename, extension))
            return load(filename, on_failure);
    }
    auto exts = utils::get_firsts<std::vector<std::string>>(kGraphLoaders);
    const std::string message = fmt::format("Cannot load graph from '{}': unrecognized file "
                                            "extension (expected one of {})",
                                            filename, fmt::join(exts, "/"));
    if (on_failure == LoadFailure::THROW) {
        throw std::runtime_error(message);
    }
    logger->error("{}", message);
    exit(1);
}

} // namespace cli
} // namespace mtg
