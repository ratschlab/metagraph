#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <future>
#include <memory>

namespace mtg {

namespace graph {
class DeBruijnGraph;
class AnnotatedDBG;
} // namespace graph

namespace cli {

class Config;

inline constexpr size_t kDefaultMaxChunksOpen = 2000;

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config,
                         size_t max_chunks_open = kDefaultMaxChunksOpen);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);

// Kick off graph loading on a worker thread (with --pattern-build-mask, the graph's mask is
// built there too, before the future is ready: build_mask_at_load).
std::shared_future<std::shared_ptr<graph::DeBruijnGraph>>
async_load_critical_dbg(const Config &config);

/**
 * --pattern-build-mask (DESIGN-pattern-search.md §4, mask: built_at_load): when |graph| is a
 * succinct graph loaded without its dummy-edge mask, build the mask in memory as
 * `transform --mask-dummy` and `build --mask-dummy` do (mask_dummy_edges), log the counts and
 * the time, and record |graph| as one whose mask was built at load. A graph that has its mask
 * (read from its .edgemask) keeps it, and one that is not succinct has none to build (the
 * pattern search refuses it for that). Exits the process when the mask cannot be built: the
 * operator asked for it, and a server without it would answer mask_required throughout.
 * |stdout_reserved|: the caller's stdout carries its answers (`metagraph pattern`), where the
 * logger writes its info lines, so the progress is logged at trace level (stderr, with -v).
 */
void build_mask_at_load(const std::shared_ptr<graph::DeBruijnGraph> &graph,
                        bool stdout_reserved = false);

// Whether this process built |graph|'s mask at load (build_mask_at_load) rather than reading
// it from its file. |graph| is the served graph: a DBGSuccinct, or the CanonicalDBG that wraps
// a PRIMARY one.
bool mask_built_at_load(const graph::DeBruijnGraph &graph);

/**
 * Start loading the graph and the AnnotatedDBG in parallel. Returns futures
 * for both. The annotated DBG future resolves to nullptr when no annotation
 * is configured.
 */
std::pair<std::shared_future<std::shared_ptr<graph::DeBruijnGraph>>,
          std::future<std::unique_ptr<graph::AnnotatedDBG>>>
load_graph_with_async_annotation(const Config &config);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
