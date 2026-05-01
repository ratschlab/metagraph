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

// Load the graph; if `config.infbase_annotators` is non-empty, also start
// loading the annotation on a background thread. Returns the graph (ready)
// and a future that resolves to the annotated DBG. The future is invalid
// (`.valid() == false`) when no annotation is configured.
//
// Lifetime: `config` is captured by reference by the spawned task and must
// outlive the returned future, including its destructor (which blocks on the
// background thread).
std::pair<std::shared_ptr<graph::DeBruijnGraph>,
          std::future<std::unique_ptr<graph::AnnotatedDBG>>>
load_graph_with_async_annotation(const Config &config);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
