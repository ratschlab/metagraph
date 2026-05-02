#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <future>
#include <memory>

namespace mtg {

namespace graph {
class DeBruijnGraph;
class AnnotatedDBG;
} // namespace graph

namespace annot {
template <typename LabelType>
class MultiLabelAnnotation;
} // namespace annot

namespace cli {

class Config;

inline constexpr size_t kDefaultMaxChunksOpen = 2000;

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config,
                         size_t max_chunks_open = kDefaultMaxChunksOpen);

std::unique_ptr<annot::MultiLabelAnnotation<std::string>>
load_annotation(std::shared_ptr<graph::DeBruijnGraph> graph,
                const Config &config,
                size_t max_chunks_open = 2000);

std::unique_ptr<annot::MultiLabelAnnotation<std::string>>
load_annotation(std::shared_future<std::shared_ptr<graph::DeBruijnGraph>> graph,
                const Config &config,
                size_t max_chunks_open = 2000);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);

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
