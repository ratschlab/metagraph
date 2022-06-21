#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__


#include <memory>

namespace mtg {

namespace graph {

class DeBruijnGraph;
class AnnotatedDBG;

} // namespace graph

namespace cli {

class Config;

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config,
                         size_t max_chunks_open = 2000);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
