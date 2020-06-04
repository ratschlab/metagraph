#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <string>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/masked_graph.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

std::unique_ptr<AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                         const Config &config);

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config);

std::unique_ptr<MaskedDeBruijnGraph>
mask_graph(const AnnotatedDBG &anno_graph, Config *config);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
