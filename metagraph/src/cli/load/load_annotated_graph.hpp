#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <string>
#include <vector>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/masked_graph.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);

std::unique_ptr<graph::MaskedDeBruijnGraph>
mask_graph(const graph::AnnotatedDBG &anno_graph, Config *config);

void
call_masked_graphs(const graph::AnnotatedDBG &anno_graph, Config *config,
                   const std::function<void(const graph::MaskedDeBruijnGraph&,
                                            const std::string& /* header */)> &callback,
                   size_t num_parallel_graphs_masked = 1,
                   size_t num_threads_per_graph = 1);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
