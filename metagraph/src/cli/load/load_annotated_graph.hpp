#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <string>

#include "graph/annotated_dbg.hpp"
#include "cli/config/config.hpp"


namespace mtg {
namespace cli {

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);


} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
