#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

namespace mtg {

namespace graph {

class DeBruijnGraph;

namespace align {
struct DBGAlignerConfig;
} // namespace align

} // namespace graph

namespace cli {

class Config;

graph::align::DBGAlignerConfig initialize_aligner_config(const Config &config,
                                                         const graph::DeBruijnGraph &graph);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
