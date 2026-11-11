#ifndef __CLEAN_GRAPH_HPP__
#define __CLEAN_GRAPH_HPP__


#include <memory>


namespace mtg {

namespace graph {
class DeBruijnGraph;
} // namespace graph

namespace cli {

class Config;

int clean_graph(Config *config);

// Same, for a graph that is already loaded from config->fnames[0]
int clean_graph(Config *config, std::shared_ptr<graph::DeBruijnGraph> graph);

} // namespace cli
} // namespace mtg

#endif // __CLEAN_GRAPH_HPP__
