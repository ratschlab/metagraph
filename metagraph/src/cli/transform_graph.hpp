#ifndef __TRANSFORM_GRAPH_HPP__
#define __TRANSFORM_GRAPH_HPP__

#include <cstddef>

namespace mtg::graph {
class DBGSuccinct;
} // namespace mtg::graph

namespace mtg {
namespace cli {

class Config;

int transform_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __TRANSFORM_GRAPH_HPP__
