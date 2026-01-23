#ifndef __TRANSFORM_GRAPH_HPP__
#define __TRANSFORM_GRAPH_HPP__

#include <cstddef>

namespace mtg {
namespace graph {
class DBGSuccinct;
} // namespace graph
} // namespace mtg

namespace mtg {
namespace cli {

class Config;

int transform_graph(Config *config);
void index_suffix_ranges(size_t suffix_length,
                         size_t num_threads,
                         graph::DBGSuccinct *dbg_succ);

} // namespace cli
} // namespace mtg

#endif // __TRANSFORM_GRAPH_HPP__
