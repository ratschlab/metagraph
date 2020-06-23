#ifndef __BUILD_GRAPH_HPP__
#define __BUILD_GRAPH_HPP__


namespace mtg {
namespace cli {

class Config;

int build_graph(Config *config);

int concatenate_graph_chunks(Config *config);

} // namespace cli
} // namespace mtg

#endif // __BUILD_GRAPH_HPP__
