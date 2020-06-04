#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

class DeBruijnGraph;
class IDBGAligner;

namespace mtg {
namespace cli {

class Config;

std::unique_ptr<IDBGAligner>
build_aligner(const DeBruijnGraph &graph, const Config &config);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
