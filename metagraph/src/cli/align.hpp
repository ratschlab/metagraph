#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

namespace mtg {

namespace graph {
    class DeBruijnGraph;
    namespace align {
        class IDBGAligner;
    }
}

namespace cli {

class Config;

std::unique_ptr<graph::align::IDBGAligner>
build_aligner(const graph::DeBruijnGraph &graph, const Config &config);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
