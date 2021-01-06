#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

namespace mtg {

namespace graph {
    class DeBruijnGraph;
    namespace align {
        class IDBGAligner;
        class DBGAlignerConfig;
    }

    class AnnotatedDBG;
}

namespace cli {

class Config;

graph::align::DBGAlignerConfig
initialize_aligner_config(size_t k, const Config &config);

std::unique_ptr<graph::align::IDBGAligner>
build_aligner(const graph::DeBruijnGraph &graph, const Config &config);

std::unique_ptr<graph::align::IDBGAligner>
build_aligner(const graph::DeBruijnGraph &graph, const graph::align::DBGAlignerConfig &aligner_config);

std::unique_ptr<graph::align::IDBGAligner>
build_labeled_aligner(const graph::AnnotatedDBG &anno_graph, const graph::align::DBGAlignerConfig &aligner_config);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
