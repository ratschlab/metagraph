#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

class IDBGAligner;
class DBGAlignerConfig;
class DeBruijnGraph;
class AnnotatedDBG;

namespace mtg {
namespace cli {

class Config;

DBGAlignerConfig
initialize_aligner_config(const DeBruijnGraph &graph, const Config &config);

std::unique_ptr<IDBGAligner>
build_aligner(const DeBruijnGraph &graph, const Config &config);

std::unique_ptr<IDBGAligner>
build_aligner(const DeBruijnGraph &graph, const DBGAlignerConfig &aligner_config);

std::unique_ptr<IDBGAligner>
build_masked_aligner(const AnnotatedDBG &anno_graph, const Config &config);

std::unique_ptr<IDBGAligner>
build_masked_aligner(const AnnotatedDBG &anno_graph,
                     const DBGAlignerConfig &aligner_config);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
