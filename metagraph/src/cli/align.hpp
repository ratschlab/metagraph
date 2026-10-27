#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>

namespace mtg {

namespace graph {

class DeBruijnGraph;

namespace align {
struct DBGAlignerConfig;
} // namespace align
namespace align_redone {
struct DBGAlignerConfig;
} // namespace align_redone

} // namespace graph

namespace cli {

class Config;

graph::align::DBGAlignerConfig initialize_aligner_config(const Config &config,
                                                         const graph::DeBruijnGraph &graph);

// Options for the new alignment engine (graph/alignment_redone), which is not
// yet used by the align command.
graph::align_redone::DBGAlignerConfig
initialize_redone_aligner_config(const Config &config, const graph::DeBruijnGraph &graph);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
