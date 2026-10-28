#ifndef __ALIGN_GRAPH_HPP__
#define __ALIGN_GRAPH_HPP__

#include <memory>
#include <optional>
#include <string_view>

#include "graph/alignment_redone/aln_match.hpp"
#include "graph/annotated_dbg.hpp"

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

// Options for the opt-in alignment engine (graph/alignment_redone), used by
// `align --align-engine redone`.
graph::align_redone::DBGAlignerConfig
initialize_redone_aligner_config(const Config &config, const graph::DeBruijnGraph &graph);

struct RedoneBestAlignment {
    graph::align_redone::Alignment alignment;
    // labels of the alignment's first label class (empty without annotation)
    std::vector<std::string> labels;
};

// Align |seq| with the opt-in engine and return its best-scoring alignment.
// |graph| must be the graph to align to directly (e.g., a PRIMARY graph
// wrapped into CanonicalDBG), |annotator| may be null.
std::optional<RedoneBestAlignment>
align_redone_best(const graph::DeBruijnGraph &graph,
                  const graph::AnnotatedDBG::Annotator *annotator,
                  const graph::align_redone::DBGAlignerConfig &config,
                  std::string_view seq,
                  bool connect_anchors,
                  bool extend_chains);

int align_to_graph(Config *config);

} // namespace cli
} // namespace mtg

#endif // __ALIGN_GRAPH_HPP__
