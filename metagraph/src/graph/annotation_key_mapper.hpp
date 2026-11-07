#pragma once

#include <span>
#include <string_view>
#include <vector>

#include "representation/base/sequence_graph.hpp"

namespace mtg::graph {

class CanonicalDBG;
class DBGSSHash;

// Map oriented graph nodes to the node IDs whose rows carry their annotation.
// IDs remain one-based here; AnnotatedDBG::graph_to_anno_index converts to rows.
class AnnotationKeyMapper {
  public:
    using node_index = DeBruijnGraph::node_index;
    explicit AnnotationKeyMapper(const DeBruijnGraph &graph);

    bool needs_sequence() const {
        return mode_ == DeBruijnGraph::CANONICAL && !canonical_ && !sshash_;
    }
    // For a nonzero node, supply its k-mer spelling when needs_sequence() is true.
    node_index key(node_index node, std::string_view kmer = {}) const {
        if (!node || mode_ == DeBruijnGraph::BASIC)
            return node;
        return canonical_key(node, kmer);
    }
    // Supply the path's spelling when needs_sequence() is true. Missing nodes
    // retain npos, and returned keys correspond to the input positions.
    std::vector<node_index> keys(std::span<const node_index> nodes,
                                 std::string_view sequence = {}) const;

  private:
    node_index canonical_key(node_index node, std::string_view kmer) const;
    const DeBruijnGraph &graph_;
    DeBruijnGraph::Mode mode_;
    const CanonicalDBG *canonical_;
    const DBGSSHash *sshash_;
};

} // namespace mtg::graph
