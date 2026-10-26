#include "annotation_key_mapper.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "representation/canonical_dbg.hpp"
#include "representation/hash/dbg_sshash.hpp"

namespace mtg::graph {

AnnotationKeyMapper::AnnotationKeyMapper(const DeBruijnGraph &graph)
      : graph_(graph), mode_(graph.get_mode()),
        canonical_(dynamic_cast<const CanonicalDBG*>(&graph)),
        sshash_(dynamic_cast<const DBGSSHash*>(&graph)) {
    if (mode_ == DeBruijnGraph::PRIMARY)
        throw std::invalid_argument("Wrap PRIMARY graphs in CanonicalDBG before mapping annotations");
}

auto AnnotationKeyMapper::canonical_key(node_index node, std::string_view kmer) const -> node_index {
    if (canonical_)
        return canonical_->get_base_node(node);
    if (sshash_)
        return std::min(node, sshash_->reverse_complement(node));
    assert(kmer.size() == graph_.get_k());
    node_index result = DeBruijnGraph::npos;
    graph_.map_to_nodes(kmer, [&](node_index key) { result = key; });
    return result;
}

auto AnnotationKeyMapper::keys(std::span<const node_index> nodes,
                                std::string_view sequence) const -> std::vector<node_index> {
    std::vector<node_index> result(nodes.size());
    if (needs_sequence()) {
        assert(sequence.size() == nodes.size() + graph_.get_k() - 1);
        size_t i = 0;
        graph_.map_to_nodes(sequence, [&](node_index key) { result[i++] = key; });
        assert(i == nodes.size());
    } else {
        std::transform(nodes.begin(), nodes.end(), result.begin(),
                       [&](node_index node) { return key(node); });
    }
    return result;
}

} // namespace mtg::graph
