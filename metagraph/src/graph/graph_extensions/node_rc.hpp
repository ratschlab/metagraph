#ifndef __NODE_RC_HPP__
#define __NODE_RC_HPP__

#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace graph {

class NodeFirstCache;

// Maps each node in a PRIMARY-mode DeBruijnGraph to nodes adjacent to its reverse
// complement.
class NodeRC : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;

    NodeRC(const DeBruijnGraph &graph);

    void adjacent_outgoing_from_rc(node_index node,
                                   const std::function<void(node_index)> &callback,
                                   const NodeFirstCache *cache = nullptr,
                                   const std::string &spelling_hint = "") const;
    void adjacent_incoming_from_rc(node_index node,
                                   const std::function<void(node_index)> &callback,
                                   const NodeFirstCache *cache = nullptr,
                                   const std::string &spelling_hint = "") const;

    void call_outgoing_from_rc(node_index node,
                               const std::function<void(node_index, char)> &callback,
                               const NodeFirstCache *cache = nullptr,
                               const std::string &spelling_hint = "") const;
    void call_incoming_from_rc(node_index node,
                               const std::function<void(node_index, char)> &callback,
                               const NodeFirstCache *cache = nullptr,
                               const std::string &spelling_hint = "") const;

    bool load(const std::string &) { throw std::runtime_error("Loading NodeRC not possible"); }
    void serialize(const std::string &) const { throw std::runtime_error("Serializing NodeRC not possible"); }

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    const DeBruijnGraph *graph_;
};

} // namespace graph
} // namespace mtg

#endif // __NODE_RC_HPP__
