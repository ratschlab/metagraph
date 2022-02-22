#ifndef __NODE_RC_HPP__
#define __NODE_RC_HPP__

#include <sdsl/dac_vector.hpp>
#include <sdsl/rrr_vector.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace graph {

class DBGSuccinct;

// Maps each node in a PRIMARY-mode DeBruijnGraph to nodes adjacent to its reverse
// complement. When using the index construct constructor on DBGSuccinct, or when
// loading such an index, this stores a map from each DBGSuccinct to the BOSS nodes
// corresponding to the reverse complements of its k-1 prefix and k-1 suffix.
template <class Indicator = bit_vector_smart,
          class Mapping = sdsl::dac_vector_dp<sdsl::rrr_vector<>>>
class NodeRC : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;

    // Use set_graph to set the graph pointer after using the default constructor
    NodeRC() : graph_(nullptr) {};

    // Construct a NodeRC index
    NodeRC(const DBGSuccinct &graph);

    void set_graph(const DeBruijnGraph &graph) {
        assert(!rc_.size() || is_compatible(graph));
        graph_ = &graph;
    }

    void call_outgoing_nodes_from_rc(node_index node, const std::function<void(node_index)> &callback) const;
    void call_incoming_nodes_from_rc(node_index node, const std::function<void(node_index)> &callback) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    const DeBruijnGraph *graph_;

    // Indicates whether a reverse complement prefix/suffix exists for each graph node.
    Indicator rc_;

    // Vector storing the indices of the reverse complements of the node prefix and
    // suffix (if they exist) in an interleaved form.
    // mapping_.size() == rc_->num_set_bits() * 2
    Mapping mapping_;

    static constexpr auto kRCExtension = ".rc";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_RC_HPP__
