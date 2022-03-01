#ifndef __NODE_RC_HPP__
#define __NODE_RC_HPP__

#include <sdsl/dac_vector.hpp>
#include <sdsl/rrr_vector.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace graph {

// Maps each node in a PRIMARY-mode DeBruijnGraph to nodes adjacent to its reverse
// complement. When using the index construct constructor on DBGSuccinct, or when
// loading such an index, this stores a map from each DBGSuccinct to the BOSS nodes
// corresponding to the reverse complements of its k-1 prefix and k-1 suffix.
class NodeRC : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;

    // If construct_index is false, then an index must be loaded to take advantage
    // of the index for call_*_from_rc calls
    NodeRC(const DeBruijnGraph &graph, bool construct_index = false);

    void call_outgoing_from_rc(node_index node, const std::function<void(node_index)> &callback) const;
    void call_incoming_from_rc(node_index node, const std::function<void(node_index)> &callback) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    const DeBruijnGraph *graph_;

    typedef bit_vector_smart Indicator;
    typedef sdsl::dac_vector_dp<sdsl::rrr_vector<>> Mapping;

    // Indicates whether a reverse complement prefix/suffix exists for each graph node.
    // If the index is initialized, then rc_.size() == graph_->max_index() + 1
    Indicator rc_;

    // Vector storing the indices of the reverse complements of the node prefix and
    // suffix (if they exist) in an interleaved form.
    // If the index is initialized, then mapping_.size() == rc_->num_set_bits() * 2
    Mapping mapping_;

    static constexpr auto kRCExtension = ".rc_adj";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_RC_HPP__
