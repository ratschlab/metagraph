#ifndef __NODE_RC_HPP__
#define __NODE_RC_HPP__

#include <sdsl/dac_vector.hpp>
#include <sdsl/rrr_vector.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"

namespace mtg {
namespace graph {

class DBGSuccinct;

class INodeRC : public SequenceGraph::GraphExtension {
  public:
    virtual uint64_t get_prefix_rc(SequenceGraph::node_index node) const = 0;
    virtual uint64_t get_suffix_rc(SequenceGraph::node_index node) const = 0;
};

// Maps each node in a PRIMARY-mode DBGSuccinct to the BOSS nodes corresponding
// to its k-1 prefix or suffix.
template <class Indicator = bit_vector_smart,
          class Mapping = sdsl::dac_vector_dp<sdsl::rrr_vector<>>>
class NodeRC : public INodeRC {
  public:
    using node_index = SequenceGraph::node_index;
    using edge_index = uint64_t;

    NodeRC() {};
    NodeRC(const DBGSuccinct &graph);

    edge_index get_prefix_rc(node_index node) const {
        if (auto rank = rc_.conditional_rank1(node))
            return mapping_[(rank - 1) * 2];

        return 0;
    }

    edge_index get_suffix_rc(node_index node) const {
        if (auto rank = rc_.conditional_rank1(node))
            return mapping_[(rank - 1) * 2 + 1];

        return 0;
    }

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
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
