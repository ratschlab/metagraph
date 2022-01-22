#ifndef __NODE_RC_HPP__
#define __NODE_RC_HPP__

#include <sdsl/dac_vector.hpp>
#include <sdsl/rrr_vector.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bit_vector.hpp"


namespace mtg {
namespace graph {

class NodeRC : public SequenceGraph::GraphExtension {
  public:
    using node_index = SequenceGraph::node_index;
    using edge_index = uint64_t;
    using mapping_type = sdsl::dac_vector_dp<sdsl::rrr_vector<>>;

    NodeRC() {};
    NodeRC(const DeBruijnGraph &graph);

    edge_index get_prefix_rc(node_index node) const {
        if (auto rank = rc_->conditional_rank1(node)) {
            return (*mapping_)[(rank - 1) * 2];
        } else {
            return 0;
        }
    }

    edge_index get_suffix_rc(node_index node) const {
        if (auto rank = rc_->conditional_rank1(node)) {
            return (*mapping_)[(rank - 1) * 2 + 1];
        } else {
            return 0;
        }
    }

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    std::shared_ptr<const bit_vector> rc_;
    std::shared_ptr<const mapping_type> mapping_;

    static constexpr auto kRCExtension = ".rc";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_RC_HPP__
