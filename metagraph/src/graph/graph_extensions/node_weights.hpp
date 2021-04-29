#ifndef __NODE_WEIGHTS_HPP__
#define __NODE_WEIGHTS_HPP__

#include <string>
#include <memory>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "graph/representation/base/sequence_graph.hpp"
#include "common/vectors/bitmap.hpp"


namespace mtg {
namespace graph {

class NodeWeights : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename SequenceGraph::node_index;
    using weight = typename sdsl::int_vector<>::value_type;

    NodeWeights() {}
    // initialize zero weights
    NodeWeights(uint64_t max_index, size_t bits_per_count);
    // initialize weights from existing vector
    NodeWeights(sdsl::int_vector<>&& weights);

    inline weight operator[](node_index i) const {
        assert(i < weights_.size());
        return weights_[i];
    }

    void add_weight(node_index i, weight w);

    void insert_nodes(const bitmap &nodes_inserted);
    void remove_nodes(const bitmap &nodes_removed);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;
    // serialize weights from buffer without loading all in RAM
    static void serialize(sdsl::int_vector_buffer<>&& weights,
                          const std::string &filename_base);

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    sdsl::int_vector<>& get_data() { return weights_; }

  private:
    sdsl::int_vector<> weights_;
    uint64_t max_weight_;

    static constexpr auto kWeightsExtension = ".weights";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_WEIGHTS_HPP__
