#ifndef __NODE_WEIGHTS_HPP__
#define __NODE_WEIGHTS_HPP__

#include <string>
#include <memory>

#include <sdsl/int_vector.hpp>

#include "sequence_graph.hpp"
#include "bitmap.hpp"


class DBGWeights : public DeBruijnGraph::GraphExtension {
  public:
    using node_index = typename DeBruijnGraph::node_index;
    using weight = typename sdsl::int_vector<>::value_type;

    DBGWeights() {}
    // initialize zero weights
    DBGWeights(uint64_t num_nodes, size_t bits_per_count);
    // initialize weights from existing vector
    DBGWeights(sdsl::int_vector<>&& weights);

    inline weight operator[](node_index i) const {
        assert(i < weights_.size());
        return weights_[i];
    }

    void add_weight(node_index i, weight w);

    void insert_nodes(bitmap *nodes_inserted);
    // remove all weight elements except those marked with 1 in |mask|
    void remove_unmasked_weights(const bitmap &mask);

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    sdsl::int_vector<> weights_;
    uint64_t max_weight_;

    static constexpr auto kWeightsExtension = ".weights";
};

#endif // __NODE_WEIGHTS_HPP__
