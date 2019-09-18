#ifndef __NODE_WEIGHTS_HPP__
#define __NODE_WEIGHTS_HPP__

#include <string>
#include <memory>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "sequence_graph.hpp"
#include "bitmap.hpp"


class DBGWeights : public DeBruijnGraph::GraphExtension {
  public:
    using node_index = typename DeBruijnGraph::node_index;
    using weight = typename sdsl::int_vector<>::value_type;

    DBGWeights(const DeBruijnGraph &graph);
    DBGWeights(const DeBruijnGraph &graph, size_t bits_per_count);
    DBGWeights(const DeBruijnGraph &graph, sdsl::int_vector<>&& weights);

    void add_kmer(const std::string&& kmer, uint32_t count);
    void add_weight(node_index i, weight w);

    void add_sequence(const std::string&& sequence,
                      bitmap *nodes_inserted = nullptr);

    void insert_node(node_index i);
    void insert_nodes(bitmap *nodes_inserted);

    template <class Bitmap>
    void remove_masked_weights(const Bitmap &mask) {
        assert(mask.size() == weights_.size());

        node_index curpos = 1;
        for (node_index i = 1; i < mask.size(); ++i) {
            if (mask[i])
                weights_[curpos++] = weights_[i];
        }
        weights_.resize(curpos);
    }

    void set_weights(sdsl::int_vector<>&& weights);

    inline weight operator[](node_index i) const {
        assert(i < weights_.size());
        return weights_[i];
    }

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;
    bool is_compatible(bool verbose = true) const;

  private:
    const DeBruijnGraph &graph_;
    sdsl::int_vector<> weights_;
    uint64_t max_weight_;

    static constexpr auto kWeightsExtension = ".weights";
};

#endif // __NODE_WEIGHTS_HPP__
