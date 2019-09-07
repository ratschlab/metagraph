#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include <string>
#include <memory>
#include <vector>

#include "sequence_graph.hpp"
#include "utils.hpp"
#include "int_vector.hpp"


template <typename Weights = sdsl::int_vector<>>
class DBGWeights : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename DeBruijnGraph::node_index;
    using weight = typename Weights::value_type;

    DBGWeights(Weights&& weights = {})
          : weights_(std::move(weights)),
            max_weight_(~uint64_t(0) >> (64 - weights_.width())) {}

    virtual void add_kmer(const DeBruijnGraph &graph,
                          const std::string&& kmer,
                          uint32_t count) {
        auto node = graph.kmer_to_node(kmer);
        add_weight(node, count);
    }

    virtual void add_sequence(const DeBruijnGraph &graph,
                              const std::string&& sequence,
                              bit_vector_dyn *nodes_inserted = nullptr) {
        if (nodes_inserted)
            utils::insert(&weights_, *nodes_inserted, 0);

        graph.map_to_nodes(sequence, [&](auto node) { add_weight(node, 1); });
    }

    virtual void insert_node(node_index i) {
        if (weights_.empty()) {
            weights_.resize(1);
            weights_[0] = 0;
        }

        assert(i <= weights_.size());
        weights_.resize(weights_.size() + 1);
        node_index j = weights_.size() - 1;

        std::copy_backward(weights_.begin() + i, weights_.begin() + j, weights_.end());

        weights_[i] = 0;
    }

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

    virtual void set_weights(Weights&& weights) {
        weights_ = std::move(weights);
        max_weight_ = ~uint64_t(0) >> (64 - weights_.width());
    }

    virtual weight operator[](node_index i) const {
        assert(i < weights_.size());
        return weights_[i];
    }

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;
    virtual bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    Weights weights_;
    uint64_t max_weight_;

    virtual void add_weight(node_index i, weight w) {
        assert(i < weights_.size());
        assert(old_weight <= max_weight_);

        uint64_t old_weight = weights_[i];

        if (old_weight < max_weight_)
            weights_[i] = old_weight + w < max_weight_
                        ? old_weight + w
                        : max_weight_;
    }

    static constexpr auto kWeightsExtension = ".weights";
};

template <typename Weights>
bool DBGWeights<Weights>::load(const std::string &filename_base) {

    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;
    try {
        std::ifstream instream(weights_filename, std::ios::binary);
        weights_.load(instream);
        max_weight_ = (~uint64_t(0) >> (64 - weights_.width()));
        return true;
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

template <typename Weights>
void DBGWeights<Weights>::serialize(const std::string &filename_base) const {

    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;

    std::ofstream outstream(weights_filename, std::ios::binary);
    weights_.serialize(outstream);
}

template <typename Weights>
bool DBGWeights<Weights>::is_compatible(const SequenceGraph &graph,
                                        bool verbose) const {
    if (graph.num_nodes() + 1 == weights_.size())
        return true;

    if (verbose)
        std::cerr << "ERROR: weights file does not match number of nodes in graph"
                  << std::endl;
    return false;
}

#endif // __WEIGHTED_GRAPH_HPP__
