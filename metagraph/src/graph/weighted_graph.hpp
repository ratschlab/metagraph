#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include <string>
#include <memory>
#include <vector>

#include "sequence_graph.hpp"
#include "utils.hpp"
#include "int_vector.hpp"


template <typename Weights = sdsl::int_vector<>>
class DBGWeights : public DBGExtension<DeBruijnGraph> {
  public:
    using node_index = typename DeBruijnGraph::node_index;
    using weight = typename Weights::value_type;

    DBGWeights() = default;
    DBGWeights(Weights&& weights) : weights_(std::move(weights)) {};
    DBGWeights(const DeBruijnGraph &graph, const std::string &filename_base) { load(graph, filename_base); };

    virtual void add_kmer(const DeBruijnGraph &graph, const std::string&& kmer, uint32_t count) {
        auto node = graph.kmer_to_node(kmer);
        add_weight(node, count);
    };

    virtual void add_sequence(const DeBruijnGraph &graph,
                              const std::string&& sequence,
                              bit_vector_dyn *nodes_inserted = nullptr) {

        if (nodes_inserted)
            utils::insert(&weights_, *nodes_inserted, 0);

        auto k = graph.get_k();
        for (size_t i = 0; i + k <= sequence.size(); ++i) {
            add_kmer(graph, sequence.substr(i, k), 1);
        }
    };

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
    };

    template <class Vector>
    void remove_masked_weights(const Vector &mask) {
        node_index curpos = 1;
        assert(mask.size() == weights_.size());
        for (node_index i = 1; i < mask.size(); ++i) {
            if (mask[i])
                weights_[curpos++] = weights_[i];
        }
        weights_.resize(curpos);
    };

    virtual void set_weights(Weights&& weights) { weights_ = std::move(weights); };

    virtual weight get_weight(node_index i) const {
        assert(i < weights_.size());
        return weights_[i];
    };

    virtual bool load(const DeBruijnGraph &graph, const std::string &filename_base);
    virtual void serialize(const DeBruijnGraph &graph, const std::string &filename_base) const;

    static bool has_file(const DeBruijnGraph &graph, const std::string &filename_base);

  private:
    Weights weights_;

    virtual void add_weight(node_index i, weight w) {
        assert(i < weights_.size());
        weight prev = weights_[i];
        weights_[i] += w;
        if (weights_[i] < prev) {
            std::cerr << "Warning: weighted graph entry " << i << " has"
                      << " overflowed:" << prev << " + " << w << " = "
                      << weights_[i] << std::endl;
        }
    };

    static constexpr auto kWeightsExtension = ".weights";
};

template <typename Weights>
bool DBGWeights<Weights>::load(const DeBruijnGraph &graph, const std::string &filename_base) {

    const auto weights_filename = utils::remove_suffix(filename_base, graph.file_extension())
                                        + graph.file_extension()
                                        + kWeightsExtension;
    try {
        std::ifstream instream(weights_filename, std::ios::binary);
        this->weights_.load(instream);

        if (graph.num_nodes() + 1 != this->weights_.size()) {
            std::cerr << "ERROR: weights file does not match number of nodes in graph "
                      << weights_filename << std::endl;
            return false;
        }

        return true;
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

template <typename Weights>
void DBGWeights<Weights>::serialize(const DeBruijnGraph &graph, const std::string &filename_base) const {

    const auto weights_filename = utils::remove_suffix(filename_base, graph.file_extension())
                                        + graph.file_extension()
                                        + kWeightsExtension;

    std::ofstream outstream(weights_filename, std::ios::binary);
    this->weights_.serialize(outstream);
}

template <typename Weights>
bool DBGWeights<Weights>::has_file(const DeBruijnGraph &graph, const std::string &filename_base) {
    const auto weights_filename = utils::remove_suffix(filename_base, graph.file_extension())
                                        + graph.file_extension()
                                        + kWeightsExtension;
    std::ifstream instream(weights_filename, std::ios::binary);
    return instream.good();
}

#endif // __WEIGHTED_GRAPH_HPP__
