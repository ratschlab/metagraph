#include "node_weights.hpp"
#include "utils.hpp"


DBGWeights::DBGWeights(const DeBruijnGraph &graph) : graph_(graph) {}

DBGWeights::DBGWeights(const DeBruijnGraph &graph, size_t bits_per_count)
      : graph_(graph),
        weights_(sdsl::int_vector<>(graph.num_nodes() + 1, 0, bits_per_count)),
        max_weight_(~uint64_t(0) >> (64 - weights_.width())) {}

DBGWeights::DBGWeights(const DeBruijnGraph &graph, sdsl::int_vector<>&& weights)
      : graph_(graph),
        weights_(std::move(weights)),
        max_weight_(~uint64_t(0) >> (64 - weights_.width())) {}

void DBGWeights::insert_node(node_index i) {
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

void DBGWeights::insert_nodes(bitmap *nodes_inserted) {
    utils::insert(&weights_, *nodes_inserted, 0);
}

void DBGWeights::set_weights(sdsl::int_vector<>&& weights) {
    weights_ = std::move(weights);
    max_weight_ = ~uint64_t(0) >> (64 - weights_.width());
}

bool DBGWeights::DBGWeights::load(const std::string &filename_base) {
    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;
    try {
        std::ifstream instream(weights_filename, std::ios::binary);
        if (!instream.good())
            return false;
        weights_.load(instream);
        max_weight_ = (~uint64_t(0) >> (64 - weights_.width()));
        return true;
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

void DBGWeights::serialize(const std::string &filename_base) const {
    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;

    std::ofstream outstream(weights_filename, std::ios::binary);
    weights_.serialize(outstream);
}

bool DBGWeights::is_compatible(bool verbose) const {
    if (graph_.num_nodes() + 1 == weights_.size())
        return true;

    if (verbose)
        std::cerr << "ERROR: weights file does not match number of nodes in graph"
                  << std::endl;
    return false;
}

void DBGWeights::add_weight(node_index i, weight w) {
    assert(i < weights_.size());
    uint64_t old_weight = weights_[i];
    assert(old_weight <= max_weight_);

    if (old_weight < max_weight_)
        weights_[i] = old_weight + w < max_weight_
                    ? old_weight + w
                    : max_weight_;
}
