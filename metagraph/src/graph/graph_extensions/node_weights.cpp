#include "node_weights.hpp"
#include "utils.hpp"


DBGWeights::DBGWeights(uint64_t num_nodes, size_t bits_per_count)
      : weights_(sdsl::int_vector<>(num_nodes, 0, bits_per_count)),
        max_weight_(~uint64_t(0) >> (64 - weights_.width())) {}

DBGWeights::DBGWeights(sdsl::int_vector<>&& weights)
      : weights_(std::move(weights)),
        max_weight_(~uint64_t(0) >> (64 - weights_.width())) {}

void DBGWeights::insert_nodes(const bitmap &nodes_inserted) {
    utils::insert(&weights_, nodes_inserted, 0);
}

void DBGWeights::remove_nodes(const bitmap &nodes_removed) {
    assert(nodes_removed.size() == weights_.size());

    node_index curpos = 1;
    for (node_index i = 1; i < nodes_removed.size(); ++i) {
        if (!nodes_removed[i])
            weights_[curpos++] = weights_[i];
    }
    weights_.resize(curpos);
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

bool DBGWeights::is_compatible(const SequenceGraph &graph, bool verbose) const {
    if (!dynamic_cast<const DeBruijnGraph*>(&graph)) {
        std::cerr << "ERROR: DBGWeights can be used only with de Bruijn graph"
                  << std::endl;
    }

    if (graph.num_nodes() + 1 == weights_.size())
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
