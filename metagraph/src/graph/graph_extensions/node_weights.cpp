#include "node_weights.hpp"
#include "utils.hpp"


NodeWeights::NodeWeights(uint64_t num_nodes, size_t bits_per_count)
      : weights_(num_nodes, 0, bits_per_count),
        max_weight_(utils::max_ull(weights_.width())) {}

NodeWeights::NodeWeights(sdsl::int_vector<>&& weights)
      : weights_(std::move(weights)),
        max_weight_(utils::max_ull(weights_.width())) {}

void NodeWeights::insert_nodes(const bitmap &nodes_inserted) {
    utils::insert(&weights_, nodes_inserted, 0);
}

void NodeWeights::remove_nodes(const bitmap &nodes_removed) {
    assert(nodes_removed.size() == weights_.size());

    node_index curpos = 1;
    for (node_index i = 1; i < nodes_removed.size(); ++i) {
        if (!nodes_removed[i])
            weights_[curpos++] = weights_[i];
    }
    weights_.resize(curpos);
}

bool NodeWeights::NodeWeights::load(const std::string &filename_base) {
    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;
    try {
        std::ifstream instream(weights_filename, std::ios::binary);
        if (!instream.good())
            return false;

        weights_.load(instream);
        max_weight_ = utils::max_ull(weights_.width());
        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

void NodeWeights::serialize(const std::string &filename_base) const {
    const auto weights_filename
        = utils::remove_suffix(filename_base, kWeightsExtension)
                                        + kWeightsExtension;

    std::ofstream outstream(weights_filename, std::ios::binary);
    weights_.serialize(outstream);
}

bool NodeWeights::is_compatible(const SequenceGraph &graph, bool verbose) const {
    // nodes plus dummy npos
    // TODO: fix this by implementing SequenceGraph::max_index()
    if (graph.num_nodes() + 1 == weights_.size())
        return true;

    if (verbose)
        std::cerr << "ERROR: weights file does not match number of nodes in graph"
                  << std::endl;
    return false;
}

void NodeWeights::add_weight(node_index i, weight w) {
    assert(i < weights_.size());
    uint64_t old_weight = weights_[i];
    assert(old_weight <= max_weight_);

    if (old_weight < max_weight_)
        weights_[i] = w < max_weight_ - old_weight
                    ? old_weight + w
                    : max_weight_;
}
