#include "node_weights.hpp"

#include <filesystem>

#include "common/algorithms.hpp"
#include "common/utils/file_utils.hpp"

namespace mtg {
namespace graph {

namespace fs = std::filesystem;

NodeWeights::NodeWeights(uint64_t max_index, size_t bits_per_count)
      : weights_(max_index, 0, bits_per_count),
        max_weight_(sdsl::bits::lo_set[weights_.width()]) {}

NodeWeights::NodeWeights(sdsl::int_vector<>&& weights)
      : weights_(std::move(weights)),
        max_weight_(sdsl::bits::lo_set[weights_.width()]) {}

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
            = utils::make_suffix(filename_base, kWeightsExtension);
    try {
        std::unique_ptr<std::ifstream> in
                = utils::open_ifstream(weights_filename, utils::with_mmap());
        if (!in->good())
            return false;

        weights_.load(*in);
        max_weight_ = sdsl::bits::lo_set[weights_.width()];
        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

void NodeWeights::serialize(const std::string &filename_base) const {
    const auto fname = utils::make_suffix(filename_base, kWeightsExtension);

    std::ofstream outstream(fname, std::ios::binary);
    weights_.serialize(outstream);
}

void NodeWeights::serialize(sdsl::int_vector_buffer<>&& weights,
                            const std::string &filename_base) {
    const auto fname = utils::make_suffix(filename_base, kWeightsExtension);
    const std::string old_fname = weights.filename();
    weights.close(false); // close without removing the file
    fs::rename(old_fname, fname);
}

bool NodeWeights::is_compatible(const SequenceGraph &graph, bool verbose) const {
    // nodes plus dummy npos
    if (graph.max_index() + 1 == weights_.size())
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

} // namespace graph
} // namespace mtg
