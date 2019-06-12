#include "weighted_graph.hpp"

#include <cassert>

#include "utils.hpp"
#include "dbg_succinct.hpp"


template <typename Weights>
WeightedDBG<Weights>::WeightedDBG(std::shared_ptr<DeBruijnGraph> graph, Weights&& weights)
      : graph_(std::move(graph)), weights_(std::move(weights)) {
    assert(graph_.get());
    assert(graph_->num_nodes() + 1 == weights_.size());
}

template <typename Weights>
WeightedDBG<Weights>::WeightedDBG(std::shared_ptr<DeBruijnGraph> graph)
      : graph_(std::move(graph)), weights_(Weights()) {
    assert(graph_.get());
}

template <typename Weights>
bool WeightedDBG<Weights>::load(const std::string &filename) {
    assert(graph_.get());

    if (dynamic_cast<DBGSuccinct*>(graph_.get())) {
        if (!dynamic_cast<DBGSuccinct&>(*graph_).load_without_mask(filename))
            return false;
    } else if (!graph_->load(filename)) {
        return false;
    }

    try {
        std::ifstream instream(utils::remove_suffix(filename, graph_->file_extension())
                                    + graph_->file_extension()
                                    + kWeightsExtension,
                               std::ios::binary);
        weights_.load(instream);
        return graph_->num_nodes() + 1 == weights_.size();
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << filename + kWeightsExtension << std::endl;
        return false;
    }
}

template <typename Weights>
void WeightedDBG<Weights>::serialize(const std::string &filename) const {
    graph_->serialize(filename);

    std::ofstream outstream(utils::remove_suffix(filename, graph_->file_extension())
                                + graph_->file_extension()
                                + kWeightsExtension,
                            std::ios::binary);

    weights_.serialize(outstream);
}

template class WeightedDBG<sdsl::int_vector<>>;
