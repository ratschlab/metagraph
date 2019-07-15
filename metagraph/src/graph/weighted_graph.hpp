#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include <string>
#include <memory>
#include <vector>

#include "sequence_graph.hpp"
#include "utils.hpp"


template <typename Weights = sdsl::int_vector<>>
class DBGWeights : public DBGExtension<DeBruijnGraph> {
  public:
    using node_index = typename DeBruijnGraph::node_index;
    using weight = typename Weights::value_type;

    DBGWeights() = default;
    DBGWeights(Weights&& weights) : weights_(std::move(weights)) {};
    DBGWeights(const DeBruijnGraph &graph, const std::string &filename_base) { load(graph, filename_base); };

    virtual void set_weights(Weights&& weights) { weights_ = std::move(weights); };
    virtual weight get_weight(node_index i) const { return weights_[i]; };

    virtual bool load(const DeBruijnGraph &graph, const std::string &filename_base);
    virtual void serialize(const DeBruijnGraph &graph, const std::string &filename_base) const;

    static bool has_file(const DeBruijnGraph &graph, const std::string &filename_base);

  private:
    Weights weights_;

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
        return graph.num_nodes() + 1 == this->weights_.size();
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}

template <typename Weights>
void DBGWeights<Weights>::serialize(const DeBruijnGraph &graph, const std::string &filename_base) const {

    std::ofstream outstream(utils::remove_suffix(filename_base, graph.file_extension())
                                + graph.file_extension()
                                + kWeightsExtension,
                            std::ios::binary);

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
