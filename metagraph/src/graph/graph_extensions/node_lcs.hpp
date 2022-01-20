#ifndef __NODE_LCS_HPP__
#define __NODE_LCS_HPP__

#include "common/vectors/wavelet_tree.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include <sdsl/int_vector_buffer.hpp>


namespace mtg {
namespace graph {

class NodeLCS : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename SequenceGraph::node_index;
    using edge_index = uint64_t;
    using lcs = uint64_t;

    NodeLCS() {};
    NodeLCS(const DeBruijnGraph &graph);

    inline lcs operator[](node_index i) const { return (*lcs_)[i]; }

    std::pair<edge_index, edge_index> expand(edge_index first, edge_index last, size_t width, size_t step = 1) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    static void serialize(sdsl::int_vector_buffer<>&& lcs,
                          const std::string &filename_base);

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    wavelet_tree& data() { return *lcs_; }

  private:
    std::unique_ptr<wavelet_tree> lcs_;

    static constexpr auto kLCSExtension = ".lcs";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_LCS_HPP__
