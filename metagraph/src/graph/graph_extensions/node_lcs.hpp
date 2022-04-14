#ifndef __NODE_LCS_HPP__
#define __NODE_LCS_HPP__

#include "common/vectors/wavelet_tree.hpp"
#include "common/vector.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

class NodeLCS : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename SequenceGraph::node_index;
    using edge_index = boss::BOSS::edge_index;
    using lcs = wavelet_tree::TAlphabet;
    using wavelet_tree_t = wavelet_tree_small;

    NodeLCS() {}
    NodeLCS(const DeBruijnGraph &graph);

    inline lcs operator[](node_index i) const { return (*lcs_)[i]; }

    std::pair<edge_index, edge_index> expand(edge_index first, edge_index last, size_t width, size_t step = 1) const;
    void call_contractions(edge_index first,
                           edge_index last,
                           size_t width,
                           const std::function<void(edge_index, edge_index)> &callback) const;

    std::pair<edge_index, edge_index>
    call_parallel_edges(node_index node,
                        size_t suffix_length,
                        const std::function<void(node_index, node_index)> &callback) const;

    std::pair<edge_index, edge_index>
    call_edges(const std::pair<edge_index, edge_index> &node_range,
               char c,
               const std::function<void(node_index, node_index)> &callback) const;

    std::vector<std::pair<edge_index, edge_index>>
    call_edges(const std::pair<edge_index, edge_index> &node_range,
               const std::function<void(node_index, const SmallVector<node_index>&)> &callback) const;

    void set_graph(const DBGSuccinct &graph);
    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    const wavelet_tree* data() const { return lcs_.get(); }

  private:
    const DBGSuccinct *graph_;
    std::unique_ptr<wavelet_tree> lcs_;

    static constexpr auto kLCSExtension = ".lcs";
};

} // namespace graph
} // namespace mtg

#endif // __NODE_LCS_HPP__
