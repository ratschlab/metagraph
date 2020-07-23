#ifndef __CANONICAL_DBG_HPP__
#define __CANONICAL_DBG_HPP__

#include <mutex>
#include <cassert>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class CanonicalDBG : public DeBruijnGraph {
  public:
    CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph, size_t cache_size_ = 100'000);
    CanonicalDBG(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size_ = 100'000);

    virtual ~CanonicalDBG() {}

    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion = [](node_index) {}) override;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override;

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override;

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override;

    virtual uint64_t num_nodes() const override;
    virtual uint64_t max_index() const override { return graph_.max_index() * 2; }

    virtual bool load(const std::string &) override {
        throw std::runtime_error("Not implemented");
    }

    virtual void serialize(const std::string &) const override {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string file_extension() const override { return graph_.file_extension(); }

    virtual const std::string& alphabet() const override { return graph_.alphabet(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override;

    virtual size_t get_k() const override { return graph_.get_k(); }

    virtual bool is_canonical_mode() const override { return true; }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const override;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const override;

    virtual size_t outdegree(node_index) const override;
    virtual size_t indegree(node_index) const override;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const override;

    virtual const DeBruijnGraph& get_graph() const { return graph_; }

    virtual bool operator==(const CanonicalDBG &other) const {
        return graph_ == other.graph_;
    }

    virtual bool operator==(const DeBruijnGraph &other) const override;

    inline node_index reverse_complement(node_index node) const {
        assert(node);
        assert(node <= offset_ * 2);
        return node > offset_ ? node - offset_ : node + offset_;
    }

    inline node_index get_base_node(node_index node) const {
        assert(node);
        assert(node <= offset_ * 2);
        return node > offset_ ? node - offset_ : node;
    }

    inline bool is_rev_comp(node_index node) const {
        assert(node);
        assert(node <= offset_ * 2);
        return node > offset_;
    }

  private:
    std::shared_ptr<const DeBruijnGraph> const_graph_ptr_;
    const DeBruijnGraph &graph_ = *const_graph_ptr_;
    size_t offset_;

    std::shared_ptr<DeBruijnGraph> graph_ptr_;

    std::array<size_t, 256> alph_map_;

    inline node_index set_offset(node_index node) const {
        assert(node <= offset_);
        assert(node);
        assert(is_rev_comp(node + offset_));
        return node + offset_;
    }

    inline node_index unset_offset(node_index node) const {
        assert(node > offset_);
        assert(node <= offset_ * 2);
        assert(!is_rev_comp(node - offset_));
        return node - offset_;
    }

    mutable caches::fixed_sized_cache<node_index, std::vector<node_index>,
                                      caches::LRUCachePolicy<node_index>> child_node_cache_;
    mutable caches::fixed_sized_cache<node_index, std::vector<node_index>,
                                      caches::LRUCachePolicy<node_index>> parent_node_cache_;
    mutable std::mutex cache_mutex_;
};

} // namespace graph
} // namespace mtg

#endif // __CANONICAL_DBG_HPP__
