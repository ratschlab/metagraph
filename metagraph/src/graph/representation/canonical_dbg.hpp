#ifndef __CANONICAL_DBG_HPP__
#define __CANONICAL_DBG_HPP__

#include <cassert>
#include <array>

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "common/vector.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"


namespace mtg {
namespace graph {

/**
 * CanonicalDBG is a wrapper which acts like a canonical-mode DeBruijnGraph, but
 * uses a PRIMARY DeBruijnGraph (constructed from primary contigs).
 */
class CanonicalDBG : public DBGWrapper<DeBruijnGraph> {
  public:
    explicit CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph,
                          size_t cache_size = 100'000);

    CanonicalDBG(const CanonicalDBG &canonical)
          : CanonicalDBG(canonical.graph_, canonical.cache_size_) {}

    /**
     * Added methods
     */
    bool operator==(const CanonicalDBG &other) const { return *graph_ == *other.graph_; }

    void reverse_complement(std::string &seq, std::vector<node_index> &path) const;
    node_index reverse_complement(node_index node) const;

    inline node_index get_base_node(node_index node) const {
        assert(node <= offset_ * 2);
        return node <= offset_ ? node : node - offset_;
    }

    /**
     * Methods from DeBruijnGraph
     */
    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override final;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override final;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override final {
        call_outgoing_kmers(kmer, get_node_sequence(kmer), callback);
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override final {
        call_incoming_kmers(kmer, get_node_sequence(kmer), callback);
    }

    void call_outgoing_kmers(node_index kmer,
                             const std::string &spelling_hint,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index kmer,
                             const std::string &spelling_hint,
                             const IncomingEdgeCallback &callback) const;

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override final;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final;

    virtual uint64_t num_nodes() const override final { return graph_->num_nodes() * 2; }
    virtual uint64_t max_index() const override final { return graph_->max_index() * 2; }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override final;

    virtual Mode get_mode() const override final { return CANONICAL; }
    virtual node_index traverse(node_index node, char next_char) const override final;
    virtual node_index traverse_back(node_index node, char prev_char) const override final;

    virtual size_t outdegree(node_index) const override final;
    virtual size_t indegree(node_index) const override final;
    virtual bool has_multiple_outgoing(node_index node) const override final;
    virtual bool has_single_incoming(node_index node) const override final;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override final;
    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const override final;

    virtual bool operator==(const DeBruijnGraph &other) const override final;

  private:
    const size_t cache_size_;

    // cache whether a given node is a palindrome (it's equal to its reverse complement)
    mutable caches::fixed_sized_cache<node_index, bool,
                                      caches::LRUCachePolicy<node_index>> is_palindrome_cache_;

    const size_t offset_;
    const bool k_odd_;
    bool has_sentinel_;

    std::array<size_t, 256> alphabet_encoder_;

    // a thread-safe 0-size NodeFirstCache to use as a fallback if another NodeFirstCache
    // extension is not available
    mutable std::unique_ptr<NodeFirstCache> fallback_cache_;

    // Find incoming nodes that are on the reverse complement strand of node
    void adjacent_incoming_rc_strand(node_index node,
                                     const std::string &spelling_hint,
                                     const std::function<void(node_index, char)> &callback) const;

    // Find outgoing nodes that are on the reverse complement strand of node
    void adjacent_outgoing_rc_strand(node_index node,
                                     const std::string &spelling_hint,
                                     const std::function<void(node_index, char)> &callback) const;

    // Can only be called when the base graph is DBGSuccinct
    NodeFirstCache& get_cache() const;
};

} // namespace graph
} // namespace mtg

#endif // __CANONICAL_DBG_HPP__
