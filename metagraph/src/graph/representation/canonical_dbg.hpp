#ifndef __CANONICAL_DBG_HPP__
#define __CANONICAL_DBG_HPP__

#include <cassert>
#include <array>

#include "graph/representation/base/dbg_wrapper.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/caches.hpp"

namespace mtg {
namespace graph {

/**
 * CanonicalDBG is a wrapper which acts like a canonical-mode DeBruijnGraph, but
 * uses a non-canonical DeBruijnGraph as the underlying storage.
 */
class CanonicalDBG : public DBGWrapper<DeBruijnGraph> {
  public:
    template <typename Graph>
    explicit CanonicalDBG(Graph&& graph, size_t cache_size = 100'000);

    // caches cannot be resized or moved, so disable these constructors
    CanonicalDBG& operator=(const CanonicalDBG &canonical) = delete;
    CanonicalDBG(CanonicalDBG&&) = delete;
    CanonicalDBG& operator=(CanonicalDBG&&) = delete;

    virtual ~CanonicalDBG() {}

    /**
     * Added methods
     */
    bool operator==(const CanonicalDBG &other) const { return *graph_ == *other.graph_; }

    void reverse_complement(std::string &seq, std::vector<node_index> &path) const;
    node_index reverse_complement(node_index node) const;

    inline node_index get_base_node(node_index node) const {
        assert(node);
        assert(node <= offset_ * 2);
        return node > offset_ ? node - offset_ : node;
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
    // Guarantees that nodes are called in the same order as the input sequence.
    // All k-mers for which the skip condition is satisfied are skipped.
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; },
                                           const std::function<bool()> &skip = [](){ return false; }) const override final;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override final;

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override final;

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;

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

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override final;
    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const override final;

    virtual bool operator==(const DeBruijnGraph &other) const override final;

    virtual void serialize(const std::string &filename_base) const override final {
        graph_->serialize(filename_base);
    }

    virtual bool load(const std::string &filename_base) override final {
        if (const_cast<DeBruijnGraph*>(graph_.get())->load(filename_base)) {
            flush();
            return true;
        }

        return false;
    }

    virtual void set_graph(std::shared_ptr<const DeBruijnGraph> graph) override final {
        graph_ = graph;
        flush();
    }

    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion) override final {
        const_cast<DeBruijnGraph*>(graph_.get())->add_sequence(sequence, on_insertion);
        flush();
    }

  private:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    // cache whether a given node is a palindrome (it's equal to its reverse complement)
    mutable common::LRUCache<node_index, bool> is_palindrome_cache_;

    // cache the BOSS node corresponding to the reverse complement of the k - 1 prefix,
    // and the number of matching characters
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_prev_cache_;

    // cache the BOSS node corresponding to the reverse complement of the k - 1 suffix,
    // and the number of matching characters
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_next_cache_;

    size_t offset_;
    bool k_odd_;
    bool has_sentinel_;

    std::array<size_t, 256> alphabet_encoder_;

    // reset all caches
    void flush();

    // find all parent nodes of node in the CanonicalDBG which are represented
    // in the reverse complement orientation in the underlying primary graph
    void append_prev_rc_nodes(node_index node, std::vector<node_index> &parents) const;

    // find all child nodes of node in the CanonicalDBG which are represented
    // in the reverse complement orientation in the underlying primary graph
    void append_next_rc_nodes(node_index node, std::vector<node_index> &children) const;

    edge_index get_rev_comp_boss_next_node(node_index node) const;
    edge_index get_rev_comp_boss_prev_node(node_index node) const;

    void call_outgoing_from_rev_comp(node_index node,
                                     const std::function<void(node_index, TAlphabet)> &callback) const;

    void call_incoming_to_rev_comp(node_index node,
                                   const std::function<void(node_index, TAlphabet)> &callback) const;
};

} // namespace graph
} // namespace mtg

#endif // __CANONICAL_DBG_HPP__
