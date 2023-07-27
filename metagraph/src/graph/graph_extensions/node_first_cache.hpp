#ifndef __NODE_FIRST_CACHE_HPP__
#define __NODE_FIRST_CACHE_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

// This cache stores intermediate results from BOSS bwd calls to speed up calls
// to call_incoming_kmers in DBGSuccinct. In addition, it caches the indices
// corresponding to the reverse complements of the k-1 prefixes and suffixes of
// each traversed node.
class NodeFirstCache : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename SequenceGraph::node_index;
    using IncomingEdgeCallback = DeBruijnGraph::IncomingEdgeCallback;
    using edge_index = boss::BOSS::edge_index;

    NodeFirstCache(const DBGSuccinct &graph, size_t cache_size = 100'000)
          : dbg_succ_(graph), cache_size_(cache_size), first_cache_(cache_size),
            prefix_rc_cache_(cache_size), suffix_rc_cache_(cache_size) {}

    // Returns the first character of the node's sequence
    char get_first_char(edge_index node, edge_index child_hint = 0) const;

    void call_incoming_edges(edge_index node, const boss::BOSS::Call<edge_index> &callback) const;
    void call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const;

    edge_index get_prefix_rc(edge_index node, const std::string &spelling) const;
    edge_index get_suffix_rc(edge_index node, const std::string &spelling) const;

    bool load(const std::string &) { throw std::runtime_error("Not implemented"); }
    void serialize(const std::string &) const { throw std::runtime_error("Not implemented"); }

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

    size_t max_size() const { return cache_size_; }

  private:
    const DBGSuccinct &dbg_succ_;
    size_t cache_size_;

    // Maps a BOSS edge e to the pair (bwd(e), bwd^(k-1)(e)), where k is the node
    // size in a BOSS graph.
    // Thus, first_cache_[e] == boss.get_minus_k_value(e, boss.get_k() - 1).first
    mutable caches::fixed_sized_cache<edge_index, std::pair<edge_index, edge_index>,
                                      caches::LRUCachePolicy<edge_index>> first_cache_;

    // Fetch or compute (bwd(edge), bwd^(k-1)(edge)), then cache the value.
    // If child_hint != 0, then check first_cache_ if child_hint_'s corresponding
    // pair is cached and compute the pair for edge using that result.
    // FYI: can only be called when the caches are initialized
    std::pair<edge_index, edge_index>
    get_parent_pair(edge_index edge, edge_index child_hint = 0) const;

    mutable caches::fixed_sized_cache<edge_index, edge_index,
                                      caches::LRUCachePolicy<edge_index>> prefix_rc_cache_;
    mutable caches::fixed_sized_cache<edge_index, edge_index,
                                      caches::LRUCachePolicy<edge_index>> suffix_rc_cache_;
};

} // namespace graph
} // namespace mtg

#endif // __NODE_FIRST_CACHE_HPP__
