#ifndef __NODE_FIRST_CACHE_HPP__
#define __NODE_FIRST_CACHE_HPP__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

// This cache stores intermediate results from BOSS bwd calls to speed up calls
// to call_incoming_kmers in DBGSuccinct.
// This stores 80 bytes per cached 8-byte node index.
class NodeFirstCache : public SequenceGraph::GraphExtension {
  public:
    using node_index = typename SequenceGraph::node_index;
    using IncomingEdgeCallback = DeBruijnGraph::IncomingEdgeCallback;

    NodeFirstCache() : dbg_succ_(nullptr), first_cache_(1) {}
    NodeFirstCache(const DBGSuccinct &graph, size_t cache_size = 100'000)
          : dbg_succ_(&graph), first_cache_(cache_size) {}

    // Returns the first character of the node's sequence
    char get_first_char(node_index node) const;

    void call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const;

    bool load(const std::string &) { throw std::runtime_error("Not implemented"); }
    void serialize(const std::string &) const { throw std::runtime_error("Not implemented"); }

    bool is_compatible(const SequenceGraph &graph, bool verbose = true) const;

  private:
    const DBGSuccinct *dbg_succ_;

    using edge_index = boss::BOSS::edge_index;

    // Maps a BOSS edge e to the pair (bwd(e), bwd^(k-1)(e)), where k is the node
    // size in a BOSS graph.
    // Thus, first_cache_[e] == boss.get_minus_k_value(e, boss.get_k() - 1).first
    mutable caches::fixed_sized_cache<edge_index, std::pair<edge_index, edge_index>,
                                      caches::LRUCachePolicy<node_index>> first_cache_;

    // Fetch or compute (bwd(edge), bwd^(k-1)(edge)), then cache the value.
    // If child_hint != 0, then check first_cache_ if child_hint_'s corresponding
    // pair is cached and compute the pair for edge using that result.
    std::pair<edge_index, edge_index>
    get_parent_pair(edge_index edge, edge_index child_hint = 0) const;
};

} // namespace graph
} // namespace mtg

#endif // __NODE_FIRST_CACHE_HPP__
