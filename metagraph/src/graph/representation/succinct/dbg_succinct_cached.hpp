#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include "dbg_succinct.hpp"
#include "kmer/kmer_extractor.hpp"

namespace mtg {
namespace graph {


template <typename KmerType>
class DBGSuccinctCachedViewImpl : public DBGSuccinct::CachedView {
  public:
    template <typename Graph>
    DBGSuccinctCachedViewImpl(Graph&& graph, size_t cache_size)
          : DBGSuccinct::CachedView(std::forward<Graph>(graph)),
            boss_(&graph_->get_boss()),
            decoded_cache_(cache_size) {}

    virtual TAlphabet get_first_value(edge_index i) const override final {
        assert(i);
        return std::get<0>(get_kmer_tuple(i))[1];
    }

    /**
     * Methods from DBGSuccinct
     */
    static constexpr auto ALWAYS_FALSE = []() { return false; };

    virtual void
    map_to_nodes_sequentially(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = ALWAYS_FALSE,
                              const std::function<bool()> &skip = ALWAYS_FALSE) const override final;

    /**
     * Methods from DeBruijnGraph
     */
    virtual std::string get_node_sequence(node_index node) const override final;

    virtual void call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const override final;

    virtual void call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const override final;

    virtual bool load(const std::string &filename_base) override final {
        if (const_cast<DBGSuccinct*>(graph_.get())->load(filename_base)) {
            decoded_cache_.Clear();
            return true;
        }

        return false;
    }

    virtual void set_graph(std::shared_ptr<const DBGSuccinct> graph) override final {
        graph_ = graph;
        decoded_cache_.Clear();
    }

    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion) override final {
        const_cast<DBGSuccinct*>(graph_.get())->add_sequence(sequence, on_insertion);
        decoded_cache_.Clear();
    }

  private:
    const boss::BOSS *boss_;

    // KmerType is the encoded k-mer
    // edge_index is the boss node whose last character is the first character of the k-mer
    typedef std::tuple<KmerType, edge_index, edge_index> CacheValue;

    mutable common::LRUCache<edge_index, CacheValue> decoded_cache_;

    // cache a computed result
    void put_kmer(edge_index key, CacheValue value) const;

    CacheValue get_kmer_tuple(edge_index i, edge_index child_hint = 0, bool extract_first = false) const;

    inline static constexpr TAlphabet encode(char c) { return kmer::KmerExtractorBOSS::encode(c); }
};


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
