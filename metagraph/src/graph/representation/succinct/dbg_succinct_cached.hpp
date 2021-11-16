#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"
#include "kmer/kmer_extractor.hpp"
#include "common/caches.hpp"

namespace mtg {
namespace graph {

class DBGSuccinct;

/**
 * A wrapper which caches computed node sequences for DBGSuccinct graphs.
 * This allows for faster get_node_sequence and call_incoming_kmers calls.
 * Traversal operations and map_to_nodes_sequentially are overridden to fill
 * the cache alongside their normal functions.
 * This wrapper stores ~144 bytes per cached node.
 */
class DBGSuccinctCachedView : public DBGWrapper<DBGSuccinct> {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    virtual ~DBGSuccinctCachedView() {}

    // if the sequence of a BOSS edge and edge label has been constructed externally,
    // cache the result
    virtual void put_decoded_edge(edge_index edge, std::string_view seq) const = 0;

    void call_and_cache_outgoing_from_rev_comp(node_index node,
                                               const std::function<void(node_index, TAlphabet)> &callback) const;

    void call_and_cache_incoming_to_rev_comp(node_index node,
                                             const std::function<void(node_index, TAlphabet)> &callback) const;


    /**
     * Methods from DeBruijnGraph
     */
    virtual bool operator==(const DeBruijnGraph &other) const override final;

    #define DELEGATE_METHOD(RETURN_TYPE, METHOD, ARG_TYPE, ARG_NAME) \
    virtual RETURN_TYPE METHOD(ARG_TYPE ARG_NAME) const override final;

    DELEGATE_METHOD(void, call_kmers, const std::function<void(node_index, const std::string&)> &, callback)
    DELEGATE_METHOD(void, call_source_nodes, const std::function<void(node_index)> &, callback)

    DELEGATE_METHOD(uint64_t, num_nodes, , )
    DELEGATE_METHOD(uint64_t, max_index, , )
    DELEGATE_METHOD(size_t, outdegree, node_index, node)
    DELEGATE_METHOD(bool, has_single_outgoing, node_index, node)
    DELEGATE_METHOD(bool, has_multiple_outgoing, node_index, node)
    DELEGATE_METHOD(size_t, indegree, node_index, node)
    DELEGATE_METHOD(bool, has_no_incoming, node_index, node)
    DELEGATE_METHOD(bool, has_single_incoming, node_index, node)
    DELEGATE_METHOD(node_index, kmer_to_node, std::string_view, kmer)

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override final;

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override final;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final;

    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override final;

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override final;

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final;

    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final;

    // TODO: these can be overloaded to cache values, but this functionality is not
    //       needed now.
    virtual node_index traverse(node_index node, char next_char) const override final;
    virtual node_index traverse_back(node_index node, char prev_char) const override final;
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate
                              = [](){ return false; }) const override final;

  protected:
    const boss::BOSS *boss_;
    size_t cache_size_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_prev_cache_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_next_cache_;

    template <typename Graph>
    explicit DBGSuccinctCachedView(Graph&& graph, size_t cache_size);

    virtual TAlphabet complement(TAlphabet c) const = 0;
    virtual std::string decode(const std::vector<TAlphabet> &v) const = 0;

    edge_index get_rev_comp_boss_next_node(node_index node) const;
    edge_index get_rev_comp_boss_prev_node(node_index node) const;

    // get the encoding of the first character of this node's sequence
    virtual TAlphabet get_first_value(edge_index i) const = 0;
};


template <typename KmerType>
class DBGSuccinctCachedViewImpl : public DBGSuccinctCachedView {
  public:
    template <typename Graph>
    DBGSuccinctCachedViewImpl(Graph&& graph, size_t cache_size);

    /**
     * Methods from DeBruijnGraph
     */
    virtual std::string get_node_sequence(node_index node) const override final;

    virtual void call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const override final;

    virtual void call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const override final;

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override final;

    /**
     * Methods from DBGSuccinctCachedView
     */
    virtual void put_decoded_edge(edge_index edge, std::string_view seq) const override final;

  private:
    // KmerType is the encoded k-mer
    // edge_index is the boss node whose last character is the first character of the k-mer
    typedef std::pair<KmerType, std::optional<edge_index>> CacheValue;

    mutable common::LRUCache<edge_index, CacheValue> decoded_cache_;

    virtual TAlphabet get_first_value(edge_index i) const override final {
        assert(i);

        // Take k - 1 traversal steps to construct the node sequence. This way,
        // the last node visited and its parent can be cached for subsequence calls.
        auto kmer_pair = get_kmer_pair(i, true);
        assert(kmer_pair.second);

        TAlphabet c = kmer_pair.first[1];

        edge_index bwd_i = boss_->bwd(i);
        auto fetch_bwd = decoded_cache_.TryGet(bwd_i);
        if (!fetch_bwd) {
            kmer_pair.second = boss_->bwd(*kmer_pair.second);
            kmer_pair.first.to_prev(get_k(), boss_->get_node_last_value(*kmer_pair.second));
            put_kmer(bwd_i, std::move(kmer_pair));
        }

        return c;
    }

    // cache a computed result
    inline void put_kmer(edge_index key, CacheValue value) const {
        assert(key);
        assert(!value.second || value.first[1] == boss_->get_node_last_value(*value.second));
        assert(value.first[1] == boss_->get_minus_k_value(key, get_k() - 2).first);
        assert(!value.second || *value.second
            == boss_->get_minus_k_value(key, get_k() - 2).second);
        decoded_cache_.Put(key, std::move(value));
    }

    inline CacheValue get_kmer_pair(edge_index i, bool check_second = false) const {
        assert(i);
        std::optional<CacheValue> kmer = decoded_cache_.TryGet(i);

        // if the kmer is not cached, or if the result of bwd^{k - 1}(i) is desired
        // and not cached, then recompute both
        if (!kmer || (check_second && !kmer->second)) {
            // get the node sequence and possibly the last node accessed while computing
            auto [seq, last_traversed_node, steps_remaining]
                = boss_->get_node_seq_with_last_traversed_node(i);

            std::optional<edge_index> last_node;

            if (steps_remaining) {
                if (check_second) {
                    last_node = boss_->get_minus_k_value(last_traversed_node,
                                                         steps_remaining).second;
                }
            } else {
                last_node = last_traversed_node;
            }

            // append the last character and cache the result
            seq.push_back(boss_->get_W(i) % boss_->alph_size);

            // cache the result
            kmer = CacheValue{ to_kmer(seq), std::move(last_node) };
            put_kmer(i, *kmer);
        }

        return *kmer;
    }

    inline void update_node_next(CacheValue kmer, edge_index next, char c) const {
        assert(c != boss::BOSS::kSentinel);

        // update the k-mer with the next character
        kmer.first.to_next(get_k(), encode(c));
        kmer.second = std::nullopt;
        put_kmer(next, std::move(kmer));
    }

    virtual std::string decode(const std::vector<TAlphabet> &v) const override final {
        return kmer::KmerExtractorBOSS::decode(v);
    }

    inline static constexpr TAlphabet encode(char c) { return kmer::KmerExtractorBOSS::encode(c); }
    inline static constexpr KmerType to_kmer(std::string_view seq) {
        return kmer::KmerExtractorBOSS::sequence_to_kmer<KmerType>(seq);
    }
    inline static constexpr KmerType to_kmer(const std::vector<TAlphabet> &seq) {
        // TODO: avoid decode step
        return to_kmer(kmer::KmerExtractorBOSS::decode(seq));
    }

    virtual TAlphabet complement(TAlphabet c) const override final {
        return kmer::KmerExtractorBOSS::complement(c);
    }
};


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
