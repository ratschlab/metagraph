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

// A wrapper which caches computed node sequences for DBGSuccinct graphs.
// This allows for faster get_node_sequence and call_incoming_kmers calls.
class DBGSuccinctCachedView : public DBGWrapper<DBGSuccinct> {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    virtual ~DBGSuccinctCachedView() {}

    // if the sequence of a node has been constructed externally, cache the result
    virtual void put_decoded_node(node_index node, std::string_view seq) const = 0;

    // get the encoding of the first character of this node's sequence
    virtual TAlphabet get_first_value(edge_index i) const = 0;

    edge_index get_rev_comp_boss_next_node(node_index node) const;
    edge_index get_rev_comp_boss_prev_node(node_index node) const;

    virtual bool operator==(const DeBruijnGraph &other) const override final;

    /**
     * Delegated methods
     */
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

    virtual node_index traverse_back(node_index node, char prev_char) const override final;

  protected:
    const boss::BOSS *boss_;
    size_t cache_size_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_prev_cache_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_next_cache_;

    template <typename Graph>
    explicit DBGSuccinctCachedView(Graph&& graph, size_t cache_size = 1024);

    virtual TAlphabet complement(TAlphabet c) const = 0;
    virtual std::string decode(const std::vector<TAlphabet> &v) const = 0;
};


template <typename KmerType>
class DBGSuccinctCachedViewImpl : public DBGSuccinctCachedView {
    typedef kmer::KmerExtractorBOSS KmerExtractor;

    // KmerType is the encoded k-mer
    // edge_index is the boss node whose last character is the first character of the k-mer
    typedef std::pair<KmerType, std::optional<edge_index>> CacheValue;

  public:
    template <typename Graph>
    DBGSuccinctCachedViewImpl(Graph&& graph, size_t cache_size);

    virtual void put_decoded_node(node_index node, std::string_view seq) const override final;

    virtual std::string get_node_sequence(node_index node) const override final;

    virtual TAlphabet get_first_value(edge_index i) const override final;

    virtual void call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const override final;

    virtual node_index traverse(node_index node, char next_char) const override final;

    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate
                              = [](){ return false; }) const override final;

    virtual void call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const override final;

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override final;

  private:
    mutable common::LRUCache<edge_index, CacheValue> decoded_cache_;

    // cache a computed result
    void put_kmer(edge_index key, CacheValue value) const {
        assert(key);
        assert(!value.second || value.first[1] == boss_->get_node_last_value(*value.second));
        assert(value.first[1] == boss_->get_minus_k_value(key, get_k() - 2).first);
        assert(!value.second || *value.second
            == boss_->get_minus_k_value(key, get_k() - 2).second);
        decoded_cache_.Put(key, std::move(value));
    }

    CacheValue get_kmer_pair(edge_index i, bool check_second = false) const {
        assert(i);
        std::optional<CacheValue> kmer = decoded_cache_.TryGet(i);

        // if the result of bwd^{k - 1}(i) is desired and not cached, then
        // clear the k-mer so it can be recomputed below
        if (check_second && kmer && !kmer->second)
            kmer = std::nullopt;

        if (!kmer) {
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

    void update_node_next(CacheValue kmer, edge_index next, char c) const {
        assert(c != boss::BOSS::kSentinel);

        // update the k-mer with the next character
        kmer.first.to_next(get_k(), encode(c));
        kmer.second = std::nullopt;
        put_kmer(next, std::move(kmer));
    }

    virtual std::string decode(const std::vector<TAlphabet> &v) const override final {
        return KmerExtractor::decode(v);
    }

    inline static constexpr TAlphabet encode(char c) { return KmerExtractor::encode(c); }
    inline static constexpr KmerType to_kmer(std::string_view seq) {
        return KmerExtractor::sequence_to_kmer<KmerType>(seq);
    }
    inline static constexpr KmerType to_kmer(const std::vector<TAlphabet> &seq) {
        // TODO: avoid decode step
        return to_kmer(KmerExtractor::decode(seq));
    }

    virtual TAlphabet complement(TAlphabet c) const override final {
        return KmerExtractor::complement(c);
    }
};


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
