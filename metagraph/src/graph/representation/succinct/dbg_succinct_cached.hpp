#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"
#include "kmer/kmer_extractor.hpp"
#include "common/caches.hpp"

namespace mtg {
namespace graph {

// A wrapper which caches computed node sequences for DBGSuccinct graphs.
// This allows for faster get_node_sequence and call_incoming_kmers calls.
class DBGSuccinctCached;

std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(std::shared_ptr<const DBGSuccinct> graph, size_t cache_size = 1024);

class DBGSuccinctCached : public DBGWrapper<DBGSuccinct> {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    virtual ~DBGSuccinctCached() {}

    // if the sequence of a node has been constructed externally, cache the result
    virtual void put_decoded_node(node_index node, std::string_view seq) const = 0;

    // get the encoding of the first character of this node's sequence
    virtual TAlphabet get_first_value(edge_index i) const = 0;

    edge_index get_rev_comp_boss_next_node(node_index node) const;
    edge_index get_rev_comp_boss_prev_node(node_index node) const;

    virtual bool operator==(const DeBruijnGraph &other) const override final {
        if (get_k() != other.get_k()
                || num_nodes() != other.num_nodes()
                || get_mode() != other.get_mode())
            return false;

        const auto *other_succ = dynamic_cast<const DBGSuccinct*>(&other);
        if (!other_succ) {
            if (const auto *other_cached = dynamic_cast<const DBGSuccinctCached*>(&other))
                other_succ = &other_cached->get_graph();
        }

        if (other_succ)
            return graph_->operator==(*other_succ);

        throw std::runtime_error("Not implemented");
        return false;
    }

    /**
     * Delegated methods
     */
    #define DELEGATE_METHOD(RETURN_TYPE, METHOD, ARG_TYPE, ARG_NAME) \
    virtual RETURN_TYPE METHOD(ARG_TYPE ARG_NAME) const override final { \
        return graph_->METHOD(ARG_NAME); \
    } \

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
                                = [](){ return false; }) const override final {
        graph_->call_nodes(callback, stop_early);
    }

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override final {
        graph_->call_sequences(callback, num_threads, kmers_in_single_form);
    }

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final {
        graph_->call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override final {
        return graph_->find(sequence, discovery_fraction);
    }

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override final {
        graph_->map_to_nodes(sequence, callback, terminate);
    }

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_outgoing_nodes(node, callback);
    }

    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_incoming_nodes(node, callback);
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        // TODO: implemented a cached version of this later if needed
        return graph_->traverse_back(node, prev_char);
    }

  protected:
    const boss::BOSS *boss_;
    size_t cache_size_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_prev_cache_;
    mutable common::LRUCache<node_index, std::pair<edge_index, size_t>> rev_comp_next_cache_;

    template <typename Graph>
    explicit DBGSuccinctCached(Graph&& graph, size_t cache_size = 1024)
          : DBGWrapper(std::forward<Graph>(graph)), boss_(&graph_->get_boss()),
            cache_size_(cache_size), rev_comp_prev_cache_(cache_size_),
            rev_comp_next_cache_(cache_size_) {}

    virtual TAlphabet complement(TAlphabet c) const = 0;
    virtual std::string decode(const std::vector<TAlphabet> &v) const = 0;
};


template <typename KmerType>
class DBGSuccinctCachedImpl : public DBGSuccinctCached {
    typedef kmer::KmerExtractorBOSS KmerExtractor;

    // KmerType is the encoded k-mer
    // edge_index is the boss node whose last character is the first character of the k-mer
    typedef std::pair<KmerType, std::optional<edge_index>> CacheValue;

  public:
    template <typename Graph>
    DBGSuccinctCachedImpl(Graph&& graph, size_t cache_size)
          : DBGSuccinctCached(std::forward<Graph>(graph), cache_size), decoded_cache_(cache_size_) {}

    virtual void put_decoded_node(node_index node, std::string_view seq) const override final {
        assert(node > 0 && node <= num_nodes());
        assert(seq.size() == graph_->get_k());
        assert(graph_->get_node_sequence(node) == seq);

        put_kmer(graph_->kmer_to_boss_index(node), CacheValue{ to_kmer(seq), std::nullopt });
    }

    virtual std::string get_node_sequence(node_index node) const override final {
        assert(node > 0 && node <= num_nodes());

        // get the sequence from either the cache, or the underlying graph
        auto ret_val = KmerExtractor::kmer_to_sequence<KmerType>(
            get_kmer_pair(graph_->kmer_to_boss_index(node)).first, get_k()
        );

        assert(ret_val == graph_->get_node_sequence(node));
        return ret_val;
    }

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

    virtual void call_outgoing_kmers(node_index node,
                                     const OutgoingEdgeCallback &callback) const override final {
        assert(node > 0 && node <= num_nodes());

        edge_index i = graph_->kmer_to_boss_index(node);
        auto kmer = decoded_cache_.TryGet(i);
        graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                if (kmer)
                    update_node_next(*kmer, graph_->kmer_to_boss_index(next), c);

                callback(next, c);
            }
        });
    }

    virtual node_index traverse(node_index node, char next_char) const override final {
        assert(node > 0 && node <= num_nodes());

        if (node_index next = graph_->traverse(node, next_char)) {
            edge_index i = graph_->kmer_to_boss_index(node);
            if (auto kmer = decoded_cache_.TryGet(i))
                update_node_next(*kmer, graph_->kmer_to_boss_index(next), next_char);

            return next;

        } else {
            return npos;
        }
    }

    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate
                              = [](){ return false; }) const override final {
        graph_->traverse(start, begin, end, callback, terminate);
    }

    virtual void call_incoming_kmers(node_index node,
                                     const IncomingEdgeCallback &callback) const override final {
        assert(node > 0 && node <= num_nodes());

        edge_index edge = graph_->kmer_to_boss_index(node);

        boss_->call_incoming_to_target(boss_->bwd(edge), boss_->get_node_last_value(edge),
            [&](edge_index incoming_boss_edge) {
                assert(boss_->get_W(incoming_boss_edge) % boss_->alph_size
                        == boss_->get_node_last_value(edge));

                auto prev = graph_->boss_to_kmer_index(incoming_boss_edge);

                // get the first character from either the cache, or the graph
                TAlphabet s = get_first_value(incoming_boss_edge);

                if (prev != npos && s != boss::BOSS::kSentinelCode)
                    callback(prev, boss_->decode(s));
            }
        );
    }

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override final {
        size_t k = get_k();
        const char *begin = sequence.data();
        const char *last = begin + k - 1;
        std::optional<KmerType> prev{std::nullopt};
        graph_->map_to_nodes_sequentially(sequence,
            [&](node_index i) {
                if (i) {
                    if (prev) {
                        // update the previous k-mer
                        prev->to_next(k, encode(*last), encode(*(last - 1)));
                    } else {
                        // initialize a new k-mer
                        prev = to_kmer(std::string_view(begin, k));
                    }

                    // cache the result
                    put_kmer(graph_->kmer_to_boss_index(i),
                             CacheValue{ *prev, std::nullopt });
                } else {
                    // reset the k-mer to null
                    prev = std::nullopt;
                }

                callback(i);
                ++begin;
                ++last;
            },
            terminate
        );
    }

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
        auto kmer = decoded_cache_.TryGet(i);

        // if the result of bwd^{k - 1}(i) is desired and not cached, then
        // clear the k-mer so it can be recomputed below
        if (check_second && kmer && !kmer->second)
            kmer = std::nullopt;

        if (!kmer) {
            // get the node sequence and possibly the last node accessed while computing
            auto [seq, last_node_opt] = boss_->get_node_seq_with_end_node(i, check_second);

            // append the last character and cache the result
            seq.push_back(boss_->get_W(i) % boss_->alph_size);

            // cache the result
            kmer = CacheValue{ to_kmer(seq), last_node_opt };
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


inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(std::shared_ptr<const DBGSuccinct> graph, size_t cache_size) {
    if (graph->get_k() * kmer::KmerExtractorBOSS::bits_per_char <= 64) {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer64>>(graph, cache_size);
    } else if (graph->get_k() * kmer::KmerExtractorBOSS::bits_per_char <= 128) {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer128>>(graph, cache_size);
    } else {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer256>>(graph, cache_size);
    }
}


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
