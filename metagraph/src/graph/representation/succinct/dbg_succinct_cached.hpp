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
    typedef std::pair<KmerType, std::optional<edge_index>> CacheValue;

    mutable common::LRUCache<edge_index, CacheValue> decoded_cache_;

    // cache a computed result
    inline void put_kmer(edge_index key, CacheValue value) const {
        assert(key);
        assert(!value.second || value.first[1] == boss_->get_node_last_value(*value.second));
        assert(value.first[1] == boss_->get_minus_k_value(key, get_k() - 2).first);
        assert(!value.second || *value.second
            == boss_->get_minus_k_value(key, get_k() - 2).second);

        if (value.second) {
            if (std::optional<CacheValue> kmer = decoded_cache_.TryGet(key)) {
                if (kmer->second) {
                    assert(kmer->second == value.second);
                    return;
                }
            }
        }

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
            kmer = CacheValue{ KmerType(seq), std::move(last_node) };
            put_kmer(i, *kmer);
        }

        return *kmer;
    }

    inline static constexpr TAlphabet encode(char c) { return kmer::KmerExtractorBOSS::encode(c); }
};


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
