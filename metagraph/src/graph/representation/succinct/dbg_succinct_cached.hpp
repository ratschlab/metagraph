#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "kmer/kmer_extractor.hpp"
#include "common/caches.hpp"

namespace mtg {
namespace graph {


// A wrapper which caches computed node sequences for DBGSuccinct graphs.
// This allows for faster get_node_sequence and call_incoming_kmers calls.
class DBGSuccinctCached;

std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size = 1024);

inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(DeBruijnGraph &graph, size_t cache_size = 1024) {
    return make_cached_dbgsuccinct({ std::shared_ptr<DeBruijnGraph>{}, &graph }, cache_size);
}


class DBGSuccinctCached : public DeBruijnGraph {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    explicit DBGSuccinctCached(std::shared_ptr<DBGSuccinct> graph, size_t cache_size)
          : graph_ptr_(graph), boss_(&graph_ptr_->get_boss()), cache_size_(cache_size) {
        assert(graph_ptr_ && "Only DBGSuccinct can be cached");
    }

    explicit DBGSuccinctCached(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size)
          : DBGSuccinctCached(std::dynamic_pointer_cast<DBGSuccinct>(graph), cache_size) {}

    explicit DBGSuccinctCached(DBGSuccinct &dbg, size_t cache_size)
          : DBGSuccinctCached(std::shared_ptr<DBGSuccinct>(
                std::shared_ptr<DBGSuccinct>{}, &dbg
            ), cache_size) {}

    explicit DBGSuccinctCached(DeBruijnGraph &graph, size_t cache_size)
          : DBGSuccinctCached(dynamic_cast<DBGSuccinct&>(graph), cache_size) {}

    virtual ~DBGSuccinctCached() {}

    // if the sequence of a node has been constructed externally, cache the result
    virtual void put_decoded_node(node_index node, std::string_view seq) const = 0;

    // get the encoding of the first character of this node's sequence
    virtual TAlphabet get_first_value(edge_index i) const = 0;

    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_ptr_->get_base_graph();
    }

    const DBGSuccinct& get_dbg_succ() const { return *graph_ptr_; }

    virtual size_t get_k() const override final { return graph_ptr_->get_k(); }

    virtual Mode get_mode() const override final { return graph_ptr_->get_mode(); }

    virtual size_t outdegree(node_index node) const override final {
        return graph_ptr_->outdegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override final {
        return graph_ptr_->has_single_outgoing(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override final {
        return graph_ptr_->has_multiple_outgoing(node);
    }

    virtual size_t indegree(node_index node) const override final {
        return graph_ptr_->indegree(node);
    }

    virtual bool has_no_incoming(node_index node) const override final {
        return graph_ptr_->has_no_incoming(node);
    }

    virtual bool has_single_incoming(node_index node) const override final {
        return graph_ptr_->has_single_incoming(node);
    }

    virtual const std::string& alphabet() const override final { return graph_ptr_->alphabet(); }

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override final {
        graph_ptr_->map_to_nodes(sequence, callback, terminate);
    }

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_ptr_->adjacent_outgoing_nodes(node, callback);
    }
    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_ptr_->adjacent_incoming_nodes(node, callback);
    }

    virtual uint64_t num_nodes() const override final { return graph_ptr_->num_nodes(); }

    virtual void serialize(const std::string &filename) const override final {
        graph_ptr_->serialize(filename);
    }

    virtual std::string file_extension() const override final {
        return graph_ptr_->file_extension();
    }

    virtual bool operator==(const DeBruijnGraph &other) const override final {
        if (get_k() != other.get_k()
                || num_nodes() != other.num_nodes()
                || get_mode() != other.get_mode())
            return false;

        const auto *other_succ = dynamic_cast<const DBGSuccinct*>(&other);
        if (!other_succ) {
            if (const auto *other_cached = dynamic_cast<const DBGSuccinctCached*>(&other))
                other_succ = &other_cached->get_dbg_succ();
        }

        if (other_succ)
            return graph_ptr_->operator==(*other_succ);

        throw std::runtime_error("Not implemented");
        return false;
    }

  protected:
    std::shared_ptr<DBGSuccinct> graph_ptr_;
    const boss::BOSS *boss_;
    size_t cache_size_;
};


template <typename KmerType>
class DBGSuccinctCachedImpl : public DBGSuccinctCached {
    typedef kmer::KmerExtractorBOSS KmerExtractor;

  public:
    template <typename... Args>
    DBGSuccinctCachedImpl(Args&&... args)
          : DBGSuccinctCached(std::forward<Args>(args)...),
            decoded_cache_(cache_size_) {}

    void put_decoded_node(node_index node, std::string_view seq) const {
        assert(node > 0 && node <= num_nodes());
        assert(seq.size() == graph_ptr_->get_k());
        assert(graph_ptr_->get_node_sequence(node) == seq);

        put_kmer(graph_ptr_->kmer_to_boss_index(node), std::make_pair(to_kmer(seq), 0));
    }

    std::string get_node_sequence(node_index node) const {
        assert(node > 0 && node <= num_nodes());

        // get the sequence from either the cache, or the underlying graph
        auto ret_val = KmerExtractor::kmer_to_sequence<KmerType>(
            get_kmer_pair(graph_ptr_->kmer_to_boss_index(node)).first,
            graph_ptr_->get_k()
        );

        assert(ret_val == graph_ptr_->get_node_sequence(node));
        return ret_val;
    }

    TAlphabet get_first_value(edge_index i) const {
        assert(i);

        // Take k - 1 traversal steps to construct the node sequence. This way,
        // the last node visited and its parent can be cached for subsequence calls.
        auto kmer_pair = get_kmer_pair(i, true);
        assert(kmer_pair.second);

        TAlphabet c = kmer_pair.first[1];

        edge_index bwd_i = boss_->bwd(i);
        auto fetch_bwd = decoded_cache_.TryGet(bwd_i);
        if (!fetch_bwd) {
            kmer_pair.second = boss_->bwd(kmer_pair.second);
            kmer_pair.first.to_prev(graph_ptr_->get_k(),
                                    boss_->get_node_last_value(kmer_pair.second));
            put_kmer(bwd_i, std::move(kmer_pair));
        }

        return c;
    }

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const {
        assert(node > 0 && node <= num_nodes());

        edge_index i = graph_ptr_->kmer_to_boss_index(node);
        auto kmer = decoded_cache_.TryGet(i);
        graph_ptr_->call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                if (kmer)
                    update_node_next(*kmer, graph_ptr_->kmer_to_boss_index(next), c);

                callback(next, c);
            }
        });
    }

    node_index traverse(node_index node, char next_char) const {
        assert(node > 0 && node <= num_nodes());

        if (node_index next = graph_ptr_->traverse(node, next_char)) {
            edge_index i = graph_ptr_->kmer_to_boss_index(node);
            if (auto kmer = decoded_cache_.TryGet(i))
                update_node_next(*kmer, graph_ptr_->kmer_to_boss_index(next), next_char);

            return next;

        } else {
            return npos;
        }
    }

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const {
        assert(node > 0 && node <= num_nodes());

        edge_index edge = graph_ptr_->kmer_to_boss_index(node);

        boss_->call_incoming_to_target(boss_->bwd(edge), boss_->get_node_last_value(edge),
            [&](edge_index incoming_boss_edge) {
                assert(boss_->get_W(incoming_boss_edge) % boss_->alph_size
                        == boss_->get_node_last_value(edge));

                auto prev = graph_ptr_->boss_to_kmer_index(incoming_boss_edge);

                // get the first character from either the cache, or the graph
                TAlphabet s = get_first_value(incoming_boss_edge);

                if (prev != npos && s != boss::BOSS::kSentinelCode)
                    callback(prev, boss_->decode(s));
            }
        );
    }

    node_index traverse_back(node_index node, char prev_char) const {
        // TODO: handle caching at some point. For now, this method isn't really used
        return graph_ptr_->traverse_back(node, prev_char);
    }

    void map_to_nodes_sequentially(std::string_view sequence,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate
                                       = [](){ return false; }) const {
        size_t k = graph_ptr_->get_k();
        const char *begin = sequence.data();
        const char *last = begin + k - 1;
        std::optional<KmerType> prev{std::nullopt};
        graph_ptr_->map_to_nodes_sequentially(sequence,
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
                    put_kmer(graph_ptr_->kmer_to_boss_index(i), std::make_pair(*prev, 0));
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

    void add_sequence(std::string_view sequence,
                      const std::function<void(node_index)> &on_insertion
                          = [](uint64_t) {}) {
        graph_ptr_->add_sequence(sequence, on_insertion);
        decoded_cache_.Clear();
    }

    bool load(const std::string &filename) {
        auto loaded = std::make_shared<DBGSuccinct>(graph_ptr_->get_k());
        if (!loaded->load(filename))
            return false;

        graph_ptr_ = loaded;
        boss_ = &graph_ptr_->get_boss();
        decoded_cache_.Clear();
        return true;
    }

  private:
    // KmerType is the encoded k-mer
    // edge_index is the boss node whose last character is the first character of the k-mer
    typedef std::pair<KmerType, edge_index> CacheValue;

    mutable common::ThreadUnsafeLRUCache<edge_index, CacheValue> decoded_cache_;

    // cache a computed result
    void put_kmer(edge_index key, CacheValue value) const {
        assert(key);
        assert(!value.second || value.first[1] == boss_->get_node_last_value(value.second));
        assert(value.first[1] == boss_->get_minus_k_value(key, get_k() - 2).first);
        assert(!value.second || value.second
            == boss_->get_minus_k_value(key, get_k() - 2).second);
        decoded_cache_.Put(key, std::move(value));
    }

    CacheValue get_kmer_pair(edge_index i, bool check_second = false) const {
        assert(i);
        auto kmer = decoded_cache_.TryGet(i);

        // if the result of bwd^{k - 1}(i) is desired and not cached, then
        // clear the k-mer so it can be recomputed below
        if (check_second && !kmer->second)
            kmer = std::nullopt;

        if (!kmer) {
            std::vector<TAlphabet> seq;

            // initialize
            kmer = CacheValue{ KmerType(), 0 };

            if (check_second) {
                // we need the result of bwd^{k - 1}(i), so take all traversal steps
                std::tie(seq, kmer->second) = boss_->get_node_seq_with_end_node(i);
            } else {
                // We don't need that value, so use the built-in method to compute
                // the sequence. This uses the BOSS suffix index to speed it up.
                seq = boss_->get_node_seq(i);
            }

            // append the last character and cache the result
            seq.push_back(boss_->get_W(i) % boss_->alph_size);
            kmer->first = to_kmer(seq);

            put_kmer(i, *kmer);
        }

        return *kmer;
    }

    void update_node_next(CacheValue kmer, edge_index next, char c) const {
        assert(c != boss::BOSS::kSentinel);

        // update the k-mer with the next character
        kmer.first.to_next(graph_ptr_->get_k(), encode(c));
        kmer.second = 0;
        put_kmer(next, std::move(kmer));
    }

    inline static std::string decode(const std::vector<TAlphabet> &v) {
        return KmerExtractor::decode(v);
    }

    inline static constexpr TAlphabet encode(char c) { return KmerExtractor::encode(c); }
    inline static constexpr KmerType to_kmer(std::string_view seq) {
        return KmerExtractor::sequence_to_kmer<KmerType>(seq);
    }
    inline static constexpr KmerType to_kmer(const std::vector<TAlphabet> &seq) {
        // TODO: avoid decode step
        return to_kmer(decode(seq));
    }
};


inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size) {
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
