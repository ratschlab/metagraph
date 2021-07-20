#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include <tsl/hopscotch_map.h>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "kmer/kmer_extractor.hpp"

namespace mtg {
namespace graph {


// A wrapper which caches the results of get_node_sequence, get_minus_k_value,
// and call_*_kmers calls for DBGSuccinct graphs.
class DBGSuccinctCached;

// A fast LRU cache
template <typename Key, typename Value>
class ThreadUnsafeLRUCache {
  public:
    ThreadUnsafeLRUCache(size_t max_size) : max_size_(max_size) {}

    std::optional<Value> TryGet(const Key &key);
    void Put(const Key &key, Value value);
    void Clear() {
        values_.clear();
        finder_.clear();
    }

  private:
    size_t max_size_;
    std::list<std::pair<Key, Value>> values_;
    tsl::hopscotch_map<Key, typename std::list<std::pair<Key, Value>>::iterator> finder_;
};

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

    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_ptr_->get_base_graph();
    }

    const DBGSuccinct& get_dbg_succ() const { return *graph_ptr_; }

    virtual void put_decoded_node(node_index node, std::string_view seq) const = 0;
    virtual TAlphabet get_first_value(node_index node) const = 0;

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

    virtual bool load(const std::string &filename) override {
        auto loaded = std::make_shared<DBGSuccinct>(graph_ptr_->get_k());
        if (!loaded->load(filename))
            return false;

        graph_ptr_ = loaded;
        boss_ = &graph_ptr_->get_boss();
        return true;
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


// A fast LRU cache which memoizes computation results. More precisely, if a key
// is already present, it assumes that the corresponding value will be the same.
template <typename KmerType>
class DBGSuccinctCachedImpl : public DBGSuccinctCached {
    typedef kmer::KmerExtractorBOSS KmerExtractor;

  public:
    template <typename... Args>
    DBGSuccinctCachedImpl(Args&&... args)
          : DBGSuccinctCached(std::forward<Args>(args)...),
            decoded_cache_(cache_size_) {}

    void put_decoded_node(node_index node, std::string_view seq) const {
        assert(seq.size() == graph_ptr_->get_k());
        put_kmer(graph_ptr_->kmer_to_boss_index(node), std::make_pair(to_kmer(seq), 0));
    }

    std::string get_node_sequence(node_index node) const {
        auto ret_val = KmerExtractor::kmer_to_sequence<KmerType>(
            get_kmer_pair(graph_ptr_->kmer_to_boss_index(node)).first,
            graph_ptr_->get_k()
        );

        assert(ret_val == graph_ptr_->get_node_sequence(node));
        return ret_val;
    }

    TAlphabet get_first_value(node_index node) const {
        edge_index i = graph_ptr_->kmer_to_boss_index(node);
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
                    update_node_next_from_kmer(*kmer, next, c);

                callback(next, c);
            }
        });
    }

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const {
        assert(node > 0 && node <= num_nodes());

        edge_index i = graph_ptr_->kmer_to_boss_index(node);
        auto kmer = decoded_cache_.TryGet(i);
        edge_index edge = graph_ptr_->kmer_to_boss_index(node);

        boss_->call_incoming_to_target(boss_->bwd(edge), boss_->get_node_last_value(edge),
            [&](edge_index incoming_boss_edge) {
                assert(boss_->get_W(incoming_boss_edge) % boss_->alph_size
                        == boss_->get_node_last_value(edge));

                auto prev = graph_ptr_->boss_to_kmer_index(incoming_boss_edge);
                TAlphabet s = get_first_value(prev);
                if (prev != npos && s != boss::BOSS::kSentinelCode) {
                    char c = boss_->decode(s);

                    if (kmer)
                        update_node_prev_from_kmer(*kmer, prev, c);

                    callback(prev, c);
                }
            }
        );
    }

    node_index traverse(node_index node, char next_char) const {
        if (node_index next = graph_ptr_->traverse(node, next_char)) {
            edge_index i = graph_ptr_->kmer_to_boss_index(node);
            if (auto kmer = decoded_cache_.TryGet(i))
                update_node_next_from_kmer(*kmer, next, next_char);

            return next;

        } else {
            return npos;
        }
    }

    node_index traverse_back(node_index node, char prev_char) const {
        if (node_index prev = graph_ptr_->traverse_back(node, prev_char)) {
            edge_index i = graph_ptr_->kmer_to_boss_index(node);
            if (auto kmer = decoded_cache_.TryGet(i))
                update_node_prev_from_kmer(*kmer, prev, prev_char);

            return prev;

        } else {
            return npos;
        }
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
        if (!DBGSuccinctCached::load(filename))
            return false;

        decoded_cache_.Clear();
        return true;
    }

  private:
    typedef std::pair<KmerType, edge_index> CacheNode;

    // TODO: switch to a more efficient LRU cache implementation?
    mutable ThreadUnsafeLRUCache<edge_index, CacheNode> decoded_cache_;

    inline static constexpr TAlphabet encode(char c) { return KmerExtractor::encode(c); }
    inline static constexpr std::string decode(const std::vector<TAlphabet> &v) {
        return KmerExtractor::decode(v);
    }
    inline static constexpr KmerType to_kmer(std::string_view seq) {
        return KmerExtractor::sequence_to_kmer<KmerType>(seq);
    }
    inline static constexpr KmerType to_kmer(const std::vector<TAlphabet> &seq) {
        // TODO: avoid decode step
        return KmerExtractor::sequence_to_kmer<KmerType>(decode(seq));
    }

    void put_kmer(edge_index key, CacheNode value) const {
        assert(!value.second || value.first[1] == boss_->get_node_last_value(value.second));
        assert(value.first[1] == boss_->get_minus_k_value(key, get_k() - 2).first);
        assert(to_kmer(graph_ptr_->get_node_sequence(graph_ptr_->boss_to_kmer_index(key)))
            == value.first);
        decoded_cache_.Put(key, std::move(value));
    }

    CacheNode get_kmer_pair(edge_index i, bool check_second = false) const {
        auto kmer = decoded_cache_.TryGet(i);
        if (!kmer || (check_second && !kmer->second)) {
            auto [seq, first] = get_boss_node_seq(i);

            kmer = !kmer ? std::make_pair(to_kmer(seq), first)
                         : std::make_pair(kmer->first, first);

            put_kmer(i, *kmer);
        }

        return *kmer;
    }

    void update_node_next_from_kmer(CacheNode kmer, node_index next, char c) const {
        kmer.first.to_next(graph_ptr_->get_k(), encode(c));
        kmer.second = 0;
        put_kmer(graph_ptr_->kmer_to_boss_index(next), std::move(kmer));
    }

    void update_node_prev_from_kmer(CacheNode kmer, node_index prev, char c) const {
        kmer.first.to_prev(graph_ptr_->get_k(), encode(c));
        kmer.second = 0;
        put_kmer(graph_ptr_->kmer_to_boss_index(prev), std::move(kmer));
    }

    std::pair<std::vector<TAlphabet>, edge_index> get_boss_node_seq(edge_index i) const {
        auto ret_val = boss_->get_node_seq_with_node(i);
        assert(ret_val.first[0] == boss_->get_node_last_value(ret_val.second));
        ret_val.first.push_back(boss_->get_W(i) % boss_->alph_size);
        return ret_val;
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


template <typename Key, typename Value>
inline std::optional<Value> ThreadUnsafeLRUCache<Key, Value>::TryGet(const Key &key) {
    auto find = finder_.find(key);
    if (find == finder_.end())
        return std::nullopt;

    values_.splice(values_.begin(), values_, find->second);
    return find->second->second;
}

template <typename Key, typename Value>
inline void ThreadUnsafeLRUCache<Key, Value>::Put(const Key &key, Value value) {
    auto find = finder_.find(key);
    if (find != finder_.end()) {
        find->second->second = std::move(value);
        values_.splice(values_.begin(), values_, find->second);
    } else {
        values_.insert(values_.begin(), std::make_pair(key, std::move(value)));
        finder_[key] = values_.begin();
    }

    if (values_.size() > max_size_) {
        finder_.erase(values_.back().first);
        values_.pop_back();
    }
}

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
