#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "kmer/kmer_extractor.hpp"

namespace mtg {
namespace graph {


class DBGSuccinctCached : public DeBruijnGraph {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    // TODO: switch to a more efficient LRU cache implementation?
    template <typename Key, typename Value>
    using LRUCache = caches::fixed_sized_cache<Key, Value, caches::LRUCachePolicy<Key>>;

    explicit DBGSuccinctCached(std::shared_ptr<DBGSuccinct> graph, size_t cache_size)
          : graph_ptr_(graph), boss_(&graph_ptr_->get_boss()),
            cache_size_(cache_size), bwd_first_cache_(cache_size_) {
        assert(graph_ptr_);
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
    virtual std::optional<std::string> get_decoded_node(node_index node) const = 0;

    inline TAlphabet get_first_value(edge_index i) const {
        TAlphabet c;
        edge_index first;

        if (auto fetch = bwd_first_cache_.TryGet(i)) {
            c = boss_->get_node_last_value(*fetch);
            first = *fetch;
        } else {
            std::tie(c, first) = boss_->get_minus_k_value(i, get_k() - 2);
            bwd_first_cache_.Put(i, first);
        }

        edge_index bwd_i = boss_->bwd(i);
        auto fetch_bwd = bwd_first_cache_.TryGet(bwd_i);
        if (!fetch_bwd)
            bwd_first_cache_.Put(bwd_i, boss_->bwd(first));

        return c;
    }

    virtual std::string get_node_sequence(node_index node) const override final {
        if (auto fetch = get_decoded_node(node)) {
            return *fetch;
        } else {
            std::string seq = graph_ptr_->get_node_sequence(node);
            put_decoded_node(node, seq);
            return seq;
        }
    }

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


    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion
                                  = [](uint64_t) {}) override {
        graph_ptr_->add_sequence(sequence, on_insertion);
        bwd_first_cache_.Clear();
    }

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
        bwd_first_cache_.Clear();
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

  private:
    mutable LRUCache<edge_index, edge_index> bwd_first_cache_;
};

template <typename KmerType>
class DBGSuccinctCachedImpl : public DBGSuccinctCached {
  public:
    typedef kmer::KmerExtractorBOSS KmerExtractor;

    template <typename... Args>
    DBGSuccinctCachedImpl(Args&&... args)
          : DBGSuccinctCached(std::forward<Args>(args)...),
            decoded_cache_(cache_size_ * graph_ptr_->alphabet().size()) {}

    void put_decoded_node(node_index node, std::string_view seq) const {
        assert(seq.size() == graph_ptr_->get_k());
        decoded_cache_.Put(node, to_kmer(seq));
    }

    std::optional<std::string> get_decoded_node(node_index node) const {
        if (auto kmer = decoded_cache_.TryGet(node)) {
            auto seq = KmerExtractor::kmer_to_sequence<KmerType>(*kmer, graph_ptr_->get_k());
            assert(seq.size() == graph_ptr_->get_k());
            return seq;
        } else {
            return std::nullopt;
        }
    }

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const {
        auto kmer = decoded_cache_.TryGet(node);
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

        auto kmer = decoded_cache_.TryGet(node);
        edge_index edge = graph_ptr_->kmer_to_boss_index(node);

        boss_->call_incoming_to_target(boss_->bwd(edge), boss_->get_node_last_value(edge),
            [&](edge_index incoming_boss_edge) {
                assert(boss_->get_W(incoming_boss_edge) % boss_->alph_size
                        == boss_->get_node_last_value(edge));

                auto prev = graph_ptr_->boss_to_kmer_index(incoming_boss_edge);
                TAlphabet s = get_first_value(incoming_boss_edge);
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
            if (auto kmer = decoded_cache_.TryGet(node))
                update_node_next_from_kmer(*kmer, next, next_char);

            return next;

        } else {
            return npos;
        }
    }

    node_index traverse_back(node_index node, char prev_char) const {
        if (node_index prev = graph_ptr_->traverse_back(node, prev_char)) {
            if (auto kmer = decoded_cache_.TryGet(node))
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
                    decoded_cache_.Put(i, *prev);
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
        DBGSuccinctCached::add_sequence(sequence, on_insertion);
        decoded_cache_.Clear();
    }

    bool load(const std::string &filename) {
        if (!DBGSuccinctCached::load(filename))
            return false;

        decoded_cache_.Clear();
        return true;
    }

  private:
    mutable LRUCache<node_index, KmerType> decoded_cache_;

    inline static constexpr TAlphabet encode(char c) { return KmerExtractor::encode(c); }
    inline static constexpr KmerType to_kmer(std::string_view seq) {
        return KmerExtractor::sequence_to_kmer<KmerType>(seq);
    }

    KmerType update_node_next_from_kmer(KmerType kmer, node_index next, char c) const {
        kmer.to_next(graph_ptr_->get_k(), encode(c));
        decoded_cache_.Put(next, kmer);
        return kmer;
    }

    KmerType update_node_prev_from_kmer(KmerType kmer, node_index prev, char c) const {
        kmer.to_prev(graph_ptr_->get_k(), encode(c));
        decoded_cache_.Put(prev, kmer);
        return kmer;
    }
};

inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size = 1000) {
    if (graph->get_k() * kmer::KmerExtractorBOSS::bits_per_char <= 64) {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer64>>(graph, cache_size);
    } else if (graph->get_k() * kmer::KmerExtractorBOSS::bits_per_char <= 128) {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer128>>(graph, cache_size);
    } else {
        return std::make_shared<DBGSuccinctCachedImpl<kmer::KmerExtractorBOSS::Kmer256>>(graph, cache_size);
    }
}

inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(DeBruijnGraph &graph, size_t cache_size = 1000) {
    return make_cached_dbgsuccinct({ std::shared_ptr<DeBruijnGraph>{}, &graph }, cache_size);
}

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
