#ifndef __DBG_SUCCINCT_CACHED__
#define __DBG_SUCCINCT_CACHED__

#include <cache.hpp>
#include <lru_cache_policy.hpp>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {


class DBGSuccinctCached : public DeBruijnGraph {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    template <typename Key, typename Value>
    using LRUCache = caches::fixed_sized_cache<Key, Value, caches::LRUCachePolicy<Key>>;

    typedef LRUCache<edge_index, edge_index> BwdFirstCache;

    DBGSuccinctCached(const DBGSuccinct &dbg, size_t cache_size)
          : graph_(dbg), boss_(graph_.get_boss()), cache_size_(cache_size),
            bwd_first_cache_(cache_size_) {}

    DBGSuccinctCached(const DeBruijnGraph &graph, size_t cache_size)
          : DBGSuccinctCached(dynamic_cast<const DBGSuccinct&>(graph), cache_size) {}

    virtual ~DBGSuccinctCached() {}

    virtual const DeBruijnGraph& get_base_graph() const override final { return graph_.get_base_graph(); }

    const DBGSuccinct& get_dbg_succ() const { return graph_; }

    virtual void put_decoded_node(node_index node, const std::string &seq) const = 0;
    virtual std::optional<std::string> get_decoded_node(node_index node) const = 0;

    inline TAlphabet get_first_value(edge_index i) const {
        TAlphabet c;
        edge_index first;

        if (auto fetch = bwd_first_cache_.TryGet(i)) {
            c = boss_.get_node_last_value(*fetch);
            first = *fetch;
        } else {
            std::tie(c, first) = boss_.get_minus_k_value(i, get_k() - 2);
            bwd_first_cache_.Put(i, first);
        }

        edge_index bwd_i = boss_.bwd(i);
        auto fetch_bwd = bwd_first_cache_.TryGet(bwd_i);
        if (!fetch_bwd)
            bwd_first_cache_.Put(bwd_i, boss_.bwd(first));

        return c;
    }

    virtual std::string get_node_sequence(node_index node) const override final {
        if (auto fetch = get_decoded_node(node)) {
            return *fetch;
        } else {
            std::string seq = graph_.get_node_sequence(node);
            put_decoded_node(node, seq);
            return seq;
        }
    }

    virtual void call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const override final {
        assert(node > 0 && node <= num_nodes());
        std::string seq;
        if (auto fetch = get_decoded_node(node))
            seq = std::string(1, boss::BOSS::kSentinel) + fetch->substr(0, get_k() - 1);

        edge_index edge = graph_.kmer_to_boss_index(node);

        boss_.call_incoming_to_target(boss_.bwd(edge), boss_.get_node_last_value(edge),
            [&](edge_index incoming_boss_edge) {
                assert(boss_.get_W(incoming_boss_edge) % boss_.alph_size
                        == boss_.get_node_last_value(edge));

                auto prev = graph_.boss_to_kmer_index(incoming_boss_edge);
                if (prev != npos) {
                    char c = boss_.decode(get_first_value(incoming_boss_edge));
                    if (seq.size()) {
                        seq[0] = c;
                        put_decoded_node(prev, seq);
                    }

                    callback(prev, c);
                }
            }
        );
    }

    virtual void call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const override final {
        std::string seq;
        if (auto fetch = get_decoded_node(node))
            seq = fetch->substr(1) + std::string(1, boss::BOSS::kSentinel);

        graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (seq.size()) {
                seq.back() = c;
                put_decoded_node(next, seq);
            }

            callback(next, c);
        });
    }

    virtual node_index traverse(node_index node, char next_char) const override final {
        if (node_index ret = graph_.traverse(node, next_char)) {
            if (auto fetch = get_decoded_node(node))
                put_decoded_node(ret, fetch->substr(1) + next_char);

            return ret;
        } else {
            return npos;
        }
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        if (node_index ret = graph_.traverse_back(node, prev_char)) {
            if (auto fetch = get_decoded_node(node)) {
                put_decoded_node(ret,
                    std::string(1, prev_char) + fetch->substr(0, get_k() - 1)
                );
            }

            return ret;
        } else {
            return npos;
        }
    }

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override final {
        size_t k = get_k();
        auto it = sequence.begin();
        graph_.map_to_nodes_sequentially(sequence,
            [&](node_index i) {
                if (i)
                    put_decoded_node(i, std::string(it, it + k));

                callback(i);
                ++it;
            },
            terminate
        );
    }

    virtual size_t get_k() const override final { return graph_.get_k(); }

    virtual Mode get_mode() const override final { return graph_.get_mode(); }

    virtual size_t outdegree(node_index node) const override final {
        return graph_.outdegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override final {
        return graph_.has_single_outgoing(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override final {
        return graph_.has_multiple_outgoing(node);
    }

    virtual size_t indegree(node_index node) const override final {
        return graph_.indegree(node);
    }

    virtual bool has_no_incoming(node_index node) const override final {
        return graph_.has_no_incoming(node);
    }

    virtual bool has_single_incoming(node_index node) const override final {
        return graph_.has_single_incoming(node);
    }

    virtual const std::string& alphabet() const override final { return graph_.alphabet(); }


    virtual void add_sequence(std::string_view /* sequence */,
                              const std::function<void(node_index)> & /* on_insertion */ = [](uint64_t) {}) override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override final {
        graph_.map_to_nodes(sequence, callback, terminate);
    }

    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final {
        graph_.adjacent_outgoing_nodes(node, callback);
    }
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final {
        graph_.adjacent_incoming_nodes(node, callback);
    }

    virtual uint64_t num_nodes() const override final { return graph_.num_nodes(); }

    virtual bool load(const std::string &) override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void serialize(const std::string &) const override final {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string file_extension() const override final {
        throw std::runtime_error("Not implemented");
    }

  protected:
    const DBGSuccinct &graph_;
    const boss::BOSS &boss_;
    size_t cache_size_;

  private:
    mutable BwdFirstCache bwd_first_cache_;
};

template <typename KmerType>
class DBGSuccinctCachedImpl : public DBGSuccinctCached {
  public:
    typedef LRUCache<node_index, KmerType> DecodedCache;

    template <typename... Args>
    DBGSuccinctCachedImpl(Args&&... args)
          : DBGSuccinctCached(std::forward<Args>(args)...),
            decoded_cache_(cache_size_ * graph_.alphabet().size()) {}

    void put_decoded_node(node_index node, const std::string &seq) const {
        decoded_cache_.Put(node, seq);
    }

    std::optional<std::string> get_decoded_node(node_index node) const {
        return decoded_cache_.TryGet(node);
    }

  private:
    mutable DecodedCache decoded_cache_;
};

inline std::shared_ptr<DBGSuccinctCached>
make_cached_dbgsuccinct(const DeBruijnGraph &graph, size_t cache_size = 1000) {
    return std::make_shared<DBGSuccinctCachedImpl<std::string>>(graph, cache_size);
}

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_CACHED__
