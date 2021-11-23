#include "dbg_succinct_cached.hpp"

#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {

/**
 * Overridden methods
 */
bool DBGSuccinct::CachedView::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || num_nodes() != other.num_nodes()
            || get_mode() != other.get_mode())
        return false;

    const auto *other_succ = dynamic_cast<const DBGSuccinct*>(&other);
    if (!other_succ) {
        if (const auto *other_cached = dynamic_cast<const DBGSuccinct::CachedView*>(&other))
            other_succ = &other_cached->get_graph();
    }

    if (other_succ)
        return graph_->operator==(*other_succ);

    throw std::runtime_error("Not implemented");
    return false;
}

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>::put_decoded_edge(edge_index edge,
                                                           std::string_view seq) const {
    assert(edge > 0 && edge <= boss_->num_edges());
    assert(seq.size() == graph_->get_k());
    assert(graph_->get_node_sequence(graph_->boss_to_kmer_index(edge)) == seq);

    put_kmer(edge, CacheValue{ to_kmer(seq), std::nullopt });
}

template <typename KmerType>
std::string DBGSuccinctCachedViewImpl<KmerType>::get_node_sequence(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    // get the sequence from either the cache, or the underlying graph
    auto ret_val = kmer::KmerExtractorBOSS::kmer_to_sequence<KmerType>(
        get_kmer_pair(graph_->kmer_to_boss_index(node)).first, get_k()
    );

    assert(ret_val == graph_->get_node_sequence(node));
    return ret_val;
}

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>
::call_outgoing_kmers(node_index node,
                      const OutgoingEdgeCallback &callback) const {
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

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>
::call_incoming_kmers(node_index node,
                      const IncomingEdgeCallback &callback) const {
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

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate,
                            const std::function<bool()> &skip) const {
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
        terminate,
        skip
    );
}


template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer64>;
template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer128>;
template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer256>;

} // namespace graph
} // namespace mtg
