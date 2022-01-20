#include "dbg_succinct_cached.hpp"

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
std::string DBGSuccinctCachedViewImpl<KmerType>::get_node_sequence(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    // get the sequence from either the cache, or the underlying graph
    auto ret_val = kmer::KmerExtractorBOSS::kmer_to_sequence<KmerType>(
        std::get<0>(get_kmer_tuple(graph_->kmer_to_boss_index(node))), get_k()
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

    // TODO: we don't need this yet, but at some point update the stored second
    // value while traversing forwards

    edge_index bwd = 0;
    if (kmer)
        std::get<2>(*kmer) = 0;

    graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
        if (c != boss::BOSS::kSentinel) {
            if (kmer) {
                // update the decoded k-mer with the next character, then cache it
                CacheValue next_kmer = *kmer;
                edge_index next_edge = graph_->kmer_to_boss_index(next);

                if (!bwd)
                    bwd = boss_->pred_W(i, boss_->get_node_last_value(next_edge));

                std::get<0>(next_kmer).to_next(get_k(), encode(c));
                std::get<1>(next_kmer) = bwd;
                put_kmer(next_edge, std::move(next_kmer));
            }

            callback(next, c);
        }
    });
}

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>
::put_kmer(edge_index key, CacheValue value) const {
    assert(key);
    assert(!std::get<2>(value) || std::get<0>(value)[1] == boss_->get_node_last_value(std::get<2>(value)));
    assert(std::get<0>(value)[1] == boss_->get_minus_k_value(key, get_k() - 2).first);

    if (std::get<2>(value)) {
        if (std::optional<CacheValue> kmer = decoded_cache_.TryGet(key)) {
            if (std::get<2>(*kmer)) {
                assert(std::get<2>(*kmer) == std::get<2>(value));
                return;
            }
        }
    }

    decoded_cache_.Put(key, std::move(value));
}

template <typename KmerType>
auto DBGSuccinctCachedViewImpl<KmerType>
::get_kmer_tuple(edge_index i, edge_index child_hint, bool extract_first) const -> CacheValue {
    assert(i);
    std::optional<CacheValue> kmer = decoded_cache_.TryGet(i);
    if (kmer) {
        if (!extract_first || std::get<2>(*kmer))
            return *kmer;
    }

    if (child_hint) {
        if (std::optional<CacheValue> child_fetch = decoded_cache_.TryGet(child_hint)) {
            auto &[child_kmer, child_bwd, child_first] = *child_fetch;
            if (child_bwd == i && child_first) {
                edge_index first = boss_->bwd(child_first);
                assert(boss_->get_node_last_value(first)
                    == boss_->get_minus_k_value(i, get_k() - 2).first);
                child_kmer.to_prev(get_k(), boss_->get_node_last_value(first));
                CacheValue kmer_tuple{ child_kmer, boss_->bwd(i), first };
                put_kmer(i, kmer_tuple);
                return kmer_tuple;
            }
        }
    }

    // if the kmer is not cached, or if the result of bwd^{k - 1}(i) is desired
    // and not cached, then recompute both
    auto [seq, last_traversed_node, steps_remaining]
        = boss_->get_node_seq_with_last_traversed_node(i);

    edge_index last_node = 0;

    if (!steps_remaining) {
        last_node = last_traversed_node;
    } else if (extract_first) {
        last_node = boss_->get_minus_k_value(last_traversed_node,
                                             steps_remaining).second;
    }

    // append the last character and cache the result
    seq.push_back(boss_->get_W(i) % boss_->alph_size);

    // cache the result
    CacheValue kmer_tuple{ KmerType(seq), boss_->bwd(i), last_node };
    put_kmer(i, kmer_tuple);
    return kmer_tuple;
}

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>
::call_incoming_kmers(node_index node,
                      const IncomingEdgeCallback &callback) const {
    assert(node > 0 && node <= num_nodes());

    edge_index edge = graph_->kmer_to_boss_index(node);
    auto [kmer, bwd, first] = get_kmer_tuple(edge);
    TAlphabet last = boss_->get_node_last_value(edge);
    if (!bwd) {
        bwd = boss_->bwd(edge);
        put_kmer(edge, CacheValue{ kmer, bwd, first });
    }

    assert(bwd == boss_->bwd(edge));

    boss_->call_incoming_to_target(bwd, last, [&](edge_index incoming_boss_edge) {
        assert(boss_->get_W(incoming_boss_edge) % boss_->alph_size
            == boss_->get_node_last_value(edge));

        auto prev = graph_->boss_to_kmer_index(incoming_boss_edge);
        if (prev == npos)
            return;

        auto [prev_kmer, prev_parent, prev_first] = get_kmer_tuple(incoming_boss_edge, edge, true);
        TAlphabet s = prev_kmer[1];

        if (s != boss::BOSS::kSentinelCode)
            callback(prev, boss_->decode(s));
    });
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
    edge_index prev_edge = 0;
    std::deque<edge_index> first_buffer;
    size_t num_contiguous = 0;
    graph_->map_to_nodes_sequentially(sequence,
        [&](node_index i) {
            if (i) {
                edge_index edge = graph_->kmer_to_boss_index(i);
                if (prev) {
                    // update the previous k-mer
                    prev->to_next(k, encode(*last), encode(*(last - 1)));
                } else {
                    // initialize a new k-mer
                    prev = kmer::KmerExtractorBOSS::sequence_to_kmer<KmerType>(
                        std::string_view(begin, k)
                    );
                }

                assert(!prev_edge || boss_->pred_W(prev_edge, boss_->get_node_last_value(edge))
                    == boss_->bwd(edge));

                edge_index first = 0;

                if (prev_edge) {
                    edge_index actual_prev = boss_->pred_W(prev_edge, boss_->get_node_last_value(edge));
                    if (actual_prev != prev_edge) {
                        prev_edge = 0;
                        first_buffer.push_back(0);
                        num_contiguous = 0;
                    } else if (!boss_->is_single_incoming(prev_edge, boss_->get_node_last_value(edge))) {
                        prev_edge = 0;
                        first_buffer.push_back(0);
                        num_contiguous = 0;
                    } else {
                        first_buffer.push_back(prev_edge);
                        if (prev_edge) {
                            ++num_contiguous;
                        } else {
                            num_contiguous = 0;
                        }
                        if (first_buffer.size() == k - 2) {
                            if (num_contiguous >= k - 2) {
                                first = first_buffer.front();
                                // std::cerr << boss_->get_node_str(first) << " " << boss_->get_node_str(boss_->get_minus_k_value(edge, k - 2).second) << "\n";
                                // assert(first == boss_->get_minus_k_value(edge, k - 2).second);
                            }

                            first_buffer.pop_front();
                        }
                    }
                } else {
                    first_buffer.push_back(0);
                    num_contiguous = 0;
                }



                // cache the result
                put_kmer(edge, CacheValue{*prev, prev_edge, first });

                prev_edge = edge;

            } else {
                // reset the k-mer to null
                prev = std::nullopt;
                prev_edge = 0;
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
