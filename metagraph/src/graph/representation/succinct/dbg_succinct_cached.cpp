#include "dbg_succinct_cached.hpp"

#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {


template <typename Graph>
DBGSuccinctCachedView::DBGSuccinctCachedView(Graph&& graph, size_t cache_size)
      : DBGWrapper(std::forward<Graph>(graph)), boss_(&graph_->get_boss()),
        cache_size_(cache_size), rev_comp_prev_cache_(cache_size_),
        rev_comp_next_cache_(cache_size_) {}

template DBGSuccinctCachedView::DBGSuccinctCachedView(const DBGSuccinct&, size_t);
template DBGSuccinctCachedView::DBGSuccinctCachedView(const std::shared_ptr<const DBGSuccinct>&, size_t);
template DBGSuccinctCachedView::DBGSuccinctCachedView(std::shared_ptr<const DBGSuccinct>&, size_t);


/**
 * Delegated methods
 */
#define DELEGATE_METHOD_IMPL(RETURN_TYPE, METHOD, ARG_TYPE, ARG_NAME) \
auto DBGSuccinctCachedView::METHOD(ARG_TYPE ARG_NAME) const -> RETURN_TYPE { \
    return graph_->METHOD(ARG_NAME); \
} \

DELEGATE_METHOD_IMPL(void, call_kmers, const std::function<void(node_index, const std::string&)> &, callback)
DELEGATE_METHOD_IMPL(void, call_source_nodes, const std::function<void(node_index)> &, callback)

DELEGATE_METHOD_IMPL(uint64_t, num_nodes, , )
DELEGATE_METHOD_IMPL(uint64_t, max_index, , )
DELEGATE_METHOD_IMPL(size_t, outdegree, node_index, node)
DELEGATE_METHOD_IMPL(bool, has_single_outgoing, node_index, node)
DELEGATE_METHOD_IMPL(bool, has_multiple_outgoing, node_index, node)
DELEGATE_METHOD_IMPL(size_t, indegree, node_index, node)
DELEGATE_METHOD_IMPL(bool, has_no_incoming, node_index, node)
DELEGATE_METHOD_IMPL(bool, has_single_incoming, node_index, node)
DELEGATE_METHOD_IMPL(node_index, kmer_to_node, std::string_view, kmer)

void DBGSuccinctCachedView::call_nodes(const std::function<void(node_index)> &callback,
                                       const std::function<bool()> &stop_early) const {
    graph_->call_nodes(callback, stop_early);
}

void DBGSuccinctCachedView::call_sequences(const CallPath &callback,
                                           size_t num_threads,
                                           bool kmers_in_single_form ) const {
    graph_->call_sequences(callback, num_threads, kmers_in_single_form);
}

void DBGSuccinctCachedView::call_unitigs(const CallPath &callback,
                                         size_t num_threads,
                                         size_t min_tip_size,
                                         bool kmers_in_single_form) const {
    graph_->call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
}

bool DBGSuccinctCachedView
::find(std::string_view sequence, double discovery_fraction) const {
    return graph_->find(sequence, discovery_fraction);
}

void DBGSuccinctCachedView::map_to_nodes(std::string_view sequence,
                                         const std::function<void(node_index)> &callback,
                                         const std::function<bool()> &terminate) const {
    graph_->map_to_nodes(sequence, callback, terminate);
}

void DBGSuccinctCachedView
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    graph_->adjacent_outgoing_nodes(node, callback);
}

void DBGSuccinctCachedView
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    graph_->adjacent_incoming_nodes(node, callback);
}

auto DBGSuccinctCachedView
::traverse_back(node_index node, char prev_char) const -> node_index {
    // TODO: implemented a cached version of this later if needed
    return graph_->traverse_back(node, prev_char);
}

void DBGSuccinctCachedView::traverse(node_index start,
                                     const char *begin,
                                     const char *end,
                                     const std::function<void(node_index)> &callback,
                                     const std::function<bool()> &terminate) const {
    graph_->traverse(start, begin, end, callback, terminate);
}

bool DBGSuccinctCachedView::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || num_nodes() != other.num_nodes()
            || get_mode() != other.get_mode())
        return false;

    const auto *other_succ = dynamic_cast<const DBGSuccinct*>(&other);
    if (!other_succ) {
        if (const auto *other_cached = dynamic_cast<const DBGSuccinctCachedView*>(&other))
            other_succ = &other_cached->get_graph();
    }

    if (other_succ)
        return graph_->operator==(*other_succ);

    throw std::runtime_error("Not implemented");
    return false;
}

auto DBGSuccinctCachedView
::get_rev_comp_boss_next_node(node_index node) const -> edge_index {
    // 78% effective
    edge_index ret_val = 0;
    size_t chars_unmatched = 0;
    if (auto fetch = rev_comp_next_cache_.TryGet(node)) {
        std::tie(ret_val, chars_unmatched) = *fetch;
        if (!ret_val && !chars_unmatched)
            return 0;
    } else {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  GAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  CTCAAGCAGAAGACGGCATACGAGATCCTC
        std::string rev_seq = get_node_sequence(node).substr(1, get_k() - 1);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_next_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss_->encode(rev_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge == edge_2);

        if (end == encoded.end()) {
            ret_val = edge;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        edge_index start = boss_->pred_W(graph_->kmer_to_boss_index(node),
                                         complement(encoded.front()));
        boss_->call_incoming_to_target(start, complement(encoded.front()),
                                       [&](edge_index adj) {
            if (node_index a = graph_->boss_to_kmer_index(adj))
                rev_comp_next_cache_.Put(a, stored_val);
        });
    }

    if (!ret_val && chars_unmatched > 1) {
        // we only matched a suffix
        // $$CTCAAGCAGAAGACGGCATACGAGATCC
        // so if we try to match the next k-mer
        // GAGGATCTCGTATGCCGTCTTCTGCTTGAGX
        // corresponding to
        //  xCTCAAGCAGAAGACGGCATACGAGATCCT
        // we may find the subrange
        // $xCTCAAGCAGAAGACGGCATACGAGATCC
        // or something shorter
        size_t next_chars_unmatched = chars_unmatched - 1;
        auto stored_val = std::make_pair(0, next_chars_unmatched);
        adjacent_outgoing_nodes(node, [&](node_index next) {
            rev_comp_next_cache_.Put(next, stored_val);
#ifndef NDEBUG
            std::string test_seq = graph_->get_node_sequence(next).substr(1, get_k() - 1);
            ::reverse_complement(test_seq.begin(), test_seq.end());
            auto encoded = boss_->encode(test_seq);
            auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
            assert(end != encoded.end());
            if (next_chars_unmatched > static_cast<size_t>(encoded.end() - end)) {
                common::logger->error("Unmatched chars incorrectly stored: {} vs. {}",
                                      next_chars_unmatched, encoded.end() - end);
                assert(false);
            }
#endif
        });
    } else if (ret_val) {
        // traverse node forward by X, then ret_val back by comp(X)
        std::vector<std::pair<node_index, edge_index>> parents(alphabet().size());
        size_t parents_count = 0;
        edge_index start = 0;
        call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                parents[boss_->encode(c)].first = next;
                edge_index next_edge = graph_->kmer_to_boss_index(next);
                if (boss_->get_last(next_edge))
                    start = next_edge;

                ++parents_count;
            }
        });

        if (parents_count) {
            edge_index ret_prev = boss_->bwd(ret_val);
            TAlphabet w = boss_->get_node_last_value(ret_val);
            boss_->call_incoming_to_target(ret_prev, w, [&](edge_index prev_edge) {
                if (boss_->get_last(prev_edge)) {
                    TAlphabet s = complement(get_first_value(prev_edge));
                    if (parents[s].first)
                        rev_comp_next_cache_.Put(parents[s].first, std::make_pair(prev_edge, 0));
                }
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = graph_->get_node_sequence(node).substr(1, get_k() - 1);
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss_->get_node_str(ret_val));
    } else {
        auto encoded = boss_->encode(test_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end());
        if (chars_unmatched > static_cast<size_t>(encoded.end() - end)) {
            common::logger->error("Unmatched chars incorrectly stored: {} vs. {}",
                                  chars_unmatched, encoded.end() - end);
            assert(false);
        }
    }
#endif

    return ret_val;
}

auto DBGSuccinctCachedView
::get_rev_comp_boss_prev_node(node_index node) const -> edge_index {
    // 8% effective
    edge_index ret_val = 0;
    size_t chars_unmatched = 0;
    if (auto fetch = rev_comp_prev_cache_.TryGet(node)) {
        std::tie(ret_val, chars_unmatched) = *fetch;
        if (!ret_val && !chars_unmatched)
            return 0;
    } else {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //-> AGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //-> TCAAGCAGAAGACGGCATACGAGATCCTCT
        std::string rev_seq = get_node_sequence(node);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_prev_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        rev_seq.pop_back();
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss_->encode(rev_seq);
        auto [edge_1, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge_1 == edge_2);
        if (end == encoded.end()) {
            ret_val = edge_1;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        boss_->call_outgoing(boss_->succ_last(graph_->kmer_to_boss_index(node)),
                             [&](edge_index adj) {
            if (node_index a = graph_->boss_to_kmer_index(adj))
                rev_comp_prev_cache_.Put(a, stored_val);
        });
    }

    // TODO: given an LCS array, we can do the same check as get_rev_comp_boss_next_node
    //       to find neighbouring nodes with rev comp matches
    //       e.g., we looked for
    //       AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
    //       TCAAGCAGAAGACGGCATACGAGATCCTCT
    //       and matched a suffix
    //       XXXXXXXXTCAAGCAGAAGACGGCATACGA
    //       so if we try to match the previous k-mer
    //       XAGAGGATCTCGTATGCCGTCTTCTGCTTGA
    //       corresponding to
    //        CAAGCAGAAGACGGCATACGAGATCCTCTx
    //       we can use the LCS array to delete the first character from the match
    //       and work from there

    if (ret_val) {
        // traverse ret_val forward by X, then node back by comp(X)
        std::vector<std::pair<node_index, edge_index>> parents(alphabet().size());
        size_t parents_count = 0;
        if (graph_->boss_to_kmer_index(ret_val)) {
            boss_->call_outgoing(ret_val, [&](edge_index next_edge) {
                if (TAlphabet w = boss_->get_W(next_edge) % boss_->alph_size) {
                    parents[w].second = boss_->fwd(boss_->pred_W(ret_val, w), w);
                    ++parents_count;
                }
            });
        }

        if (parents_count) {
            call_incoming_kmers(node, [&](node_index prev, char c) {
                TAlphabet s = complement(boss_->encode(c));
                if (parents[s].second)
                    rev_comp_prev_cache_.Put(prev, std::make_pair(parents[s].second, 0));
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = graph_->get_node_sequence(node);
    test_seq.pop_back();
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss_->get_node_str(ret_val));
    } else {
        auto encoded = boss_->encode(test_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end());
        assert(chars_unmatched <= static_cast<size_t>(encoded.end() - end));
    }
#endif

    return ret_val;
}


template <typename KmerType>
template <typename Graph>
DBGSuccinctCachedViewImpl<KmerType>
::DBGSuccinctCachedViewImpl(Graph&& graph, size_t cache_size)
      : DBGSuccinctCachedView(std::forward<Graph>(graph), cache_size),
        decoded_cache_(cache_size_) {}

template <typename KmerType>
void DBGSuccinctCachedViewImpl<KmerType>::put_decoded_node(node_index node,
                                                           std::string_view seq) const {
    assert(node > 0 && node <= num_nodes());
    assert(seq.size() == graph_->get_k());
    assert(graph_->get_node_sequence(node) == seq);

    put_kmer(graph_->kmer_to_boss_index(node), CacheValue{ to_kmer(seq), std::nullopt });
}

template <typename KmerType>
std::string DBGSuccinctCachedViewImpl<KmerType>::get_node_sequence(node_index node) const {
    assert(node > 0 && node <= num_nodes());

    // get the sequence from either the cache, or the underlying graph
    auto ret_val = KmerExtractor::kmer_to_sequence<KmerType>(
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
auto DBGSuccinctCachedViewImpl<KmerType>
::traverse(node_index node, char next_char) const -> node_index {
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
                            const std::function<bool()> &terminate) const {
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

template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer64>;
template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer128>;
template class DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer256>;

template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer64>::DBGSuccinctCachedViewImpl(const DBGSuccinct&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer128>::DBGSuccinctCachedViewImpl(const DBGSuccinct&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer256>::DBGSuccinctCachedViewImpl(const DBGSuccinct&, size_t);

template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer64>::DBGSuccinctCachedViewImpl(const std::shared_ptr<const DBGSuccinct>&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer128>::DBGSuccinctCachedViewImpl(const std::shared_ptr<const DBGSuccinct>&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer256>::DBGSuccinctCachedViewImpl(const std::shared_ptr<const DBGSuccinct>&, size_t);

template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer64>::DBGSuccinctCachedViewImpl(std::shared_ptr<const DBGSuccinct>&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer128>::DBGSuccinctCachedViewImpl(std::shared_ptr<const DBGSuccinct>&, size_t);
template DBGSuccinctCachedViewImpl<kmer::KmerExtractorBOSS::Kmer256>::DBGSuccinctCachedViewImpl(std::shared_ptr<const DBGSuccinct>&, size_t);

} // namespace graph
} // namespace mtg
