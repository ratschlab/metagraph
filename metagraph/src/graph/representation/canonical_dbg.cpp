#include "canonical_dbg.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/kmer_extractor.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;

inline const DBGSuccinct* get_dbg_succ(const DeBruijnGraph &graph) {
    return dynamic_cast<const DBGSuccinct*>(&graph.get_base_graph());
}

template <typename Graph>
CanonicalDBG::CanonicalDBG(Graph&& graph, size_t cache_size)
      : DBGWrapper<DeBruijnGraph>(std::forward<Graph>(graph)),
        is_palindrome_cache_(cache_size),
        rev_comp_prev_cache_(cache_size),
        rev_comp_next_cache_(cache_size) {
    static_assert(!std::is_same_v<Graph, std::shared_ptr<CanonicalDBG>>);
    static_assert(!std::is_same_v<Graph, std::shared_ptr<const CanonicalDBG>>);
    flush();
}

template CanonicalDBG::CanonicalDBG(std::shared_ptr<DeBruijnGraph>&&, size_t);
template CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph>&&, size_t);
template CanonicalDBG::CanonicalDBG(std::shared_ptr<DeBruijnGraph>&, size_t);
template CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph>&, size_t);
template CanonicalDBG::CanonicalDBG(std::shared_ptr<DBGSuccinct>&, size_t);

void CanonicalDBG::flush() {
    if (graph_->get_mode() != DeBruijnGraph::PRIMARY) {
        logger->error("Only primary graphs can be wrapped in CanonicalDBG");
        exit(1);
    }

    is_palindrome_cache_.Clear();
    rev_comp_next_cache_.Clear();
    rev_comp_prev_cache_.Clear();

    offset_ = graph_->max_index();
    k_odd_ = (graph_->get_k() % 2);
    has_sentinel_ = false;
    alphabet_encoder_.fill(graph_->alphabet().size());

    for (size_t i = 0; i < graph_->alphabet().size(); ++i) {
        alphabet_encoder_[graph_->alphabet()[i]] = i;
        if (graph_->alphabet()[i] == boss::BOSS::kSentinel)
            has_sentinel_ = true;
    }
}

void CanonicalDBG
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate,
                            const std::function<bool()> &skip) const {
    if (sequence.size() < get_k())
        return;

    std::vector<node_index> path;
    path.reserve(sequence.size() - get_k() + 1);

    // map until the first mismatch
    bool stop = false;
    graph_->map_to_nodes_sequentially(sequence,
        [&](node_index node) {
            if (node) {
                path.push_back(node);
            } else {
                stop = true;
            }
        },
        [&]() { return stop; },
        skip
    );

    for (node_index node : path) {
        if (terminate())
            return;

        callback(node);
    }

    // trim the mapped prefix
    sequence = sequence.substr(path.size());
    path.resize(0);
    if (sequence.size() < get_k())
        return;

    std::string rev_seq(sequence);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    // map the reverse-complement
    std::vector<node_index> rev_path = map_sequence_to_nodes(*graph_, rev_seq);

    // map the forward
    const auto *dbg_succ = get_dbg_succ(*graph_);
    if (dbg_succ && k_odd_) {
        // if it's a boss table with odd k (without palindromic k-mers),
        // we can skip k-mers that have been found in the rev-compl sequence

        // the initial forward mapping stopped on this k-mer,
        // hence it's missing and we skip it
        path.push_back(npos);
        auto it = rev_path.rbegin() + 1;

        auto map_remaining = [&](const auto &graph) {
            graph.map_to_nodes_sequentially(sequence.substr(1),
                [&](node_index next) {
                    path.push_back(next);
                    ++it;
                },
                []() { return false; },
                [&]() { return *it; }
            );
        };

        if (const auto *cached = dynamic_cast<const DBGSuccinct::CachedView*>(graph_.get())) {
            map_remaining(*cached);
        } else {
            map_remaining(*dbg_succ);
        }

        assert(it == rev_path.rend());

    } else {
        path = map_sequence_to_nodes(*graph_, sequence);
    }

    assert(path.size() == rev_path.size());

    auto it = rev_path.rbegin();
    for (auto jt = path.begin(); jt != path.end(); ++jt, ++it) {
        if (terminate())
            return;

        if (skip()) {
            callback(npos);
            continue;
        }

        if (*jt != npos) {
            callback(*jt);
        } else if (*it != npos) {
            callback(*it + offset_);
        } else {
            callback(npos);
        }
    }
}

void CanonicalDBG::map_to_nodes(std::string_view sequence,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    map_to_nodes_sequentially(sequence, [&](node_index i) {
        callback(i != npos ? get_base_node(i) : i);
    }, terminate);
}

void CanonicalDBG::append_next_rc_nodes(node_index node,
                                        std::vector<node_index> &children) const {
    /**
     *
     * find children of node by searching for parents of its reverse complement
     * e.g., node = ATGGCT. Find TGGCTA and TGGCTT by looking for  TAGCCA and AAGCCA
     *         TGGCTA      TAGCCA
     *        /                  \
     *  ATGGCT         ->         AGCCAT
     *        \                  /
     *         TGGCTT      AAGCCA
     */

    const auto &alphabet = graph_->alphabet();

    // for each n, check for nAGCCA. If found, define and store the index for
    // TGGCTrc(n) as index(nAGCCA) + offset_

    // callback for DBGSuccinct and DBGSuccinct::CachedView
    auto next_callback_succ = [&](node_index next, boss::BOSS::TAlphabet c) {
        if (c == boss::BOSS::kSentinelCode)
            return;

        c = kmer::KmerExtractorBOSS::complement(c);

        if (children[c] == npos) {
            children[c] = next + offset_;

            return;
        }

        if (k_odd_) {
            logger->error(
                "Primary graph contains both forward and reverse complement: {} {} -> {} {}\t{} {}",
                node, graph_->get_node_sequence(node),
                children[c], graph_->get_node_sequence(children[c]),
                next, graph_->get_node_sequence(next));
            exit(1);
        }

        is_palindrome_cache_.Put(next, true);
    };

    if (dynamic_cast<const DBGSuccinct::CachedView*>(graph_.get())) {
        call_outgoing_from_rev_comp(node, next_callback_succ);
        return;
    }

    //        rshift    rc
    // ATGGCT -> TGGCT* -> *AGCCA
    std::string rev_seq = get_node_sequence(node).substr(1) + std::string(1, '\0');
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    assert(rev_seq[0] == '\0');

    if (const auto *dbg_succ = get_dbg_succ(*graph_)) {
        dbg_succ->call_nodes_with_suffix_matching_longest_prefix(
            std::string_view(&rev_seq[1], get_k() - 1),
            [&](node_index next, uint64_t /* match_length */) {
                next_callback_succ(next,
                                   dbg_succ->get_boss().get_minus_k_value(
                                       dbg_succ->kmer_to_boss_index(next), get_k() - 2
                                   ).first);
            },
            get_k() - 1
        );

        return;
    }

    for (size_t c = 0; c < alphabet.size(); ++c) {
        // Do the checks by directly mapping the sequences of the desired k-mers.
        // For non-DBGSuccinct graphs, this should be fast enough.
        if (alphabet[c] != boss::BOSS::kSentinel && children[c] == npos) {
            rev_seq[0] = complement(alphabet[c]);
            node_index next = graph_->kmer_to_node(rev_seq);
            if (next != npos)
                children[c] = next + offset_;
        }
    }
}

void CanonicalDBG
::call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (node > offset_) {
        call_incoming_kmers(node - offset_, [&](node_index next, char c) {
            c = complement(c);
            callback(reverse_complement(next), c);
            assert(traverse(node, c) == reverse_complement(next));
        });
        return;
    }

    const auto &alphabet = graph_->alphabet();

    std::vector<node_index> children(alphabet.size(), npos);
    size_t max_num_edges_left = children.size() - has_sentinel_;

    graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
        if (c != boss::BOSS::kSentinel) {
            children[alphabet_encoder_[c]] = next;
            --max_num_edges_left;
        }
    });

    if (max_num_edges_left)
        append_next_rc_nodes(node, children);

    for (size_t c = 0; c < children.size(); ++c) {
        if (children[c] != npos) {
            callback(children[c], alphabet[c]);
            assert(traverse(node, alphabet[c]) == children[c]);
        }
    }
}

void CanonicalDBG::append_prev_rc_nodes(node_index node,
                                        std::vector<node_index> &parents) const {
    /**
     * find parents of node by searching for children of its reverse complement
     * e.g., node = AGCCAT. Find TAGCCA and AAGCCA by looking for TGGCTA and TGGCTT.
     *  TAGCCA                    TGGCTA
     *        \                  /
     *         AGCCAT  ->  ATGGCT
     *        /                  \
     *  AAGCCA                    TGGCTT
     */

    const auto &alphabet = graph_->alphabet();

    // for each n, check for TGGCTn. If found, define and store the index for
    // rc(n)AGCCA as index(TGGCTn) + offset_

    // callback for DBGSuccinct and DBGSuccinct::CachedView
    auto prev_callback_succ = [&](node_index prev, boss::BOSS::TAlphabet c) {
        if (prev == npos || c == boss::BOSS::kSentinelCode)
            return;

        c = kmer::KmerExtractorBOSS::complement(c);

        if (parents[c] == npos) {
            parents[c] = prev + offset_;

            return;
        }

        if (k_odd_) {
            logger->error(
                "Primary graph contains both forward and reverse complement: {} {} -> {} {}\t{} {}",
                node, graph_->get_node_sequence(node),
                parents[c], graph_->get_node_sequence(parents[c]),
                prev, graph_->get_node_sequence(prev)
            );
            exit(1);
        }

        is_palindrome_cache_.Put(prev, true);
    };

    if (dynamic_cast<const DBGSuccinct::CachedView*>(graph_.get())) {
        call_incoming_to_rev_comp(node, prev_callback_succ);
        return;
    }

    //        lshift    rc
    // AGCCAT -> *AGCCA -> TGGCT*
    std::string rev_seq = std::string(1, '\0') + get_node_sequence(node).substr(0, get_k() - 1);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    assert(rev_seq.back() == '\0');

    if (const auto *dbg_succ = get_dbg_succ(*graph_)) {
        // Find the BOSS node TGGCT and iterate through all of its outdoing edges.
        // Then, convert the edge indices to get the DBGSuccinct node indices
        const auto &boss = dbg_succ->get_boss();
        auto encoded = boss.encode(std::string_view(rev_seq.data(), get_k() - 1));
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge == edge_2);
        if (end == encoded.end()) {
            boss.call_outgoing(edge, [&](boss::BOSS::edge_index adjacent_edge) {
                prev_callback_succ(dbg_succ->boss_to_kmer_index(adjacent_edge),
                                   boss.get_W(adjacent_edge) % boss.alph_size);
            });
        }

        return;
    }

    for (size_t c = 0; c < alphabet.size(); ++c) {
        // Do the checks by directly mapping the sequences of the desired k-mers.
        // For non-DBGSuccinct graphs, this should be fast enough.
        if (alphabet[c] != boss::BOSS::kSentinel && parents[c] == npos) {
            rev_seq.back() = complement(alphabet[c]);
            node_index prev = graph_->kmer_to_node(rev_seq);
            if (prev != npos)
                parents[c] = prev + offset_;
        }
    }
}

void CanonicalDBG
::call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (node > offset_) {
        call_outgoing_kmers(node - offset_, [&](node_index prev, char c) {
            c = complement(c);
            callback(reverse_complement(prev), c);
            assert(traverse_back(node, c) == reverse_complement(prev));
        });
        return;
    }

    const auto &alphabet = graph_->alphabet();

    std::vector<node_index> parents(alphabet.size(), npos);
    size_t max_num_edges_left = parents.size() - has_sentinel_;

    graph_->call_incoming_kmers(node, [&](node_index prev, char c) {
        if (c != boss::BOSS::kSentinel) {
            parents[alphabet_encoder_[c]] = prev;
            --max_num_edges_left;
        }
    });

    if (max_num_edges_left)
        append_prev_rc_nodes(node, parents);

    for (size_t c = 0; c < parents.size(); ++c) {
        if (parents[c] != npos) {
            callback(parents[c], alphabet[c]);
            assert(traverse_back(node, alphabet[c]) == parents[c]);
        }
    }
}

void CanonicalDBG
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    // TODO: better implementation
    call_outgoing_kmers(node, [&](node_index i, char) { callback(i); });
}

void CanonicalDBG
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    // TODO: better implementation
    call_incoming_kmers(node, [&](node_index i, char) { callback(i); });
}

size_t CanonicalDBG::outdegree(node_index node) const {
    // TODO: better implementation
    size_t outdegree = 0;
    adjacent_outgoing_nodes(node, [&](node_index) { ++outdegree; });
    return outdegree;
}

size_t CanonicalDBG::indegree(node_index node) const {
    // TODO: better implementation
    size_t indegree = 0;
    adjacent_incoming_nodes(node, [&](node_index) { ++indegree; });
    return indegree;
}

void CanonicalDBG::call_sequences(const CallPath &callback,
                                  size_t num_threads,
                                  bool kmers_in_single_form) const {
    if (kmers_in_single_form) {
        graph_->call_sequences(callback, num_threads, false);
    } else {
        // TODO: port over implementation from DBGSuccinct to DeBruijnGraph
        DeBruijnGraph::call_sequences(callback, num_threads, false);
    }
}

void CanonicalDBG::call_unitigs(const CallPath &callback,
                                size_t num_threads,
                                size_t min_tip_size,
                                bool kmers_in_single_form) const {
    // TODO: port over implementation from DBGSuccinct to DeBruijnGraph
    DeBruijnGraph::call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
}

std::string CanonicalDBG::get_node_sequence(node_index index) const {
    assert(index <= offset_ * 2);
    node_index node = get_base_node(index);
    std::string seq = graph_->get_node_sequence(node);

    if (node != index)
        ::reverse_complement(seq.begin(), seq.end());

    return seq;
}

DeBruijnGraph::node_index CanonicalDBG::traverse(node_index node, char next_char) const {
    assert(node <= offset_ * 2);
    if (node > offset_) {
        node = traverse_back(node - offset_, complement(next_char));
        return node != npos ? reverse_complement(node) : npos;
    } else {
        node_index next = graph_->traverse(node, next_char);
        if (next != npos)
            return next;

        std::string rev_seq = get_node_sequence(node).substr(1) + next_char;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        next = graph_->kmer_to_node(rev_seq);
        return next != npos ? reverse_complement(next) : next;
    }
}

DeBruijnGraph::node_index CanonicalDBG::traverse_back(node_index node,
                                                      char prev_char) const {
    assert(node <= offset_ * 2);
    if (node > offset_) {
        node = traverse(node - offset_, complement(prev_char));
        return node != npos ? reverse_complement(node) : npos;
    } else {
        node_index prev = graph_->traverse_back(node, prev_char);
        if (prev != npos)
            return prev;

        std::string rev_seq = std::string(1, prev_char)
            + get_node_sequence(node).substr(0, get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        prev = graph_->kmer_to_node(rev_seq);
        return prev != npos ? reverse_complement(prev) : prev;
    }
}

void CanonicalDBG::call_nodes(const std::function<void(node_index)> &callback,
                              const std::function<bool()> &stop_early) const {
    graph_->call_nodes(
        [&](node_index i) {
            callback(i);
            if (!stop_early()) {
                node_index j = reverse_complement(i);
                if (j != i)
                    callback(j);
            }
        },
        stop_early
    );
}

void CanonicalDBG
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    graph_->call_kmers([&](node_index i, const std::string &seq) {
        callback(i, seq);
        node_index j = reverse_complement(i);
        if (j != i) {
            std::string rseq(seq);
            ::reverse_complement(rseq.begin(), rseq.end());
            callback(j, rseq);
        }
    });
}

bool CanonicalDBG::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const CanonicalDBG*>(&other)) {
        return *this == dynamic_cast<const CanonicalDBG&>(other);
    }

    return DeBruijnGraph::operator==(other);
}

DeBruijnGraph::node_index CanonicalDBG::reverse_complement(node_index node) const {
    assert(node);
    assert(node <= offset_ * 2);

    if (node > offset_) {
        // we know that this node is definitely not present in the base graph

        if (!k_odd_)
            is_palindrome_cache_.Put(node - offset_, false);

        return node - offset_;

    } else if (k_odd_) {
        // if k is odd, then we know that the reverse complement doesn't exist,
        // so we apply the offset
        return node + offset_;

    }

    if (auto fetch = is_palindrome_cache_.TryGet(node))
        return *fetch ? node : node + offset_;

    std::string seq = graph_->get_node_sequence(node);
    std::string rev_seq = seq;
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    bool palindrome = (rev_seq == seq);

    assert(palindrome || graph_->kmer_to_node(rev_seq) == npos);

    is_palindrome_cache_.Put(node, palindrome);
    return palindrome ? node : node + offset_;
}

void CanonicalDBG::reverse_complement(std::string &seq,
                                      std::vector<node_index> &path) const {
    ::reverse_complement(seq.begin(), seq.end());

    std::vector<node_index> rev_path(path.size());
    std::transform(path.begin(), path.end(), rev_path.rbegin(), [&](node_index i) {
        return i ? reverse_complement(i) : i;
    });
    std::swap(path, rev_path);
}



void CanonicalDBG
::call_outgoing_from_rev_comp(node_index node,
                              const std::function<void(node_index, TAlphabet)> &callback) const {
    if (edge_index e = get_rev_comp_boss_next_node(node)) {
        const auto &cached = dynamic_cast<const DBGSuccinct::CachedView&>(*graph_);
        const DBGSuccinct &dbg_succ = cached.get_graph();
        const boss::BOSS &boss = dbg_succ.get_boss();
        //        rshift    rc
        // ATGGCT -> TGGCT* -> *AGCCA
        std::string rev_comp_suffix = cached.get_node_sequence(node).substr(1) + std::string(1, '\0');
        ::reverse_complement(rev_comp_suffix.begin(), rev_comp_suffix.end());
        boss.call_incoming_to_target(boss.bwd(e), boss.get_node_last_value(e),
                                     [&](edge_index incoming_edge) {
            TAlphabet c = cached.get_first_value(incoming_edge);
            rev_comp_suffix[0] = boss.decode(c);
            cached.put_decoded_edge(incoming_edge, rev_comp_suffix);
            auto kmer_index = dbg_succ.boss_to_kmer_index(incoming_edge);
            if (kmer_index != npos)
                callback(kmer_index, c);
        });
    }
}

void CanonicalDBG
::call_incoming_to_rev_comp(node_index node,
                            const std::function<void(node_index, TAlphabet)> &callback) const {
    if (edge_index e = get_rev_comp_boss_prev_node(node)) {
        const auto &cached = dynamic_cast<const DBGSuccinct::CachedView&>(*graph_);
        const DBGSuccinct &dbg_succ = cached.get_graph();
        const boss::BOSS &boss = dbg_succ.get_boss();
        //        lshift    rc
        // AGCCAT -> *AGCCA -> TGGCT*
        std::string rev_comp_prefix = std::string(1, '\0') + cached.get_node_sequence(node).substr(0, get_k() - 1);
        ::reverse_complement(rev_comp_prefix.begin(), rev_comp_prefix.end());
        boss.call_outgoing(e, [&](edge_index adjacent_edge) {
            TAlphabet c = boss.get_W(adjacent_edge) % boss.alph_size;
            rev_comp_prefix.back() = boss.decode(c);
            cached.put_decoded_edge(adjacent_edge, rev_comp_prefix);
            auto kmer_index = dbg_succ.boss_to_kmer_index(adjacent_edge);
            if (kmer_index != npos)
                callback(kmer_index, c);
        });
    }
}

auto CanonicalDBG
::get_rev_comp_boss_next_node(node_index node) const -> edge_index {
    const auto &cached = dynamic_cast<const DBGSuccinct::CachedView&>(*graph_);
    const DBGSuccinct &dbg_succ = cached.get_graph();
    const boss::BOSS &boss = dbg_succ.get_boss();
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
        std::string rev_seq = cached.get_node_sequence(node).substr(1, get_k() - 1);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_next_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss.encode(rev_seq);
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge == edge_2);

        if (end == encoded.end()) {
            ret_val = edge;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        edge_index start = boss.pred_W(dbg_succ.kmer_to_boss_index(node),
                                       cached.complement(encoded.front()));
        boss.call_incoming_to_target(start, cached.complement(encoded.front()),
                                     [&](edge_index adj) {
            if (node_index a = dbg_succ.boss_to_kmer_index(adj))
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
        cached.adjacent_outgoing_nodes(node, [&](node_index next) {
            rev_comp_next_cache_.Put(next, stored_val);
#ifndef NDEBUG
            std::string test_seq = dbg_succ.get_node_sequence(next).substr(1, get_k() - 1);
            ::reverse_complement(test_seq.begin(), test_seq.end());
            auto encoded = boss.encode(test_seq);
            auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
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
        cached.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                parents[boss.encode(c)].first = next;
                edge_index next_edge = dbg_succ.kmer_to_boss_index(next);
                if (boss.get_last(next_edge))
                    start = next_edge;

                ++parents_count;
            }
        });

        if (parents_count) {
            edge_index ret_prev = boss.bwd(ret_val);
            TAlphabet w = boss.get_node_last_value(ret_val);
            boss.call_incoming_to_target(ret_prev, w, [&](edge_index prev_edge) {
                if (boss.get_last(prev_edge)) {
                    TAlphabet s = cached.complement(cached.get_first_value(prev_edge));
                    if (parents[s].first)
                        rev_comp_next_cache_.Put(parents[s].first, std::make_pair(prev_edge, 0));
                }
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = dbg_succ.get_node_sequence(node).substr(1, get_k() - 1);
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss.get_node_str(ret_val));
    } else {
        auto encoded = boss.encode(test_seq);
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
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

auto CanonicalDBG
::get_rev_comp_boss_prev_node(node_index node) const -> edge_index {
    // 8% effective
    const auto &cached = dynamic_cast<const DBGSuccinct::CachedView&>(*graph_);
    const DBGSuccinct &dbg_succ = cached.get_graph();
    const boss::BOSS &boss = dbg_succ.get_boss();
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
        std::string rev_seq = cached.get_node_sequence(node);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_prev_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        rev_seq.pop_back();
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss.encode(rev_seq);
        auto [edge_1, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge_1 == edge_2);
        if (end == encoded.end()) {
            ret_val = edge_1;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        boss.call_outgoing(boss.succ_last(dbg_succ.kmer_to_boss_index(node)),
                           [&](edge_index adj) {
            if (node_index a = dbg_succ.boss_to_kmer_index(adj))
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
        if (dbg_succ.boss_to_kmer_index(ret_val)) {
            boss.call_outgoing(ret_val, [&](edge_index next_edge) {
                if (TAlphabet w = boss.get_W(next_edge) % boss.alph_size) {
                    parents[w].second = boss.fwd(boss.pred_W(ret_val, w), w);
                    ++parents_count;
                }
            });
        }

        if (parents_count) {
            cached.call_incoming_kmers(node, [&](node_index prev, char c) {
                TAlphabet s = cached.complement(boss.encode(c));
                if (parents[s].second)
                    rev_comp_prev_cache_.Put(prev, std::make_pair(parents[s].second, 0));
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = dbg_succ.get_node_sequence(node);
    test_seq.pop_back();
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss.get_node_str(ret_val));
    } else {
        auto encoded = boss.encode(test_seq);
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end());
        assert(chars_unmatched <= static_cast<size_t>(encoded.end() - end));
    }
#endif

    return ret_val;
}

} // namespace graph
} // namespace mtg
