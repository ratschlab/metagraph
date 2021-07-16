#include "canonical_dbg.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/dbg_succinct_cached.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;

inline const DBGSuccinct* get_dbg_succ(const DeBruijnGraph &graph) {
    return dynamic_cast<const DBGSuccinct*>(&graph.get_base_graph());
}


CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph, size_t cache_size)
      : const_graph_ptr_(graph),
        offset_(graph_.max_index()),
        k_odd_(graph_.get_k() % 2),
        has_sentinel_(false),
        alphabet_encoder_({ graph_.alphabet().size() }),
        child_node_cache_(cache_size),
        parent_node_cache_(cache_size),
        is_palindrome_cache_(k_odd_ ? 0 : cache_size) {
    if (graph->get_mode() != DeBruijnGraph::PRIMARY)
        throw std::runtime_error("Only primary graphs can be wrapped in CanonicalDBG");

    for (size_t i = 0; i < graph_.alphabet().size(); ++i) {
        alphabet_encoder_[graph_.alphabet()[i]] = i;
        if (graph_.alphabet()[i] == boss::BOSS::kSentinel)
            has_sentinel_ = true;
    }
}

CanonicalDBG::CanonicalDBG(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size)
      : CanonicalDBG(std::dynamic_pointer_cast<const DeBruijnGraph>(graph), cache_size) {
    graph_ptr_ = graph;
}

CanonicalDBG::CanonicalDBG(const DeBruijnGraph &graph, size_t cache_size)
      : CanonicalDBG(std::shared_ptr<const DeBruijnGraph>(&graph, [](const auto*) {}),
                     cache_size) {}

CanonicalDBG::CanonicalDBG(DeBruijnGraph &graph, size_t cache_size)
      : CanonicalDBG(std::shared_ptr<DeBruijnGraph>(&graph, [](const auto*) {}),
                     cache_size) {}

uint64_t CanonicalDBG::num_nodes() const {
    logger->trace("Number of nodes may be overestimated if k is even or reverse complements are present in the graph");
    return graph_.num_nodes() * 2;
}

void CanonicalDBG
::add_sequence(std::string_view sequence,
               const std::function<void(node_index)> &on_insertion) {
    if (!graph_ptr_)
        throw std::runtime_error("add_sequence only supported for non-const graphs.");

    graph_ptr_->add_sequence(sequence, on_insertion);
    offset_ = graph_.max_index();
    child_node_cache_.Clear();
    parent_node_cache_.Clear();
}

void CanonicalDBG
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    if (sequence.size() < get_k())
        return;

    std::vector<node_index> path;
    path.reserve(sequence.size() - get_k() + 1);

    // map until the first mismatch
    bool stop = false;
    graph_.map_to_nodes_sequentially(sequence,
        [&](node_index node) {
            if (node) {
                path.push_back(node);
            } else {
                stop = true;
            }
        },
        [&]() { return stop; }
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
    std::vector<node_index> rev_path = map_sequence_to_nodes(graph_, rev_seq);

    // map the forward
    const auto *dbg_succ = get_dbg_succ(graph_);

    if (dbg_succ && get_k() % 2) {
        const auto *cached = dynamic_cast<const DBGSuccinctCached*>(&graph_);

        // if it's a boss table with odd k (without palindromic k-mers),
        // we can skip k-mers that have been found in the rev-compl sequence
        const auto &boss = dbg_succ->get_boss();
        // the initial forward mapping stopped on this k-mer,
        // hence it's missing and we skip it
        path.push_back(npos);
        auto it = rev_path.rbegin() + 1;
        auto is_missing = get_missing_kmer_skipper(dbg_succ->get_bloom_filter(),
                                                   sequence.substr(1));
        auto jt = sequence.begin() + 1;
        boss.map_to_edges(sequence.substr(1),
            [&](boss::BOSS::edge_index edge) {
                path.push_back(dbg_succ->boss_to_kmer_index(edge));
                if (cached && path.back())
                    cached->put_decoded_node(path.back(), std::string(jt, jt + get_k()));

                ++it;
                ++jt;
            },
            []() { return false; },
            [&]() {
                if (is_missing() || *it) {
                    path.push_back(npos);
                    ++it;
                    ++jt;
                    return true;
                } else {
                    return false;
                }
            }
        );

        assert(it == rev_path.rend());

    } else {
        path = map_sequence_to_nodes(graph_, sequence);
    }

    assert(path.size() == rev_path.size());

    auto it = rev_path.rbegin();
    for (auto jt = path.begin(); jt != path.end(); ++jt, ++it) {
        if (terminate())
            return;

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
        callback(get_base_node(i));
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

    const auto &alphabet = graph_.alphabet();

    //        rshift    rc
    // ATGGCT -> TGGCT* -> *AGCCA
    std::string rev_seq = get_node_sequence(node).substr(1) + std::string(1, '\0');
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    assert(rev_seq[0] == '\0');

    // for each n, check for nAGCCA. If found, define and store the index for
    // TGGCTrc(n) as index(nAGCCA) + offset_
    if (const auto *dbg_succ = get_dbg_succ(graph_)) {
        const auto *cached = dynamic_cast<const DBGSuccinctCached*>(&graph_);
        const auto &boss = dbg_succ->get_boss();
        dbg_succ->call_nodes_with_suffix_matching_longest_prefix(
            std::string_view(&rev_seq[1], get_k() - 1),
            [&](node_index next, uint64_t /* match length */) {
                auto edge = dbg_succ->kmer_to_boss_index(next);
                boss::BOSS::TAlphabet c = cached
                    ? cached->get_first_value(edge)
                    : boss.get_minus_k_value(edge, get_k() - 2).first;
                if (c == boss::BOSS::kSentinelCode)
                    return;

                rev_seq[0] = boss.decode(c);
                c = kmer::KmerExtractorBOSS::complement(c);

                if (children[c] != npos) {
                    if (!k_odd_) {
                        is_palindrome_cache_.Put(next, true);
                    } else {
                        throw std::runtime_error(fmt::format(
                            "Primary graph contains both forward and reverse complement: {} {} -> {} {}\t{} {}",
                            node, graph_.get_node_sequence(node),
                            children[c], graph_.get_node_sequence(children[c]),
                            next, graph_.get_node_sequence(next)
                        ));
                    }

                } else {
                    children[c] = next + offset_;
                    if (cached)
                        cached->put_decoded_node(next, rev_seq);
                }
            },
            get_k() - 1
        );

    } else {
        for (size_t c = 0; c < alphabet.size(); ++c) {
            if (alphabet[c] != boss::BOSS::kSentinel && children[c] == npos) {
                rev_seq[0] = complement(alphabet[c]);
                node_index next = graph_.kmer_to_node(rev_seq);
                if (next != npos)
                    children[c] = next + offset_;
            }
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

    const auto &alphabet = graph_.alphabet();

    std::vector<node_index> children;
    if (auto fetch = child_node_cache_.TryGet(node)) {
        children = std::move(*fetch);

    } else {
        children.resize(alphabet.size(), npos);
        size_t max_num_edges_left = children.size() - has_sentinel_;

        graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                children[alphabet_encoder_[c]] = next;
                --max_num_edges_left;
            }
        });

        if (max_num_edges_left)
            append_next_rc_nodes(node, children);

        child_node_cache_.Put(node, children);
    }

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

    const auto &alphabet = graph_.alphabet();

    //        lshift    rc
    // AGCCAT -> *AGCCA -> TGGCT*
    std::string rev_seq = std::string(1, '\0') + get_node_sequence(node).substr(0, get_k() - 1);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    assert(rev_seq.back() == '\0');

    // for each n, check for TGGCTn. If found, define and store the index for
    // rc(n)AGCCA as index(TGGCTn) + offset_
    if (const auto *dbg_succ = get_dbg_succ(graph_)) {
        // Find the BOSS node TGGCT and iterate through all of its outdoing edges.
        // Then, convert the edge indices to get the DBGSuccinct node indices
        const auto *cached = dynamic_cast<const DBGSuccinctCached*>(&graph_);
        const auto &boss = dbg_succ->get_boss();

        auto encoded = boss.encode(std::string_view(rev_seq.data(), get_k() - 1));
        auto [edge, edge_2, end] = boss.index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge == edge_2);

        if (end == encoded.end()) {
            boss.call_outgoing(edge, [&](auto adjacent_edge) {
                node_index prev = dbg_succ->boss_to_kmer_index(adjacent_edge);
                if (prev) {
                    boss::BOSS::TAlphabet c = boss.get_W(adjacent_edge) % boss.alph_size;
                    if (c == boss::BOSS::kSentinelCode)
                        return;

                    rev_seq.back() = boss.decode(c);
                    c = kmer::KmerExtractorBOSS::complement(c);

                    if (parents[c] != npos) {
                        if (!k_odd_) {
                            is_palindrome_cache_.Put(prev, true);
                        } else {
                            throw std::runtime_error(fmt::format(
                                "Primary graph contains both forward and reverse complement: {} {} -> {} {}\t{} {}",
                                node, graph_.get_node_sequence(node),
                                parents[c], graph_.get_node_sequence(parents[c]),
                                prev, graph_.get_node_sequence(prev)
                            ));
                        }

                    } else {
                        parents[c] = prev + offset_;
                        if (cached)
                            cached->put_decoded_node(prev, rev_seq);
                    }
                }
            });
        }

    } else {
        for (size_t c = 0; c < alphabet.size(); ++c) {
            if (alphabet[c] != boss::BOSS::kSentinel && parents[c] == npos) {
                rev_seq.back() = complement(alphabet[c]);
                node_index prev = graph_.kmer_to_node(rev_seq);
                if (prev != npos)
                    parents[c] = prev + offset_;
            }
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

    const auto &alphabet = graph_.alphabet();

    std::vector<node_index> parents;
    if (auto fetch = parent_node_cache_.TryGet(node)) {
        parents = std::move(*fetch);

    } else {
        parents.resize(alphabet.size(), npos);
        size_t max_num_edges_left = parents.size() - has_sentinel_;

        graph_.call_incoming_kmers(node, [&](node_index prev, char c) {
            if (c != boss::BOSS::kSentinel) {
                parents[alphabet_encoder_[c]] = prev;
                --max_num_edges_left;
            }
        });

        if (max_num_edges_left)
            append_prev_rc_nodes(node, parents);

        parent_node_cache_.Put(node, parents);
    }

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
        graph_.call_sequences(callback, num_threads, false);
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

    std::string seq = graph_.get_node_sequence(node);

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
        node_index next = graph_.traverse(node, next_char);
        if (next != npos)
            return next;

        std::string rev_seq = get_node_sequence(node).substr(1) + next_char;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        next = graph_.kmer_to_node(rev_seq);
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
        node_index prev = graph_.traverse_back(node, prev_char);
        if (prev != npos)
            return prev;

        std::string rev_seq = std::string(1, prev_char)
            + get_node_sequence(node).substr(0, get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        prev = graph_.kmer_to_node(rev_seq);
        return prev != npos ? reverse_complement(prev) : prev;
    }
}

void CanonicalDBG::call_nodes(const std::function<void(node_index)> &callback,
                              const std::function<bool()> &stop_early) const {
    graph_.call_nodes([&](node_index i) {
                          callback(i);
                          if (!stop_early()) {
                              node_index j = reverse_complement(i);
                              if (j != i)
                                  callback(j);
                          }
                      }, stop_early);
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

    std::string seq = graph_.get_node_sequence(node);
    std::string rev_seq = seq;
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    bool palindrome = (rev_seq == seq);

    assert(palindrome || graph_.kmer_to_node(rev_seq) == npos);

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

auto CanonicalDBG::get_base_path(const std::vector<node_index> &path,
                                 const std::string &sequence) const
        -> std::pair<std::vector<node_index>, bool /* is reversed */> {
    auto ret_val = graph_.get_base_path(path, sequence);

    if (ret_val.first.size() && get_base_node(ret_val.first[0]) != ret_val.first[0]) {
        ret_val.second = !ret_val.second;
        std::reverse(ret_val.first.begin(), ret_val.first.end());
    }

    std::transform(ret_val.first.begin(), ret_val.first.end(), ret_val.first.begin(),
                   [this](node_index node) { return get_base_node(node); });

    return ret_val;
}

} // namespace graph
} // namespace mtg
