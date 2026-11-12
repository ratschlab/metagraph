#include "canonical_dbg.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "common/logger.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using edge_index = boss::BOSS::edge_index;

inline const DBGSuccinct* get_dbg_succ(const DeBruijnGraph &graph) {
    const DeBruijnGraph *base_graph = &graph;
    while (const auto *wrapper = dynamic_cast<const DBGWrapper<>*>(base_graph)) {
        base_graph = &wrapper->get_graph();
    }

    return dynamic_cast<const DBGSuccinct*>(base_graph);
}

CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph, size_t cache_size)
      : DBGWrapper<DeBruijnGraph>(std::move(graph)),
        cache_size_(cache_size), is_palindrome_cache_(cache_size_),
        offset_(graph_->max_index()),
        k_odd_(graph_->get_k() % 2) {
    assert(graph_->get_mode() == DeBruijnGraph::PRIMARY
                && "Only primary graphs can be wrapped in CanonicalDBG");

    has_sentinel_ = false;
    alphabet_encoder_.fill(graph_->alphabet().size());

    for (size_t i = 0; i < graph_->alphabet().size(); ++i) {
        alphabet_encoder_[graph_->alphabet()[i]] = i;
        if (graph_->alphabet()[i] == boss::BOSS::kSentinel)
            has_sentinel_ = true;
    }

    if (const DBGSuccinct *dbg_succ = get_dbg_succ(*graph_)) {
        has_sentinel_ = !dbg_succ->get_mask();
        fallback_cache_ = std::make_unique<NodeFirstCache>(*dbg_succ, 0);
    }
}

std::shared_ptr<CanonicalDBG>
CanonicalDBG::with_node_first_cache(std::shared_ptr<const DeBruijnGraph> graph) {
    auto canonical = std::make_shared<CanonicalDBG>(std::move(graph));
    if (const DBGSuccinct *dbg_succ = get_dbg_succ(canonical->get_graph()))
        canonical->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));
    return canonical;
}

NodeFirstCache& CanonicalDBG::get_cache() const {
    if (auto *cache = get_extension_threadsafe<NodeFirstCache>())
        return *cache;

    assert(fallback_cache_);

    return *fallback_cache_;
}

std::unique_ptr<NodeStream>
CanonicalDBG::node_stream(std::string_view sequence, bool sequential) const {
    class Stream final : public NodeStream {
      public:
        Stream(const CanonicalDBG &graph, std::string_view sequence)
              : graph_(graph), sequence_(sequence),
                prefix_(graph.graph_->node_stream(sequence, true)) {}

        size_t next(std::span<node_index> out) override {
            size_t n = 0;
            // map until the first mismatch, one node at a time: lookups in
            // the base graph may be slow and the suffix maps them differently
            while (prefix_ && n < out.size()) {
                node_index node;
                if (!prefix_->next({ &node, 1 })) {
                    prefix_.reset();
                } else if (node) {
                    out[n++] = node;
                    num_mapped_++;
                } else {
                    prefix_.reset();
                    suffix_ = map_suffix(sequence_.substr(num_mapped_));
                }
            }
            size_t m = std::min(out.size() - n, suffix_.size() - pos_);
            std::copy_n(suffix_.begin() + pos_, m, out.begin() + n);
            pos_ += m;
            return n + m;
        }

      private:
        const CanonicalDBG &graph_;
        std::string_view sequence_;
        std::unique_ptr<NodeStream> prefix_;
        size_t num_mapped_ = 0;
        std::vector<node_index> suffix_;
        size_t pos_ = 0;

        // Map a sequence whose first k-mer is missing in the base graph
        std::vector<node_index> map_suffix(std::string_view sequence) const {
            const DeBruijnGraph &base = *graph_.graph_;

            std::string rev_seq(sequence);
            ::reverse_complement(rev_seq.begin(), rev_seq.end());
            // map the reverse-complement
            std::vector<node_index> rev_path = graph::map_to_nodes_sequentially(base, rev_seq);

            // map the forward
            std::vector<node_index> path;
            const auto *dbg_succ = get_dbg_succ(base);
            if (dbg_succ && graph_.k_odd_) {
                // if it's a boss table with odd k (without palindromic k-mers),
                // we can skip k-mers that have been found in the rev-compl sequence
                const auto &boss = dbg_succ->get_boss();
                path.reserve(rev_path.size());
                // the first k-mer is missing and we skip it
                path.push_back(npos);
                auto it = rev_path.rbegin() + 1;
                auto is_missing = get_missing_kmer_skipper(dbg_succ->get_bloom_filter(),
                                                           sequence.substr(1));
                boss.map_to_edges(sequence.substr(1),
                    [&](boss::BOSS::edge_index edge) {
                        path.push_back(dbg_succ->validate_edge(edge));
                        ++it;
                    },
                    []() { return false; },
                    [&]() { return is_missing() || *it; }
                );

                assert(it == rev_path.rend());

            } else {
                path = graph::map_to_nodes_sequentially(base, sequence);
            }

            assert(path.size() == rev_path.size());

            auto it = rev_path.rbegin();
            for (auto jt = path.begin(); jt != path.end(); ++jt, ++it) {
                if (*jt == npos && *it != npos)
                    *jt = *it + graph_.offset_;
            }
            return path;
        }
    };

    auto stream = std::make_unique<Stream>(*this, sequence);
    if (sequential)
        return stream;

    return map_node_stream(std::move(stream),
                           [this](node_index node) { return get_base_node(node); });
}

template <class Callback>
void CanonicalDBG::call_outgoing_kmers(node_index node,
                                       const std::string &spelling_hint,
                                       const Callback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);

    if (node > offset_) {
        std::string rc_hint = spelling_hint;
        ::reverse_complement(rc_hint.begin(), rc_hint.end());
        // type-erased, so that the two templates do not instantiate each other forever
        const std::function<void(node_index, char)> rc_callback = [&](node_index next, char c) {
            c = complement(c);
            assert(c == boss::BOSS::kSentinel
                || traverse(node, c) == reverse_complement(next));
            callback(reverse_complement(next), c);
        };
        call_incoming_kmers(node - offset_, rc_hint, rc_callback);
        return;
    }

    // includes `$` for DBGSuccinct
    const auto &alphabet = graph_->alphabet();

    SmallVector<node_index> children(alphabet.size(), npos);
    // "- has_sentinel_" because there can't be a dummy sink with another non-dummy edge
    size_t max_num_edges_left = children.size() - has_sentinel_;

    for (auto [next, c] : graph_->outgoing_kmers(node)) {
        assert(has_sentinel_ || c != boss::BOSS::kSentinel);
        assert(c == boss::BOSS::kSentinel || traverse(node, c) == next);

        // don't call dummy sink, it will be handled in the adjacent strand below
        if (c != boss::BOSS::kSentinel) {
            callback(next, c);
            --max_num_edges_left;
        }

        children[alphabet_encoder_[c]] = next;
    }

    if (!max_num_edges_left)
        return;

    // jump to the other strand and check incoming k-mers and jump back
    adjacent_outgoing_rc_strand(node, spelling_hint, [&](node_index next, char c) {
        assert(next <= offset_);

        // transform `c`, `next` will be transformed later
        c = complement(c);
        auto s = alphabet_encoder_[c];

        if (children[s] != npos && c != boss::BOSS::kSentinel) {
            // `next` must be a palindrome
            is_palindrome_cache_.insert(next, true);
            if (k_odd_) { // palindromes exist only for even k
                logger->error(
                    "Forward traversal: Primary graph contains both forward and reverse complement: {} {} -> {} {}\t{} {}",
                    node, get_node_sequence(node),
                    children[s], get_node_sequence(children[s]),
                    next, get_node_sequence(next));
                throw std::runtime_error("");
            }
            return;
        }

        is_palindrome_cache_.insert(next, false);

        // transform `next` to the adj strand
        next = reverse_complement(next);

        assert(c == get_node_sequence(next).back());
        // traverse works only for real edges
        assert(c == boss::BOSS::kSentinel || traverse(node, c) == next);

        if (c != boss::BOSS::kSentinel) {
            callback(next, c);
            children[s] = next;
            --max_num_edges_left;
        }
    });

    if (has_sentinel_ && children[alphabet_encoder_[boss::BOSS::kSentinel]]
            && max_num_edges_left + 1 == alphabet.size()) {
        callback(children[alphabet_encoder_[boss::BOSS::kSentinel]],
                 boss::BOSS::kSentinel);
    }
}

template <class Callback>
void CanonicalDBG::call_incoming_kmers(node_index node,
                                       const std::string &spelling_hint,
                                       const Callback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);

    if (node > offset_) {
        std::string rc_hint = spelling_hint;
        ::reverse_complement(rc_hint.begin(), rc_hint.end());
        // type-erased, so that the two templates do not instantiate each other forever
        const std::function<void(node_index, char)> rc_callback = [&](node_index prev, char c) {
            c = complement(c);
            assert(c == boss::BOSS::kSentinel
                || traverse_back(node, c) == reverse_complement(prev));
            callback(reverse_complement(prev), c);
        };
        call_outgoing_kmers(node - offset_, rc_hint, rc_callback);
        return;
    }

    // includes `$` for DBGSuccinct
    const auto &alphabet = graph_->alphabet();

    SmallVector<node_index> parents(alphabet.size(), npos);
    // "- has_sentinel_" because there can't be a dummy sink with another non-dummy edge
    size_t max_num_edges_left = parents.size() - has_sentinel_;
    for (auto [prev, c] : get_dbg_succ(*graph_) ? get_cache().incoming_kmers(node)
                                                : graph_->incoming_kmers(node)) {
        assert(has_sentinel_ || c != boss::BOSS::kSentinel);
        assert(c == boss::BOSS::kSentinel || traverse_back(node, c) == prev);

        // don't call dummy source, it will be handled in the adjacent strand below
        if (c != boss::BOSS::kSentinel) {
            callback(prev, c);
            --max_num_edges_left;
        }

        parents[alphabet_encoder_[c]] = prev;
    }

    if (!max_num_edges_left)
        return;

    adjacent_incoming_rc_strand(node, spelling_hint, [&](node_index prev, char c) {
        assert(prev <= offset_);

        // transform `c`, `prev` will be transformed later
        c = complement(c);
        auto s = alphabet_encoder_[c];

        if (parents[s] != npos && c != boss::BOSS::kSentinel) {
            // `prev` must be a palindrome
            is_palindrome_cache_.insert(prev, true);
            if (k_odd_) {
                logger->error(
                    "Backward traversal: Primary graph contains both forward and reverse complement: {} {} <- {} {}\t{} {}",
                    node, get_node_sequence(node),
                    parents[s], get_node_sequence(parents[s]),
                    prev, get_node_sequence(prev));
                throw std::runtime_error("");
            }
            return;
        }

        is_palindrome_cache_.insert(prev, false);

        // transform `prev` to the adj strand
        prev = reverse_complement(prev);

        assert(c == get_node_sequence(prev)[0]);
        assert(c == boss::BOSS::kSentinel || traverse_back(node, c) == prev);

        if (c != boss::BOSS::kSentinel) {
            callback(prev, c);
            parents[s] = prev;
            --max_num_edges_left;
        }
    });

    if (has_sentinel_ && parents[alphabet_encoder_[boss::BOSS::kSentinel]]
            && max_num_edges_left + 1 == alphabet.size()) {
        callback(parents[alphabet_encoder_[boss::BOSS::kSentinel]],
                 boss::BOSS::kSentinel);
    }
}

DeBruijnGraph::Edges CanonicalDBG::outgoing_kmers(node_index node,
                                                  const std::string &spelling_hint) const {
    Edges edges;
    call_outgoing_kmers(node, spelling_hint, [&](node_index next, char c) {
        edges.emplace_back(next, c);
    });
    return edges;
}

DeBruijnGraph::Edges CanonicalDBG::incoming_kmers(node_index node,
                                                  const std::string &spelling_hint) const {
    Edges edges;
    call_incoming_kmers(node, spelling_hint, [&](node_index prev, char c) {
        edges.emplace_back(prev, c);
    });
    return edges;
}

template <class Callback>
void CanonicalDBG::call_incoming_nodes(node_index node, const Callback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);

    if (node > offset_) {
        const std::function<void(node_index)> rc_callback = [&](node_index prev) {
            callback(reverse_complement(prev));
        };
        call_outgoing_nodes(node - offset_, rc_callback);
    } else {
        call_incoming_kmers(node, get_node_sequence(node),
                            [&](node_index prev, char) { callback(prev); });
    }
}

template <class Callback>
void CanonicalDBG::call_outgoing_nodes(node_index node, const Callback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);

    if (node > offset_) {
        const std::function<void(node_index)> rc_callback = [&](node_index next) {
            callback(reverse_complement(next));
        };
        call_incoming_nodes(node - offset_, rc_callback);
    } else {
        call_outgoing_kmers(node, get_node_sequence(node),
                            [&](node_index next, char) { callback(next); });
    }
}

DeBruijnGraph::Nodes CanonicalDBG::adjacent_incoming_nodes(node_index node) const {
    Nodes nodes;
    call_incoming_nodes(node, [&](node_index prev) { nodes.push_back(prev); });
    return nodes;
}

DeBruijnGraph::Nodes CanonicalDBG::adjacent_outgoing_nodes(node_index node) const {
    Nodes nodes;
    call_outgoing_nodes(node, [&](node_index next) { nodes.push_back(next); });
    return nodes;
}

bool CanonicalDBG::has_multiple_outgoing(node_index node) const {
    // stop at the second node: this can skip the lookup on the other strand
    bool found = false;
    try {
        call_outgoing_nodes(node, [&](node_index) {
            if (found)
                throw early_term();

            found = true;
        });
    } catch (const early_term&) {
        return true;
    }

    return false;
}

bool CanonicalDBG::has_single_incoming(node_index node) const {
    size_t count = 0;
    try {
        call_incoming_nodes(node, [&](node_index) {
            if (++count > 1)
                throw early_term();
        });
    } catch (const early_term&) {}

    return count == 1;
}

size_t CanonicalDBG::outdegree(node_index node) const {
    // TODO: better implementation
    return adjacent_outgoing_nodes(node).size();
}

size_t CanonicalDBG::indegree(node_index node) const {
    // TODO: better implementation
    return adjacent_incoming_nodes(node).size();
}

void CanonicalDBG::call_sequences(const CallPath &callback,
                                  size_t num_threads,
                                  bool kmers_in_single_form,
                                  bool verbose) const {
    if (kmers_in_single_form) {
        graph_->call_sequences(callback, num_threads, true, verbose);
    } else {
        // TODO: port over implementation from DBGSuccinct to DeBruijnGraph
        DeBruijnGraph::call_sequences(callback, num_threads, false, verbose);
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
        // TODO: use the cache here?
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
                              const std::function<bool()> &stop_early,
                              size_t num_threads,
                              size_t batch_size) const {
    graph_->call_nodes(
        [&](node_index i) {
            callback(i);
            if (!stop_early()) {
                node_index j = reverse_complement(i);
                if (j != i)
                    callback(j);
            }
        },
        stop_early,
        num_threads,
        batch_size
    );
}

void CanonicalDBG
::call_kmers(const std::function<void(node_index, const std::string&)> &callback,
             const std::function<bool()> &stop_early) const {
    graph_->call_kmers([&](node_index i, const std::string &seq) {
        callback(i, seq);
        if (stop_early())
            return;

        node_index j = reverse_complement(i);
        if (j != i) {
            std::string rseq(seq);
            ::reverse_complement(rseq.begin(), rseq.end());
            callback(j, rseq);
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
            is_palindrome_cache_.insert(node - offset_, false);

        return node - offset_;

    } else if (k_odd_) {
        // if k is odd, then we know that the reverse complement doesn't exist,
        // so we apply the offset
        return node + offset_;

    }

    assert(!k_odd_ && node <= offset_);

    if (auto fetch = is_palindrome_cache_.find(node))
        return *fetch ? node : node + offset_;

    std::string seq = get_node_sequence(node);
    std::string rev_seq = seq;
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    bool palindrome = (rev_seq == seq);

    assert(palindrome || graph_->kmer_to_node(rev_seq) == npos);

    is_palindrome_cache_.insert(node, palindrome);

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

/**
 * find parents of node by searching for children of its reverse complement
 * e.g., node = AGCCAT. Find TAGCCA and AAGCCA by looking for TGGCTA and TGGCTT.
 *  TAGCCA                    TGGCTA
 *        \                  /
 *         AGCCAT  ->  ATGGCT
 *        /                  \
 *  AAGCCA                    TGGCTT
 */

// for each n, check for TGGCTn. If found, define and the index for
// rc(n)AGCCA as index(TGGCTn) + offset_
void CanonicalDBG
::adjacent_incoming_rc_strand(node_index node,
                              const std::string &spelling_hint,
                              const std::function<void(node_index, char)> &callback) const {
    //        lshift    rc
    // AGCCAT -> *AGCCA -> TGGCT*
    if (const auto *dbg_succ_ = get_dbg_succ(*graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //-> AGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //-> TCAAGCAGAAGACGGCATACGAGATCCTCT
        const boss::BOSS &boss = dbg_succ_->get_boss();

        boss::BOSS::edge_index rc_edge = get_cache().get_prefix_rc(node, spelling_hint);

        if (!rc_edge)
            return;

        boss.call_outgoing(rc_edge, [&](boss::BOSS::edge_index adjacent_edge) {
            assert(dbg_succ_);
            node_index prev = adjacent_edge;
            if (!dbg_succ_->in_graph(prev))
                return;

            char c = boss.decode(boss.get_W(adjacent_edge) % boss.alph_size);
            // make sure that the following type of backward traversal doesn't happen
            // TCCTGCGCTTCGTACATATTCCCGCCGCAC$ -> rc($TCCTGCGCTTCGTACATATTCCCGCCGCAC)
            if (spelling_hint.back() == boss::BOSS::kSentinel && c == boss::BOSS::kSentinel)
                return;

            callback(prev, c);
        });
    } else {
        // Do the checks by directly mapping the sequences of the desired k-mers.
        // For non-DBGSuccinct graphs, this should be fast enough.
        std::string rev_seq = graph_->get_node_sequence(node).substr(0, graph_->get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        rev_seq.push_back('\0');

        for (char c : graph_->alphabet()) {
            rev_seq.back() = c;
            node_index prev = graph_->kmer_to_node(rev_seq);
            if (prev != npos)
                callback(prev, c);
        }
    }
}

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

// for each n, check for nAGCCA. If found, define the index for
// TGGCTrc(n) as index(nAGCCA) + offset_
void CanonicalDBG
::adjacent_outgoing_rc_strand(node_index node,
                              const std::string &spelling_hint,
                              const std::function<void(node_index, char)> &callback) const {
    //        rshift    rc
    // ATGGCT -> TGGCT* -> *AGCCA
    if (const auto *dbg_succ_ = get_dbg_succ(*graph_)) {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  GAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  CTCAAGCAGAAGACGGCATACGAGATCCTC

        auto &cache = get_cache();

        boss::BOSS::edge_index rc_edge = cache.get_suffix_rc(node, spelling_hint);

        if (!rc_edge)
            return;

        cache.call_incoming_edges(rc_edge,
            [&](edge_index prev_edge) {
                node_index prev = prev_edge;
                if (!dbg_succ_->in_graph(prev))
                    return;

                char c = cache.get_first_char(prev_edge, rc_edge);
                // meant to make sure that the following type of traversal doesn't happen
                // $TGTGCGGCGGGAATATGTACGAAGCGCAGG -> rc(TGTGCGGCGGGAATATGTACGAAGCGCAGG$)
                if (spelling_hint[0] == boss::BOSS::kSentinel && c == boss::BOSS::kSentinel)
                    return;

                callback(prev, c);
            }
        );
    } else {
        // Do the checks by directly mapping the sequences of the desired k-mers.
        // For non-DBGSuccinct graphs, this should be fast enough.
        std::string rev_seq = graph_->get_node_sequence(node).substr(1) + std::string(1, '\0');
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        assert(rev_seq[0] == '\0');

        for (char c : graph_->alphabet()) {
            rev_seq[0] = c;
            node_index prev = graph_->kmer_to_node(rev_seq);
            if (prev != npos)
                callback(prev, c);
        }
    }
}

} // namespace graph
} // namespace mtg
