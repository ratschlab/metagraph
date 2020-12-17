#include "canonical_dbg.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/dbg_succinct_range.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;


// If the graph has even order and is primary, then only is_palindrome_cache_
// is required. rev_comp_cache_ is only required if the underlying graph is not
// primary.
CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph,
                           bool primary,
                           size_t cache_size)
      : const_graph_ptr_(graph),
        offset_(graph_.max_index()),
        alphabet_encoder_({ graph_.alphabet().size() }),
        primary_(primary && !graph->is_canonical_mode()
                     && !dynamic_cast<const DBGSuccinctRange*>(graph.get())),
        child_node_cache_(cache_size),
        parent_node_cache_(cache_size),
        rev_comp_cache_(primary_ ? 0 : cache_size),
        is_palindrome_cache_(graph_.get_k() % 2 || !primary_ ? 0 : cache_size) {
    for (size_t i = 0; i < graph_.alphabet().size(); ++i) {
        alphabet_encoder_[graph_.alphabet()[i]] = i;
    }
}

CanonicalDBG::CanonicalDBG(std::shared_ptr<DeBruijnGraph> graph,
                           bool primary,
                           size_t cache_size)
      : CanonicalDBG(std::dynamic_pointer_cast<const DeBruijnGraph>(graph),
                     primary, cache_size) { graph_ptr_ = graph; }

CanonicalDBG::CanonicalDBG(const DeBruijnGraph &graph, bool primary, size_t cache_size)
      : CanonicalDBG(std::shared_ptr<const DeBruijnGraph>(&graph, [](const auto*) {}),
                     primary,
                     cache_size) {}

CanonicalDBG::CanonicalDBG(DeBruijnGraph &graph, bool primary, size_t cache_size)
      : CanonicalDBG(std::shared_ptr<DeBruijnGraph>(&graph, [](const auto*) {}),
                     primary,
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
    rev_comp_cache_.Clear();
}


void CanonicalDBG
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    std::vector<node_index> path = map_sequence_to_nodes(graph_, sequence);
    auto first_not_found = graph_.is_canonical_mode()
        ? path.end()
        : std::find_if(path.begin(), path.end(),
                       [&](node_index node) {
                           return graph_.get_node_length(node) < graph_.get_k();
                       });

    for (auto jt = path.begin(); jt != first_not_found; ++jt) {
        if (terminate())
            return;

        callback(*jt);
    }

    if (first_not_found == path.end())
        return;

    size_t start = first_not_found - path.begin();

    std::string rev_seq(sequence.begin() + start, sequence.end());
    std::vector<node_index> rev_path(path.begin() + start, path.end());
    reverse_complement_seq_path(graph_, rev_seq, rev_path);
    std::reverse(rev_path.begin(), rev_path.end());

    if (path.size() > rev_path.size() + start) {
        rev_path.resize(path.size() - start);
    } else if (path.size() < rev_path.size() + start) {
        path.resize(rev_path.size() + start);
    }

    assert(rev_path.size() + start == path.size());

    auto it = rev_path.begin();
    for (auto jt = first_not_found; jt != path.end(); ++jt) {
        assert(it != rev_path.end());

        if (terminate())
            return;

        if (*it && *jt) {
            size_t fwd_length = graph_.get_node_length(*jt);
            size_t rev_length = graph_.get_node_length(*it);
            if (fwd_length > rev_length) {
                *it = DeBruijnGraph::npos;
            } else if (fwd_length < rev_length) {
                *jt = DeBruijnGraph::npos;
            }
        }

        if (*jt != DeBruijnGraph::npos) {
            if (!primary_)
                rev_comp_cache_.Put(*jt, *it != DeBruijnGraph::npos ? *it : *jt + offset_);

            callback(*jt);

        } else if (*it != DeBruijnGraph::npos) {
            if (!primary_)
                rev_comp_cache_.Put(*it, *it + offset_);

            callback(*it + offset_);
        } else {
            callback(DeBruijnGraph::npos);
        }
        ++it;
    }
}

void CanonicalDBG::map_to_nodes(std::string_view sequence,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    if (graph_.is_canonical_mode()) {
        graph_.map_to_nodes(sequence, callback, terminate);
    } else {
        map_to_nodes_sequentially(sequence, [&](node_index i) {
            callback(i != DeBruijnGraph::npos ? get_base_node(i) : i);
        }, terminate);
    }
}

void CanonicalDBG::get_kmers_from_suffix(node_index node,
                                         std::vector<node_index> &children) const {
    std::string rev_seq = get_node_sequence(node).substr(1)
        + std::string(1, boss::BOSS::kSentinel);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    const auto &alphabet = graph_.alphabet();

    std::shared_ptr<const DBGSuccinctRange> range_graph;
    const auto *range_cast = dynamic_cast<const DBGSuccinctRange*>(&graph_);
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph_);

    if (range_cast) {
        range_graph = std::shared_ptr<const DBGSuccinctRange>(range_graph, range_cast);
    } else if (dbg_succ) {
        range_graph = std::make_shared<const DBGSuccinctRange>(*dbg_succ);
    }

    if (range_graph) {
        node_index next_range = range_graph->toggle_node_sink_source(
            range_graph->kmer_to_node(rev_seq)
        );

        if (next_range) {
            assert(range_graph->is_sink(next_range));
            assert(range_graph->get_offset(next_range) == 1);
            range_graph->call_incoming_kmers(next_range, [&](node_index next, char c) {
                assert(range_graph->get_offset(next) == 0);
                ::reverse_complement(&c, &c + 1);
                auto i = alphabet_encoder_[c];
                if (children[i] != DeBruijnGraph::npos || c == boss::BOSS::kSentinel)
                    return;

                assert(!graph_.traverse(node, c));

                if (!primary_)
                    rev_comp_cache_.Put(next, next + offset_);

                children[i] = next + offset_;
            });
        }

        return;
    }

    node_index next;
    for (size_t i = 0; i < alphabet.size(); ++i) {
        char c = alphabet[i];

        if (children[i] != DeBruijnGraph::npos || c == boss::BOSS::kSentinel)
            continue;

        ::reverse_complement(&c, &c + 1);

        rev_seq[0] = c;
        next = graph_.kmer_to_node(rev_seq);
        if (next != DeBruijnGraph::npos) {
            assert(!graph_.traverse(node, alphabet[i]));
            if (!primary_)
                rev_comp_cache_.Put(next, next + offset_);

            children[i] = next + offset_;
        }
    }
}

void CanonicalDBG
::call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (node > offset_) {
        call_incoming_kmers(node - offset_, [&](node_index next, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(next), c);
            assert(traverse(node, c) == reverse_complement(next));
        });
        return;
    }

    const auto &alphabet = graph_.alphabet();

    try {
        auto children = child_node_cache_.Get(node);
        for (size_t i = 0; i < alphabet.size(); ++i) {
            if (children[i] != DeBruijnGraph::npos)
                callback(children[i], alphabet[i]);
        }

    } catch (...) {
        size_t max_num_edges_left = alphabet.size();
        std::vector<node_index> children(max_num_edges_left);

        graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                children[alphabet_encoder_[c]] = next;

            --max_num_edges_left;
        });

        if (!graph_.is_canonical_mode() && max_num_edges_left)
            get_kmers_from_suffix(node, children);

        child_node_cache_.Put(node, children);
        for (size_t i = 0; i < children.size(); ++i) {
            if (children[i] != DeBruijnGraph::npos) {
                callback(children[i], alphabet[i]);
                assert(traverse(node, alphabet[i]) == children[i]);
            }
        }
    }
}

void CanonicalDBG::get_kmers_from_prefix(node_index node,
                                         std::vector<node_index> &parents) const {
    std::string rev_seq = std::string(1, boss::BOSS::kSentinel)
        + get_node_sequence(node).substr(0, get_k() - 1);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());
    const auto &alphabet = graph_.alphabet();

    std::shared_ptr<const DBGSuccinctRange> range_graph;
    const auto *range_cast = dynamic_cast<const DBGSuccinctRange*>(&graph_);
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph_);

    if (range_cast) {
        range_graph = std::shared_ptr<const DBGSuccinctRange>(range_graph, range_cast);
    } else if (dbg_succ) {
        range_graph = std::make_shared<const DBGSuccinctRange>(*dbg_succ);
    }

    if (range_graph) {
        node_index prev_range = range_graph->toggle_node_sink_source(
            range_graph->kmer_to_node(rev_seq)
        );

        if (prev_range) {
            assert(!range_graph->is_sink(prev_range));
            assert(range_graph->get_offset(prev_range) == 1);
            range_graph->call_outgoing_kmers(prev_range, [&](node_index prev, char c) {
                assert(range_graph->get_offset(prev) == 0);
                ::reverse_complement(&c, &c + 1);
                auto i = alphabet_encoder_[c];
                if (parents[i] != DeBruijnGraph::npos || c == boss::BOSS::kSentinel)
                    return;

                assert(!graph_.traverse_back(node, c));

                if (!primary_)
                    rev_comp_cache_.Put(prev, prev + offset_);

                parents[i] = prev + offset_;
            });
        }

        return;
    }

    node_index prev;
    for (size_t i = 0; i < alphabet.size(); ++i) {
        char c = alphabet[i];

        if (parents[i] != DeBruijnGraph::npos || c == boss::BOSS::kSentinel)
            continue;

        ::reverse_complement(&c, &c + 1);

        rev_seq.back() = c;
        prev = graph_.kmer_to_node(rev_seq);
        if (prev != DeBruijnGraph::npos) {
            assert(!graph_.traverse_back(node, alphabet[i]));
            if (!primary_)
                rev_comp_cache_.Put(prev, prev + offset_);

            parents[i] = prev + offset_;
        }
    }
}

void CanonicalDBG
::call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (node > offset_) {
        call_outgoing_kmers(node - offset_, [&](node_index prev, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(prev), c);
            assert(traverse_back(node, c) == reverse_complement(prev));
        });
        return;
    }

    const auto &alphabet = graph_.alphabet();

    try {
        auto parents = parent_node_cache_.Get(node);
        for (size_t i = 0; i < alphabet.size(); ++i) {
            if (parents[i] != DeBruijnGraph::npos)
                callback(parents[i], alphabet[i]);
        }

    } catch (...) {
        size_t max_num_edges_left = alphabet.size();
        std::vector<node_index> parents(max_num_edges_left);

        graph_.call_incoming_kmers(node, [&](node_index prev, char c) {
            if (c != boss::BOSS::kSentinel)
                parents[alphabet_encoder_[c]] = prev;

            --max_num_edges_left;
        });

        if (!graph_.is_canonical_mode() && max_num_edges_left)
            get_kmers_from_prefix(node, parents);

        parent_node_cache_.Put(node, parents);
        for (size_t i = 0; i < parents.size(); ++i) {
            if (parents[i] != DeBruijnGraph::npos) {
                callback(parents[i], alphabet[i]);
                assert(traverse_back(node, alphabet[i]) == parents[i]);
            }
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
        ::reverse_complement(&next_char, &next_char + 1);
        node = traverse_back(node - offset_, next_char);
        return node != DeBruijnGraph::npos ? reverse_complement(node) : DeBruijnGraph::npos;
    } else {
        node_index next = graph_.traverse(node, next_char);
        if (next != DeBruijnGraph::npos || graph_.is_canonical_mode())
            return next;

        std::string rev_seq = get_node_sequence(node).substr(1) + next_char;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        next = graph_.kmer_to_node(rev_seq);
        return next != DeBruijnGraph::npos ? reverse_complement(next) : next;
    }
}

DeBruijnGraph::node_index CanonicalDBG::traverse_back(node_index node,
                                                      char prev_char) const {
    assert(node <= offset_ * 2);
    if (node > offset_) {
        ::reverse_complement(&prev_char, &prev_char + 1);
        node = traverse(node - offset_, prev_char);
        return node != DeBruijnGraph::npos ? reverse_complement(node) : DeBruijnGraph::npos;
    } else {
        node_index prev = graph_.traverse_back(node, prev_char);
        if (prev != DeBruijnGraph::npos || graph_.is_canonical_mode())
            return prev;

        std::string rev_seq = std::string(1, prev_char)
            + get_node_sequence(node).substr(0, get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        prev = graph_.kmer_to_node(rev_seq);
        return prev != DeBruijnGraph::npos ? reverse_complement(prev) : prev;
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

        if (!primary_) {
            rev_comp_cache_.Put(node - offset_, node);
        } else if (!(graph_.get_k() & 1)) {
            is_palindrome_cache_.Put(node - offset_, false);
        }

        return node - offset_;

    } else if ((graph_.get_k() & 1) && primary_) {
        // if k is odd and the underlying graph is primary, then we know that
        // the reverse complement doesn't exist, so we apply the offset
        return node + offset_;

    }

    if (!primary_) {
        try {
            return rev_comp_cache_.Get(node);

        } catch (...) {
            std::string rev_seq = graph_.get_node_sequence(node);
            ::reverse_complement(rev_seq.begin(), rev_seq.end());

            node_index rev_node = graph_.kmer_to_node(rev_seq);
            rev_node = rev_node != DeBruijnGraph::npos ? rev_node : node + offset_;
            rev_comp_cache_.Put(node, rev_node);
            return rev_node;
        }
    }

    try {
        return is_palindrome_cache_.Get(node) ? node : node + offset_;

    } catch (...) {
        std::string seq = graph_.get_node_sequence(node);
        std::string rev_seq = seq;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        bool palindrome = (rev_seq == seq);

        assert(palindrome || graph_.kmer_to_node(rev_seq) == DeBruijnGraph::npos);

        is_palindrome_cache_.Put(node, palindrome);
        return palindrome ? node : node + offset_;
    }

    assert(false && "All cases should have been captured until now.");
    return DeBruijnGraph::npos;
}

void CanonicalDBG::reverse_complement(std::string &seq,
                                      std::vector<node_index> &path) const {
    ::reverse_complement(seq.begin(), seq.end());

    if (primary_) {
        std::vector<node_index> rev_path(path.size());
        std::transform(path.begin(), path.end(), rev_path.rbegin(), [&](node_index i) {
            return reverse_complement(i);
        });
        std::swap(path, rev_path);

    } else {
        size_t old_size = path.size();
        path = map_sequence_to_nodes(*this, seq);
        assert(path.size() == old_size || dynamic_cast<const DBGSuccinctRange*>(&graph_));
        path.resize(old_size);
    }
}

} // namespace graph
} // namespace mtg
