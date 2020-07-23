#include "canonical_dbg.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;

CanonicalDBG::CanonicalDBG(std::shared_ptr<const DeBruijnGraph> graph, size_t cache_size)
      : const_graph_ptr_(graph),
        offset_(graph_.max_index()),
        alph_map_({ graph_.alphabet().size() }),
        child_node_cache_(cache_size),
        parent_node_cache_(cache_size) {
    if (graph_.is_canonical_mode())
        throw std::runtime_error("CanonicalDBG should not be used as a wrapper for a canonical graph");

    for (size_t i = 0; i < graph_.alphabet().size(); ++i) {
        alph_map_[graph_.alphabet()[i]] = i;
    }
}

CanonicalDBG::CanonicalDBG(std::shared_ptr<DeBruijnGraph> graph, size_t cache_size)
      : CanonicalDBG(std::dynamic_pointer_cast<const DeBruijnGraph>(graph),
                             cache_size) { graph_ptr_ = graph; }

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
    std::string rev_seq(sequence);
    ::reverse_complement(rev_seq.begin(), rev_seq.end());

    std::vector<node_index> rev_path = map_sequence_to_nodes(graph_, rev_seq);

    auto it = rev_path.rbegin();
    graph_.map_to_nodes_sequentially(sequence, [&](node_index i) {
        assert(it != rev_path.rend());

        // if k is even, then some k-mers may be palindromes
        assert(!(get_k() & 1) || i == DeBruijnGraph::npos || *it == DeBruijnGraph::npos);

        if (i != DeBruijnGraph::npos) {
            assert(i <= offset_ * 2);
            callback(i);
        } else if (*it != DeBruijnGraph::npos) {
            assert(set_offset(*it) <= offset_ * 2);
            callback(set_offset(*it));
        } else {
            callback(DeBruijnGraph::npos);
        }
        ++it;
    }, terminate);
}

void CanonicalDBG::map_to_nodes(std::string_view sequence,
                                const std::function<void(node_index)> &callback,
                                const std::function<bool()> &terminate) const {
    map_to_nodes_sequentially(sequence, [&](node_index i) {
        callback(i != DeBruijnGraph::npos ? get_base_node(i) : i);
    }, terminate);
}

void CanonicalDBG
::call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (is_rev_comp(node)) {
        call_incoming_kmers(unset_offset(node), [&](node_index next, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(next), c);
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
        size_t alph_size = alphabet.size();
        std::vector<node_index> children(alph_size);

        graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
            children[alph_map_[c]] = next;
            --alph_size;
        });

        if (alph_size) {
            node_index next;
            std::string rev_seq = get_node_sequence(node).substr(1) + std::string(1, '\0');
            ::reverse_complement(rev_seq.begin(), rev_seq.end());
            assert(rev_seq[0] == '\0');
            for (char c : graph_.alphabet()) {
                char d = c;
                ::reverse_complement(&d, &d + 1);
                if (children[alph_map_[d]] != DeBruijnGraph::npos)
                    continue;

                rev_seq[0] = c;
                next = map_sequence_to_nodes(graph_, rev_seq)[0];
                if (next != DeBruijnGraph::npos)
                    children[alph_map_[d]] = set_offset(next);
            }
        }

        child_node_cache_.Put(node, children);
        for (size_t i = 0; i < children.size(); ++i) {
            if (children[i] != DeBruijnGraph::npos)
                callback(children[i], alphabet[i]);
        }
    }
}

void CanonicalDBG
::call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (is_rev_comp(node)) {
        call_outgoing_kmers(unset_offset(node), [&](node_index prev, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(prev), c);
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
        size_t alph_size = alphabet.size();
        std::vector<node_index> parents(alph_size);

        graph_.call_incoming_kmers(node, [&](node_index prev, char c) {
            parents[alph_map_[c]] = prev;
            --alph_size;
        });

        if (alph_size) {
            node_index prev;
            std::string rev_seq = std::string(1, '\0') + get_node_sequence(node).substr(0, get_k() - 1);
            ::reverse_complement(rev_seq.begin(), rev_seq.end());
            assert(rev_seq.back() == '\0');
            for (char c : graph_.alphabet()) {
                char d = c;
                ::reverse_complement(&d, &d + 1);
                if (parents[alph_map_[d]] != DeBruijnGraph::npos)
                    continue;

                rev_seq.back() = c;
                prev = map_sequence_to_nodes(graph_, rev_seq)[0];
                if (prev != DeBruijnGraph::npos)
                    parents[alph_map_[d]] = set_offset(prev);
            }
        }

        parent_node_cache_.Put(node, parents);
        for (size_t i = 0; i < parents.size(); ++i) {
            if (parents[i] != DeBruijnGraph::npos)
                callback(parents[i], alphabet[i]);
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
    if (is_rev_comp(node)) {
        ::reverse_complement(&next_char, &next_char + 1);
        node = traverse_back(unset_offset(node), next_char);
        return node != DeBruijnGraph::npos ? reverse_complement(node) : DeBruijnGraph::npos;
    } else {
        node_index next = graph_.traverse(node, next_char);
        if (next != DeBruijnGraph::npos)
            return next;

        std::string rev_seq = get_node_sequence(node).substr(1) + next_char;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        next = map_sequence_to_nodes(graph_, rev_seq)[0];
        return next != DeBruijnGraph::npos ? reverse_complement(next) : next;
    }
}

DeBruijnGraph::node_index CanonicalDBG::traverse_back(node_index node,
                                                      char prev_char) const {
    assert(node <= offset_ * 2);
    if (is_rev_comp(node)) {
        ::reverse_complement(&prev_char, &prev_char + 1);
        node = traverse(unset_offset(node), prev_char);
        return node != DeBruijnGraph::npos ? reverse_complement(node) : DeBruijnGraph::npos;
    } else {
        node_index prev = graph_.traverse_back(node, prev_char);
        if (prev != DeBruijnGraph::npos)
            return prev;

        std::string rev_seq = std::string(1, prev_char)
            + get_node_sequence(node).substr(0, get_k() - 1);
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        prev = map_sequence_to_nodes(graph_, rev_seq)[0];
        return prev != DeBruijnGraph::npos ? reverse_complement(prev) : prev;
    }
}

void CanonicalDBG::call_nodes(const std::function<void(node_index)> &callback,
                              const std::function<bool()> &stop_early) const {
    graph_.call_nodes([&](node_index i) {
                          callback(i);
                          if (!stop_early())
                              callback(set_offset(i));
                      }, stop_early);
}

bool CanonicalDBG::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const CanonicalDBG*>(&other)) {
        return *this == dynamic_cast<const CanonicalDBG&>(other);
    }

    return DeBruijnGraph::operator==(other);
}

} // namespace graph
} // namespace mtg
