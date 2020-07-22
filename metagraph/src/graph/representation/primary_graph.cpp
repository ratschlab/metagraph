#include "primary_graph.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/seq_tools/reverse_complement.hpp"


PrimaryDeBruijnGraph::PrimaryDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                                           size_t num_seqs_cached)
      : graph_ptr_(graph), offset_(graph_.max_index()), seq_cache_(num_seqs_cached) {
    if (!(get_k() & 1))
        throw std::runtime_error("PrimaryDeBruijnGraph not supported for even k.");
}


void PrimaryDeBruijnGraph::map_to_nodes(std::string_view sequence,
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

void PrimaryDeBruijnGraph
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    map_to_nodes(sequence, callback, terminate);
}

void PrimaryDeBruijnGraph
::call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (is_rev_comp(node)) {
        call_incoming_kmers(unset_offset(node), [&](node_index next, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(next), c);
        });
        return;
    } else {
        size_t alph_size = graph_.alphabet().size();
        graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
            callback(next, c);
            --alph_size;
        });

        if (!alph_size)
            return;

        node_index next;
        std::string rev_seq = get_node_sequence(node).substr(1) + std::string(1, '\0');
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        assert(rev_seq[0] == '\0');
        for (char c : graph_.alphabet()) {
            rev_seq[0] = c;
            next = map_sequence_to_nodes(graph_, rev_seq)[0];
            if (next != DeBruijnGraph::npos) {
                ::reverse_complement(&c, &c + 1);
                callback(set_offset(next), c);
            }
        }
    }
}

void PrimaryDeBruijnGraph
::call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const {
    assert(node);
    assert(node <= offset_ * 2);
    if (is_rev_comp(node)) {
        call_outgoing_kmers(unset_offset(node), [&](node_index prev, char c) {
            ::reverse_complement(&c, &c + 1);
            callback(reverse_complement(prev), c);
        });
        return;
    } else {
        size_t alph_size = graph_.alphabet().size();
        graph_.call_incoming_kmers(node, [&](node_index next, char c) {
            callback(next, c);
            --alph_size;
        });

        if (!alph_size)
            return;

        std::string rev_seq = std::string(1, '\0') + get_node_sequence(node).substr(0, get_k() - 1);
        node_index prev;
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        assert(rev_seq.back() == '\0');
        for (char c : graph_.alphabet()) {
            rev_seq.back() = c;
            prev = map_sequence_to_nodes(graph_, rev_seq)[0];
            if (prev != DeBruijnGraph::npos) {
                ::reverse_complement(&c, &c + 1);
                callback(set_offset(prev), c);
            }
        }
    }
}

void PrimaryDeBruijnGraph
::adjacent_outgoing_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    // TODO: better implementation
    call_outgoing_kmers(node, [&](node_index i, char) { callback(i); });
}

void PrimaryDeBruijnGraph
::adjacent_incoming_nodes(node_index node,
                          const std::function<void(node_index)> &callback) const {
    // TODO: better implementation
    call_incoming_kmers(node, [&](node_index i, char) { callback(i); });
}

size_t PrimaryDeBruijnGraph::outdegree(node_index node) const {
    // TODO: better implementation
    size_t outdegree = 0;
    adjacent_outgoing_nodes(node, [&](node_index) { ++outdegree; });
    return outdegree;
}

size_t PrimaryDeBruijnGraph::indegree(node_index node) const {
    // TODO: better implementation
    size_t indegree = 0;
    adjacent_incoming_nodes(node, [&](node_index) { ++indegree; });
    return indegree;
}

void PrimaryDeBruijnGraph::call_sequences(const CallPath &callback,
                                          size_t num_threads,
                                          bool kmers_in_single_form) const {
    if (kmers_in_single_form) {
        graph_.call_sequences(callback, num_threads, false);
    } else {
        // TODO: port over implementation from DBGSuccinct to DeBruijnGraph
        DeBruijnGraph::call_sequences(callback, num_threads, false);
    }
}

void PrimaryDeBruijnGraph::call_unitigs(const CallPath &callback,
                                        size_t num_threads,
                                        size_t min_tip_size,
                                        bool kmers_in_single_form) const {
    // TODO: port over implementation from DBGSuccinct to DeBruijnGraph
    DeBruijnGraph::call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
}

std::string PrimaryDeBruijnGraph::get_node_sequence(node_index index) const {
    assert(index <= offset_ * 2);
    std::string seq;
    node_index node = get_base_node(index);
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        try {
            seq = seq_cache_.Get(node);
        } catch (...) {
            seq = graph_.get_node_sequence(node);
            seq_cache_.Put(node, seq);
        }
    }

    if (node != index)
        ::reverse_complement(seq.begin(), seq.end());

    return seq;
}

DeBruijnGraph::node_index PrimaryDeBruijnGraph
::traverse(node_index node, char next_char) const {
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

DeBruijnGraph::node_index PrimaryDeBruijnGraph
::traverse_back(node_index node, char prev_char) const {
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

void PrimaryDeBruijnGraph::call_nodes(const std::function<void(node_index)> &callback,
                                      const std::function<bool()> &stop_early) const {
    graph_.call_nodes([&](node_index i) {
                          callback(i);
                          if (!stop_early())
                              callback(set_offset(i));
                      }, stop_early);
}

bool PrimaryDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const PrimaryDeBruijnGraph*>(&other)) {
        return *this == dynamic_cast<const PrimaryDeBruijnGraph&>(other);
    }

    return DeBruijnGraph::operator==(other);
}
