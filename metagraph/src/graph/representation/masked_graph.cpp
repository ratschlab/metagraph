#include "masked_graph.hpp"

#include "common/serialization.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace graph {

// Traverse the outgoing edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse(node_index node, char next_char) const {
    assert(in_subgraph(node));

    auto index = graph_->traverse(node, next_char);
    return index && in_subgraph(index) ? index : npos;
}
// Traverse the incoming edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse_back(node_index node, char prev_char) const {
    assert(in_subgraph(node));

    auto index = graph_->traverse_back(node, prev_char);
    return index && in_subgraph(index) ? index : npos;
}

size_t MaskedDeBruijnGraph::outdegree(node_index node) const {
    assert(in_subgraph(node));

    size_t outdegree = 0;
    graph_->adjacent_outgoing_nodes(node, [&](auto index) {
        outdegree += in_subgraph(index);
    });
    return outdegree;
}

size_t MaskedDeBruijnGraph::indegree(node_index node) const {
    assert(in_subgraph(node));

    size_t indegree = 0;
    graph_->adjacent_incoming_nodes(node, [&](auto index) {
        indegree += in_subgraph(index);
    });
    return indegree;
}

void MaskedDeBruijnGraph
::adjacent_outgoing_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_subgraph(node));

    graph_->adjacent_outgoing_nodes(node, [&](auto node) {
        if (in_subgraph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::adjacent_incoming_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_subgraph(node));

    graph_->adjacent_incoming_nodes(node, [&](auto node) {
        if (in_subgraph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::call_outgoing_kmers(node_index kmer, const OutgoingEdgeCallback &callback) const {
    assert(in_subgraph(kmer));

    graph_->call_outgoing_kmers(kmer, [&](const auto &index, auto c) {
        if (in_subgraph(index))
            callback(index, c);
    });
}

void MaskedDeBruijnGraph
::call_incoming_kmers(node_index kmer, const IncomingEdgeCallback &callback) const {
    assert(in_subgraph(kmer));

    graph_->call_incoming_kmers(kmer, [&](const auto &index, auto c) {
        if (in_subgraph(index))
            callback(index, c);
    });
}

bit_vector_stat get_boss_mask(const DBGSuccinct &dbg_succ,
                              const bitmap &kmers_in_graph,
                              bool only_valid_nodes_in_mask) {
    sdsl::bit_vector mask_bv(dbg_succ.get_boss().num_edges() + 1, false);
    if (only_valid_nodes_in_mask) {
        kmers_in_graph.call_ones([&](auto i) {
            assert(dbg_succ.kmer_to_boss_index(i));
            mask_bv[dbg_succ.kmer_to_boss_index(i)] = true;
        });
    } else {
        dbg_succ.call_nodes([&](auto i) {
            assert(dbg_succ.kmer_to_boss_index(i));
            if (kmers_in_graph[i])
                mask_bv[dbg_succ.kmer_to_boss_index(i)] = true;
        });
    }
    return bit_vector_stat(std::move(mask_bv));
}

void MaskedDeBruijnGraph::call_sequences(const CallPath &callback,
                                         size_t num_threads,
                                         bool kmers_in_single_form) const {
    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {
        bit_vector_stat mask = get_boss_mask(*dbg_succ, *kmers_in_graph_,
                                             only_valid_nodes_in_mask_);

        dbg_succ->get_boss().call_sequences([&](std::string&& sequence, auto&& path) {
            for (auto &node : path) {
                node = dbg_succ->boss_to_kmer_index(node);
            }
            callback(sequence, path);

        }, num_threads, kmers_in_single_form, &mask);

    } else {
        DeBruijnGraph::call_sequences(callback, num_threads, kmers_in_single_form);
    }
}

void MaskedDeBruijnGraph::call_unitigs(const CallPath &callback,
                                       size_t num_threads,
                                       size_t min_tip_size,
                                       bool kmers_in_single_form) const {
    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {
        bit_vector_stat mask = get_boss_mask(*dbg_succ, *kmers_in_graph_,
                                             only_valid_nodes_in_mask_);

        dbg_succ->get_boss().call_unitigs([&](std::string&& sequence, auto&& path) {
            for (auto &node : path) {
                node = dbg_succ->boss_to_kmer_index(node);
            }
            callback(sequence, path);

        }, num_threads, min_tip_size, kmers_in_single_form, &mask);

    } else {
        DeBruijnGraph::call_unitigs(callback,
                                    num_threads,
                                    min_tip_size,
                                    kmers_in_single_form);
    }
}

void MaskedDeBruijnGraph
::call_nodes(const std::function<void(node_index)> &callback,
             const std::function<bool()> &stop_early) const {
    assert(max_index() + 1 == kmers_in_graph_->size());

    bool stop = false;

    if (only_valid_nodes_in_mask_) {
        // iterate only through the nodes marked in the mask
        // TODO: add terminate<bool(void)> to call_ones
        kmers_in_graph_->call_ones([&](auto index) {
            if (stop || !index)
                return;

            assert(in_subgraph(index));

            if (stop_early()) {
                stop = true;
            } else {
                callback(index);
            }
        });
    } else {
        // call all nodes in the base graph and check the mask
        graph_->call_nodes(
            [&](auto index) {
                if (in_subgraph(index))
                    callback(index);
            },
            stop_early
        );
    }
}

void MaskedDeBruijnGraph
::call_kmers(const std::function<void(node_index, const std::string&)> &callback) const {
    assert(max_index() + 1 == kmers_in_graph_->size());

    if (only_valid_nodes_in_mask_) {
        // iterate only through the nodes marked in the mask
        kmers_in_graph_->call_ones([&](node_index index) {
            if (index) {
                assert(in_subgraph(index));
                // TODO: make this more efficient
                callback(index, get_node_sequence(index));
            }
        });
    } else {
        // call all nodes in the base graph and check the mask
        graph_->call_kmers([&](node_index index, const std::string &seq) {
            if (in_subgraph(index))
                callback(index, seq);
        });
    }
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
void MaskedDeBruijnGraph::map_to_nodes(std::string_view sequence,
                                       const std::function<void(node_index)> &callback,
                                       const std::function<bool()> &terminate) const {
    graph_->map_to_nodes(
        sequence,
        [&](const node_index &index) {
            callback(index && in_subgraph(index) ? index : npos);
        },
        terminate
    );
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence
void MaskedDeBruijnGraph
::map_to_nodes_sequentially(std::string_view sequence,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate,
                            const std::function<bool()> &skip) const {
    graph_->map_to_nodes_sequentially(
        sequence,
        [&](const node_index &index) {
            callback(index && in_subgraph(index) ? index : npos);
        },
        terminate,
        skip
    );
}

// Get string corresponding to |node_index|.
// Note: Not efficient if sequences in nodes overlap. Use sparingly.
std::string MaskedDeBruijnGraph::get_node_sequence(node_index index) const {
    assert(in_subgraph(index));

    return graph_->get_node_sequence(index);
}

bool MaskedDeBruijnGraph::operator==(const MaskedDeBruijnGraph &other) const {
    return get_k() == other.get_k()
            && get_mode() == other.get_mode()
            && num_nodes() == other.num_nodes()
            && *kmers_in_graph_ == *other.kmers_in_graph_
            && *graph_ == *other.graph_;
}

bool MaskedDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || get_mode() != other.get_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (dynamic_cast<const MaskedDeBruijnGraph*>(&other))
        return operator==(dynamic_cast<const MaskedDeBruijnGraph&>(other));

    return DeBruijnGraph::operator==(other);
}

} // namespace graph
} // namespace mtg
