#include "masked_graph.hpp"

#include <stack>

#include "sequence_graph.hpp"
#include "serialization.hpp"
#include "dbg_succinct.hpp"


MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                      std::unique_ptr<bitmap>&& kmers_in_graph)
      : graph_(graph),
        kmers_in_graph_(std::move(kmers_in_graph)) {
    assert(kmers_in_graph_.get());
    assert(kmers_in_graph_->size() == graph->num_nodes() + 1);
}

MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                      std::function<bool(const DeBruijnGraph::node_index&)>&& callback,
                      size_t num_set_bits)
      : MaskedDeBruijnGraph(graph,
                            std::make_unique<bitmap_lazy>(
                                std::move(callback),
                                graph->num_nodes() + 1,
                                num_set_bits)
                            ) {}

// Traverse the outgoing edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse(node_index node, char next_char) const {
    assert(in_graph(node));

    auto index = graph_->traverse(node, next_char);
    return index && in_graph(index) ? index : DeBruijnGraph::npos;
}
// Traverse the incoming edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse_back(node_index node, char prev_char) const {
    assert(in_graph(node));

    auto index = graph_->traverse_back(node, prev_char);
    return index && in_graph(index) ? index : DeBruijnGraph::npos;
}

size_t MaskedDeBruijnGraph::outdegree(node_index node) const {
    assert(in_graph(node));

    size_t outdegree = 0;
    graph_->adjacent_outgoing_nodes(node, [&](auto index) { outdegree += in_graph(index); });
    return outdegree;
}

size_t MaskedDeBruijnGraph::indegree(node_index node) const {
    assert(in_graph(node));

    size_t indegree = 0;
    graph_->adjacent_incoming_nodes(node, [&](auto index) { indegree += in_graph(index); });
    return indegree;
}

void MaskedDeBruijnGraph
::adjacent_outgoing_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    graph_->adjacent_outgoing_nodes(node, [&](auto node) {
        if (in_graph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::adjacent_incoming_nodes(node_index node, const std::function<void(node_index)> &callback) const {
    assert(in_graph(node));

    graph_->adjacent_incoming_nodes(node, [&](auto node) {
        if (in_graph(node))
            callback(node);
    });
}

void MaskedDeBruijnGraph
::call_outgoing_kmers(node_index kmer,
                      const OutgoingEdgeCallback &callback) const {
    assert(in_graph(kmer));

    graph_->call_outgoing_kmers(
        kmer,
        [&](const auto &index, auto c) {
            if (in_graph(index))
                callback(index, c);
        }
    );
}

void MaskedDeBruijnGraph
::call_incoming_kmers(node_index kmer,
                      const IncomingEdgeCallback &callback) const {
    assert(in_graph(kmer));

    graph_->call_incoming_kmers(
        kmer,
        [&](const auto &index, auto c) {
            if (in_graph(index))
                callback(index, c);
        }
    );
}

void MaskedDeBruijnGraph
::call_sequences(const std::function<void(const std::string&)> &callback) const {

    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {

        bitmap_lazy mask(
            [&](const auto &i) {
                auto node = dbg_succ->boss_to_kmer_index(i);
                return node && in_graph(node);
            },
            dbg_succ->get_boss().num_edges() + 1
        );

        // TODO: call_sequences will eventually call all indexes in bitmap_lazy
        // this is inefficient! pass sdsl::bit_vector instead
        dbg_succ->get_boss().call_sequences(callback, &mask);

    } else {
        DeBruijnGraph::call_sequences(callback);
    }
}

void MaskedDeBruijnGraph
::call_unitigs(const std::function<void(const std::string&)> &callback,
               size_t min_tip_size) const {

    if (auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get())) {
        bitmap_lazy mask(
            [&](const auto &i) {
                auto node = dbg_succ->boss_to_kmer_index(i);
                return node && in_graph(node);
            },
            dbg_succ->get_boss().num_edges() + 1
        );

        // TODO: call_sequences will eventually call all indexes in bitmap_lazy
        // this is inefficient! pass sdsl::bit_vector instead
        dbg_succ->get_boss().call_unitigs(callback, min_tip_size, &mask);

    } else {
        DeBruijnGraph::call_unitigs(callback, min_tip_size);
    }
}

void MaskedDeBruijnGraph
::call_nodes(const std::function<void(node_index)> &callback,
             const std::function<bool()> &stop_early) const {
    kmers_in_graph_->call_ones(
        [&](auto index) {
            if (!index)
                return;

            assert(in_graph(index));

            if (!stop_early())
                callback(index);
        }
    );
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied
void MaskedDeBruijnGraph
::map_to_nodes(const std::string &sequence,
               const std::function<void(node_index)> &callback,
               const std::function<bool()> &terminate) const {
    graph_->map_to_nodes(
        sequence,
        [&](const node_index &index) {
            callback(index && in_graph(index) ? index : npos);
        },
        terminate
    );
}

// Traverse graph mapping sequence to the graph nodes
// and run callback for each node until the termination condition is satisfied.
// Guarantees that nodes are called in the same order as the input sequence
void MaskedDeBruijnGraph
::map_to_nodes_sequentially(std::string::const_iterator begin,
                            std::string::const_iterator end,
                            const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate) const {
    graph_->map_to_nodes_sequentially(
        begin, end,
        [&](const node_index &index) {
            callback(index && in_graph(index) ? index : npos);
        },
        terminate
    );
}

// TODO: This should ideally return kmers_in_graph_->num_set_bits(), but the
//       API assumes that all node indices from 1 to num_nodes are valid. Until
//       we remove that restriction that, this will stay like this.
uint64_t MaskedDeBruijnGraph::num_nodes() const {
    return graph_->num_nodes();
}

// Get string corresponding to |node_index|.
// Note: Not efficient if sequences in nodes overlap. Use sparingly.
std::string MaskedDeBruijnGraph::get_node_sequence(node_index index) const {
    assert(in_graph(index));

    return graph_->get_node_sequence(index);
}

bool MaskedDeBruijnGraph::operator==(const MaskedDeBruijnGraph &other) const {
    return get_k() == other.get_k()
            && is_canonical_mode() == other.is_canonical_mode()
            && num_nodes() == other.num_nodes()
            && *kmers_in_graph_ == *other.kmers_in_graph_
            && *graph_ == *other.graph_;
}

bool MaskedDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (get_k() != other.get_k()
            || is_canonical_mode() != other.is_canonical_mode()
            || num_nodes() != other.num_nodes())
        return false;

    if (dynamic_cast<const MaskedDeBruijnGraph*>(&other))
        return operator==(dynamic_cast<const MaskedDeBruijnGraph&>(other));

    return DeBruijnGraph::operator==(other);
}
