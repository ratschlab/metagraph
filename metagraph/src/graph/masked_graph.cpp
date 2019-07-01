#include "masked_graph.hpp"

#include <stack>

#include "sequence_graph.hpp"
#include "serialization.hpp"


MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph, bitmap *mask)
      : graph_(graph),
        is_target_mask_(mask) {
    if (!is_target_mask_.get())
        is_target_mask_.reset(new bitmap_lazy(
            [&](const auto &i) { return i != DeBruijnGraph::npos; },
            graph->num_nodes() + 1,
            graph->num_nodes()
        ));

    assert(is_target_mask_->size() == graph->num_nodes() + 1);
    assert(!(*is_target_mask_)[DeBruijnGraph::npos]);
}

MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                      std::function<bool(const DeBruijnGraph::node_index&)>&& callback,
                      size_t num_set_bits)
      : MaskedDeBruijnGraph(graph,
                            new bitmap_lazy(std::move(callback),
                                            graph->num_nodes() + 1,
                                            num_set_bits)) {}

// Traverse the outgoing edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse(node_index node, char next_char) const {
    auto index = graph_->traverse(node, next_char);
    return in_graph(index) ? index : DeBruijnGraph::npos;
}
// Traverse the incoming edge
MaskedDeBruijnGraph::node_index MaskedDeBruijnGraph
::traverse_back(node_index node, char prev_char) const {
    auto index = graph_->traverse_back(node, prev_char);
    return in_graph(index) ? index : DeBruijnGraph::npos;
}

size_t MaskedDeBruijnGraph::outdegree(node_index node) const {
    std::vector<node_index> outgoing;
    adjacent_outgoing_nodes(node, &outgoing);
    return std::count_if(outgoing.begin(), outgoing.end(),
                         [&](const auto &index) { return in_graph(index); });
}

size_t MaskedDeBruijnGraph::indegree(node_index node) const {
    std::vector<node_index> incoming;
    adjacent_incoming_nodes(node, &incoming);
    return std::count_if(incoming.begin(), incoming.end(),
                         [&](const auto &index) { return in_graph(index); });
}

// Given a node index and a pointer to a vector of node indices, iterates
// over all the outgoing edges and pushes back indices of their target nodes.
void MaskedDeBruijnGraph
::adjacent_outgoing_nodes(node_index node, std::vector<node_index>* target_nodes) const {
    assert(target_nodes);

    graph_->adjacent_outgoing_nodes(node, target_nodes);
    target_nodes->erase(
        std::remove_if(target_nodes->begin(), target_nodes->end(),
                       [&](const auto &index) { return !in_graph(index); }),
        target_nodes->end()
    );
}

// Given a node index and a pointer to a vector of node indices, iterates
// over all the incoming edges and pushes back indices of their source nodes.
void MaskedDeBruijnGraph
::adjacent_incoming_nodes(node_index node,
                          std::vector<node_index> *source_nodes) const {
    assert(source_nodes);

    graph_->adjacent_incoming_nodes(node, source_nodes);
    source_nodes->erase(
        std::remove_if(source_nodes->begin(), source_nodes->end(),
                       [&](const auto &index) { return !in_graph(index); }),
        source_nodes->end()
    );
}

void MaskedDeBruijnGraph
::call_outgoing_kmers(node_index kmer,
                      const OutgoingEdgeCallback &callback) const {
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
    graph_->call_incoming_kmers(
        kmer,
        [&](const auto &index, auto c) {
            if (in_graph(index))
                callback(index, c);
        }
    );
}

uint64_t MaskedDeBruijnGraph
::unmasked_outdegree(node_index node) const {
    return graph_->outdegree(node);
}

uint64_t MaskedDeBruijnGraph
::unmasked_indegree(node_index node) const {
    return graph_->indegree(node);
}

void MaskedDeBruijnGraph
::call_nodes(const std::function<void(node_index)> &callback,
             const std::function<bool()> &stop_early) const {
    assert(is_target_mask_.get());

    is_target_mask_->call_ones(
        [&](auto index) {
            if (!stop_early())
                callback(index);
        }
    );
}


void MaskedDeBruijnGraph
::add_sequence(const std::string &,
               bit_vector_dyn *) {
    throw std::runtime_error("Not implemented");
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
            if (in_graph(index))
                callback(index);
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
            if (in_graph(index))
                callback(index);
        },
        terminate
    );
}

uint64_t MaskedDeBruijnGraph::num_nodes() const {
    return graph_->num_nodes();
}

bool MaskedDeBruijnGraph::load(const std::string &) {
    throw std::runtime_error("Not implemented");
}

void MaskedDeBruijnGraph::serialize(const std::string &) const {
    throw std::runtime_error("Not implemented");
}

// Get string corresponding to |node_index|.
// Note: Not efficient if sequences in nodes overlap. Use sparingly.
std::string MaskedDeBruijnGraph::get_node_sequence(node_index index) const {
    return graph_->get_node_sequence(index);
}

bool MaskedDeBruijnGraph::operator==(const MaskedDeBruijnGraph &other) const {
    if (graph_ != other.graph_) {
        if (!graph_ || !other.graph_)
            return false;

        if (*graph_ != *other.graph_)
            return false;
    }

    if (is_target_mask_.get() == other.is_target_mask_.get())
        return true;

    if (!is_target_mask_.get() || !other.is_target_mask_.get())
        return false;

    if (is_target_mask_->num_set_bits() != other.is_target_mask_->num_set_bits())
        return false;

    uint64_t count = 0;
    uint64_t i = 0;
    is_target_mask_->call_ones([&](const auto &index) {
        if (i == count && (*other.is_target_mask_)[index])
            ++count;

        ++i;
    });

    return i == count;
}

bool MaskedDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const MaskedDeBruijnGraph*>(&other))
        return operator==(*dynamic_cast<const MaskedDeBruijnGraph*>(&other));

    return *graph_ == other;
}
