#include "masked_graph.hpp"

#include <stack>

#include "sequence_graph.hpp"


MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<DeBruijnGraph> graph, bit_vector *mask)
      : graph_(graph),
        is_target_mask_(mask) {
    assert(!is_target_mask_.get() || !(*mask)[DeBruijnGraph::npos]);
}

MaskedDeBruijnGraph
::MaskedDeBruijnGraph(std::shared_ptr<DeBruijnGraph> graph,
                      const std::function<bool(const node_index&)> &is_target_callback)
      : graph_(graph),
        is_target_callback_(new std::function<bool(const node_index&)>(is_target_callback)) {}

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
    adjacent_outgoing_nodes(node, &incoming);
    return std::count_if(incoming.begin(), incoming.end(),
                         [&](const auto &index) { return in_graph(index); });
}


// TODO: rewrite to be more like call_sequences
void MaskedDeBruijnGraph
::call_outgoing_paths(node_index node,
                      const std::function<void(const std::string&)> &callback,
                      size_t max_num_steps) const {
    if (!max_num_steps)
        return;

    max_num_steps += get_k() - 1;
    std::stack<std::pair<node_index, std::string>> paths;
    paths.emplace(node, get_node_sequence(node));
    while (paths.size()) {
        auto path = std::move(paths.top());
        paths.pop();
        while (path.second.size() < max_num_steps) {
            size_t cur_size = path.second.size();
            call_outgoing_kmers(
                path.first,
                [&](const auto &index, char out_char) {
                    if (path.second.size() == cur_size) {
                        path.second.push_back(out_char);
                        path.first = index;
                    } else {
                        paths.emplace(index,
                                      path.second.substr(0, path.second.length() - 1)
                                          + out_char);
                    }
                }
            );
            if (path.second.size() == cur_size)
                break;
        }
        if (path.second.size() == max_num_steps)
            callback(path.second);
    }
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

void MaskedDeBruijnGraph
::call_nodes(const std::function<void(const node_index&)> &callback,
             const std::function<bool()> &stop_early) const {
    if (is_target_callback_.get()) {
        for (size_t i = 1; i <= graph_->num_nodes(); ++i) {
            assert(i != npos);
            if ((*is_target_callback_)(i)) {
                callback(i);

                if (stop_early())
                    return;
            }
        }
    } else {
        if (is_target_mask_.get()) {
            uint64_t num_set_bits = is_target_mask_->num_set_bits();
            for (size_t i = 1; i <= num_set_bits; ++i) {
                assert(is_target_mask_->select1(i) != npos);
                callback(is_target_mask_->select1(i));

                if (stop_early())
                    return;
            }
        } else {
            for (DeBruijnGraph::node_index i = 1; i <= graph_->num_nodes(); ++i) {
                callback(i);

                if (stop_early())
                    return;
            }
        }
    }
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

bool MaskedDeBruijnGraph::load(const std::string &filename_base) {
    try {
        const_cast<DeBruijnGraph*>(graph_.get())->load(filename_base);

        std::ifstream in(filename_base + ".edgemask", std::ios::binary);
        is_target_mask_->load(in);
    } catch (...) {
        return false;
    }

    return true;
}

void MaskedDeBruijnGraph::serialize(const std::string &filename_base) const {
    if (is_target_callback_.get())
        throw std::runtime_error("Not implemented");

    graph_->serialize(filename_base);

    std::ofstream out(filename_base + ".edgemask", std::ios::binary);
    is_target_mask_->serialize(out);
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

    if (is_target_callback_.get() != other.is_target_callback_.get())
        return false;

    if (is_target_mask_.get() == other.is_target_mask_.get())
        return true;

    if (!is_target_mask_.get() || !other.is_target_mask_.get())
        return false;

    if (is_target_mask_->num_set_bits() != other.is_target_mask_->num_set_bits())
        return false;

    uint64_t cur_bit = is_target_mask_->num_set_bits();
    while (cur_bit) {
        if (is_target_mask_->select1(cur_bit) != other.is_target_mask_->select1(cur_bit))
            return false;

        cur_bit--;
    }

    return true;
}

bool MaskedDeBruijnGraph::operator==(const DeBruijnGraph &other) const {
    if (dynamic_cast<const MaskedDeBruijnGraph*>(&other))
        return operator==(*dynamic_cast<const MaskedDeBruijnGraph*>(&other));

    return *graph_ == other;
}
