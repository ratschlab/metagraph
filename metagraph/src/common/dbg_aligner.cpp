#include "dbg_aligner.hpp"

#include <string>
#include <random>

DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_threads) : AnnotatedDBG(dbg, annotation, num_threads) {
    sub_loss = {
        {'a', {{'t', 2}, {'c', 2}, {'g', 1}}},
        {'g', {{'t', 2}, {'c', 2}, {'a', 1}}},
        {'c', {{'a', 2}, {'g', 2}, {'t', 1}}},
        {'t', {{'a', 2}, {'g', 2}, {'c', 1}}}};
}

DBGAligner::AlignedPath DBGAligner::align(const std::string& sequence) const {
    // The current strategy is to make choice in point
    // rather than exploring the possible solutions.
    if (sequence.size() < graph_->get_k())
        return {};

    std::string mapped_sequence(sequence);

    // TODO: What is the best container for the queue?
    std::priority_queue<AlignedPath, std::vector<AlignedPath>, AlignedPathCompare> queue;
    queue.push(AlignedPath());

    while (!queue.empty()) {
        AlignedPath path = queue.top();
        queue.pop();
        std::string::iterator sequence_it = std::begin(mapped_sequence) + path.size();

        if (sequence_it + graph_->get_k() > std::end(mapped_sequence))
            return path;

        node_index seed = path.size() > 0 ? path.back() : graph_->npos;
        graph_->extend_from_seed(sequence_it, std::end(mapped_sequence),
                             [&](node_index node) {
                                path.push_back(node, 0, annotator_->get(node));
                                ++ sequence_it; },
                             [&]() {
                                return (path.size() > 0 && graph_->outdegree(path.back()) != 1); },
                                seed);
        if (sequence_it + graph_->get_k() > std::end(mapped_sequence))
            return path;

        inexact_map(sequence_it, std::end(mapped_sequence), path.back(), path.get_labels(),
                    [&](node_index node, std::string::const_iterator last_mapped_position) {
                        AlignedPath alternative_path(path);
                        alternative_path.push_back(node, single_node_loss(node, sequence_it,
                            last_mapped_position), annotator_->get(node));
                        queue.push(alternative_path);
                        } );
        ++ sequence_it;
    }
    return {};
}

float DBGAligner::single_node_loss(const node_index& node, std::string::const_iterator begin,
                                   std::string::const_iterator end) const {
    auto node_sequence_it = std::begin(graph_->get_node_sequence(node));
    float loss = 0;
    while (begin != end) {
        if (std::tolower(*node_sequence_it) != std::tolower(*begin))
            loss += sub_loss.at(std::tolower(*node_sequence_it)).at(std::tolower(*begin));
        // TODO: compute indel loss as well.
        ++ begin;
        ++ node_sequence_it;
    }
    return loss;
}

// Strategy: Randomly choose between possible outgoing neighbors.
void DBGAligner::randomly_pick_strategy(std::vector<node_index> out_neighbors,
                                        const std::function<void(node_index)> &callback) const {
    std::random_device r;
    std::default_random_engine engine(r());
    std::uniform_int_distribution<int> uniform_dist(0, out_neighbors.size() - 1);
    callback(out_neighbors[uniform_dist(engine)]);
}

// Strategy: Call callback for all edges. This is equivalent to exhaustive search.
void DBGAligner::pick_all_strategy(std::vector<node_index> out_neighbors,
                                  const std::function<void(node_index)> &callback) const {
    for (auto &out_neighbor : out_neighbors) {
        callback(out_neighbor);
    }
}

//TODO: Write a strategy to handle indels.
void DBGAligner::inexact_map(std::string::const_iterator begin,
                             std::string::const_iterator end,
                             const node_index &last_mapped_node,
                             const Annotator::VLabels &,
                             const std::function<void(node_index,
                                std::string::const_iterator)> &callback,
                             const std::function<bool()> &terminate) const {
    if (terminate())
        return;
    if (begin + graph_->get_k() > end)
        return;
    // TODO: Do inexact_seeding in case last_mapped_node is npos.
    std::vector<node_index> out_neighbors;
    graph_->adjacent_outgoing_nodes(last_mapped_node, &out_neighbors);

    pick_all_strategy(out_neighbors, [&](node_index node) {
                           callback(node, begin + graph_->get_k()); });
}

std::string DBGAligner::get_path_sequence(const std::vector<node_index>& path) const {
    if (path.size() < 1)
        return "";
    std::string sequence = graph_->get_node_sequence(path.front());
    for (auto path_it = std::next(path.begin()); path_it != path.end(); ++path_it) {
        sequence.push_back(graph_->get_node_sequence(*path_it).back());
    }
    return sequence;
}
