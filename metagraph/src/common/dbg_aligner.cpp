#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "bounded_priority_queue.hpp"

DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_threads, uint64_t search_space_size) :
                            AnnotatedDBG(dbg, annotation, num_threads),
                            search_space_size_(search_space_size) {
    // Substitution loss for each pair of nucleotides.
    // Transition and transversion mutations have different loss values.
    sub_loss_ = {
        {'a', {{'t', 2}, {'c', 2}, {'g', 1}}},
        {'g', {{'t', 2}, {'c', 2}, {'a', 1}}},
        {'c', {{'a', 2}, {'g', 2}, {'t', 1}}},
        {'t', {{'a', 2}, {'g', 2}, {'c', 1}}}};
}

DBGAligner::AlignedPath DBGAligner::align(const std::string& sequence) const {
    if (sequence.size() < graph_->get_k())
        return AlignedPath(sequence.end());

    BoundedPriorityQueue<AlignedPath> queue (search_space_size_);
    queue.push(AlignedPath(std::begin(sequence)));

    while (!queue.empty()) {
        AlignedPath path = queue.top();
        queue.pop();

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        node_index seed = path.size() > 0 ? path.back() : graph_->npos;
        graph_->extend_from_seed(path.get_sequence_it(), std::end(sequence),
                                       [&](node_index node) {
                                            path.push_back(node, 0, annotator_->get(node));},
                                       [&]() { return (path.size() > 0 &&
                                                graph_->outdegree(path.back()) != 1); },
                                       seed);

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        inexact_map(path, std::end(sequence),
                    [&](node_index node, std::string::const_iterator last_mapped_position) {
                        AlignedPath alternative_path(path);
                        alternative_path.set_sequence_it(last_mapped_position);
                        alternative_path.push_back(node, single_node_loss(node,
                            alternative_path.get_sequence_it()), annotator_->get(node));
                        queue.push(alternative_path);
                        } );
    }
    return AlignedPath(sequence.end());
}

float DBGAligner::single_node_loss(node_index node, std::string::const_iterator begin) const {
    // TODO: Compute indel loss as well.
    std::string node_sequence = graph_->get_node_sequence(node);
    auto node_sequence_it = std::begin(node_sequence);
    std::string::const_iterator end = begin + graph_->get_k();
    float loss = 0;
    while (begin != end) {
        if (std::tolower(*node_sequence_it) != std::tolower(*begin))
            loss += sub_loss_.at(std::tolower(*node_sequence_it)).at(std::tolower(*begin));
        ++ begin;
        ++ node_sequence_it;
    }
    return loss;
}

void DBGAligner::randomly_pick_strategy(std::vector<node_index> out_neighbors,
                                        const std::function<void(node_index)> &callback) const {
    std::random_device r;
    std::default_random_engine engine(r());
    std::uniform_int_distribution<int> uniform_dist(0, out_neighbors.size() - 1);
    callback(out_neighbors[uniform_dist(engine)]);
}

void DBGAligner::pick_all_strategy(std::vector<node_index> out_neighbors,
                                  const std::function<void(node_index)> &callback) const {
    for (auto &out_neighbor : out_neighbors) {
        callback(out_neighbor);
    }
}

void DBGAligner::inexact_map(const AlignedPath &path,
                             std::string::const_iterator end,
                             const std::function<void(node_index,
                                std::string::const_iterator)> &callback) const {
    if (path.get_sequence_it() + graph_->get_k() > end)
        return;
    // TODO: Do inexact_seeding in case last_mapped_node is npos.
    std::vector<node_index> out_neighbors;
    graph_->adjacent_outgoing_nodes(path.back(), &out_neighbors);

    pick_all_strategy(out_neighbors, [&](node_index node) {
                           callback(node, path.get_sequence_it()); });
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
