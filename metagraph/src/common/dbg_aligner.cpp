#include "dbg_aligner.hpp"

#include <string>
#include <random>


DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_threads) : AnnotatedDBG(dbg, annotation, num_threads) {
}

std::vector<DBGAligner::node_index> DBGAligner::align(const std::string& sequence) const {
    // The current strategy is to make choice in point
    // rather than exploring the possible solutions.
    if (sequence.size() < graph_->get_k())
        return {};

    std::string mapped_sequence(sequence);
    std::string::iterator sequence_it = mapped_sequence.begin();

    std::vector<node_index> path;
    Annotator::VLabels label_set;
    while (sequence_it + graph_->get_k() <= std::end(mapped_sequence)) {
        graph_->map_to_nodes_sequentially(sequence_it, std::end(mapped_sequence),
                             [&](node_index index) {
                                path.push_back(index);
                                auto new_labels = annotator_->get(index);
                                label_set.insert(std::end(label_set),
                                                 std::begin(new_labels),
                                                 std::end(new_labels));
                                ++ sequence_it; });

        if (sequence_it + graph_->get_k() > std::end(mapped_sequence))
            break;
        inexact_map(sequence_it, std::end(mapped_sequence), path.back(), label_set,
                    [&](node_index index) {
                        path.push_back(index);
                        auto new_labels = annotator_->get(index);
                        label_set.insert(std::end(label_set),
                                         std::begin(new_labels),
                                         std::end(new_labels));
                        // Update alterations to the original sequence
                        // to apply for later alignments.
                        *(sequence_it + graph_->get_k() - 1) = graph_->get_node_sequence(index).back();
                        ++ sequence_it; } );
    }
    return path;
}

void DBGAligner::inexact_map(std::string::const_iterator begin,
                             std::string::const_iterator end,
                             const node_index &last_mapped_node,
                             const Annotator::VLabels &,
                             const std::function<void(node_index)> &callback,
                             const std::function<bool()> &terminate) const {
    if (terminate())
        return;

    if (begin + graph_->get_k() > end)
        return;

    std::vector<node_index> out_neighbors;

    // TODO: do inexact_seeding in case last_mapped_node is npos.

    graph_->adjacent_outgoing_nodes(last_mapped_node, &out_neighbors);

    // Strategy: Randomly choose between possible outgoing neighbors.
    std::random_device r;
    std::default_random_engine engine(r());
    std::uniform_int_distribution<int> uniform_dist(0, out_neighbors.size() - 1);
    node_index next_node = out_neighbors[uniform_dist(engine)];

    //callback(graph_->traverse(last_mapped_node, graph_->get_kmer_last_char(next_node)));
    callback(next_node);
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
