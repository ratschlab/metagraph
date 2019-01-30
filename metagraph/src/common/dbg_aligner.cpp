#include "dbg_aligner.hpp"

#include <string>


DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_threads) : AnnotatedDBG(dbg, annotation, num_threads) {
}

std::vector<DBGAligner::node_index> DBGAligner::align(const std::string& sequence) const {
    if (sequence.size() < graph_->get_k())
        return {};

    std::vector<node_index> path;
    Annotator::VLabels label_set;
    auto sequence_it = std::begin(sequence);
    while (sequence_it + graph_->get_k() <= std::end(sequence)) {
        graph_->map_to_nodes_sequentially(sequence_it, std::end(sequence),
                             [&](node_index index) {
                                path.push_back(index);
                                auto new_labels = annotator_->get(index);
                                label_set.insert(std::end(label_set),
                                                 std::begin(new_labels),
                                                 std::end(new_labels));
                                ++ sequence_it; });

        if (sequence_it + graph_->get_k() > sequence.end())
            break;
        // TODO: Perform inexact mapping.
    }
    return path;
}

std::string DBGAligner::get_path_sequence(const std::vector<node_index>& path) const {
    if(path.size() < 1)
        return "";
    std::string sequence = graph_->get_node_sequence(path.front());
    for (auto path_it = std::next(path.begin()); path_it != path.end(); ++path_it) {
        sequence.push_back(graph_->get_node_sequence(*path_it).back());
    }
    return sequence;
}
