#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "bounded_priority_queue.hpp"

DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_threads, size_t search_space_size,
                       float path_loss_weak_threshold,
                       float insertion_penalty, float deletion_penalty, size_t max_sw_table_size) :
                            AnnotatedDBG(dbg, annotation, num_threads),
                            search_space_size_(search_space_size),
                            path_loss_weak_threshold_(path_loss_weak_threshold),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty),
                            max_sw_table_size_(max_sw_table_size) {
    // Substitution loss for each pair of nucleotides.
    // Transition and transversion mutations have different loss values.
    sub_loss_ = {
        {'a', {{'a', 0}, {'t', 2}, {'c', 2}, {'g', 1}}},
        {'g', {{'g', 0}, {'t', 2}, {'c', 2}, {'a', 1}}},
        {'c', {{'c', 0}, {'a', 2}, {'g', 2}, {'t', 1}}},
        {'t', {{'t', 0}, {'a', 2}, {'g', 2}, {'c', 1}}},
        {'$', {{'$', 2}}}};
}

DBGAligner::AlignedPath DBGAligner::align(const std::string& sequence) const {
    if (sequence.size() < graph_->get_k())
        return AlignedPath(sequence.end());

    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath> queue(search_space_size_);
    queue.push(std::move(AlignedPath(std::begin(sequence))));
    while (!queue.empty()) {
        auto path = std::move(queue.top());
        queue.pop();
        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        node_index seed = path.size() > 0 ? path.back() : graph_->npos;
        graph_->extend_from_seed(path.get_sequence_it(), std::end(sequence),
                                       [&](node_index node) {
                                            path.push_back(node, annotator_->get(node));},
                                       [&]() { return (path.size() > 0 &&
                                                (graph_->outdegree(path.back()) > 1)); },
                                       seed);

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        inexact_map(path, std::end(sequence),
                    [&](node_index node, std::string::const_iterator last_mapped_position) {
                        AlignedPath alternative_path(path);
                        alternative_path.set_sequence_it(last_mapped_position);
                        size_t sw_table_size = (alternative_path.size() + graph_->get_k() - 1)
                                               * (alternative_path.get_sequence_it() - std::begin(sequence));
                        // If the loss is less than the threshold or the Smith-Waterman table size is too large,
                        // continue with single char loss computation.
                        if (alternative_path.get_total_loss() < path_loss_weak_threshold_
                            || sw_table_size > max_sw_table_size_) {
                            // TODO: Construct sequence last char more efficiently.
                            alternative_path.push_back(node, annotator_->get(node),
                                single_char_loss(*(alternative_path.get_sequence_it() +
                                    graph_->get_k() - 1), graph_->get_node_sequence(node).back()));
                        } else {
                            alternative_path.push_back(node, annotator_->get(node));
                            alternative_path.update_total_loss(
                                whole_path_loss(alternative_path, std::begin(sequence)));
                        }

                        DPAlignmentKey alternative_key{.node = alternative_path.back(),
                                           .query_it = alternative_path.get_sequence_it()};
                        DPAlignmentValue alternative_value = {.parent = alternative_path.last_parent(),
                                                              .loss = alternative_path.get_total_loss()};
                        auto dp_alignment_it = dp_alignment.find(alternative_key);
                        if (dp_alignment_it == dp_alignment.end()) {
                            dp_alignment[alternative_key] = alternative_value;
                            queue.push(std::move(alternative_path));
                        }
                        else if (dp_alignment_it->second.loss > alternative_path.get_total_loss()) {
                            dp_alignment_it->second = alternative_value;
                            queue.push(std::move(alternative_path));
                        }
                        });
    }
    return AlignedPath(sequence.end());
}

float DBGAligner::single_char_loss(char char_in_query, char char_in_graph) const {
    try {
        return sub_loss_.at(std::tolower(char_in_query)).at(std::tolower(char_in_graph));
    }
    catch (const std::out_of_range&) {
        return sub_loss_.at('$').at('$');
    }
}

float DBGAligner::whole_path_loss(const AlignedPath& path, std::string::const_iterator begin) const {
    auto path_sequence = get_path_sequence(path.get_nodes());
    std::string query_sequence(begin, path.get_sequence_it() + graph_->get_k() - 1);
    float alignment_loss[query_sequence.size() + 1][path_sequence.size() + 1] = {};
    for (size_t i = 1; i <= query_sequence.size(); ++i) {
        for (size_t j = 1; j <= path_sequence.size(); ++j) {
            alignment_loss[i][j] = std::min({alignment_loss[i-1][j-1] +
                    single_char_loss(query_sequence[i-1], path_sequence[j-1]),
                alignment_loss[i][j-1] + insertion_penalty_, alignment_loss[i-1][j] + deletion_penalty_});
        }
    }
    return alignment_loss[query_sequence.size()][path_sequence.size()];
}

void DBGAligner::randomly_pick_strategy(std::vector<node_index> out_neighbors,
                                        const std::function<void(node_index)> &callback) const {
    // random_device creates a different seed everytime this function is called.
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
    // TODO: Do inexact_seeding in case last mapped node is npos.
    if (path.size() == 0)
        return;
    std::vector<node_index> out_neighbors;
   // TODO: char is not used currently.
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char)
                                             { out_neighbors.push_back(node); });
    pick_all_strategy(out_neighbors, [&](node_index node) {
                           callback(node, path.get_sequence_it()); });
}

std::string DBGAligner::get_path_sequence(const std::vector<node_index>& path) const {
    if (path.size() < 1)
        return "";
    // TODO: Construct sequence more efficiently.
    std::string sequence = graph_->get_node_sequence(path.front());
    for (auto path_it = std::next(path.begin()); path_it != path.end(); ++path_it) {
        sequence.push_back(graph_->get_node_sequence(*path_it).back());
    }
    return sequence;
}
