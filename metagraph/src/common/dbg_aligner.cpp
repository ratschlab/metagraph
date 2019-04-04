#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "bounded_priority_queue.hpp"
#include "ssw_cpp.h"

DBGAligner::DBGAligner(DeBruijnGraph *graph, size_t num_top_paths,
                       float sw_threshold, bool verbose,
                       float insertion_penalty, float deletion_penalty) :
                            graph_(graph), num_top_paths_(num_top_paths),
                            sw_threshold_(sw_threshold), verbose_(verbose),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty) {
    // Substitution loss for each pair of nucleotides.
    // Transition and transversion mutations have different loss values.
    sub_loss_ = {
        {'a', {{'a', 0}, {'t', 2}, {'c', 2}, {'g', 1}}},
        {'g', {{'g', 0}, {'t', 2}, {'c', 2}, {'a', 1}}},
        {'c', {{'c', 0}, {'a', 2}, {'g', 2}, {'t', 1}}},
        {'t', {{'t', 0}, {'a', 2}, {'g', 2}, {'c', 1}}},
        {'$', {{'$', 2}}}};
}

namespace {

template <typename T>
float std_dev (const std::vector<T>& list, T mean) {
    float std_dev = 0;
    for (auto it = list.begin(); it != list.end(); ++it)
        std_dev += std::pow((*it) - mean, 2.0);
    return std::pow(std_dev / list.size(), 0.5);
}

}  // namespace

DBGAligner::AlignedPath DBGAligner::align(const std::string& sequence) const {
    if (sequence.size() < graph_->get_k())
        return AlignedPath(sequence.end());

    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath> queue(num_top_paths_);
    queue.push(std::move(AlignedPath(std::begin(sequence))));
    while (!queue.empty()) {
        auto path = std::move(queue.top());
        queue.pop();
        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        node_index seed = path.size() > 0 ? path.back() : graph_->npos;
        graph_->extend_from_seed(path.get_sequence_it(), std::end(sequence),
                                       [&](node_index node) {
                                        // TODO: Check if this is necessary.
                                        if (node != graph_->npos)
                                            path.push_back(node, {});
                                        },
                                       [&]() { return (path.size() > 0 &&
                                                (graph_->outdegree(path.back()) > 1)); },
                                       seed);

        if (path.get_sequence_it() + graph_->get_k() > std::end(sequence))
            return path;

        inexact_map(path, std::end(sequence),
                    [&](node_index node, std::string::const_iterator last_mapped_position) {
                        AlignedPath alternative_path(path);
                        alternative_path.set_sequence_it(last_mapped_position);
                        // If the loss is less than the threshold continue with single char loss computation.
                        if (alternative_path.get_total_loss() < sw_threshold_ * sequence.size()) {
                            // TODO: Construct sequence last char more efficiently.
                            alternative_path.push_back(node, {},
                                single_char_loss(*(alternative_path.get_sequence_it() +
                                    graph_->get_k() - 1), graph_->get_node_sequence(node).back()));
                        } else {
                            alternative_path.push_back(node, {});
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

        if (verbose_ && path.size() % 50 == 0) {
            auto shadow_queue(queue);
            std::cout << "Loss mean and std of paths (top path size: "
                      << path.size() << "): ";
            std::vector<float> loss_values;
            loss_values.reserve(queue.size());
            while (!shadow_queue.empty()) {
                loss_values.push_back(shadow_queue.top().get_total_loss());
                shadow_queue.pop();
            }
            float mean = std::accumulate(loss_values.begin(), loss_values.end(), 0);
            mean /= loss_values.size();
            std::cout << mean << ", " << std_dev<float>(loss_values, mean) << std::endl;
        }
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

namespace {
// helper function to print an alignment based on SSW library.
static void PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
    std::cout << "===== SSW result =====" << std::endl;
    std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
       << "Number of mismatches:\t" << alignment.mismatches << std::endl
       << "Cigar: " << alignment.cigar_string << std::endl;
    std::cout << "======================" << std::endl;
}

}  // namespace

float DBGAligner::ssw_loss(const AlignedPath& path, std::string::const_iterator begin) const {
    auto ref = get_path_sequence(path.get_nodes());
    std::string query(begin, path.get_sequence_it() + graph_->get_k() - 1);

    int32_t maskLen = strlen(query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
    if (verbose_) {
        PrintAlignment(alignment);
    }
    return alignment.sw_score;
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
