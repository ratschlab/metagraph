#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "reverse_complement.hpp"

#define NUMERICAL_DIGITS "0123456789"


DBGAligner::DBGAligner(DeBruijnGraph *graph,
                       size_t num_top_paths, bool verbose,
                       float sw_threshold, float re_seeding_threshold,
                       float insertion_penalty, float deletion_penalty) :
                            graph_(graph),
                            num_top_paths_(num_top_paths),
                            verbose_(verbose),
                            sw_threshold_(sw_threshold),
                            re_seeding_threshold_(re_seeding_threshold),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty),
                            merged_paths_counter_(0) {
    match_score_ = 2;
    k_ = graph_->get_k();
    // Substitution score for each pair of nucleotides.
    // Transition and transversion mutations have different score values.
    sub_score_ = {
        {'a', {{'a', match_score_}, {'t', -2}, {'c', -2}, {'g', -1}}},
        {'g', {{'g', match_score_}, {'t', -2}, {'c', -2}, {'a', -1}}},
        {'c', {{'c', match_score_}, {'a', -2}, {'g', -2}, {'t', -1}}},
        {'t', {{'t', match_score_}, {'a', -2}, {'g', -2}, {'c', -1}}},
        {'$', {{'$', -2}}}};

    // score for the following chars respectively:
    // A  C  G  T   N (or other ambiguous code)
    int8_t cssw_score_matrix_[5*5] = {
    sub_score_['a']['a'], sub_score_['a']['c'], sub_score_['a']['g'], sub_score_['a']['t'], sub_score_['$']['$'],
    sub_score_['c']['a'], sub_score_['c']['c'], sub_score_['c']['g'], sub_score_['c']['t'], sub_score_['$']['$'],
    sub_score_['g']['a'], sub_score_['g']['c'], sub_score_['g']['g'], sub_score_['g']['t'], sub_score_['$']['$'],
    sub_score_['t']['a'], sub_score_['t']['c'], sub_score_['t']['g'], sub_score_['t']['t'], sub_score_['$']['$'],
    sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$']};

    int8_t cssw_translation_matrix_[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4};

    cssw_aligner_.ReBuild(cssw_score_matrix_, 5, cssw_translation_matrix_, 128);
    cssw_aligner_.SetGapPenalty(std::abs(deletion_penalty_), std::abs(deletion_penalty_));
}

namespace { // Helper function.
template <typename T>
float std_dev(const std::vector<T> &list, T mean) {
    float std_dev = 0;
    for (auto it = list.begin(); it != list.end(); ++it)
        std_dev += std::pow((*it) - mean, 2.0);
    return std::pow(std_dev / list.size(), 0.5);
}
}  // namespace

DBGAligner::AlignedPath DBGAligner::map_to_nodes_forward_reverse_complement(const std::string &sequence) {
    auto reverse_complement_sequence(sequence);
    reverse_complement(std::begin(reverse_complement_sequence), std::end(reverse_complement_sequence));
    auto path = map_to_nodes(sequence);
    auto reverse_path = map_to_nodes(reverse_complement_sequence);
    return path.get_total_score() > reverse_path.get_total_score() ? path : reverse_path;
}

DBGAligner::AlignedPath DBGAligner::map_to_nodes(const std::string &sequence) {
    assert(sequence.size() >= k_);

    std::vector<node_index> nodes;
    graph_->map_to_nodes(sequence, [&] (node_index node) { nodes.push_back(node); });

    if (verbose_) {
        std::cout << "Aligning read to reference graph using map_to_nodes." << std::endl;
        std::cout << "Number of zeros in initial mapping: "
                  << std::count(std::begin(nodes), std::end(nodes), 0) << std::endl;
    }

    // TODO: perform inexact seeding.

    // If nodes is all 0, return.
    auto non_zero_it = std::find_if(std::begin(nodes), std::end(nodes),
                                    [] (node_index node) { return node != 0; });
    if (non_zero_it == std::end(nodes))
        return AlignedPath(k_, std::begin(sequence), std::end(sequence));
    // Adjust nodes and sequence starting positions.
    auto delta = non_zero_it - std::begin(nodes);
    auto nodes_begin = std::begin(nodes) + delta;
    auto sequence_begin = std::begin(sequence) + delta;

    auto last_mapped_it = std::find(nodes_begin, std::end(nodes), 0) - 1;
    auto target_node_it = std::find_if(last_mapped_it + 1, std::end(nodes),
                                       [] (node_index node) { return node != 0; });
    // Seed the path and extend exactly as much as possible.
    AlignedPath stitched_path(k_, sequence_begin, sequence_begin);
    stitched_path.seed(*nodes_begin, {}, graph_->get_node_sequence(*nodes_begin), k_ * match_score_);
    for (auto it = nodes_begin + 1; it <= last_mapped_it; ++it) {
        stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
        if (graph_->outdegree(*it) > 1) {
            last_mapped_it = it;
            break;
        }
    }
    // Aligned kmers with no exact map.
    while (last_mapped_it + 1 != std::end(nodes)) {
        auto unmapped_string_begin_it = sequence_begin + (last_mapped_it - nodes_begin);
        auto unmapped_string = std::string(unmapped_string_begin_it, std::end(sequence));

        auto re_mapped_paths = align(unmapped_string, [&] (node_index node, const std::string::const_iterator& query_it) {
                // Terminate if an exact mapped kmer with same query_it exists.
                auto index_in_original_seq = query_it - std::begin(unmapped_string) +
                                             unmapped_string_begin_it - sequence_begin;
                assert(index_in_original_seq > 0 && index_in_original_seq < nodes.size());
                // if a node is exactly mapped at the same index (except from the first node or seed), stop exploring the graph.
                if (nodes[index_in_original_seq] == node && node != 0 && graph_->outdegree(node) <= 1 &&
                    index_in_original_seq != last_mapped_it - nodes_begin) {
                    return true;
                }
                return false;
            });
        for (const AlignedPath& partial_path : re_mapped_paths) {
            // Append a path to the current object. Fills in 0 in case of spaced nodes.
            // Avoid adding duplicate nodes and chars in case of positive overlap.
            int64_t overlap_length = stitched_path.num_kmers_in_query()
                                     - (partial_path.get_query_begin_it() - std::begin(unmapped_string)
                                        + unmapped_string_begin_it - std::begin(sequence));
            auto updated_score = overlap_length <= 0 ? 0 :
                                 stitched_path.get_total_score() + partial_path.get_total_score() - overlap_length * match_score_;
            stitched_path.append_path(partial_path, overlap_length, updated_score);
        }
        target_node_it = nodes_begin + stitched_path.num_kmers_in_query();
        last_mapped_it = std::find(target_node_it + 1, std::end(nodes), 0) - 1;

        // Append the exactly mapped regions to the final path.
        for (auto it = target_node_it; it <= last_mapped_it; ++it) {
            stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
            if (graph_->outdegree(*it) > 1) {
                last_mapped_it = it;
                break;
            }
        }
        target_node_it = std::find_if(last_mapped_it + 1, std::end(nodes),
                                      [] (node_index node) { return node != 0; });
    }
    return stitched_path;
}

std::vector<DBGAligner::AlignedPath> DBGAligner::align(const std::string &sequence,
                                                       const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    if (sequence.size() < k_)
        return std::vector<AlignedPath>();

    std::vector<AlignedPath> partial_paths;
    uint32_t print_frequency = 50;
    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath> queue(num_top_paths_);
    AlignedPath initial_path(k_, std::begin(sequence), std::begin(sequence));
    queue.push(std::move(initial_path));
    while (!queue.empty()) {
        auto path = std::move(queue.top());
        queue.pop();

        if (path.get_total_score() < 0) {
            std::cerr << "Trimming and storing partial path\n." << std::endl;
            partial_paths.push_back(path);
            queue = BoundedPriorityQueue<AlignedPath>(num_top_paths_);
            dp_alignment.clear();
            path = AlignedPath(k_, path.get_query_it() + 1, path.get_query_it() + 1);
        }

        if (path.get_query_it() + k_ > std::end(sequence)) {
            partial_paths.push_back(path);
            break;
        }

        if (!exact_map(path, sequence, dp_alignment, terminate) || path.get_query_it() + k_ > std::end(sequence)) {
            partial_paths.push_back(path);
            break;
        }

        if (path.get_similar())
            continue;

        if (!inexact_map(path, queue, dp_alignment, terminate) || queue.empty()) {
            partial_paths.push_back(path);
            break;
        }

        if (verbose_ && path.size() % print_frequency == 0) {
            auto shadow_queue(queue);
            std::cout << "Loss mean and std of paths (top path size: "
                      << path.size() << "): ";
            std::vector<float> score_values;
            score_values.reserve(queue.size());
            while (!shadow_queue.empty()) {
                score_values.push_back(shadow_queue.top().get_total_score());
                shadow_queue.pop();
            }
            float mean = std::accumulate(score_values.begin(), score_values.end(), 0);
            mean /= score_values.size();
            std::cout << mean << ", " << std_dev<float>(score_values, mean) << std::endl;
        }
    }
    if (partial_paths.empty())
        return std::vector<AlignedPath>();

    std::vector<AlignedPath> complete_path;
    for (auto partial_path_it = std::begin(partial_paths);
         partial_path_it != std::end(partial_paths); ++ partial_path_it) {
        // TODO: Find shortest path between the end node of a partial path
        // and the first node of the next partial path.
        StripedSmithWaterman::Alignment alignment;
        if (!cssw_align(*partial_path_it, alignment)) {
            std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
            break;
        }
        partial_path_it->set_cigar(alignment.cigar_string);
        trim(*partial_path_it);
        complete_path.push_back(*partial_path_it);
    }

    // TODO: Stitch partial paths together and return the corresponding path.
    // Insertions is the parts in query that aren't matched and for deletion,
    // we have to find approximation or exact distance using labeling schemes.
    return complete_path;
}

bool DBGAligner::exact_map(AlignedPath &path, const std::string &sequence,
                           std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                           const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    node_index seed = path.size() > 0 ? path.back() : graph_->npos;
    bool continue_alignment = true;
    graph_->extend_from_seed(path.get_query_it(), std::end(sequence),
            [&](node_index node) {
            if (node == graph_->npos)
                return;
            if (path.size() > 0) {
                if (terminate(node, path.get_query_it())) {
                    continue_alignment = false;
                    return;
                }
                path.extend(node, {},
                            *(path.get_query_it() + k_ - 1), match_score_);
            } else {
                path.seed(node, {}, graph_->get_node_sequence(node),
                          match_score_ * k_);
            }

        // Merge Similar paths.
        DPAlignmentKey key{.node = path.back(),
                           .query_begin_it = path.get_query_begin_it(),
                           .query_it = path.get_query_it()};
        DPAlignmentValue value = {.score = path.get_total_score()};
        auto dp_alignment_it = dp_alignment.find(key);
        if (dp_alignment_it == dp_alignment.end()) {
            dp_alignment[key] = value;
        }
        else if (dp_alignment_it->second.score < path.get_total_score()) {
            dp_alignment_it->second = value;
        }
        else {
            merged_paths_counter_ ++;
            std::cerr << "Similar path: " << path.get_sequence() << " detected.\n";
            path.set_similar();
        }},
            [&]() { return (path.get_similar() || !continue_alignment ||
                        (path.size() > 0 && (graph_->outdegree(path.back()) > 1))); },
            seed);
    return continue_alignment;
}

bool DBGAligner::inexact_map(AlignedPath& path,
                             BoundedPriorityQueue<AlignedPath> &queue,
                             std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                             const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    bool continue_alignment = true;
    if (path.size() == 0) {
        // To do inexact seeding, a path with incremented query_it is pushed to the queue.
        path.set_query_begin_it(path.get_query_it() + 1);
        path.set_query_it(path.get_query_it() + 1);
        queue.push(std::move(path));
        return continue_alignment;
    }
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char extension) {
        AlignedPath alternative_path(path);
        alternative_path.set_query_it(path.get_query_it());

        if (!continue_alignment || terminate(node, alternative_path.get_query_it())) {
            continue_alignment = false;
            return;
        }

        if (alternative_path.get_total_score() > sw_threshold_ * path.size()) {
            alternative_path.extend(node, {}, extension,
                single_char_score(*(alternative_path.get_query_it() +
                    k_ - 1), extension));
        } else {
            alternative_path.extend(node, {}, extension);
            StripedSmithWaterman::Alignment alignment;
            if (!cssw_align(alternative_path, alignment))
                std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
            alternative_path.update_total_score(alignment.sw_score);

            alternative_path.set_cigar(alignment.cigar_string);
        }
        // Merge Similar paths.
        DPAlignmentKey alternative_key{.node = alternative_path.back(),
                           .query_begin_it = alternative_path.get_query_begin_it(),
                           .query_it = alternative_path.get_query_it()};
        DPAlignmentValue alternative_value = {.score = alternative_path.get_total_score()};
        auto dp_alignment_it = dp_alignment.find(alternative_key);
        if (dp_alignment_it == dp_alignment.end()) {
            dp_alignment[alternative_key] = alternative_value;
            queue.push(std::move(alternative_path));
        }
        else if (dp_alignment_it->second.score < alternative_path.get_total_score()) {
            dp_alignment_it->second = alternative_value;
            queue.push(std::move(alternative_path));
        }
        else {
            merged_paths_counter_ ++;
            path.set_similar();
        }
        });
    return continue_alignment;
}

float DBGAligner::single_char_score(char char_in_query, char char_in_graph) const {
    try {
        return sub_score_.at(std::tolower(char_in_query)).at(std::tolower(char_in_graph));
    }
    catch (const std::out_of_range&) {
        return sub_score_.at('$').at('$');
    }
}

namespace {
// helper function to print an alignment based on SSW library.
static void PrintAlignment(const StripedSmithWaterman::Alignment &alignment){
    std::cout << "===== SSW result =====" << std::endl;
    std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
       << "Reference start:\t" << alignment.ref_begin << std::endl
       << "Reference end:\t" << alignment.ref_end << std::endl
       << "Query start:\t" << alignment.query_begin << std::endl
       << "Query end:\t" << alignment.query_end << std::endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
       << "Number of mismatches:\t" << alignment.mismatches << std::endl
       << "Cigar: " << alignment.cigar_string << std::endl;
    std::cout << "======================" << std::endl;
}
}  // namespace

bool DBGAligner::cssw_align(const AlignedPath &path, StripedSmithWaterman::Alignment& alignment) const {
    auto ref = path.get_sequence();
    std::string query(path.get_query_begin_it(), path.get_query_it() + k_ - 1);

    int32_t maskLen = strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    StripedSmithWaterman::Filter filter;
    return cssw_aligner_.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
}

float DBGAligner::ssw_score(const AlignedPath &path) const {
    StripedSmithWaterman::Alignment alignment;
    if (verbose_) {
        PrintAlignment(alignment);
    }
    if (!cssw_align(path, alignment)) {
        std::cout << "Failure in SSW calculation. Computing SSW inefficiently.\n";
        return whole_path_score(path);
    }
    return alignment.sw_score;
}

float DBGAligner::whole_path_score(const AlignedPath &path) const {
    auto path_sequence = path.get_sequence();
    std::string query_sequence(path.get_query_begin_it(), path.get_query_it() + k_ - 1);
    float alignment_score[query_sequence.size() + 1][path_sequence.size() + 1] = {};
    for (size_t i = 1; i <= query_sequence.size(); ++i) {
        for (size_t j = 1; j <= path_sequence.size(); ++j) {
            alignment_score[i][j] = std::max({alignment_score[i-1][j-1] +
                    single_char_score(query_sequence[i-1], path_sequence[j-1]),
                alignment_score[i][j-1] - insertion_penalty_, alignment_score[i-1][j] - deletion_penalty_});
        }
    }
    return alignment_score[query_sequence.size()][path_sequence.size()];
}

void DBGAligner::trim(AlignedPath &path) const {
    StripedSmithWaterman::Alignment alignment;
    if (!cssw_align(path, alignment)) {
        std::cout << "Failure in SSW calculation. Computing SSW inefficiently.\n";
        return;
    }
    if (verbose_) {
        PrintAlignment(alignment);
    }
    // Returns in case of no clipping in CSSW cigar.
    if (alignment.cigar_string.back() != 'S')
        return;
    auto num_starting_pos = alignment.cigar_string.find_last_not_of(NUMERICAL_DIGITS,
                                alignment.cigar_string.size() - 2) + 1;
    auto trim_length = std::atoi(alignment.cigar_string.substr(num_starting_pos,
                                 alignment.cigar_string.size() - num_starting_pos - 1).c_str());
    path.trim_with_score(trim_length, alignment.sw_score);
    path.set_cigar(alignment.cigar_string);
}
