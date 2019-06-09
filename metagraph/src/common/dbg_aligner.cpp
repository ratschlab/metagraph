#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "reverse_complement.hpp"

typedef std::pair<uint64_t, std::vector<DBGAligner::AlignedPath>::iterator> ScoredPathIt;

DBGAligner::DBGAligner(std::shared_ptr<DeBruijnGraph> graph, size_t num_top_paths,
                       size_t num_alternative_paths, uint8_t path_comparison_code,
                       bool verbose, bool discard_similar_paths,
                       float insertion_penalty, float deletion_penalty,
                       float gap_openning_penalty, float gap_extension_penalty) :
                            graph_(graph),
                            num_top_paths_(num_top_paths),
                            num_alternative_paths_(num_alternative_paths),
                            path_comparison_code_(path_comparison_code), verbose_(verbose),
                            discard_similar_paths_(discard_similar_paths),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty),
                            gap_openning_penalty_(gap_openning_penalty),
                            gap_extension_penalty_(gap_extension_penalty),
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
    cssw_aligner_.SetGapPenalty(std::abs(gap_openning_penalty_), std::abs(gap_extension_penalty_));
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

namespace { // helper functions
// Adjust nodes and sequence starting positions to start from a mapped node.
bool adjust_positions(std::vector<DBGAligner::node_index>::iterator& nodes_begin,
                      const std::vector<DBGAligner::node_index>::iterator& nodes_end,
                      std::string::const_iterator& sequence_begin,
                      std::vector<DBGAligner::node_index>::iterator& last_mapped_it,
                      std::vector<DBGAligner::node_index>::iterator& target_node_it) {
    // If nodes is all 0, return.
    auto non_zero_it = std::find_if(nodes_begin, nodes_end,
                                    [] (DBGAligner::node_index node) { return node != 0; });
    if (non_zero_it == nodes_end)
        return false;

    auto delta = non_zero_it - nodes_begin;
    nodes_begin += delta;
    sequence_begin += delta;

    last_mapped_it = std::find(nodes_begin, nodes_end, 0) - 1;
    target_node_it = std::find_if(last_mapped_it + 1, nodes_end,
                                       [] (DBGAligner::node_index node) { return node != 0; });
    return true;
}
}  // namespace

DBGAligner::AlignedPath DBGAligner::map_to_nodes(const std::string &sequence) {
    assert(sequence.size() >= k_);

    std::vector<node_index> nodes;
    graph_->map_to_nodes(sequence, [&] (node_index node) { nodes.push_back(node); });

    // TODO: perform inexact seeding.
    // TODO: enable alternative paths when calling map_to_nodes.

    auto nodes_begin = std::begin(nodes);
    auto sequence_begin = std::begin(sequence);
    std::vector<node_index>::iterator last_mapped_it, target_node_it;
    if (!adjust_positions(nodes_begin, std::end(nodes), sequence_begin, last_mapped_it, target_node_it))
        return AlignedPath(k_, std::begin(sequence), std::end(sequence), path_comparison_code_);

    // Seed the path and extend exactly as long as there is no branching points
    // and nodes are exactly mapped.
    AlignedPath stitched_path(k_, sequence_begin, sequence_begin, path_comparison_code_);
    stitched_path.seed(*nodes_begin, {}, graph_->get_node_sequence(*nodes_begin), k_ * match_score_);
    for (auto it = nodes_begin + 1; it <= last_mapped_it; ++it) {
        stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
        if (graph_->outdegree(*it) > 1) {
            last_mapped_it = it;
            break;
        }
    }

    BoundedPriorityQueue<AlignedPath> alternative_paths(num_alternative_paths_);
//    std::cerr << "map_to_nodes called for " << sequence <<std::endl;
    // Continue alignment and fill in the zeros in nodes.
    while (last_mapped_it + 1 != std::end(nodes)) {
        auto unmapped_string_begin_it = sequence_begin + (last_mapped_it - nodes_begin);
        // Explore to find alignments.
//        std::cerr << "Calling align for seq of length " << std::end(sequence) - unmapped_string_begin_it << std::endl;
//        std::cerr << "First node index " << *(unmapped_string_begin_it - sequence_begin + nodes_begin) << std::endl;
//        std::cerr << "First node index kmer" << graph_->get_node_sequence(*(unmapped_string_begin_it - sequence_begin + nodes_begin)) << std::endl;
//        std::cerr << std::string(unmapped_string_begin_it, end(sequence)) << std::endl;
        auto re_mapped_paths = align(unmapped_string_begin_it, std::end(sequence),
                                     [&] (node_index node,
                                          const std::string::const_iterator& query_it) {
                // Terminate if an exact mapped kmer with same query_it exists.
                auto index_in_original_seq = query_it - sequence_begin;
                assert(index_in_original_seq > 0 && uint64_t(index_in_original_seq) < nodes.size());
                // If a node is exactly mapped at the same index (except from the first node or seed),
                // and there is no branching point, stop exploring the graph.
                if (nodes[index_in_original_seq] == node && node != 0 &&
                    graph_->outdegree(node) <= 1 &&
                    index_in_original_seq != last_mapped_it - nodes_begin) {
                    return true;
                }
                return false;
            });

//        std::cerr << "Aligned with " << re_mapped_paths.front().get_sequence() << std::endl;
        // Stitch partial paths resulting from exploration to the final path.
        for (const AlignedPath& alternative_partial_path : re_mapped_paths) {
            // Append a path to the current result. Fills in zeros in original nodes vector.
            // Avoid adding duplicate nodes and chars in case of positive overlap.
            // Negative overlap is a sign of a gap in query.
            int64_t overlap_length = stitched_path.num_kmers_in_query()
                                     - (alternative_partial_path.get_query_begin_it() - sequence_begin);
            auto updated_score = stitched_path.get_total_score()
                                 + alternative_partial_path.get_total_score()
                                 - overlap_length * match_score_;

            if (overlap_length <= 0) {
                std::cerr << "Negative overlap when appending a re_mapped path" << std::endl;
                // Push back the stitched path until now to partial paths and set stitched_path to start from current
                // re_mapped_path.
                auto next_kmer_it = stitched_path.get_query_it() - overlap_length;

                StripedSmithWaterman::Alignment alignment;
                if (!cssw_align(stitched_path, alignment)) {
                    std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
                }
                trim(stitched_path, alignment);
                alternative_paths.push(std::move(stitched_path));

                stitched_path = AlignedPath(k_, next_kmer_it, next_kmer_it, path_comparison_code_);
                updated_score += overlap_length * match_score_;

                nodes_begin += next_kmer_it - sequence_begin;
                sequence_begin = next_kmer_it;

            }
            stitched_path.append_path(alternative_partial_path, overlap_length, updated_score);
        }
        target_node_it = nodes_begin + stitched_path.num_kmers_in_query();
        last_mapped_it = std::find(target_node_it + 1, std::end(nodes), 0) - 1;
//        std::cerr << "Stitched after align " << stitched_path.get_sequence() << std::endl;
        // If re_mapped_path wasn't able to fill in all the zeros in nodes, report a partial path
        // and continue from the next exact map.
        if (*target_node_it == 0) {
            std::cerr << "Gap after re_mapping" << std::endl;

            StripedSmithWaterman::Alignment alignment;
            if (!cssw_align(stitched_path, alignment)) {
                std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
            }
            trim(stitched_path, alignment);
//            std::cout << "stitched_path:" << stitched_path.get_sequence() << std::endl;
//            std::cout << "stitched_path it: " << stitched_path.get_query_it() - stitched_path.get_query_begin_it() << std::endl;
//            for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++ nodes_it)
//                std::cout << *nodes_it << std::endl;
            alternative_paths.push(std::move(stitched_path));

            auto non_zero_target_node_it = std::find_if(target_node_it + 1, std::end(nodes),
                                      [] (node_index node) { return node != 0; });
            if (non_zero_target_node_it == std::end(nodes))
                break;
//            std::cout << "target it in nodes vector: " << target_node_it - nodes_begin << std::endl;
//            std::cout << "non_zero target it in nodes vector: " << non_zero_target_node_it - nodes_begin << std::endl;
            auto next_kmer_it = non_zero_target_node_it - nodes_begin + sequence_begin;
//            std::cout << "next_kmer ind: " << next_kmer_it - sequence_begin << std::endl;
            stitched_path = AlignedPath(k_, next_kmer_it, next_kmer_it, path_comparison_code_);
            stitched_path.seed(*non_zero_target_node_it, {},
                               graph_->get_node_sequence(*non_zero_target_node_it),
                               k_ * match_score_);

            nodes_begin += next_kmer_it - sequence_begin;
            sequence_begin = next_kmer_it;
            // first exact mapped kmer is seeded;
            target_node_it = non_zero_target_node_it + 1;
        }

        // Append the exactly mapped regions to the final path from nodes vector.
        for (auto it = target_node_it; it <= last_mapped_it; ++it) {
            stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
            if (graph_->outdegree(*it) > 1) {
                last_mapped_it = it;
                break;
            }
        }

//        std::cerr << "stitched after exact nodes " << stitched_path.get_sequence() << std::endl;
        // Update last_mapped_it, target_node_it to point to the beginning and end of the inexactly mapped region.
        target_node_it = std::find_if(last_mapped_it + 1, std::end(nodes),
                                      [] (node_index node) { return node != 0; });
    }
    StripedSmithWaterman::Alignment alignment;
    if (!cssw_align(stitched_path, alignment)) {
        std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
    }
    trim(stitched_path, alignment);
    alternative_paths.push(std::move(stitched_path));
    return alternative_paths.top();
}

std::vector<DBGAligner::AlignedPath> DBGAligner::align(const std::string::const_iterator &sequence_begin, const std::string::const_iterator &sequence_end,
                                                       const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    if (sequence_begin + k_ > sequence_end)
        return std::vector<AlignedPath>();

    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    std::vector<AlignedPath> alternative_paths;

    BoundedPriorityQueue<AlignedPath> queue(num_top_paths_);
    AlignedPath path(k_, sequence_begin, sequence_begin, path_comparison_code_);
    queue.push(path);

    while (!queue.empty() && alternative_paths.size() < num_alternative_paths_) {
        path = std::move(queue.top());
        queue.pop();
//        std::cerr << "popping path (#matches, size, seq) : " << path.get_num_matches()
//                  << " " << path.size() + k_ - 1
//                  << " " << path.get_sequence() << std::endl;
        // Early cutoff if the top path has poor alignment.
//        if (path.size() > 10 && 2 * path.get_num_matches() < path.size() + k_ - 1) {
//            alternative_paths.push_back(path);
//            break;
//        }

        if (path.get_query_it() + k_ > sequence_end) {
            if (path.size() != 0)
                alternative_paths.push_back(path);
            continue;
        }
//        std::cerr << "calling exact map" << std::endl;

        if (!exact_map(path, sequence_end, dp_alignment, terminate)) {
            if (path.size() != 0)
                alternative_paths.push_back(path);
            break;
        }

        if (path.get_query_it() + k_ > sequence_end) {
            if (path.size() != 0)
                alternative_paths.push_back(path);
            continue;
        }

        if (path.get_similar())
            continue;

//        std::cerr << "calling inexact map" << std::endl;

        if (!inexact_map(path, queue, dp_alignment, terminate) || queue.empty()) {
            if (path.size() != 0)
                alternative_paths.push_back(path);
            break;
        }
    }
    
    for (auto alternative_path_it = std::begin(alternative_paths);
            alternative_path_it != std::end(alternative_paths); ++ alternative_path_it) {
        StripedSmithWaterman::Alignment alignment;
        if (!cssw_align(*alternative_path_it, alignment)) {
            std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
            break;
        }
        trim(*alternative_path_it, alignment);
    }
    return alternative_paths;
}

bool DBGAligner::exact_map(AlignedPath &path, const std::string::const_iterator &sequence_end,
                           std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                           const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    node_index seed = path.size() > 0 ? path.back() : graph_->npos;
    bool continue_alignment = true;
    graph_->extend_from_seed(path.get_query_it(), sequence_end,
            [&](node_index node) {
            if (node == graph_->npos)
                return;
            if (path.size() > 0) {
                if (terminate(node, path.get_query_it())) {
                    continue_alignment = false;
                    return;
                }
//                std::cerr << "Exact extension, terminate: " << !continue_alignment << std::endl;
                path.extend(node, {},
                            *(path.get_query_it() + k_ - 1), match_score_);
            } else {
                path.seed(node, {}, graph_->get_node_sequence(node),
                          match_score_ * k_);
            }

            // Merge Similar paths.
            if (discard_similar_paths_) {
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
                    path.set_similar();
                }
            }
            },
            [&]() {
//                    std::cerr << "path.get_similar(): " << path.get_similar()
//                              << ", !continue_alignment: " << !continue_alignment
//                              << ", outdegree > 1: " << (path.size() > 0 && (graph_->outdegree(path.back()) > 1)) << std::endl;
                    return (path.get_similar() || !continue_alignment ||
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
        queue.push(path);
        return continue_alignment;
    }
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char extension) {
        AlignedPath alternative_path(path);
        alternative_path.set_query_it(path.get_query_it());
        // Stop early in case of early termination condition by the caller function.
        if (!continue_alignment || terminate(node, alternative_path.get_query_it())) {
            continue_alignment = false;
            return;
        }
        // Approximate Smith-Waterman score by assuming either match/mismatch for the next character.
        alternative_path.extend(node, {}, extension,
            single_char_score(*(alternative_path.get_query_it() +
                k_ - 1), extension));
        // Merge Similar paths.
        if (discard_similar_paths_) {
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
        }
        else {
            queue.push(std::move(alternative_path));
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

bool DBGAligner::cssw_align(AlignedPath &path, StripedSmithWaterman::Alignment& alignment) const {
    if (path.is_score_updated()) {
        alignment = path.get_alignment();
        return true;
    }

    auto ref = path.get_sequence();
    std::string query(path.get_query_begin_it(), path.get_query_it() + k_ - 1);

    int32_t maskLen = strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    StripedSmithWaterman::Filter filter;
    if (!cssw_aligner_.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen))
        return false;

    path.update_alignment(alignment);

    return true;
}

//void DBGAligner::whole_path_score(const AlignedPath &path, GlobalSW& global_sw) const {
//    auto path_sequence = path.get_sequence();
//    std::string query_sequence(path.get_query_begin_it(), path.get_query_it() + k_ - 1);
//    SWDpCell alignment_score[query_sequence.size() + 1][path_sequence.size() + 1] = {};
//    for (size_t i = 1; i <= query_sequence.size(); ++i) {
//        for (size_t j = 1; j <= path_sequence.size(); ++j) {
//            std::vector<SWDpCell> options = {
//                {.value = alignment_score[i-1][j-1] + single_char_score(query_sequence[i-1], path_sequence[j-1]),
//                 .parent_row = i-1,
//                 .parent_column = j-1},
//                {.value = alignment_score[i][j-1] - insertion_penalty_,
//                 .parent_row = i,
//                 .parent_column = j-1},
//                {.value = alignment_score[i-1][j] - deletion_penalty_,
//                 .parent_row = i-1,
//                 .parent_column = j};
//            auto max_valued_option = std::max(options);
//            alignment_score[i][j] = {.value = max_valued_option.value,
//                                     .parent_row = max_valued_option.parent_row,
//                                     .parent_column = max_valued_option.parent_column};
//        }
//    }
//    auto last_i = query_sequence.size();
//    auto last_j = path_sequence.size();
//    auto optimum_alignment_cell = alignment_score[last_i][last_j];
//    global_sw.value = optimum_alignment_cell.value;
//    while (optimum_alignment_cell.parent_row >= 0 &&
//           optimum_alignment_cell.parent_column >= 0) {
//        optimum_alignment_cell =
//            alignment_score[optimum_alignment_cell.parent_row]
//                           [optimum_alignment_cell.parent_column];
//    }
//}

// Note: Assume Cigar is updated.
void DBGAligner::trim(AlignedPath &path, const StripedSmithWaterman::Alignment& alignment) const {
    auto cigar = path.get_cigar();
    // Trim path if reference end in cigar is not according to reference end in path.
    uint64_t path_trim_length = path.get_sequence().size() - alignment.ref_end - 1;
    // Trim query if marked by S in cigar.
    uint64_t query_trim_length = 0;
    if (cigar.back() == 'S') {
        std::smatch match;
        std::regex_search(cigar, match, std::regex("[0-9]+S$"));
        std::string position_str = match.str();
        position_str.pop_back();
        query_trim_length += atoi(position_str.c_str());
    }
    path.trim(query_trim_length, path_trim_length);
}
