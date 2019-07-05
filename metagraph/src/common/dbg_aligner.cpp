#include "dbg_aligner.hpp"

#include <string>
#include <random>

#include "reverse_complement.hpp"

typedef std::pair<uint64_t, std::vector<DBGAligner::AlignedPath>::iterator> ScoredPathIt;

std::function<bool(const DBGAligner::AlignedPath&, const DBGAligner::AlignedPath&)> DBGAligner::priority_function_;

DBGAligner::DBGAligner(std::shared_ptr<DeBruijnGraph> graph, DBGAlignerConfig dbg_aligner_config) :
                            graph_(graph),
                            num_top_paths_(dbg_aligner_config.num_top_paths),
                            num_alternative_paths_(dbg_aligner_config.num_alternative_paths),
                            path_comparison_code_(dbg_aligner_config.path_comparison_code),
                            verbose_(dbg_aligner_config.verbose),
                            discard_similar_paths_(dbg_aligner_config.discard_similar_paths),
                            use_cssw_lib_(dbg_aligner_config.use_cssw_lib),
                            sw_threshold_for_stitched_path_(dbg_aligner_config.sw_threshold_for_stitched_path),
                            insertion_penalty_(dbg_aligner_config.insertion_penalty),
                            deletion_penalty_(dbg_aligner_config.deletion_penalty),
                            gap_opening_penalty_(dbg_aligner_config.gap_opening_penalty),
                            gap_extension_penalty_(dbg_aligner_config.gap_extension_penalty),
                            merged_paths_counter_(0), explored_paths_counter_(0), explored_seeds_counter_(0) {
    k_ = graph_->get_k();
    match_score_ = 2;
    mm_transition_ = -1;
    mm_transversion_ = -2;
    sw_prev_row_ = new std::vector<SWDpCell>(100);
    sw_cur_row_ = new std::vector<SWDpCell>(100);
    // Substitution score for each pair of nucleotides.
    // Transition and transversion mutations have different score values.
    fill_score_matrix(128*128);

    assert(path_comparison_code_ >= 0 && path_comparison_code_ <= 2);
    switch (path_comparison_code_) {
        case 0:
            DBGAligner::priority_function_ = AlignedPath::total_score_less_than;
            break;
        case 1:
            DBGAligner::priority_function_ = AlignedPath::normalized_score_less_than;
            break;
        case 2:
            DBGAligner::priority_function_ = AlignedPath::num_matches_less_than;
            break;
        default:
            std::cout << "Warning: Priority function is not set properly in path construction."
                      << std::endl << "Setting it to default (total score comparison)." << std::endl;
            DBGAligner::priority_function_ = AlignedPath::total_score_less_than;
    }

    // score for the following chars respectively:
    // A  C  G  T   N (or other ambiguous code)
    int8_t cssw_score_matrix_[5*5] = {
    score_matrix_['a']['a'], score_matrix_['a']['c'], score_matrix_['a']['g'], score_matrix_['a']['t'], score_matrix_['$']['$'],
    score_matrix_['c']['a'], score_matrix_['c']['c'], score_matrix_['c']['g'], score_matrix_['c']['t'], score_matrix_['$']['$'],
    score_matrix_['g']['a'], score_matrix_['g']['c'], score_matrix_['g']['g'], score_matrix_['g']['t'], score_matrix_['$']['$'],
    score_matrix_['t']['a'], score_matrix_['t']['c'], score_matrix_['t']['g'], score_matrix_['t']['t'], score_matrix_['$']['$'],
    score_matrix_['$']['$'], score_matrix_['$']['$'], score_matrix_['$']['$'], score_matrix_['$']['$'], score_matrix_['$']['$']};

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
    cssw_aligner_.SetGapPenalty(std::abs(gap_opening_penalty_), std::abs(gap_extension_penalty_));
}

DBGAligner::AlignedPath DBGAligner::map_to_nodes_forward_reverse_complement(const std::string &sequence) {
    auto reverse_complement_sequence(sequence);
    reverse_complement(std::begin(reverse_complement_sequence),
                       std::end(reverse_complement_sequence));
    std::vector<AlignedPath> alternative_paths;
    AlignedPath path(k_, std::begin(sequence), std::end(sequence));
    if (map_to_nodes(sequence, alternative_paths))
        path = std::move(alternative_paths.front());
    alternative_paths.clear();
    AlignedPath reverse_path(k_, std::begin(sequence), std::end(sequence));
    if  (map_to_nodes(reverse_complement_sequence, alternative_paths, path.get_total_score()))
        reverse_path = std::move(alternative_paths.front());
    return path.get_total_score() > reverse_path.get_total_score() ? path : reverse_path;
}

bool DBGAligner::map_to_nodes(const std::string &sequence,
                              std::vector<DBGAligner::AlignedPath>& alternative_paths_vec,
                              int64_t alternative_alignment_score) {
    assert(sequence.size() >= k_);

    std::vector<node_index> nodes;
    graph_->map_to_nodes(sequence, [&] (node_index node) { nodes.push_back(node); });

    auto nodes_begin = std::begin(nodes);
    auto sequence_begin = std::begin(sequence);
    auto nodes_end = std::end(nodes);
    std::vector<node_index>::iterator last_mapped_it, target_node_it;
    last_mapped_it = (*nodes_begin) == 0 ? nodes_begin : std::find(nodes_begin, nodes_end, 0) - 1;
    target_node_it = std::find_if(last_mapped_it + 1, nodes_end,
                                       [] (DBGAligner::node_index node) { return node != 0; });

    // If no kmer has mapped.
    if (std::find_if(nodes_begin, nodes_end, [] (DBGAligner::node_index node) { return node != 0; }) == nodes_end) {
        std::cerr << "No kmer not mapped for " << sequence << std::endl;
    }
    // Seed the path and extend exactly as long as there is no branching points
    // and nodes are exactly mapped.
    AlignedPath stitched_path(k_, sequence_begin, sequence_begin);
    if (*last_mapped_it != 0) {
        stitched_path.seed(*nodes_begin, {}, graph_->get_node_sequence(*nodes_begin), k_ * match_score_);
        for (auto it = nodes_begin + 1; it <= last_mapped_it; ++it) {
            stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
            if (graph_->outdegree(*it) > 1) {
                last_mapped_it = it;
                break;
            }
        }
    }
    else {
        std::cerr << "First kmer not mapped for " << sequence << std::endl;
    }

    BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>,
                         decltype(&AlignedPath::total_score_less_than)>
                            alternative_paths(AlignedPath::total_score_less_than, num_alternative_paths_);
//    std::cerr << "map_to_nodes called for " << sequence << std::endl;
    // Continue alignment and fill in the zeros in nodes.
    while (last_mapped_it + 1 < nodes_end) {
        auto unmapped_string_begin_it = sequence_begin + (last_mapped_it - nodes_begin);
        // Explore to find alignments.
//        std::cerr << "Calling align for " << std::string(unmapped_string_begin_it, end(sequence)) << std::endl;
        auto re_mapped_paths = align_by_graph_exploration(unmapped_string_begin_it, std::end(sequence),
                                     [&] (AlignedPath& path) {
                return (stitched_path.get_total_score() + path.get_total_score()
                    + (std::end(sequence) - path.get_query_it()) * match_score_ < alternative_alignment_score);
                },
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
        // In case no path was found while exploring the graph.
        // This can happen if all paths result in low scoring paths compared to alternative_alignment_score
        // which is the score from either forward or reverse_complement of the sequence. Or in the case that
        // seeding was not possible.
        if (re_mapped_paths.size() == 0) {
            //std::cerr << "Terminating early because of poor alignment compared to forward sequence" << std::endl;
            break;
        }
        auto path_before_graph_exploration(stitched_path);
        // Stitch partial paths resulting from exploration to the final path.
        for (const AlignedPath& alternative_partial_path : re_mapped_paths) {
//            std::cerr << "Aligned with " << alternative_partial_path.get_sequence()
//                      << " cigar " << alternative_partial_path.get_cigar() << std::endl;
            // Append a path to the current result. Fills in zeros in original nodes vector.
            // Avoid adding duplicate nodes and chars in case of positive overlap.
            unsigned stitched_path_offset = stitched_path.get_query_begin_it() - sequence_begin;
            int64_t overlap_length = stitched_path.num_kmers_in_query() + stitched_path_offset
                                     - (alternative_partial_path.get_query_begin_it() - sequence_begin);
            int64_t updated_score = stitched_path.get_total_score()
                                 + alternative_partial_path.get_total_score()
                                 - (overlap_length + k_ - 1)  * match_score_;
            // This case shouldn't normally happen.
            if (overlap_length < 0) {
                stitched_path = AlignedPath(k_, sequence_begin, sequence_begin);
                std::cerr << "Negative overlap when appending a re_mapped path" << std::endl;
                break;
            }
            stitched_path.append_path(alternative_partial_path, overlap_length, updated_score);
            target_node_it = nodes_begin + stitched_path.num_kmers_in_query() + stitched_path_offset;
            last_mapped_it = std::find(target_node_it + 1, nodes_end, 0) - 1;
            // If re_mapped_path wasn't able to fill in all the zeros in nodes, report a partial path
            // and continue from the next exact map.
            if (stitched_path.size() + alternative_partial_path.size() < sw_threshold_for_stitched_path_) {
                smith_waterman_score(stitched_path);
                trim(stitched_path);
            }
            // In case all zeros were filled with the graph exploration call, continue with a single nice long
            // path and discard alternative alignments for the zero region for now.
            if (*target_node_it != 0)
                break;
            // If all zeros were not filled with graph exploration call, report stitched path so far
            // as an alternative path.
//            std::cerr << "Pushing stitched path" << std::endl;
            alternative_paths.push(std::move(stitched_path));
//            std::cerr << "Pushed stitched path" << std::endl;
            stitched_path = path_before_graph_exploration;
        }
        if (stitched_path.size() == 0) {
            std::cerr << "Terminating early because of negative overlap when appending a re_mapped path" << std::endl;
            break;
        }
        // Skip the remainder of the zero region and adjust last_mapped_it and target_node_it accordingly.
        // Reset stitched_path to start at the first exactly mapped after the current kmer.
        if (*target_node_it == 0) {
            auto non_zero_target_node_it = std::find_if(target_node_it + 1, nodes_end,
                                      [] (node_index node) { return node != 0; });
            if (non_zero_target_node_it == nodes_end)
                break;
            auto next_kmer_it = non_zero_target_node_it - nodes_begin + sequence_begin;
            // first exact mapped kmer is added to the path as seed.
            stitched_path = AlignedPath(k_, next_kmer_it, next_kmer_it);
            stitched_path.seed(*non_zero_target_node_it, {},
                               graph_->get_node_sequence(*non_zero_target_node_it),
                               k_ * match_score_);

            target_node_it = std::find_if(non_zero_target_node_it + 1, nodes_end,
                                      [] (node_index node) { return node != 0; });
            last_mapped_it = std::find(std::min(target_node_it + 1, nodes_end), nodes_end, 0) - 1;
        }
//        std::cerr << "Stitched after filling remaining gaps from re-mapped " << stitched_path.get_sequence()
//                      << " cigar " << stitched_path.get_cigar() << std::endl;
        if (last_mapped_it == nodes_end || target_node_it == nodes_end)
            break;
//        std::cerr << "appending exactly_mapped regions" << std::endl;
        // Append the exactly mapped regions to the final path from nodes vector.
        for (auto it = target_node_it; it <= last_mapped_it; ++it) {
            stitched_path.extend(*it, {}, *(sequence_begin + (it - nodes_begin) + k_ - 1), match_score_);
            if (graph_->outdegree(*it) > 1) {
                last_mapped_it = it;
                break;
            }
        }
//        std::cerr << "stitched after exact nodes " << stitched_path.get_sequence()
//                      << " cigar " << stitched_path.get_cigar() << std::endl;
        // Update last_mapped_it, target_node_it to point to the beginning and end of the inexactly mapped region.
        target_node_it = std::find_if(last_mapped_it + 1, nodes_end,
                                      [] (node_index node) { return node != 0; });
//        std::cerr << "target it in nodes vector: " << target_node_it - nodes_begin << std::endl;
    }

    if (stitched_path.size() == 0)
        return false;

    if (stitched_path.size() < sw_threshold_for_stitched_path_) {
        smith_waterman_score(stitched_path);
        trim(stitched_path);
    }
    alternative_paths.push(std::move(stitched_path));
//    std::cerr << "Pushed stitched path" << std::endl;

    while(!alternative_paths.empty() && alternative_paths_vec.size() < num_alternative_paths_) {
        alternative_paths_vec.push_back(std::move(alternative_paths.top()));
        alternative_paths.pop();
    }
    return alternative_paths_vec.size() > 0;
}

std::vector<DBGAligner::AlignedPath> DBGAligner::align_by_graph_exploration(
            const std::string::const_iterator &sequence_begin, const std::string::const_iterator &sequence_end,
            const std::function<bool (AlignedPath&)>& early_discard_path,
            const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate_mapping) {

    if (sequence_begin + k_ > sequence_end)
        return std::vector<AlignedPath>();

    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>,
                         decltype(&AlignedPath::total_score_less_than)>
                            alternative_paths(AlignedPath::total_score_less_than, num_alternative_paths_);

    BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>, decltype(priority_function_)> queue(priority_function_, num_top_paths_);
    AlignedPath path(k_, sequence_begin, sequence_begin);

    suffix_seed(queue, sequence_begin, sequence_end);
    explored_seeds_counter_ += queue.get_num_items_pushed();
    // Graph exploration process.
    while (!queue.empty() && alternative_paths.size() < num_alternative_paths_) {
        path = std::move(queue.top());
        queue.pop();
//        std::cerr << "popping path (#matches, size, cigar, seq) : " << path.get_num_matches()
//                  << " " << path.size() + k_ - 1
//                  << " " << path.get_cigar()
//                  << " " << path.get_sequence() << std::endl;
        if (path.get_query_it() + k_ > sequence_end) {
            if (path.size() != 0 && path.size() < sw_threshold_for_stitched_path_) {
                smith_waterman_score(path);
                trim(path);
            }
            alternative_paths.push(std::move(path));
            continue;
        }
        if (early_discard_path(path))
            continue;
//        std::cerr << "calling exact map" << std::endl;

        if (!exact_map(path, sequence_end, dp_alignment, terminate_mapping)) {
            if (path.size() != 0 && path.size() < sw_threshold_for_stitched_path_) {
                smith_waterman_score(path);
                trim(path);
            }
            alternative_paths.push(std::move(path));
            break;
        }

        if (path.get_query_it() + k_ > sequence_end) {
            if (path.size() != 0 && path.size() < sw_threshold_for_stitched_path_) {
                smith_waterman_score(path);
                trim(path);
            }
            alternative_paths.push(std::move(path));
            continue;
        }

        if (path.get_similar())
            continue;

//        std::cerr << "calling inexact map" << std::endl;

        if (!inexact_map(path, queue, dp_alignment, terminate_mapping) || queue.empty()) {
//            std::cerr << "Aligning and trim after inexact map: " << path.get_sequence() << " " << path.get_cigar() << std::endl;
            if (path.size() != 0 && path.size() < sw_threshold_for_stitched_path_) {
                smith_waterman_score(path);
                trim(path);
            }
            alternative_paths.push(std::move(path));
            break;
        }
    }
    if (alternative_paths.size() == 0)
        return std::vector<AlignedPath>();
    std::vector<AlignedPath> alternative_paths_vec;
    while(!alternative_paths.empty()) {
        alternative_paths_vec.push_back(std::move(alternative_paths.top()));
        alternative_paths.pop();
    }
    explored_paths_counter_ += queue.get_num_items_pushed();
    return alternative_paths_vec;
}

void DBGAligner::suffix_seed(BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>, decltype(priority_function_)>& queue,
                             const std::string::const_iterator &sequence_begin,
                             const std::string::const_iterator &sequence_end) const {

    AlignedPath path(k_, sequence_begin, sequence_begin);

    // TODO: perform inexact seeding in case of DBGs other than BOSS by skipping ahead.
    auto seed = graph_->kmer_to_node(std::string(sequence_begin, sequence_begin + k_));
    if (seed != graph_->npos) {
        path.seed(seed, {}, graph_->get_node_sequence(seed), match_score_ * k_);
        queue.push(std::move(path));
        return;
    }
    // Inexact seeding.
    auto begin = sequence_begin;
    // A while loop helps to pick a starting position where there aren't too many suffix seeds found.
    // Too many suffix seeds are expensive to explore and may lead to discarding desirable seeds incorrectly.
    while (begin + k_ < sequence_end) {
        size_t k_for_seeding_l = k_ > 12 ? 12 : 1;
        size_t k_for_seeding_r = k_;
//        std::cerr << "inexact seeding " << std::endl;

        // Perform binary search to find the maximum k_for_seeding for which a seed can be found.
        while (k_for_seeding_l + 1 < k_for_seeding_r) {
            size_t k_for_seeding_m = (k_for_seeding_l + k_for_seeding_r) / 2;
            bool seed_found = false;
            graph_->suffix_seeding(begin + k_ - k_for_seeding_m, begin + k_,
                    [&] (node_index) { seed_found = true; },
                    [&] () { return seed_found; });
            if (seed_found)
                k_for_seeding_l = k_for_seeding_m;
            else
                k_for_seeding_r = k_for_seeding_m;
        }
        // Seed by suffix with the k range after the binary search.
        graph_->suffix_seeding(begin + k_ - k_for_seeding_l, begin + k_, [&] (node_index node) {
            auto alternative_path(path);
//            std::cerr << "seeding alternative: " << graph_->get_node_sequence(node) << std::endl;
            alternative_path.seed(node, {}, graph_->get_node_sequence(node));
            queue.push(std::move(alternative_path));
            });
        if (queue.size() > 0)
            return;
//        std::cerr << "Continue alignment with l and r: " << k_for_seeding_l << ", " << k_for_seeding_r << std::endl;
        begin += k_ - k_for_seeding_r + 1;
    }
    std::cout << "Warning! no seed found for path" << std::endl;
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
                           (path.size() != 0 && graph_->outdegree(path.back()) > 1)); },
            seed);
    return continue_alignment;
}

bool DBGAligner::inexact_map(AlignedPath& path,
                             BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>, decltype(priority_function_)> &queue,
                             std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                             const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate) {
    bool continue_alignment = true;
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char extension) {
        AlignedPath alternative_path(path);
        alternative_path.set_query_it(path.get_query_it());
//        std::cerr << "Alternative path: " << alternative_path.get_sequence() << std::endl;
        // Stop early in case of early termination condition by the caller function.
        if (!continue_alignment || terminate(node, alternative_path.get_query_it())) {
            continue_alignment = false;
            return;
        }
//        std::cerr << "Computing sw score" << std::endl;
        // Approximate Smith-Waterman score by assuming either match/mismatch for the next character.
        if (use_cssw_lib_ //|| (queue.size() > 0 && queue.back() < alternative_path)
            || path.size() > sw_threshold_for_stitched_path_) {
        alternative_path.extend(node, {}, extension,
            single_char_score(*(alternative_path.get_query_it() +
                k_ - 1), extension));
        } else {
            alternative_path.extend(node, {}, extension);
            smith_waterman_score(alternative_path);
        }
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
//                std::cerr << "similar path" << std::endl;
                path.set_similar();
            }
        }
        else {
            queue.push(std::move(alternative_path));
        }
        });
    return continue_alignment;
}

int8_t DBGAligner::single_char_score(char char_in_query, char char_in_graph) const {
    return score_matrix_[int8_t(char_in_query)][int8_t(char_in_graph)];
}

bool DBGAligner::smith_waterman_score(AlignedPath &path) {
    if (path.is_score_updated())
        return true;

    StripedSmithWaterman::Alignment alignment;
    // Calculate and set the alignment properties according to the smith_waterman_core function.
    if (!use_cssw_lib_) {
        SWDpCell optimum_alignment;
        smith_waterman_core(path, optimum_alignment);
        if (optimum_alignment.score >= 0) {
            alignment.cigar_string = optimum_alignment.cigar.to_string();
            alignment.sw_score = optimum_alignment.score;
            alignment.ref_begin = optimum_alignment.cigar.get_ref_begin();
            alignment.query_begin = optimum_alignment.cigar.get_query_begin();
            alignment.ref_end = path.get_sequence().size() - 1 - optimum_alignment.cigar.get_ref_end();
            alignment.query_end = path.get_query_sequence().size() - 1 - optimum_alignment.cigar.get_ref_end();
        } else {
            std::cerr << "Warning: Negative Smith Waterman score!" << std::endl;
            alignment.cigar_string = "";
            alignment.sw_score = 0u;
            alignment.ref_begin = optimum_alignment.cigar.get_ref_begin();
            alignment.query_begin = optimum_alignment.cigar.get_query_begin();
            alignment.ref_end = path.get_sequence().size() - 1;
            alignment.query_end = path.get_query_sequence().size() - 1;
        }
        path.update_alignment(alignment);
//        std::cerr << "path " << path.get_sequence() << ", cigar: " << path.get_cigar()
//                  << ", size " << path.get_sequence().size() << ", ref_end: " << alignment.ref_end
 //                 << ", query_end: " << alignment.query_end << std::endl;
        return true;
    }

    auto ref = path.get_sequence();
    std::string query(path.get_query_begin_it(), path.get_query_it() + k_ - 1);
    int32_t maskLen = strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    StripedSmithWaterman::Filter filter;
    if (!cssw_aligner_.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen)){
        std::cout << "Failure in SSW calculation. Cigar string might be missing.\n";
        return false;
    }

    path.update_alignment(alignment);
    return true;
}

void DBGAligner::smith_waterman_dp_step(size_t i, size_t j, const std::string& path_sequence, const std::string& query_sequence) {
    std::vector<SWDpCell> options(3);

    // Option for insertion operation
    const auto& parent_cell_for_in_op = sw_prev_row_->operator[](j);
    options[0].cigar = parent_cell_for_in_op.cigar;
    options[0].cigar.append(Cigar::Operator::INSERTION);
    if (parent_cell_for_in_op.cigar.size() == 0
        || parent_cell_for_in_op.cigar.back() != Cigar::Operator::INSERTION) // Gap opening
        options[0].score = parent_cell_for_in_op.score - gap_opening_penalty_;
    else // Gap extension
        options[0].score = parent_cell_for_in_op.score - gap_extension_penalty_;

    // Option for deletion operation
    const auto& parent_cell_for_del_op = sw_cur_row_->operator[](j-1);
    options[1].cigar = parent_cell_for_del_op.cigar;
    options[1].cigar.append(Cigar::Operator::DELETION);
    if (parent_cell_for_del_op.cigar.size() == 0
        || parent_cell_for_del_op.cigar.back() != Cigar::Operator::DELETION) // Gap opening
        options[1].score = parent_cell_for_del_op.score - gap_opening_penalty_;
    else // Gap extension
        options[1].score = parent_cell_for_del_op.score - gap_extension_penalty_;

    // Option for match/mismatch operation
   const  auto& parent_cell_for_sub_op = sw_prev_row_->operator[](j-1);
    auto match_mismatch_score = single_char_score(query_sequence[i-1], path_sequence[j-1]);
    options[2].score = parent_cell_for_sub_op.score + match_mismatch_score;
    options[2].cigar = parent_cell_for_sub_op.cigar;
    if (match_mismatch_score > 0) // Match
        options[2].cigar.append(Cigar::Operator::MATCH);
    else if (match_mismatch_score == -1) // Mismatch with transition mutation
        options[2].cigar.append(Cigar::Operator::MISMATCH_TRANSITION);
    else // Mismatch with tranversion mutation
        options[2].cigar.append(Cigar::Operator::MISMATCH_TRANSVERSION);

    auto max_option = std::max_element(std::begin(options), std::end(options));
    sw_cur_row_->operator[](j).score = max_option->score;
    sw_cur_row_->operator[](j).cigar = std::move(max_option->cigar);
//    std::cerr << "max: " << sw_cur_row_->operator[](j).cigar.to_string() << ", " << sw_cur_row_->operator[](j).score << std::endl;
}

void DBGAligner::smith_waterman_core(AlignedPath& path, SWDpCell& optimum_alignment) {
    auto path_sequence = path.get_sequence();
    std::string query_sequence(path.get_query_begin_it(), path.get_query_it() + k_ - 1);
    auto num_cols = path_sequence.size() + 1;
    auto num_rows =  query_sequence.size() + 1;
    // Keep only two rows of the table at a time to save space.
    // Keep the last column in path to continue
    // filling the table when path grows from that column.
    // Check if SW related vectors are allocated and otherwise, allocate them.
    if (sw_prev_row_->size() < num_cols)
        sw_prev_row_->resize(2 * num_cols);
    if (sw_cur_row_->size() < num_cols)
        sw_cur_row_->resize(2 * num_cols);
    if (sw_last_column_.size() < num_rows)
        sw_last_column_.resize(2 * num_rows);

    auto num_rows_stored = path.get_sw_num_rows_stored();
    auto num_cols_stored = path.get_sw_num_cols_stored();

    assert(sw_cur_row_->size() == sw_prev_row_->size());
    assert(num_rows_stored < sw_cur_row_->size());
    // Fill the first row with trivial values. Note the swap action in the beginning of the loops.
    for (size_t i = 0; i < num_cols; ++i)
        sw_cur_row_->operator[](i) = SWDpCell{};
//    std::cerr << "smith waterman score for path: " << path_sequence
//              << " sequence: " << query_sequence << std::endl;
    for (uint64_t i = 1; i < num_rows; ++i) {
        size_t j = 1;
        auto tmp = sw_prev_row_;
        sw_prev_row_ = sw_cur_row_;
        sw_cur_row_ = tmp;

        if (i < num_rows_stored) {
            j = num_cols_stored - 1;
            // prev_row is filled with correct values for the rows after the first row by construction.
            sw_cur_row_->operator[](j) = path.get_sw_last_column(i);
            ++ j;
        } else if (i == num_rows_stored) {
            std::copy(path.get_sw_last_row_it(), path.get_sw_last_row_it() + num_cols_stored, sw_prev_row_->begin());
//            for (size_t k = 0; k < num_cols_stored; ++k)
//                sw_prev_row_->operator[](k) = path.get_sw_last_row(k);
        }
        for (; j < num_cols; ++j) {
//            std::cerr << "Smith waterman step for i, j: " << i << ", " << j << std::endl;
            smith_waterman_dp_step(i, j, path_sequence, query_sequence);
        }
        sw_last_column_[i] = sw_cur_row_->operator[](path_sequence.size());
    }
    optimum_alignment = sw_cur_row_->operator[](path_sequence.size());
//    std::cerr << "optimum_alignment with score " << optimum_alignment.score << " and cigar: " << optimum_alignment.cigar.to_string() << std::endl;
    auto updated_score_delta = optimum_alignment.cigar.clip(gap_opening_penalty_, gap_extension_penalty_, -1 * mm_transition_, -1 * mm_transversion_);
    optimum_alignment.score += updated_score_delta;
    path.store_sw_intermediate_info(sw_cur_row_->begin(), path_sequence.size() + 1, std::begin(sw_last_column_), query_sequence.size() + 1);
//    std::cerr << "optimum alinment after clipping: " << optimum_alignment.score << " " << optimum_alignment.cigar.to_string() << std::endl;
}

// Note: Assume Cigar is updated.
void DBGAligner::trim(AlignedPath &path) const {
    auto cigar = path.get_cigar();
    // Trim path if reference end in cigar is not according to reference end in path.
    uint64_t path_trim_length = path.get_sequence().size() - path.get_alignment().ref_end - 1;
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

void DBGAligner::fill_score_matrix(size_t num_elements) {
    std::fill_n(&(score_matrix_[0][0]), num_elements, std::min(mm_transversion_, mm_transition_));
    // Match scores
    score_matrix_['a']['a'] = score_matrix_['A']['a'] = score_matrix_['a']['A'] = score_matrix_['A']['A'] = match_score_;
    score_matrix_['c']['c'] = score_matrix_['C']['c'] = score_matrix_['c']['C'] = score_matrix_['C']['C'] = match_score_;
    score_matrix_['g']['g'] = score_matrix_['G']['g'] = score_matrix_['g']['G'] = score_matrix_['G']['G'] = match_score_;
    score_matrix_['t']['t'] = score_matrix_['T']['t'] = score_matrix_['t']['T'] = score_matrix_['T']['T'] = match_score_;
    // Mismatch transition scores
    score_matrix_['a']['g'] = score_matrix_['A']['g'] = score_matrix_['a']['G'] = score_matrix_['A']['G'] = mm_transition_;
    score_matrix_['g']['a'] = score_matrix_['G']['a'] = score_matrix_['g']['A'] = score_matrix_['G']['A'] = mm_transition_;
    score_matrix_['c']['t'] = score_matrix_['C']['t'] = score_matrix_['c']['T'] = score_matrix_['C']['T'] = mm_transition_;
    score_matrix_['t']['c'] = score_matrix_['T']['c'] = score_matrix_['t']['C'] = score_matrix_['T']['C'] = mm_transition_;
    // Mismatch transversion scores
    score_matrix_['a']['c'] = score_matrix_['A']['c'] = score_matrix_['a']['C'] = score_matrix_['A']['C'] = mm_transversion_;
    score_matrix_['a']['t'] = score_matrix_['A']['t'] = score_matrix_['a']['T'] = score_matrix_['A']['T'] = mm_transversion_;

    score_matrix_['g']['c'] = score_matrix_['G']['c'] = score_matrix_['g']['C'] = score_matrix_['G']['C'] = mm_transversion_;
    score_matrix_['g']['t'] = score_matrix_['G']['t'] = score_matrix_['g']['T'] = score_matrix_['G']['T'] = mm_transversion_;

    score_matrix_['c']['a'] = score_matrix_['C']['a'] = score_matrix_['c']['A'] = score_matrix_['C']['A'] = mm_transversion_;
    score_matrix_['c']['g'] = score_matrix_['C']['g'] = score_matrix_['c']['G'] = score_matrix_['C']['G'] = mm_transversion_;

    score_matrix_['t']['a'] = score_matrix_['T']['a'] = score_matrix_['t']['A'] = score_matrix_['T']['A'] = mm_transversion_;
    score_matrix_['t']['g'] = score_matrix_['T']['g'] = score_matrix_['t']['G'] = score_matrix_['T']['G'] = mm_transversion_;
}
