#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cstring>

#include "annotated_dbg.hpp"
#include "bounded_priority_queue.hpp"
#include "path.hpp"
#include "sequence_graph.hpp"
#include "ssw_cpp.h"
#include "aligner_helper.hpp"

struct DBGAlignerConfig {
    size_t num_top_paths = 10;
    size_t num_alternative_paths = 1;
    uint8_t path_comparison_code = 0;
    bool verbose = false;
    bool discard_similar_paths = true;
    bool use_cssw_lib = false;
    bool disable_cssw_speedup = false;
    size_t sw_threshold_for_stitched_path = 200;
    float gap_opening_penalty = 3;
    float gap_extension_penalty = 1;
    float  mm_transition = -1;
    float mm_transversion = -2;
    float match_score = 2;
};

class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, AnnotatedDBG::Annotator::VLabels> AlignedPath;

    DBGAligner(std::shared_ptr<DeBruijnGraph> graph, DBGAlignerConfig dbg_aligner_config = default_config);

    DBGAligner() = delete;
    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Detecting similar paths.
    struct DPAlignmentKey {
        // The node in the underlying graph.
        node_index node;
        std::string::const_iterator query_begin_it;
        std::string::const_iterator query_it;

        bool operator< (const DPAlignmentKey &other) const {
            return (query_begin_it == other.query_begin_it) ?
            ((query_it == other.query_it) ? (node < other.node) :
            (query_it < other.query_it)) : query_begin_it < other.query_begin_it;
        }
    };
    struct DPAlignmentValue {
        float score;
    };

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    std::vector<AlignedPath> align_by_graph_exploration(const std::string::const_iterator &sequence_begin,
        const std::string::const_iterator &sequence_end,
        const std::function<bool (AlignedPath&)>& early_discard_path = [](AlignedPath&) { return false; },
        const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate_mapping
            = [](node_index, const std::string::const_iterator&) { return false; });

    // Align a sequence to the underlying graph using map_to_nodes.
    bool map_to_nodes(const std::string &sequence, std::vector<AlignedPath>& alternative_paths,
                      int64_t alternative_alignment_score = 0);

    // Align to nodes for both sequence and reverse complement sequence and return the higher scoring one.
    AlignedPath map_to_nodes_forward_reverse_complement(const std::string &sequence);

    float get_match_score() const { return match_score_; }
    long long get_num_merge_paths() const { return merged_paths_counter_; }
    long long get_num_explored_paths() const { return explored_paths_counter_; }
    long long get_num_explored_seeds() const { return explored_seeds_counter_; }

    // Priority function for path extension phase.
    static std::function<bool(const AlignedPath&, const AlignedPath&)> priority_function_;

  private:
    std::shared_ptr<DeBruijnGraph> graph_;

   constexpr static DBGAlignerConfig default_config = {
    .num_top_paths = 10,
    .num_alternative_paths = 1,
    .path_comparison_code = 0,
    .verbose = false,
    .discard_similar_paths = true,
    .use_cssw_lib = false,
    .disable_cssw_speedup = false,
    .sw_threshold_for_stitched_path = 200,
    .gap_opening_penalty = 3,
    .gap_extension_penalty = 1,
    .mm_transition = -1,
    .mm_transversion = -2,
    .match_score = 2};

    // Substitution score for each pair of nucleotides.
    //std::map<char, std::map<char, int8_t>> sub_score_;
    int8_t score_matrix_[128][128];
    int8_t mm_transition_;
    int8_t mm_transversion_;
    // Maximum number of paths to explore at the same time.
    size_t num_top_paths_;
    size_t num_alternative_paths_;

    uint8_t path_comparison_code_;
    bool verbose_;
    bool discard_similar_paths_;
    bool use_cssw_lib_;
    bool disable_cssw_speedup_;
    size_t sw_threshold_for_stitched_path_;

    int8_t match_score_;
    float gap_opening_penalty_;
    float gap_extension_penalty_;

    StripedSmithWaterman::Aligner cssw_aligner_;

    long long merged_paths_counter_;
    long long explored_paths_counter_;
    long long explored_seeds_counter_;
    size_t k_;

    // Smith-Waterman method related fields.
    // Keep only two rows of the table at a time to save space.
    std::vector<SWDpCell>* sw_prev_row_;
    std::vector<SWDpCell>* sw_cur_row_;
    // Keep the last column in path to continue
    // filling the table when path grows from that column.
    std::vector<SWDpCell> sw_last_column_;

    // Seed the paths. In case of no exact seed, run suffix seeding. To find the size of suffix to seed based on,
    // Perform binary search and report a value of suffix seeding lead to more than 0, but limited number of seeds.
    // If too many seeds are extractable, it is possible that some seeds are discarded which leads to possible poor alignment.
    void suffix_seed(BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>, decltype(priority_function_)> &queue,
                     const std::string::const_iterator &sequence_begin,
                     const std::string::const_iterator &sequence_end);

    // Align part of a sequence to the graph in the case of no exact map
    // based on internal strategy. Calls callback for every possible alternative path.
    bool inexact_map(AlignedPath &path, BoundedPriorityQueue<AlignedPath, std::vector<AlignedPath>, decltype(priority_function_)> &queue,
                     std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                     const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate
                         = [](node_index, const std::string::const_iterator&) { return false; });

    // Align the path to the graph based on the query until either exact alignment is not
    // possible or there exists a branching point in the graph.
    bool exact_map(AlignedPath &path, const std::string::const_iterator &sequence_end,
                   std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                   const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate
                        = [](node_index, const std::string::const_iterator&) { return false; });

    // Return the score of substitution. If not in sub_score_ return a fixed maximized score value.
    int8_t single_char_score(char char_in_query, char char_in_graph) const;

    // Compute the edit distance between the query sequence and the aligned path
    // according to score parameters in this class.
    void smith_waterman_core(AlignedPath& path, SWDpCell& global_sw, bool clip);

    // Perform one dynamic programming step consisting of making a choice among insertion, deletion, map and mismatch
    // for a single cell in the dynamic programming table.
    void smith_waterman_dp_step(size_t i, size_t j, const std::string& path_sequence, const std::string& query_sequence);

    // Compute Smith-Waterman score based on either CSSW library or our implementation
    // and set cigar and score in path accordingly.
    // Clip can only be true for the case of using CSSW.
    bool smith_waterman_score(AlignedPath &path, bool clip=true);

    // Trim the path to remove unmapped regions from the tail using CSSW lib clipping.
    void trim(AlignedPath &path) const;

    // Fill the score matrix with total number of num_elements cells in the matrix.
    void fill_score_matrix(size_t num_elements);
};

#endif // __DBG_ALIGNER_HPP__
