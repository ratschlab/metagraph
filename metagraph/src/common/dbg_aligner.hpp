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

class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, AnnotatedDBG::Annotator::VLabels> AlignedPath;

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

    DBGAligner(std::shared_ptr<DeBruijnGraph> graph,
               size_t num_top_paths = 10,
               size_t num_alternative_paths = 1,
               uint8_t path_comparison_code = 0,
               bool verbose = false,
               bool discard_similar_paths = false,
               float insertion_penalty = 3,
               float deletion_penalty = 3,
               float gap_openning_penalty = 3,
               float gap_extension_penalty = 1);

    DBGAligner() = delete;
    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    std::vector<AlignedPath> align(const std::string::const_iterator &sequence_begin,
        const std::string::const_iterator &sequence_end,
        const std::function<bool(node_index, const std::string::const_iterator& query_it)>& terminate
            = [](node_index, const std::string::const_iterator&) { return false; });

    // Align a sequence to the underlying graph using map_to_nodes.
    AlignedPath map_to_nodes(const std::string &sequence);

    // Align to nodes for both sequence and reverse complement sequence and return the higher scoring one.
    AlignedPath map_to_nodes_forward_reverse_complement(const std::string &sequence);

    float get_match_score() const { return match_score_; }
    long long get_num_merge_paths() const { return merged_paths_counter_; }

  private:
    std::shared_ptr<DeBruijnGraph> graph_;
    // Substitution score for each pair of nucleotides.
    std::map<char, std::map<char, int8_t>> sub_score_;
    // Maximum number of paths to explore at the same time.
    size_t num_top_paths_;
    size_t num_alternative_paths_;
    uint8_t path_comparison_code_;
    bool verbose_;
    bool discard_similar_paths_;
    int8_t match_score_;
    float insertion_penalty_;
    float deletion_penalty_;
    float gap_openning_penalty_;
    float gap_extension_penalty_;
    StripedSmithWaterman::Aligner cssw_aligner_;
    long long merged_paths_counter_;
    size_t k_;

    // Align part of a sequence to the graph in the case of no exact map
    // based on internal strategy. Calls callback for every possible alternative path.
    bool inexact_map(AlignedPath &path, BoundedPriorityQueue<AlignedPath> &queue,
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
    float single_char_score(char char_in_query, char char_in_graph) const;

    // Compute the edit distance between the query sequence and the aligned path
    // according to score parameters in this class.
//    void whole_path_score(const AlignedPath &path, GlobalSW& global_sw) const;

    // Compute Smith-Waterman score based on CSSW library and set cigar and score in path accordingly.
    bool cssw_align(AlignedPath &path, StripedSmithWaterman::Alignment& alignment) const;

    // Trim the path to remove unmapped regions from the tail using CSSW lib clipping.
    void trim(AlignedPath &path, const StripedSmithWaterman::Alignment& alignment) const;
};

#endif // __DBG_ALIGNER_HPP__
