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

class DBGAligner : public AnnotatedDBG {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, Annotator::VLabels> AlignedPath;

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

    DBGAligner(DeBruijnGraph *dbg,
               Annotator *annotation,
               size_t num_top_paths = 10,
               bool verbose = false,
               float sw_threshold = 0.8,
               float re_seeding_threshold = 0.6,
               float insertion_penalty = 3,
               float deletion_penalty = 3,
               size_t num_threads = 0);

    DBGAligner(const DBGAligner&) = default;
    DBGAligner(DBGAligner&&) = default;
    DBGAligner& operator= (const DBGAligner&) = default;
    DBGAligner& operator= (DBGAligner&&) = default;

    // Align a sequence to the underlying graph based on the strategy defined in the graph.
    std::vector<AlignedPath> align(const std::string &sequence);

    float get_match_score() const { return match_score_; }
    long long get_num_merge_paths() const { return merged_paths_counter_; }

  private:
    // Substitution score for each pair of nucleotides.
    std::map<char, std::map<char, int8_t>> sub_score_;
    // Maximum number of paths to explore at the same time.
    size_t num_top_paths_;
    bool verbose_;
    float sw_threshold_;
    float re_seeding_threshold_;
    int8_t match_score_;
    float insertion_penalty_;
    float deletion_penalty_;
    StripedSmithWaterman::Aligner cssw_aligner_;
    long long merged_paths_counter_;
    size_t k_;

    // Align part of a sequence to the graph in the case of no exact map
    // based on internal strategy. Calls callback for every possible alternative path.
    void inexact_map(AlignedPath &path, BoundedPriorityQueue<AlignedPath> &queue,
                     std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment);

    // Align the path to the graph based on the query until either exact alignment is not
    // possible or there exists a branching point in the graph.
    void exact_map(AlignedPath &path, const std::string &sequence,
                   std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment);

    // Return the score of substitution. If not in sub_score_ return a fixed maximized score value.
    float single_char_score(char char_in_query, char char_in_graph) const;

    // Compute the edit distance between the query sequence and the aligned path
    // according to score parameters in this class.
    float whole_path_score(const AlignedPath &path) const;

    // Compute the distance between the query sequence and the aligned path sequence
    // according to the CSSW library.
    float ssw_score(const AlignedPath &path) const;

    // Compute Smith-Waterman score based on CSSW library.
    bool cssw_align(const AlignedPath &path, StripedSmithWaterman::Alignment& alignment) const;

    // Trim the path to remove unmapped regions from the tail using CSSW lib clipping.
    void trim(AlignedPath &path) const;
};

#endif // __DBG_ALIGNER_HPP__
