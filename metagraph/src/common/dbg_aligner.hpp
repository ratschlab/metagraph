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


class Cigar {
  public:

    enum Operator {
        MATCH,
        MISMATCH_TRANSITION,
        MISMATCH_TRANSVERSION,
        INSERTION,
        DELETION,
        CLIPPED };

    size_t size() { return cigar_.size(); }
    Operator back() { return cigar_.back().first; }
    uint64_t get_ref_begin() const { return ref_begin_; }
    uint64_t get_ref_end() const { return ref_end_; }
    uint64_t get_query_begin() const { return query_begin_; }
    uint64_t get_query_end() const { return query_end_; }

    std::string to_string() const {
        std::string str = "";
        for (auto cigar_part = std::begin(cigar_);
             cigar_part != std::end(cigar_); ++ cigar_part) {
            str += std::to_string(cigar_part->second);
            str += opt_to_str(cigar_part->first);
        }
        return str;
    }
    void append(const Operator& op) {
        if (cigar_.size() == 0 || cigar_.back().first != op)
            cigar_.push_back(std::make_pair(op, 1));
        else
            cigar_.back().second += 1;
    }
    // Clip mismatches, insertions and deletions from the begin and end.
    // Adjust query and ref, begin and end positions accordingly.
    // returns a positive score saved by not calculating clipped penalties.
    uint64_t clip(uint64_t gap_opening_penalty, uint64_t gap_extension_penalty,
                  uint64_t mismatch_penalty_transition, uint64_t mismatch_penalty_transversion);

  private:
    uint64_t ref_begin_, ref_end_, query_begin_, query_end_;
    std::vector<std::pair<Operator, uint64_t>> cigar_;

    static std::string opt_to_str(const Operator& op) {
        switch (op) {
            case MATCH:
                return "=";
            case MISMATCH_TRANSITION:
                return "X";
            case MISMATCH_TRANSVERSION:
                return "X";
            case INSERTION:
                return "I";
            case DELETION:
                return "D";
            case CLIPPED:
                return "S";
            default:
                assert(false);
                return "";
        }
    }
};

// Smith-Waterman table cell type
struct SWDpCell {
    int64_t score;
    Cigar cigar;
    bool operator< (const SWDpCell& other) const {
        return score < other.score;
    }
};


class DBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Path<node_index, AnnotatedDBG::Annotator::VLabels> AlignedPath;


    DBGAligner(std::shared_ptr<DeBruijnGraph> graph,
               size_t num_top_paths = 10,
               size_t num_alternative_paths = 1,
               uint8_t path_comparison_code = 0,
               bool verbose = false,
               bool discard_similar_paths = false,
               bool use_cssw_lib = false,
               float insertion_penalty = 3,
               float deletion_penalty = 3,
               float gap_opening_penalty = 3,
               float gap_extension_penalty = 1);

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
    int8_t mm_transition_;
    int8_t mm_transversion_;
    // Maximum number of paths to explore at the same time.
    size_t num_top_paths_;
    size_t num_alternative_paths_;
    uint8_t path_comparison_code_;
    bool verbose_;
    bool discard_similar_paths_;
    bool use_cssw_lib_;
    int8_t match_score_;
    float insertion_penalty_;
    float deletion_penalty_;
    float gap_opening_penalty_;
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
    int8_t single_char_score(char char_in_query, char char_in_graph) const;

    // Compute the edit distance between the query sequence and the aligned path
    // according to score parameters in this class.
    void whole_path_score(AlignedPath& path, SWDpCell& global_sw) const;

    // Compute Smith-Waterman score based on either CSSW library or our implementation
    // and set cigar and score in path accordingly.
    bool smith_waterman_score(AlignedPath &path) const;

    // Trim the path to remove unmapped regions from the tail using CSSW lib clipping.
    void trim(AlignedPath &path) const;
};

#endif // __DBG_ALIGNER_HPP__
