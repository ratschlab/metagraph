#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <vector>
#include <regex>
#include <cstdlib>
#include <iostream>

#include "ssw_cpp.h"
#include "aligner_helper.hpp"

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path(uint64_t k, std::string::const_iterator query_begin_it,
         std::string::const_iterator query_it,
         uint8_t prioritizing_function_code) :
            k_(k),
            query_begin_it_(query_begin_it),
            query_it_(query_it) {
        score_ = 0.0;
        num_matches_ = 0;
        is_similar_ = false;
        is_score_updated_ = false;

        assert(prioritizing_function_code >= 0 && prioritizing_function_code <= 2);
        switch (prioritizing_function_code) {
            case 0:
                prioritizing_function = total_score_less_than;
                break;
            case 1:
                prioritizing_function = normalized_score_less_than;
                break;
            case 2:
                prioritizing_function = num_matches_less_than;
                break;
            default:
                std::cout << "Warning: Prioritizing function is not set properly in path construction."
                          << std::endl << "Setting it to default (total score comparison)." << std::endl;
                prioritizing_function = total_score_less_than;
            }
        }

    Path(const Path&) = default;
    Path(Path&&) = default;
    Path& operator= (const Path&) = default;
    Path& operator= (Path&&) = default;

    void seed(NodeType node, const VLabels &labels, std::string sequence, float score=0) {
        push_back(node, labels, score);
        num_matches_ += sequence.size();
        sequence_ += sequence;
    }
    void extend(NodeType node, const VLabels &labels, char extention, float score=0) {
        push_back(node, labels, score);
        //In case of a match, increment the counter.
        if (score > 0)
            num_matches_ += 1;
        is_score_updated_ = false;
        sequence_ += extention;
    }
    // Append a path to the current object. Fills in 0 in case of spaced nodes.
    // Avoid adding duplicate nodes and chars in case of positive overlap.
    void append_path(const Path& path, int64_t overlap_length, float score) {
        is_score_updated_ = false;
        assert(path.nodes_.size() > overlap_length);
        assert(overlap_length >= 0);
        nodes_.insert(std::end(nodes_), std::begin(path.nodes_) + overlap_length, std::end(path.nodes_));
        sequence_ += path.sequence_.substr(k_ - 1 + overlap_length);
        label_set_.insert(std::end(label_set_), std::begin(path.label_set_),
                          std::end(path.label_set_));
        query_it_ = path.query_it_;
        score_ = score;
    }
    // Note: Assume score is updated by user accordingly.
    void trim(uint64_t query_trim_length, uint64_t path_trim_length) {
        auto last_op_ind = cigar_.find_last_of("MX=ID");
        if (last_op_ind != std::string::npos)
            cigar_ = cigar_.substr(0, last_op_ind + 1);

        query_it_ -= query_trim_length;
        sequence_.resize(sequence_.size() - path_trim_length);
        nodes_.resize(nodes_.size() - path_trim_length);
        if (sw_last_row_.size() >= path_trim_length
                && sw_last_column_.size() >= query_trim_length) {
            sw_last_row_.resize(sw_last_row_.size() - path_trim_length);
            sw_last_column_.resize(sw_last_column_.size() - query_trim_length);
        }
    }

    void set_query_begin_it(std::string::const_iterator query_begin_it) {
        query_begin_it_ = query_begin_it; }
    void set_query_it(std::string::const_iterator query_it) {
        query_it_ = query_it; }
    void update_alignment(const StripedSmithWaterman::Alignment& alignment) {
        alignment_ = alignment;
        score_ = alignment.sw_score;
        is_score_updated_ = true;
        cigar_ = alignment.cigar_string;
        // Extracting the number of matches.
        num_matches_ = 0;
        std::regex match_regex("[0-9]+[M=]");
        std::smatch match;
            std::string cigar_cpy = cigar_;
            while (std::regex_search(cigar_cpy, match, match_regex)) {
                std::string position_str = match.str();
                position_str.pop_back();
                num_matches_ += atoi(position_str.c_str());
                cigar_cpy = match.suffix();
            }
    }

    void store_sw_intermediate_info(std::vector<SWDpCell>&& sw_last_row, std::vector<SWDpCell>&& sw_last_column) {
        sw_last_row_ = std::move(sw_last_row);
        sw_last_column_ = std::move(sw_last_column);
    }

    std::vector<SWDpCell> get_sw_last_column() const { return sw_last_column_; }
    std::vector<SWDpCell> get_sw_last_row() const { return sw_last_row_; }
    void set_similar() { is_similar_ = true; }
    NodeType back() const { return nodes_.back(); }
    NodeType front() const { return nodes_.front(); }
    NodeType last_parent() const { return nodes_.at(nodes_.size() - 1); }
    size_t size() const { return nodes_.size(); }
    size_t query_size() const { return query_it_ - query_begin_it_ + k_ - 1; }
    size_t num_kmers_in_query() const { return query_it_ - query_begin_it_; }

    float get_total_score() const { return score_; }
    VLabels get_labels() const { return label_set_; }
    std::vector<NodeType> get_nodes() const { return nodes_; }

    std::string::const_iterator get_query_it() const { return query_it_; }
    std::string::const_iterator get_query_begin_it() const { return query_begin_it_; }
    std::string get_sequence() const { return sequence_; }
    std::string get_query_sequence() const { return std::string(query_begin_it_, query_it_ + k_ - 1); }
    std::string get_cigar() const { return cigar_; }
    uint64_t get_num_matches() const { return num_matches_; }
    bool get_similar() const { return is_similar_; }
    bool is_score_updated() const { return is_score_updated_; }
    StripedSmithWaterman::Alignment get_alignment() { return alignment_; }

    bool operator< (const Path &other) const {
        return prioritizing_function(*this, other);
    }

  private:
    float score_;
    uint64_t k_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
    std::string::const_iterator query_begin_it_;
    std::string::const_iterator query_it_;
    std::string sequence_;
    std::string cigar_;
    bool is_similar_;
    uint64_t num_matches_;
    bool is_score_updated_;
    StripedSmithWaterman::Alignment alignment_;
    // Store intermediate values in Smith-Waterman table to save computations
    // when Smith-Waterman is called multiple times.
    std::vector<SWDpCell> sw_last_row_;
    std::vector<SWDpCell> sw_last_column_;
    std::function<bool(const Path&, const Path&)> prioritizing_function;

    void push_back(NodeType node, const VLabels &labels, float score=0) {
        nodes_.push_back(node);
        score_ += score;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
        ++ query_it_;
    }

    // Functions to Prioritize paths.
    static bool total_score_less_than (const Path& first, const Path& second) {
        return (first.score_ < second.score_); }
    static bool normalized_score_less_than (const Path& first, const Path& second) {
        return (first.score_/float(first.size()) < second.score_/float(second.size())); }
    static bool num_matches_less_than (const Path& first, const Path& second) {
        return (first.num_matches_ < second.num_matches_); }

};

#endif  // __PATH_HPP__
