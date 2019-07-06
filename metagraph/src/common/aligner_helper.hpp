#ifndef __ALIGNER_HELPER_HPP__
#define __ALIGNER_HELPER_HPP__

#include <vector>
#include <string>
#include <cassert>

class Cigar {

  public:
    enum Operator {
        MATCH,
        MISMATCH_TRANSITION,
        MISMATCH_TRANSVERSION,
        INSERTION,
        DELETION,
        CLIPPED };

    Cigar() {
        cigar_.reserve(10);
    }

    size_t size() const { return cigar_.size(); }
    Operator back() const { return cigar_.back().first; }
    uint64_t get_ref_begin() const { return ref_begin_; }
    uint64_t get_ref_end_offset() const { return ref_end_offset_; }
    uint64_t get_query_begin() const { return query_begin_; }
    uint64_t get_query_end_offset() const { return query_end_offset_; }

    std::string to_string() const;

    void append(const Operator& op, uint64_t num = 1);

    // Clip mismatches, insertions and deletions from the begin and end.
    // Adjust query and ref, begin and end positions accordingly.
    // returns a positive score saved by not calculating clipped penalties.
    uint64_t clip(uint64_t match_score, uint64_t gap_opening_penalty, uint64_t gap_extension_penalty,
              uint64_t mismatch_penalty_transition, uint64_t mismatch_penalty_transversion);

  private:
    uint64_t ref_begin_, ref_end_offset_, query_begin_, query_end_offset_;
    std::vector<std::pair<Operator, uint64_t>> cigar_;

    static std::string opt_to_str(const Operator& op);
};

// Smith-Waterman table cell type
struct SWDpCell {
    int64_t score;
    Cigar cigar;
    bool operator< (const SWDpCell& other) const {
        return score < other.score;
    }
};

#endif  // __ALIGNER_HELPER_HPP__
