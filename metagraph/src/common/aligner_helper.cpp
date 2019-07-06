#include "aligner_helper.hpp"

#include <iostream>


std::string Cigar::to_string() const {
    std::string str = "";
    for (auto cigar_part = std::begin(cigar_);
         cigar_part != std::end(cigar_); ++ cigar_part) {
        str += std::to_string(cigar_part->second);
        str += opt_to_str(cigar_part->first);
    }
    return str;
}

void Cigar::append(const Operator& op, uint64_t num) {
    if (cigar_.size() == 0 || cigar_.back().first != op)
        cigar_.push_back(std::make_pair(op, num));
    else
        cigar_.back().second += num;
}

uint64_t Cigar::clip(uint64_t match_score, uint64_t gap_opening_penalty, uint64_t gap_extension_penalty,
          uint64_t mismatch_penalty_transition, uint64_t mismatch_penalty_transversion) {
    if (cigar_.size() == 0)
        return 0;
    query_begin_ = 0;
    ref_begin_ = 0;
    query_end_offset_ = 0;
    ref_end_offset_ = 0;
    auto clipped_cigar(cigar_);
    int64_t updated_score_delta = 0;
    Cigar max_scoring_cigar(*this);
    int64_t max_score_delta = 0;
//    std::cerr << "Clipping from front: " << std::endl;
    // Detect mismatches, insertions and deletions from the beginning
    // and adjust cigar, query and ref begin values.
    for (auto cigar_op = std::begin(clipped_cigar);
         cigar_op != std::end(clipped_cigar); ++ cigar_op) {
        if (cigar_op->first == Operator::INSERTION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            query_begin_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSITION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transition * cigar_op->second;
            query_begin_ += cigar_op->second;
            ref_begin_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSVERSION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transversion * cigar_op->second;
            query_begin_ += cigar_op->second;
            ref_begin_ += cigar_op->second;
        } else if (cigar_op->first == Operator::DELETION) {
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            ref_begin_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MATCH) {
            cigar_op->first = CLIPPED;
            updated_score_delta -= match_score * cigar_op->second;
            query_begin_ += cigar_op->second;
            ref_begin_ += cigar_op->second;
        }
        // Update cigar with maximum score after clipping.
        if (updated_score_delta >= max_score_delta) {
            max_score_delta = updated_score_delta;
            max_scoring_cigar.cigar_ = clipped_cigar;
            max_scoring_cigar.query_begin_ = query_begin_;
            max_scoring_cigar.ref_begin_ = ref_begin_;
        }
    }
//    std::cerr << "Max scoring cigar " << max_scoring_cigar.to_string() << std::endl;
    clipped_cigar = max_scoring_cigar.cigar_;
    updated_score_delta = max_score_delta;
//    std::cerr << "Clipping from back: " << std::endl;
    // Detect mismatches, insertions and deletions from the end
    // and adjust cigar, query and ref end values.
    for (auto cigar_op = std::end(clipped_cigar) - 1;
             cigar_op != std::begin(clipped_cigar); -- cigar_op) {
        if (cigar_op->first == Operator::INSERTION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            query_end_offset_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSITION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transition * cigar_op->second;
            query_end_offset_ += cigar_op->second;
            ref_end_offset_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSVERSION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transversion * cigar_op->second;
            query_end_offset_ += cigar_op->second;
            ref_end_offset_ += cigar_op->second;
        } else if (cigar_op->first == Operator::DELETION) {
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            ref_end_offset_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MATCH) {
            cigar_op->first = CLIPPED;
            updated_score_delta -= match_score * cigar_op->second;
            query_end_offset_ += cigar_op->second;
            ref_end_offset_ += cigar_op->second;
        } 
        // Update cigar with maximum score after clipping.
        if (updated_score_delta >= max_score_delta) {
            max_score_delta = updated_score_delta;
            max_scoring_cigar.cigar_ = clipped_cigar;
            max_scoring_cigar.query_end_offset_ = query_end_offset_;
            max_scoring_cigar.ref_end_offset_ = ref_end_offset_;
        }
    }
//    std::cerr << "Max scoring cigar " << max_scoring_cigar.to_string() << std::endl;

    ref_begin_ = max_scoring_cigar.ref_begin_;
    ref_end_offset_ = max_scoring_cigar.ref_end_offset_;
    query_begin_ = max_scoring_cigar.query_begin_;
    query_end_offset_ = max_scoring_cigar.query_end_offset_;
    clipped_cigar = max_scoring_cigar.cigar_;
//    std::cerr << "query_end_offset_: " << query_end_offset_ << std::endl;
    updated_score_delta = max_score_delta;

//    std::cerr << "Updated score delta: " << updated_score_delta
//              << " cigar before: " << to_string() << std::endl;
    cigar_ = std::vector<std::pair<Operator, uint64_t>>();
    cigar_.reserve(clipped_cigar.size());
    // In case of adjacent clipped operations in cigar, merge them.
    for (auto cigar_op = std::begin(clipped_cigar);
             cigar_op < std::end(clipped_cigar); ++ cigar_op) {
        this->append(cigar_op->first, cigar_op->second);
    }
//    std::cerr << " cigar after: " << to_string() << std::endl;
    return updated_score_delta;
}

std::string Cigar::opt_to_str(const Operator& op) {
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
