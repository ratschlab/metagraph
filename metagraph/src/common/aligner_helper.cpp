#include "aligner_helper.hpp"


std::string Cigar::to_string() const {
    std::string str = "";
    for (auto cigar_part = std::begin(cigar_);
         cigar_part != std::end(cigar_); ++ cigar_part) {
        str += std::to_string(cigar_part->second);
        str += opt_to_str(cigar_part->first);
    }
    return str;
}

void Cigar::append(const Operator& op) {
    if (cigar_.size() == 0 || cigar_.back().first != op)
        cigar_.push_back(std::make_pair(op, 1));
    else
        cigar_.back().second += 1;
}

uint64_t Cigar::clip(uint64_t gap_opening_penalty, uint64_t gap_extension_penalty,
          uint64_t mismatch_penalty_transition, uint64_t mismatch_penalty_transversion) {
    assert(cigar_.size() == 0);
    query_begin_ = 0;
    ref_begin_ = 0;
    query_end_ = 0;
    ref_end_ = 0;
    int64_t updated_score_delta = 0;
    // Detect mismatches, insertions and deletions from the beginning
    // and adjust cigar, query and ref begin values.
    for (auto cigar_op = std::begin(cigar_);
         cigar_op != std::end(cigar_)
            && cigar_op->first != Operator::MATCH
            && cigar_op->first != Operator::CLIPPED; ++ cigar_op) {
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
        }
    }
    // Detect mismatches, insertions and deletions from the end
    // and adjust cigar, query and ref end values.
    for (auto cigar_op = std::end(cigar_) - 1;
             cigar_op != std::begin(cigar_)
                && cigar_op->first != Operator::MATCH
                && cigar_op->first != Operator::CLIPPED; -- cigar_op) {
        if (cigar_op->first == Operator::INSERTION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            query_end_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSITION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transition * cigar_op->second;
            query_end_ += cigar_op->second;
            ref_end_ += cigar_op->second;
        } else if (cigar_op->first == Operator::MISMATCH_TRANSVERSION) {
            cigar_op->first = CLIPPED;
            updated_score_delta += mismatch_penalty_transversion * cigar_op->second;
            query_end_ += cigar_op->second;
            ref_end_ += cigar_op->second;
        } else if (cigar_op->first == Operator::DELETION) {
            updated_score_delta += gap_opening_penalty + gap_extension_penalty * (cigar_op->second - 1);
            ref_end_ += cigar_op->second;
        }
    }
    // In case of adjacent clipped operations in cigar, merge them.
    for (auto cigar_op = std::begin(cigar_);
             cigar_op < std::end(cigar_) - 1; ++ cigar_op) {
        if (cigar_op->first == Operator::CLIPPED
                && (cigar_op + 1)->first == Operator::CLIPPED) {
            cigar_op->second += (cigar_op + 1)->second;
            cigar_op = cigar_.erase(cigar_op + 1) - 2;
        }
    }
//    std::cerr << "Updated score delta: " << updated_score_delta
//              << " for cigar: " << to_string() << std::endl;
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
