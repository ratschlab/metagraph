#include "aln_cigar.hpp"

#include "graph/representation/succinct/boss.hpp"

namespace mtg::graph::align_redone {


Cigar::Cigar(std::string_view cigar_str) {
    std::string op_count;
    for (char c : cigar_str) {
        switch (c) {
            case '=':
                cigar_.emplace_back(Cigar::MATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'X':
                cigar_.emplace_back(Cigar::MISMATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'I':
                cigar_.emplace_back(Cigar::INSERTION, std::stol(op_count));
                op_count.clear();
                break;
            case 'D':
                cigar_.emplace_back(Cigar::DELETION, std::stol(op_count));
                op_count.clear();
                break;
            case 'S':
                cigar_.emplace_back(Cigar::CLIPPED, std::stol(op_count));
                op_count.clear();
                break;
            case 'G':
                cigar_.emplace_back(Cigar::NODE_INSERTION, std::stol(op_count));
                op_count.clear();
                break;
            default:
                op_count += c;
        }
    }
}

std::string Cigar::to_string() const {
    std::string cigar_string;

    for (const auto &pair : cigar_) {
        cigar_string += std::to_string(pair.second) + opt_to_char(pair.first);
    }

    return cigar_string;
}

void Cigar::append(Operator op, LengthType num) {
    if (!num)
        return;

    if (cigar_.empty() || cigar_.back().first != op) {
        cigar_.emplace_back(op, num);
    } else {
        cigar_.back().second += num;
    }
}

void Cigar::append(Cigar&& other) {
    if (other.empty())
        return;

    append(other.cigar_.front().first, other.cigar_.front().second);
    cigar_.insert(cigar_.end(), std::next(other.cigar_.begin()), other.cigar_.end());
}

bool Cigar::is_valid(std::string_view reference, std::string_view query) const {
    size_t ref_pos = 0;
    size_t alt_pos = 0;

    for (size_t i = 0; i < cigar_.size(); ++i) {
        const auto &op = cigar_[i];
        if (!op.second) {
            std::cerr << "Empty operation found in CIGAR" << std::endl
                      << to_string() << std::endl
                      << reference << std::endl
                      << query << std::endl;
            return false;
        }

        switch (op.first) {
            case CLIPPED: {
                if (op.second > query.size() - alt_pos) {
                    std::cerr << "Query too short after "
                                << Cigar::opt_to_char(op.first) << std::endl
                                << to_string() << std::endl
                                << reference << std::endl
                                << query << std::endl;
                    return false;
                }

                alt_pos += op.second;
            } break;
            case MATCH:
            case MISMATCH: {
                if (op.second > reference.size() - ref_pos) {
                    std::cerr << "Reference too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (op.second > query.size() - alt_pos) {
                    std::cerr << "Query too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (std::equal(reference.begin() + ref_pos, reference.begin() + ref_pos + op.second,
                               query.begin() + alt_pos)
                        == (op.first != Cigar::MATCH)) {
                    std::cerr << "Mismatch despite MATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_pos += op.second;
                alt_pos += op.second;
            } break;
            case INSERTION: {
                if (i && cigar_[i - 1].first == DELETION) {
                    std::cerr << "INSERTION after DELETION" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (op.second > query.size() - alt_pos) {
                    std::cerr << "Query too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                alt_pos += op.second;
            } break;
            case DELETION: {
                if (i && cigar_[i - 1].first == INSERTION
                        && (ref_pos == reference.size() || reference[ref_pos] != boss::BOSS::kSentinel)) {
                    std::cerr << "DELETION after INSERTION" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (op.second > reference.size() - ref_pos) {
                    std::cerr << "Reference too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_pos += op.second;
            } break;
            case NODE_INSERTION: {
                // do nothing
            } break;
        }
    }

    if (ref_pos != reference.size()) {
        std::cerr << "Reference end not reached" << std::endl
                  << to_string() << std::endl
                  << reference << std::endl
                  << query << std::endl;
        return false;
    }

    if (alt_pos != query.size()) {
        std::cerr << "Query end not reached" << std::endl
                  << to_string() << std::endl
                  << reference << std::endl
                  << query << std::endl;
        return false;
    }

    return true;
}

} // namespace mtg::graph::align
