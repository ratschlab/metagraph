#include "aligner_helper.hpp"


Cigar::Cigar(const std::string &cigar_str) {
    std::string op_count;
    for (auto c : cigar_str) {
        switch (c) {
            case '=':
                cigar_.emplace_back(Cigar::Operator::MATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'X':
                cigar_.emplace_back(Cigar::Operator::MISMATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'I':
                cigar_.emplace_back(Cigar::Operator::INSERTION, std::stol(op_count));
                op_count.clear();
                break;
            case 'D':
                cigar_.emplace_back(Cigar::Operator::DELETION, std::stol(op_count));
                op_count.clear();
                break;
            case 'S':
                cigar_.emplace_back(Cigar::Operator::CLIPPED, std::stol(op_count));
                op_count.clear();
                break;
            default:
                op_count += c;
        }
    }
}

Cigar::OperatorTable Cigar::char_to_op;

void Cigar::initialize_opt_table(const std::string &alphabet, const uint8_t *encoding) {
    for (auto& row : char_to_op) {
        row.fill(Cigar::Operator::MISMATCH);
    }

    for (uint8_t c : alphabet) {
        if (encoding[c] == encoding[0])
            continue;

        char upper = toupper(c);
        char lower = tolower(c);

        char_to_op[upper][upper]
            = char_to_op[upper][lower]
            = char_to_op[lower][upper]
            = char_to_op[lower][lower] = Cigar::Operator::MATCH;
    }
}

char opt_to_char(Cigar::Operator op) {
    switch (op) {
        case Cigar::Operator::MATCH: return '=';
        case Cigar::Operator::MISMATCH: return 'X';
        case Cigar::Operator::INSERTION: return 'I';
        case Cigar::Operator::DELETION: return 'D';
        case Cigar::Operator::CLIPPED: return 'S';
    }

    assert(false);
    return '\0';
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

bool Cigar::is_valid(const std::string_view reference,
                     const std::string_view query) const {
    auto ref_it = reference.begin();
    auto alt_it = query.begin();

    for (const auto &op : cigar_) {
        if (!op.second) {
            std::cerr << "Empty operation found in CIGAR" << std::endl
                      << to_string() << std::endl
                      << reference << std::endl
                      << query << std::endl;
            return false;
        }

        switch (op.first) {
            case Operator::CLIPPED: {
                if ((ref_it != reference.begin() || alt_it != query.begin())
                        && (ref_it != reference.end() || alt_it != query.end())) {
                    std::cerr << "Internal clipping found in CIGAR" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }
            } break;
            case Operator::MATCH: {
                if (ref_it > reference.end() - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (alt_it > query.end() - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (strncmp(ref_it, alt_it, op.second)) {
                    std::cerr << "Mismatch despite MATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Operator::MISMATCH: {
                if (ref_it > reference.end() - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (alt_it > query.end() - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (std::mismatch(ref_it, ref_it + op.second,
                                  alt_it, alt_it + op.second).first != ref_it) {
                    std::cerr << "Match despite MISMATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Operator::INSERTION: {
                if (alt_it > query.end() - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                alt_it += op.second;
            } break;
            case Operator::DELETION: {
                if (ref_it > reference.end() - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
            } break;
        }
    }

    if (ref_it != reference.end()) {
        std::cerr << "Reference end not reached" << std::endl
                  << to_string() << std::endl
                  << reference << std::endl
                  << query << std::endl;
        return false;
    }

    if (alt_it != query.end()) {
        std::cerr << "Query end not reached" << std::endl
                  << to_string() << std::endl
                  << reference << std::endl
                  << query << std::endl;
        return false;
    }

    return true;
}
