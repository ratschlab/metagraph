#include "aligner_cigar.hpp"

#include "kmer/alphabets.hpp"

namespace mtg {
namespace graph {
namespace align {


OperatorTable initialize_opt_table() {
    OperatorTable char_to_op;

    // TODO: fix this when alphabets are no longer set at compile time
    #if _PROTEIN_GRAPH
        const auto *alphabet = kmer::alphabets::kAlphabetProtein;
        const auto *alphabet_encoding = kmer::alphabets::kCharToProtein;
    #elif _DNA_CASE_SENSITIVE_GRAPH
        const auto *alphabet = kmer::alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = kmer::alphabets::kCharToDNA;
    #elif _DNA5_GRAPH
        const auto *alphabet = kmer::alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = kmer::alphabets::kCharToDNA;
    #elif _DNA_GRAPH
        const auto *alphabet = kmer::alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = kmer::alphabets::kCharToDNA;
    #else
        static_assert(false,
            "Define an alphabet: either "
            "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
        );
    #endif

    for (auto &row : char_to_op) {
        row.fill(Cigar::MISMATCH);
    }

    for (uint8_t c : std::string(alphabet)) {
        if (alphabet_encoding[c] == alphabet_encoding[0])
            continue;

        char upper = toupper(c);
        char lower = tolower(c);

        char_to_op[upper][upper]
            = char_to_op[upper][lower]
            = char_to_op[lower][upper]
            = char_to_op[lower][lower] = Cigar::MATCH;
    }

    return char_to_op;
}

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
    auto ref_it = reference.begin();
    auto alt_it = query.begin();

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
            case Operator::CLIPPED: {
                if ((ref_it != reference.begin() || alt_it != query.begin())
                        && (ref_it != reference.end() || alt_it != query.end())) {
                    if (alt_it > query.end() - op.second) {
                        std::cerr << "Query too short after "
                                  << Cigar::opt_to_char(op.first) << std::endl
                                  << to_string() << std::endl
                                  << reference << std::endl
                                  << query << std::endl;
                        return false;
                    }

                    alt_it += op.second;
                }
            } break;
            case Operator::MATCH:
            case Operator::MISMATCH: {
                if (ref_it > reference.end() - op.second) {
                    std::cerr << "Reference too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (alt_it > query.end() - op.second) {
                    std::cerr << "Query too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (std::equal(ref_it, ref_it + op.second, alt_it)
                        == (op.first != Cigar::MATCH)) {
                    std::cerr << "Mismatch despite MATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Operator::INSERTION: {
                if (i && cigar_[i - 1].first == Operator::DELETION) {
                    std::cerr << "INSERTION after DELETION" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (alt_it > query.end() - op.second) {
                    std::cerr << "Query too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                alt_it += op.second;
            } break;
            case Operator::DELETION: {
                if (i && cigar_[i - 1].first == Operator::INSERTION) {
                    std::cerr << "DELETION after INSERTION" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                if (ref_it > reference.end() - op.second) {
                    std::cerr << "Reference too short after "
                              << Cigar::opt_to_char(op.first) << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
            } break;
            case Operator::NODE_INSERTION: {
                // do nothing
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

} // namespace align
} // namespace graph
} // namespace mtg
