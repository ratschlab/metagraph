#include "aligner_cigar.hpp"

#include <iterator>
#include <limits>
#include <stdexcept>

#include "kmer/alphabets.hpp"

namespace mtg {
namespace graph {
namespace detail {

CigarCore::CigarCore(std::string_view cigar_str) {
    LengthType count = 0;
    bool have_count = false;
    for (char c : cigar_str) {
        if (c >= '0' && c <= '9') {
            const auto digit = static_cast<LengthType>(c - '0');
            if (count > (std::numeric_limits<LengthType>::max() - digit) / 10)
                throw std::out_of_range("CIGAR operation length is too large");
            count = count * 10 + digit;
            have_count = true;
            continue;
        }
        const auto op = std::string_view(op_str_).find(c);
        if (op == std::string_view::npos || !have_count || !count)
            throw std::invalid_argument("Invalid CIGAR operation");
        cigar_.emplace_back(static_cast<Operator>(op), count);
        count = 0;
        have_count = false;
    }
    if (have_count)
        throw std::invalid_argument("CIGAR length has no operation");
}

std::string CigarCore::to_string() const {
    std::string cigar_string;

    for (const auto &pair : cigar_) {
        cigar_string += std::to_string(pair.second) + opt_to_char(pair.first);
    }

    return cigar_string;
}

void CigarCore::append(Operator op, LengthType num) {
    if (!num)
        return;

    if (cigar_.empty() || cigar_.back().first != op) {
        cigar_.emplace_back(op, num);
    } else {
        cigar_.back().second += num;
    }
}

void CigarCore::append_cigar(CigarCore&& other) {
    if (other.empty())
        return;

    append(other.cigar_.front().first, other.cigar_.front().second);
    cigar_.insert(cigar_.end(), std::next(other.cigar_.begin()), other.cigar_.end());
}

} // namespace detail

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
        const auto *alphabet = kmer::alphabets::kAlphabetDNA5;
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
        if (alphabet_encoding[c] == alphabet_encoding[0]
#if _DNA5_GRAPH
                && c != 'N'
#endif
                )
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
            case CLIPPED: {
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
            case MATCH:
            case MISMATCH: {
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

                if (!std::equal(ref_it, ref_it + op.second, alt_it, [&](char ref, char alt) {
                        return (ref == alt) == (op.first == Cigar::MATCH);
                    })) {
                    std::cerr << "Character disagrees with CIGAR operation" << std::endl
                              << to_string() << std::endl
                              << reference << std::endl
                              << query << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case INSERTION: {
                if (i && cigar_[i - 1].first == DELETION) {
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
            case DELETION: {
                if (i && cigar_[i - 1].first == INSERTION) {
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
            case NODE_INSERTION: {
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
