#include "aln_cigar.hpp"

#include "graph/representation/succinct/boss.hpp"

namespace mtg::graph::align_redone {


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
                if (i && i + 1 != cigar_.size()) {
                    // Post-chaining skips unaligned query bases beside a
                    // graph jump, represented by G or a deleted sentinel.
                    const auto [previous, previous_length] = cigar_[i - 1];
                    const auto [next, next_length] = cigar_[i + 1];
                    const bool graph_jump = previous == NODE_INSERTION || next == NODE_INSERTION
                        || (previous == DELETION && previous_length <= ref_pos
                            && reference.substr(ref_pos - previous_length, previous_length)
                                .find(boss::BOSS::kSentinel) != std::string_view::npos)
                        || (next == DELETION && next_length <= reference.size() - ref_pos
                            && reference.substr(ref_pos, next_length)
                                .find(boss::BOSS::kSentinel) != std::string_view::npos);
                    if (!graph_jump) {
                        std::cerr << "Internal clipping without a graph jump in CIGAR" << std::endl;
                        return false;
                    }
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

                if (!std::equal(reference.begin() + ref_pos, reference.begin() + ref_pos + op.second,
                                query.begin() + alt_pos, [&](char ref, char alt) {
                                    return (ref == alt) == (op.first == Cigar::MATCH);
                                })) {
                    std::cerr << "Character disagrees with CIGAR operation" << std::endl
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

} // namespace mtg::graph::align_redone
