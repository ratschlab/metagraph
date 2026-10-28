#include "aligner_config.hpp"

#include "aln_cigar.hpp"


namespace mtg::graph::align_redone {

DBGAlignerConfig::score_t DBGAlignerConfig
::score_cigar(std::string_view reference,
              std::string_view query,
              const Cigar &cigar) const {
    assert(cigar.is_valid(reference, query));

    if (cigar.empty())
        return 0;

    score_t score = (!cigar.get_clipping() ? left_end_bonus : 0)
                    + (!cigar.get_end_clipping() ? right_end_bonus : 0);

    auto ref_it = reference.begin();
    auto alt_it = query.begin();
    auto it = cigar.data().begin();
    if (it->first == Cigar::CLIPPED) {
        alt_it += it->second;
        ++it;
    }

    for ( ; it != cigar.data().end(); ++it) {
        const auto &op = *it;
        switch (op.first) {
            case Cigar::CLIPPED: {
                alt_it += op.second;
            } break;
            case Cigar::MATCH: {
                score += match_score(std::string_view(ref_it, op.second));
                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Cigar::MISMATCH: {
                score += score_sequences(std::string_view(ref_it, op.second),
                                         std::string_view(alt_it, op.second));
                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Cigar::INSERTION: {
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
                alt_it += op.second;
            } break;
            case Cigar::DELETION: {
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
                ref_it += op.second;
                if (it >= cigar.data().begin() + 2 && (it - 2)->first == Cigar::DELETION
                        && (it - 1)->first == Cigar::NODE_INSERTION) {
                    score -= gap_opening_penalty - gap_extension_penalty;
                }
            } break;
            case Cigar::NODE_INSERTION: {
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
            } break;
        }
    }

    assert(ref_it == reference.end());
    assert(alt_it == query.end());

    return score;
}

} // namespace mtg::graph::align_redone
