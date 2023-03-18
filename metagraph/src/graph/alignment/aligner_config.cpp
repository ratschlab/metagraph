#include "aligner_config.hpp"

#include "aligner_cigar.hpp"
#include "kmer/alphabets.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

void DBGAlignerConfig::print_summary() const {
    logger->trace("Alignment settings:");
    logger->trace("\t Alignments to report: {}", num_alternative_paths);
    logger->trace("\t Min seed length: {}", min_seed_length);
    logger->trace("\t Max seed length: {}", max_seed_length);
    logger->trace("\t Max num seeds per locus: {}", max_num_seeds_per_locus);
    logger->trace("\t Max num nodes per sequence char: {}", max_nodes_per_seq_char);
    logger->trace("\t Max RAM per alignment: {}", max_ram_per_alignment);
    logger->trace("\t Gap opening penalty: {}", (int64_t)gap_opening_penalty);
    logger->trace("\t Gap extension penalty: {}", (int64_t)gap_extension_penalty);
    logger->trace("\t Min alignment score: {}", min_path_score);
    logger->trace("\t X drop-off: {}", xdrop);
    logger->trace("\t Allow node insertion: {}", allow_jump);
    logger->trace("\t Path insertion penalty: {}", node_insertion_penalty);
    logger->trace("\t Exact nucleotide match threshold: {}", min_exact_match);
    logger->trace("\t Chain alignments: {}", chain_alignments || allow_jump || allow_label_change);

    logger->trace("\t Scoring matrix: {}", alignment_edit_distance ? "unit costs" : "matrix");
    if (!alignment_edit_distance) {
        logger->trace("\t\t Match score: {}", (int64_t)alignment_match_score);
        logger->trace("\t\t (DNA) Transition score: {}", (int64_t)alignment_mm_transition_score);
        logger->trace("\t\t (DNA) Transversion score: {}", (int64_t)alignment_mm_transversion_score);
    }
}

// check to make sure the current scoring system won't underflow
bool DBGAlignerConfig::check_config_scores() const {
    int8_t min_penalty_score = std::numeric_limits<int8_t>::max();
    for (const auto &row : score_matrix) {
        min_penalty_score = std::min(min_penalty_score,
                                     *std::min_element(row.begin(), row.end()));
    }

    if (gap_opening_penalty * 2 >= min_penalty_score) {
        common::logger->error(
            "|gap_opening_penalty| * 2 should be greater than greatest mismatch penalty: "
            "{} >= {}",
            gap_opening_penalty * 2, (score_t)min_penalty_score
        );
        return false;
    }

    min_penalty_score = std::min({ min_penalty_score,
                                   gap_opening_penalty,
                                   gap_extension_penalty });

    assert(min_penalty_score < 0 && "min scores must be negative");

    if (min_cell_score >= std::numeric_limits<score_t>::min() - min_penalty_score)
        return true;

    std::cerr << "min_cell_score is too small: " << min_cell_score << std::endl;
    return false;
}

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
    if (it->first == Cigar::CLIPPED)
        ++it;

    for ( ; it != cigar.data().end(); ++it) {
        const auto &op = *it;
        switch (op.first) {
            case Cigar::CLIPPED: {
                if (it + 1 != cigar.data().end())
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
            case Cigar::NODE_INSERTION: { score += node_insertion_penalty; } break;
        }
    }

    assert(ref_it == reference.end());
    assert(alt_it == query.end());

    return score;
}

std::vector<DBGAlignerConfig::score_t> DBGAlignerConfig
::get_per_char_scores(std::string_view reference,
                      std::string_view query,
                      const Cigar &cigar) const {
    assert(cigar.is_valid(reference, query));

    if (cigar.empty())
        return {};

    std::vector<score_t> scores;
    scores.reserve(query.size() - cigar.get_clipping() - cigar.get_end_clipping());

    score_t score = !cigar.get_clipping() ? left_end_bonus : 0;
    score_t last_score = 0;

    auto ref_it = reference.begin();
    auto alt_it = query.begin();
    auto it = cigar.data().begin();
    if (it->first == Cigar::CLIPPED)
        ++it;

    for ( ; it != cigar.data().end(); ++it) {
        const auto &op = *it;
        switch (op.first) {
            case Cigar::CLIPPED: {
                if (it + 1 != cigar.data().end())
                    alt_it += op.second;
            } break;
            case Cigar::MATCH: {
                for (size_t i = 0; i < op.second; ++i) {
                    score += match_score(std::string_view(ref_it, 1));
                    scores.emplace_back(score - last_score);
                    last_score = score;
                    ++ref_it;
                    ++alt_it;
                }
            } break;
            case Cigar::MISMATCH: {
                for (size_t i = 0; i < op.second; ++i) {
                    score += score_sequences(std::string_view(ref_it, 1),
                                             std::string_view(alt_it, 1));
                    scores.emplace_back(score - last_score);
                    last_score = score;
                    ++ref_it;
                    ++alt_it;
                }
            } break;
            case Cigar::INSERTION: {
                scores.emplace_back(gap_opening_penalty);
                for (size_t i = 1; i < op.second; ++i) {
                    scores.emplace_back(gap_extension_penalty);
                }
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
                last_score = score;
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
            case Cigar::NODE_INSERTION: { score += node_insertion_penalty; } break;
        }
    }

    scores.back() += !cigar.get_end_clipping() ? right_end_bonus : 0;

    assert(ref_it == reference.end());
    assert(alt_it == query.end());
    assert(scores.size() == query.size() - cigar.get_clipping() - cigar.get_end_clipping());

    return scores;
}

void DBGAlignerConfig::set_scoring_matrix() {
    if (alignment_edit_distance) {
        // TODO: REPLACE THIS
        #if _PROTEIN_GRAPH
            const auto *alphabet = kmer::alphabets::kAlphabetProtein;
            const auto *alphabet_encoding = kmer::alphabets::kCharToProtein;
        #elif _DNA_GRAPH || _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
            const auto *alphabet = kmer::alphabets::kAlphabetDNA;
            const auto *alphabet_encoding = kmer::alphabets::kCharToDNA;
        #else
            static_assert(false,
                "Define an alphabet: either "
                "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
            );
        #endif

        score_matrix = unit_scoring_matrix(1, alphabet, alphabet_encoding);
        left_end_bonus = 0;
        right_end_bonus = 0;

    } else {
        #if _PROTEIN_GRAPH
            score_matrix = score_matrix_blosum62;
        #elif _DNA_GRAPH || _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
            score_matrix = dna_scoring_matrix(alignment_match_score,
                                              -alignment_mm_transition_score,
                                              -alignment_mm_transversion_score);
        #else
            static_assert(false,
                "Define an alphabet: either "
                "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
            );
        #endif
    }
}

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::dna_scoring_matrix(int8_t match_score,
                     int8_t mm_transition_score,
                     int8_t mm_transversion_score) {
    ScoreMatrix score_matrix;
    for (auto& row : score_matrix) {
        row.fill(mm_transversion_score);
    }

    score_matrix['A']['G'] = mm_transition_score;
    score_matrix['G']['A'] = mm_transition_score;
    score_matrix['C']['T'] = mm_transition_score;
    score_matrix['T']['C'] = mm_transition_score;
    score_matrix['A']['A'] = match_score;
    score_matrix['C']['C'] = match_score;
    score_matrix['G']['G'] = match_score;
    score_matrix['T']['T'] = match_score;

    return score_matrix;
}

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::unit_scoring_matrix(int8_t match_score,
                      const std::string &alphabet,
                      const uint8_t *encoding) {
    ScoreMatrix score_matrix;
    for (auto& row : score_matrix) {
        row.fill(-match_score);
    }

    for (uint8_t c : alphabet) {
        // if a character is invalid, don't count matches of that character
        if (encoding[c] == encoding[0])
            continue;

        char upper = toupper(c);

        score_matrix[upper][upper] = match_score;
    }

    return score_matrix;
}

DBGAlignerConfig::ScoreMatrix blosum62_scoring_matrix() {
    std::string alphabet = "ARNDCQEGHILKMFPSTWYVBZX";

    std::vector<std::vector<int8_t>> scores = {
        {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0 },
        { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1 },
        { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1 },
        { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1 },
        {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2 },
        { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1 },
        { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },
        {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1 },
        { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1 },
        { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1 },
        { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1 },
        { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1 },
        { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1 },
        { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1 },
        { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2 },
        {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0 },
        {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0 },
        { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2 },
        { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1 },
        {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1 },
        { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1 },
        { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },
        {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1 }
    };

    DBGAlignerConfig::ScoreMatrix score_matrix;

    for (size_t i = 0; i < score_matrix.size(); ++i) {
        score_matrix[i].fill(-4);

        // meant to handle the letters J, O, U
        score_matrix[i][i] = 1;
    }

    for (size_t i = 0; i < alphabet.size(); ++i) {
        for (size_t j = 0; j < alphabet.size(); ++j) {
            score_matrix[alphabet[i]][alphabet[j]] = scores[i][j];
        }
    }

    return score_matrix;
}

const DBGAlignerConfig::ScoreMatrix DBGAlignerConfig::score_matrix_blosum62
    = blosum62_scoring_matrix();

} // namespace align
} // namespace graph
} // namespace mtg
