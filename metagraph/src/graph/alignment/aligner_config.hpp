#ifndef __ALIGNER_CONFIG_HPP__
#define __ALIGNER_CONFIG_HPP__

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string_view>


namespace mtg {
namespace graph {
namespace align {

class Cigar;

struct DBGAlignerConfig {
    typedef int32_t score_t;
    typedef std::array<int8_t, 128> ScoreMatrixRow;
    typedef std::array<ScoreMatrixRow, 128> ScoreMatrix;

    size_t num_alternative_paths = 1;
    size_t min_seed_length = 0;
    size_t max_seed_length = 0;
    size_t max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
    size_t row_batch_size = std::numeric_limits<size_t>::max();

    // Lowest possible score. 100 is added to prevent underflow during operations.
    // For this to work, all penalties should be less than 100.
    // This is checked whenever an aligner is initialized.
    constexpr static score_t ninf = std::numeric_limits<score_t>::min() + 100;

    // thresholds for scores
    score_t min_cell_score = ninf;
    score_t min_path_score = 0;
    score_t xdrop = std::numeric_limits<score_t>::max();

    double min_exact_match = 0.0;
    double max_nodes_per_seq_char = std::numeric_limits<double>::max();
    double max_ram_per_alignment = std::numeric_limits<double>::max();
    double rel_score_cutoff = 0.0;

    int8_t gap_opening_penalty = -5;
    int8_t gap_extension_penalty = -2;
    int8_t left_end_bonus = 0;
    int8_t right_end_bonus = 0;

    bool forward_and_reverse_complement = true;
    bool chain_alignments = false;
    bool post_chain_alignments = false;
    bool global_xdrop = true;
    bool allow_left_trim = true;
    bool no_backtrack = false;

    bool alignment_edit_distance;
    int8_t alignment_match_score;
    int8_t alignment_mm_transition_score;
    int8_t alignment_mm_transversion_score;

    ScoreMatrix score_matrix;

    void print_summary() const;

    score_t score_sequences(std::string_view a, std::string_view b) const {
        return std::inner_product(
            a.begin(), a.end(), b.begin(), score_t(0), std::plus<score_t>(),
            [&](char a, char b) -> score_t { return score_matrix[a][b]; }
        );
    }

    score_t match_score(std::string_view query) const {
        return score_sequences(query, query);
    }

    score_t score_cigar(std::string_view reference,
                        std::string_view query,
                        const Cigar &cigar) const;

    bool check_config_scores() const;

    void set_scoring_matrix();

    // Protein matrices
    static const ScoreMatrix score_matrix_blosum62;

    static ScoreMatrix dna_scoring_matrix(int8_t match_score,
                                          int8_t mm_transition_score,
                                          int8_t mm_transversion_score);

    static ScoreMatrix unit_scoring_matrix(int8_t match_score,
                                           const std::string &alphabet,
                                           const uint8_t *encoding);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_CONFIG_HPP__
