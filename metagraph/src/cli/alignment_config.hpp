#pragma once

#include <cstddef>
#include <cstdint>

#include "config/config.hpp"


namespace mtg::cli::detail {

template <class AlignerConfig>
AlignerConfig make_alignment_config(const Config &config, size_t min_seed_length) {
    return {
        .num_alternative_paths = config.alignment_num_alternative_paths,
        .min_seed_length = min_seed_length,
        .max_seed_length = config.alignment_max_seed_length,
        .max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus,
        .min_path_score = config.alignment_min_path_score,
        .xdrop = config.alignment_xdrop,
        .min_exact_match = config.alignment_min_exact_match,
        .max_nodes_per_seq_char = config.alignment_max_nodes_per_seq_char,
        .max_ram_per_alignment = config.alignment_max_ram,
        .rel_score_cutoff = config.alignment_rel_score_cutoff,
        .gap_opening_penalty = static_cast<int8_t>(-config.alignment_gap_opening_penalty),
        .gap_extension_penalty = static_cast<int8_t>(-config.alignment_gap_extension_penalty),
        .left_end_bonus = config.alignment_end_bonus,
        .right_end_bonus = config.alignment_end_bonus,
        .forward_and_reverse_complement = !config.align_only_forwards,
        .chain_alignments = config.alignment_chain,
        .two_sided = config.alignment_two_sided,
        .post_chain_alignments = config.alignment_post_chain,
        .seed_complexity_filter = config.alignment_seed_complexity_filter,
        .alignment_edit_distance = config.alignment_edit_distance,
        .alignment_match_score = config.alignment_match_score,
        .alignment_mm_transition_score = config.alignment_mm_transition_score,
        .alignment_mm_transversion_score = config.alignment_mm_transversion_score,
        .score_matrix = typename AlignerConfig::ScoreMatrix{},
    };
}

} // namespace mtg::cli::detail
