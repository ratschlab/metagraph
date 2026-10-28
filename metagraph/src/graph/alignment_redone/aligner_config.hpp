#pragma once

#include <string_view>

#include "graph/alignment/aligner_config.hpp"


namespace mtg::graph::align_redone {

class Cigar;

// The new engine shares the scoring and search options of the existing aligner
// and adds the parameters of its anchor chaining.
struct DBGAlignerConfig : public align::DBGAlignerConfig {
    size_t max_dist_between_seeds = 400;
    size_t gap_shrinking_factor = 4;
    double chaining_algorithm_switch_cutoff = 1.0;
    size_t chaining_bandwidth = 65;

    score_t score_cigar(std::string_view reference,
                        std::string_view query,
                        const Cigar &cigar) const;
};

} // namespace mtg::graph::align_redone
