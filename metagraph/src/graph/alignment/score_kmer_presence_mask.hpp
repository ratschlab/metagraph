#ifndef __SCORE_KMER_PRESENCE_MASK_HPP__
#define __SCORE_KMER_PRESENCE_MASK_HPP__

#include <cstdint>
#include <sdsl/int_vector.hpp>


namespace mtg {
namespace graph {
namespace align {

/**
 * Score a k-mer presence mask using the BIGSI penalty function.
 *
 * Treats each run of present k-mers as a match block and each run of absent
 * k-mers as a putative-SNP penalty, then normalizes by sequence length.
 *
 * @param k                  k-mer length used to recover sequence length from the mask.
 * @param kmer_presence_mask bit per k-mer position; 1 = present in graph, 0 = absent.
 * @param match_score        score added per matching k-mer.
 * @param mismatch_score     penalty per inferred SNP in a non-matching run.
 * @return                   non-negative score, scaled to sequence length.
 */
int32_t score_kmer_presence_mask(size_t k,
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score = 1,
                                 int32_t mismatch_score = 2);

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __SCORE_KMER_PRESENCE_MASK_HPP__
