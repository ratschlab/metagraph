#ifndef __SCORE_KMER_PRESENCE_MASK_HPP__
#define __SCORE_KMER_PRESENCE_MASK_HPP__

#include <cstdint>
#include <sdsl/int_vector.hpp>


namespace mtg {
namespace common {

int32_t score_kmer_presence_mask(size_t k,
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score = 1,
                                 int32_t mismatch_score = 2);

} // namespace common
} // namespace mtg

#endif // __SCORE_KMER_PRESENCE_MASK_HPP__
