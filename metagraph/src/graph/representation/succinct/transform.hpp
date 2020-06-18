#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "common/utils/template_utils.hpp"
#include "kmer/kmer_extractor.hpp"

namespace mtg {
namespace kmer {

using TAlphabet = KmerExtractorBOSS::TAlphabet;
template <typename T>
inline T
reverse_complement(size_t k, const T &v, const std::vector<TAlphabet> &complement_code) {
    using KMER = utils::get_first_type_t<T>;
    using INT = typename KMER::WordType;
    INT kmer = utils::get_first(v).data();
    constexpr uint64_t mask = KMER::kFirstCharMask;
    INT last_two_chars = complement_code[static_cast<TAlphabet>(kmer & mask)];
    kmer >>= KMER::kBitsPerChar;
    last_two_chars = (last_two_chars << KMER::kBitsPerChar)
            | complement_code[static_cast<TAlphabet>(kmer & mask)];
    kmer >>= KMER::kBitsPerChar;
    INT result = 0;
    for (uint32_t i = 2; i < k; ++i) {
        TAlphabet next_char = kmer & mask;
        assert(next_char >= 0 && next_char < complement_code.size());
        result = (result << KMER::kBitsPerChar) | complement_code[next_char];
        kmer >>= KMER::kBitsPerChar;
    }
    result = (result << 2 * KMER::kBitsPerChar) | last_two_chars;
    if constexpr (utils::is_pair_v<T>) {
        return T(KMER(result), v.second);
    } else {
        return KMER(result);
    }
}
} // namespace kmer
} // namespace mtg
