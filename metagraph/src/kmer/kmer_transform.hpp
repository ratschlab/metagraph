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
    return KMER(result);
}

// construct word `...1001001001` with k ones for lifting k-mers
template <typename WordType>
inline WordType get_sentinel_delta(size_t char_width, size_t k) {
    assert(char_width * k <= sizeof(WordType) * 8);
    WordType word = 0;
    for (size_t i = 0; i < k; ++i) {
        word <<= char_width;
        word |= 1;
    }
    return word;
}

// shift to the next dummy sink and add +1 to each character of the k-mer
template <typename KMER_TO, typename KMER_FROM>
inline typename KMER_TO::WordType get_sink_and_lift(const KMER_FROM &kmer, size_t k) {
    static constexpr int L1 = KMER_FROM::kBitsPerChar;
    static constexpr int L2 = KMER_TO::kBitsPerChar;
    static_assert(L2 >= L1);
    static_assert(L2 <= L1 + 1);
    assert(sizeof(typename KMER_TO::WordType)
                   >= sizeof(typename KMER_FROM::WordType));

    static constexpr uint64_t first_char_mask_1 = (1ull << L1) - 1;

    typename KMER_TO::WordType word = (kmer.data() & first_char_mask_1) + 1;

    for (int pos = L1 * (k - 1); pos >= L1 * 2; pos -= L1) {
        word <<= L2;
        assert(kmer[pos / L1] + 1 <= sdsl::bits::lo_set[L2]);
        word |= (static_cast<uint64_t>(kmer.data() >> pos) & first_char_mask_1) + 1;
    }

    word <<= L2;

    return word;
}

// transforms k-mer to the new character width
template <typename KMER_TO, typename KMER_FROM>
inline __attribute__((always_inline))
typename KMER_TO::WordType transform(const KMER_FROM &kmer, size_t k) {
    static constexpr size_t L1 = KMER_FROM::kBitsPerChar;
    static constexpr size_t L2 = KMER_TO::kBitsPerChar;
    static_assert(L2 >= L1);
    static_assert(L2 <= L1 + 1);
    assert(sizeof(typename KMER_TO::WordType) >= sizeof(typename KMER_FROM::WordType));

    if constexpr(L1 == L2) {
        return kmer.data();

    } else {
        typename KMER_TO::WordType word = 0;

        static constexpr uint64_t char_mask = (1ull << L1) - 1;

        for (int pos = L1 * (k - 1); pos >= 0; pos -= L1) {
            word <<= L2;
            assert(kmer[pos / L1] + 1 <= sdsl::bits::lo_set[L2]);
            word |= static_cast<uint64_t>(kmer.data() >> pos) & char_mask;
        }

        return word;
    }
}

// Expand each of a byte's four 2-bit characters into a 3-bit field.
// A raw array avoids GCC's missed optimization for std::array's size_t index.
// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=127883.
static constexpr auto lookup_2_3_storage = [] {
    struct { uint16_t values[256]; } table{};
    for (size_t value = 0; value < 256; ++value) {
        for (size_t character = 0; character < 4; ++character)
            table.values[value] |= ((value >> (2 * character)) & 3) << (3 * character);
    }
    return table;
}();
static constexpr auto &lookup_2_3 = lookup_2_3_storage.values;

template <>
inline __attribute__((always_inline)) sdsl::uint128_t
transform<kmer::KMerBOSS<sdsl::uint128_t, 3>, kmer::KMerBOSS<uint64_t, 2>>(
        const kmer::KMerBOSS<uint64_t, 2> &kmer, size_t /*k*/) {
    const uint8_t *kmer_char = reinterpret_cast<const uint8_t *>(&kmer);
    // transform 64-bit kmer to 128 bits
    return static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[7]]) << 84
            | static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[6]]) << 72
            | static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[5]]) << 60
            | static_cast<uint64_t>(lookup_2_3[kmer_char[4]]) << 48
            | static_cast<uint64_t>(lookup_2_3[kmer_char[3]]) << 36
            | static_cast<uint64_t>(lookup_2_3[kmer_char[2]]) << 24
            | static_cast<uint64_t>(lookup_2_3[kmer_char[1]]) << 12
            | lookup_2_3[kmer_char[0]];
}

template <>
inline __attribute__((always_inline)) sdsl::uint256_t
transform<kmer::KMerBOSS<sdsl::uint256_t, 3>, kmer::KMerBOSS<sdsl::uint128_t, 2>>(
        const kmer::KMerBOSS<sdsl::uint128_t, 2> &kmer, size_t /*k*/) {
    const uint8_t *kmer_char = reinterpret_cast<const uint8_t *>(&kmer);
    // transform 128-bit kmer to 256 bits
    return sdsl::uint256_t(0, static_cast<sdsl::uint128_t>(lookup_2_3[kmer_char[15]]) << 52)
            | sdsl::uint256_t(static_cast<uint64_t>(lookup_2_3[kmer_char[14]]) << 48
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[13]]) << 36
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[12]]) << 24
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[11]]) << 12
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[10]])
            ) << 120
            | sdsl::uint128_t(static_cast<uint64_t>(lookup_2_3[kmer_char[9]]) << 48
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[8]]) << 36
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[7]]) << 24
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[6]]) << 12
                              | static_cast<uint64_t>(lookup_2_3[kmer_char[5]])
            ) << 60
            | (static_cast<uint64_t>(lookup_2_3[kmer_char[4]]) << 48
            | static_cast<uint64_t>(lookup_2_3[kmer_char[3]]) << 36
            | static_cast<uint64_t>(lookup_2_3[kmer_char[2]]) << 24
            | static_cast<uint64_t>(lookup_2_3[kmer_char[1]]) << 12
            | lookup_2_3[kmer_char[0]]);
}

} // namespace kmer
} // namespace mtg
