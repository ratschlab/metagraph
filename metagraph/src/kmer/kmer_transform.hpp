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

static const uint16_t lookup_2_3[256] = {
        0,    1,    2,    3,    8,    9,    10,   11,   16,   17,   18,   19,   24,
        25,   26,   27,   64,   65,   66,   67,   72,   73,   74,   75,   80,   81,
        82,   83,   88,   89,   90,   91,   128,  129,  130,  131,  136,  137,  138,
        139,  144,  145,  146,  147,  152,  153,  154,  155,  192,  193,  194,  195,
        200,  201,  202,  203,  208,  209,  210,  211,  216,  217,  218,  219,  512,
        513,  514,  515,  520,  521,  522,  523,  528,  529,  530,  531,  536,  537,
        538,  539,  576,  577,  578,  579,  584,  585,  586,  587,  592,  593,  594,
        595,  600,  601,  602,  603,  640,  641,  642,  643,  648,  649,  650,  651,
        656,  657,  658,  659,  664,  665,  666,  667,  704,  705,  706,  707,  712,
        713,  714,  715,  720,  721,  722,  723,  728,  729,  730,  731,  1024, 1025,
        1026, 1027, 1032, 1033, 1034, 1035, 1040, 1041, 1042, 1043, 1048, 1049, 1050,
        1051, 1088, 1089, 1090, 1091, 1096, 1097, 1098, 1099, 1104, 1105, 1106, 1107,
        1112, 1113, 1114, 1115, 1152, 1153, 1154, 1155, 1160, 1161, 1162, 1163, 1168,
        1169, 1170, 1171, 1176, 1177, 1178, 1179, 1216, 1217, 1218, 1219, 1224, 1225,
        1226, 1227, 1232, 1233, 1234, 1235, 1240, 1241, 1242, 1243, 1536, 1537, 1538,
        1539, 1544, 1545, 1546, 1547, 1552, 1553, 1554, 1555, 1560, 1561, 1562, 1563,
        1600, 1601, 1602, 1603, 1608, 1609, 1610, 1611, 1616, 1617, 1618, 1619, 1624,
        1625, 1626, 1627, 1664, 1665, 1666, 1667, 1672, 1673, 1674, 1675, 1680, 1681,
        1682, 1683, 1688, 1689, 1690, 1691, 1728, 1729, 1730, 1731, 1736, 1737, 1738,
        1739, 1744, 1745, 1746, 1747, 1752, 1753, 1754, 1755
};

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
