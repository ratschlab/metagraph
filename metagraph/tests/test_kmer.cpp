#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#define private public
#define protected public

#include "kmer/kmer.hpp"
#include "kmer/kmer_extractor.hpp"


namespace {

using namespace mtg::kmer;

typedef uint8_t TAlphabet;

const size_t kBitsPerChar = KmerExtractor2Bit::bits_per_char;

const KmerExtractor2Bit kmer_extractor;


template <typename IntType>
std::string kmer_packed_codec(const std::string &test_kmer) {
    std::vector<TAlphabet> kmer(test_kmer.size());
    std::transform(test_kmer.begin(), test_kmer.end(), kmer.begin(),
        [](char c) { return kmer_extractor.encode(c); }
    );
    return KMer<IntType, kBitsPerChar>(kmer).to_string(test_kmer.length(),
                                                       kmer_extractor.alphabet);
}

template <typename IntType>
void test_kmer_packed_codec(const std::string &test_kmer,
                            const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_packed_codec<IntType>(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_packed_codec<IntType>(test_kmer));
}

// Workaround to avoid the issues with the uint128_t type traits undefined in AppleClang
struct UINT64_ { typedef uint64_t type; };
struct UINT128_ { typedef sdsl::uint128_t type; };
struct UINT256_ { typedef sdsl::uint256_t type; };
typedef ::testing::Types<UINT64_,
                         UINT128_,
                         UINT256_> IntTypes;
template <class IntType>
class Kmer : public ::testing::Test { };
TYPED_TEST_SUITE(Kmer, IntTypes);

template <class IntType>
class Integer : public ::testing::Test { };
TYPED_TEST_SUITE(Integer, IntTypes);

TYPED_TEST(Integer, compare0) {
    typename TypeParam::type value = 0;

    ASSERT_TRUE(value < 7ull);
    ASSERT_TRUE(value < 6ull);
    ASSERT_TRUE(value < 3ull);
    ASSERT_TRUE(value < 1ull);
    ASSERT_FALSE(value < 0ull);

    ASSERT_FALSE(value == 7ull);
    ASSERT_FALSE(value == 6ull);
    ASSERT_FALSE(value == 3ull);
    ASSERT_FALSE(value == 1ull);
    ASSERT_TRUE(value == 0ull);

    ASSERT_FALSE(value > 7ull);
    ASSERT_FALSE(value > 6ull);
    ASSERT_FALSE(value > 3ull);
    ASSERT_FALSE(value > 1ull);
    ASSERT_FALSE(value > 0ull);
}

TYPED_TEST(Integer, compare7) {
    typename TypeParam::type value = 7;

    ASSERT_TRUE(value >= 7ull);
    ASSERT_TRUE(value >= 6ull);
    ASSERT_TRUE(value >= 3ull);
    ASSERT_TRUE(value >= 1ull);
    ASSERT_TRUE(value >= 0ull);
}

TYPED_TEST(Integer, int_compare0) {
    typename TypeParam::type value = 0;

    ASSERT_TRUE(7ull > value);
    ASSERT_TRUE(6ull > value);
    ASSERT_TRUE(3ull > value);
    ASSERT_TRUE(1ull > value);
    ASSERT_FALSE(0ull > value);

    ASSERT_FALSE(7ull == value);
    ASSERT_FALSE(6ull == value);
    ASSERT_FALSE(3ull == value);
    ASSERT_FALSE(1ull == value);
    ASSERT_TRUE(0ull == value);

    ASSERT_FALSE(7ull < value);
    ASSERT_FALSE(6ull < value);
    ASSERT_FALSE(3ull < value);
    ASSERT_FALSE(1ull < value);
    ASSERT_FALSE(0ull < value);
}

TYPED_TEST(Integer, int_compare7) {
    typename TypeParam::type value = 7;

    ASSERT_TRUE(7ull <= value);
    ASSERT_TRUE(6ull <= value);
    ASSERT_TRUE(3ull <= value);
    ASSERT_TRUE(1ull <= value);
    ASSERT_TRUE(0ull <= value);
}

TYPED_TEST(Integer, Minus) {
    typename TypeParam::type value = 0;
    value -= 1llu;
    ASSERT_TRUE(value == (typename TypeParam::type(0) - 1llu));
    ASSERT_TRUE(typename TypeParam::type(0) == ~value);
}

TYPED_TEST(Integer, Plus) {
    typename TypeParam::type value = 0;
    value += sdsl::bits::lo_set[64];
    ASSERT_TRUE(typename TypeParam::type(0) + sdsl::bits::lo_set[64] == value);
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[64]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[64] == value);
}

TYPED_TEST(Integer, Plus1) {
    typename TypeParam::type value = sdsl::bits::lo_set[31];
    value <<= 32;
    value += sdsl::bits::lo_set[32];
    value += 1;
    value >>= 63;
    ASSERT_TRUE(1llu == value);
}

TYPED_TEST(Integer, And1) {
    typename TypeParam::type value = sdsl::bits::lo_set[32];
    ASSERT_TRUE(sdsl::bits::lo_set[32] == (value & sdsl::bits::lo_set[32]));
    value &= sdsl::bits::lo_set[32];
    ASSERT_TRUE(sdsl::bits::lo_set[32] == value);
    value <<= 32;
    ASSERT_TRUE(sdsl::bits::lo_unset[32] == (value & sdsl::bits::lo_set[64]));
    value &= sdsl::bits::lo_set[64];
    ASSERT_TRUE(sdsl::bits::lo_unset[32] == value);

    ASSERT_TRUE(0llu == ((value << (sizeof(typename TypeParam::type) * 8 - 64)) & 0));
    ASSERT_TRUE(0llu != value);
    value <<= sizeof(typename TypeParam::type) * 8 - 64;
    ASSERT_TRUE(0llu != value);
    value &= 0llu;
    ASSERT_TRUE(0llu == value);
}

TYPED_TEST(Integer, LeftShift) {
    typename TypeParam::type value = sdsl::bits::lo_set[32];
    value <<= 32;
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[32]) << 32 == value);
    ASSERT_TRUE(typename TypeParam::type(uint64_t(sdsl::bits::lo_set[32]) << 32) == value);
    ASSERT_TRUE((uint64_t(sdsl::bits::lo_set[32]) << 32) == value);

    // Avoid warning "shift count >= width of type"
    ASSERT_TRUE(((value << (sizeof(typename TypeParam::type) * 8 - 1)) << 1) == 0llu);
    // Avoid warning "shift count >= width of type"
    value <<= sizeof(typename TypeParam::type) * 8 - 1;
    value <<= 1;
    ASSERT_TRUE(0llu == value);

    value = sdsl::bits::lo_unset[32];
    ASSERT_TRUE(0llu != (value << (sizeof(typename TypeParam::type) * 8 - 64)));
    value <<= sizeof(typename TypeParam::type) * 8 - 64;
    ASSERT_TRUE(0llu != value);
}

TYPED_TEST(Integer, RightShift) {
    typename TypeParam::type value = uint64_t(sdsl::bits::lo_set[32]) << 32;
    value >>= 32;
    ASSERT_TRUE((typename TypeParam::type(uint64_t(sdsl::bits::lo_set[32]) << 32) >> 32) == value);
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[32]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[32] == value);
    value >>= 32;
    ASSERT_TRUE(0llu == value);
}

TYPED_TEST(Integer, LeftRightShift) {
    typename TypeParam::type value = sdsl::bits::lo_set[32];

    value <<= sizeof(typename TypeParam::type) * 8 - 32;
    value >>= sizeof(typename TypeParam::type) * 8 - 32;
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[32]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[32] == value);

    value <<= sizeof(typename TypeParam::type) * 8 - 31;
    value >>= sizeof(typename TypeParam::type) * 8 - 31;
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[31]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[31] == value);

    value <<= sizeof(typename TypeParam::type) * 8 - 28;
    value >>= sizeof(typename TypeParam::type) * 8 - 28;
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[28]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[28] == value);

    value <<= sizeof(typename TypeParam::type) * 8 - 3;
    value >>= sizeof(typename TypeParam::type) * 8 - 3;
    ASSERT_TRUE(typename TypeParam::type(sdsl::bits::lo_set[3]) == value);
    ASSERT_TRUE(sdsl::bits::lo_set[3] == value);
}

TYPED_TEST(Kmer, Simple) {
    KMer<typename TypeParam::type, kBitsPerChar> kmer(
        std::vector<uint64_t>(sizeof(typename TypeParam::type) * 8 / kBitsPerChar, 0)
    );
    for (size_t i = 0; i < sizeof(typename TypeParam::type) * 8 / kBitsPerChar; ++i) {
        ASSERT_EQ(0llu, kmer[i]);
    }
}

TYPED_TEST(Kmer, Invertible) {
    test_kmer_packed_codec<typename TypeParam::type>("ATGG", "ATGG");
}

TYPED_TEST(Kmer, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < sizeof(typename TypeParam::type) * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, sizeof(typename TypeParam::type) * 8 / kBitsPerChar);
    //test bit shifting
    KMer<typename TypeParam::type, kBitsPerChar> kmer_builtup(0u);
    size_t k = 0;
    //ASSERT_EQ(k, kmer_builtup.get_k());
    ASSERT_EQ(k * kBitsPerChar, sdsl::bits::hi(kmer_builtup.seq_));
    for (int i = long_seq.length() - 1; i >= 0; --i) {
        kmer_builtup.seq_ <<= static_cast<uint64_t>(kBitsPerChar);
        kmer_builtup.seq_ |= kmer_extractor.encode(long_seq[i]);
        ++k;
    }
    std::string dec = kmer_builtup.to_string(long_seq.length(),
                                             kmer_extractor.alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_packed_codec<typename TypeParam::type>(long_seq, long_seq);
}

TYPED_TEST(Kmer, UpdateKmer) {
    KMer<typename TypeParam::type, kBitsPerChar> kmer[2] = {
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode("ATGC")),
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode("TGCT"))
    };
    KMer<typename TypeParam::type, kBitsPerChar> updated = kmer[0];
    updated.to_next(4, kmer_extractor.encode('T'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, kmer_extractor.encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TYPED_TEST(Kmer, NextPrevKmer) {
    KMer<typename TypeParam::type, kBitsPerChar> kmer[2] = {
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode("ATGC")),
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, kmer_extractor.encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, kmer_extractor.encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TYPED_TEST(Kmer, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < sizeof(typename TypeParam::type) * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, sizeof(typename TypeParam::type) * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMer<typename TypeParam::type, kBitsPerChar> kmer[2] = {
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode(long_seq)),
        KMer<typename TypeParam::type, kBitsPerChar>(kmer_extractor.encode(long_seq_alt))
    };

    kmer[0].to_next(long_seq.length(), kmer_extractor.encode('T'));

    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), kmer_extractor.alphabet),
              kmer[0].to_string(long_seq.length(), kmer_extractor.alphabet));
}

TYPED_TEST(Kmer, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA";
    long_seq0.resize(std::min(sizeof(typename TypeParam::type) * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAT";
    long_seq1.resize(std::min(sizeof(typename TypeParam::type) * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = kmer_extractor.encode(long_seq0);
    KMer<typename TypeParam::type, kBitsPerChar> kmer0(seq0.begin(), seq0.size());

    kmer0.to_next(long_seq0.length(), kmer_extractor.encode(long_seq1.back()));

    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               kmer_extractor.alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(kmer_extractor.encode(long_seq1.back()));
    KMer<typename TypeParam::type, kBitsPerChar> kmer1(seq0.begin() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               kmer_extractor.alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TYPED_TEST(Kmer, InvertibleEndDol) {
#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("ATG$", "ATGN");
#elif _PROTEIN_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("ATG$", "ATGX");
#elif _DNA_GRAPH
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("ATG$"), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

TYPED_TEST(Kmer, InvertibleStartDol) {
#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("$ATGG", "NATGG");
#elif _PROTEIN_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("$ATGG", "XATGG");
#elif _DNA_GRAPH
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("$ATGG"), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

TYPED_TEST(Kmer, InvertibleBothDol) {
#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("$ATG$", "NATGN");
#elif _PROTEIN_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("$ATG$", "XATGX");
#elif _DNA_GRAPH
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("$ATG$"), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

TYPED_TEST(Kmer, InvalidChars) {
    KMer<typename TypeParam::type, kBitsPerChar> kmer(kmer_extractor.encode("ATGC"));

#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("ATGH", "ATGN");
    test_kmer_packed_codec<typename TypeParam::type>("ATGЯ", "ATGNN"); // cyrillic
    test_kmer_packed_codec<typename TypeParam::type>("ATGН", "ATGNN"); // cyrillic
#elif _PROTEIN_GRAPH
    test_kmer_packed_codec<typename TypeParam::type>("ATGH", "ATGH");
    test_kmer_packed_codec<typename TypeParam::type>("ATGЯ", "ATGXX"); // cyrillic
    test_kmer_packed_codec<typename TypeParam::type>("ATGН", "ATGXX"); // cyrillic
#elif _DNA_GRAPH
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("ATGH"), "");
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("ATGЯ"), ""); // cyrillic
    ASSERT_DEBUG_DEATH(kmer_packed_codec<typename TypeParam::type>("ATGН"), ""); // cyrillic
    ASSERT_DEBUG_DEATH(kmer.to_next(4, kmer_extractor.encode('N')), "");
    ASSERT_DEBUG_DEATH(kmer.to_next(4, kmer_extractor.encode("Я")[0]), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

template <typename IntType>
void test_kmer_packed_less(const std::string &k1,
                           const std::string &k2, bool truth) {
    KMer<IntType, kBitsPerChar> kmer[2] = {
        KMer<IntType, kBitsPerChar>(kmer_extractor.encode(k1)),
        KMer<IntType, kBitsPerChar>(kmer_extractor.encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TYPED_TEST(Kmer, LessEdge) {
    test_kmer_packed_less<typename TypeParam::type>("ATGC", "ATGG", true);
}

TYPED_TEST(Kmer, Less) {
    test_kmer_packed_less<typename TypeParam::type>("ACTG", "GCTG", true);
}

TYPED_TEST(Kmer, LessLong) {
    test_kmer_packed_less<typename TypeParam::type>(
        std::string(sizeof(typename TypeParam::type) * 8 / kBitsPerChar - 1, 'A') +  "C",
        std::string(sizeof(typename TypeParam::type) * 8 / kBitsPerChar - 1, 'A') +  "T",
        true
    );

    test_kmer_packed_less<typename TypeParam::type>(
        std::string(sizeof(typename TypeParam::type) * 8 / kBitsPerChar - 2, 'A') + "CA",
        std::string(sizeof(typename TypeParam::type) * 8 / kBitsPerChar - 2, 'A') + "TA",
        true
    );
}

TEST(Kmer, TestPrint64) {
    size_t size = sizeof(uint64_t) * 8 / kBitsPerChar;
    KMer<uint64_t, kBitsPerChar> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _PROTEIN_GRAPH
    EXPECT_EQ("0000000000000000000000000000000000000000000000000084210842108421", out);
#else
    EXPECT_EQ("0000000000000000000000000000000000000000000000005555555555555555", out);
#endif
}

TEST(Kmer, TestPrint128) {
    size_t size = sizeof(sdsl::uint128_t) * 8 / kBitsPerChar;
    KMer<sdsl::uint128_t, kBitsPerChar> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _PROTEIN_GRAPH
    EXPECT_EQ("0000000000000000000000000000000001084210842108421084210842108421", out);
#else
    EXPECT_EQ("0000000000000000000000000000000055555555555555555555555555555555", out);
#endif
}

TEST(Kmer, TestPrint256) {
    size_t size = sizeof(sdsl::uint256_t) * 8 / kBitsPerChar;
    KMer<sdsl::uint256_t, kBitsPerChar> kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _PROTEIN_GRAPH
    EXPECT_EQ("0421084210842108421084210842108421084210842108421084210842108421", out);
#else
    EXPECT_EQ("5555555555555555555555555555555555555555555555555555555555555555", out);
#endif
}


template <typename T>
T encode_c(char c, const T *char_map) {
    assert(static_cast<size_t>(c) < 128);
    return char_map[static_cast<size_t>(c)];
}

template <typename T>
std::vector<T> encode_c(const std::string &sequence, const T *char_map) {
    std::vector<T> encoded;
    std::transform(sequence.begin(), sequence.end(),
                   std::back_inserter(encoded),
                   [&](char c) { return encode_c(c, char_map); });
    assert(encoded.size() == sequence.size());

    return encoded;
}


// Nucleotide 2 bit
std::vector<TAlphabet> encode(const std::string &sequence) {
    const TAlphabet kCharToNucleotide[128] = {
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 1,  0, 0, 0, 2,  0, 0, 0, 0,  0, 0, 0, 0,
        0, 0, 0, 0,  3, 3, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0
    };

    return encode_c(sequence, kCharToNucleotide);
}

template <typename G, int L>
void test_kmer_codec(const std::string &sequence) {
    const auto encoded = encode(sequence);

    for (uint64_t k = 2; k < sizeof(G) * 8 / L; ++k) {
        if (k > encoded.size())
            break;

        ASSERT_LE(k, encoded.size());
        std::vector<KMer<G, L>> kmers;
        KMer<G, L> kmer_packed(encoded.data(), k);

        for (uint64_t i = 0; i + k <= encoded.size(); ++i) {
            kmers.emplace_back(std::vector<TAlphabet>(encoded.begin() + i,
                                                      encoded.begin() + i + k));
            assert(i + k <= sequence.size());
            ASSERT_EQ(sequence.substr(i, k), kmers.back().to_string(k, "ACGT"))
                 << sequence.substr(i, k) << " " << k << " " << i
                 << "\n" << kmers.back()
                 << "\n" << kmers.back()[0]
                    << " " << kmers.back()[1]
                    << " " << kmers.back()[2]
                    << " " << kmers.back()[3]
                 << "\n" << kmers.back().to_string(k, "ACGT");

            KMer<G, L> kmer_alt(encoded.data() + i, k);
            ASSERT_EQ(kmers.back(), kmer_alt) << sequence.substr(i, k) << " " << k << " " << i;

            ASSERT_EQ(kmers.back(), kmer_packed) << sequence.substr(i, k) << " " << k << " " << i;

            if (i + k < encoded.size())
                kmer_packed.to_next(k, encoded[i + k]);
        }
    }
}

TYPED_TEST(Kmer, nucleotide_alphabet_pack_6_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<typename TypeParam::type, 2>(sequence);
}

TYPED_TEST(Kmer, nucleotide_alphabet_pack_6) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<typename TypeParam::type, 3>(sequence);
}

} // namespace
