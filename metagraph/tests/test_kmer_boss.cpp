#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#define private public
#define protected public

#include "kmer/kmer_boss.hpp"
#include "kmer/kmer_extractor.hpp"


namespace {

using namespace mtg::kmer;

template <typename T>
T encode_c(char c, const T *char_map) {
    assert(static_cast<size_t>(c) < 128);
    return char_map[static_cast<size_t>(c)];
}

template <typename T>
char decode_c(T a, const std::string &alphabet) {
    assert(a < alphabet.size());
    return alphabet[a];
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

template <typename T>
std::string decode_c(const std::vector<T> &encoded, const std::string &alphabet) {
    std::string decoded;
    std::transform(encoded.begin(), encoded.end(),
                   std::back_inserter(decoded),
                   [&](T a) { return decode_c(a, alphabet); });
    assert(decoded.size() == encoded.size());

    return decoded;
}

typedef uint8_t TAlphabet;
template std::vector<TAlphabet> encode_c(const std::string &sequence, const TAlphabet *char_map);


 // Nucleotide
std::vector<TAlphabet> encode_nucleotide(const std::string &sequence) {
    const TAlphabet kCharToNucleotide[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    return encode_c(sequence, kCharToNucleotide);
}

template <typename G, int L>
std::string decode_nucleotide(const KMerBOSS<G, L> &kmer, size_t k) {
    return kmer.to_string(k, "ACGTN");
}

// Workaround to avoid the issues with the uint128_t type traits undefined in AppleClang
struct UINT64_ { typedef uint64_t type; };
struct UINT128_ { typedef sdsl::uint128_t type; };
struct UINT256_ { typedef sdsl::uint256_t type; };
typedef ::testing::Types<UINT64_,
                         UINT128_,
                         UINT256_> IntTypes;
template <class KMER>
class KmerBOSS : public ::testing::Test { };
TYPED_TEST_SUITE(KmerBOSS, IntTypes);

TYPED_TEST(KmerBOSS, nucleotide_alphabet_pack) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    const auto encoded = encode_nucleotide(sequence);

    for (uint64_t k = 2; k < sizeof(typename TypeParam::type) * 8 / 3; ++k) {
        if (k > encoded.size())
            break;

        ASSERT_LE(k, encoded.size());
        std::vector<KMerBOSS<typename TypeParam::type, 3>> kmers;
        KMerBOSS<typename TypeParam::type, 3> kmer_packed(encoded.data(), k);
        for (uint64_t i = 0; i + k <= encoded.size(); ++i) {
            kmers.emplace_back(std::vector<TAlphabet>(encoded.begin() + i,
                                                      encoded.begin() + i + k));
            auto decoded = decode_nucleotide<typename TypeParam::type, 3>(kmers.back(), k);
            EXPECT_EQ(sequence.substr(i, k), decoded) << k << " " << i;

            KMerBOSS<typename TypeParam::type, 3> kmer_alt(encoded.data() + i, k);
            EXPECT_EQ(kmers.back(), kmer_alt) << k << " " << i;

            EXPECT_EQ(kmers.back(), kmer_packed) << k << " " << i;

            if (i + k < encoded.size())
                kmer_packed.to_next(k, encoded[i + k], encoded[i + k - 1]);
        }
    }
}

//typedef uint64_t KMerBaseType;
const size_t kBitsPerChar = KmerExtractorBOSS::bits_per_char;

template <typename T>
using KMER = KMerBOSS<T, kBitsPerChar>;

#define kSizeOfKmer ( sizeof(typename TypeParam::type) )
//typedef KMerBOSS<KMerBaseType, kBitsPerChar> KMER;
//const size_t kSizeOfKmer = sizeof(KMerBaseType);

template <typename KMER>
std::string kmer_codec(const std::string &test_kmer) {
    std::vector<uint64_t> kmer(test_kmer.size());
    std::transform(test_kmer.begin(), test_kmer.end(), kmer.begin(),
        [](char c) {
            return c == KmerExtractorBOSS::alphabet[0]
                        ? 0
                        : KmerExtractorBOSS::encode(c);
        }
    );
    return KMER(kmer).to_string(test_kmer.length(), KmerExtractorBOSS::alphabet);
}

template <typename T>
void test_kmer_codec(const std::string &test_kmer,
                     const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_codec<KMER<T>>(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_codec<KMER<T>>(test_kmer));
}

TYPED_TEST(KmerBOSS, Invertible) {
    test_kmer_codec<typename TypeParam::type>("ATGG", "ATGG");
}

TYPED_TEST(KmerBOSS, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    //test bit shifting
    KMER<typename TypeParam::type> kmer_builtup(KmerExtractorBOSS::encode(
        std::string(long_seq.rbegin() + 1,
                    long_seq.rbegin() + 3)
    ));
    for (int i = long_seq.length() - 4; i >= 0; --i) {
        kmer_builtup.seq_ <<= static_cast<uint64_t>(kBitsPerChar);
        kmer_builtup.seq_ |= KmerExtractorBOSS::encode(long_seq[i]);
    }
    kmer_builtup.seq_ <<= static_cast<uint64_t>(kBitsPerChar);
    kmer_builtup.seq_ |= KmerExtractorBOSS::encode(long_seq[long_seq.length() - 1]);
    std::string dec = kmer_builtup.to_string(long_seq.length(), KmerExtractorBOSS::alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_codec<typename TypeParam::type>(long_seq, long_seq);
}

TYPED_TEST(KmerBOSS, UpdateKmer) {
    KMER<typename TypeParam::type> kmer[2] = {
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode("ATGC")),
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode("TGCT"))
    };
    KMER<typename TypeParam::type> updated = kmer[0];
    updated.to_next(4, KmerExtractorBOSS::encode('T'), KmerExtractorBOSS::encode('C'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractorBOSS::encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TYPED_TEST(KmerBOSS, NextPrevKmer) {
    KMER<typename TypeParam::type> kmer[2] = {
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode("ATGC")),
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractorBOSS::encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, KmerExtractorBOSS::encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TYPED_TEST(KmerBOSS, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER<typename TypeParam::type> kmer[2] = {
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode(long_seq)),
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode(long_seq_alt))
    };
    kmer[0].to_next(long_seq.length(), KmerExtractorBOSS::encode('T'),
                                       KmerExtractorBOSS::encode(long_seq.back()));
    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), KmerExtractorBOSS::alphabet),
              kmer[0].to_string(long_seq.length(), KmerExtractorBOSS::alphabet));
}

TYPED_TEST(KmerBOSS, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = KmerExtractorBOSS::encode(long_seq0);
    KMER<typename TypeParam::type> kmer0(seq0.data(), seq0.size());
    kmer0.to_next(long_seq0.length(), KmerExtractorBOSS::encode(long_seq1.back()),
                                      KmerExtractorBOSS::encode(long_seq0.back()));
    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               KmerExtractorBOSS::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(KmerExtractorBOSS::encode(long_seq1.back()));
    KMER<typename TypeParam::type> kmer1(seq0.data() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               KmerExtractorBOSS::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TYPED_TEST(KmerBOSS, InvertibleEndDol) {
    test_kmer_codec<typename TypeParam::type>("ATG$", "ATG$");
}

TYPED_TEST(KmerBOSS, InvertibleStartDol) {
    test_kmer_codec<typename TypeParam::type>("$ATGG", "$ATGG");
}

TYPED_TEST(KmerBOSS, InvertibleBothDol) {
    test_kmer_codec<typename TypeParam::type>("$ATG$", "$ATG$");
}

TYPED_TEST(KmerBOSS, InvalidChars) {
#if _DNA5_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_codec<typename TypeParam::type>("ATGH", "ATGN");
    test_kmer_codec<typename TypeParam::type>("ATGЯ", "ATGNN");
#elif _PROTEIN_GRAPH
    test_kmer_codec<typename TypeParam::type>("ATGH", "ATGH");
    test_kmer_codec<typename TypeParam::type>("ATGЯ", "ATGXX");
#elif _DNA_GRAPH
    ASSERT_DEBUG_DEATH(kmer_codec<KMER<typename TypeParam::type>>("ATGH"), "");
    ASSERT_DEBUG_DEATH(kmer_codec<KMER<typename TypeParam::type>>("ATGЯ"), "");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

template <typename T>
void test_kmer_less(const std::string &k1,
                    const std::string &k2, bool truth) {
    KMER<T> kmer[2] = {
        KMER<T>(KmerExtractorBOSS::encode(k1)),
        KMER<T>(KmerExtractorBOSS::encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TYPED_TEST(KmerBOSS, LessEdge) {
    test_kmer_less<typename TypeParam::type>("ATGC", "ATGG", true);
}

TYPED_TEST(KmerBOSS, Less) {
    test_kmer_less<typename TypeParam::type>("ACTG", "GCTG", true);
}

TYPED_TEST(KmerBOSS, LessLong) {
    test_kmer_less<typename TypeParam::type>(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "C",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "T",
        true
    );

    test_kmer_less<typename TypeParam::type>(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "CA",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "TA",
        true
    );
}

template <typename T>
void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    KMER<T> kmer[2] = {
        KMER<T>(KmerExtractorBOSS::encode(k1)),
        KMER<T>(KmerExtractorBOSS::encode(k2))
    };
    ASSERT_EQ(truth, KMER<T>::compare_suffix(kmer[0], kmer[1], 1));
}

TYPED_TEST(KmerBOSS, CompareSuffixTrue) {
    test_kmer_suffix<typename TypeParam::type>("ACTG", "GCTG", true);
}

TYPED_TEST(KmerBOSS, CompareSuffixFalse) {
    test_kmer_suffix<typename TypeParam::type>("ATTG", "ACTG", false);
}

TYPED_TEST(KmerBOSS, CompareSuffixTrueLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[0] = 'T';
    KMER<typename TypeParam::type> kmer[2] = {
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode(long_seq)),
        KMER<typename TypeParam::type>(KmerExtractorBOSS::encode(long_seq_alt))
    };
    ASSERT_TRUE(KMER<typename TypeParam::type>::compare_suffix(kmer[0], kmer[1], 1));

    //shift, then compare
    long_seq_alt[kSizeOfKmer * 8 / kBitsPerChar - 2] = 'T';

    kmer[0].seq_
        = kmer[0].seq_ >> static_cast<int>((kSizeOfKmer * 8 / kBitsPerChar - 2)
                                                * kBitsPerChar);

    kmer[1] = KMER<typename TypeParam::type>(KmerExtractorBOSS::encode(
        long_seq_alt.substr(kSizeOfKmer * 8 / kBitsPerChar - 2)
    ));

    ASSERT_TRUE(KMER<typename TypeParam::type>::compare_suffix(kmer[0], kmer[1], 1));
}

TYPED_TEST(KmerBOSS, CompareSuffixFalseLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[1] = 'T';

    test_kmer_suffix<typename TypeParam::type>(long_seq, long_seq_alt, false);
}

TYPED_TEST(KmerBOSS, SizeOfClass) {
    EXPECT_EQ(kSizeOfKmer, sizeof(KMER<typename TypeParam::type>));
}


TEST(KmerBOSS, TestPrint64) {
    size_t size = sizeof(uint64_t) * 8 / kBitsPerChar;
    KMER<uint64_t> kmer(std::vector<uint64_t>(size, 1));
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("0000000000000000000000000000000000000000000000001249249249249249", out);
#endif
}

TEST(KmerBOSS, TestPrint128) {
    size_t size = sizeof(sdsl::uint128_t) * 8 / kBitsPerChar;
    KMER<sdsl::uint128_t> kmer(std::vector<uint64_t>(size, 1));
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("0000000000000000000000000000000009249249249249249249249249249249", out);
#endif
}

TEST(KmerBOSS, TestPrint256) {
    size_t size = sizeof(sdsl::uint256_t) * 8 / kBitsPerChar;
    KMER<sdsl::uint256_t> kmer(std::vector<uint64_t>(size, 1));
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA_GRAPH || _DNA5_GRAPH
    EXPECT_EQ("1249249249249249249249249249249249249249249249249249249249249249", out);
#endif
}

} // namespace
