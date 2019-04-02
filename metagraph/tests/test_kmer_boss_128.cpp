#include <stdio.h>

#include "gtest/gtest.h"

#define private public
#define protected public

#include "dbg_succinct.hpp"
#include "kmer_boss.hpp"
#include "kmer_extractor.hpp"
#include "utils.hpp"

typedef sdsl::uint128_t KMerBaseType;
const size_t kBitsPerChar = KmerExtractor::kLogSigma;
typedef KMerBOSS<KMerBaseType, kBitsPerChar> KMER;
const size_t kSizeOfKmer = sizeof(KMerBaseType);


template <typename KMER>
std::string kmer_codec(const std::string &test_kmer);

std::string kmer_codec_128(const std::string &test_kmer) {
    return kmer_codec<KMER>(test_kmer);
}

void test_kmer_codec_128(const std::string &test_kmer,
                         const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_codec_128(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_codec_128(test_kmer));
}

TEST(KmerEncodeTest_128, Invertible) {
    test_kmer_codec_128("ATGG", "ATGG");
}

TEST(KmerEncodeTest_128, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    //test bit shifting
    KMER kmer_builtup(KmerExtractor::encode(
        std::string(long_seq.rbegin() + 1,
                    long_seq.rbegin() + 3)
    ));
    for (int i = long_seq.length() - 4; i >= 0; --i) {
        kmer_builtup.seq_ = kmer_builtup.seq_.operator<<(kBitsPerChar);
        kmer_builtup.seq_ |= KmerExtractor::encode(long_seq[i]);
    }
    kmer_builtup.seq_ = kmer_builtup.seq_.operator<<(kBitsPerChar);
    kmer_builtup.seq_ |= KmerExtractor::encode(long_seq[long_seq.length() - 1]);
    std::string dec = kmer_builtup.to_string(long_seq.length(), KmerExtractor::alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_codec_128(long_seq, long_seq);
}

TEST(KmerEncodeTest_128, UpdateKmer) {
    KMER kmer[2] = {
        KMER(KmerExtractor::encode("ATGC")),
        KMER(KmerExtractor::encode("TGCT"))
    };
    KMER updated = kmer[0];
    updated.to_next(4, KmerExtractor::encode('T'), KmerExtractor::encode('C'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractor::encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TEST(KmerEncodeTest_128, NextPrevKmer) {
    KMER kmer[2] = {
        KMER(KmerExtractor::encode("ATGC")),
        KMER(KmerExtractor::encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractor::encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, KmerExtractor::encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TEST(KmerEncodeTest_128, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER kmer[2] = {
        KMER(KmerExtractor::encode(long_seq)),
        KMER(KmerExtractor::encode(long_seq_alt))
    };
    kmer[0].to_next(long_seq.length(), KmerExtractor::encode('T'),
                                       KmerExtractor::encode(long_seq.back()));
    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), KmerExtractor::alphabet),
              kmer[0].to_string(long_seq.length(), KmerExtractor::alphabet));
}

TEST(KmerEncodeTest_128, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAANT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = KmerExtractor::encode(long_seq0);
    KMER kmer0(seq0.begin(), seq0.size());
    kmer0.to_next(long_seq0.length(), KmerExtractor::encode(long_seq1.back()),
                                      KmerExtractor::encode(long_seq0.back()));
    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               KmerExtractor::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(KmerExtractor::encode(long_seq1.back()));
    KMER kmer1(seq0.begin() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               KmerExtractor::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TEST(KmerEncodeTest_128, InvertibleEndDol) {
    test_kmer_codec_128("ATG$", "ATG$");
}

TEST(KmerEncodeTest_128, InvertibleStartDol) {
    test_kmer_codec_128("$ATGG", "$ATGG");
}

TEST(KmerEncodeTest_128, InvertibleBothDol) {
    test_kmer_codec_128("$ATG$", "$ATG$");
}

TEST(KmerEncodeTest_128, InvalidChars) {
#if _DNA_GRAPH || _DNA_CASE_SENSITIVE_GRAPH
    test_kmer_codec_128("ATGH", "ATGN");
    test_kmer_codec_128("ATGЯ", "ATGNN");
#elif _PROTEIN_GRAPH
    test_kmer_codec_128("ATGH", "ATGH");
    test_kmer_codec_128("ATGЯ", "ATGXX");
#elif _DNA4_GRAPH
    test_kmer_codec_128("ATGH", "ATGA");
    test_kmer_codec_128("ATGЯ", "ATGAA");
#else
    static_assert(false,
        "Add a unit test for checking behavior with invalid characters"
    );
#endif
}

void test_kmer_less_128(const std::string &k1,
                        const std::string &k2, bool truth) {
    KMER kmer[2] = {
        KMER(KmerExtractor::encode(k1)),
        KMER(KmerExtractor::encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TEST(KmerEncodeTest_128, LessEdge) {
    test_kmer_less_128("ATGC", "ATGG", true);
}

TEST(KmerEncodeTest_128, Less) {
    test_kmer_less_128("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest_128, LessLong) {
    test_kmer_less_128(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "C",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "T",
        true
    );

    test_kmer_less_128(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "CA",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "TA",
        true
    );
}

void test_kmer_suffix_128(std::string k1, std::string k2, bool truth) {
    KMER kmer[2] = {
        KMER(KmerExtractor::encode(k1)),
        KMER(KmerExtractor::encode(k2))
    };
    ASSERT_EQ(truth, KMER::compare_suffix(kmer[0], kmer[1], 1));
}

TEST(KmerEncodeTest_128, CompareSuffixTrue) {
    test_kmer_suffix_128("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest_128, CompareSuffixFalse) {
    test_kmer_suffix_128("ATTG", "ACTG", false);
}

TEST(KmerEncodeTest_128, CompareSuffixTrueLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[0] = 'T';
    KMER kmer[2] = {
        KMER(KmerExtractor::encode(long_seq)),
        KMER(KmerExtractor::encode(long_seq_alt))
    };
    ASSERT_TRUE(KMER::compare_suffix(kmer[0], kmer[1], 1));

    //shift, then compare
    long_seq_alt[kSizeOfKmer * 8 / kBitsPerChar - 2] = 'T';

    kmer[0].seq_
        = kmer[0].seq_ >> static_cast<int>((kSizeOfKmer * 8 / kBitsPerChar - 2)
                                                * kBitsPerChar);

    kmer[1] = KMER(KmerExtractor::encode(
        long_seq_alt.substr(kSizeOfKmer * 8 / kBitsPerChar - 2)
    ));

    ASSERT_TRUE(KMER::compare_suffix(kmer[0], kmer[1], 1));
}

TEST(KmerEncodeTest_128, CompareSuffixFalseLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[1] = 'T';

    test_kmer_suffix_128(long_seq, long_seq_alt, false);
}

TEST(KmerEncodeTest_128, SizeOfClass) {
    EXPECT_EQ(kSizeOfKmer, sizeof(KMER));
}

TEST(KmerEncodeTest_128, TestPrint) {
    size_t size = sizeof(KMerBaseType) * 8 / kBitsPerChar;
    KMER kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
#if _DNA4_GRAPH
    EXPECT_EQ("0000000000000000000000000000000009249249249249249249249249249249", out);
#elif _DNA_GRAPH
    EXPECT_EQ("0000000000000000000000000000000009249249249249249249249249249249", out);
#endif
}
