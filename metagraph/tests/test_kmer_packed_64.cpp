#include <stdio.h>
#include <sdsl/bits.hpp>

#include "gtest/gtest.h"

#define private public
#define protected public

#include "dbg_succinct.hpp"
#include "kmer_packed.hpp"
#include "kmer_extractor.hpp"
#include "utils.hpp"

typedef uint64_t KMerPackedBaseType;
const size_t kBitsPerChar = KmerExtractor2Bit::kLogSigma;
typedef KMerPacked<KMerPackedBaseType, kBitsPerChar> KMER;
const size_t kSizeOfKmer = sizeof(KMerPackedBaseType);


template <typename KMER>
std::string kmer_packed_codec(const std::string &test_kmer) {
    std::vector<uint64_t> kmer(test_kmer.size());
    std::transform(test_kmer.begin(), test_kmer.end(), kmer.begin(),
        [](char c) { return KmerExtractor2Bit::encode(c); }
    );
    return KMER(kmer).to_string(test_kmer.length(), KmerExtractor2Bit::alphabet);
}

template std::string kmer_packed_codec<KMerPacked<uint64_t, KmerExtractor2Bit::kLogSigma>>(const std::string &test_kmer);
template std::string kmer_packed_codec<KMerPacked<sdsl::uint256_t, KmerExtractor2Bit::kLogSigma>>(const std::string &test_kmer);
template std::string kmer_packed_codec<KMerPacked<sdsl::uint128_t, KmerExtractor2Bit::kLogSigma>>(const std::string &test_kmer);

std::string kmer_packed_codec_64(const std::string &test_kmer) {
    return kmer_packed_codec<KMER>(test_kmer);
}

void test_kmer_packed_codec_64(const std::string &test_kmer,
                              const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_packed_codec_64(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_packed_codec_64(test_kmer));
}

TEST(KmerPackedEncodeTest_64, Invertible) {
    test_kmer_packed_codec_64("ATGG", "ATGG");
}

/*
TEST(KmerPackedEncodeTest_64, Operations) {
    for (uint8_t j = 1; j <= kMax; ++j) {
        char curchar = KmerExtractor2Bit::decode(j - 1);
        std::string long_seq = std::string(2, curchar);
        KMER kmer(long_seq, KmerExtractor2Bit::encode);
        int shift = kSizeOfKmer * 8 / kBitsPerChar;
        for (int i = 3; i <= shift; ++i) {
            kmer <<= kBitsPerChar;
            kmer |= j;
            ASSERT_EQ(kmer.to_string(long_seq.length(), KmerExtractor2Bit::alphabet),
                      long_seq + std::string(i - 2, curchar));
        }
        while (shift--) {
            kmer >>= kBitsPerChar;
            ASSERT_EQ(kmer.to_string(long_seq.length(), KmerExtractor2Bit::alphabet),
                      long_seq + std::string(shift - 2, curchar));
        }
    }
}
*/

TEST(KmerPackedEncodeTest_64, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() + 1 < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar - 1);
    //test bit shifting
    KMER kmer_builtup(0u);
    size_t k = 0;
    //ASSERT_EQ(k, kmer_builtup.get_k());
    ASSERT_EQ(k * kBitsPerChar, sdsl::bits::hi(kmer_builtup.seq_));
    for (int i = long_seq.length() - 1; i >= 0; --i) {
        kmer_builtup.seq_ = kmer_builtup.seq_ << kBitsPerChar;
        kmer_builtup.seq_ |= KmerExtractor2Bit::encode(long_seq[i]);
        ++k;
    }
    std::string dec = kmer_builtup.to_string(long_seq.length(),
                                             KmerExtractor2Bit::alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_packed_codec_64(long_seq, long_seq);
}

TEST(KmerPackedEncodeTest_64, UpdateKmer) {
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode("ATGC")),
        KMER(KmerExtractor2Bit::encode("TGCT"))
    };

    KMER updated = kmer[0];
    KMER::update_kmer(4, KmerExtractor2Bit::encode('T'),
                         KmerExtractor2Bit::encode('C'), &updated.seq_);
    EXPECT_EQ(kmer[1], updated);
}

TEST(KmerPackedEncodeTest_64, NextPrevKmer) {
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode("ATGC")),
        KMER(KmerExtractor2Bit::encode("TGCT"))
    };

    EXPECT_EQ(kmer[0], kmer[1].prev_kmer(4, KmerExtractor2Bit::encode('A')));
    kmer[0].next_kmer(4, KmerExtractor2Bit::encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TEST(KmerPackedEncodeTest_64, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() + 1 < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar - 1);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode(long_seq)),
        KMER(KmerExtractor2Bit::encode(long_seq_alt))
    };
    KMER::update_kmer(long_seq.length(),
                KmerExtractor2Bit::encode('T'),
                KmerExtractor2Bit::encode(long_seq.back()),
                &kmer[0].seq_);
    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), KmerExtractor2Bit::alphabet),
              kmer[0].to_string(long_seq.length(), KmerExtractor2Bit::alphabet));
}

TEST(KmerPackedEncodeTest_64, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar - 1,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAANT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar - 1,
                              long_seq0.size()));
    auto seq0 = KmerExtractor2Bit::encode(long_seq0);
    KMER kmer0(KMER::pack_kmer(seq0.begin(), seq0.size()));
    KMER::update_kmer(
            long_seq0.length(),
            KmerExtractor2Bit::encode(long_seq1.back()),
            KmerExtractor2Bit::encode(long_seq0.back()),
            reinterpret_cast<KMerPackedBaseType*>(&kmer0));
    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               KmerExtractor2Bit::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(KmerExtractor2Bit::encode(long_seq1.back()));
    KMER kmer1(KMER::pack_kmer(seq0.begin() + 1, seq0.size() - 1));
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               KmerExtractor2Bit::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TEST(KmerPackedEncodeTest_64, InvertibleEndDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("ATG$", "ATGA"), "");
}

TEST(KmerPackedEncodeTest_64, InvertibleStartDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("$ATGG", "AATGG"), "");
}

TEST(KmerPackedEncodeTest_64, InvertibleBothDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("$ATG$", "AATGA"), "");
}

/*
#ifndef _PROTEIN_GRAPH
TEST(KmerPackedEncodeTest_64, InvalidChars) {
    test_kmer_packed_codec_64("ATGH", "ATGA");
}
#endif
*/

void test_kmer_packed_less_64(const std::string &k1,
                       const std::string &k2, bool truth) {
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode(k1)),
        KMER(KmerExtractor2Bit::encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TEST(KmerPackedEncodeTest_64, LessEdge) {
    test_kmer_packed_less_64("ATGC", "ATGG", true);
}

TEST(KmerPackedEncodeTest_64, Less) {
    test_kmer_packed_less_64("ACTG", "GCTG", true);
}

TEST(KmerPackedEncodeTest_64, LessLong) {
    test_kmer_packed_less_64(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') +  "C",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') +  "T",
        true
    );

    test_kmer_packed_less_64(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 3, 'A') + "CA",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 3, 'A') + "TA",
        true
    );
}

void test_kmer_packed_suffix_64(std::string k1, std::string k2, bool truth) {
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode(k1)),
        KMER(KmerExtractor2Bit::encode(k2))
    };
    ASSERT_EQ(truth, KMER::compare_suffix(kmer[0], kmer[1], 1));
}

TEST(KmerPackedEncodeTest_64, CompareSuffixTrue) {
    test_kmer_packed_suffix_64("ACTG", "GCTG", true);
}

TEST(KmerPackedEncodeTest_64, CompareSuffixFalse) {
    test_kmer_packed_suffix_64("ATTG", "ACTG", false);
}

TEST(KmerPackedEncodeTest_64, CompareSuffixTrueLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[0] = 'T';
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode(long_seq)),
        KMER(KmerExtractor2Bit::encode(long_seq_alt))
    };
    ASSERT_TRUE(KMER::compare_suffix(kmer[0], kmer[1], 1));

    //shift, then compare
    long_seq_alt[kSizeOfKmer * 8 / kBitsPerChar - 3] = 'T';

    kmer[0].seq_
        = kmer[0].seq_ >> static_cast<int>((kSizeOfKmer * 8 / kBitsPerChar - 3)
                                                * kBitsPerChar);

    kmer[1] = KMER(KmerExtractor2Bit::encode(
        long_seq_alt.substr(kSizeOfKmer * 8 / kBitsPerChar - 3)
    ));

    ASSERT_TRUE(KMER::compare_suffix(kmer[0], kmer[1], 1));
}

TEST(KmerPackedEncodeTest_64, CompareSuffixFalseLong) {
    std::string long_seq(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A');

    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';

    std::string long_seq_alt(long_seq);

    long_seq_alt[1] = 'T';

    test_kmer_packed_suffix_64(long_seq, long_seq_alt, false);
}

TEST(KmerPackedEncodeTest_64, SizeOfClass) {
    EXPECT_EQ(kSizeOfKmer, sizeof(KMER));
}

/*
#if _DNA_GRAPH
TEST(KmerPackedEncodeTest_64, TestPrint) {
    KMER kmer(std::vector<uint64_t>(sizeof(KMerPackedBaseType) * 8 / kBitsPerChar - 1, 1),
              sizeof(KMerPackedBaseType) * 8 / kBitsPerChar - 1);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
    EXPECT_EQ("0000000000000000000000000000000000000000000000002492492492492492", out);
}
#endif
*/
