#include <stdio.h>

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

TEST(KmerPackedEncodeTest_64, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
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
    updated.to_next(4, KmerExtractor2Bit::encode('T'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractor2Bit::encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TEST(KmerPackedEncodeTest_64, NextPrevKmer) {
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode("ATGC")),
        KMER(KmerExtractor2Bit::encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, KmerExtractor2Bit::encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, KmerExtractor2Bit::encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TEST(KmerPackedEncodeTest_64, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER kmer[2] = {
        KMER(KmerExtractor2Bit::encode(long_seq)),
        KMER(KmerExtractor2Bit::encode(long_seq_alt))
    };

    kmer[0].to_next(long_seq.length(), KmerExtractor2Bit::encode('T'));

    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), KmerExtractor2Bit::alphabet),
              kmer[0].to_string(long_seq.length(), KmerExtractor2Bit::alphabet));
}

TEST(KmerPackedEncodeTest_64, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAN";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAANT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = KmerExtractor2Bit::encode(long_seq0);
    KMER kmer0(seq0.begin(), seq0.size());

    kmer0.to_next(long_seq0.length(), KmerExtractor2Bit::encode(long_seq1.back()));

    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               KmerExtractor2Bit::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(KmerExtractor2Bit::encode(long_seq1.back()));
    KMER kmer1(seq0.begin() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               KmerExtractor2Bit::alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TEST(KmerPackedEncodeTest_64, InvertibleEndDol) {
    test_kmer_packed_codec_64("ATG$", "ATGA");
}

TEST(KmerPackedEncodeTest_64, InvertibleStartDol) {
    test_kmer_packed_codec_64("$ATGG", "AATGG");
}

TEST(KmerPackedEncodeTest_64, InvertibleBothDol) {
    test_kmer_packed_codec_64("$ATG$", "AATGA");
}

TEST(KmerPackedEncodeTest_64, InvalidChars) {
#ifndef _PROTEIN_GRAPH
    test_kmer_packed_codec_64("ATGH", "ATGA");
#endif
    test_kmer_packed_codec_64("ATGЯ", "ATGAA");
}

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
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "C",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 1, 'A') +  "T",
        true
    );

    test_kmer_packed_less_64(
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "CA",
        std::string(kSizeOfKmer * 8 / kBitsPerChar - 2, 'A') + "TA",
        true
    );
}

#if _DNA_GRAPH
TEST(KmerPackedEncodeTest_64, TestPrint) {
    size_t size = sizeof(KMerPackedBaseType) * 8 / kBitsPerChar;
    KMER kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
    EXPECT_EQ("0000000000000000000000000000000000000000000000005555555555555555", out);
}
#endif
