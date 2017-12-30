#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "kmer.hpp"


std::string kmer_codec(const std::string &test_kmer) {
    std::string kmer_s = KMer(
        test_kmer,
        DBG_succ::encode
    ).to_string(
        DBG_succ::alphabet
    );
    kmer_s.push_back(kmer_s[0]);
    kmer_s.erase(kmer_s.begin());
    return kmer_s;
}

void test_kmer_codec(const std::string &test_kmer, const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_codec(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_codec(test_kmer));
}

TEST(KmerEncodeTest, Invertible) {
    test_kmer_codec("ATGG", "ATGG");
}

TEST(KmerEncodeTest, InvertibleEndDol) {
    test_kmer_codec("ATG$", "ATG$");
}

TEST(KmerEncodeTest, InvertibleStartDol) {
    test_kmer_codec("$ATGG", "$ATGG");
}

TEST(KmerEncodeTest, InvertibleBothDol) {
    test_kmer_codec("$ATG$", "$ATG$");
}

TEST(KmerEncodeTest, InvalidChars) {
    test_kmer_codec("ATGH", "ATGN");
}

void test_kmer_less(std::string k1, std::string k2, bool truth) {
    KMer kmer[2] = {
        KMer(k1, DBG_succ::encode),
        KMer(k2, DBG_succ::encode)
    };
    ASSERT_EQ(kmer[0] < kmer[1], truth);
}

TEST(KmerEncodeTest, LessEdge) {
    test_kmer_less("ATGC", "ATGG", true);
}

TEST(KmerEncodeTest, Less) {
    test_kmer_less("ACTG", "GCTG", true);
}

void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    KMer kmer[2] = {
        KMer(k1, DBG_succ::encode),
        KMer(k2, DBG_succ::encode)
    };
    ASSERT_EQ(KMer::compare_kmer_suffix(kmer[0], kmer[1]), truth);
}

TEST(KmerEncodeTest, CompareSuffixTrue) {
    test_kmer_less("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest, CompareSuffixFalse) {
    test_kmer_less("ATTG", "ACTG", false);
}

TEST(KmerTest, SizeOfClass) {
    EXPECT_EQ(sizeof(ui256), sizeof(KMer));
}
