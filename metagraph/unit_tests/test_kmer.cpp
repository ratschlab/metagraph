#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "construct.hpp"
#include "kmer.hpp"

std::string kmer_codec(std::string &test_kmer) {
    const char *kmer_s = kmer_boost::kmertos(kmer_boost::stokmer(test_kmer.c_str(), test_kmer.length(), construct::nt_lookup), DBG_succ::default_alphabet, DBG_succ::default_alph_size);
    std::string kmer = std::string(kmer_s+1) + kmer_s[0];
    free((char*)kmer_s);
    return kmer;
}

void test_kmer_codec(std::string test_kmer, std::string test_compare_kmer = "") {
    if (test_compare_kmer == "")
        test_compare_kmer = test_kmer;
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    std::string kmer = kmer_codec(test_kmer);
    ASSERT_EQ(test_compare_kmer.length(), kmer.length());
    EXPECT_EQ(test_compare_kmer, kmer);
}

TEST(KmerEncodeTest, Invertible) {
    test_kmer_codec("ATGG");
}

TEST(KmerEncodeTest, InvertibleEndDol) {
    test_kmer_codec("ATG$");
}

TEST(KmerEncodeTest, InvertibleStartDol) {
    test_kmer_codec("$ATGG");
}

TEST(KmerEncodeTest, InvertibleBothDol) {
    test_kmer_codec("$ATG$");
}

TEST(KmerEncodeTest, InvalidChars) {
    test_kmer_codec("ATGH", "ATGN");
}

void test_kmer_less(std::string k1, std::string k2, bool truth) {
    ui256 kmer[2] = {\
        kmer_boost::stokmer(k1.c_str(), k1.length(), construct::nt_lookup),
        kmer_boost::stokmer(k2.c_str(), k2.length(), construct::nt_lookup)
    };
    ASSERT_EQ(kmer[0] < kmer[1], truth);
}

TEST(KmerEncodeTest, LessEdge) {
    test_kmer_less("ATGC","ATGG", true);
}

TEST(KmerEncodeTest, Less) {
    test_kmer_less("ACTG", "GCTG", true);
}

void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    ui256 kmer[2] = {\
        kmer_boost::stokmer(k1.c_str(), k1.length(), construct::nt_lookup),
        kmer_boost::stokmer(k2.c_str(), k2.length(), construct::nt_lookup)
    };
    ASSERT_EQ(kmer_boost::compare_kmer_suffix(kmer[0], kmer[1]), truth);
}

TEST(KmerEncodeTest, CompareSuffixTrue) {
    test_kmer_less("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest, CompareSuffixFalse) {
    test_kmer_less("ATTG", "ACTG", false);
}

