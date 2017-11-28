#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "kmer.hpp"


std::string kmer_codec(std::string &test_kmer) {
    const char *kmer_s = KMer::kmertos(
        KMer::stokmer(test_kmer, DBG_succ::get_alphabet_number),
        DBG_succ::default_alphabet, DBG_succ::default_alph_size
    );
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
    KMer kmer[2] = {
        KMer::stokmer(k1.c_str(), k1.length(), DBG_succ::get_alphabet_number),
        KMer::stokmer(k2.c_str(), k2.length(), DBG_succ::get_alphabet_number)
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
    KMer kmer[2] = {
        KMer::stokmer(k1.c_str(), k1.length(), DBG_succ::get_alphabet_number),
        KMer::stokmer(k2.c_str(), k2.length(), DBG_succ::get_alphabet_number)
    };
    ASSERT_EQ(KMer::compare_kmer_suffix(kmer[0], kmer[1]), truth);
}

TEST(KmerEncodeTest, CompareSuffixTrue) {
    test_kmer_less("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest, CompareSuffixFalse) {
    test_kmer_less("ATTG", "ACTG", false);
}

