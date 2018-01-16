#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "kmer.hpp"

#include <utility>

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

/*

TEST(KmerEncodeTest, Operations) {
    for (uint8_t j = 1; j <= kMax; ++j) {
        char curchar = DBG_succ::decode(j - 1);
        std::string long_seq = std::string(2, curchar);
        KMer kmer(long_seq, DBG_succ::encode);
        int shift = 256 / kBitsPerChar;
        for (int i = 3; i <= shift; ++i) {
            kmer <<= kBitsPerChar;
            kmer |= j;
            ASSERT_EQ(kmer.to_string(DBG_succ::alphabet), long_seq + std::string(i - 2, curchar));
        }
        while (shift--) {
            kmer >>= kBitsPerChar;
            ASSERT_EQ(kmer.to_string(DBG_succ::alphabet), long_seq + std::string(shift - 2, curchar));
        }
    }
}
*/

TEST(KmerEncodeTest, BitShiftBuild) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < 256 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, 256 / kBitsPerChar);
    assert(long_seq.back() != long_seq.front());
    //test bit shifting
    KMer kmer_builtup(std::string(long_seq.rbegin() + 1, long_seq.rbegin() + 3), DBG_succ::encode);
    for (int i = long_seq.length() - 4; i >= 0; --i) {
        kmer_builtup <<= kBitsPerChar;
        kmer_builtup |= DBG_succ::encode(long_seq[i]) + 1;
    }
    kmer_builtup <<= kBitsPerChar;
    kmer_builtup |= DBG_succ::encode(long_seq[long_seq.length() - 1]) + 1;
    std::string dec = kmer_builtup.to_string(DBG_succ::alphabet);
    dec.push_back(dec.front());
    dec.erase(dec.begin());
    ASSERT_EQ(long_seq, dec);

    KMer kmer(long_seq, DBG_succ::encode);
    test_kmer_codec(long_seq, long_seq);
}

TEST(KmerEncodeTest, UpdateKmer) {
    KMer kmer[2] = {
        KMer(std::string("ATGC"), DBG_succ::encode),
        KMer(std::string("TGCT"), DBG_succ::encode)
    };
    kmer[0].update(3, 4);
    EXPECT_EQ(kmer[0], kmer[1]);
}

TEST(KmerEncodeTest, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < 256 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, 256 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMer kmer[2] = {
        KMer(long_seq, DBG_succ::encode),
        KMer(long_seq_alt, DBG_succ::encode)
    };
    kmer[0].update(long_seq.length() - 1, 4);
    EXPECT_EQ(kmer[0], kmer[1]);
    EXPECT_EQ(kmer[0].to_string(DBG_succ::alphabet), kmer[1].to_string(DBG_succ::alphabet));
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

TEST(KmerEncodeTest, LessLong) {
    test_kmer_less(std::string(79, 'A') +  "C", std::string(79, 'A') +  "T", true);
    test_kmer_less(std::string(79, 'A') + "CA", std::string(79, 'A') + "TA", true);
}

void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    KMer kmer[2] = {
        KMer(k1, DBG_succ::encode),
        KMer(k2, DBG_succ::encode)
    };
    ASSERT_EQ(KMer::compare_kmer_suffix(kmer[0], kmer[1], 1), truth);
}

TEST(KmerEncodeTest, CompareSuffixTrue) {
    test_kmer_suffix("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest, CompareSuffixFalse) {
    test_kmer_suffix("ATTG", "ACTG", false);
}

TEST(KmerEncodeTest, CompareSuffixTrueLong) {
    std::string long_seq(256 / kBitsPerChar, 'A');
    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';
    std::string long_seq_alt(long_seq);

    long_seq_alt[0] = 'T';
    KMer kmer[2] = {
        KMer(long_seq, DBG_succ::encode),
        KMer(long_seq_alt, DBG_succ::encode)
    };
    ASSERT_EQ(KMer::compare_kmer_suffix(kmer[0], kmer[1], 1), true);

    //shift, then compare
    long_seq_alt[22] = 'T';
    kmer[0] >>= 66;
    kmer[1] = KMer(long_seq_alt.substr(22), DBG_succ::encode);
    ASSERT_EQ(KMer::compare_kmer_suffix(kmer[0], kmer[1], 1), true);
}

TEST(KmerEncodeTest, CompareSuffixFalseLong) {
    std::string long_seq(256 / kBitsPerChar, 'A');
    *(long_seq.rbegin()) = 'T';
    *(++long_seq.rbegin()) = 'C';
    std::string long_seq_alt(long_seq);
    long_seq_alt[1] = 'T';
    test_kmer_suffix(long_seq, long_seq_alt, false);
}

TEST(KmerTest, SizeOfClass) {
    EXPECT_EQ(sizeof(uint64_t) * 4, sizeof(KMer));
}


