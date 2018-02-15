#include <stdio.h>

#include "gtest/gtest.h"

#define private public
#define protected public

#include "dbg_succinct.hpp"
#include "kmer.hpp"
#include <sdsl/uint128_t.hpp>


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

void test_kmer_codec(const std::string &test_kmer,
                     const std::string &test_compare_kmer) {
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
            ASSERT_EQ(kmer.to_string(DBG_succ::alphabet),
                      long_seq + std::string(i - 2, curchar));
        }
        while (shift--) {
            kmer >>= kBitsPerChar;
            ASSERT_EQ(kmer.to_string(DBG_succ::alphabet),
                      long_seq + std::string(shift - 2, curchar));
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
    KMer kmer_builtup(std::string(long_seq.rbegin() + 1,
                      long_seq.rbegin() + 3), DBG_succ::encode);
    for (int i = long_seq.length() - 4; i >= 0; --i) {
        kmer_builtup.seq_ = kmer_builtup.seq_ << kBitsPerChar;
        kmer_builtup.seq_ |= DBG_succ::encode(long_seq[i]) + 1;
    }
    kmer_builtup.seq_ = kmer_builtup.seq_ << kBitsPerChar;
    kmer_builtup.seq_ |= DBG_succ::encode(long_seq[long_seq.length() - 1]) + 1;
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
    KMer::update_kmer(3, DBG_succ::encode('T'),
                   DBG_succ::encode('C'), &kmer[0].seq_);
    EXPECT_EQ(kmer[1], kmer[0]);
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
    KMer::update_kmer(long_seq.length() - 1,
                DBG_succ::encode('T'),
                DBG_succ::encode(long_seq.back()),
                &kmer[0].seq_);
    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(DBG_succ::alphabet),
              kmer[0].to_string(DBG_succ::alphabet));
}


TEST(KmerEncodeTest, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA$";
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA$T";
    std::deque<TAlphabet> seq0(long_seq0.length());
    std::transform(long_seq0.begin(), long_seq0.end(), seq0.begin(), DBG_succ::encode);
    KMer kmer0(KMer::pack_kmer(seq0.begin(), seq0.size()));
    KMer::update_kmer(
            long_seq0.length() - 1,
            DBG_succ::encode('T'),
            DBG_succ::encode('$'),
            reinterpret_cast<sdsl::uint256_t*>(&kmer0));
    std::string reconst_seq1 = kmer0.to_string(DBG_succ::alphabet);
    reconst_seq1.push_back(reconst_seq1[0]);
    reconst_seq1 = reconst_seq1.substr(1);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(DBG_succ::encode('T'));
    KMer kmer1(KMer::pack_kmer(seq0.begin() + 1, seq0.size() - 1));
    std::string reconst_seq2 = kmer1.to_string(DBG_succ::alphabet);
    reconst_seq2.push_back(reconst_seq2[0]);
    reconst_seq2 = reconst_seq2.substr(1);
    EXPECT_EQ(long_seq1, reconst_seq2);
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
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TEST(KmerEncodeTest, LessEdge) {
    test_kmer_less("ATGC", "ATGG", true);
}

TEST(KmerEncodeTest, Less) {
    test_kmer_less("ACTG", "GCTG", true);
}

TEST(KmerEncodeTest, LessLong) {
    test_kmer_less(std::string(79, 'A') +  "C",
                   std::string(79, 'A') +  "T", true);

    test_kmer_less(std::string(79, 'A') + "CA",
                   std::string(79, 'A') + "TA", true);
}

void test_kmer_suffix(std::string k1, std::string k2, bool truth) {
    KMer kmer[2] = {
        KMer(k1, DBG_succ::encode),
        KMer(k2, DBG_succ::encode)
    };
    ASSERT_EQ(truth, KMer::compare_kmer_suffix(kmer[0], kmer[1], 1));
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
    ASSERT_TRUE(KMer::compare_kmer_suffix(kmer[0], kmer[1], 1));

    //shift, then compare
    long_seq_alt[22] = 'T';

    kmer[0].seq_ = kmer[0].seq_ >> 66;
    kmer[1] = KMer(long_seq_alt.substr(22), DBG_succ::encode);

    ASSERT_TRUE(KMer::compare_kmer_suffix(kmer[0], kmer[1], 1));
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
    EXPECT_EQ(32u, sizeof(KMer));
}
