#include "gtest/gtest.h"

#include <vector>
#include <functional>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#define private public
#define protected public

#include "kmer.hpp"
#include "kmer_extractor.hpp"


typedef uint8_t TAlphabet;

typedef uint64_t KMerBaseType;
const size_t kBitsPerChar = KmerExtractor2Bit::bits_per_char;
typedef KMer<KMerBaseType, kBitsPerChar> KMER;
const size_t kSizeOfKmer = sizeof(KMerBaseType);

const KmerExtractor2Bit kmer_extractor;


template <typename KMER>
std::string kmer_packed_codec(const std::string &test_kmer) {
    std::vector<uint64_t> kmer(test_kmer.size());
    std::transform(test_kmer.begin(), test_kmer.end(), kmer.begin(),
        [](char c) { return kmer_extractor.encode(c); }
    );
    return KMER(kmer).to_string(test_kmer.length(), kmer_extractor.alphabet);
}

std::string kmer_packed_codec_64(const std::string &test_kmer) {
    return kmer_packed_codec<KMER>(test_kmer);
}

void test_kmer_packed_codec_64(const std::string &test_kmer,
                              const std::string &test_compare_kmer) {
    ASSERT_EQ(test_kmer.length(), test_compare_kmer.length());
    ASSERT_EQ(test_compare_kmer.length(), kmer_packed_codec_64(test_kmer).length());
    EXPECT_EQ(test_compare_kmer, kmer_packed_codec_64(test_kmer));
}

TEST(Kmer, Invertible) {
    test_kmer_packed_codec_64("ATGG", "ATGG");
}

TEST(Kmer, BitShiftBuild) {
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
        kmer_builtup.seq_ |= kmer_extractor.encode(long_seq[i]);
        ++k;
    }
    std::string dec = kmer_builtup.to_string(long_seq.length(),
                                             kmer_extractor.alphabet);
    ASSERT_EQ(long_seq, dec);

    test_kmer_packed_codec_64(long_seq, long_seq);
}

TEST(Kmer, UpdateKmer) {
    KMER kmer[2] = {
        KMER(kmer_extractor.encode("ATGC")),
        KMER(kmer_extractor.encode("TGCT"))
    };
    KMER updated = kmer[0];
    updated.to_next(4, kmer_extractor.encode('T'));
    EXPECT_EQ(kmer[1], updated);
    auto prev = kmer[1];
    prev.to_prev(4, kmer_extractor.encode('A'));
    EXPECT_EQ(kmer[0], prev);
}

TEST(Kmer, NextPrevKmer) {
    KMER kmer[2] = {
        KMER(kmer_extractor.encode("ATGC")),
        KMER(kmer_extractor.encode("TGCT"))
    };

    auto prev = kmer[1];
    prev.to_prev(4, kmer_extractor.encode('A'));
    EXPECT_EQ(kmer[0], prev);
    kmer[0].to_next(4, kmer_extractor.encode('T'));
    EXPECT_EQ(kmer[1], kmer[0]);
}

TEST(Kmer, UpdateKmerLong) {
    std::string long_seq = "ATGCCTGA";
    while (long_seq.length() < kSizeOfKmer * 8 / kBitsPerChar) {
        long_seq += long_seq;
    }
    long_seq = long_seq.substr(0, kSizeOfKmer * 8 / kBitsPerChar);
    std::string long_seq_alt(long_seq.substr(1));
    long_seq_alt.push_back('T');
    KMER kmer[2] = {
        KMER(kmer_extractor.encode(long_seq)),
        KMER(kmer_extractor.encode(long_seq_alt))
    };

    kmer[0].to_next(long_seq.length(), kmer_extractor.encode('T'));

    EXPECT_EQ(kmer[1], kmer[0]);
    EXPECT_EQ(kmer[1].to_string(long_seq.length(), kmer_extractor.alphabet),
              kmer[0].to_string(long_seq.length(), kmer_extractor.alphabet));
}

TEST(Kmer, UpdateKmerVsConstruct) {
    std::string long_seq0 = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAA";
    long_seq0.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    std::string long_seq1 =  "AGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAT";
    long_seq1.resize(std::min(kSizeOfKmer * 8 / kBitsPerChar,
                              long_seq0.size()));
    auto seq0 = kmer_extractor.encode(long_seq0);
    KMER kmer0(seq0.begin(), seq0.size());

    kmer0.to_next(long_seq0.length(), kmer_extractor.encode(long_seq1.back()));

    std::string reconst_seq1 = kmer0.to_string(long_seq0.length(),
                                               kmer_extractor.alphabet);
    EXPECT_EQ(long_seq1, reconst_seq1);

    seq0.emplace_back(kmer_extractor.encode(long_seq1.back()));
    KMER kmer1(seq0.begin() + 1, seq0.size() - 1);
    std::string reconst_seq2 = kmer1.to_string(long_seq1.length(),
                                               kmer_extractor.alphabet);
    EXPECT_EQ(long_seq1, reconst_seq2);
}

TEST(Kmer, InvertibleEndDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("ATG$", "ATGA"), "");
}

TEST(Kmer, InvertibleStartDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("$ATGG", "AATGG"), "");
}

TEST(Kmer, InvertibleBothDol) {
    ASSERT_DEATH(test_kmer_packed_codec_64("$ATG$", "AATGA"), "");
}

TEST(Kmer, InvalidChars) {
    KMER kmer(kmer_extractor.encode("ATGC"));

    ASSERT_DEATH(test_kmer_packed_codec_64("ATGH", "ATGA"), "");
    ASSERT_DEATH(test_kmer_packed_codec_64("ATGЯ", "ATGAA"), "");

    ASSERT_DEATH(kmer.to_next(4, kmer_extractor.encode('N')), "");
    ASSERT_DEATH(kmer.to_next(4, kmer_extractor.encode("Я")[0]), "");
}

void test_kmer_packed_less_64(const std::string &k1,
                              const std::string &k2, bool truth) {
    KMER kmer[2] = {
        KMER(kmer_extractor.encode(k1)),
        KMER(kmer_extractor.encode(k2))
    };
    ASSERT_EQ(truth, kmer[0] < kmer[1]);
}

TEST(Kmer, LessEdge) {
    test_kmer_packed_less_64("ATGC", "ATGG", true);
}

TEST(Kmer, Less) {
    test_kmer_packed_less_64("ACTG", "GCTG", true);
}

TEST(Kmer, LessLong) {
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

TEST(Kmer, TestPrint) {
    size_t size = sizeof(KMerBaseType) * 8 / kBitsPerChar;
    KMER kmer(std::vector<uint64_t>(size, 1), size);
    std::stringstream ss;
    ss << kmer;
    std::string out;
    ss >> out;
    EXPECT_EQ("0000000000000000000000000000000000000000000000005555555555555555", out);
}


template <typename T>
std::vector<T> encode_c(const std::string &sequence, const T *char_map);

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
std::string decode(const KMer<G, L> &kmer, size_t k) {
    return kmer.to_string(k, "ACGT");
}


template <typename G, int L>
void test_kmer_codec(const std::string &sequence,
                     const std::function<std::vector<TAlphabet>(const std::string&)> &encoder,
                     const std::function<std::string(const KMer<G, L>&, size_t)> &decoder) {
    const auto encoded = encoder(sequence);

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
            ASSERT_EQ(sequence.substr(i, k), decoder(kmers.back(), k))
                 << sequence.substr(i, k) << " " << k << " " << i;

            KMer<G, L> kmer_alt(encoded.data() + i, k);
            ASSERT_EQ(kmers.back(), kmer_alt) << sequence.substr(i, k) << " " << k << " " << i;

            ASSERT_EQ(kmers.back(), kmer_packed) << sequence.substr(i, k) << " " << k << " " << i;

            if (i + k < encoded.size())
                kmer_packed.to_next(k, encoded[i + k]);
        }
    }
}

TEST(Kmer, nucleotide_alphabet_pack_6_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<uint64_t, 2>(sequence, encode, decode<uint64_t, 2>);
}

TEST(Kmer, nucleotide_alphabet_pack_128_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint128_t, 2>(sequence, encode, decode<sdsl::uint128_t, 2>);
}

TEST(Kmer, nucleotide_alphabet_pack_256_2Bit) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint256_t, 2>(sequence, encode, decode<sdsl::uint256_t, 2>);
}

TEST(Kmer, nucleotide_alphabet_pack_6) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<uint64_t, 3>(sequence, encode, decode<uint64_t, 3>);
}

TEST(Kmer, nucleotide_alphabet_pack_128) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint128_t, 3>(sequence, encode, decode<sdsl::uint128_t, 3>);
}

TEST(Kmer, nucleotide_alphabet_pack_256) {
    const std::string sequence = "AAGGCAGCCTACCCCTCTGTCTCCACCTTTGAGAAACACTCATCCTCAGGCCATGCAGTGGAAA";
    test_kmer_codec<sdsl::uint256_t, 3>(sequence, encode, decode<sdsl::uint256_t, 3>);
}
