#include "gtest/gtest.h"

#define private public
#define protected public

#include "kmer_extractor.hpp"
#include "utils.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(a, b) (void)0
#endif
#endif


TEST(KmerExtractor, encode_decode) {
    KmerExtractor encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    EXPECT_EQ('N', encoder.decode(encoder.encode('N')));
#ifndef _PROTEIN_GRAPH
    EXPECT_EQ('N', encoder.decode(encoder.encode('X')));
    EXPECT_EQ('N', encoder.decode(encoder.encode(-1)));
#else
    EXPECT_EQ('X', encoder.decode(encoder.encode('X')));
    EXPECT_EQ('X', encoder.decode(encoder.encode(-1)));
#endif
}

KmerExtractor::Kmer64 to_kmer(const KmerExtractor &encoder,
                              const std::string &kmer) {
    Vector<KmerExtractor::Kmer64> kmers;
    encoder.sequence_to_kmers(encoder.encode(kmer), kmer.size(), {}, &kmers);
    return kmers[0];
}

TEST(KmerExtractor, encode_decode_kmer) {
    KmerExtractor encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer))) << kmer;
}

TEST(KmerExtractor, encode_decode_string) {
    KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor::Kmer256> kmers;

        encoder.sequence_to_kmers(encoder.encode(sequence), k, {}, &kmers);
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0]);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i])[k - 1]);
        }

        EXPECT_EQ(sequence, reconstructed);
    }
}


TEST(KmerExtractor2Bit, encode_decode) {
    KmerExtractor2Bit encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    EXPECT_EQ('A', encoder.decode(encoder.encode('N')));
    EXPECT_EQ('A', encoder.decode(encoder.encode('X')));
}

TEST(KmerExtractor2Bit, encode_decode_kmer) {
    KmerExtractor2Bit encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(std::string("AAAAAAAAAAA"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(std::string("AAAAATAAAAA"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(std::string("AAAAAAAAGAT"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
}

TEST(KmerExtractor2Bit, encode_decode_string) {
    KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGTAN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto kmers = encoder.sequence_to_kmers(sequence, k);
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0]);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i])[k - 1]);
        }

        EXPECT_EQ(std::string("AAGGCAGCCTACCCCTCTGTAA"), reconstructed);
    }
}
