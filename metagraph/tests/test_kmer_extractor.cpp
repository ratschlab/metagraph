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
    encoder.sequence_to_kmers(kmer, kmer.size(), {}, &kmers);
    return kmers.at(1);
}

TEST(KmerExtractor, encode_decode_kmer) {
    KmerExtractor encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
}

TEST(KmerExtractor, encode_decode_string) {
    KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);
        ASSERT_LT(2u, kmers.size());
        kmers.erase(kmers.begin());
        kmers.erase(kmers.end() - 1);

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
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
}

KmerExtractor2Bit::Kmer64 to_kmer(const KmerExtractor2Bit &encoder,
                              const std::string &kmer) {
    Vector<KmerExtractor2Bit::Kmer64> kmers;
    encoder.sequence_to_kmers(kmer, kmer.size(), {}, &kmers);
    return kmers.at(0);
}

TEST(KmerExtractor2Bit, encode_decode_kmer) {
    KmerExtractor2Bit encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(std::string("AAAAAAAAAAA"), encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(std::string("AAAAATAAAAA"), encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(std::string("AAAAAAAAGAT"), encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
}

TEST(KmerExtractor2Bit, encode_decode_string) {
    KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor2Bit::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
        }

        EXPECT_EQ(std::string("AAGGCAGCCTACCCCTCTGA"), reconstructed);
    }
}
