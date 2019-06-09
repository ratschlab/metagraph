#include "gtest/gtest.h"

#define private public
#define protected public

#include "kmer_extractor.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif
#endif

template <typename Kmer>
class ExtractKmers2Bit : public ::testing::Test { };

typedef ::testing::Types<KMer<uint64_t, KmerExtractor2Bit::bits_per_char>,
                         KMer<sdsl::uint128_t, KmerExtractor2Bit::bits_per_char>,
                         KMer<sdsl::uint256_t, KmerExtractor2Bit::bits_per_char>> KmerTypes;

TYPED_TEST_CASE(ExtractKmers2Bit, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractor2Bit::bits_per_char )

const KmerExtractor2Bit kmer_extractor;


TEST(KmerExtractor2Bit, encode_decode) {
    KmerExtractor2Bit encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    // #if _DNA4_GRAPH
        // N->A in 2Bit mode
        EXPECT_EQ('A', encoder.decode(encoder.encode('N')));
    // #else
    //     EXPECT_EQ('N', encoder.decode(encoder.encode('N')));
    // #endif
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
    std::string expected;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;

    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;

    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;

    kmer = "ANANANANANA";
    // #if _DNA4_GRAPH
        expected = std::string("AAAAAAAAAAA");
    // #else
    //     expected = std::string("ANANANANANA");
    // #endif
    EXPECT_EQ(expected,
              encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;

    kmer = "ANANATANANA";
    // #if _DNA4_GRAPH
        expected = std::string("AAAAATAAAAA");
    // #else
    //     expected = std::string("ANANATANANA");
    // #endif
    EXPECT_EQ(expected,
              encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;

    kmer = "ANANANANGNT";
    // #if _DNA4_GRAPH
        expected = std::string("AAAAAAAAGAT");
    // #else
    //     expected = std::string("ANANANANGNT");
    // #endif
    EXPECT_EQ(expected,
              encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
}

TEST(KmerExtractor2Bit, encode_decode_string) {
    KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACNCCCTCTG";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor2Bit::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);
        EXPECT_EQ(kmers, encoder.sequence_to_kmers<KmerExtractor2Bit::Kmer256>(sequence, k));
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
        }
        // #if _DNA4_GRAPH
            EXPECT_EQ(std::string("AAGGCAGCCTACACCCTCTG"), reconstructed);
        // #else
        //     EXPECT_EQ(std::string("AAGGCAGCCTACNCCCTCTG"), reconstructed);
        // #endif
    }
}

TEST(KmerExtractor2Bit, encode_decode_string_suffix) {
    KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTG";
    std::vector<bool> bits;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        for (size_t len = 1; len < std::min(k, uint64_t(5)); ++len) {
            bits.assign(sequence.size() + 1 - k, false);
            for (const auto &suffix : utils::generate_strings("ATGC", len)) {
                uint64_t it = k - len;
                Vector<KmerExtractor2Bit::Kmer256> kmers;
                encoder.sequence_to_kmers(sequence, k, encoder.encode(suffix), &kmers);
                for (const auto &kmer : kmers) {
                    auto jt = sequence.find(suffix, it);
                    auto kmer_str = encoder.kmer_to_sequence(kmer, k);
                    EXPECT_EQ(suffix, kmer_str.substr(kmer_str.size() - len, len));
                    ASSERT_NE(std::string::npos, jt);
                    ++jt;
                    EXPECT_EQ(
                        std::string(
                            sequence.begin() + jt + len - k - 1,
                            sequence.begin() + jt + len - 1
                        ), kmer_str
                    );
                    ASSERT_GT(bits.size(), jt + len - k - 1);
                    ASSERT_FALSE(bits[jt + len - k - 1]);
                    bits[jt + len - k - 1] = 1;
                    it = jt;
                }

                EXPECT_EQ(
                    kmers,
                    encoder.sequence_to_kmers<KmerExtractor2Bit::Kmer256>(
                        sequence, k, false, encoder.encode(suffix)
                    )
                );
            }
            EXPECT_EQ(bits.size(), std::accumulate(bits.begin(), bits.end(), 0u));
        }
    }
}

TEST(KmerExtractor2Bit, encode_decode_string_canonical_suffix) {
    KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTG";
    std::vector<bool> bits;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        for (size_t len = 1; len < std::min(k, uint64_t(5)); ++len) {
            bits.assign(sequence.size() + 1 - k, false);
            for (const auto &suffix : utils::generate_strings("ATGC", len)) {
                Vector<KmerExtractor2Bit::Kmer256> kmers;
                encoder.sequence_to_kmers(sequence, k, encoder.encode(suffix), &kmers, true);
                EXPECT_EQ(
                    kmers,
                    encoder.sequence_to_kmers<KmerExtractor2Bit::Kmer256>(
                        sequence, k, true, encoder.encode(suffix)
                    )
                );
                for (const auto &kmer : kmers) {
                    auto kmer_str = encoder.kmer_to_sequence(kmer, k);
                    EXPECT_EQ(suffix, kmer_str.substr(kmer_str.size() - len, len));
                    auto it = sequence.find(kmer_str);
                    uint64_t rev_count = 0;
                    if (it == std::string::npos) {
                        reverse_complement(kmer_str.begin(), kmer_str.end());
                        rev_count++;
                        it = sequence.find(kmer_str);
                    }
                    ASSERT_NE(std::string::npos, it);
                    ASSERT_GT(bits.size(), it);
                    while (bits[it]) {
                        it = sequence.find(kmer_str, it + 1);
                        if (it == std::string::npos) {
                            reverse_complement(kmer_str.begin(), kmer_str.end());
                            ++rev_count;
                            ASSERT_GT(2u, rev_count);
                            it = sequence.find(kmer_str);
                        }
                        ASSERT_NE(std::string::npos, it)
                            << k << " " << len << " " << suffix << " "
                            << encoder.kmer_to_sequence(kmer, k) << " " << kmer_str;
                        ASSERT_GT(bits.size(), it);
                    }
                    bits[it] = 1;
                }
            }
            EXPECT_EQ(bits.size(), std::accumulate(bits.begin(), bits.end(), 0u))
                << k << " " << len;
        }
    }
}

using TAlphabet = KmerExtractor2Bit::TAlphabet;

TYPED_TEST(ExtractKmers2Bit, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        // NNN -> $NNN$
        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, {}, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers2Bit, ExtractKmersFromStringWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }
    }

    suffix.assign({ kmer_extractor.encode('T') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, suffix, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers2Bit, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 3; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TT"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, { 0, kmer_extractor.encode('T') }, &result
            );
            ASSERT_EQ(0u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TT"), &result
            );

            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            kmer_extractor.sequence_to_kmers(
                std::string(length, 'T'), k, kmer_extractor.encode("TA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'T');
            sequence[k - 1] = 'A';

            kmer_extractor.sequence_to_kmers(
                sequence, k, kmer_extractor.encode("TA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers2Bit, ExtractKmersFromStringAppend) {
    Vector<TypeParam> result;

    kmer_extractor.sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u, result.size());

    kmer_extractor.sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u * 2, result.size());
}
