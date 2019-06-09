#include "gtest/gtest.h"

#define private public
#define protected public

#include "kmer_extractor.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(a, b) (void)0
#endif
#endif

template <typename Kmer>
class ExtractKmers : public ::testing::Test { };

typedef ::testing::Types<KMerBOSS<uint64_t, KmerExtractor::bits_per_char>,
                         KMerBOSS<sdsl::uint128_t, KmerExtractor::bits_per_char>,
                         KMerBOSS<sdsl::uint256_t, KmerExtractor::bits_per_char>> KmerTypes;

TYPED_TEST_CASE(ExtractKmers, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractor::bits_per_char )


TEST(ExtractKmers, encode_decode) {
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

TEST(ExtractKmers, encode_decode_kmer) {
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

TEST(ExtractKmers, encode_decode_string) {
    KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);
        EXPECT_EQ(kmers, encoder.sequence_to_kmers<KmerExtractor::Kmer256>(sequence, k));
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

TEST(ExtractKmers, encode_decode_string_suffix) {
    KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    std::vector<bool> bits(sequence.size(), false);
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto sequence_dummy = std::string(k - 1, '$') + sequence + "$";
        for (size_t len = 1; len < std::min(k, uint64_t(5)); ++len) {
            bits.assign(sequence.length() + 1, false);
            for (const auto &suffix : utils::generate_strings("ATGCN$", len)) {
                std::vector<uint8_t> suffix_encoded;
                std::transform(suffix.begin(), suffix.end(), std::back_inserter(suffix_encoded),
                    [&](char c) { return c == '$' ? 0 : encoder.encode(c); });
                Vector<KmerExtractor::Kmer256> kmers;
                encoder.sequence_to_kmers(sequence, k, suffix_encoded, &kmers);
                EXPECT_EQ(
                    kmers,
                    encoder.sequence_to_kmers<KmerExtractor::Kmer256>(
                        sequence, k, false, suffix_encoded
                    )
                );
                auto it = k - len - 1;
                for (const auto &kmer : kmers) {
                    auto kmer_str = encoder.kmer_to_sequence(kmer, k);
                    EXPECT_EQ(suffix, kmer_str.substr(kmer_str.size() - len - 1, len));
                    auto jt = sequence_dummy.find(suffix, it);
                    ASSERT_NE(std::string::npos, jt);
                    ++jt;
                    ASSERT_GE(jt + len + 1, k);
                    ASSERT_GE(sequence_dummy.length(), jt + len);
                    EXPECT_EQ(
                        std::string(
                            sequence_dummy.begin() + jt + len - k,
                            sequence_dummy.begin() + jt + len
                        ), kmer_str
                    );
                    ASSERT_GT(bits.size(), jt + len - k);
                    ASSERT_FALSE(bits[jt + len - k]);
                    bits[jt + len - k] = 1;
                    it = jt;
                }
            }
            EXPECT_EQ(bits.size(), std::accumulate(bits.begin(), bits.end(), 0u))
                << k << " " << len;
        }
    }
}

TEST(ExtractKmers, encode_decode_string_canonical_suffix) {
    KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    std::vector<bool> bits;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto sequence_dummy = std::string(k - 1, '$') + sequence + "$";
        for (size_t len = 1; len < std::min(k - 1, uint64_t(5)); ++len) {
            bits.assign(sequence.length() + 1, false);
            for (const auto &suffix : utils::generate_strings("ATGCN$", len)) {
                Vector<KmerExtractor::Kmer256> kmers;
                std::vector<uint8_t> encoded;
                std::transform(suffix.begin(), suffix.end(), std::back_inserter(encoded),
                    [&](auto c) { return c == '$' ? 0 : encoder.encode(c); });
                encoder.sequence_to_kmers(sequence, k, encoded, &kmers, true);
                EXPECT_EQ(
                    kmers,
                    encoder.sequence_to_kmers<KmerExtractor::Kmer256>(
                        sequence, k, true, encoded
                    )
                );
                for (const auto &kmer : kmers) {
                    auto kmer_str = encoder.kmer_to_sequence(kmer, k);
                    EXPECT_EQ(suffix, kmer_str.substr(kmer_str.size() - len - 1, len));
                    auto it = sequence_dummy.find(kmer_str);
                    uint64_t rev_count = 0;
                    if (it == std::string::npos) {
                        reverse_complement(kmer_str.begin(), kmer_str.end());
                        rev_count++;
                        it = sequence_dummy.find(kmer_str);
                    }
                    ASSERT_NE(std::string::npos, it);
                    ASSERT_GT(bits.size(), it);
                    while (bits[it]) {
                        it = sequence_dummy.find(kmer_str, it + 1);
                        if (it == std::string::npos) {
                            reverse_complement(kmer_str.begin(), kmer_str.end());
                            ++rev_count;
                            ASSERT_GT(2u, rev_count);
                            it = sequence_dummy.find(kmer_str);
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

using TAlphabet = KmerExtractor::TAlphabet;

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        // NNN -> $NNN$
        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, {}, &result
            );
            ASSERT_EQ(length - k + 3, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }

    suffix.assign({ KmerExtractor::encode('N') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        // NN   -> $NN$
        // NNN  -> $$NNN$
        // NNNN -> $$$NNNN$
        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, suffix, &result
            );
            ASSERT_EQ(length, result.size()) << "k: " << k
                                             << ", length: " << length;
        }
    }
}


TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithFilteringTwo) {
    for (size_t k = 3; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NN"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        // NNN -> $NNN$
        for (size_t length = k; length < 200; ++length) {
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, { 0, KmerExtractor::encode('N') }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NN"), &result
            );

            // NNN  -> $$NNN$
            // NNNN -> $$$NNNN$
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            KmerExtractor::sequence_to_kmers(
                std::string(length, 'N'), k, KmerExtractor::encode("NA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'N');
            sequence[k - 1] = 'A';

            KmerExtractor::sequence_to_kmers(
                sequence, k, KmerExtractor::encode("NA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringAppend) {
    Vector<TypeParam> result;

    // A...A -> $A...A$
    KmerExtractor::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(501u, result.size());

    KmerExtractor::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(501u * 2, result.size());
}
