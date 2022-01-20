#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include <numeric>

#define private public
#define protected public

#include "kmer/kmer_extractor.hpp"
#include "common/utils/string_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace {

using namespace mtg;

using mtg::kmer::KmerExtractorBOSS;

template <typename Kmer>
class ExtractKmers : public ::testing::Test { };

typedef ::testing::Types<KmerExtractorBOSS::Kmer64,
                         KmerExtractorBOSS::Kmer128,
                         KmerExtractorBOSS::Kmer256> KmerTypes;

TYPED_TEST_SUITE(ExtractKmers, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractorBOSS::bits_per_char )


TEST(ExtractKmers, encode_decode) {
    KmerExtractorBOSS encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
#if _DNA_GRAPH
    ASSERT_THROW(encoder.decode(encoder.encode('N')), std::exception);
    ASSERT_THROW(encoder.decode(encoder.encode('X')), std::exception);
    ASSERT_THROW(encoder.decode(encoder.encode(-1)), std::exception);
#elif _PROTEIN_GRAPH
    EXPECT_EQ('X', encoder.decode(encoder.encode('X')));
    EXPECT_EQ('X', encoder.decode(encoder.encode(-1)));
#else
    EXPECT_EQ('N', encoder.decode(encoder.encode('N')));
#endif
}

TEST(ExtractKmers, encode_decode_kmer) {
    KmerExtractorBOSS encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
#if _DNA_GRAPH
    kmer = "ANANANANANA";
    ASSERT_THROW(encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length()), std::exception) << kmer;
    kmer = "ANANATANANA";
    ASSERT_THROW(encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length()), std::exception) << kmer;
    kmer = "ANANANANGNT";
    ASSERT_THROW(encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length()), std::exception) << kmer;
#else
    kmer = "ANANANANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(encoder.sequence_to_kmer<KmerExtractorBOSS::Kmer64>(kmer), kmer.length())) << kmer;
#endif
}

TEST(ExtractKmers, encode_decode_string) {
    KmerExtractorBOSS encoder;
    std::string first_part = "AAGGCAGCCTAC";
    std::string last_part = "CCCTCTG";
    std::string sequence = first_part + 'N' + last_part;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractorBOSS::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);

#if _DNA_GRAPH
        ASSERT_EQ(k <= last_part.size()
                    ? sequence.size() - 2 * k + 1
                    : (k <= first_part.size() ? first_part.size() - k + 1 : 0),
                  kmers.size()) << k;

        if (!kmers.size())
            continue;

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        uint64_t i;
        for (i = 1; i <= first_part.size() - k; ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
        }
        EXPECT_EQ(first_part, reconstructed);

        if (k > last_part.size())
            continue;

        reconstructed = encoder.kmer_to_sequence(kmers[i], k);
        while (++i < kmers.size()) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
        }
        EXPECT_EQ(last_part, reconstructed);
#else
        ASSERT_LT(0u, kmers.size());
        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i], k)[k - 1]);
        }

        EXPECT_EQ(sequence, reconstructed);
#endif
    }
}

TEST(ExtractKmers, encode_decode_string_suffix) {
    KmerExtractorBOSS encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTG";
    std::vector<bool> bits(sequence.size(), false);
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto sequence_dummy = std::string(k - 1, '$') + sequence + "$";
        for (size_t len = 1; len < std::min(k, uint64_t(5)); ++len) {
            bits.assign(sequence.length() + 1, false);
            for (const auto &suffix : utils::generate_strings(encoder.alphabet, len)) {
                std::vector<uint8_t> suffix_encoded;
                std::transform(suffix.begin(), suffix.end(), std::back_inserter(suffix_encoded),
                    [&](char c) { return c == '$' ? 0 : encoder.encode(c); });
                Vector<KmerExtractorBOSS::Kmer256> kmers;
                encoder.sequence_to_kmers(sequence, k, suffix_encoded, &kmers);
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
    KmerExtractorBOSS encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTG";
    std::vector<bool> bits;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto sequence_dummy = std::string(k - 1, '$') + sequence + "$";
        for (size_t len = 1; len < std::min(k - 1, uint64_t(5)); ++len) {
            bits.assign(sequence.length() + 1, false);
            for (const auto &suffix : utils::generate_strings(encoder.alphabet, len)) {
                Vector<KmerExtractorBOSS::Kmer256> kmers;
                std::vector<uint8_t> encoded;
                std::transform(suffix.begin(), suffix.end(), std::back_inserter(encoded),
                    [&](auto c) { return c == '$' ? 0 : encoder.encode(c); });
                encoder.sequence_to_kmers(sequence, k, encoded, &kmers, true);
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

using TAlphabet = KmerExtractorBOSS::TAlphabet;

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithoutFiltering) {
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        // AAA -> $AAA$
        for (size_t length = 0; length < k; ++length) {
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, {}, &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, {}, &result
            );
            ASSERT_EQ(length - k + 1, result.size()) << "k: " << k
                                                     << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringWithFilteringOne) {
    std::vector<TAlphabet> suffix = { 0 };

    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, suffix, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }

    suffix.assign({ KmerExtractorBOSS::encode('A') });
    for (size_t k = 2; k <= kMaxK; ++k) {
        Vector<TypeParam> result;

        for (size_t length = 0; length < k; ++length) {
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, suffix, &result
            );
            ASSERT_TRUE(result.empty());
        }

        // AA   -> $AA$
        // AAA  -> $$AAA$
        // AAAA -> $$$AAAA$
        for (size_t length = k; length < 500; ++length) {
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'A'), k, suffix, &result
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
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractorBOSS::encode("TT"), &result
            );
            ASSERT_TRUE(result.empty()) << "k: " << k
                                        << ", length: " << length;
        }

        // TTT -> $TTT$
        for (size_t length = k; length < 200; ++length) {
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'T'), k, { 0, 0 }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'T'), k, { 0, KmerExtractorBOSS::encode('T') }, &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
            result.clear();
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractorBOSS::encode("TT"), &result
            );

            // TTT  -> $$TTT$
            // TTTT -> $$$TTTT$
            ASSERT_EQ(length - 1, result.size()) << "k: " << k
                                                 << ", length: " << length;
        }

        result.clear();

        for (size_t length = 0; length < k; ++length) {
            KmerExtractorBOSS::sequence_to_kmers(
                std::string(length, 'T'), k, KmerExtractorBOSS::encode("TA"), &result
            );
            ASSERT_TRUE(result.empty());
        }

        for (size_t length = k; length < 200; ++length) {
            result.clear();

            std::string sequence(length, 'T');
            sequence[k - 1] = 'A';

            KmerExtractorBOSS::sequence_to_kmers(
                sequence, k, KmerExtractorBOSS::encode("TA"), &result
            );
            ASSERT_EQ(1u, result.size()) << "k: " << k
                                         << ", length: " << length;
        }
    }
}

TYPED_TEST(ExtractKmers, ExtractKmersFromStringAppend) {
    Vector<TypeParam> result;

    KmerExtractorBOSS::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    // A...AA and A...A$
    ASSERT_EQ(499u, result.size());

    KmerExtractorBOSS::sequence_to_kmers(
        std::string(500, 'A'), 2, {}, &result
    );
    ASSERT_EQ(499u * 2, result.size());
}

} // namespace
