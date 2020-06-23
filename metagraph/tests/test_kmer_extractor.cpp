#include "gtest/gtest.h"

#include <numeric>

#include "kmer/kmer_extractor.hpp"
#include "common/utils/string_utils.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace {

using namespace mtg;
using mtg::kmer::KmerExtractor2Bit;

template <typename Kmer>
class ExtractKmers2Bit : public ::testing::Test { };

typedef ::testing::Types<kmer::KMer<uint64_t, KmerExtractor2Bit::bits_per_char>,
                         kmer::KMer<sdsl::uint128_t, KmerExtractor2Bit::bits_per_char>,
                         kmer::KMer<sdsl::uint256_t, KmerExtractor2Bit::bits_per_char>> KmerTypes;

TYPED_TEST_SUITE(ExtractKmers2Bit, KmerTypes);

#define kMaxK ( sizeof(TypeParam) * 8 / KmerExtractor2Bit::bits_per_char )

const KmerExtractor2Bit kmer_extractor;


TEST(KmerExtractor2Bit, encode_decode) {
    KmerExtractor2Bit encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    EXPECT_EQ('A', encoder.decode(encoder.encode('a')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('c')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('g')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('t')));
    ASSERT_THROW(encoder.decode(encoder.encode('N')), std::exception);
    ASSERT_THROW(encoder.decode(encoder.encode('n')), std::exception);
    ASSERT_THROW(encoder.decode(encoder.encode('y')), std::exception);
}

KmerExtractor2Bit::Kmer64 to_kmer(const KmerExtractor2Bit &encoder,
                                  const std::string &kmer) {
    Vector<KmerExtractor2Bit::Kmer64> kmers;
    encoder.sequence_to_kmers(kmer, kmer.size(), {}, &kmers);
    return kmers.at(0);
}

TEST(KmerExtractor2Bit, encode_decode_kmer) {
    KmerExtractor2Bit encoder;
    std::string kmer, kmer_lower, kmer_mixed;

    kmer = "ACGT";
    kmer_lower = "acgt";
    kmer_mixed = "AcGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_lower), kmer_lower.length())) << kmer_lower;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_mixed), kmer_mixed.length())) << kmer_mixed;

    kmer = "AAAAAAAAA";
    kmer_lower = "aaaaaaaaa";
    kmer_mixed = "AaAAAaAAa";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_lower), kmer_lower.length())) << kmer_lower;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_mixed), kmer_mixed.length())) << kmer_mixed;

    kmer = "TTTTTTTTT";
    kmer_lower = "ttttttttt";
    kmer_mixed = "TTtTTTTTt";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length())) << kmer;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_lower), kmer_lower.length())) << kmer_lower;
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(to_kmer(encoder, kmer_mixed), kmer_mixed.length())) << kmer_mixed;

    kmer = "ANANANANANA";
    ASSERT_THROW(encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length()), std::exception);

    kmer = "ANANATANANA";
    ASSERT_THROW(encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length()), std::exception);

    kmer = "ANANANANGNT";
    ASSERT_THROW(encoder.kmer_to_sequence(to_kmer(encoder, kmer), kmer.length()), std::exception);
}

TEST(KmerExtractor2Bit, encode_decode_string) {
    KmerExtractor2Bit encoder;
    std::string first_part = "AAGGCAGCCTAC";
    std::string last_part = "CCCTCTG";
    std::string sequence = first_part + 'N' + last_part;
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        Vector<KmerExtractor2Bit::Kmer256> kmers;

        encoder.sequence_to_kmers(sequence, k, {}, &kmers);
        ASSERT_EQ(k <= last_part.size()
                    ? sequence.size() - 2 * k + 1
                    : (k <= first_part.size() ? first_part.size() - k + 1 : 0),
                  kmers.size()) << k;

        if (!kmers.size())
            continue;

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0], k);
        uint64_t i;
        for (i = 1; i < first_part.size() - k + 1; ++i) {
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

} // namespace
