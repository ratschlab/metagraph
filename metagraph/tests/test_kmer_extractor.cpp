#include "gtest/gtest.h"

#define protected public
#define private public

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
