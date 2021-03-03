#include "kmer/kmer_transform.hpp"

#include <random>
#include <vector>

#include "kmer/kmer_extractor.hpp"

#include <gtest/gtest.h>
#include "tests/utils/gtest_patch.hpp"


namespace {
using namespace mtg;

#if ! _PROTEIN_GRAPH

using TAlphabet = kmer::KmerExtractorBOSS::TAlphabet;
constexpr auto &bits_per_char = kmer::KmerExtractorBOSS::bits_per_char;

using KmerTypes = ::testing::Types<kmer::KMerBOSS<uint64_t, bits_per_char>,
                                   kmer::KMerBOSS<sdsl::uint128_t, bits_per_char>,
                                   kmer::KMerBOSS<sdsl::uint256_t, bits_per_char>>;

template <typename Kmer>
class ReverseComplement : public ::testing::Test {};

TYPED_TEST_SUITE(ReverseComplement, KmerTypes);

TYPED_TEST(ReverseComplement, Palindrome) {
    std::vector<uint8_t> seq = { 1, 2, 3, 4 };
    TypeParam kmer_boss(seq); // ACGT
    EXPECT_EQ(kmer_boss,
              kmer::reverse_complement(4, kmer_boss,
                                       kmer::KmerExtractorBOSS::kComplementCode));
}

TYPED_TEST(ReverseComplement, Random) {
    const auto complement_code = kmer::KmerExtractorBOSS::kComplementCode;
    std::mt19937 gen(12345);
    std::uniform_int_distribution<uint64_t> dis(0, complement_code.size() - 1);
    for (uint32_t k = 2; k < sizeof(TypeParam) * 8 / TypeParam::kBitsPerChar; ++k) {
        for (uint32_t trial = 0; trial < 10; ++trial) {
            std::vector<uint8_t> kmer(k);
            std::vector<uint8_t> complement_kmer(k);
            for (uint32_t i = 0; i < k; ++i) {
                kmer[i] = dis(gen);
                complement_kmer[k - i - 1] = complement_code[kmer[i]];
            }
            TypeParam kmer_boss(kmer);
            TypeParam expected(complement_kmer);
            EXPECT_EQ(expected, kmer::reverse_complement(k, kmer_boss, complement_code));
        }
    }
}

#endif // ! _PROTEIN_GRAPH

} // namespace
