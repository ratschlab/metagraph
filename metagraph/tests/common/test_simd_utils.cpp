#include "common/utils/simd_utils.hpp"

#include <cmath>

#include <gtest/gtest.h>


namespace {

TEST(SIMDUtils, restrict_to) {
    // test multiples of powers of two, and powers of two minus 1
    for (uint64_t i = 0; i < 64; ++i) {
        for (uint64_t j = 0; j < 64; ++j) {
            ASSERT_EQ(i + j >= 64 ? 1llu << (i + j - 64) : 0,
                      restrict_to(1llu << i, 1llu << j));
        }
    }

    ASSERT_EQ(0UL, restrict_to(0, 0));
    ASSERT_EQ(0UL, restrict_to(0, 0xFFFFFFFFFFFFFFFF));
    ASSERT_EQ(0UL, restrict_to(0xFFFFFFFFFFFFFFFF, 0));
    ASSERT_EQ(0xFFFFFFFFFFFFFFFEUL, restrict_to(0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF));

    // test multiples of powers of 10
    constexpr double log10maxull = 19.2659197225;
    for (uint64_t i = 0; i < 11; ++i) {
        for (uint64_t j = 0; j < 11; ++j) {
            // compute upper half of product in log10 space
            ASSERT_EQ(i + j >= log10maxull
                          ? std::floor(std::pow(10, i + j - log10maxull))
                          : 0,
                      restrict_to(std::pow(10, i), std::pow(10, j)));
        }
    }
}

} // namespace
