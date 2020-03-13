#include <cstdint>
#include <filesystem>
#include <numeric>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/elias_fano.hpp"

#include "common/utils/file_utils.hpp"

namespace {

using namespace mg;

template <typename T>
class EliasFanoTest : public ::testing::Test {};

typedef ::testing::Types<uint32_t, uint64_t> ValueTypes;

TYPED_TEST_SUITE(EliasFanoTest, ValueTypes);

TYPED_TEST(EliasFanoTest, WriteEmpty) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(0, 0, out.ofstream());
    size_t file_size = encoder.finish();
    // 25 = 3*8 + 1; no data is written to the file except number of low/high bytes (8
    // bytes each), number of low bits (1 byte) and number of elements (8 bytes)
    EXPECT_EQ(25, file_size);
    EXPECT_EQ(25, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadEmpty) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(0, 0, out.ofstream());
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(out.ifstream());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, WriteOne) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(1, 1234, out.ofstream());
    encoder.add(1234);
    size_t file_size = encoder.finish();
    // 1234 is encoded in 3 bytes plus the additional 25 byte header overhead
    EXPECT_EQ(25 + 3, file_size);
    EXPECT_EQ(25 + 3, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadOne) {
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(1, 1234, file.ofstream());
    encoder.add(1234);
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    std::optional<int64_t> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, WriteTwo) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(2, 4321, out.ofstream());
    encoder.add(1234);
    encoder.add(4321);
    size_t file_size = encoder.finish();
    // 1234  and 4321 are encoded in 2 bytes plus the additional 25 byte header overhead
    EXPECT_EQ(25 + 2 * 2, file_size);
    EXPECT_EQ(25 + 2 * 2, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadTwo) {
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(2, 4321, file.ofstream());
    encoder.add(1234);
    encoder.add(4321);
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    std::optional<int64_t> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value());
    EXPECT_TRUE(decoded.has_value());
    decoded = decoder.next();
    EXPECT_EQ(static_cast<TypeParam>(4321), decoded.value());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, ReadWriteIncrementOne) {
    Vector<TypeParam> values(100);
    std::iota(values.begin(), values.end(), 0);
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
    size_t file_size = encoder.finish();
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(25 + (2 * 100) / 8, file_size);
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<int64_t> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(static_cast<TypeParam>(i), decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, ReadWriteIncrementTwo) {
    Vector<TypeParam> values(100);
    uint32_t i = 0;
    std::for_each(values.begin(), values.end(), [&i](TypeParam &v) { v = 2 * i++; });
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
    size_t file_size = encoder.finish();
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<int64_t> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(static_cast<TypeParam>(2 * i), decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, VariousSizes) {
    for (uint32_t size = 100; size < 116; ++size) {
        Vector<TypeParam> values(size);
        std::iota(values.begin(), values.end(), 0);
        utils::TempFile file;
        common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
        size_t file_size = encoder.finish();
        // each value is represented in 2 bits, plus 25 bytes overhead for the header
        EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

        common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
        for (uint32_t i = 0; i < size; ++i) {
            std::optional<int64_t> decoded = decoder.next();
            EXPECT_TRUE(decoded.has_value());
            EXPECT_EQ(static_cast<TypeParam>(i), decoded.value());
        }
        EXPECT_FALSE(decoder.next().has_value());
    }
}

/**
 * These sorted numbers are picked such that the #num_lower_bits_ will be 2, so the last
 * element in the array fills in exactly 64 bits in lower_. This wasy we test if this
 * border-case is handled correctly by the algorithm.
 */
TYPED_TEST(EliasFanoTest, ReadWriteExactly64LowBits) {
    Vector<TypeParam> values = { 1,   5,   7,   12,  16,  17,  25,  31,  32,  37,  40,
                                50,  53,  62,  71,  74,  82,  92,  97,  103, 104, 105,
                                107, 114, 122, 123, 125, 129, 130, 139, 147, 150, 153 };
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
    size_t file_size = encoder.finish();
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    for (uint32_t i = 0; i < values.size(); ++i) {
        std::optional<int64_t> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(values[i], decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, ReadWriteExactly128LowBits) {
    Vector<TypeParam> values
            = { 1,   5,   7,   12,  16,  16,  25,  31,  32,  37,  40,  50,  53,
                62,  71,  74,  82,  92,  97,  103, 104, 105, 107, 114, 122, 123,
                125, 129, 129, 139, 147, 150, 153, 161, 169, 174, 184, 194, 203,
                210, 213, 220, 230, 234, 236, 243, 247, 253, 261, 268, 275, 279,
                287, 289, 296, 298, 299, 300, 307, 311, 317, 321, 326, 333, 338 };
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
    size_t file_size = encoder.finish();
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
    for (uint32_t i = 0; i < values.size(); ++i) {
        std::optional<int64_t> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(values[i], decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, ReadWriteRandom) {
    for (uint32_t size = 1000; size < 1016; ++size) {
        Vector<TypeParam> values(size);
        for (uint32_t trial = 0; trial < 16; ++trial) {
            std::mt19937 rng(123457);
            std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

            uint32_t i = 0;
            std::for_each(values.begin(), values.end(), [&](TypeParam &v) {
                i += dist10(rng);
                v = i;
            });
            utils::TempFile file;
            common::EliasFanoEncoder<TypeParam> encoder(values, file.ofstream());
            size_t file_size = encoder.finish();
            // each value is represented in 2 bits, plus 25 bytes overhead for the header
            EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

            common::EliasFanoDecoder<TypeParam> decoder(file.ifstream());
            for (uint32_t i = 0; i < size; ++i) {
                std::optional<int64_t> decoded = decoder.next();
                EXPECT_TRUE(decoded.has_value());
                EXPECT_EQ(values[i], decoded.value());
            }
            EXPECT_FALSE(decoder.next().has_value());
        }
    }
}

} // namespace
