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
    common::EliasFanoEncoder<TypeParam> encoder(0, 0, out.name());
    size_t file_size = encoder.finish();
    // 25 = 3*8 + 1; no data is written to the file except number of low/high bytes (8
    // bytes each), number of low bits (1 byte) and number of elements (8 bytes)
    EXPECT_EQ(0U, file_size);
    EXPECT_EQ(0U, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadEmpty) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(0, 0, out.name());
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(out.name());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, WriteOne) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(1, 1234, out.name());
    encoder.add(1234);
    size_t file_size = encoder.finish();
    // 25 = 3*8 + 1; is the overhead, i.e. the number of low/high bytes (8
    // bytes each), number of low bits (1 byte) and number of elements (8 bytes)
    // 1234 is encoded in 3 bytes plus the additional 25 byte header overhead
    EXPECT_EQ(25 + 3U, file_size);
    EXPECT_EQ(25 + 3U, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadOne) {
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(1, 1234, file.name());
    encoder.add(1234);
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    std::optional<TypeParam> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTest, WriteTwo) {
    utils::TempFile out;
    common::EliasFanoEncoder<TypeParam> encoder(2, 4321, out.name());
    encoder.add(1234);
    encoder.add(4321);
    size_t file_size = encoder.finish();
    // 1234  and 4321 are encoded in 2 bytes plus the additional 25 byte header overhead
    EXPECT_EQ(25 + 2 * 2U, file_size);
    EXPECT_EQ(25 + 2 * 2U, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTest, ReadTwo) {
    utils::TempFile file;
    common::EliasFanoEncoder<TypeParam> encoder(2, 4321, file.name());
    encoder.add(1234);
    encoder.add(4321);
    encoder.finish();

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    std::optional<TypeParam> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value());
    EXPECT_TRUE(decoded.has_value());
    decoded = decoder.next();
    EXPECT_EQ(static_cast<TypeParam>(4321), decoded.value());
    EXPECT_FALSE(decoder.next().has_value());
}

template <typename T>
size_t encode(const Vector<T> &values, const std::string &file_name) {
    common::EliasFanoEncoder<T> encoder(values.size(), values.back(), file_name);
    for (const auto &v : values) {
        encoder.add(v);
    }
    return encoder.finish();
}

TYPED_TEST(EliasFanoTest, ReadWriteIncrementOne) {
    Vector<TypeParam> values(100);
    std::iota(values.begin(), values.end(), 0);
    utils::TempFile file;
    size_t file_size = encode(values, file.name());
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(25 + (2 * 100) / 8U, file_size);
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<TypeParam> decoded = decoder.next();
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
    size_t file_size = encode(values, file.name());
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<TypeParam> decoded = decoder.next();
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
        size_t file_size = encode(values, file.name());
        // each value is represented in 2 bits, plus 25 bytes overhead for the header
        EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

        common::EliasFanoDecoder<TypeParam> decoder(file.name());
        for (uint32_t i = 0; i < size; ++i) {
            std::optional<TypeParam> decoded = decoder.next();
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
    size_t file_size = encode(values, file.name());
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < values.size(); ++i) {
        std::optional<TypeParam> decoded = decoder.next();
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
    size_t file_size = encode(values, file.name());
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < values.size(); ++i) {
        std::optional<TypeParam> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(values[i], decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

/**
 * Tests that the last byte of the encoding is correctly read, even if it's at position
 * 8*n+1.
 */
TYPED_TEST(EliasFanoTest, LastByteRead) {
    Vector<TypeParam> values
            = { 0,    585,  587,  588,  601,  604,  609,  612,  650,  651,  658,
                676,  713,  715,  737,  739,  777,  780,  787,  801,  804,  1105,
                1169, 1170, 1172, 1178, 1228, 1242, 1297, 1300, 1313, 1316, 1611,
                1618, 1625, 1634, 1739, 1754, 1755, 1763, 1801, 1812, 1819, 1820,
                2121, 2124, 2129, 2145, 2147, 2148, 2196, 2201, 2203, 2210, 2275,
                2313, 2314, 2316, 2324, 2337, 2338, 2340 };
    utils::TempFile file;
    size_t file_size = encode(values, file.name());
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < values.size(); ++i) {
        std::optional<TypeParam> decoded = decoder.next();
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
            size_t file_size = encode(values, file.name());
            // each value is represented in 2 bits, plus 25 bytes overhead for the header
            EXPECT_EQ(file_size, std::filesystem::file_size(file.name()));

            common::EliasFanoDecoder<TypeParam> decoder(file.name());
            for (uint32_t i = 0; i < size; ++i) {
                std::optional<TypeParam> decoded = decoder.next();
                EXPECT_TRUE(decoded.has_value());
                EXPECT_EQ(values[i], decoded.value());
            }
            EXPECT_FALSE(decoder.next().has_value());
        }
    }
}

template <typename T>
class EliasFanoBufferedTest : public ::testing::Test {};
TYPED_TEST_SUITE(EliasFanoBufferedTest, ValueTypes);

TYPED_TEST(EliasFanoBufferedTest, Empty) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<TypeParam> under_test(file.name(), 100);
    under_test.finish();
    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedTest, InsertOne) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<TypeParam> under_test(file.name(), 100);
    under_test.add(43);
    under_test.finish();
    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    std::optional<TypeParam> decoded = decoder.next();
    ASSERT_TRUE(decoded.has_value());
    EXPECT_EQ(43U, decoded.value());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedTest, InsertFullChunk) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<TypeParam> under_test(file.name(), 3);
    for (uint32_t i = 0; i < 3; ++i) {
        under_test.add(2 * i);
    }
    under_test.finish();
    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < 3; ++i) {
        std::optional<TypeParam> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedTest, InsertTwoChunks) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<TypeParam> under_test(file.name(), 3);
    for (uint32_t i = 0; i < 4; ++i) {
        under_test.add(2 * i);
    }
    under_test.finish();
    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < 4; ++i) {
        std::optional<TypeParam> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedTest, InsertManyChunks) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<TypeParam> under_test(file.name(), 10);
    for (uint32_t i = 0; i < 75; ++i) {
        under_test.add(2 * i);
    }
    under_test.finish();
    common::EliasFanoDecoder<TypeParam> decoder(file.name());
    for (uint32_t i = 0; i < 75; ++i) {
        std::optional<TypeParam> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value());
    }
    EXPECT_FALSE(decoder.next().has_value());
}

} // namespace
