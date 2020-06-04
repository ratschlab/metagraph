#include <cstdint>
#include <filesystem>
#include <numeric>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include <sdsl/uint128_t.hpp>

#include "common/elias_fano.hpp"
#include "common/utils/file_utils.hpp"
#include "tests/utils/gtest_patch.hpp"


namespace {

using namespace mtg;

template <typename T>
class EliasFanoTestPair : public ::testing::Test {};

typedef ::testing::Types<uint64_t, sdsl::uint128_t, sdsl::uint256_t> ValueTypes;

TYPED_TEST_SUITE(EliasFanoTestPair, ValueTypes);

TYPED_TEST(EliasFanoTestPair, WriteEmpty) {
    utils::TempFile out;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(0, 0, 0, out.name());
    size_t file_size = encoder.finish();
    // 25 = 3*8 + 1; no data is written to the file except number of low/high bytes (8
    // bytes each), number of low bits (1 byte) and number of elements (8 bytes)
    EXPECT_EQ(0U, file_size);
    EXPECT_EQ(0U, std::filesystem::file_size(out.name()));
}

TYPED_TEST(EliasFanoTestPair, ReadEmpty) {
    utils::TempFile out;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(0, 0, 0, out.name());
    encoder.finish();

    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(out.name());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTestPair, WriteOne) {
    utils::TempFile out;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(1, 1234, 1234,
                                                                    out.name());
    encoder.add({ 1234, 5 });
    size_t total_size = encoder.finish();
    // 24 + 1 + size(T); is the overhead, i.e. the number of low/high bytes (8
    // bytes each), number of low bits (1 byte), the number of elements (8 bytes) and the
    // offset (sizeof(T)).
    // 1234 is encoded in 1 byte plus the additional header overhead; another byte is used
    // for the count
    EXPECT_EQ(25 + sizeof(TypeParam) + 1U + 1U, total_size);
    EXPECT_EQ(25 + sizeof(TypeParam) + 1U + 1U,
              std::filesystem::file_size(out.name())
                      + std::filesystem::file_size(out.name() + ".count")
                      + std::filesystem::file_size(out.name() + ".up"));
}

TYPED_TEST(EliasFanoTestPair, ReadOne) {
    utils::TempFile file;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(1, 1234, 1234,
                                                                    file.name());
    encoder.add({ 1234, 5 });
    encoder.finish();

    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value().first);
    EXPECT_EQ(static_cast<TypeParam>(5), decoded.value().second);
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTestPair, WriteTwo) {
    utils::TempFile out;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(2, 1234, 4321,
                                                                    out.name());
    encoder.add({ 1234, 5 });
    encoder.add({ 4321, 3 });
    size_t file_size = encoder.finish();
    // 1234  and 4321 are encoded in 4 bytes plus the additional 25 byte header overhead
    // + 2 bytes for the counts
    EXPECT_EQ(25 + sizeof(TypeParam) + 4U + 2U, file_size);
    EXPECT_EQ(25 + sizeof(TypeParam) + 4U + 2U,
              std::filesystem::file_size(out.name())
                      + std::filesystem::file_size(out.name() + ".count")
                      + std::filesystem::file_size(out.name() + ".up"));
}

TYPED_TEST(EliasFanoTestPair, ReadTwo) {
    utils::TempFile file;
    common::EliasFanoEncoder<std::pair<TypeParam, uint8_t>> encoder(2, 1234, 4321,
                                                                    file.name());
    encoder.add({ 1234, 5 });
    encoder.add({ 4321, 3 });
    encoder.finish();

    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
    EXPECT_TRUE(decoded.has_value());
    EXPECT_EQ(static_cast<TypeParam>(1234), decoded.value().first);
    EXPECT_EQ(static_cast<TypeParam>(5), decoded.value().second);
    EXPECT_TRUE(decoded.has_value());
    decoded = decoder.next();
    EXPECT_EQ(static_cast<TypeParam>(4321), decoded.value().first);
    EXPECT_EQ(static_cast<TypeParam>(3), decoded.value().second);
    EXPECT_FALSE(decoder.next().has_value());
}


template <typename T, typename C>
size_t encode(const std::vector<std::pair<T, C>> &values, const std::string &file_name) {
    common::EliasFanoEncoder<std::pair<T, C>> encoder(values.size(), values.front().first,
                                                      values.back().first, file_name);
    for (const auto &v : values) {
        encoder.add(v);
    }
    return encoder.finish();
}

TYPED_TEST(EliasFanoTestPair, ReadWriteIncrementOne) {
    std::vector<std::pair<TypeParam, uint8_t>> values(100);
    for (uint32_t i = 0; i < 100; ++i) {
        values[i] = { i, 100 + i };
    }
    utils::TempFile file;
    size_t file_size = encode(values, file.name());
    // each value is represented in 2 bits, plus 25 bytes overhead for the header
    EXPECT_EQ(25 + sizeof(TypeParam) + (2 * 100) / 8U + 100, file_size);
    EXPECT_EQ(file_size,
              std::filesystem::file_size(file.name())
                      + std::filesystem::file_size(file.name() + ".count")
                      + std::filesystem::file_size(file.name() + ".up"));

    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(static_cast<TypeParam>(i), decoded.value().first);
        EXPECT_EQ(100 + i, decoded.value().second);
    }
    EXPECT_FALSE(decoder.next().has_value());
}


TYPED_TEST(EliasFanoTestPair, ReadWriteIncrementTwo) {
    std::vector<std::pair<TypeParam, uint8_t>> values(100);
    uint32_t i = 0;
    std::for_each(values.begin(), values.end(), [&i](std::pair<TypeParam, uint8_t> &v) {
        v = { 2 * i, i };
        i++;
    });
    utils::TempFile file;
    size_t file_size = encode(values, file.name());
    EXPECT_EQ(file_size,
              std::filesystem::file_size(file.name())
                      + std::filesystem::file_size(file.name() + ".count")
                      + std::filesystem::file_size(file.name() + ".up"));

    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    for (uint32_t i = 0; i < 100; ++i) {
        std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
        EXPECT_TRUE(decoded.has_value());
        EXPECT_EQ(static_cast<TypeParam>(2 * i), decoded.value().first);
        EXPECT_EQ(i, decoded.value().second);
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoTestPair, VariousSizes) {
    for (uint32_t size = 100; size < 116; ++size) {
        std::vector<std::pair<TypeParam, uint8_t>> values(size);
        for (uint32_t i = 0; i < size; ++i) {
            values[i] = { i, size + i };
        }
        utils::TempFile file;
        size_t file_size = encode(values, file.name());
        // each value is represented in 2 bits, plus 25 bytes overhead for the header
        EXPECT_EQ(file_size,
                  std::filesystem::file_size(file.name())
                          + std::filesystem::file_size(file.name() + ".count")
                          + std::filesystem::file_size(file.name() + ".up"));

        common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
        for (uint32_t i = 0; i < size; ++i) {
            std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
            EXPECT_TRUE(decoded.has_value());
            EXPECT_EQ(static_cast<TypeParam>(i), decoded.value().first);
            EXPECT_EQ(size + i, decoded.value().second);
        }
        EXPECT_FALSE(decoder.next().has_value());
    }
}

TYPED_TEST(EliasFanoTestPair, ReadWriteRandom) {
    for (uint32_t size = 1000; size < 1016; ++size) {
        std::vector<std::pair<TypeParam, uint16_t>> values(size);
        for (uint32_t trial = 0; trial < 16; ++trial) {
            std::mt19937 rng(123457);
            std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

            uint32_t i = 0;
            std::for_each(values.begin(), values.end(),
                          [&](std::pair<TypeParam, uint16_t> &v) {
                              i += dist10(rng);
                              v = { i, size };
                          });
            utils::TempFile file;
            size_t file_size = encode(values, file.name());
            // each value is represented in 2 bits, plus 25 bytes overhead for the header
            EXPECT_EQ(file_size,
                      std::filesystem::file_size(file.name())
                              + std::filesystem::file_size(file.name() + ".count")
                              + std::filesystem::file_size(file.name() + ".up"));

            common::EliasFanoDecoder<std::pair<TypeParam, uint16_t>> decoder(file.name());
            for (uint32_t i = 0; i < size; ++i) {
                std::optional<std::pair<TypeParam, uint16_t>> decoded = decoder.next();
                EXPECT_TRUE(decoded.has_value());
                EXPECT_EQ(values[i].first, decoded.value().first);
                EXPECT_EQ(values[i].second, decoded.value().second);
            }
            EXPECT_FALSE(decoder.next().has_value());
        }
    }
}


template <typename T>
class EliasFanoBufferedPairTest : public ::testing::Test {};
TYPED_TEST_SUITE(EliasFanoBufferedPairTest, ValueTypes);

TYPED_TEST(EliasFanoBufferedPairTest, Empty) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<std::pair<TypeParam, uint8_t>> under_test(file.name(),
                                                                               100);
    under_test.finish();
    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedPairTest, InsertOne) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<std::pair<TypeParam, uint8_t>> under_test(file.name(),
                                                                               100);
    under_test.add({ 43, 15 });
    under_test.finish();
    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
    ASSERT_TRUE(decoded.has_value());
    EXPECT_EQ(43U, decoded.value().first);
    EXPECT_EQ(15U, decoded.value().second);
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedPairTest, InsertFullChunk) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<std::pair<TypeParam, uint8_t>> under_test(file.name(),
                                                                               3);
    for (uint32_t i = 0; i < 3; ++i) {
        under_test.add({ 2 * i, i });
    }
    under_test.finish();
    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    for (uint32_t i = 0; i < 3; ++i) {
        std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value().first);
        EXPECT_EQ(i, decoded.value().second);
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedPairTest, InsertTwoChunks) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<std::pair<TypeParam, uint8_t>> under_test(file.name(),
                                                                               3);
    for (uint32_t i = 0; i < 4; ++i) {
        under_test.add({ 2 * i, i });
    }
    under_test.finish();
    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    for (uint32_t i = 0; i < 4; ++i) {
        std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value().first);
        EXPECT_EQ(i, decoded.value().second);
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TYPED_TEST(EliasFanoBufferedPairTest, InsertManyChunks) {
    utils::TempFile file;
    common::EliasFanoEncoderBuffered<std::pair<TypeParam, uint8_t>> under_test(file.name(),
                                                                               10);
    for (uint32_t i = 0; i < 75; ++i) {
        under_test.add({ 2 * i, i });
    }
    under_test.finish();
    common::EliasFanoDecoder<std::pair<TypeParam, uint8_t>> decoder(file.name());
    for (uint32_t i = 0; i < 75; ++i) {
        std::optional<std::pair<TypeParam, uint8_t>> decoded = decoder.next();
        ASSERT_TRUE(decoded.has_value());
        EXPECT_EQ(2 * i, decoded.value().first);
        EXPECT_EQ(i, decoded.value().second);
    }
    EXPECT_FALSE(decoder.next().has_value());
}

TEST(EliasFanoTestPair128, ReadWriteRandom) {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

    for (uint32_t size = 1000; size < 1016; ++size) {
        std::vector<std::pair<sdsl::uint128_t, uint16_t>> values(size);
        for (uint32_t trial = 0; trial < 16; ++trial) {
            uint32_t i = 0;
            std::for_each(values.begin(), values.end(),
                          [&](std::pair<sdsl::uint128_t, uint16_t> &v) {
                              i += dist10(rng);
                              v = { i, trial };
                          });
            utils::TempFile file;
            size_t file_size = encode(values, file.name());
            // each value is represented in 2 bits, plus 25 bytes overhead for the header
            EXPECT_EQ(file_size,
                      std::filesystem::file_size(file.name())
                              + std::filesystem::file_size(file.name() + ".count")
                              + std::filesystem::file_size(file.name() + ".up"));

            common::EliasFanoDecoder<std::pair<sdsl::uint128_t, uint16_t>> decoder(
                    file.name());
            for (uint32_t i = 0; i < size; ++i) {
                std::optional<std::pair<sdsl::uint128_t, uint16_t>> decoded = decoder.next();
                EXPECT_TRUE(decoded.has_value());
                EXPECT_EQ(values[i].first, decoded.value().first);
                EXPECT_EQ(values[i].second, decoded.value().second);
            }
            EXPECT_FALSE(decoder.next().has_value());
        }
    }
}

TEST(EliasFanoTestPair128, ReadWriteRandomLarge) {
    std::mt19937 rng(123457);
    std::uniform_int_distribution<std::mt19937::result_type> dist10(0, 10);

    for (uint32_t size = 1000; size < 1016; ++size) {
        std::vector<std::pair<sdsl::uint128_t, uint16_t>> values(size);
        for (uint32_t trial = 0; trial < 16; ++trial) {
            sdsl::uint128_t i = 0;
            std::for_each(values.begin(), values.end(),
                          [&](std::pair<sdsl::uint128_t, uint16_t> &v) {
                              i += dist10(rng);
                              if (dist10(rng)
                                  < 1) { // increase the hi 64 bits every ~10th element
                                  i = ((sdsl::uint128_t)((uint64_t)(i >> 64) + 1) << 64)
                                          + (uint64_t)i;
                              }
                              v = { i, size };
                          });
            utils::TempFile file;
            size_t file_size = encode(values, file.name());
            // each value is represented in 2 bits, plus 25 bytes overhead for the header
            EXPECT_EQ(file_size,
                      std::filesystem::file_size(file.name())
                              + std::filesystem::file_size(file.name() + ".count")
                              + std::filesystem::file_size(file.name() + ".up"));

            common::EliasFanoDecoder<std::pair<sdsl::uint128_t, uint16_t>> decoder(
                    file.name());
            for (uint32_t i = 0; i < size; ++i) {
                std::optional<std::pair<sdsl::uint128_t, uint16_t>> decoded = decoder.next();
                EXPECT_TRUE(decoded.has_value());
                EXPECT_EQ(values[i].first, decoded.value().first);
                EXPECT_EQ(values[i].second, decoded.value().second);
            }
            EXPECT_FALSE(decoder.next().has_value());
        }
    }
}

} // namespace
