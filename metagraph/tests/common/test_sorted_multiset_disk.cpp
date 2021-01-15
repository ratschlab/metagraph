#include "common/sorted_sets/sorted_multiset_disk.hpp"
#include "common/threads/chunked_wait_queue.hpp"

#include <gtest/gtest.h>

#include "tests/utils/gtest_patch.hpp"

#include <array>
#include <numeric>
#include <filesystem>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace {

using namespace mtg;

template <typename T>
class SortedMultisetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t,
                         sdsl::uint128_t,
                         sdsl::uint256_t> SortedDiskElementTypes;

TYPED_TEST_SUITE(SortedMultisetDiskTest, SortedDiskElementTypes);

template <typename TypeParam>
void expect_equals(common::SortedMultisetDisk<TypeParam, uint8_t> &underTest,
                   const std::vector<std::pair<TypeParam, uint8_t>> &expectedValues) {
    using Pair = std::pair<TypeParam, uint8_t>;
    uint32_t size = 0;
    common::ChunkedWaitQueue<Pair> &merge_queue = underTest.data();
    for (auto &it = merge_queue.begin(); it != merge_queue.end(); ++it) {
        EXPECT_EQ(expectedValues[size], *it);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
}

template <typename T>
common::SortedMultisetDisk<T, uint8_t> create_sorted_set_disk(size_t container_size = 8) {
    constexpr size_t thread_count = 1;
    constexpr size_t max_disk_space = 1e6;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    return common::SortedMultisetDisk<T, uint8_t>(
            thread_count, container_size, "./test_chunk_", max_disk_space);
}

TYPED_TEST(SortedMultisetDiskTest, Empty) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    expect_equals(underTest, {});
}

TYPED_TEST(SortedMultisetDiskTest, InsertOneElement) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 1> elements = { 42 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { { 42, 1 } });
}

TYPED_TEST(SortedMultisetDiskTest, InsertOneRange) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { { 42, 2 }, { 43, 2 }, { 44, 1 }, { 45, 2 } });
}

TYPED_TEST(SortedMultisetDiskTest, InsertOneRangeLargerThanBuffer) {
    common::SortedMultisetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(3);
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { { 42, 2 }, { 43, 2 }, { 44, 1 }, { 45, 2 } });
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedMultisetDiskTest, OneInsertMultipleFiles) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    for (uint32_t buf = 0; buf < 4; ++buf) {
        std::vector<TypeParam> values(underTest.buffer_size());
        std::iota(values.begin(), values.end(), 0);
        underTest.insert(values.begin(), values.end());
    }
    auto &data = underTest.data();
    uint32_t i = 0;
    for (auto &it = data.begin(); it != data.end(); ++it) {
        ASSERT_EQ(std::make_pair((TypeParam)i, (uint8_t)4), *it);
        i++;
    }
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedMultisetDiskTest, MultipleInsertMultipleFiles) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<std::pair<TypeParam, uint8_t>> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        for (const auto &el : elements) {
            expected_result.push_back(std::make_pair(el, (uint8_t)1));
        }
    }
    expect_equals(underTest, expected_result);
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedMultisetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        underTest.insert(elements.begin(), elements.end());
    }
    std::vector<std::pair<TypeParam, uint8_t>> expected_result = { { TypeParam(0), 100 },
                                                                   { TypeParam(1), 100 },
                                                                   { TypeParam(2), 100 },
                                                                   { TypeParam(3), 100 } };
    expect_equals(underTest, expected_result);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedMultisetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
    common::SortedMultisetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<std::thread> workers;
    std::vector<std::pair<TypeParam, uint8_t>> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.push_back(std::thread([&underTest, i]() {
            std::array<TypeParam, 4> elements
                    = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                        TypeParam(4 * i + 3) };
            underTest.insert(elements.begin(), elements.end());
        }));
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        for (const auto &el : elements) {
            expected_result.push_back(std::make_pair(el, 1));
        }
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(underTest, expected_result);
}

/**
 * Test that elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedMultisetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
    common::SortedMultisetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(100);
    std::vector<std::thread> workers;
    std::vector<std::pair<TypeParam, uint8_t>> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.emplace_back([&underTest, i]() {
            std::array<TypeParam, 2> elements = { TypeParam(3 * i), TypeParam(3 * i + 1) };
            underTest.insert(elements.begin(), elements.end());
            elements = { TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
            underTest.insert(elements.begin(), elements.end());
        });

        expected_result.push_back(std::make_pair(TypeParam(3 * i), 1));
        expected_result.push_back(std::make_pair(TypeParam(3 * i + 1), 2));
        expected_result.push_back(std::make_pair(TypeParam(3 * i + 2), 1));
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(underTest, expected_result);
}

/**
 * Test that count overflows are handled correctly in file_merger, i.e. the counter stays
 * at the maximum possible value rather than rolling over.
 */
TYPED_TEST(SortedMultisetDiskTest, CounterOverflowAtMergeDisk) {
    constexpr uint32_t value = 12342341;
    // make sure we correctly count up to the max value of the counter
    // the container size is 8, so we are guaranteed to generate many chunk files that
    // will have to be merged and overflow handled correctly
    {
        common::SortedMultisetDisk<TypeParam, uint8_t> underTest
                = create_sorted_set_disk<TypeParam>();
        std::vector<TypeParam> values = { TypeParam(value) };
        for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max(); ++idx) {
            underTest.insert(values.begin(), values.end());
        }
        expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
    }

    // now let's generate an overflow in the counter
    {
        common::SortedMultisetDisk<TypeParam, uint8_t> underTest
                = create_sorted_set_disk<TypeParam>();
        for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max() + 10; ++idx) {
            std::vector<TypeParam> values = { TypeParam(value) };
            underTest.insert(values.begin(), values.end());
        }
        expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
    }
}

/**
 * Test that count overflows are also handled correctly when de-duping in memory
 * (before writing to disk).
 */
TYPED_TEST(SortedMultisetDiskTest, CounterOverflowAtMergeMemory) {
    constexpr uint32_t value = 12342341;
    // make sure the container is large enough to hold all values - this way we are
    // guaranteed to test de-duping in memory rather than on disk
    constexpr uint32_t container_size = 2 * std::numeric_limits<uint8_t>::max();

    constexpr size_t thread_count = 1;
    constexpr size_t max_disk_space = 1e6;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    {
        common::SortedMultisetDisk<TypeParam, uint8_t> underTest(
                thread_count, container_size, "./test_chunk_", max_disk_space);
        // make sure we correctly count up to the max value of the counter
        for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max(); ++idx) {
            std::vector<TypeParam> values = { TypeParam(value) };
            underTest.insert(values.begin(), values.end());
        }
        expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
    }

    // now let's generate an overflow in the counter
    {
        common::SortedMultisetDisk<TypeParam, uint8_t> underTest(
                thread_count, container_size, "./test_chunk_", max_disk_space);
        for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max() + 10; ++idx) {
            std::vector<TypeParam> values = { TypeParam(value) };
            underTest.insert(values.begin(), values.end());
        }
        expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
    }
}

/**
 * Test that reaching the maximum allowed disk space is handled correctly: the data is
 * de-duped and building can continue with the reduced disk space
 */
TYPED_TEST(SortedMultisetDiskTest, ExhaustMaxAllowedDiskSpace) {
    // make sure the container is large enough to hold all values - this way we are
    // guaranteed to test de-duping in memory rather than on disk
    constexpr uint32_t container_size = 256;

    constexpr size_t thread_count = 1;
    constexpr size_t max_disk_space
            = container_size * sizeof(std::pair<uint32_t, uint8_t>) * 2;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    common::SortedMultisetDisk<TypeParam, uint8_t> underTest(
            thread_count, container_size, "./test_chunk_", max_disk_space);
    std::vector<TypeParam> values(container_size / 1.5);
    std::iota(values.begin(), values.end(), 0);
    // these values will fill the buffer and write to disk filling half the allowed space
    underTest.insert(values.begin(), values.end());
    // now we fill the other half
    underTest.insert(values.begin(), values.end());
    // now the disk space will be exhausted, but it's okay, because duplicates will be
    // reduced
    underTest.insert(values.begin(), values.end());
    underTest.insert(values.begin(), values.end());

    auto &data = underTest.data();
    uint32_t i = 0;
    for (auto &it = data.begin(); it != data.end(); ++it) {
        ASSERT_EQ(std::make_pair((TypeParam)i, (uint8_t)4), *it);
        i++;
    }
}

} // namespace
