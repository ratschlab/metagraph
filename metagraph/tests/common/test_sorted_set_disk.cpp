#include "common/sorted_sets/sorted_set_disk.hpp"
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
class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t,
                         sdsl::uint128_t,
                         sdsl::uint256_t> SortedDiskElementTypes;

TYPED_TEST_SUITE(SortedSetDiskTest, SortedDiskElementTypes);


template <typename TypeParam>
void expect_equals(common::SortedSetDisk<TypeParam> &underTest,
                   const std::vector<TypeParam> &expectedValues) {
    uint32_t size = 0;
    common::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.data();
    for (auto &it = merge_queue.begin(); it != merge_queue.end(); ++it) {
        EXPECT_EQ(expectedValues[size], *it);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
}

template <typename T>
common::SortedSetDisk<T> create_sorted_set_disk(size_t container_size = 8, size_t merge_count = 4) {
    constexpr size_t thread_count = 1;
    constexpr size_t max_disk_space = 1e6;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    return common::SortedSetDisk<T>(thread_count, container_size,
                                    "./test_chunk_", max_disk_space, merge_count);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    std::filesystem::remove_all("./test_chunk_");
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    expect_equals(underTest, {});
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    std::array<TypeParam, 1> elements = { 42 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42 });
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42, 43, 44, 45 });
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, elements);
}

/**
 * Test that we correctly deal with inserting a range larger than the buffer's size.
 */
TYPED_TEST(SortedSetDiskTest, OneInsertLargerThanBuffer) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>(3);
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, elements);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFiles) {
    constexpr size_t container_size = 8;
    for (size_t merge_count : {0,4}) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, merge_count);
        std::vector<TypeParam> expected_result;
        for (uint32_t i = 0; i < 100; ++i) {
            std::array<TypeParam, 4> elements
                    = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                        TypeParam(4 * i + 3) };
            underTest.insert(elements.begin(), elements.end());
            expected_result.insert(expected_result.end(), elements.begin(), elements.end());
        }
        expect_equals(underTest, expected_result);
    }
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
    for (size_t merge_count : {0,4}) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(8, merge_count);
        for (uint32_t i = 0; i < 100; ++i) {
            std::array<TypeParam, 4> elements
                    = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
            underTest.insert(elements.begin(), elements.end());
        }
        std::vector<TypeParam> expected_result
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        expect_equals(underTest, expected_result);
    }
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.emplace_back([&underTest, i]() {
            std::array<TypeParam, 4> elements
                = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                    TypeParam(4 * i + 3) };
            underTest.insert(elements.begin(), elements.end());
        });
    }
    for (uint32_t i = 0; i < 400; ++i) {
        expected_result.emplace_back(i);
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(underTest, expected_result);
}

/**
 * Test that elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.emplace_back([&underTest, i]() {
            std::array<TypeParam, 4> elements = { TypeParam(3 * i), TypeParam(3 * i + 1) };
            underTest.insert(elements.begin(), elements.end());
            elements = { TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
            underTest.insert(elements.begin(), elements.end());
        });
    }
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 3> elements
                = { TypeParam(3 * i), TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(underTest, expected_result);
}

/**
 * Test that exceeding the allocated disk space and then merging all data to reduce
 * space works correctly.
 */
TYPED_TEST(SortedSetDiskTest, DiskExceeded) {
    constexpr size_t thread_count = 1;
    constexpr size_t reserved_num_elements = 100;
    constexpr size_t max_disk_space = 100;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    auto underTest = common::SortedSetDisk<TypeParam>(thread_count, reserved_num_elements,
                                                      "./test_chunk_", max_disk_space);
    std::vector<TypeParam> elements(100);
    std::iota(elements.begin(), elements.end(), 0);
    for(uint32_t i = 0; i<10;++i) {
        underTest.insert(elements.begin(), elements.end());
    }
    expect_equals(underTest, elements);
}

TYPED_TEST(SortedSetDiskTest, InsertSortedOnly) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                            TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert_sorted(elements);
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expect_equals(underTest, expected_result);
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Overlap) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                            TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert_sorted(elements);
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expect_equals(underTest, expected_result);
    std::filesystem::remove_all("./test_chunk_");
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Distinct) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements1 = { TypeParam(4 * i), TypeParam(4 * i + 1) };
        std::vector<TypeParam> elements2 = { TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert_sorted(elements1);
        underTest.insert(elements2.begin(), elements2.end());
    }
    std::vector<TypeParam> expected_result(400);
    std::iota(expected_result.begin(), expected_result.end(), 0);
    expect_equals(underTest, expected_result);
    std::filesystem::remove_all("./test_chunk_");
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Intertwined) {
    constexpr size_t container_size = 8;
    common::SortedSetDisk<TypeParam> underTest
            = create_sorted_set_disk<TypeParam>(container_size);
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements1 = { TypeParam(4 * i), TypeParam(4 * i + 2) };
        std::vector<TypeParam> elements2 = { TypeParam(4 * i + 1), TypeParam(4 * i + 3) };
        underTest.insert_sorted(elements1);
        underTest.insert(elements2.begin(), elements2.end());
    }
    std::vector<TypeParam> expected_result(400);
    std::iota(expected_result.begin(), expected_result.end(), 0);
    expect_equals(underTest, expected_result);
    std::filesystem::remove_all("./test_chunk_");
}

} // namespace
