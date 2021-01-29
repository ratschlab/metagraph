#include "common/sorted_sets/sorted_set_disk.hpp"
#include "common/threads/chunked_wait_queue.hpp"
#include "common/utils/file_utils.hpp"

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
void expect_equals(common::SortedSetDisk<TypeParam> &under_test,
                   const std::vector<TypeParam> &expectedValues) {
    uint32_t size = 0;
    common::ChunkedWaitQueue<TypeParam> &merge_queue = under_test.data();
    for (auto &it = merge_queue.begin(); it != merge_queue.end(); ++it) {
        EXPECT_EQ(expectedValues[size], *it);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
}

template <typename T>
common::SortedSetDisk<T> create_sorted_set_disk(const std::string &tmp_dir,
                                                size_t container_size = 8,
                                                size_t merge_count = 4) {
    constexpr size_t thread_count = 1;
    constexpr size_t disk_cap_bytes = 1e6;
    return common::SortedSetDisk<T>(thread_count, container_size, tmp_dir,
                                    disk_cap_bytes, merge_count);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test = create_sorted_set_disk<TypeParam>(tmp_dir);
    expect_equals(under_test, {});
    std::filesystem::remove(tmp_dir);
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::array<TypeParam, 1> elements = { 42 };
    under_test.insert(elements.begin(), elements.end());
    expect_equals(under_test, { 42 });
    std::filesystem::remove(tmp_dir);
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    under_test.insert(elements.begin(), elements.end());
    expect_equals(under_test, { 42, 43, 44, 45 });
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    under_test.insert(elements.begin(), elements.end());
    expect_equals(under_test, elements);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that we correctly deal with inserting a range larger than the buffer's size.
 */
TYPED_TEST(SortedSetDiskTest, OneInsertLargerThanBuffer) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    constexpr size_t container_size = 3;
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir, container_size);
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    under_test.insert(elements.begin(), elements.end());
    expect_equals(under_test, elements);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFiles) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    constexpr size_t container_size = 8;
    for (size_t merge_count : {0,4}) {
        common::SortedSetDisk<TypeParam> under_test
                = create_sorted_set_disk<TypeParam>(tmp_dir, container_size, merge_count);
        std::vector<TypeParam> expected_result;
        for (uint32_t i = 0; i < 100; ++i) {
            std::array<TypeParam, 4> elements
                    = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                        TypeParam(4 * i + 3) };
            under_test.insert(elements.begin(), elements.end());
            expected_result.insert(expected_result.end(), elements.begin(), elements.end());
        }
        expect_equals(under_test, expected_result);
    }
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    constexpr size_t container_size = 8;
    for (size_t merge_count : {0,4}) {
        common::SortedSetDisk<TypeParam> under_test
                = create_sorted_set_disk<TypeParam>(tmp_dir, container_size, merge_count);
        for (uint32_t i = 0; i < 100; ++i) {
            std::array<TypeParam, 4> elements
                    = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
            under_test.insert(elements.begin(), elements.end());
        }
        std::vector<TypeParam> expected_result
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        expect_equals(under_test, expected_result);
    }
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.emplace_back([&under_test, i]() {
            std::array<TypeParam, 4> elements
                = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                    TypeParam(4 * i + 3) };
            under_test.insert(elements.begin(), elements.end());
        });
    }
    for (uint32_t i = 0; i < 400; ++i) {
        expected_result.emplace_back(i);
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.emplace_back([&under_test, i]() {
            std::array<TypeParam, 4> elements = { TypeParam(3 * i), TypeParam(3 * i + 1) };
            under_test.insert(elements.begin(), elements.end());
            elements = { TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
            under_test.insert(elements.begin(), elements.end());
        });
    }
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 3> elements
                = { TypeParam(3 * i), TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test that exceeding the allocated disk space and then merging all data to reduce
 * space works correctly.
 */
TYPED_TEST(SortedSetDiskTest, DiskExceeded) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    constexpr size_t thread_count = 1;
    constexpr size_t reserved_num_elements = 100;
    constexpr size_t disk_cap_bytes = 100;
    auto under_test = common::SortedSetDisk<TypeParam>(thread_count, reserved_num_elements,
                                                      tmp_dir, disk_cap_bytes);
    std::vector<TypeParam> elements(100);
    std::iota(elements.begin(), elements.end(), 0);
    for(uint32_t i = 0; i<10;++i) {
        under_test.insert(elements.begin(), elements.end());
    }
    expect_equals(under_test, elements);
    std::filesystem::remove(tmp_dir);
}

TYPED_TEST(SortedSetDiskTest, InsertSortedOnly) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                            TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        under_test.insert_sorted(elements);
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Overlap) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                            TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        under_test.insert_sorted(elements);
        under_test.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Distinct) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test
            = create_sorted_set_disk<TypeParam>(tmp_dir);
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements1 = { TypeParam(4 * i), TypeParam(4 * i + 1) };
        std::vector<TypeParam> elements2 = { TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        under_test.insert_sorted(elements1);
        under_test.insert(elements2.begin(), elements2.end());
    }
    std::vector<TypeParam> expected_result(400);
    std::iota(expected_result.begin(), expected_result.end(), 0);
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

/**
 * Test #insert_sorted combined with #insert gives correct results
 */
TYPED_TEST(SortedSetDiskTest, InsertSortedAndInsert_Intertwined) {
    std::filesystem::path tmp_dir = utils::create_temp_dir("", "test_ssd");
    common::SortedSetDisk<TypeParam> under_test = create_sorted_set_disk<TypeParam>(tmp_dir);
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<TypeParam> elements1 = { TypeParam(4 * i), TypeParam(4 * i + 2) };
        std::vector<TypeParam> elements2 = { TypeParam(4 * i + 1), TypeParam(4 * i + 3) };
        under_test.insert_sorted(elements1);
        under_test.insert(elements2.begin(), elements2.end());
    }
    std::vector<TypeParam> expected_result(400);
    std::iota(expected_result.begin(), expected_result.end(), 0);
    expect_equals(under_test, expected_result);
    std::filesystem::remove(tmp_dir);
}

} // namespace
