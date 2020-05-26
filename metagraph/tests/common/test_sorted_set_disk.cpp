#include "common/sorted_set_disk.hpp"
#include "common/threads/chunked_wait_queue.hpp"

#include <gtest/gtest.h>

#include "tests/utils/gtest_patch.hpp"

#include <array>
#include <filesystem>

namespace {
using namespace mg;

template <typename T>
class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, sdsl::uint128_t> SortedDiskElementTypes;

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
common::SortedSetDisk<T> create_sorted_set_disk(size_t container_size = 8,
                                                size_t num_elements_cached = 2) {
    constexpr size_t thread_count = 1;
    constexpr size_t max_disk_space = 1e6;
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    auto nocleanup = [](typename common::SortedSetDisk<T>::storage_type *) {};
    return common::SortedSetDisk<T>(nocleanup, thread_count, container_size,
                                    "./test_chunk_", max_disk_space,
                                    num_elements_cached);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    constexpr size_t container_size = 8;
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
        expect_equals(underTest, {});
    }
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    constexpr size_t container_size = 8;
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
        std::array<TypeParam, 1> elements = { 42 };
        underTest.insert(elements.begin(), elements.end());
        expect_equals(underTest, { 42 });
    }
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    constexpr size_t container_size = 8;
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
        std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
        underTest.insert(elements.begin(), elements.end());
        expect_equals(underTest, { 42, 43, 44, 45 });
    }
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    constexpr size_t container_size = 8;
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
        std::vector<TypeParam> elements = { 42, 43, 44, 45 };
        underTest.insert(elements.begin(), elements.end());
        expect_equals(underTest, elements);
    }
}

/**
 * Test that we correctly deal with inserting a range larger than the buffer's size.
 */
TYPED_TEST(SortedSetDiskTest, OneInsertLargerThanBuffer) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>(3, 1);
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
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        underTest.insert(elements.begin(), elements.end());
    }
    std::vector<TypeParam> expected_result
            = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
    expect_equals(underTest, expected_result);
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

TYPED_TEST(SortedSetDiskTest, IterateBackwards) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>(100, 10);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }

    uint32_t size = 0;
    using ChunkedQueueIterator = typename common::ChunkedWaitQueue<TypeParam>::Iterator;
    common::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.data();
    for (ChunkedQueueIterator &iterator = merge_queue.begin();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ((TypeParam)size, *iterator);
        for (uint32_t idx = 0; idx < 10 && size - idx > 0; ++idx) {
            EXPECT_EQ((TypeParam)(size - idx), *iterator)
                    << "Size: " << size << " Index: " << idx;
            --iterator;
        }
        TypeParam value = *iterator;
        for (uint32_t idx = 0; idx < 10 && size - idx > 0; ++idx) {
            EXPECT_EQ((TypeParam)(value + idx), *iterator);
            ++iterator;
        }
        size++;
    }
}

TYPED_TEST(SortedSetDiskTest, IterateBackwardsFromEnd) {
    constexpr size_t container_size = 8;
    for (uint32_t cached = 1; cached < container_size / 3; ++cached) {
        common::SortedSetDisk<TypeParam> underTest
                = create_sorted_set_disk<TypeParam>(container_size, cached);
        std::array<TypeParam, 2> elements = { 42, 45 };
        underTest.insert(elements.begin(), elements.end());
        using ChunkedQueueIterator = typename common::ChunkedWaitQueue<TypeParam>::Iterator;
        common::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.data();
        ChunkedQueueIterator &iterator = merge_queue.begin();
        for (; iterator != merge_queue.end(); ++iterator) {
        }
        --iterator;
        EXPECT_EQ(TypeParam(45), *iterator);
        --iterator;
        EXPECT_EQ(TypeParam(42), *iterator);
        ++iterator;
        ++iterator;
        EXPECT_EQ(merge_queue.end(), iterator);
    }
}

/**
 * Test that exceeding the allocated disk space and then merging all data to reduce
 * space works correctly.
 */
TYPED_TEST(SortedSetDiskTest, DiskExceeded) {
    constexpr size_t thread_count = 1;
    constexpr size_t reserved_num_elements = 100;
    constexpr size_t max_disk_space = 100;
    auto nocleanup = [](typename common::SortedSetDisk<TypeParam>::storage_type *) {};
    std::filesystem::create_directory("./test_chunk_");
    std::atexit([]() { std::filesystem::remove_all("./test_chunk_"); });
    auto underTest = common::SortedSetDisk<TypeParam>(nocleanup, thread_count,
                                                      reserved_num_elements,
                                                      "./test_chunk_", max_disk_space, 1);
    std::vector<TypeParam> elements(100);
    std::iota(elements.begin(), elements.end(), 0);
    for(uint32_t i = 0; i<10;++i) {
        underTest.insert(elements.begin(), elements.end());
    }
    expect_equals(underTest, elements);
}

} // namespace
