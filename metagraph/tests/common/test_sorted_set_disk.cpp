#include "common/chunked_wait_queue.hpp"
#include "common/sorted_set_disk.hpp"

#include <gtest/gtest.h>

#include <array>

template <typename T>
class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, int32_t> SortedDiskElementTypes;

TYPED_TEST_CASE(SortedSetDiskTest, SortedDiskElementTypes);

template <typename TypeParam>
void expectEquals(SortedSetDisk<TypeParam> &underTest,
                  const std::vector<TypeParam> &expectedValues) {
    uint32_t size = 0;
    using ChunkedQueueIterator = typename common::ChunkedWaitQueue<TypeParam>::Iterator;
    common::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.dataStream();
    for (ChunkedQueueIterator &iterator = merge_queue.iterator();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ(expectedValues[size], *iterator);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    SortedSetDisk<TypeParam> underTest;
    expectEquals(underTest, {});
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    SortedSetDisk<TypeParam> underTest;
    std::array<TypeParam, 1> elements = { 42 };
    underTest.insert(elements.begin(), elements.end());
    expectEquals(underTest, { 42 });
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    SortedSetDisk<TypeParam> underTest;
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expectEquals(underTest, { 42, 43, 44, 45 });
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    auto cleanup = [](typename SortedSetDisk<TypeParam>::storage_type *) {};
    SortedSetDisk<TypeParam> underTest(cleanup, 1 /*threads*/, false /*verbose*/,
                                       8 /* Max buffer size */);
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    underTest.insert(elements.begin(), elements.end());
    expectEquals(underTest, elements);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFiles) {
    auto cleanup = [](typename SortedSetDisk<TypeParam>::storage_type *) {};
    SortedSetDisk<TypeParam> underTest(cleanup, 1 /*threads*/, false /*verbose*/,
                                       8 /* Max buffer size */);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expectEquals(underTest, expected_result);
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
    auto cleanup = [](typename SortedSetDisk<TypeParam>::storage_type *) {};
    SortedSetDisk<TypeParam> underTest(cleanup, 1 /*threads*/, false /*verbose*/,
                                       8 /* Max buffer size */);
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        underTest.insert(elements.begin(), elements.end());
    }
    expectEquals(underTest, { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) });
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
    auto cleanup = [](typename SortedSetDisk<TypeParam>::storage_type *) {};
    SortedSetDisk<TypeParam> underTest(cleanup, 1 /*threads*/, false /*verbose*/,
                                       8 /* Max buffer size */);
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.push_back(std::thread([&underTest, i]() {
            std::array<TypeParam, 4> elements
                    = { TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
                        TypeParam(4 * i + 3) };
            underTest.insert(elements.begin(), elements.end());
        }));
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expectEquals(underTest, expected_result);
}

/**
 * Test that elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
    SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/, 1 /*threads*/,
                                       false
                                       /*verbose*/,
                                       8 /* Max buffer size */);
    std::vector<std::thread> workers;
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.push_back(std::thread([&underTest, i]() {
            std::array<TypeParam, 4> elements = { TypeParam(3 * i), TypeParam(3 * i + 1) };
            underTest.insert(elements.begin(), elements.end());
            elements = { TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
            underTest.insert(elements.begin(), elements.end());
        }));
        std::array<TypeParam, 3> elements
                = { TypeParam(3 * i), TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expectEquals(underTest, expected_result);
}

TYPED_TEST(SortedSetDiskTest, IterateBackwards) {
    auto cleanup = [](typename SortedSetDisk<TypeParam>::storage_type *) {};
    SortedSetDisk<TypeParam> underTest(cleanup, 1 /*threads*/, false /*verbose*/,
                                       8 /* Max buffer size */);
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }

    uint32_t size = 0;
    using ChunkedQueueIterator = typename common::ChunkedWaitQueue<TypeParam>::Iterator;
    common::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.dataStream();
    for (ChunkedQueueIterator &iterator = merge_queue.iterator();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ((TypeParam)size, *iterator);
        for (uint32_t idx = 0; idx < 10 && size - idx > 0; ++idx) {
            EXPECT_EQ((TypeParam)(size - idx), *iterator);
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
