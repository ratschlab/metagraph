#include "common/sorted_multiset_disk.hpp"
#include "common/threads/chunked_wait_queue.hpp"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>

namespace {
using namespace mg;

template <typename T>
class SortedMultisetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, int32_t> SortedDiskElementTypes;

TYPED_TEST_SUITE(SortedMultisetDiskTest, SortedDiskElementTypes);

template <typename T>
static void nocleanup(T *) {}

template <typename TypeParam>
void expect_equals(common::SortedMultisetDisk<TypeParam, uint8_t> &underTest,
                   const std::vector<std::pair<TypeParam, uint8_t>> &expectedValues) {
    using Pair = std::pair<TypeParam, uint8_t>;
    uint32_t size = 0;
    using ChunkedQueueIterator = typename common::ChunkedWaitQueue<Pair>::Iterator;
    common::ChunkedWaitQueue<Pair> &merge_queue = underTest.data();
    for (ChunkedQueueIterator &iterator = merge_queue.begin();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ(expectedValues[size], *iterator);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
    merge_queue.shutdown();
}

template <typename T>
common::SortedMultisetDisk<T, uint8_t>
create_sorted_set_disk(size_t container_size = 8, size_t num_elements_cached = 2) {
    constexpr size_t thread_count = 1;
    auto on_item_pushed = [](const std::pair<T, uint8_t> &) {};
    return common::SortedMultisetDisk<T, uint8_t>(
            nocleanup<typename common::SortedMultisetDisk<T>::storage_type>, thread_count,
            container_size, "/tmp/test_chunk_", on_item_pushed, num_elements_cached);
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
            = create_sorted_set_disk<TypeParam>(3, 1);
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
            = create_sorted_set_disk<TypeParam>(100, 3);
    std::vector<std::thread> workers;
    std::vector<std::pair<TypeParam, uint8_t>> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        workers.push_back(std::thread([&underTest, i]() {
            std::array<TypeParam, 2> elements = { TypeParam(3 * i), TypeParam(3 * i + 1) };
            underTest.insert(elements.begin(), elements.end());
            elements = { TypeParam(3 * i + 1), TypeParam(3 * i + 2) };
            underTest.insert(elements.begin(), elements.end());
        }));

        expected_result.push_back(std::make_pair(TypeParam(3 * i), 1));
        expected_result.push_back(std::make_pair(TypeParam(3 * i + 1), 2));
        expected_result.push_back(std::make_pair(TypeParam(3 * i + 2), 1));
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
    expect_equals(underTest, expected_result);
}

TYPED_TEST(SortedMultisetDiskTest, IterateBackwards) {
    common::SortedMultisetDisk<TypeParam, uint8_t> underTest
            = create_sorted_set_disk<TypeParam>(100, 10);
    using Pair = std::pair<TypeParam, uint8_t>;
    std::vector<Pair> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        for (const auto &el : elements) {
            expected_result.push_back(std::make_pair(el, 1));
        }
    }

    uint32_t size = 0;
    using ChunkedQueueIterator = typename common::ChunkedWaitQueue<Pair>::Iterator;
    common::ChunkedWaitQueue<Pair> &merge_queue = underTest.data();
    for (ChunkedQueueIterator &iterator = merge_queue.begin();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ((TypeParam)size, (*iterator).first);
        for (uint32_t idx = 0; idx < 10 && size - idx > 0; ++idx) {
            EXPECT_EQ((TypeParam)(size - idx), (*iterator).first);
            --iterator;
        }
        Pair value = *iterator;
        for (uint32_t idx = 0; idx < 10 && size - idx > 0; ++idx) {
            EXPECT_EQ((TypeParam)(value.first + idx), (*iterator).first);
            ++iterator;
        }
        size++;
    }
}

/**
 * Test that count overflows are handled correctly in file_merger, i.e. the counter stays
 * at the maximum possible value rather than rolling over.
 */
TYPED_TEST(SortedMultisetDiskTest, CounterOverflowAtMergeDisk) {
    constexpr uint32_t value = 12342341;
    // make sure we correctly count up to the max value of the counter
    common::SortedMultisetDisk<TypeParam, uint8_t> underTest
            = create_sorted_set_disk<TypeParam>();
    for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max(); ++idx) {
        std::vector<TypeParam> values = { TypeParam(value) };
        underTest.insert(values.begin(), values.end());
    }
    expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });

    // now let's generate an overflow in the counter
    underTest.clear();
    for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max() + 10; ++idx) {
        std::vector<TypeParam> values = { TypeParam(value) };
        underTest.insert(values.begin(), values.end());
    }
    expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
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
    auto on_item_pushed = [](const std::pair<TypeParam, uint8_t> &) {};
    common::SortedMultisetDisk<TypeParam, uint8_t> underTest(
            nocleanup<typename common::SortedMultisetDisk<TypeParam>::storage_type>,
            thread_count, container_size, "/tmp/test_chunk_", on_item_pushed, 2);
    // make sure we correctly count up to the max value of the counter
    for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max(); ++idx) {
        std::vector<TypeParam> values = { TypeParam(value) };
        underTest.insert(values.begin(), values.end());
    }
    expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });

    // now let's generate an overflow in the counter
    underTest.clear();
    for (uint32_t idx = 0; idx < std::numeric_limits<uint8_t>::max() + 10; ++idx) {
        std::vector<TypeParam> values = { TypeParam(value) };
        underTest.insert(values.begin(), values.end());
    }
    expect_equals(underTest, { std::make_pair(TypeParam(value), 255) });
}

} // namespace
