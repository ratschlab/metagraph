#include "common/chunked_wait_queue.hpp"
#include "common/sorted_set_disk.hpp"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>

template <typename T>
class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, int32_t> SortedDiskElementTypes;

TYPED_TEST_CASE(SortedSetDiskTest, SortedDiskElementTypes);

template <typename TypeParam>
void expect_equals(SortedSetDisk<TypeParam> &underTest,
                   const std::vector<TypeParam> &expectedValues) {
    uint32_t size = 0;
    using ChunkedQueueIterator = typename threads::ChunkedWaitQueue<TypeParam>::Iterator;
    threads::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.dataStream();
    for (ChunkedQueueIterator &iterator = merge_queue.iterator();
         iterator != merge_queue.end(); ++iterator) {
        EXPECT_EQ(expectedValues[size], *iterator);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);
}

const std::string out_file = "/tmp/out";
template <typename TypeParam>
void expect_disk_data(const std::string &file_name,
                      const std::vector<TypeParam> &expectedValues) {
    uint32_t size = 0;
    std::ifstream in(file_name, std::ios::binary);
    while (true) {
        TypeParam v;
        if (!in.read(reinterpret_cast<char *>(&v), sizeof(TypeParam))) {
            break;
        }
        EXPECT_LT(size, expectedValues.size());
        EXPECT_EQ(expectedValues[size], v);
        size++;
    }
    EXPECT_EQ(expectedValues.size(), size);

    std::filesystem::remove(out_file);
}

template <typename T>
SortedSetDisk<T> create_sorted_set_disk() {
    constexpr bool verbose = false;
    constexpr size_t thread_count = 1;
    constexpr size_t container_size = 8;
    constexpr size_t merge_queue_size = 1000;
    constexpr size_t merge_queue_backwards_count = 10;
    auto cleanup = [](typename SortedSetDisk<T>::storage_type *) {};
    return SortedSetDisk<T>(cleanup, out_file, thread_count, verbose, container_size,
                            merge_queue_size, merge_queue_backwards_count);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    expect_equals(underTest, {});
    expect_disk_data(out_file, std::vector<TypeParam>());
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 1> elements = { 42 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42 });
    expect_disk_data<TypeParam>(out_file, { 42 });
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42, 43, 44, 45 });
    expect_disk_data<TypeParam>(out_file, { 42, 43, 44, 45 });
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<TypeParam> elements = { 42, 43, 44, 45 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, elements);
    expect_disk_data(out_file, elements);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFiles) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }
    expect_equals(underTest, expected_result);
    expect_disk_data(out_file, expected_result);
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements
                = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
        underTest.insert(elements.begin(), elements.end());
    }
    std::vector<TypeParam> expected_result
            = { TypeParam(0), TypeParam(1), TypeParam(2), TypeParam(3) };
    expect_equals(underTest, expected_result);
    expect_disk_data(out_file, expected_result);
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    expect_equals(underTest, expected_result);
    expect_disk_data(out_file, expected_result);
}

/**
 * Test that elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    expect_equals(underTest, expected_result);
    expect_disk_data(out_file, expected_result);
}

TYPED_TEST(SortedSetDiskTest, IterateBackwards) {
    SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::vector<TypeParam> expected_result;
    for (uint32_t i = 0; i < 100; ++i) {
        std::array<TypeParam, 4> elements = { TypeParam(4 * i), TypeParam(4 * i + 1),
                                              TypeParam(4 * i + 2), TypeParam(4 * i + 3) };
        underTest.insert(elements.begin(), elements.end());
        expected_result.insert(expected_result.end(), elements.begin(), elements.end());
    }

    uint32_t size = 0;
    using ChunkedQueueIterator = typename threads::ChunkedWaitQueue<TypeParam>::Iterator;
    threads::ChunkedWaitQueue<TypeParam> &merge_queue = underTest.dataStream();
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
    expect_disk_data(out_file, expected_result);
}
