#include "common/chunked_wait_queue.hpp"
#include "common/sorted_set_disk.hpp"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>

namespace {
using namespace mg;
template <typename T>
class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, int32_t> SortedDiskElementTypes;

TYPED_TEST_CASE(SortedSetDiskTest, SortedDiskElementTypes);

template <typename TypeParam>
void expect_equals(common::SortedSetDisk<TypeParam> &underTest,
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
common::SortedSetDisk<T> create_sorted_set_disk() {
    constexpr bool verbose = false;
    constexpr size_t thread_count = 1;
    constexpr size_t container_size = 8;
    constexpr size_t merge_queue_size = 1000;
    constexpr size_t num_last_elements_cached = 10;
    auto cleanup = [](typename common::SortedSetDisk<T>::storage_type *) {};
    return common::SortedSetDisk<T>(cleanup, out_file, thread_count, verbose,
            container_size,
                            merge_queue_size, num_last_elements_cached);
}

TYPED_TEST(SortedSetDiskTest, Empty) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    expect_equals(underTest, {});
    expect_disk_data(out_file, std::vector<TypeParam>());
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 1> elements = { 42 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42 });
    expect_disk_data<TypeParam>(out_file, { 42 });
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
    std::array<TypeParam, 7> elements = { 43, 42, 42, 45, 44, 45, 43 };
    underTest.insert(elements.begin(), elements.end());
    expect_equals(underTest, { 42, 43, 44, 45 });
    expect_disk_data<TypeParam>(out_file, { 42, 43, 44, 45 });
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    common::SortedSetDisk<TypeParam> underTest = create_sorted_set_disk<TypeParam>();
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
    expect_disk_data(out_file, expected_result);
}
} // namespace
