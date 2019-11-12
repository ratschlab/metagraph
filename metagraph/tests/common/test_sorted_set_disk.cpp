#include "common/sorted_set_disk.hpp"

#include "kmer/kmer_boss.hpp"

#include <gtest/gtest.h>

#include <array>

template <typename T> class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, int32_t> SortedDiskElementTypes;

TYPED_TEST_CASE(SortedSetDiskTest, SortedDiskElementTypes);

TYPED_TEST(SortedSetDiskTest, Empty) {
  SortedSetDisk<TypeParam> underTest;
  EXPECT_EQ(0UL, underTest.data().size());
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
  SortedSetDisk<TypeParam> underTest;
  std::array<TypeParam, 1> elements = {42};
  underTest.insert(elements.begin(), elements.end());
  EXPECT_EQ(1U, underTest.data().size());
  EXPECT_EQ(TypeParam(42), underTest.data()[0]);
}

TYPED_TEST(SortedSetDiskTest, InsertOneRange) {
  SortedSetDisk<TypeParam> underTest;
  std::array<TypeParam, 4> elements = {42, 43, 44, 45};
  underTest.insert(elements.begin(), elements.end());
  EXPECT_EQ(4UL, underTest.data().size());
  // TODO(ddanciu) - use EXPECT_THAT(..., ::testing::ElementsAre(...)); once we
  // can use gmock
  EXPECT_EQ(underTest.data()[0], TypeParam(42));
  EXPECT_EQ(underTest.data()[3], TypeParam(45));
}

/**
 * Test that elements are correctly merged from multiple buffers
 */
TYPED_TEST(SortedSetDiskTest, OneInsertMultipleFiles) {
  SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/,
                                     1 /*threads*/, false /*verbose*/,
                                     8 /* Max buffer size */);
  std::array<TypeParam, 4> elements = {42, 43, 44, 45};
  underTest.insert(elements.begin(), elements.end());
  EXPECT_EQ(4UL, underTest.data().size());
  EXPECT_EQ(underTest.data()[0], TypeParam(42));
  EXPECT_EQ(underTest.data()[3], TypeParam(45));
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFiles) {
  SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/,
                                     1 /*threads*/, false /*verbose*/,
                                     8 /* Max buffer size */);
  for (uint32_t i = 0; i < 100; ++i) {
    std::array<TypeParam, 4> elements = {TypeParam(4 * i), TypeParam(4 * i + 1),
                                         TypeParam(4 * i + 2),
                                         TypeParam(4 * i + 3)};
    underTest.insert(elements.begin(), elements.end());
  }
  EXPECT_EQ(400UL, underTest.data().size());
  for (uint32_t i = 0; i < underTest.data().size(); ++i) {
    ASSERT_EQ(TypeParam(i), underTest.data()[i]);
  }
}

/**
 * Test that non-distinct elements are correctly merged from multiple buffers
 * across multiple inserts.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesNonDistinct) {
  SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/,
                                     1 /*threads*/, false /*verbose*/,
                                     8 /* Max buffer size */);
  for (uint32_t i = 0; i < 100; ++i) {
    std::array<TypeParam, 4> elements = {TypeParam(0), TypeParam(1),
                                         TypeParam(2), TypeParam(3)};
    underTest.insert(elements.begin(), elements.end());
  }
  EXPECT_EQ(4UL, underTest.data().size());
  EXPECT_EQ(underTest.data()[0], TypeParam(0));
  EXPECT_EQ(underTest.data()[3], TypeParam(3));
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreads) {
  SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/,
                                     1 /*threads*/, false /*verbose*/,
                                     8 /* Max buffer size */);
  std::vector<std::thread> workers;
  for (uint32_t i = 0; i < 100; ++i) {
    workers.push_back(std::thread([&underTest, i]() {
      std::array<TypeParam, 4> elements = {
          TypeParam(4 * i), TypeParam(4 * i + 1), TypeParam(4 * i + 2),
          TypeParam(4 * i + 3)};
      underTest.insert(elements.begin(), elements.end());
    }));
  }
  std::for_each(workers.begin(), workers.end(),
                [](std::thread &t) { t.join(); });
  EXPECT_EQ(400UL, underTest.data().size());
  for (uint32_t i = 0; i < underTest.data().size(); ++i) {
    ASSERT_EQ(TypeParam(i), underTest.data()[i]);
  }
}

/**
 * Test that distinct elements are correctly merged from multiple buffers across
 * multiple inserts across multiple threads. Each insert will have dupes.
 */
TYPED_TEST(SortedSetDiskTest, MultipleInsertMultipleFilesMultipleThreadsDupes) {
  SortedSetDisk<TypeParam> underTest([](auto *) {} /*cleanup*/,
                                     1 /*threads*/, false /*verbose*/,
                                     8 /* Max buffer size */);
  std::vector<std::thread> workers;
  for (uint32_t i = 0; i < 100; ++i) {
    workers.push_back(std::thread([&underTest, i]() {
      std::array<TypeParam, 4> elements = {TypeParam(3 * i),
                                           TypeParam(3 * i + 1)};
      underTest.insert(elements.begin(), elements.end());
      elements = {TypeParam(3 * i + 1), TypeParam(3 * i + 2)};
      underTest.insert(elements.begin(), elements.end());
    }));
  }
  std::for_each(workers.begin(), workers.end(),
                [](std::thread &t) { t.join(); });
  EXPECT_EQ(300UL, underTest.data().size());
  for (uint32_t i = 0; i < underTest.data().size(); ++i) {
    ASSERT_EQ(TypeParam(i), underTest.data()[i]);
  }
}
