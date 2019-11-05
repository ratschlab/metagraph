#include "kmer_boss.hpp"
#include "sorted_set_disk.hpp"

#include <gtest/gtest.h>

template <typename T> class SortedSetDiskTest : public ::testing::Test {};

typedef ::testing::Types<uint64_t, KMerBOSS<uint64_t, 4>>
    SortedDiskElementTypes;

TYPED_TEST_CASE(SortedSetDiskTest, SortedDiskElementTypes);

TYPED_TEST(SortedSetDiskTest, Empty) {
  SortedSetDisk<TypeParam> underTest;
  EXPECT_EQ(0UL, underTest.data().size());
}

TYPED_TEST(SortedSetDiskTest, InsertOneElement) {
  SortedSetDisk<TypeParam> underTest;
  std::array<TypeParam, 1> elements;
  underTest.insert(elements.begin(), elements.end());
  EXPECT_EQ(1UL, underTest.data().size());
}
