#include <stdio.h>

#include "gtest/gtest.h"
#include "dbg_succinct_libmaus.hpp"
#include "construct.hpp"

int Sum(int a, int b) {
  return a + b;
}

TEST(SumTest, HandlesZeroInput) {
  // Yoda condition.
  ASSERT_EQ(3, Sum(1, 2));
  ASSERT_EQ(4, Sum(2, 2));
  printf("--- After ASSERT_EQ\n");
  EXPECT_EQ(3, Sum(1, 2));
  printf("--- After EXPECT_EQ\n");
  ASSERT_EQ(4, Sum(2, 2)) << "Something is wrong, 2+2 isn't equal to 4 !";
}

TEST(SumTest, HandlesNegativeInput) {
  ASSERT_EQ(-3, Sum(-1, -2));
}

