#ifndef __TEST_HELPERS_HPP__
#define __TEST_HELPERS_HPP__


// Disable death tests
#ifdef _NO_DEATH_TEST

#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif

#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(a, b) (void)0
#endif

#endif // _NO_DEATH_TEST

#endif // __TEST_HELPERS_HPP__
