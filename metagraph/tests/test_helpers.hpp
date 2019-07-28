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


// Colored logging
// https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework
namespace testing
{
 namespace internal
 {
  enum GTestColor {
      COLOR_DEFAULT,
      COLOR_RED,
      COLOR_GREEN,
      COLOR_YELLOW
  };

  extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
 }
}

#define GTEST_LOG_COLOR testing::internal::COLOR_GREEN

#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(GTEST_LOG_COLOR, __VA_ARGS__); } while(0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
    ~TestCout()
    {
        PRINTF("%s\n",str().c_str());
    }
};

#define TEST_COUT  TestCout()


#endif // __TEST_HELPERS_HPP__
