#ifndef __TEST_HELPERS_HPP__
#define __TEST_HELPERS_HPP__

#include <set>


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
namespace testing {
    namespace internal {
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

#define PRINTF(...) do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(GTEST_LOG_COLOR, __VA_ARGS__); } while(0)

// C++ stream interface
class TestCout : public std::stringstream {
  public:
    ~TestCout() {
        PRINTF("%s\n", str().c_str());
    }
};

#define TEST_COUT TestCout()

template <typename T>
inline std::set<T> convert_to_set(const std::vector<T> &vector) {
    return std::set<T>(vector.begin(), vector.end());
}

template <class Container>
inline std::set<typename Container::value_type> convert_to_set(const Container &vector) {
    return std::set<typename Container::value_type>(vector.begin(), vector.end());
}

// To support calls 'convert_to_set({ "string_1", "string_2" });'
inline std::set<std::string> convert_to_set(const std::vector<std::string> &vector) {
    return std::set<std::string>(vector.begin(), vector.end());
}

#endif // __TEST_HELPERS_HPP__
