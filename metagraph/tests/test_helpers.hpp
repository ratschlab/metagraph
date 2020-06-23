#ifndef __TEST_HELPERS_HPP__
#define __TEST_HELPERS_HPP__

#include <set>

#include <gtest/gtest.h>

#include "common/logger.hpp"


namespace {

using namespace mtg;

#ifdef _NO_DEATH_TEST
// Disable death tests

#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(...) (void)0
#define ASSERT_DEATH_SILENT(...) (void)0
#endif

#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(...) (void)0
#define EXPECT_DEATH_SILENT(...) (void)0
#endif

#else
// Define death tests without error messages in log
#define ASSERT_DEATH_SILENT(...) \
    mtg::common::logger->set_level(spdlog::level::critical); \
    ASSERT_DEATH(__VA_ARGS__)
#define EXPECT_DEATH_SILENT(...) \
    mtg::common::logger->set_level(spdlog::level::critical); \
    EXPECT_DEATH(__VA_ARGS__)

#endif // _NO_DEATH_TEST

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

} // namespace

#endif // __TEST_HELPERS_HPP__
