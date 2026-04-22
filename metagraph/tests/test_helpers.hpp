#ifndef __TEST_HELPERS_HPP__
#define __TEST_HELPERS_HPP__

#include <cstdlib>
#include <filesystem>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <unistd.h>

#include <gtest/gtest.h>

#include "common/logger.hpp"


// Per-process scratch directory for unit-test dump outputs. Created on first
// call under $TEST_DUMP_DIR (or the system tempdir) and removed at exit.
// Must live outside the anonymous namespace so the static local is shared
// across translation units — so a serialize in one TU and load in another hit
// the same path.
//
// Note: uses mkdtemp directly rather than utils::create_temp_dir because
// test_dump_basename-style constants at namespace scope in the test TUs
// trigger this during static initialization, and create_temp_dir logs via
// mtg::common::logger — which sits in another TU and may not be initialized
// yet (classic SIOF). Direct mkdtemp has no such dependency.
//
// The atexit handler guards on pid to avoid wiping the directory from forked
// children — gtest's threadsafe death-test style forks, and a naive cleanup
// from the child would pull the dir out from under the parent's later tests.
inline const std::string& test_dump_dir() {
    static const pid_t owner_pid = getpid();
    static const std::string dir = [] {
        namespace fs = std::filesystem;
        if (const char* env = std::getenv("TEST_DUMP_DIR"); env && *env) {
            fs::create_directories(env);
            return std::string(env);
        }
        std::string tmpl = (fs::temp_directory_path() / "metagraph_tests_XXXXXX").string();
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (!mkdtemp(buf.data()))
            throw std::runtime_error("mkdtemp failed: " + tmpl);
        std::atexit([] {
            if (getpid() != owner_pid)
                return;
            std::error_code ec;
            fs::remove_all(test_dump_dir(), ec);
        });
        return std::string(buf.data());
    }();
    return dir;
}


namespace {

using namespace mtg;

// run debug death tests only for debug to check assertions
#ifdef NDEBUG
#undef ASSERT_DEBUG_DEATH
#define ASSERT_DEBUG_DEATH(...) (void)0
#endif

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

// C++ stream interface
class TestCout : public std::stringstream {
  public:
    ~TestCout() {
        fmt::print("[          ] {:s}\n", str().c_str());
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
