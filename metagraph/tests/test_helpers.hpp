#ifndef __TEST_HELPERS_HPP__
#define __TEST_HELPERS_HPP__

#include <cstdlib>
#include <filesystem>
#include <set>
#include <string>

#include <gtest/gtest.h>

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "graph/representation/succinct/boss_construct.hpp"


// Per-process scratch directory for unit-test dump outputs. Shared across
// translation units via the function-local static, so a serialize in one TU
// and a load in another hit the same path.
//
// When $TEST_DUMP_DIR is set (e.g. by scripts/unit_tests_parallel.sh), the
// caller owns the dir's lifecycle. Otherwise delegate to utils::create_temp_dir,
// which installs the PID-guarded atexit + signal cleanup.
inline const std::string& test_dump_dir() {
    static const std::string dir = [] {
        namespace fs = std::filesystem;
        if (const char* env = std::getenv("TEST_DUMP_DIR"); env && *env) {
            fs::create_directories(env);
            mtg::common::logger->trace("Using TEST_DUMP_DIR={}", env);
            return std::string(env);
        }
        return utils::create_temp_dir(fs::temp_directory_path(),
                                      "metagraph_tests").string();
    }();
    return dir;
}

// BOSSConstructor factory that pins swap_dir to test_dump_dir(). Every
// BOSS::Chunk created during construction scribbles buffered int_vectors
// (W/last/weights) into that swap_dir — including empty Chunks — because
// those members are sdsl::int_vector_buffer and always disk-backed. Using
// the production default "/tmp/" leaks dirs to /tmp on any abnormal exit;
// routing through test_dump_dir() lets the atexit sweep pick them up.
inline mtg::graph::boss::BOSSConstructor make_test_boss_constructor(
        size_t k,
        bool both_strands = false,
        uint8_t bits_per_count = 0,
        const std::string &filter_suffix = "",
        size_t indexed_suffix_length = 0,
        size_t num_threads = 1,
        double memory_preallocated = 0,
        mtg::kmer::ContainerType container = mtg::kmer::ContainerType::VECTOR) {
    return mtg::graph::boss::BOSSConstructor(
        k, both_strands, bits_per_count, filter_suffix,
        indexed_suffix_length, num_threads, memory_preallocated,
        container, test_dump_dir());
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
