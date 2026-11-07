#ifndef __TESTS_SOURCE_PATH_HPP__
#define __TESTS_SOURCE_PATH_HPP__

#include <string>

// A file or directory of the source tree, given relative to metagraph/ ("tests/data/...",
// "docs/...", "api/..."), as an absolute path: a test reads it from whatever directory it runs
// in. unit_tests_parallel runs each shard in a directory of its own, where a path relative to
// the build directory ("../tests/data/...") does not resolve. METAGRAPH_SOURCE_DIR (metagraph/)
// is compiled into unit_tests by CMake
inline std::string source_path(const std::string &path) {
    return std::string(METAGRAPH_SOURCE_DIR) + "/" + path;
}

#endif // __TESTS_SOURCE_PATH_HPP__
