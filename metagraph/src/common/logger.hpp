#pragma once
#include <spdlog/fmt/ostr.h> // for logging custom classes
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace mtg {
namespace common {

bool get_verbose();
void set_verbose(bool verbose);

extern std::shared_ptr<spdlog::logger> logger;

} // namespace common

#ifndef NDEBUG
#define DEBUG_LOG(...) common::logger->trace(__VA_ARGS__)
#else
#define DEBUG_LOG(...) (void)0
#endif

} // namespace mtg
