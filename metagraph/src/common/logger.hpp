#pragma once
#include <memory>

#include <spdlog/fmt/ostr.h> // for logging custom classes
#include <spdlog/fmt/std.h>
#include <spdlog/fmt/fmt.h>
#include <spdlog/spdlog.h>

namespace mtg {
namespace common {

bool get_verbose();
void set_verbose(bool verbose);

// Function-local static ("Meyer's singleton"): initialized on first call,
// which sidesteps static-initialization-order issues for other TUs that
// reach for the logger during their own dynamic initialization.
std::shared_ptr<spdlog::logger>& get_logger();

// Thin proxy that forwards to get_logger(). Zero-sized and aggregate-default
// constructible, so `inline LoggerProxy logger{};` is constant-initialized
// and safe to reference from any other TU's dynamic initializer — every
// access routes through get_logger() which triggers the Meyer's singleton.
struct LoggerProxy {
    spdlog::logger *operator->() const { return get_logger().get(); }
    spdlog::logger &operator*() const { return *get_logger(); }
    operator std::shared_ptr<spdlog::logger>&() const { return get_logger(); }
};

inline LoggerProxy logger{};

} // namespace common

#ifndef NDEBUG
#define DEBUG_LOG(...) common::logger->trace(__VA_ARGS__)
#else
#define DEBUG_LOG(...) (void)0
#endif

} // namespace mtg
