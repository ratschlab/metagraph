#pragma once
#include <spdlog/fmt/ostr.h> // for logging custom classes
#include <spdlog/spdlog.h>

namespace mtg {
namespace common {
extern std::shared_ptr<spdlog::logger> logger;
} // namespace common
} // namespace mtg
