#include "logger.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace mg {
namespace common {
std::shared_ptr<spdlog::logger> logger = spdlog::stdout_color_mt("mg");
} // namespace common
} // namespace mg
