#include "logger.hpp"

#include <spdlog/sinks/stdout_color_sinks.h>

namespace mtg {
namespace common {

static bool VERBOSE = false;

bool get_verbose() { return VERBOSE; }
void set_verbose(bool verbose) { VERBOSE = verbose; }

std::shared_ptr<spdlog::logger> logger = spdlog::default_logger();

} // namespace common
} // namespace mtg
