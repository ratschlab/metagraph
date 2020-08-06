#include "logger.hpp"

#include <spdlog/sinks/stdout_sinks.h>

namespace mtg {
namespace common {

static bool VERBOSE = false;

bool get_verbose() { return VERBOSE; }
void set_verbose(bool verbose) { VERBOSE = verbose; }

class split_sink : public spdlog::sinks::sink {
    void log(const spdlog::details::log_msg &msg) {
        if (msg.level == spdlog::level::info) {
            out_.log(msg);
        } else {
            err_.log(msg);
        }
    }

    void flush() {
        out_.flush();
        err_.flush();
    }

    void set_pattern(const std::string &pattern) {
        out_.set_pattern(pattern);
        err_.set_pattern(pattern);
    }

    void set_formatter(std::unique_ptr<spdlog::formatter> sink_formatter) {
        out_.set_formatter(sink_formatter->clone());
        err_.set_formatter(std::move(sink_formatter));
    }

    spdlog::sinks::stdout_sink_mt out_;
    spdlog::sinks::stderr_sink_mt err_;
};

std::shared_ptr<spdlog::logger> make_logger() {
    spdlog::set_automatic_registration(false);
    auto sink = std::make_shared<split_sink>();
    return std::make_shared<spdlog::logger>("", sink);
}

std::shared_ptr<spdlog::logger> logger = make_logger();

} // namespace common
} // namespace mtg
