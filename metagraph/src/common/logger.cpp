#include "logger.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <string>

#include <spdlog/sinks/stdout_color_sinks.h>

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

    spdlog::sinks::stdout_color_sink_mt out_;
    spdlog::sinks::stderr_color_sink_mt err_;
};

static std::shared_ptr<spdlog::logger> make_logger() {
    spdlog::set_automatic_registration(false);
    auto sink = std::make_shared<split_sink>();
    auto logger = std::make_shared<spdlog::logger>("", sink);
    // Honor SPDLOG_LEVEL at construction time so any log calls made during
    // static/dynamic init of other TUs (e.g. tests creating temp dirs) are
    // filtered at the intended level, not the default of `info`.
    if (const char* env = std::getenv("SPDLOG_LEVEL"); env && *env) {
        // spdlog::level::from_str is case-sensitive; normalize the env value
        // so users can set SPDLOG_LEVEL=TRACE or SPDLOG_LEVEL=trace.
        std::string level_env(env);
        std::transform(level_env.begin(), level_env.end(), level_env.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        auto lvl = spdlog::level::from_str(level_env);
        // from_str returns `off` for unrecognized names, so distinguish a real
        // "off" from a parse failure.
        bool recognized = (lvl != spdlog::level::off) || (level_env == "off");
        if (recognized) {
            logger->set_level(lvl);
        } else {
            std::fprintf(stderr, "Warning: failed to parse SPDLOG_LEVEL='%s'; "
                         "expected one of trace|debug|info|warn|error|critical|off\n", env);
        }
    }
    return logger;
}

std::shared_ptr<spdlog::logger>& get_logger() {
    static std::shared_ptr<spdlog::logger> instance = make_logger();
    return instance;
}

} // namespace common
} // namespace mtg
