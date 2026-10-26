#pragma once

#include <chrono>
#include <functional>
#include <utility>

namespace mtg::common {

// A request's work deadline reserves time for serializing its answer. The
// injected clock lets callers test both boundaries without sleeping.
class Deadline {
  public:
    using Clock = std::chrono::steady_clock;

    Deadline(Clock::time_point start, double time_budget_ms, double finalize_reserve_ms,
             std::function<Clock::time_point()> clock = &Clock::now)
          : start_(start), time_budget_ms_(time_budget_ms),
            finalize_reserve_ms_(finalize_reserve_ms),
            clock_(clock ? std::move(clock) : std::function<Clock::time_point()>(&Clock::now)) {}

    bool work_expired() const { return elapsed_ms() >= time_budget_ms_ - finalize_reserve_ms_; }
    bool respond_expired() const { return elapsed_ms() >= time_budget_ms_; }
    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(clock_() - start_).count();
    }
    double time_budget_ms() const { return time_budget_ms_; }
    double finalize_reserve_ms() const { return finalize_reserve_ms_; }

  private:
    Clock::time_point start_;
    double time_budget_ms_;
    double finalize_reserve_ms_;
    std::function<Clock::time_point()> clock_;
};

} // namespace mtg::common
