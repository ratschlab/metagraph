#include "traverse_attempts.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <random>

#include <spdlog/fmt/fmt.h>

#include "json_helpers.hpp"
#include "traverse.hpp"


namespace mtg {
namespace cli {

using graph::traversal::ExternalStop;
using Clock = Attempt::Clock;

namespace {

// The rules of the capabilities' attempts block, stated by reference to the SPEC section that
// holds each (the documents a service returns in one piece have a ceiling of 32 KiB; the
// numbers a ledger computes with are the block's fields). ASCII only
constexpr char kBoundRule[] = "SPEC-labeled-traversal-core.md section 6.8, the attempt's bound";
constexpr char kDeliveryReserveRule[]
        = "SPEC-labeled-traversal-core.md section 6.8, the delivery reserve";
constexpr char kCalibrationRule[]
        = "SPEC-labeled-traversal-core.md section 6.8, the delivery reserve: calibration";
constexpr char kNotAfterRule[] = "SPEC-labeled-traversal-core.md section 5, not_after_ms";
constexpr char kInstanceRule[]
        = "SPEC-labeled-traversal-core.md section 5, expect_server_instance";
constexpr char kSuppressionRule[]
        = "SPEC-labeled-traversal-core.md section 10.3, POST /traverse/cancel: the tombstone";
constexpr char kReleaseRule[] = "SPEC-labeled-traversal-core.md section 10.3, the release rule";

// whole milliseconds, rounded up (a duration a ledger leases by is never understated); null
// when there is no bound (an infinite time budget without the HTTP server's cap: the CLI)
Json::Value ceil_ms_json(double ms) {
    if (!std::isfinite(ms))
        return Json::Value();
    return uint_json(static_cast<uint64_t>(std::ceil(std::max(0.0, ms))));
}

// Writes the ids of |ids| into |j|: attempt_id (always with |always_attempt_id|, else when
// given), budget_id and locus_id when given, and not_after_ms when given
void put_attempt_ids(Json::Value *j, const AttemptIds &ids, bool always_attempt_id) {
    if (always_attempt_id || !ids.attempt_id.empty())
        (*j)["attempt_id"] = ids.attempt_id;
    if (!ids.budget_id.empty())
        (*j)["budget_id"] = ids.budget_id;
    if (!ids.locus_id.empty())
        (*j)["locus_id"] = ids.locus_id;
    if (ids.not_after_ms)
        (*j)["not_after_ms"] = uint_json(*ids.not_after_ms);
}

// Appends |x| to |window|, which keeps the last kRateWindow values
void push_window(std::deque<double> *window, double x) {
    window->push_back(x);
    if (window->size() > kRateWindow)
        window->pop_front();
}

// the smaller of two measurements where both were made (> 0), else the one that was made
double prefer_measured(double a, double b) {
    return a > 0 && b > 0 ? std::min(a, b) : std::max(a, b);
}

// |unknown| (AttemptRegistry::unknown_json of |id|) as the answer for an id a cancel
// tombstoned before any request with it arrived
Json::Value tombstoned_json(Json::Value unknown, const std::string &id) {
    unknown["error"] = "unknown attempt_id '" + id + "': a cancel named it before any request "
                       "with it arrived; a request with it will be refused (409)";
    unknown["tombstone"] = true;
    return unknown;
}

// |unknown| as the 429 of a cancel that could not tombstone the id: |error|, and |reason|
Json::Value not_tombstoned_json(Json::Value unknown, std::string error, const char *reason) {
    unknown["error"] = std::move(error);
    unknown["cancelled"] = false;
    unknown["tombstone"] = false;
    unknown["reason"] = reason;
    return unknown;
}

const char* stop_name(ExternalStop s) {
    switch (s) {
        case ExternalStop::NONE: return nullptr;
        case ExternalStop::CANCELLED: return "cancel";
        case ExternalStop::ATTEMPT_DEADLINE: return "deadline";
        case ExternalStop::CLIENT_GONE: return "client_gone";
    }
    return nullptr;
}

std::string random_instance() {
    std::random_device device;
    std::mt19937_64 rng((static_cast<uint64_t>(device()) << 32) ^ device()
                        ^ static_cast<uint64_t>(std::chrono::system_clock::now()
                                                        .time_since_epoch().count()));
    return fmt::format("{:016x}", rng());
}

} // namespace

bool valid_attempt_id(const std::string &id) {
    if (id.empty() || id.size() > 128)
        return false;
    for (char c : id) {
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')
                     || c == '.' || c == '_' || c == ':' || c == '-';
        if (!ok)
            return false;
    }
    return true;
}

std::string iso_utc(std::chrono::system_clock::time_point t) {
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t.time_since_epoch()).count();
    const std::time_t seconds = static_cast<std::time_t>(ms / 1000);
    std::tm tm {};
    gmtime_r(&seconds, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return fmt::format("{}.{:03d}Z", buf, static_cast<int>(ms % 1000));
}

static std::string iso_of_ms(uint64_t ms) {
    return iso_utc(std::chrono::system_clock::time_point(std::chrono::milliseconds(ms)));
}

bool not_after_passed(const AttemptIds &ids, uint64_t now_ms) {
    return ids.not_after_ms && now_ms > *ids.not_after_ms;
}

Json::Value expired_json(const AttemptIds &ids, uint64_t now_ms,
                         const std::string &server_instance) {
    assert(ids.not_after_ms);
    Json::Value j;
    j["error"] = fmt::format("not_after_ms {} ({}) has passed on this server's clock ({}, {}) "
                             "when the request's handler started: it was not started",
                             *ids.not_after_ms, iso_of_ms(*ids.not_after_ms), now_ms,
                             iso_of_ms(now_ms));
    // told apart from the other 409 (a duplicate id, which carries `attempt`) by its state
    j["state"] = "expired";
    j["server_time_ms"] = uint_json(now_ms);
    put_attempt_ids(&j, ids, false);
    j["server_instance"] = server_instance;
    return j;
}

std::optional<uint64_t> read_not_after_ms(const Json::Value &request, const std::string &where) {
    if (!request.isObject() || !request.isMember("not_after_ms"))
        return std::nullopt;
    // an integer in [0, 2^53 - 1]: a fraction, a sign, a string or a value no JSON reader keeps
    // exactly is refused, never rounded into an instant the ledger did not mean
    const Json::Value &v = request["not_after_ms"];
    if (!v.isIntegral() || (v.isInt64() && v.asInt64() < 0) || v.asUInt64() > kMaxNotAfterMs) {
        throw InvalidRequest(where + ".not_after_ms: expected an integer in [0, "
                             + std::to_string(kMaxNotAfterMs) + "] (Unix epoch, ms)");
    }
    return v.asUInt64();
}

bool instance_mismatch(const AttemptIds &ids, const std::string &instance) {
    return !ids.expect_server_instance.empty() && ids.expect_server_instance != instance;
}

Json::Value instance_mismatch_json(const AttemptIds &ids, const std::string &instance) {
    Json::Value j;
    j["error"] = "expect_server_instance '" + ids.expect_server_instance + "' is not this "
                 "server's instance ('" + instance + "'; the process restarted, or the request "
                 "reached another server): it was not started";
    j["state"] = "instance_mismatch";
    j["expect_server_instance"] = ids.expect_server_instance;
    j["server_instance"] = instance;
    put_attempt_ids(&j, ids, false);
    return j;
}

AttemptIds attempt_ids(const Json::Value &request) {
    AttemptIds ids;
    if (!request.isObject())
        return ids;   // parse_traverse_request names the problem
    ids.not_after_ms = read_not_after_ms(request, "request");
    auto read = [&](const char *field, std::string *out) {
        if (!request.isMember(field))
            return;
        const Json::Value &v = request[field];
        if (!v.isString() || !valid_attempt_id(v.asString())) {
            throw InvalidRequest(std::string("request.") + field + ": expected a string of 1 to "
                                 "128 characters from [A-Za-z0-9._:-]");
        }
        *out = v.asString();
    };
    read("attempt_id", &ids.attempt_id);
    read("budget_id", &ids.budget_id);
    read("locus_id", &ids.locus_id);
    // the server_instance is 16 hex digits; any token is compared as given (another one is
    // refused, 409), so a future format does not turn a mismatch into a 400
    read("expect_server_instance", &ids.expect_server_instance);
    // an echo without the attempt it belongs to would be usage nobody can reconcile
    if (ids.attempt_id.empty() && !ids.budget_id.empty())
        throw InvalidRequest("request.budget_id: given without request.attempt_id");
    if (ids.attempt_id.empty() && !ids.locus_id.empty())
        throw InvalidRequest("request.locus_id: given without request.attempt_id");
    if (ids.attempt_id.empty() && !ids.expect_server_instance.empty())
        throw InvalidRequest("request.expect_server_instance: given without request.attempt_id");
    return ids;
}


// ---------------------------------------------------------------- Attempt

Attempt::Attempt(size_t request_id, std::chrono::system_clock::time_point header_read,
                 const AttemptSettings &settings, std::string server_instance, bool enforced,
                 std::function<bool()> peer_gone)
      : request_id_(request_id), settings_(settings), server_instance_(std::move(server_instance)),
        enforced_(enforced), peer_gone_(std::move(peer_gone)) {
    if (!settings_.clock)
        settings_.clock = [] { return Clock::now(); };
    const auto wall = std::chrono::system_clock::now();
    const Clock::time_point steady = now();
    // the clock starts when the header was read (the HTTP server's content timeout runs from
    // then), which the handler learns only later; never in the future
    if (header_read.time_since_epoch().count() == 0 || header_read > wall)
        header_read = wall;
    received_wall_ = header_read;
    received_ = steady - std::chrono::duration_cast<Clock::duration>(wall - header_read);
    next_peer_check_ = steady;
    // until the request is parsed, the HTTP server's cap alone bounds it
    bound_ms_ = settings_.hard_cap_ms > 0 ? settings_.hard_cap_ms
                                          : std::numeric_limits<double>::infinity();
    walk_until_ms_ = std::max(0.0, bound_ms_ - settings_.allowance_ms / 2);
    capped_ = settings_.hard_cap_ms > 0;
}

Clock::time_point Attempt::now() const { return settings_.clock(); }

double Attempt::ms_since_received(Clock::time_point t) const {
    return std::chrono::duration<double, std::milli>(t - received_).count();
}

double Attempt::elapsed_ms() const { return ms_since_received(now()); }

std::string Attempt::iso(Clock::time_point t) const {
    return iso_utc(received_wall_ + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                                            t - received_));
}

void Attempt::set_bound(size_t seeds, double time_budget_ms, uint64_t memory_budget) {
    std::lock_guard<std::mutex> lock(mutex_);
    memory_budget_ = memory_budget;
    // a non-positive budget means "no extension" (the walk of such a seed ends at once); the
    // allowance covers reading it
    const double t = time_budget_ms > 0 ? time_budget_ms : 0;
    const double raw = static_cast<double>(seeds) * t + settings_.allowance_ms;
    bound_seeds_ = seeds;
    bound_time_budget_ms_ = t;
    capped_ = settings_.hard_cap_ms > 0 && !(raw <= settings_.hard_cap_ms);
    bound_ms_ = capped_ ? settings_.hard_cap_ms : raw;
    bound_set_ = true;
    checked_walk_until_ms_.store(std::numeric_limits<double>::infinity(),
                                 std::memory_order_relaxed);
    update_walk_until_locked();
    seeds_.assign(seeds, SeedUsage());
    seeds_requested_ = seeds;
}

double Attempt::ms_left() const {
    if (!enforced())
        return std::numeric_limits<double>::infinity();
    return walk_until_ms_ - elapsed_ms();
}

double Attempt::ratio_locked() const {
    // the measured ratios replace the configured one, the smaller (more text per account) of
    // this server's and this attempt's own
    if (server_ratio_ > 0 || own_ratio_ > 0)
        return prefer_measured(server_ratio_, own_ratio_);
    return configured_ratio_;
}

double Attempt::reserve_ms() const {
    // the measured rates replace the configured ones: this server's (the slowest of its recent
    // responses) and, for building, this attempt's own seeds, the slower of the two
    double build_mbps = settings_.delivery_build_mbps;
    if (server_.build_mbps > 0 || measured_build_mbps_ > 0)
        build_mbps = prefer_measured(server_.build_mbps, measured_build_mbps_);
    const double compress_mbps = server_.compress_mbps > 0 ? server_.compress_mbps
                                                           : settings_.delivery_compress_mbps;
    // the walked seed's text, estimated from its account: its record coordinates' share at
    // their own bound (kCoordinateAccountPerTextByte), the rest at the ratio in use. Without
    // coordinates the share is 0 and the estimate is the one before the split
    const uint64_t coordinates = std::min(walking_coordinates_, walking_account_);
    const double walking
        = std::ceil(static_cast<double>(walking_account_ - coordinates) / ratio_locked())
        + std::ceil(static_cast<double>(coordinates)
                    / static_cast<double>(kCoordinateAccountPerTextByte));
    // MB/s are bytes per microsecond: x 1000 bytes per ms. With a margin for rates that vary
    // between responses, and the time the walk takes to end once its walk-until passed (it
    // stops at its next poll, then finalises the stopped seed), the configured one or the
    // longest this server measured recently
    const double model = (static_cast<double>(delivered_bytes_) + walking) / (compress_mbps * 1000)
                       + walking / (build_mbps * 1000);
    return kReserveMargin * model + std::max(settings_.delivery_stop_ms, server_.stop_ms);
}

void Attempt::set_measured(const DeliveryMeasurements &measured) {
    std::lock_guard<std::mutex> lock(mutex_);
    server_ = measured;
    if (auto it = server_.account_per_text_byte.find(detail_);
            it != server_.account_per_text_byte.end()) {
        server_ratio_ = it->second;
    }
    update_walk_until_locked();
}

double Attempt::own_build_mbps() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return measured_build_mbps_;
}

double Attempt::own_account_per_text_byte() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return own_ratio_;
}

double Attempt::own_stop_ms() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return own_stop_ms_;
}

void Attempt::update_walk_until_locked() {
    // the reserve keeps time back for a delivery the server bounds: an attempt whose bound is
    // not enforced (the CLI) states the floor
    walk_until_ms_ = std::max(0.0, bound_ms_ - std::max(settings_.allowance_ms / 2,
                                                        enforced() ? reserve_ms() : 0.0));
}

void Attempt::set_delivery_detail(const std::string &detail) {
    std::lock_guard<std::mutex> lock(mutex_);
    detail_ = detail;
    configured_ratio_ = detail == "graphlet" ? settings_.account_per_text_byte_graphlet
                                             : settings_.account_per_text_byte_json;
    auto it = server_.account_per_text_byte.find(detail);
    server_ratio_ = it != server_.account_per_text_byte.end() ? it->second : 0;
    update_walk_until_locked();
}

void Attempt::progress(uint64_t account, uint64_t coordinates) {
    std::lock_guard<std::mutex> lock(mutex_);
    walking_account_ = account;
    walking_coordinates_ = coordinates;
    update_walk_until_locked();
}

void Attempt::note_delivered(uint64_t text_bytes, double build_seconds, uint64_t account,
                             uint64_t coordinate_account, uint64_t coordinate_text) {
    std::lock_guard<std::mutex> lock(mutex_);
    delivered_bytes_ += text_bytes;
    walking_account_ = 0;
    walking_coordinates_ = 0;
    // the build rate over the whole text: the coordinates' is built like the rest
    if (text_bytes >= kMeasuredTextBytes && build_seconds > 0) {
        const double mbps = static_cast<double>(text_bytes) / build_seconds / 1e6;
        measured_build_mbps_ = measured_build_mbps_ > 0 ? std::min(measured_build_mbps_, mbps)
                                                        : mbps;
    }
    // The ratio of the rest of the output: the coordinate share's account and its exact text
    // both left out, so that the sample is the one the same walk without coordinates gives,
    // and measured only on a rest of at least kMeasuredTextBytes (a seed whose text is
    // mostly coordinates measures no ratio, like a small seed). Without coordinates both are
    // 0 and the sample is the whole output's
    if (account > coordinate_account && text_bytes > coordinate_text
            && text_bytes - coordinate_text >= kMeasuredTextBytes) {
        const double ratio = static_cast<double>(account - coordinate_account)
                           / static_cast<double>(text_bytes - coordinate_text);
        own_ratio_ = own_ratio_ > 0 ? std::min(own_ratio_, ratio) : ratio;
    }
    update_walk_until_locked();
}

void Attempt::note_max_uninterruptible_ms(double ms) {
    std::lock_guard<std::mutex> lock(mutex_);
    max_uninterruptible_ms_ = std::max(max_uninterruptible_ms_, ms);
}

double Attempt::max_uninterruptible_ms() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return max_uninterruptible_ms_;
}

void Attempt::note_delivery_gap_ms(double ms) {
    std::lock_guard<std::mutex> lock(mutex_);
    max_delivery_gap_ms_ = std::max(max_delivery_gap_ms_, ms);
}

double Attempt::max_delivery_gap_ms() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return max_delivery_gap_ms_;
}

bool Attempt::request_stop(ExternalStop reason) {
    uint8_t none = 0;
    if (!stop_.compare_exchange_strong(none, static_cast<uint8_t>(reason)))
        return false;
    std::lock_guard<std::mutex> lock(mutex_);
    stop_requested_at_ = now();
    return true;
}

bool Attempt::cancel() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        // checked with finish() under one lock: a finished attempt is not cancelled (404)
        if (finished_at_)
            return false;
        if (!cancel_requested_at_)
            cancel_requested_at_ = now();
    }
    request_stop(ExternalStop::CANCELLED);
    return true;
}

ExternalStop Attempt::poll(bool force) {
    uint8_t stop = stop_.load(std::memory_order_acquire);
    // the clock only every poll_stride-th poll: a poll precedes every head, and reading the
    // clock there cost the walk about one percent
    if (!stop && (force || ++polls_ >= settings_.poll_stride)) {
        polls_ = 0;
        if (enforced() && elapsed_ms() >= walk_until_ms_) {
            if (request_stop(ExternalStop::ATTEMPT_DEADLINE)) {
                // where the walk was stopped: the walk-until in force now (it moves with the
                // delivery reserve, up again once a seed is delivered)
                std::lock_guard<std::mutex> lock(mutex_);
                tripped_walk_until_ms_ = walk_until_ms_;
            }
        } else if (enforced()
                   && walk_until_ms_ < checked_walk_until_ms_.load(std::memory_order_relaxed)) {
            // the walk was bounded by it here (usage.bound.walk_until_ms)
            checked_walk_until_ms_.store(walk_until_ms_, std::memory_order_relaxed);
        }
        if (!stop_.load(std::memory_order_acquire) && peer_gone_) {
            const Clock::time_point t = now();
            if (t >= next_peer_check_) {
                next_peer_check_ = t + std::chrono::milliseconds(settings_.client_check_ms);
                if (peer_gone_())
                    request_stop(ExternalStop::CLIENT_GONE);
            }
        }
        stop = stop_.load(std::memory_order_acquire);
    }
    if (stop && !stop_seen_) {
        // the first poll that hands the stop to the walk: when it stopped (the walker acts on
        // a stop at the checkpoint that polled it)
        std::lock_guard<std::mutex> lock(mutex_);
        stop_seen_ = true;
        if (!stopped_at_)
            stopped_at_ = now();
    }
    return static_cast<ExternalStop>(stop);
}

void Attempt::check_delivery() {
    if (!client_gone_ && peer_gone_) {
        const Clock::time_point t = now();
        if (t >= next_peer_check_) {
            next_peer_check_ = t + std::chrono::milliseconds(settings_.client_check_ms);
            client_gone_ = peer_gone_();
        }
    }
    // a client gone after a cancel or the deadline stopped the walk is still gone: nothing can
    // be delivered (the stop keeps its first reason)
    if (client_gone_
            || stop_.load(std::memory_order_acquire) == static_cast<uint8_t>(ExternalStop::CLIENT_GONE)) {
        client_gone_ = true;
        request_stop(ExternalStop::CLIENT_GONE);
        throw graph::traversal::AttemptAborted("the client closed its connection while the "
                                               "response was built or written");
    }
    if (enforced() && elapsed_ms() >= bound_ms_) {
        request_stop(ExternalStop::ATTEMPT_DEADLINE);
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_seen_ = true;
            if (!stopped_at_)
                stopped_at_ = now();
        }
        // the bound as enforced: capped at the content timeout less one second
        throw AttemptAtBound(fmt::format(
                "the attempt reached the duration bound the server enforces for it ({} ms: the "
                "seeds' time budgets plus the server's allowance, at most hard_cap_ms, the "
                "content timeout less one second; see usage.bound) while its response was built "
                "or written: nothing of it is delivered",
                static_cast<uint64_t>(std::ceil(bound_ms_))));
    }
}

void Attempt::walk_ended() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!stopped_at_)
        stopped_at_ = now();
}

const char* Attempt::walk_reason() const {
    if (!stop_seen_)
        return "completed";
    switch (static_cast<ExternalStop>(stop_.load(std::memory_order_acquire))) {
        case ExternalStop::CANCELLED: return "cancelled";
        case ExternalStop::ATTEMPT_DEADLINE: return "deadline";
        case ExternalStop::CLIENT_GONE: return "client_gone";
        case ExternalStop::NONE: break;
    }
    return "completed";
}

void Attempt::seed_started(size_t index) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (index >= seeds_.size())
        seeds_.resize(index + 1);
    seeds_requested_ = std::max(seeds_requested_, seeds_.size());
    seeds_[index].started = true;
    ++seeds_started_;
}

void Attempt::seed_walked(size_t index, const SeedUsage &usage) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (index >= seeds_.size())
        seeds_.resize(index + 1);
    seeds_requested_ = std::max(seeds_requested_, seeds_.size());
    // once per seed: a seed abandoned while its result was built, or refused after its walk,
    // was already counted. A walk the client's departure cut did not finish: it is counted as
    // abandoned, not as finished
    const bool first = !seeds_[index].walked;
    if (first) {
        if (usage.outcome == "abandoned") {
            seeds_abandoned_++;
        } else {
            seeds_walked_++;
        }
    }
    seeds_[index] = usage;
    seeds_[index].started = true;
    seeds_[index].walked = true;
    // A walk that ended after the walk-until (stopped at a poll after it, or ended between two
    // of them) measures what the reserve must keep beyond building and compressing — at the
    // walk's end only, the first call for the seed. A second call comes after the walk, while
    // its result was built (the client left: abandoned; a writer or the seed's refusal: failed),
    // and the time since the walk-until then includes building, which the reserve prices on its
    // own: measured as stop latency it would become the server's measured_stop_ms (the longest
    // of its last rate_window attempts), and every later attempt whose bound is no longer than
    // that build time would walk nothing (10.5 s of building recorded after a 400 ms stop, for
    // example)
    if (first && enforced() && bound_set_) {
        const double late = elapsed_ms() - tripped_walk_until_ms_.value_or(walk_until_ms_);
        if (late > 0)
            own_stop_ms_ = std::max(own_stop_ms_, late);
    }
}

void Attempt::seed_not_started(size_t index, const SeedUsage &usage) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (index >= seeds_.size())
        seeds_.resize(index + 1);
    seeds_requested_ = std::max(seeds_requested_, seeds_.size());
    seeds_[index] = usage;
    seeds_[index].outcome = "not_started";
    seeds_[index].started = false;
    seeds_[index].walked = false;
}

void Attempt::seed_delivered(size_t index, const std::string &outcome, double elapsed_ms) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (index >= seeds_.size())
        return;
    seeds_[index].outcome = outcome;
    seeds_[index].elapsed_ms = elapsed_ms;
}

Json::Value Attempt::bound_json() const {
    Json::Value b;
    b["seeds"] = uint_json(bound_seeds_);
    b["time_budget_ms"] = bound_set_ ? Json::Value(bound_time_budget_ms_) : Json::Value();
    b["allowance_ms"] = ceil_ms_json(settings_.allowance_ms);
    // When the walk-until stopped the walk, the walk-until in force then: where the seeds
    // stopped being walked (the walk stopped at its first poll that read the clock after it,
    // usage.stopped_at). Otherwise the lowest walk-until the walk's clock-reading polls checked
    // (one in poll_stride, every forced one): the floor (bound - allowance /
    // 2) unless the delivery reserve moved it before a check. A value no poll read bounds no
    // walk: the lowest computed would read 14905 ms for walks stopped near 16000, and 20174 ms
    // (from the last seed's text, written after its walk) for a walk its own 30 s budget ended
    const double checked = checked_walk_until_ms_.load(std::memory_order_relaxed);
    b["walk_until_ms"] = ceil_ms_json(tripped_walk_until_ms_ ? *tripped_walk_until_ms_
                               : std::isfinite(checked) ? checked : walk_until_ms_);
    b["capped_by"] = capped_ ? Json::Value("content_timeout") : Json::Value();
    b["enforced"] = enforced();
    return b;
}

Json::Value Attempt::seeds_json() const {
    Json::Value seeds;
    seeds["requested"] = uint_json(seeds_requested_);
    seeds["started"] = uint_json(seeds_started_);
    // a seed whose walk ended (complete, partial, failed), delivered or not
    seeds["finished"] = uint_json(seeds_walked_);
    // a seed whose walk the client's departure cut: neither finished nor delivered
    seeds["abandoned"] = uint_json(seeds_abandoned_);
    return seeds;
}

Json::Value Attempt::usage_json(const std::string &reason, bool per_seed) const {
    std::lock_guard<std::mutex> lock(mutex_);
    Json::Value u;
    // not_after_ms echoed only when the request gave it, so that every other usage keeps its
    // keys
    put_attempt_ids(&u, ids_, true);
    u["server_instance"] = server_instance_;
    u["reason"] = reason;
    u["received_at"] = iso(received_);
    u["stopped_at"] = stopped_at_ ? Json::Value(iso(*stopped_at_)) : Json::Value();
    u["elapsed_ms"] = ceil_ms_json(ms_since_received(finished_at_ ? *finished_at_ : now()));
    u["bound_ms"] = ceil_ms_json(bound_ms_);
    u["bound"] = bound_json();
    u["seeds"] = seeds_json();
    // The request: work adds up; memory is held at once by the seed being walked and the
    // results of the seeds before it, which stay until the response is written — a walked
    // result its account at the end, a failed or never started one its priced echo
    // (final_bytes). Under a memory budget a seed's walk holds at most the budget and its soft
    // excess (as observed): what the request held at once is bounded by that, the results
    // before it included
    uint64_t work = 0, held = 0, peak = 0, soft = 0, bound = 0;
    Json::Value list(Json::arrayValue);
    for (size_t i = 0; i < seeds_.size(); ++i) {
        const SeedUsage &s = seeds_[i];
        const graph::traversal::AttemptMeter &m = s.meter;
        work += m.work_units;
        if (s.walked) {
            peak = std::max(peak, held + m.memory_peak);
            bound = std::max(bound, held + memory_budget_ + m.soft_excess);
        }
        soft = std::max(soft, m.soft_excess);
        // an abandoned walk has no result (nothing is written)
        if (s.outcome != "abandoned")
            held += m.memory_final;
        peak = std::max(peak, held);
        bound = std::max(bound, held);
        if (!per_seed)
            continue;
        Json::Value e;
        e["index"] = uint_json(i);
        e["outcome"] = s.outcome;
        e["stopped_by"] = s.stopped_by.empty() ? Json::Value() : Json::Value(s.stopped_by);
        e["work_units"] = uint_json(m.work_units);
        e["work_seed"] = uint_json(m.work_seed);
        e["peak_admitted_bytes"] = uint_json(m.memory_peak);
        e["final_bytes"] = uint_json(m.memory_final);
        e["soft_excess_bytes"] = memory_budget_ ? uint_json(m.soft_excess) : Json::Value();
        e["refused_bytes"] = s.refused_bytes ? uint_json(*s.refused_bytes) : Json::Value();
        e["elapsed_ms"] = ceil_ms_json(s.elapsed_ms);
        list.append(std::move(e));
    }
    u["work_units"] = uint_json(work);
    // the longest read or head piece of the request's walks: how late a stop could be seen
    // in them (an observation of this attempt, not a bound; see deadline_check)
    u["observed_max_uninterruptible_ms"] = ceil_ms_json(max_uninterruptible_ms_);
    Json::Value memory;
    memory["peak_admitted_bytes"] = uint_json(peak);
    // the soft part is observed only under a memory budget, and per seed (each seed's excess
    // over its own budget): the request states the largest
    memory["soft_excess_bytes"] = memory_budget_ ? uint_json(soft) : Json::Value();
    // Without a memory budget the account leaves out the label caches, the lookahead and a
    // level's decoded rows, so nothing bounds what was held (a walk whose account states
    // 0.79 MB can raise the RSS by 111 MB): null, never a number that reads as a bound
    memory["held_bound_bytes"] = memory_budget_ ? uint_json(bound) : Json::Value();
    u["memory"] = std::move(memory);
    if (per_seed)
        u["per_seed"] = std::move(list);
    return u;
}

Json::Value Attempt::state_json() const {
    std::lock_guard<std::mutex> lock(mutex_);
    Json::Value j;
    put_attempt_ids(&j, ids_, true);
    j["server_instance"] = server_instance_;
    if (tombstone_) {
        // a cancel that came first: no request with this id ran here, and none will
        j["state"] = "unknown";
        j["tombstone"] = true;
        j["cancel_requested_at"] = cancel_requested_at_ ? Json::Value(iso(*cancel_requested_at_))
                                                        : Json::Value();
        return j;
    }
    const auto stop = static_cast<ExternalStop>(stop_.load(std::memory_order_acquire));
    j["state"] = finished_at_ ? "finished" : stop != ExternalStop::NONE ? "stopping" : "running";
    j["reason"] = finished_at_ ? Json::Value(reason_) : Json::Value();
    j["stop_requested_by"] = stop_name(stop) ? Json::Value(stop_name(stop)) : Json::Value();
    auto at = [&](const std::optional<Clock::time_point> &t) {
        return t ? Json::Value(iso(*t)) : Json::Value();
    };
    j["received_at"] = iso(received_);
    j["cancel_requested_at"] = at(cancel_requested_at_);
    j["stop_requested_at"] = at(stop_requested_at_);
    j["stopped_at"] = at(stopped_at_);
    j["finished_at"] = at(finished_at_);
    j["deadline_at"] = std::isfinite(bound_ms_)
        ? Json::Value(iso(received_ + std::chrono::duration_cast<Clock::duration>(
                                          std::chrono::duration<double, std::milli>(bound_ms_))))
        : Json::Value();
    j["bound_ms"] = ceil_ms_json(bound_ms_);
    j["elapsed_ms"] = ceil_ms_json(ms_since_received(finished_at_ ? *finished_at_ : now()));
    j["seeds"] = seeds_json();
    Json::Value response;
    response["written"] = finished_at_ ? Json::Value(status_ != 0) : Json::Value();
    response["status"] = status_ ? Json::Value(status_) : Json::Value();
    response["bytes"] = bytes_ ? uint_json(*bytes_) : Json::Value();
    j["response"] = std::move(response);
    j["usage"] = finished_at_ ? final_usage_ : Json::Value();
    return j;
}

std::string Attempt::stop_summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::string out = stopped_at_ ? "walk stopped at " + iso(*stopped_at_) : "walk not stopped";
    const auto stop = static_cast<ExternalStop>(stop_.load(std::memory_order_acquire));
    if (stopped_at_ && stop_requested_at_ && stop != ExternalStop::NONE) {
        const auto asked = cancel_requested_at_ ? *cancel_requested_at_ : *stop_requested_at_;
        out += fmt::format(" ({:.0f} ms after the stop was requested by {}{})",
                           std::max(0.0, std::chrono::duration<double, std::milli>(
                                                 *stopped_at_ - asked).count()),
                           stop_name(stop), stop_seen_ ? "" : "; the walk had already ended");
    }
    out += fmt::format(", {} of {} seed(s) walked{}, bound {} ms", seeds_walked_, seeds_requested_,
                       seeds_abandoned_ ? fmt::format(" ({} abandoned in its walk)", seeds_abandoned_)
                                        : std::string(),
                       std::isfinite(bound_ms_) ? fmt::format("{:.0f}", std::ceil(bound_ms_))
                                                : std::string("none"));
    return out;
}

void Attempt::finish(const std::string &reason, int status, std::optional<size_t> bytes) {
    // the usage kept for GET: as the response stated it, without per_seed
    Json::Value usage = managed() && !tombstone_ ? usage_json(reason, false) : Json::Value();
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (finished_at_)
            return;
        finished_at_ = now();
        if (!stopped_at_)
            stopped_at_ = finished_at_;
        reason_ = reason;
        status_ = status;
        bytes_ = bytes;
        if (!usage.isNull()) {
            usage["elapsed_ms"] = ceil_ms_json(ms_since_received(*finished_at_));
            usage["stopped_at"] = iso(*stopped_at_);
        }
        final_usage_ = std::move(usage);
        // a finished attempt is kept for its state only (retention_count of them): what it
        // held for its walk goes — the per-seed records (the usage kept has the totals) and the
        // client check, which holds the request and its body
        std::vector<SeedUsage>().swap(seeds_);
    }
    // the handler's thread (on_written), the only one that calls it
    peer_gone_ = nullptr;
    finished_cv_.notify_all();
}

bool Attempt::finished() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return finished_at_.has_value();
}

Clock::time_point Attempt::finished_at() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return finished_at_ ? *finished_at_ : now();
}

bool Attempt::wait_finished(uint64_t ms) const {
    std::unique_lock<std::mutex> lock(mutex_);
    return finished_cv_.wait_for(lock, std::chrono::milliseconds(ms),
                                 [&] { return finished_at_.has_value(); });
}

std::shared_ptr<Attempt> Attempt::make_tombstone(const std::string &id,
                                                 const AttemptSettings &settings,
                                                 std::string server_instance) {
    auto a = std::make_shared<Attempt>(0, std::chrono::system_clock::now(), settings,
                                       std::move(server_instance), false);
    AttemptIds ids;
    ids.attempt_id = id;
    a->set_ids(std::move(ids));
    a->tombstone_ = true;
    a->cancel_requested_at_ = a->now();
    a->finished_at_ = a->cancel_requested_at_;
    a->reason_ = "cancelled";
    return a;
}


// ---------------------------------------------------------------- AttemptRegistry

AttemptRegistry::AttemptRegistry(AttemptSettings settings)
      : settings_(std::move(settings)), instance_(random_instance()) {
    if (!settings_.clock)
        settings_.clock = [] { return Clock::now(); };
}

Clock::time_point AttemptRegistry::now() const { return settings_.clock(); }

std::string AttemptRegistry::retention_text() const {
    if (!settings_.retention_s) {
        return fmt::format("attempts are not kept after they finish (retention 0 s); no cancel "
                           "of an unknown id is tombstoned (refused, 429: nothing is promised "
                           "for it)");
    }
    return fmt::format("attempts are kept {} s after they finish, at most the last {} (oldest "
                       "dropped first), and one sent with not_after_ms until that + {} ms, at "
                       "most {} s after it finished or after its latest refused copy (held "
                       "among the tombstones); a cancel of an "
                       "unknown id tombstones it at least {} s, and until the not_after_ms it "
                       "names + {} ms, at most {} s, while fewer than {} tombstones and held "
                       "attempts are held (else the cancel is refused, 429)",
                       settings_.retention_s, settings_.retention_count, settings_.clock_skew_ms,
                       tombstone_cap_ms() / 1000, settings_.retention_s,
                       settings_.clock_skew_ms, tombstone_cap_ms() / 1000,
                       settings_.retention_count);
}

uint64_t AttemptRegistry::tombstone_cap_ms() const {
    return std::max(settings_.tombstone_max_s, settings_.retention_s) * 1000;
}

void AttemptRegistry::note_uninterruptible(double ms) {
    if (!(ms > 0))
        return;
    const uint64_t us = static_cast<uint64_t>(std::ceil(std::min(ms, 1e12) * 1000));
    uint64_t seen = max_uninterruptible_us_.load(std::memory_order_relaxed);
    while (us > seen && !max_uninterruptible_us_.compare_exchange_weak(seen, us)) {}
}

uint64_t AttemptRegistry::observed_max_uninterruptible_ms() const {
    return (max_uninterruptible_us_.load(std::memory_order_relaxed) + 999) / 1000;
}

void AttemptRegistry::note_build_rate(double mbps) {
    if (!(mbps > 0) || !std::isfinite(mbps))
        return;
    std::lock_guard<std::mutex> lock(rates_mutex_);
    push_window(&build_rates_, mbps);
}

void AttemptRegistry::note_compress_rate(double mbps) {
    if (!(mbps > 0) || !std::isfinite(mbps))
        return;
    std::lock_guard<std::mutex> lock(rates_mutex_);
    push_window(&compress_rates_, mbps);
}

void AttemptRegistry::note_account_per_text_byte(const std::string &detail, double ratio) {
    if (!(ratio > 0) || !std::isfinite(ratio) || detail.empty())
        return;
    std::lock_guard<std::mutex> lock(rates_mutex_);
    push_window(&ratios_[detail], ratio);
}

void AttemptRegistry::note_stop_latency(double ms) {
    if (!(ms > 0) || !std::isfinite(ms))
        return;
    std::lock_guard<std::mutex> lock(rates_mutex_);
    push_window(&stop_latencies_, ms);
}

DeliveryMeasurements AttemptRegistry::measured() const {
    std::lock_guard<std::mutex> lock(rates_mutex_);
    auto least = [](const std::deque<double> &values) {
        return values.empty() ? 0.0 : *std::min_element(values.begin(), values.end());
    };
    DeliveryMeasurements m;
    m.build_mbps = least(build_rates_);
    m.compress_mbps = least(compress_rates_);
    m.stop_ms = stop_latencies_.empty()
        ? 0.0 : *std::max_element(stop_latencies_.begin(), stop_latencies_.end());
    for (const auto &[detail, ratios] : ratios_) {
        if (!ratios.empty())
            m.account_per_text_byte[detail] = least(ratios);
    }
    return m;
}

bool AttemptRegistry::hold_live_locked(const Attempt &a, Clock::time_point t,
                                       uint64_t wall) const {
    return a.tomb_steady_until_ > t || a.tomb_wall_until_ms_ >= wall;
}

void AttemptRegistry::expire_locked() {
    const Clock::time_point t = now();
    const auto keep = std::chrono::seconds(settings_.retention_s);
    uint64_t wall = 0;
    bool wall_read = false;
    auto wall_now = [&]() {
        if (!wall_read) {
            wall = now_ms();
            wall_read = true;
        }
        return wall;
    };
    while (!retained_.empty()
            && (retained_.size() > settings_.retention_count
                || t - retained_.front()->finished_at() >= keep)) {
        const std::shared_ptr<Attempt> old = retained_.front();
        retained_.pop_front();
        auto it = attempts_.find(old->ids().attempt_id);
        if (it == attempts_.end() || it->second != old)
            continue;
        if (hold_live_locked(*old, t, wall_now())) {
            // Its retention ended, by age or by count, while a copy of the request could still
            // be admitted: the id stays refused, among the tombstones, until its hold ends (a
            // ledger releases on the finished state, and a replay after retention_s, or after
            // retention_count later finishes, would run again). Never dropped early: it ran,
            // so it cannot be refused like a cancel
            old->held_ = true;
            tomb_expiry_.emplace(old->tomb_steady_until_, old->ids().attempt_id);
            continue;
        }
        attempts_.erase(it);
    }
    // by its hold only: a tombstone (or a held finished attempt) promises that its id does not
    // run here while it is live, and it is live while either clock says so
    while (!tomb_expiry_.empty() && tomb_expiry_.begin()->first <= t) {
        const auto [key, id] = *tomb_expiry_.begin();
        tomb_expiry_.erase(tomb_expiry_.begin());
        auto it = attempts_.find(id);
        if (it == attempts_.end())
            continue;
        Attempt &held = *it->second;
        if (!(held.tombstone() || held.held_) || held.tomb_steady_until_ != key)
            continue;   // not this entry's (hold_locked re-keys its own)
        if (held.tomb_wall_until_ms_ >= wall_now()) {
            // the wall clock stepped back since the hold was set: kept until it reads past the
            // expiry too (re-keyed strictly later, so the loop ends)
            held.tomb_steady_until_
                = t + std::chrono::milliseconds(held.tomb_wall_until_ms_ - wall + 1);
            tomb_expiry_.emplace(held.tomb_steady_until_, id);
            continue;
        }
        attempts_.erase(it);
    }
}

void AttemptRegistry::hold_locked(Attempt &held, std::optional<uint64_t> not_after_ms,
                                  bool at_least_retention) {
    const uint64_t wall = now_ms();
    uint64_t hold_ms = at_least_retention ? settings_.retention_s * 1000 : 0;
    if (not_after_ms) {
        held.tomb_not_after_ms_ = std::max(held.tomb_not_after_ms_.value_or(0), *not_after_ms);
        // covered once this server's clock passed not_after_ms + the skew a ledger adds (both
        // at most 2^53 - 1: no overflow)
        const uint64_t target = *not_after_ms + settings_.clock_skew_ms;
        hold_ms = std::max(hold_ms, std::min(target > wall ? target - wall : 0,
                                             tombstone_cap_ms()));
    } else if (!at_least_retention) {
        return;   // a finished attempt's hold is that of a not_after_ms alone
    }
    // live through wall + hold_ms inclusive; the steady hold one ms longer to match
    held.tomb_wall_until_ms_ = std::max(held.tomb_wall_until_ms_, wall + hold_ms);
    const Clock::time_point until = now() + std::chrono::milliseconds(hold_ms + 1);
    if (until > held.tomb_steady_until_) {
        const std::string &id = held.ids().attempt_id;
        const bool keyed = held.tombstone() || held.held_;
        if (keyed)
            tomb_expiry_.erase({ held.tomb_steady_until_, id });
        held.tomb_steady_until_ = until;
        if (keyed)
            tomb_expiry_.emplace(until, id);
    }
}

void AttemptRegistry::add_suppression_locked(Json::Value *j, const Attempt &tomb,
                                             std::optional<uint64_t> against,
                                             bool fallback) const {
    if (!against && fallback)
        against = tomb.tomb_not_after_ms_;
    // the hold's last instant on the wall clock (inclusive): the tombstone is live whenever
    // this server's clock reads at most it — the steady hold lasts as long from when it was
    // set, so a forward step of the clock does not end it earlier — and once it is gone the
    // clock has read later
    (*j)["suppressed_until_ms"] = uint_json(tomb.tomb_wall_until_ms_);
    if (against)
        (*j)["not_after_ms"] = uint_json(*against);
    const bool covers = against
                     && tomb.tomb_wall_until_ms_ >= *against + settings_.clock_skew_ms;
    (*j)["covers_admission"] = covers;
    if (!covers) {
        // no not_after_ms: nothing bounds when a copy of the request could still be admitted;
        // beyond the cap: the tombstone expires while one could
        (*j)["covers_admission_reason"] = against ? "beyond_tombstone_max" : "no_not_after_ms";
    }
}

Json::Value AttemptRegistry::unknown_json(const std::string &id) const {
    Json::Value j;
    j["error"] = "unknown attempt_id '" + id + "': no attempt with it is running on this server "
                 "and none finished within the retention period (" + retention_text() + ")";
    j["attempt_id"] = id;
    j["state"] = "unknown";
    j["server_instance"] = instance_;
    return j;
}

uint64_t AttemptRegistry::now_ms() const {
    const auto t = settings_.wall_clock ? settings_.wall_clock() : std::chrono::system_clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t.time_since_epoch()).count();
    return ms > 0 ? static_cast<uint64_t>(ms) : 0;
}

std::optional<AttemptRegistry::StartRefusal>
AttemptRegistry::start(const std::shared_ptr<Attempt> &attempt) {
    std::lock_guard<std::mutex> lock(mutex_);
    expire_locked();
    // meant for another process (this one restarted, or the request reached another server):
    // its tombstones, if any, are not here, so nothing here may run it
    if (instance_mismatch(attempt->ids(), instance_))
        return StartRefusal { false, instance_mismatch_json(attempt->ids(), instance_), true };
    const std::string &id = attempt->ids().attempt_id;
    auto it = attempts_.find(id);
    if (it != attempts_.end()) {
        Json::Value body = it->second->state_json();
        if (it->second->tombstone()) {
            // A copy of a cancelled request arrived (a proxy's replay carries the same
            // not_after_ms): the tombstone is extended to cover its admission, so that every
            // later copy is refused too while it could still be admitted, and the refusal
            // states the suppression judged against this request's not_after_ms
            hold_locked(*it->second, attempt->ids().not_after_ms);
            // only this request's own not_after_ms decides whether a copy of THIS request is
            // covered, never a cancel's: a copy without one would otherwise be told
            // covers_admission and run once re-sent after the tombstone expired
            add_suppression_locked(&body, *it->second, attempt->ids().not_after_ms, false);
        } else if (settings_.retention_s && attempt->ids().not_after_ms) {
            // a copy of a running or finished request (a replay, or another request with its
            // id): the id's hold covers the copy's not_after_ms too (a running attempt's is
            // applied again from its finish)
            hold_locked(*it->second, attempt->ids().not_after_ms, false);
        }
        return StartRefusal { false, std::move(body) };
    }
    // not run, and not kept: a later copy of the request is expired as well, so nothing has
    // to remember it (GET /traverse/attempt answers 404 for it)
    const uint64_t now = now_ms();
    if (not_after_passed(attempt->ids(), now))
        return StartRefusal { true, expired_json(attempt->ids(), now, instance_) };
    attempts_.emplace(id, attempt);
    return std::nullopt;
}

Json::Value AttemptRegistry::capabilities_json() const {
    Json::Value att;
    Json::Value fields(Json::arrayValue);
    for (const char *f : kAttemptFields) {
        fields.append(f);
    }
    att["fields"] = std::move(fields);
    att["id_pattern"] = "^[A-Za-z0-9._:-]{1,128}$";
    att["cancel"] = "POST /traverse/cancel";
    Json::Value cancel_fields(Json::arrayValue);
    for (const char *f : { "attempt_id", "wait_ms", "not_after_ms" }) {
        cancel_fields.append(f);
    }
    att["cancel_fields"] = std::move(cancel_fields);
    att["state"] = "GET /traverse/attempt/{attempt_id}";
    att["server_instance"] = instance_;
    att["retention_s"] = uint_json(settings_.retention_s);
    att["retention_count"] = uint_json(settings_.retention_count);
    // the cap of a tombstone's hold, as applied: max(tombstone_max_s, retention_s)
    att["tombstone_max_s"] = uint_json(tombstone_cap_ms() / 1000);
    // integers (ms), as usage.bound states them: a ledger compares them with its own
    att["allowance_ms"] = ceil_ms_json(settings_.allowance_ms);
    att["hard_cap_ms"] = ceil_ms_json(settings_.hard_cap_ms);
    att["content_timeout_s"] = uint_json(settings_.content_timeout_s);
    att["client_check_ms"] = uint_json(settings_.client_check_ms);
    att["clock_skew_allowance_ms"] = uint_json(settings_.clock_skew_ms);
    // what the bound is and how it is enforced (the cap, the walk-until stated, what runs past
    // it up to the next delivery check): the SPEC's rule, the numbers in the fields above and
    // in deadline_check (poll_stride)
    att["bound"] = kBoundRule;
    // the time kept back from the walk to build and compress what was walked
    Json::Value reserve;
    reserve["compress_mbps"] = settings_.delivery_compress_mbps;
    reserve["build_mbps"] = settings_.delivery_build_mbps;
    Json::Value per;
    per["json"] = settings_.account_per_text_byte_json;
    per["graphlet"] = settings_.account_per_text_byte_graphlet;
    reserve["account_per_text_byte"] = std::move(per);
    // the fixed bound the record coordinates' share is estimated with (feature level 6): a
    // number, so that a ledger can reproduce the reserve
    reserve["coordinate_account_per_text_byte"] = uint_json(kCoordinateAccountPerTextByte);
    reserve["measured_text_bytes"] = uint_json(kMeasuredTextBytes);
    const DeliveryMeasurements m = measured();
    reserve["measured_build_mbps"] = m.build_mbps > 0 ? Json::Value(m.build_mbps) : Json::Value();
    reserve["measured_compress_mbps"] = m.compress_mbps > 0 ? Json::Value(m.compress_mbps)
                                                            : Json::Value();
    Json::Value ratios;
    for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
        auto it = m.account_per_text_byte.find(detail);
        ratios[detail] = it != m.account_per_text_byte.end() ? Json::Value(it->second)
                                                             : Json::Value();
    }
    reserve["measured_account_per_text_byte"] = std::move(ratios);
    reserve["rate_window"] = uint_json(kRateWindow);
    reserve["margin"] = kReserveMargin;
    reserve["stop_ms"] = ceil_ms_json(settings_.delivery_stop_ms);
    // where the configured starting estimates come from: the SPEC's table of the measurements
    reserve["calibration"] = kCalibrationRule;
    reserve["measured_stop_ms"] = m.stop_ms > 0 ? ceil_ms_json(m.stop_ms) : Json::Value();
    // the reserve's formula and when a response is still not delivered: the SPEC's rule, its
    // numbers the fields of this block
    reserve["rule"] = kDeliveryReserveRule;
    att["delivery_reserve"] = std::move(reserve);
    // not_after_ms and expect_server_instance (their checks, their 409s and the order of the
    // refusals), the tombstone of a cancel and what it suppresses (retention_s 0: none), and
    // when a ledger may release an attempt's capacity: the SPEC's rules, by reference
    att["not_after"] = kNotAfterRule;
    att["instance"] = kInstanceRule;
    att["suppression"] = kSuppressionRule;
    att["release_rule"] = kReleaseRule;
    return att;
}

std::pair<int, Json::Value> AttemptRegistry::cancel(const std::string &id, uint64_t wait_ms,
                                                    std::optional<uint64_t> not_after_ms) {
    std::shared_ptr<Attempt> attempt;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        expire_locked();
        auto it = attempts_.find(id);
        if (it == attempts_.end()) {
            if (!settings_.retention_s) {
                // Retention 0: this server suppresses nothing, and says so (a tombstone held
                // 0 s would expire at once); a retry gives the same answer
                return { 429, not_tombstoned_json(unknown_json(id),
                    "unknown attempt_id '" + id + "', and it was NOT tombstoned: this server "
                    "keeps no tombstones (--traverse-attempt-retention-s 0), so a request with it "
                    "that arrives later may still run here", "no_suppression") };
            }
            if (tomb_expiry_.size() >= settings_.retention_count) {
                // No tombstone can be kept for its whole hold without evicting another one
                // early, which would break that one's promise: nothing is promised for this
                // id, said so, and the cancel can be retried
                return { 429, not_tombstoned_json(unknown_json(id),
                    "unknown attempt_id '" + id + "', and it was NOT tombstoned: "
                    + std::to_string(tomb_expiry_.size()) + " cancels of unknown ids are held "
                    "already (" + retention_text() + "), so a request with it that arrives later "
                    "may still run here; retry the cancel later", "tombstones_full") };
            }
            // A cancel can overtake its request (still queued, or on the wire, its body half
            // uploaded): the id is tombstoned, so a request that arrives later is refused
            // (409) and never runs here — for retention_s, and until the not_after_ms the
            // cancel names + the clock skew allowance, after which the request's own
            // not_after check refuses it. Only then (covers_admission) does a 404 from this
            // server_instance mean that the attempt will not run on it at all (with retention
            // 1 s alone, a half-uploaded request completing at 1.18 s would run)
            auto tomb = Attempt::make_tombstone(id, settings_, instance_);
            attempts_.emplace(id, tomb);
            hold_locked(*tomb, not_after_ms);
            Json::Value j = unknown_json(id);
            j["cancelled"] = false;
            j["tombstone"] = true;
            add_suppression_locked(&j, *tomb, not_after_ms, true);
            return { 404, j };
        }
        attempt = it->second;
        if (attempt->tombstone()) {
            // repeated: never shortened, extended to this cancel's not_after_ms, and at least
            // retention_s from now
            hold_locked(*attempt, not_after_ms);
            Json::Value j = tombstoned_json(unknown_json(id), id);
            j["cancelled"] = false;
            add_suppression_locked(&j, *attempt, not_after_ms, true);
            return { 404, j };
        }
    }
    Json::Value j;
    j["attempt_id"] = id;
    j["server_instance"] = instance_;
    if (!attempt->cancel()) {
        j["error"] = "attempt '" + id + "' has finished: nothing to cancel";
        j["state"] = "finished";
        j["cancelled"] = false;
        j["attempt"] = attempt->state_json();
        return { 404, j };
    }
    if (wait_ms)
        attempt->wait_finished(wait_ms);
    j["cancelled"] = true;
    j["state"] = attempt->finished() ? "finished" : "stopping";
    j["attempt"] = attempt->state_json();
    return { 200, j };
}

std::pair<int, Json::Value> AttemptRegistry::state(const std::string &id) {
    std::shared_ptr<Attempt> attempt;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        expire_locked();
        auto it = attempts_.find(id);
        if (it == attempts_.end())
            return { 404, unknown_json(id) };
        attempt = it->second;
        if (attempt->tombstone()) {
            // read only: the hold is not extended by looking at it
            Json::Value j = tombstoned_json(unknown_json(id), id);
            add_suppression_locked(&j, *attempt, std::nullopt, true);
            return { 404, j };
        }
    }
    return { 200, attempt->state_json() };
}

void AttemptRegistry::finish(const std::shared_ptr<Attempt> &attempt, const std::string &reason,
                             int status, std::optional<size_t> bytes) {
    attempt->finish(reason, status, bytes);
    if (!attempt->managed())
        return;
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = attempts_.find(attempt->ids().attempt_id);
    if (it != attempts_.end() && it->second == attempt) {
        retained_.push_back(attempt);
        // Its id stays refused until no copy of the request can be admitted any more: until
        // its not_after_ms (and any a refused copy named) + clock_skew_ms, within the cap from
        // now — beyond its retention if need be (expire_locked). Nothing is kept with
        // retention 0 (stated)
        std::optional<uint64_t> not_after = attempt->ids().not_after_ms;
        if (attempt->tomb_not_after_ms_)
            not_after = std::max(not_after.value_or(0), *attempt->tomb_not_after_ms_);
        if (settings_.retention_s && not_after)
            hold_locked(*attempt, not_after, false);
    }
    expire_locked();
}

} // namespace cli
} // namespace mtg
