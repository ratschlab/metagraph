#include "traverse_attempts.hpp"

#include <cmath>
#include <ctime>
#include <random>

#include <fmt/format.h>

#include "traverse.hpp"


namespace mtg {
namespace cli {

using graph::traversal::ExternalStop;
using Clock = Attempt::Clock;

namespace {

Json::Value uint_value(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

// whole milliseconds, rounded up (a duration a ledger leases by is never understated); null
// when there is no bound (an infinite time budget without the HTTP server's cap: the CLI)
Json::Value ms_json(double ms) {
    if (!std::isfinite(ms))
        return Json::Value();
    return uint_value(static_cast<uint64_t>(std::ceil(std::max(0.0, ms))));
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

AttemptIds attempt_ids(const Json::Value &request) {
    AttemptIds ids;
    if (!request.isObject())
        return ids;   // parse_traverse_request names the problem
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
    // an echo without the attempt it belongs to would be usage nobody can reconcile
    if (ids.attempt_id.empty() && !ids.budget_id.empty())
        throw InvalidRequest("request.budget_id: given without request.attempt_id");
    if (ids.attempt_id.empty() && !ids.locus_id.empty())
        throw InvalidRequest("request.locus_id: given without request.attempt_id");
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
    const auto wall = received_wall_ + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                                               t - received_);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(wall.time_since_epoch()).count();
    const std::time_t seconds = static_cast<std::time_t>(ms / 1000);
    std::tm tm {};
    gmtime_r(&seconds, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return fmt::format("{}.{:03d}Z", buf, static_cast<int>(ms % 1000));
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
    walk_until_ms_ = std::max(0.0, bound_ms_ - settings_.allowance_ms / 2);
    bound_set_ = true;
    seeds_.assign(seeds, SeedUsage());
    seeds_requested_ = seeds;
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
            request_stop(ExternalStop::ATTEMPT_DEADLINE);
        } else if (peer_gone_) {
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
        throw AttemptAtBound(fmt::format(
                "the attempt reached the duration bound the server enforces for it ({} ms: the "
                "seeds' time budgets plus the server's allowance, see usage.bound) while its "
                "response was built or written: nothing of it is delivered",
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
    // abandoned, not as finished (review of the stage-4 backend, F5)
    if (!seeds_[index].walked) {
        if (usage.outcome == "abandoned") {
            seeds_abandoned_++;
        } else {
            seeds_walked_++;
        }
    }
    seeds_[index] = usage;
    seeds_[index].started = true;
    seeds_[index].walked = true;
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
    b["seeds"] = uint_value(bound_seeds_);
    b["time_budget_ms"] = bound_set_ ? Json::Value(bound_time_budget_ms_) : Json::Value();
    b["allowance_ms"] = ms_json(settings_.allowance_ms);
    b["walk_until_ms"] = ms_json(walk_until_ms_);
    b["capped_by"] = capped_ ? Json::Value("content_timeout") : Json::Value();
    b["enforced"] = enforced();
    return b;
}

Json::Value Attempt::seeds_json() const {
    Json::Value seeds;
    seeds["requested"] = uint_value(seeds_requested_);
    seeds["started"] = uint_value(seeds_started_);
    // a seed whose walk ended (complete, partial, failed), delivered or not
    seeds["finished"] = uint_value(seeds_walked_);
    // a seed whose walk the client's departure cut: neither finished nor delivered
    seeds["abandoned"] = uint_value(seeds_abandoned_);
    return seeds;
}

Json::Value Attempt::usage_json(const std::string &reason, bool per_seed) const {
    std::lock_guard<std::mutex> lock(mutex_);
    Json::Value u;
    u["attempt_id"] = ids_.attempt_id;
    if (!ids_.budget_id.empty())
        u["budget_id"] = ids_.budget_id;
    if (!ids_.locus_id.empty())
        u["locus_id"] = ids_.locus_id;
    u["server_instance"] = server_instance_;
    u["reason"] = reason;
    u["received_at"] = iso(received_);
    u["stopped_at"] = stopped_at_ ? Json::Value(iso(*stopped_at_)) : Json::Value();
    u["elapsed_ms"] = ms_json(ms_since_received(finished_at_ ? *finished_at_ : now()));
    u["bound_ms"] = ms_json(bound_ms_);
    u["bound"] = bound_json();
    u["seeds"] = seeds_json();
    // The request: work adds up; memory is held at once by the seed being walked and the
    // results of the seeds before it, which stay until the response is written — a walked
    // result its account at the end, a failed or never started one its priced echo
    // (final_bytes; review of the stage-4 backend, F1: they were counted as 0). Under a memory
    // budget a seed's walk holds at most the budget and its soft excess (as observed): what
    // the request held at once is bounded by that, the results before it included
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
        e["index"] = uint_value(i);
        e["outcome"] = s.outcome;
        e["stopped_by"] = s.stopped_by.empty() ? Json::Value() : Json::Value(s.stopped_by);
        e["work_units"] = uint_value(m.work_units);
        e["work_seed"] = uint_value(m.work_seed);
        e["peak_admitted_bytes"] = uint_value(m.memory_peak);
        e["final_bytes"] = uint_value(m.memory_final);
        e["soft_excess_bytes"] = memory_budget_ ? uint_value(m.soft_excess) : Json::Value();
        e["refused_bytes"] = s.refused_bytes ? uint_value(*s.refused_bytes) : Json::Value();
        e["elapsed_ms"] = ms_json(s.elapsed_ms);
        list.append(std::move(e));
    }
    u["work_units"] = uint_value(work);
    Json::Value memory;
    memory["peak_admitted_bytes"] = uint_value(peak);
    // the soft part is observed only under a memory budget, and per seed (each seed's excess
    // over its own budget): the request states the largest
    memory["soft_excess_bytes"] = memory_budget_ ? uint_value(soft) : Json::Value();
    // Without a memory budget the account leaves out the label caches, the lookahead and a
    // level's decoded rows, so nothing bounds what was held (review of the stage-4 backend,
    // F2: 0.79 MB stated for a walk that raised the RSS by 111 MB): null, never a number that
    // reads as a bound
    memory["held_bound_bytes"] = memory_budget_ ? uint_value(bound) : Json::Value();
    u["memory"] = std::move(memory);
    if (per_seed)
        u["per_seed"] = std::move(list);
    return u;
}

Json::Value Attempt::state_json() const {
    std::lock_guard<std::mutex> lock(mutex_);
    Json::Value j;
    j["attempt_id"] = ids_.attempt_id;
    if (!ids_.budget_id.empty())
        j["budget_id"] = ids_.budget_id;
    if (!ids_.locus_id.empty())
        j["locus_id"] = ids_.locus_id;
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
    j["bound_ms"] = ms_json(bound_ms_);
    j["elapsed_ms"] = ms_json(ms_since_received(finished_at_ ? *finished_at_ : now()));
    j["seeds"] = seeds_json();
    Json::Value response;
    response["written"] = finished_at_ ? Json::Value(status_ != 0) : Json::Value();
    response["status"] = status_ ? Json::Value(status_) : Json::Value();
    response["bytes"] = bytes_ ? uint_value(*bytes_) : Json::Value();
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
            usage["elapsed_ms"] = ms_json(ms_since_received(*finished_at_));
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
    return fmt::format("attempts are kept {} s after they finish, at most the last {} (oldest "
                       "dropped first); a cancel of an unknown id tombstones it for {} s, at most "
                       "{} tombstones at once (beyond that the cancel is refused, 429)",
                       settings_.retention_s, settings_.retention_count, settings_.retention_s,
                       settings_.retention_count);
}

void AttemptRegistry::expire_locked() {
    const Clock::time_point t = now();
    const auto keep = std::chrono::seconds(settings_.retention_s);
    auto drop = [&](std::deque<std::shared_ptr<Attempt>> *queue) {
        const std::shared_ptr<Attempt> &old = queue->front();
        auto it = attempts_.find(old->ids().attempt_id);
        if (it != attempts_.end() && it->second == old)
            attempts_.erase(it);
        queue->pop_front();
    };
    while (!retained_.empty()
            && (retained_.size() > settings_.retention_count
                || t - retained_.front()->finished_at() >= keep)) {
        drop(&retained_);
    }
    // by age only: a tombstone promises that its id does not run here for retention_s
    while (!tombstones_.empty() && t - tombstones_.front()->finished_at() >= keep) {
        drop(&tombstones_);
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

std::optional<Json::Value> AttemptRegistry::start(const std::shared_ptr<Attempt> &attempt) {
    std::lock_guard<std::mutex> lock(mutex_);
    expire_locked();
    const std::string &id = attempt->ids().attempt_id;
    auto it = attempts_.find(id);
    if (it != attempts_.end())
        return it->second->state_json();
    attempts_.emplace(id, attempt);
    return std::nullopt;
}

std::pair<int, Json::Value> AttemptRegistry::cancel(const std::string &id, uint64_t wait_ms) {
    std::shared_ptr<Attempt> attempt;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        expire_locked();
        auto it = attempts_.find(id);
        if (it == attempts_.end()) {
            if (tombstones_.size() >= settings_.retention_count) {
                // No tombstone can be kept for the whole retention period without evicting
                // another one early, which would break that one's promise: nothing is promised
                // for this id, said so, and the cancel can be retried
                Json::Value j = unknown_json(id);
                j["error"] = "unknown attempt_id '" + id + "', and it was NOT tombstoned: "
                             + std::to_string(tombstones_.size()) + " cancels of unknown ids are "
                             "held already (" + retention_text() + "), so a request with it "
                             "that arrives later may still run here; retry the cancel later";
                j["cancelled"] = false;
                j["tombstone"] = false;
                return { 429, j };
            }
            // A cancel can overtake its request (still queued, or on the wire): the id is
            // tombstoned for the retention period, so a request that arrives later is refused
            // (409) and never runs here. A 404 from this server_instance therefore means the
            // attempt will not run on it for retention_s, and a ledger can release its
            // reservation
            auto tomb = Attempt::make_tombstone(id, settings_, instance_);
            attempts_.emplace(id, tomb);
            tombstones_.push_back(tomb);
            Json::Value j = unknown_json(id);
            j["cancelled"] = false;
            j["tombstone"] = true;
            return { 404, j };
        }
        attempt = it->second;
    }
    Json::Value j;
    j["attempt_id"] = id;
    j["server_instance"] = instance_;
    if (attempt->tombstone()) {
        j["error"] = "unknown attempt_id '" + id + "': a cancel named it before any request with "
                     "it arrived; a request with it will be refused (409)";
        j["state"] = "unknown";
        j["cancelled"] = false;
        j["tombstone"] = true;
        return { 404, j };
    }
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
    }
    if (attempt->tombstone()) {
        Json::Value j = unknown_json(id);
        j["error"] = "unknown attempt_id '" + id + "': a cancel named it before any request with "
                     "it arrived; a request with it will be refused (409)";
        j["tombstone"] = true;
        return { 404, j };
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
    if (it != attempts_.end() && it->second == attempt)
        retained_.push_back(attempt);
    expire_locked();
}

} // namespace cli
} // namespace mtg
