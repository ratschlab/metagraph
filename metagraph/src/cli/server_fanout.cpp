#include "server_fanout.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include <spdlog/fmt/fmt.h>

#include "cli/json_helpers.hpp"
#include "config/config.hpp"
#include "server_envelope.hpp"
#include "traverse.hpp"
#include "traverse_attempts.hpp"


namespace mtg {
namespace cli {

namespace {

using MilliDouble = std::chrono::duration<double, std::milli>;

SteadyClock::duration steady_ms(double ms) {
    return std::chrono::duration_cast<SteadyClock::duration>(MilliDouble(ms));
}

WallClock::time_point wall_of_ms(uint64_t ms) {
    return WallClock::time_point(std::chrono::milliseconds(ms));
}

uint64_t ms_of_wall(WallClock::time_point t) {
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            t.time_since_epoch()).count();
    return ms > 0 ? static_cast<uint64_t>(ms) : 0;
}

double ms_between(SteadyClock::time_point from, SteadyClock::time_point to) {
    return std::max(0.0, MilliDouble(to - from).count());
}

// days since 1970-01-01 of a civil date (H. Hinnant's days_from_civil): the instants the
// attempt registry states are read back without the platform's timegm
int64_t days_from_civil(int64_t y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int64_t era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + static_cast<int64_t>(doe) - 719468;
}

// an instant iso_utc wrote ("2026-10-31T09:20:58.000Z"), none for anything else
std::optional<WallClock::time_point> parse_iso_utc(const Json::Value &iso) {
    if (!iso.isString())
        return std::nullopt;
    const std::string s = iso.asString();
    if (s.size() != 24 || s[4] != '-' || s[10] != 'T' || s[19] != '.' || s[23] != 'Z')
        return std::nullopt;
    try {
        const int64_t days = days_from_civil(std::stoll(s.substr(0, 4)),
                                             static_cast<unsigned>(std::stoul(s.substr(5, 2))),
                                             static_cast<unsigned>(std::stoul(s.substr(8, 2))));
        const int64_t ms = ((days * 24 + std::stoll(s.substr(11, 2))) * 60
                                + std::stoll(s.substr(14, 2))) * 60'000
                         + std::stoll(s.substr(17, 2)) * 1000 + std::stoll(s.substr(20, 3));
        return WallClock::time_point(std::chrono::milliseconds(ms));
    } catch (const std::exception &) {
        return std::nullopt;
    }
}

Json::Value null_pair_ref() {
    Json::Value ref(Json::objectValue);
    ref["graph"] = Json::Value();
    ref["graph_path"] = Json::Value();
    ref["annotation_path"] = Json::Value();
    return ref;
}

} // namespace

// ---------------------------------------------------------------- the request's handles

RequestHandles read_request_handles(const Json::Value &request) {
    RequestHandles handles;
    if (!request.isObject())
        return handles;   // the route's parser names the problem
    if (request.isMember("request_id")) {
        const Json::Value &v = request["request_id"];
        if (!v.isString() || !valid_attempt_id(v.asString())) {
            throw InvalidRequest("request.request_id: expected a string of 1 to 128 characters "
                                 "from [A-Za-z0-9._:-] (attempt_id's rule)");
        }
        handles.request_id = v.asString();
    }
    handles.not_after_ms = read_not_after_ms(request, "request");
    handles.finish_by_ms = read_finish_by_ms(request, "request");
    return handles;
}

// ---------------------------------------------------------------- the deadlines

SteadyClock::time_point RequestDeadlines::steady_of(WallClock::time_point wall) const {
    return arrival + std::chrono::duration_cast<SteadyClock::duration>(wall - arrival_wall);
}

WallClock::time_point RequestDeadlines::wall_of(SteadyClock::time_point steady) const {
    return arrival_wall + std::chrono::duration_cast<WallClock::duration>(steady - arrival);
}

Json::Value RequestDeadlines::json() const {
    const bool given = not_after_ms || finish_by_ms;
    return deadline_json(not_after_ms, finish_by_ms,
                         given ? std::optional<WallClock::time_point>(effective_wall)
                               : std::nullopt,
                         given ? std::optional<WallClock::time_point>(cutoff_wall)
                               : std::nullopt);
}

RequestDeadlines request_deadlines(WallClock::time_point arrival_wall,
                                   WallClock::time_point now_wall,
                                   SteadyClock::time_point now,
                                   double budget_ms,
                                   const RequestHandles &handles,
                                   uint64_t skew_ms,
                                   double reserve_ms,
                                   double transport_ms) {
    RequestDeadlines d;
    d.not_after_ms = handles.not_after_ms;
    d.finish_by_ms = handles.finish_by_ms;
    // the request arrived when its header was read (the transport's timeout runs from then),
    // which its handler learns later; never in the future
    if (arrival_wall.time_since_epoch().count() == 0 || arrival_wall > now_wall)
        arrival_wall = now_wall;
    d.arrival_wall = arrival_wall;
    d.arrival = now - std::chrono::duration_cast<SteadyClock::duration>(now_wall - arrival_wall);
    // a budget past the transport's end binds nothing either (the cap of /pattern's budget is
    // the operator's, up to 2^53 - 1 ms): it is saturated just past it before it meets the
    // steady clock's nanoseconds, which it would overflow
    const double transport_end_ms = std::max(0.0, transport_ms);
    const SteadyClock::time_point budget_end
            = d.arrival + steady_ms(std::min(std::max(0.0, budget_ms), transport_end_ms + 1));
    const SteadyClock::time_point transport_end = d.arrival + steady_ms(transport_end_ms);
    d.transport_bound = !(budget_end < transport_end);
    d.effective = std::min(budget_end, transport_end);
    if (handles.finish_by_ms) {
        uint64_t finish_ms = *handles.finish_by_ms > skew_ms
                ? *handles.finish_by_ms - skew_ms : 0;
        // A finish after the transport's end cannot bind (the effective deadline is at most
        // that end): saturated just past it, so that a far-future finish_by_ms (read up to
        // 2^53 - 1) is no instant the clocks cannot hold (on Linux the system clock counts
        // nanoseconds, which overflow past the year 2262, and the steady difference past 292
        // years)
        const uint64_t past_transport = ms_of_wall(arrival_wall)
                + static_cast<uint64_t>(std::ceil(std::max(0.0, transport_ms))) + 1;
        finish_ms = std::min(finish_ms, past_transport);
        d.finish_wall = wall_of_ms(finish_ms);
        d.finish = d.steady_of(*d.finish_wall);
        if (*d.finish < d.effective) {
            d.effective = *d.finish;
            d.transport_bound = false;
        }
    }
    d.cutoff = d.effective - steady_ms(reserve_ms);
    d.effective_wall = d.wall_of(d.effective);
    d.cutoff_wall = d.wall_of(d.cutoff);
    return d;
}

// ---------------------------------------------------------------- the request registry

const char* to_string(PairState state) {
    switch (state) {
        case PairState::QUEUED: return "queued";
        case PairState::WAITING_MEMORY: return "waiting_memory";
        case PairState::LOADING: return "loading";
        case PairState::RUNNING: return "running";
        case PairState::ENDED: return "ended";
        case PairState::NEVER_STARTED: return "never_started";
    }
    return "unknown";
}

RequestEntry::RequestEntry(std::optional<std::string> request_id, std::string route,
                           std::string server_instance, WallClock::time_point received_wall,
                           std::vector<Json::Value> pairs, Clocks clocks)
      : request_id_(std::move(request_id)), route_(std::move(route)),
        server_instance_(std::move(server_instance)), received_wall_(received_wall),
        clocks_(std::move(clocks)) {
    if (!clocks_.steady)
        clocks_.steady = []() { return SteadyClock::now(); };
    if (!clocks_.wall)
        clocks_.wall = []() { return WallClock::now(); };
    // the steady instant of the arrival, through one reading of both clocks
    const WallClock::time_point wall = clocks_.wall();
    received_ = clocks_.steady()
              - std::chrono::duration_cast<SteadyClock::duration>(
                        wall > received_wall_ ? wall - received_wall_ : WallClock::duration(0));
    pairs_.resize(pairs.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        pairs_[i].ref = std::move(pairs[i]);
    }
}

WallClock::time_point RequestEntry::wall_of_locked(SteadyClock::time_point t) const {
    return received_wall_ + std::chrono::duration_cast<WallClock::duration>(t - received_);
}

void RequestEntry::set_on_ended(std::function<void()> on_ended) {
    std::lock_guard<std::mutex> lock(mutex_);
    on_ended_ = std::move(on_ended);
}

void RequestEntry::pair_waiting_memory(size_t i) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (pairs_[i].state == PairState::QUEUED)
        pairs_[i].state = PairState::WAITING_MEMORY;
}

void RequestEntry::pair_loading(size_t i, uint64_t bytes) {
    std::lock_guard<std::mutex> lock(mutex_);
    Pair &p = pairs_[i];
    p.state = PairState::LOADING;
    p.memory_bytes = bytes;
    if (!p.started) {
        p.started = clocks_.steady();
        p.started_wall = wall_of_locked(*p.started);
    }
}

void RequestEntry::pair_running(size_t i) {
    std::lock_guard<std::mutex> lock(mutex_);
    Pair &p = pairs_[i];
    p.state = PairState::RUNNING;
    if (!p.started) {
        p.started = clocks_.steady();
        p.started_wall = wall_of_locked(*p.started);
    }
}

bool RequestEntry::pair_started(size_t i) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pairs_[i].started.has_value();
}

PairState RequestEntry::pair_state(size_t i) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pairs_[i].state;
}

void RequestEntry::pair_cpu(size_t i, double cpu_ms) {
    std::lock_guard<std::mutex> lock(mutex_);
    pairs_[i].cpu_ms = cpu_ms;
    pairs_[i].has_cpu = true;
}

void RequestEntry::pair_outcome(size_t i, const std::string &outcome) {
    std::lock_guard<std::mutex> lock(mutex_);
    pairs_[i].outcome = outcome;
}

void RequestEntry::pair_ended(size_t i, std::optional<uint64_t> memory_bytes,
                              std::optional<double> cpu_ms) {
    std::function<void()> ended;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        Pair &p = pairs_[i];
        if (p.state == PairState::ENDED)
            return;
        if (!p.started) {
            p.state = PairState::NEVER_STARTED;
        } else {
            p.state = PairState::ENDED;
            p.finished = clocks_.steady();
            if (memory_bytes)
                p.memory_bytes = *memory_bytes;
            if (cpu_ms) {
                p.cpu_ms = *cpu_ms;
                p.has_cpu = true;
            }
        }
        ended = end_locked();
    }
    if (ended)
        ended();
}

void RequestEntry::pair_never_started(size_t i) {
    std::function<void()> ended;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        Pair &p = pairs_[i];
        if (p.started || p.state == PairState::NEVER_STARTED)
            return;
        p.state = PairState::NEVER_STARTED;
        ended = end_locked();
    }
    if (ended)
        ended();
}

Json::Value RequestEntry::pair_usage_now(size_t i, const std::string &state) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const Pair &p = pairs_[i];
    Json::Value u;
    u["state"] = state;
    u["started_at"] = p.started_wall ? Json::Value(iso_utc(*p.started_wall)) : Json::Value();
    u["elapsed_ms"] = p.started ? number_json(std::round(ms_between(*p.started, clocks_.steady())))
                                : Json::Value();
    u["cpu_ms"] = p.has_cpu ? Json::Value(p.cpu_ms) : Json::Value();
    u["memory_bytes"] = uint_json(p.memory_bytes);
    return u;
}

void RequestEntry::response(int status, std::optional<size_t> bytes) {
    std::function<void()> ended;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (responded_)
            return;
        responded_ = true;
        status_ = status;
        bytes_ = bytes;
        ended = end_locked();
    }
    if (ended)
        ended();
}

std::function<void()> RequestEntry::end_locked() {
    if (!responded_ || finished_wall_)
        return nullptr;
    for (const Pair &p : pairs_) {
        if (p.state != PairState::ENDED && p.state != PairState::NEVER_STARTED)
            return nullptr;
    }
    finished_wall_ = wall_of_locked(clocks_.steady());
    // what the request ended with: its client gone (nothing written, whatever its workers
    // did), its answer, the whole request's deadline, or another failure
    reason_ = !status_ ? "client_gone"
            : status_ == 200 ? "answered"
            : status_ == 503 ? "deadline"
            : "error";
    if (ended_called_ || !on_ended_)
        return nullptr;
    ended_called_ = true;
    return on_ended_;
}

bool RequestEntry::ended() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return finished_wall_.has_value();
}

std::optional<WallClock::time_point> RequestEntry::finished_at() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return finished_wall_;
}

Json::Value RequestEntry::usage_locked(const Pair &p) const {
    if (!p.started)
        return Json::Value();
    Json::Value u;
    u["started_at"] = iso_utc(*p.started_wall);
    u["finished_at"] = p.finished ? Json::Value(iso_utc(wall_of_locked(*p.finished)))
                                  : Json::Value();
    u["elapsed_ms"] = number_json(std::round(
            ms_between(*p.started, p.finished ? *p.finished : clocks_.steady())));
    u["cpu_ms"] = p.has_cpu ? Json::Value(p.cpu_ms) : Json::Value();
    u["memory_bytes"] = uint_json(p.memory_bytes);
    return u;
}

Json::Value RequestEntry::json(uint64_t retention_s) const {
    std::lock_guard<std::mutex> lock(mutex_);
    Json::Value j;
    j["request_id"] = request_id_ ? Json::Value(*request_id_) : Json::Value();
    j["server_instance"] = server_instance_;
    j["route"] = route_;
    j["state"] = finished_wall_ ? "ended" : "running";
    j["received_at"] = iso_utc(received_wall_);
    j["finished_at"] = finished_wall_ ? Json::Value(iso_utc(*finished_wall_)) : Json::Value();
    j["reason"] = finished_wall_ ? Json::Value(reason_) : Json::Value();
    if (responded_) {
        Json::Value r;
        r["status"] = status_ ? Json::Value(status_) : Json::Value();
        r["bytes"] = bytes_ ? uint_json(*bytes_) : Json::Value();
        r["written"] = status_ != 0;
        j["response"] = std::move(r);
    } else {
        j["response"] = Json::Value();
    }
    j["retained_until"] = finished_wall_
            ? Json::Value(iso_utc(*finished_wall_ + std::chrono::seconds(retention_s)))
            : Json::Value();
    Json::Value pairs(Json::arrayValue);
    for (const Pair &p : pairs_) {
        Json::Value pj = p.ref;
        pj["state"] = to_string(p.state);
        pj["outcome"] = p.outcome.empty() ? Json::Value() : Json::Value(p.outcome);
        pj["usage"] = usage_locked(p);
        pairs.append(std::move(pj));
    }
    j["pairs"] = std::move(pairs);
    return j;
}

SteadyClock::time_point RequestRegistry::Shared::now() const {
    return settings.clocks.steady ? settings.clocks.steady() : SteadyClock::now();
}

void RequestRegistry::Shared::expire_locked() {
    const SteadyClock::time_point t = now();
    const auto keep = std::chrono::seconds(settings.retention_s);
    auto expired = [&](const auto &queue) {
        return !queue.empty()
            && (queue.size() > settings.retention_count || t - queue.front().first >= keep);
    };
    while (expired(ended)) {
        const std::string id = ended.front().second;
        ended.pop_front();
        auto it = entries.find(id);
        // (an id that ran again after it was dropped is a newer entry, still held)
        if (it != entries.end() && it->second->ended()) {
            bool newer = false;
            for (const auto &[_, later] : ended) {
                newer |= later == id;
            }
            if (!newer)
                entries.erase(it);
        }
    }
}

RequestRegistry::RequestRegistry(RequestSettings settings)
      : shared_(std::make_shared<Shared>()) {
    shared_->settings = std::move(settings);
}

std::string RequestRegistry::retention_text() const {
    const RequestSettings &s = shared_->settings;
    if (!s.retention_s || !s.retention_count)
        return "requests are not kept after they finish (--request-retention-s 0)";
    return fmt::format("requests are kept {} s after they finish, at most {}", s.retention_s,
                       s.retention_count);
}

Json::Value RequestRegistry::unknown_json(const std::string &request_id) const {
    Json::Value j;
    j["error"] = "request_id '" + request_id + "' is not known on this server: never received, "
                 "or no longer retained (" + retention_text() + ")";
    j["request_id"] = request_id;
    j["server_instance"] = shared_->settings.server_instance;
    return j;
}

std::string RequestRegistry::repeated_text(const std::string &request_id) const {
    return "request_id '" + request_id + "' is running, or retained after it finished, on this "
           "server (" + retention_text() + "): it is not run again; its state is this body";
}

std::optional<Json::Value> RequestRegistry::start(
        const std::shared_ptr<RequestEntry> &entry,
        const std::function<std::optional<Json::Value>()> &admission) {
    assert(entry->request_id());
    const std::string &id = *entry->request_id();
    std::lock_guard<std::mutex> lock(shared_->mutex);
    shared_->expire_locked();
    auto it = shared_->entries.find(id);
    if (it != shared_->entries.end()) {
        Json::Value body = it->second->json(shared_->settings.retention_s);
        body["error"] = repeated_text(id);
        return body;
    }
    if (admission) {
        if (auto refused = admission())
            return refused;
    }
    shared_->entries.emplace(id, entry);
    // the retention starts when the request ended (its response recorded, its workers ended);
    // the registry is held weakly, so that an entry kept by a late worker does not keep it
    std::weak_ptr<Shared> weak = shared_;
    entry->set_on_ended([weak, id]() {
        auto shared = weak.lock();
        if (!shared)
            return;
        std::lock_guard<std::mutex> lock(shared->mutex);
        shared->ended.emplace_back(shared->now(), id);
        shared->expire_locked();
    });
    return std::nullopt;
}

void RequestRegistry::forget(const std::shared_ptr<RequestEntry> &entry) {
    if (!entry->request_id())
        return;
    // no retention, whenever its workers end (it holds no record in |ended| yet: it has not
    // recorded its response)
    entry->set_on_ended(nullptr);
    std::lock_guard<std::mutex> lock(shared_->mutex);
    auto it = shared_->entries.find(*entry->request_id());
    if (it != shared_->entries.end() && it->second == entry)
        shared_->entries.erase(it);
}

std::pair<int, Json::Value> RequestRegistry::state(const std::string &request_id) const {
    std::lock_guard<std::mutex> lock(shared_->mutex);
    shared_->expire_locked();
    auto it = shared_->entries.find(request_id);
    if (it == shared_->entries.end())
        return { 404, unknown_json(request_id) };
    return { 200, it->second->json(shared_->settings.retention_s) };
}

size_t RequestRegistry::size() const {
    std::lock_guard<std::mutex> lock(shared_->mutex);
    shared_->expire_locked();
    return shared_->entries.size();
}

Json::Value traverse_request_json(const Json::Value &attempt, const Json::Value &pair,
                                  uint64_t retention_s) {
    const bool ended = attempt["state"].asString() == "finished";
    Json::Value j;
    j["request_id"] = attempt["attempt_id"];
    j["server_instance"] = attempt["server_instance"];
    j["route"] = "/traverse";
    j["state"] = ended ? "ended" : "running";
    j["received_at"] = attempt["received_at"];
    j["finished_at"] = attempt["finished_at"];
    j["reason"] = attempt["reason"];
    j["response"] = ended ? attempt["response"] : Json::Value();
    const auto finished = parse_iso_utc(attempt["finished_at"]);
    j["retained_until"] = ended && finished
            ? Json::Value(iso_utc(*finished + std::chrono::seconds(retention_s)))
            : Json::Value();
    Json::Value p = pair.isObject() ? pair : null_pair_ref();
    p["state"] = ended ? "ended" : "running";
    const Json::Value &status = attempt["response"]["status"];
    p["outcome"] = !ended || !status.isIntegral() ? Json::Value()
                 : status.asInt() == 200 ? Json::Value("answered") : Json::Value("refused");
    // the attempt's usage is the request's (attempt.usage)
    p["usage"] = Json::Value();
    j["pairs"] = Json::Value(Json::arrayValue);
    j["pairs"].append(std::move(p));
    j["attempt"] = attempt;
    return j;
}

Json::Value request_as_attempt_json(const Json::Value &entry) {
    Json::Value attempt = entry;
    attempt["attempt_id"] = entry["request_id"];
    attempt["state"] = entry["state"].asString() == "ended" ? "finished" : "running";
    return attempt;
}

std::optional<Json::Value> IdAdmission::start_request(
        const std::shared_ptr<RequestEntry> &entry,
        const std::function<std::optional<Json::Value>()> &admission) {
    assert(entry->request_id());
    std::lock_guard<std::mutex> lock(mutex_);
    return requests_.start(entry, [&]() -> std::optional<Json::Value> {
        const std::string &id = *entry->request_id();
        Json::Value pair;
        auto [status, attempt] = attempts_.state(id, &pair);
        if (status == 200 && attempt.get("state", "").asString() != "unknown") {
            Json::Value body = traverse_request_json(attempt, pair,
                                                     attempts_.settings().retention_s);
            body["error"] = requests_.repeated_text(id);
            return body;
        }
        return admission ? admission() : std::nullopt;
    });
}

std::optional<AttemptRegistry::StartRefusal> IdAdmission::start_attempt(
        const std::shared_ptr<Attempt> &attempt) {
    std::lock_guard<std::mutex> lock(mutex_);
    return attempts_.start(attempt, [&]() -> std::optional<Json::Value> {
        auto [status, entry] = requests_.state(attempt->ids().attempt_id);
        if (status != 200)
            return std::nullopt;
        return request_as_attempt_json(entry);
    });
}

// ---------------------------------------------------------------- the fan-out

namespace {

std::atomic<size_t> g_fanouts_alive { 0 };

} // namespace

/**
 * The COMPLETE state of one fan-out request, shared by its handler and every worker (each holds
 * a shared_ptr): the settings and the job (which holds the parsed request), the slots and the
 * stop flag. A worker that outlives the HTTP answer reads nothing else. |mutex| guards the
 * slots, |closed| and the transitions of the registry entry that decide an outcome (a pair
 * starting, the handler closing), so that a pair is either started or not_reached, never both.
 */
class FanoutState {
  public:
    struct Slot {
        enum Phase { QUEUED, WAITING, STARTED, DONE };
        Phase phase = QUEUED;
        bool loading = false;
        // the answer's entry is decided (by the worker, or by the handler at the cutoff)
        bool resolved = false;
        Json::Value entry;
        bool aborted = false;
        common::ThreadUsage usage = common::ThreadUsage::zero();
        bool usage_final = false;
        // the worker's thread: its CPU time at the job's start, and the next sample's time
        common::ThreadUsage start_usage = common::ThreadUsage::zero();
        SteadyClock::time_point next_sample;
    };

    FanoutState(FanoutSettings s, std::shared_ptr<RequestEntry> e, PairJob j)
          : settings(std::move(s)), entry(std::move(e)), job(std::move(j)),
            inline_run(settings.pairs.size() == 1), slots(settings.pairs.size()) {
        if (!settings.clock)
            settings.clock = []() { return SteadyClock::now(); };
        if (!settings.wall_clock)
            settings.wall_clock = []() { return WallClock::now(); };
        // the handler's check is the handler's alone (it reads the HTTP request, which a
        // worker may outlive): kept apart, given to the inline pair only while it runs
        inline_gone = std::move(settings.gone);
        settings.gone = nullptr;
        g_fanouts_alive++;
    }
    ~FanoutState() { g_fanouts_alive--; }

    SteadyClock::time_point now() const { return settings.clock(); }
    bool not_after_passed() const {
        return settings.deadlines.not_after_ms
            && ms_of_wall(settings.wall_clock()) > *settings.deadlines.not_after_ms;
    }
    std::string route_name() const {
        return settings.route.size() && settings.route[0] == '/' ? settings.route.substr(1)
                                                                 : settings.route;
    }
    std::string label(size_t i) const {
        const SelectedPair &p = settings.pairs[i];
        return "(" + p.name + ", " + p.pair.first + ", " + p.pair.second + ")";
    }
    std::string fp(size_t i) const {
        return i < settings.index_fps.size() ? settings.index_fps[i] : std::string();
    }
    std::string usage_route() const {
        return entry->request_id() ? "its final usage is in GET /requests/" + *entry->request_id()
                                   : "the request has no request_id: its final usage is not "
                                     "retained";
    }

    // why a pair not started by now is not reached (section 4): its not_after_ms passed
    // (|expired|), else the request's work cutoff, the transport's when the transport's cap
    // is its effective deadline — one rule, whichever thread decides the pair
    std::string not_reached_reason(bool expired) const {
        return expired ? "not_after_ms"
             : settings.deadlines.transport_bound ? "transport" : "cutoff";
    }

    // the not_reached entry of pair |i| (under |mutex|): |reason| cutoff | not_after_ms |
    // transport, |why| what the pair was doing
    void not_reached_locked(size_t i, const std::string &reason, const std::string &why) {
        Slot &slot = slots[i];
        Json::Value body;
        body["code"] = "not_reached";
        body["reason"] = reason;
        if (reason == "not_after_ms") {
            body["error"] = route_name() + ": pair " + label(i) + " was not started: not_after_ms "
                          + std::to_string(*settings.deadlines.not_after_ms) + " had passed "
                          + why;
        } else if (reason == "transport") {
            body["error"] = route_name() + ": pair " + label(i) + " was not started before the "
                            "request's work cutoff (" + iso_utc(settings.deadlines.cutoff_wall)
                          + ", the transport's cap of " + std::to_string(kServerMaxDeadlineMs)
                          + " ms after its arrival less the finalisation reserve): " + why;
        } else {
            body["error"] = route_name() + ": pair " + label(i) + " was not started before the "
                            "request's work cutoff (" + iso_utc(settings.deadlines.cutoff_wall)
                          + "): " + why;
        }
        slot.entry = envelope_entry_refused(settings.pairs[i], fp(i), 503, body);
        slot.resolved = true;
        entry->pair_outcome(i, "not_reached");
        entry->pair_never_started(i);
    }

    // the unfinished entry of started pair |i| at the cutoff (under |mutex|)
    void unfinished_locked(size_t i) {
        Slot &slot = slots[i];
        const bool loading = slot.loading;
        Json::Value body;
        body["code"] = "unfinished";
        body["reason"] = loading ? "loading" : "cutoff";
        body["usage"] = entry->pair_usage_now(i, loading ? "loading" : "stopping");
        const std::string cutoff = iso_utc(settings.deadlines.cutoff_wall);
        const std::string started = body["usage"]["started_at"].asString();
        body["error"] = loading
            ? route_name() + ": pair " + label(i) + " was loading its index into RAM for the "
              "request (since " + started + ") at the request's work cutoff (" + cutoff + "): a "
              "load is not interrupted, its work is told to stop once it ends; " + usage_route()
            : route_name() + ": pair " + label(i) + " started at " + started + " and had not "
              "finished at the request's work cutoff (" + cutoff + "): told to stop; "
              + usage_route();
        slot.entry = envelope_entry_refused(settings.pairs[i], fp(i), 503, body);
        slot.resolved = true;
        entry->pair_outcome(i, "unfinished");
    }

    // The fan-out failed while it handed its pairs over (under |mutex|; the exception leaves
    // run_fanout, so nothing of the entries is written, and nothing here allocates): the
    // request closes, the workers handed over are told to stop, and every pair not started
    // ends now as not_reached, so that the registry entry ends with the response and its
    // workers
    void abandon_locked() {
        closed = true;
        stop.store(true, std::memory_order_release);
        for (size_t i = 0; i < slots.size(); ++i) {
            Slot &slot = slots[i];
            if (slot.resolved)
                continue;
            slot.resolved = true;
            if (slot.phase == Slot::STARTED) {
                entry->pair_outcome(i, "unfinished");
            } else {
                entry->pair_outcome(i, "not_reached");
                entry->pair_never_started(i);
            }
        }
        changed.notify_all();
    }

    // The handler stops waiting (under |mutex|): the cutoff passed (|gone| false) or the
    // client left. Every pair not done is decided now, and the work told to stop
    void close_locked(bool gone) {
        closed = true;
        stop.store(true, std::memory_order_release);
        const bool expired = not_after_passed();
        for (size_t i = 0; i < slots.size(); ++i) {
            Slot &slot = slots[i];
            if (slot.resolved)
                continue;
            if (slot.phase == Slot::STARTED) {
                unfinished_locked(i);
                continue;
            }
            const std::string what = slot.phase == Slot::WAITING
                    ? "waited for the memory to load its index"
                    : "waited for a pool thread";
            if (gone) {
                // nothing is answered; the registry states what became of the pair
                slot.resolved = true;
                entry->pair_outcome(i, "not_reached");
                entry->pair_never_started(i);
            } else {
                not_reached_locked(i, not_reached_reason(expired),
                                   (expired ? "while it " : "it ") + what);
            }
        }
        changed.notify_all();
    }

    // samples the worker's CPU time for the registry (the worker's thread, at most every
    // 100 ms): an unfinished pair states its CPU time so far
    void sample(size_t i) {
        Slot &slot = slots[i];
        const SteadyClock::time_point t = now();
        if (t < slot.next_sample)
            return;
        slot.next_sample = t + std::chrono::milliseconds(100);
        const common::ThreadUsage u = common::ThreadUsage::now().since(slot.start_usage);
        if (u.has_cpu)
            entry->pair_cpu(i, static_cast<double>(u.cpu_ns) / 1e6);
    }

    bool stopped(size_t i) {
        // (the slot's sampling fields are the worker's own: no lock)
        sample(i);
        if (stop.load(std::memory_order_acquire))
            return true;
        if (inline_run && inline_gone && inline_gone()) {
            stop.store(true, std::memory_order_release);
            return true;
        }
        return false;
    }

    // the limit of an in_ram memory wait: the cutoff in a pooled fan-out (the handler stops
    // waiting there), finish_by_ms - the skew allowance or the transport's cap inline (a
    // one-pair request's budgets start after its load, so its own budget does not bound the
    // wait)
    SteadyClock::time_point wait_limit() const {
        const RequestDeadlines &d = settings.deadlines;
        if (!inline_run)
            return d.cutoff;
        SteadyClock::time_point limit = d.arrival + steady_ms(kServerMaxDeadlineMs);
        if (d.finish)
            limit = std::min(limit, *d.finish);
        return limit;
    }

    // a pair starts (the worker's thread): false when it may not. Every limit is read under
    // |mutex| with the transition, whichever thread reaches the cutoff first: a pooled pair
    // starts nothing past the cutoff, also when the handler has not closed the request yet (a
    // load cannot be interrupted once started). An inline pair is its request (section 4): it
    // is not cut, and only a load is not started once its client left or its wait limit passed
    // (the in_ram memory wait is not entered once finish_by_ms has passed, and a load would
    // only hold memory for an answer already due): the pair is then not reached
    bool start(size_t i, bool loading, uint64_t bytes) {
        // (its client's check reads the socket, outside the lock: no other thread decides an
        // inline pair)
        const bool inline_stop = inline_run && loading && stopped(i);
        std::lock_guard<std::mutex> lock(mutex);
        Slot &slot = slots[i];
        if (slot.phase == Slot::STARTED) {
            // the load ended: the work begins
            slot.loading = false;
            entry->pair_running(i);
            return true;
        }
        if (slot.resolved || closed || inline_stop)
            return false;
        if (not_after_passed()) {
            not_reached_locked(i, "not_after_ms", loading ? "when the memory to load its index "
                                                            "was granted"
                                                          : "when it was to start");
            changed.notify_all();
            return false;
        }
        if (inline_run) {
            if (loading && now() >= wait_limit())
                return false;
        } else if (now() >= wait_limit()) {
            not_reached_locked(i, not_reached_reason(false),
                               slot.phase == Slot::WAITING
                                   ? "it waited for the memory to load its index"
                                   : "it waited for a pool thread");
            changed.notify_all();
            return false;
        }
        slot.phase = Slot::STARTED;
        slot.loading = loading;
        if (loading) {
            entry->pair_loading(i, bytes);
        } else {
            entry->pair_running(i);
        }
        return true;
    }

    // the worker's result (the worker's thread): published once, unless the handler decided
    // the pair already (it is then only recorded in the registry)
    void publish(size_t i, PairResult result) {
        std::lock_guard<std::mutex> lock(mutex);
        Slot &slot = slots[i];
        const bool started = slot.phase == Slot::STARTED;
        if (!slot.resolved) {
            if (result.not_started || (!started && result.aborted)) {
                if (inline_run && stop.load(std::memory_order_acquire)) {
                    // an inline pair whose client left while it waited for memory
                    slot.aborted = true;
                    slot.resolved = true;
                    entry->pair_outcome(i, "not_reached");
                } else {
                    // (a worker whose wait gave up at the cutoff may decide the pair before
                    // the handler does: the same reason either way)
                    const bool expired = not_after_passed();
                    not_reached_locked(i, not_reached_reason(expired),
                                       expired ? "while it waited for the memory to load its "
                                                 "index"
                                               : "it waited for the memory to load its index");
                }
            } else if (result.aborted) {
                slot.aborted = true;
                slot.resolved = true;
                entry->pair_outcome(i, "unfinished");
            } else {
                slot.entry = std::move(result.entry);
                slot.resolved = true;
                slot.usage = result.usage;
                slot.usage_final = true;
                entry->pair_outcome(i, result.answered ? "answered" : "refused");
            }
        }
        slot.phase = Slot::DONE;
        if (started) {
            entry->pair_ended(i, result.memory_bytes,
                              result.usage.has_cpu
                                  ? std::optional<double>(static_cast<double>(result.usage.cpu_ns)
                                                          / 1e6)
                                  : std::nullopt);
        } else {
            entry->pair_never_started(i);
        }
        changed.notify_all();
    }

    // one pair, on the thread that took it
    static void run_pair(const std::shared_ptr<FanoutState> &self, size_t i) {
        {
            std::lock_guard<std::mutex> lock(self->mutex);
            Slot &slot = self->slots[i];
            if (slot.resolved) {
                // the handler decided it while it waited for a thread: nothing to do
                slot.phase = Slot::DONE;
                self->entry->pair_never_started(i);
                self->changed.notify_all();
                return;
            }
            if (self->not_after_passed()) {
                self->not_reached_locked(i, "not_after_ms",
                                         self->inline_run ? "when its handler was to start it"
                                                          : "when a pool thread took it");
                slot.phase = Slot::DONE;
                self->changed.notify_all();
                return;
            }
            slot.start_usage = common::ThreadUsage::now();
            slot.next_sample = self->now();
        }
        const common::ThreadUsage start = common::ThreadUsage::now();
        PairResult result;
        try {
            PairControl control(self, i);
            result = self->job(control);
        } catch (...) {
            // a job answers its failures in its entry; anything that escapes it is the pair's
            // internal error, never the process's (a pool thread must not throw)
            Json::Value body;
            body["error"] = "Internal server error";
            result = PairResult();
            result.entry = envelope_entry_refused(self->settings.pairs[i], self->fp(i), 500,
                                                  body);
        }
        result.usage = common::ThreadUsage::now().since(start);
        self->publish(i, std::move(result));
    }

    FanoutSettings settings;
    std::shared_ptr<RequestEntry> entry;
    PairJob job;
    const bool inline_run;
    std::function<bool()> inline_gone;
    std::atomic<bool> stop { false };
    std::mutex mutex;
    std::condition_variable changed;
    std::vector<Slot> slots;
    bool closed = false;
};

PairControl::PairControl(std::shared_ptr<FanoutState> state, size_t index)
      : state_(std::move(state)), index_(index) {}

const SelectedPair& PairControl::pair() const { return state_->settings.pairs[index_]; }

const RequestDeadlines& PairControl::deadlines() const { return state_->settings.deadlines; }

bool PairControl::inline_run() const { return state_->inline_run; }

bool PairControl::stopped() { return state_->stopped(index_); }

std::function<bool()> PairControl::stop_predicate() {
    return [state = state_, i = index_]() { return state->stopped(i); };
}

std::function<bool()> PairControl::wait_predicate() {
    return [state = state_, i = index_]() {
        return state->stopped(i) || state->not_after_passed()
            || state->now() >= state->wait_limit();
    };
}

void PairControl::waiting_memory() {
    std::lock_guard<std::mutex> lock(state_->mutex);
    FanoutState::Slot &slot = state_->slots[index_];
    if (slot.phase == FanoutState::Slot::QUEUED)
        slot.phase = FanoutState::Slot::WAITING;
    state_->entry->pair_waiting_memory(index_);
}

bool PairControl::start_loading(uint64_t bytes) { return state_->start(index_, true, bytes); }

bool PairControl::start_running() { return state_->start(index_, false, 0); }

std::optional<SteadyClock::time_point> PairControl::deadline_cap() const {
    // finish_by_ms (less the skew allowance): each pair keeps its own budget, which starts when
    // it starts (section 4), capped by it as a single-graph request is
    return state_->settings.deadlines.finish;
}

std::optional<SteadyClock::time_point> PairControl::cutoff_cap() const {
    // a pooled pair's budget ends by the cutoff, its time stop answered rather than cut
    // unfinished (the owner's decision of 2026-11-01); an inline pair is the request
    if (state_->inline_run)
        return std::nullopt;
    return state_->settings.deadlines.cutoff;
}

FanoutOutcome run_fanout(const FanoutSettings &settings,
                         const std::shared_ptr<RequestEntry> &entry,
                         PairJob job,
                         const Launch &launch) {
    auto state = std::make_shared<FanoutState>(settings, entry, std::move(job));
    FanoutOutcome outcome;
    const size_t n = state->slots.size();
    if (state->inline_run) {
        // one pair: on the handler's thread, as a single-graph request (no pool queue)
        FanoutState::run_pair(state, 0);
        state->inline_gone = nullptr;
        std::lock_guard<std::mutex> lock(state->mutex);
        FanoutState::Slot &slot = state->slots[0];
        outcome.client_gone = slot.aborted;
        outcome.entries.push_back(std::move(slot.entry));
        return outcome;
    }
    const RequestDeadlines &d = state->settings.deadlines;
    const std::function<bool()> &gone = state->inline_gone;
    std::unique_lock<std::mutex> lock(state->mutex);
    auto all_resolved = [&]() {
        return std::all_of(state->slots.begin(), state->slots.end(),
                           [](const FanoutState::Slot &s) { return s.resolved; });
    };
    // the pairs handed to the pool so far, in their order
    size_t next = 0;
    // hands over the next pairs (without the lock: |launch| may wait for room until |until|);
    // each task holds the state: it may run after the handler answered. A pair decided
    // meanwhile is not handed over
    auto hand_over = [&](SteadyClock::time_point until) {
        for ( ; next < n; ++next) {
            if (state->slots[next].resolved)
                continue;
            lock.unlock();
            const bool taken = launch(
                    [state, i = next]() { FanoutState::run_pair(state, i); }, until);
            lock.lock();
            if (!taken)
                return;
        }
    };
    while (!all_resolved()) {
        const SteadyClock::time_point t = state->now();
        if (t >= d.cutoff) {
            state->close_locked(false);
            break;
        }
        if (state->not_after_passed()) {
            // no pair still queued (or waiting for memory) can start any more: decided now
            // rather than waited for behind the pool's queue
            for (size_t i = 0; i < n; ++i) {
                FanoutState::Slot &slot = state->slots[i];
                if (!slot.resolved && slot.phase != FanoutState::Slot::STARTED) {
                    state->not_reached_locked(
                            i, "not_after_ms",
                            slot.phase == FanoutState::Slot::WAITING
                                ? "while it waited for the memory to load its index"
                                : "while it waited for a pool thread");
                }
            }
            if (all_resolved())
                break;
        }
        lock.unlock();
        const bool left = gone && gone();
        lock.lock();
        if (left) {
            state->close_locked(true);
            outcome.client_gone = true;
            break;
        }
        SteadyClock::time_point until = std::min(d.cutoff,
                                                 t + std::chrono::milliseconds(settings.poll_ms));
        try {
            hand_over(until);
        } catch (...) {
            // a pair handed over runs on; the pairs after it never will: without this, they
            // would hold the registry entry running for good
            if (!lock.owns_lock())
                lock.lock();
            state->abandon_locked();
            state->inline_gone = nullptr;
            throw;
        }
        // (a full queue was the wait: the handler's checks come first). A pair may have ended
        // while the lock was released to hand it over, its notice sent to nobody: the wait
        // reads the slots first
        if (next == n)
            state->changed.wait_until(lock, until, all_resolved);
    }
    // a pair decided by the handler is told to stop (its worker ends at its next check)
    if (!state->closed && std::any_of(state->slots.begin(), state->slots.end(),
                                      [](const FanoutState::Slot &s) {
                                          return s.phase != FanoutState::Slot::DONE;
                                      })) {
        state->closed = true;
        state->stop.store(true, std::memory_order_release);
    }
    for (FanoutState::Slot &slot : state->slots) {
        if (slot.usage_final)
            outcome.pairs_usage.add(slot.usage);
        outcome.entries.push_back(std::move(slot.entry));
        slot.entry = Json::Value();
    }
    state->inline_gone = nullptr;
    return outcome;
}

size_t fanouts_alive() { return g_fanouts_alive.load(); }

} // namespace cli
} // namespace mtg
