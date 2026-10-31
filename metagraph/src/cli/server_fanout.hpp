#ifndef __METAGRAPH_SERVER_FANOUT_HPP__
#define __METAGRAPH_SERVER_FANOUT_HPP__

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>

#include "common/threads/thread_usage.hpp"
#include "server_selection.hpp"


namespace mtg {
namespace cli {

// The fan-out of a multi-graph /resolve or /pattern over its selected pairs, its deadlines and
// the registry of requests by request_id (SPEC-multigraph-envelope.md sections 4 and 5), apart
// from the server and without the HTTP library, so that the lifecycle of the workers is tested
// on its own: a worker may outlive the HTTP answer of its request (it was cut at the work
// cutoff, or its client left), so everything it reads is held by it, never borrowed from the
// handler's stack

using SteadyClock = std::chrono::steady_clock;
using WallClock = std::chrono::system_clock;

// ---------------------------------------------------------------- the request's handles

/**
 * The request-level fields of section 4 and 5 that /resolve and /pattern read before anything
 * runs: request_id (attempt_id's token rule, valid_attempt_id), not_after_ms (the latest START,
 * as /traverse's) and finish_by_ms (the completion deadline), each absent or as given. Read by
 * the routes' strict parsers too, which refuse any other spelling.
 */
struct RequestHandles {
    std::optional<std::string> request_id;
    std::optional<uint64_t> not_after_ms;
    std::optional<uint64_t> finish_by_ms;
};
// The three fields of |request| (an object; anything else: none of them), checked: throws
// InvalidRequest (traverse.hpp) naming request.<field> for a request_id that is not a token, or
// a not_after_ms or finish_by_ms that is not an integer in [0, kMaxNotAfterMs]
RequestHandles read_request_handles(const Json::Value &request);

// ---------------------------------------------------------------- the deadlines

/**
 * The deadlines of one request (section 4), fixed when its handler starts:
 *  - |effective|: min(arrival + the route's request budget, finish_by_ms - the clock skew
 *    allowance, arrival + the transport cap): the deadline the envelope's writer enforces;
 *  - |cutoff|: |effective| less the route's finalisation reserve: the fan-out stops admitting
 *    and waiting for pairs there;
 *  - |finish|: finish_by_ms - the skew allowance (none without finish_by_ms): the bound of each
 *    pair's own deadline and of an in_ram memory wait.
 * Each instant on both clocks: the steady one the waits read, the wall one the answer states.
 */
struct RequestDeadlines {
    std::optional<uint64_t> not_after_ms;
    std::optional<uint64_t> finish_by_ms;
    SteadyClock::time_point arrival;
    WallClock::time_point arrival_wall;
    SteadyClock::time_point effective;
    WallClock::time_point effective_wall;
    SteadyClock::time_point cutoff;
    WallClock::time_point cutoff_wall;
    std::optional<SteadyClock::time_point> finish;
    std::optional<WallClock::time_point> finish_wall;
    // the effective deadline is the transport's cap (no shorter budget, no finish_by_ms): a pair
    // not started by the cutoff is not_reached with reason "transport"
    bool transport_bound = false;

    // the steady instant of |wall| (through the arrival's two readings)
    SteadyClock::time_point steady_of(WallClock::time_point wall) const;
    WallClock::time_point wall_of(SteadyClock::time_point steady) const;
    // the envelope's `deadline` (deadline_json): not_after_ms and finish_by_ms as given, and
    // effective_deadline_at and cutoff_at stated when the request gave either of them (null
    // otherwise, as the reference samples state a request without them: no deadline given)
    Json::Value json() const;
};

/**
 * The deadlines of a request that arrived at |arrival_wall| (its header read) as |now_wall| and
 * |now| read the clocks when its handler started: |budget_ms|, the route's request budget
 * (/pattern's time_budget_ms, /resolve's bounds.time_budget_ms or the server's maximum), |
 * |skew_ms|, the clock skew allowance, |reserve_ms|, the route's finalisation reserve, and
 * |transport_ms|, the transport's cap from the arrival (kServerMaxDeadlineMs).
 */
RequestDeadlines request_deadlines(WallClock::time_point arrival_wall,
                                   WallClock::time_point now_wall,
                                   SteadyClock::time_point now,
                                   double budget_ms,
                                   const RequestHandles &handles,
                                   uint64_t skew_ms,
                                   double reserve_ms,
                                   double transport_ms);

// ---------------------------------------------------------------- the request registry

// Where a pair of a request is (GET /requests/{request_id}'s pairs[].state)
enum class PairState { QUEUED, WAITING_MEMORY, LOADING, RUNNING, ENDED, NEVER_STARTED };
const char* to_string(PairState state);

/**
 * One request as GET /requests/{request_id} states it (section 5): its pairs' states and
 * usage, and its response once written. Thread-safe: the handler and every worker of the
 * request update it, a GET reads it. It never holds a result (results are not retained). The
 * request has ENDED once its response was recorded (written, or nothing written) and every
 * pair's worker ended (or the pair never started); then |on_ended| is called once.
 */
class RequestEntry {
  public:
    struct Clocks {
        std::function<SteadyClock::time_point()> steady;
        std::function<WallClock::time_point()> wall;
    };

    // |pairs|: {graph, graph_path, annotation_path} each (null tags on a single-graph server)
    RequestEntry(std::optional<std::string> request_id, std::string route,
                 std::string server_instance, WallClock::time_point received_wall,
                 std::vector<Json::Value> pairs, Clocks clocks = {});

    const std::optional<std::string>& request_id() const { return request_id_; }
    const std::string& route() const { return route_; }
    size_t num_pairs() const { return pairs_.size(); }

    // ---- the pairs (the fan-out's workers, or the handler of a single-graph request)
    // waiting for the memory to load its index (not started)
    void pair_waiting_memory(size_t i);
    // started: loading its index into RAM for the request (|bytes| reserved), or working
    void pair_loading(size_t i, uint64_t bytes);
    void pair_running(size_t i);
    bool pair_started(size_t i) const;
    PairState pair_state(size_t i) const;
    // the pair's CPU time so far (ms), sampled by its worker while it works
    void pair_cpu(size_t i, double cpu_ms);
    // answered | refused | not_reached | unfinished
    void pair_outcome(size_t i, const std::string &outcome);
    // the pair's worker ended: ENDED with its final usage (|memory_bytes|, |cpu_ms| when
    // known), or NEVER_STARTED when it never started
    void pair_ended(size_t i, std::optional<uint64_t> memory_bytes = std::nullopt,
                    std::optional<double> cpu_ms = std::nullopt);
    // a pair that will never start (the request closed while it was queued or waited for
    // memory): NEVER_STARTED now, whether or not its worker still runs (it ends without work)
    void pair_never_started(size_t i);
    // the usage of a started pair so far, as an `unfinished` refusal states it: {state
    // (|state|: "stopping" or "loading"), started_at, elapsed_ms, cpu_ms, memory_bytes}
    Json::Value pair_usage_now(size_t i, const std::string &state) const;

    // ---- the response: |status| (0: nothing written, the client is gone) and its bytes
    void response(int status, std::optional<size_t> bytes);
    bool ended() const;
    std::optional<WallClock::time_point> finished_at() const;

    // GET /requests/{request_id}: {request_id, server_instance, route, state, received_at,
    // finished_at, reason, response, retained_until, pairs}; |retention_s| sets
    // retained_until (finished_at + it) once ended
    Json::Value json(uint64_t retention_s) const;

    // called once, without the entry's lock held, when the request ended
    void set_on_ended(std::function<void()> on_ended);

  private:
    struct Pair {
        Json::Value ref;
        PairState state = PairState::QUEUED;
        std::string outcome;
        std::optional<SteadyClock::time_point> started;
        std::optional<WallClock::time_point> started_wall;
        std::optional<SteadyClock::time_point> finished;
        double cpu_ms = 0;
        bool has_cpu = false;
        uint64_t memory_bytes = 0;
    };
    // the end, when it is reached now (under |mutex_|): the on_ended to call after unlocking
    std::function<void()> end_locked();
    Json::Value usage_locked(const Pair &p) const;
    WallClock::time_point wall_of_locked(SteadyClock::time_point t) const;

    const std::optional<std::string> request_id_;
    const std::string route_;
    const std::string server_instance_;
    const WallClock::time_point received_wall_;
    SteadyClock::time_point received_;
    Clocks clocks_;
    mutable std::mutex mutex_;
    std::vector<Pair> pairs_;
    bool responded_ = false;
    int status_ = 0;
    std::optional<size_t> bytes_;
    std::optional<WallClock::time_point> finished_wall_;
    std::string reason_;
    std::function<void()> on_ended_;
    bool ended_called_ = false;
};

// The registry's settings (--request-retention-s, --request-retention)
struct RequestSettings {
    // ended requests are kept this long (s), at most this many (oldest dropped first); 0:
    // none kept after it ended (its id runs again)
    uint64_t retention_s = 3600;
    size_t retention_count = 10'000;
    std::string server_instance;
    RequestEntry::Clocks clocks;
};

/**
 * The requests of one server process by request_id (section 5): running ones, and ended ones
 * for the retention period. A request_id runs once while it is running or retained: a repeat
 * starts no work (409, the entry and an error). For /traverse the attempt_id is the handle and
 * the attempt registry answers; this registry keeps only the pair each attempt ran on (for
 * the pairs of GET /requests), as long as it keeps a request.
 */
class RequestRegistry {
  public:
    explicit RequestRegistry(RequestSettings settings);

    const RequestSettings& settings() const { return shared_->settings; }

    // Registers |entry| (with a request_id) unless its id is running or retained: then the
    // 409 body, the other entry's json with `error`, and nothing is registered. |admission|,
    // when given, is asked after that check and before the entry is registered (under the
    // registry's lock): a body it returns (an expired not_after_ms) refuses the request
    // instead, nothing registered, so that such a refusal always means that no request with
    // the id is held here (the attempt registry's order)
    std::optional<Json::Value> start(
            const std::shared_ptr<RequestEntry> &entry,
            const std::function<std::optional<Json::Value>()> &admission = nullptr);
    // GET /requests/{request_id}: (200, the entry) or (404, {error, request_id,
    // server_instance}) for an id that is not held (never received, or past retention)
    std::pair<int, Json::Value> state(const std::string &request_id) const;
    // the retention rule, stated by the 404 and the 409
    std::string retention_text() const;
    // the 404 body of an unknown id
    Json::Value unknown_json(const std::string &request_id) const;

    // /traverse: the pair the attempt |attempt_id| ran on ({graph, graph_path,
    // annotation_path}), kept as requests are (by the same retention, from now)
    void note_traverse_pair(const std::string &attempt_id, Json::Value pair);
    std::optional<Json::Value> traverse_pair(const std::string &attempt_id) const;

    // the requests held now (running and retained), for tests
    size_t size() const;

  private:
    struct Shared {
        RequestSettings settings;
        mutable std::mutex mutex;
        std::map<std::string, std::shared_ptr<RequestEntry>> entries;
        // ended entries, oldest first: (when it ended, id)
        std::deque<std::pair<SteadyClock::time_point, std::string>> ended;
        std::map<std::string, Json::Value> traverse_pairs;
        std::deque<std::pair<SteadyClock::time_point, std::string>> traverse_order;
        SteadyClock::time_point now() const;
        void expire_locked();
    };
    std::shared_ptr<Shared> shared_;
};

// GET /requests/{attempt_id} of a /traverse attempt (sample 21): the request object built from
// the attempt's state (GET /traverse/attempt/{attempt_id}, |attempt|), its one pair (|pair|,
// null tags when unknown: a single-graph server) and |retention_s| for retained_until
Json::Value traverse_request_json(const Json::Value &attempt, const Json::Value &pair,
                                  uint64_t retention_s);

// ---------------------------------------------------------------- the fan-out

// What a pair's job hands back (PairJob)
struct PairResult {
    // the envelope entry: envelope_entry_answered or envelope_entry_refused
    Json::Value entry;
    bool answered = false;
    // the job never started the pair (PairControl::start_* refused: the request closed, or
    // not_after_ms passed): the fan-out states the not_reached refusal
    bool not_started = false;
    // the work was abandoned (the request's stop: its client left or the cutoff passed): no
    // entry; the fan-out states what it decided
    bool aborted = false;
    // the pair's final account (the answer's work.memory_bytes, or what it held)
    std::optional<uint64_t> memory_bytes;
    // the job's CPU time and page faults on its thread (from the job's start)
    common::ThreadUsage usage = common::ThreadUsage::zero();
};

class FanoutState;

/**
 * What a pair's job reads from and reports to the fan-out, on the job's thread, for the job's
 * duration. Its predicates hold the request's shared state, so a copy given to the work stays
 * valid however long the work runs.
 */
class PairControl {
  public:
    PairControl(std::shared_ptr<FanoutState> state, size_t index);

    size_t index() const { return index_; }
    const SelectedPair& pair() const;
    const RequestDeadlines& deadlines() const;
    bool inline_run() const;

    // whether the work must stop: the request's stop flag (the cutoff passed, the client is
    // gone, the server stops), or for an inline pair its client check. Samples the thread's
    // CPU time for the registry (at most every 100 ms)
    bool stopped();
    // stopped() as a predicate for the work (holds the shared state)
    std::function<bool()> stop_predicate();
    // the predicate an in_ram memory wait polls: stopped(), not_after_ms passed (a pair not
    // started by then never starts), or the wait's limit passed (the cutoff in a pooled
    // fan-out; finish_by_ms - the skew allowance, or the transport's cap, inline)
    std::function<bool()> wait_predicate();

    // the pair waits for the memory to load its index into RAM: not started yet
    void waiting_memory();
    // The pair starts: loading its index (|bytes| reserved for it), or its work. False when it
    // may not (the request closed, or not_after_ms passed): the job releases what it holds
    // and returns not_started. Once started, the pair is unfinished if it is still running at
    // the cutoff. start_running after start_loading only notes the end of the load
    bool start_loading(uint64_t bytes);
    bool start_running();

    // The cap of the pair's own deadline (its route's deadline, which starts when it starts,
    // becomes the minimum of its own and this): finish_by_ms - the skew allowance, none
    // without it. The work cutoff caps nothing: a pair still working there is unfinished
    std::optional<SteadyClock::time_point> deadline_cap() const;

  private:
    std::shared_ptr<FanoutState> state_;
    size_t index_;
    common::ThreadUsage start_usage_;
    SteadyClock::time_point next_sample_;
};

// The job of one pair: answers it as one single-graph request would be (its entry), on the
// thread it is called on. Must not throw (the fan-out turns an escaping exception into a 500
// refusal); must not read the handler's stack (it may outlive it): it holds what it reads
using PairJob = std::function<PairResult(PairControl &)>;
// Runs a task on another thread (the graphs' pool)
using Launch = std::function<void(std::function<void()>)>;

struct FanoutSettings {
    // "/resolve" or "/pattern": the messages name the route ("resolve: pair (...) ...")
    std::string route;
    std::vector<SelectedPair> pairs;
    // per pair: its index_fp ("" none)
    std::vector<std::string> index_fps;
    RequestDeadlines deadlines;
    // the handler's check while it waits (and an inline pair's): the client is gone or the
    // server stops; asked on the handler's thread only, every |poll_ms|
    std::function<bool()> gone;
    uint64_t poll_ms = 100;
    // the clocks (tests inject them)
    std::function<SteadyClock::time_point()> clock;
    std::function<WallClock::time_point()> wall_clock;
};

// What the fan-out answered: the entries in pair order, and the request's state
struct FanoutOutcome {
    std::vector<Json::Value> entries;
    // the client left (or the server stops) before the answer: nothing is written
    bool client_gone = false;
    // the completed pairs' CPU time and page faults
    common::ThreadUsage pairs_usage = common::ThreadUsage::zero();
};

/**
 * Runs |job| on every pair of |settings| (section 4): one pair inline on the calling thread,
 * two or more each as a task of |launch|. The request's complete state is one object shared by
 * the handler and every worker; each pair's result is published once under its lock, and the
 * handler reads a snapshot. The handler waits until every pair is done, its client is gone, or
 * the work cutoff: then pairs never started are not_reached (reason cutoff, not_after_ms or
 * transport), pairs started and not done are unfinished with their usage so far and told to
 * stop; their workers end later and |entry| (the registry's, always given) records their
 * final usage. A pair that a worker takes after not_after_ms is not_reached (not_after_ms).
 */
FanoutOutcome run_fanout(const FanoutSettings &settings,
                         const std::shared_ptr<RequestEntry> &entry,
                         PairJob job,
                         const Launch &launch);

// For tests: the shared state of the fan-outs whose workers have not all ended
size_t fanouts_alive();

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_FANOUT_HPP__
