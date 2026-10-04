#ifndef __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__
#define __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>

#include "graph/traversal/walker.hpp"


namespace mtg {
namespace cli {

/**
 * The backend half of stage 4 of DESIGN-traverse-graphlet.md §14 (the ledger itself lives in
 * the search service): a /traverse request that carries an `attempt_id` is an ATTEMPT the
 * service's ledger reserved an allowance for. The server registers it while it runs, lets it
 * be cancelled by id (POST /traverse/cancel), reports its state (GET /traverse/attempt/{id})
 * for a retention period after it finished, states its usage in every response to it, and
 * enforces a duration bound on it itself — §14 v5.1: a lease is a valid release point only
 * because the backend enforces the same deadline, so an attempt cannot outlive its lease, and
 * a cancellation the backend has not acknowledged releases nothing.
 *
 * Every /traverse (with or without attempt_id) is also stopped when its client is gone (the
 * connection closed): the walk is abandoned and nothing is written.
 */

// attempt_id, budget_id, locus_id: ^[A-Za-z0-9._:-]{1,128}$ (a token in logs, URLs and keys)
bool valid_attempt_id(const std::string &id);

struct AttemptIds {
    std::string attempt_id;     // "" = the request is not ledger-managed
    std::string budget_id;      // echoed only
    std::string locus_id;       // echoed only
};

// The three request fields of /traverse (top level). Throws InvalidRequest when one is not a
// string matching the pattern, or when budget_id or locus_id is given without attempt_id.
AttemptIds attempt_ids(const Json::Value &request);

// The server's settings for attempts (Config::traverse_attempt_*)
struct AttemptSettings {
    // added to n_seeds x the effective per-seed time budget: the time to read the seeds,
    // build the response and write it (ms)
    double allowance_ms = 10'000;
    // the bound never exceeds this (ms, on the attempt's clock, which starts when the server
    // read the request's header): the content timeout of the HTTP server, after which it shuts
    // the connection, less one second for the response to leave. 0: none (the CLI)
    double hard_cap_ms = 899'000;
    // finished attempts are kept this long, and at most this many (oldest dropped first)
    uint64_t retention_s = 3600;
    size_t retention_count = 10'000;
    // the client's connection is checked at most this often (ms)
    uint64_t client_check_ms = 100;
    // the clock (the bound, the client check's interval) is read at every |poll_stride|-th poll
    // of the walk only: a poll comes before every head, and the stop flag alone is read at every
    // one (tests set 1)
    uint32_t poll_stride = 8;
    // the clock (tests inject one)
    std::function<std::chrono::steady_clock::time_point()> clock;
};

// What one seed of the attempt consumed (the usage block's per_seed)
struct SeedUsage {
    // complete | partial | failed | not_started, or abandoned (the client went away before its
    // walk ended: no response is written, so no per_seed states it; the attempt counts it in
    // seeds.abandoned, not in seeds.finished)
    std::string outcome = "not_started";
    // what stopped it, when a resource did: memory | work | time | cancelled | attempt_deadline
    std::string stopped_by;
    // The walker's meter, as the result states it: for a failed or never started seed, whose
    // result is not the walk's, |memory_final| is what that result holds until the response is
    // written (its fixed part and its echo of seed_id, as failed_soft prices it) and
    // |soft_excess| what its memory_bound_soft states (the echo included; review of the stage-4
    // backend, F1); |memory_peak| is at most the budget (what the budget admitted)
    graph::traversal::AttemptMeter meter;
    // the demand a memory budget refused when it stopped or failed the seed (bytes at that
    // budget: a head, a read, or the seed's depth-0 state); none when no memory budget did
    std::optional<uint64_t> refused_bytes;
    double elapsed_ms = 0;
    bool started = false;
    bool walked = false;
};

// The attempt reached its bound while its response was built or written: nothing more is
// written than what the bound allows (the server answers 503 with the usage)
class AttemptAtBound : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

/**
 * One /traverse request as the server controls it. The handler's thread drives it (poll,
 * seeds, bound, usage); POST /traverse/cancel and GET /traverse/attempt/{id} reach it from
 * other threads through the registry: the stop flag is atomic (the first reason wins), the
 * state they read is under the attempt's mutex.
 */
class Attempt {
  public:
    using Clock = std::chrono::steady_clock;

    // |header_read|: when the server read the request's header (the attempt's clock starts
    // there: the HTTP server's content timeout runs from then). |enforced|: the server enforces
    // the bound (the CLI states it only). |peer_gone|: whether the client is gone (null: never).
    Attempt(size_t request_id, std::chrono::system_clock::time_point header_read,
            const AttemptSettings &settings, std::string server_instance, bool enforced,
            std::function<bool()> peer_gone = nullptr);

    // ---- identity (set before the attempt is registered, then fixed)
    void set_ids(AttemptIds ids) { ids_ = std::move(ids); }
    const AttemptIds& ids() const { return ids_; }
    // the request carries an attempt_id: it is registered, its bound enforced (by the server)
    // and every response to it carries usage
    bool managed() const { return !ids_.attempt_id.empty(); }
    size_t request_id() const { return request_id_; }

    // ---- the attempt's clock: ms since the server read the request's header
    double elapsed_ms() const;

    // ---- the bound: min(n_seeds x T + allowance, hard cap), T the effective per-seed time
    // budget (0 when it is not positive). The seeds stop being walked at
    // bound - allowance / 2, leaving the rest of the allowance to deliver what was walked.
    // Set once the request is parsed; before that the hard cap alone bounds the attempt.
    // |memory_budget|: the request's bounds.max_memory_mb in bytes, 0 for none (the soft
    // excess is observed under one, and what was held is bounded by it)
    void set_bound(size_t seeds, double time_budget_ms, uint64_t memory_budget = 0);
    double bound_ms() const { return bound_ms_; }
    double walk_until_ms() const { return walk_until_ms_; }
    bool enforced() const { return enforced_ && managed(); }

    // ---- stops
    // asks the walk to stop; false when another reason was first (the first one wins)
    bool request_stop(graph::traversal::ExternalStop reason);
    // POST /traverse/cancel: records when it was asked, then request_stop(CANCELLED)
    bool cancel();
    // the handler's poll, at the walker's checkpoints and between seeds: the stop flag, and at
    // every poll_stride-th poll (or with |force|: between seeds) the bound (enforced only) and
    // the client's connection (at most every client_check_ms). Records the first poll that
    // returns a stop: when the walk stopped.
    graph::traversal::ExternalStop poll(bool force = false);
    // while the response is built and written: the client's connection (AttemptAborted) and
    // the bound itself (AttemptAtBound); a cancel does not discard a finished walk
    void check_delivery();
    // the walk of every seed is over (normally, or at a stop): stopped_at, if not yet
    void walk_ended();
    // "cancelled" | "deadline" when a stop reached the walk, "completed" otherwise
    const char* walk_reason() const;

    // ---- seeds (the handler's thread)
    void seed_started(size_t index);
    // the seed's walk is over: its meter, its outcome so far and what stopped it (called again
    // when its result turns out not to be the walk's: a refused seed); |usage.outcome|
    // abandoned: the client went away before the walk ended (seeds.abandoned)
    void seed_walked(size_t index, const SeedUsage &usage);
    // a seed the attempt never started (a stop came first): what its failed result states
    // (|usage.meter|'s memory_final and soft_excess); neither started nor finished
    void seed_not_started(size_t index, const SeedUsage &usage);
    // the seed's outcome once its result is built, and its elapsed time with the building
    void seed_delivered(size_t index, const std::string &outcome, double elapsed_ms);

    // ---- statements
    // the response-level `usage` block, |reason| completed | cancelled | deadline | error
    Json::Value usage_json(const std::string &reason, bool per_seed = true) const;
    // the attempt's state (GET /traverse/attempt/{id})
    Json::Value state_json() const;
    // for the server's log: when the walk stopped, how long after the stop was requested, and
    // the seeds walked
    std::string stop_summary() const;

    // ---- the end (the registry's; also when the handler returned without registering)
    // |reason| completed | cancelled | client_gone | deadline | error; |status| 0 and |bytes|
    // nullopt: nothing was written
    void finish(const std::string &reason, int status, std::optional<size_t> bytes);
    bool finished() const;
    // waits until the attempt finished, at most |ms|; true when it did
    bool wait_finished(uint64_t ms) const;
    // tombstone: a cancel named this id before any request with it arrived
    bool tombstone() const { return tombstone_; }
    static std::shared_ptr<Attempt> make_tombstone(const std::string &id,
                                                   const AttemptSettings &settings,
                                                   std::string server_instance);
    Clock::time_point finished_at() const;

  private:
    Clock::time_point now() const;
    std::string iso(Clock::time_point t) const;
    double ms_since_received(Clock::time_point t) const;
    Json::Value bound_json() const;
    // {requested, started, finished, abandoned}, under |mutex_|
    Json::Value seeds_json() const;

    const size_t request_id_;
    AttemptSettings settings_;
    const std::string server_instance_;
    const bool enforced_;
    std::function<bool()> peer_gone_;
    AttemptIds ids_;
    // the attempt's clock: |received_| on the steady clock, |received_wall_| when that was
    Clock::time_point received_;
    std::chrono::system_clock::time_point received_wall_;

    // the bound (written by the handler's thread under |mutex_|, read by it without)
    size_t bound_seeds_ = 0;
    double bound_time_budget_ms_ = 0;
    double bound_ms_ = 0;
    double walk_until_ms_ = 0;
    bool bound_set_ = false;
    bool capped_ = false;
    uint64_t memory_budget_ = 0;

    std::atomic<uint8_t> stop_ { 0 };
    // the handler's thread only: the next instant the client's connection is checked, and
    // whether it was found gone (whatever stopped the walk first)
    Clock::time_point next_peer_check_;
    bool client_gone_ = false;
    uint32_t polls_ = 0;

    mutable std::mutex mutex_;
    mutable std::condition_variable finished_cv_;
    std::optional<Clock::time_point> cancel_requested_at_;
    std::optional<Clock::time_point> stop_requested_at_;
    std::optional<Clock::time_point> stopped_at_;     // when the walk stopped
    std::optional<Clock::time_point> finished_at_;    // when the handler returned
    bool stop_seen_ = false;                         // a poll returned the stop to the walk
    std::vector<SeedUsage> seeds_;                    // freed when the attempt finishes
    size_t seeds_requested_ = 0;
    size_t seeds_started_ = 0;
    size_t seeds_walked_ = 0;
    size_t seeds_abandoned_ = 0;                      // walks the client's departure cut
    std::string reason_;                              // set by finish()
    int status_ = 0;
    std::optional<size_t> bytes_;
    Json::Value final_usage_;                         // the usage without per_seed, at finish
    bool tombstone_ = false;
};

/**
 * The attempts of one server process: the running ones by attempt_id, and the finished ones
 * (and tombstones) for the retention period. An attempt_id is unique per server process while
 * it is running or retained: a second request with it is refused (409), never run twice.
 */
class AttemptRegistry {
  public:
    explicit AttemptRegistry(AttemptSettings settings);

    const AttemptSettings& settings() const { return settings_; }
    // random per process (16 hex digits): a ledger tells a restarted backend (whose attempts
    // are gone) from one that forgot an attempt after the retention period
    const std::string& server_instance() const { return instance_; }

    // registers a managed attempt; when its id is running or retained (or tombstoned),
    // nothing is registered and the other attempt's state is returned (the 409 body's)
    std::optional<Json::Value> start(const std::shared_ptr<Attempt> &attempt);
    // POST /traverse/cancel {attempt_id, wait_ms}: (HTTP status, body). An unknown id is
    // tombstoned for retention_s (404, tombstone: true); when retention_count tombstones are
    // held already it is not, and the cancel is refused (429, tombstone: false): the promise
    // that the id will not run here is never given for less than retention_s
    std::pair<int, Json::Value> cancel(const std::string &id, uint64_t wait_ms);
    // GET /traverse/attempt/{id}: (HTTP status, body)
    std::pair<int, Json::Value> state(const std::string &id);
    // the handler returned: the attempt is finished, and kept for the retention period
    void finish(const std::shared_ptr<Attempt> &attempt, const std::string &reason, int status,
                std::optional<size_t> bytes);
    // the retention rule, stated by 404s and GET /traverse/capabilities
    std::string retention_text() const;

  private:
    void expire_locked();
    Attempt::Clock::time_point now() const;
    Json::Value unknown_json(const std::string &id) const;

    AttemptSettings settings_;
    std::string instance_;
    mutable std::mutex mutex_;
    std::map<std::string, std::shared_ptr<Attempt>> attempts_;
    // finished attempts, oldest first: kept retention_s, at most retention_count of them
    std::deque<std::shared_ptr<Attempt>> retained_;
    // tombstones, oldest first, apart from the finished attempts: each is kept the full
    // retention_s, whatever finishes after it (review of the stage-4 backend, F6: sharing
    // retention_count, later finishes evicted a tombstone early and the cancelled id ran)
    std::deque<std::shared_ptr<Attempt>> tombstones_;
};

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__
