#ifndef __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__
#define __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <limits>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>

#include "graph/traversal/walker.hpp"


namespace mtg {
namespace cli {

/**
 * The backend half of the attempt ledger of DESIGN-traverse-graphlet.md §14 (the ledger
 * itself lives in the search service): a /traverse request that carries an `attempt_id` is an
 * ATTEMPT the service's ledger reserved an allowance for. The server registers it while it
 * runs, lets it be cancelled by id (POST /traverse/cancel), reports its state (GET
 * /traverse/attempt/{id}) for a retention period after it finished, states its usage in every
 * response to it, and enforces a duration bound on it itself — §14: a lease is a valid
 * release point only because the backend enforces the same deadline, so an attempt cannot
 * outlive its lease, and a cancellation the backend has not acknowledged releases nothing.
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
    // The instant (Unix epoch, ms, on the server's clock) after which the request must not be
    // started: a ledger that lost track of an attempt treats it as one that cannot start
    // subsequently once its own clock passes this plus the server's stated clock skew
    // allowance (an unanswered request may already be running: the attempt's bound covers
    // that), so the server refuses it at handler start once it has passed (409, state
    // expired) rather than start work the ledger has already released. With or without
    // attempt_id; echoed in usage
    std::optional<uint64_t> not_after_ms;
    // The server_instance the request is meant for (with attempt_id only): a process with
    // another one refuses it before anything runs (409, state instance_mismatch). Tombstones
    // live in memory, so after a restart a delayed copy of a cancelled request would find
    // none and run: with this field it is refused instead
    std::string expect_server_instance;
};

// The names of AttemptIds' request fields, in the order the capabilities list them
constexpr const char *kAttemptFields[] = { "attempt_id", "budget_id", "locus_id",
                                           "not_after_ms", "expect_server_instance" };

// The largest not_after_ms accepted: 2^53 - 1, the largest integer every JSON reader keeps
// exactly (a ledger in JavaScript or Python writes and reads it unchanged)
constexpr uint64_t kMaxNotAfterMs = (uint64_t(1) << 53) - 1;

// The request fields of /traverse (top level) the server reads before the request is parsed.
// Throws InvalidRequest when an id is not a string matching the pattern, when budget_id,
// locus_id or expect_server_instance is given without attempt_id, or when not_after_ms is not
// an integer in [0, kMaxNotAfterMs].
AttemptIds attempt_ids(const Json::Value &request);
// not_after_ms as /traverse and POST /traverse/cancel read it: absent (nullopt), or an integer
// in [0, kMaxNotAfterMs]; throws InvalidRequest naming |where| otherwise
std::optional<uint64_t> read_not_after_ms(const Json::Value &request, const std::string &where);

// Whether a request carrying |ids| names another server_instance than |instance| (it must then
// be refused before anything runs), and the 409 body: {error, state: "instance_mismatch",
// expect_server_instance, server_instance, attempt_id/budget_id/locus_id as given,
// not_after_ms when given}
bool instance_mismatch(const AttemptIds &ids, const std::string &instance);
Json::Value instance_mismatch_json(const AttemptIds &ids, const std::string &instance);

// Whether a request carrying |ids| must be refused at handler start: its not_after_ms is
// earlier than |now_ms| (the server's clock, Unix epoch ms; strict: no allowance is added,
// the allowance is the ledger's to add). The 409 body when it must: {error, state: "expired",
// not_after_ms, server_time_ms, attempt_id/budget_id/locus_id as given, server_instance}
bool not_after_passed(const AttemptIds &ids, uint64_t now_ms);
Json::Value expired_json(const AttemptIds &ids, uint64_t now_ms,
                         const std::string &server_instance);
// |t| as ISO-8601 UTC with milliseconds (the instants of usage and the attempt's state)
std::string iso_utc(std::chrono::system_clock::time_point t);

// The server's settings for attempts (Config::traverse_attempt_*)
struct AttemptSettings {
    // added to n_seeds x the effective per-seed time budget: the time to read the seeds,
    // build the response and write it (ms)
    double allowance_ms = 10'000;
    // the bound never exceeds this (ms, on the attempt's clock, which starts when the server
    // read the request's header): the content timeout of the HTTP server, after which it shuts
    // the connection, less one second for the response to leave. 0: none (the CLI)
    double hard_cap_ms = 899'000;
    // finished attempts are kept this long, and at most this many (oldest dropped first) —
    // one sent with not_after_ms is held longer, until not_after_ms + clock_skew_ms (within
    // the cap), among the tombstones —; a cancel of an unknown id tombstones it at least
    // retention_s, at most retention_count tombstones at once. 0: nothing is kept (no
    // tombstone: such a cancel is refused, 429; no finished attempt is held)
    uint64_t retention_s = 3600;
    size_t retention_count = 10'000;
    // the longest a tombstone is held to cover a not_after_ms (+ clock_skew_ms) named by a
    // cancel or by a refused copy of the request, and a finished attempt's id from its finish
    // (s; never less than retention_s)
    uint64_t tombstone_max_s = 86'400;
    // the client's connection is checked at most this often (ms)
    uint64_t client_check_ms = 100;
    // the clock (the walk-until, the client check's interval) is read at every
    // |poll_stride|-th poll of the walk only, and at every forced poll: a poll comes before
    // every head, and the stop flag alone is read at every one (tests set 1). The bound itself
    // is compared only by check_delivery
    uint32_t poll_stride = 8;
    // the clock (tests inject one)
    std::function<std::chrono::steady_clock::time_point()> clock;
    // the wall clock not_after_ms is compared with (tests inject one; system_clock by default)
    std::function<std::chrono::system_clock::time_point()> wall_clock;
    // What a ledger must add to not_after_ms before it treats an unanswered attempt as one that
    // cannot start subsequently (stated in capabilities; the server's own not_after check is
    // strict, and a tombstone covering a not_after_ms is held to it plus this): how far this
    // server's clock may be from the ledger's (NTP keeps them closer; ms)
    uint64_t clock_skew_ms = 2000;
    // the HTTP server's content timeout (s), from which hard_cap_ms is derived: stated
    uint64_t content_timeout_s = 900;
    // The delivery reserve: the rates (MB/s) at which the response is assumed to be compressed
    // and its results' text built — the latter replaced by the slowest rate measured on the
    // attempt's own seeds of at least kMeasuredTextBytes — and how many bytes of the walker's
    // modelled account one byte of a seed's text is at most taken for (the model prices every
    // delivered byte several times over). The ratios are starting estimates just below the
    // smallest measured on real responses — 33.5 (UHGG full) to 1,344 (SRA summary) for the
    // JSON details, 58.4 and more for a graphlet —, which a server replaces by its own
    // measurements. Conservative: a server's first large attempt of a detail can still be cut
    // early until it measured that detail (a warm SRA server cut the first tree and full
    // attempts of a 16S beam at 3-5 s of 40 s with 20/40 and with 30/50 alike, the tree
    // measuring 115-129)
    double delivery_compress_mbps = 50;
    double delivery_build_mbps = 10;
    double account_per_text_byte_json = 30;
    double account_per_text_byte_graphlet = 50;
    // The walk does not stop at its walk-until but at the first poll that reads the clock after
    // it — after a chunk of an annotation read, a lookahead's poll, the heads between two
    // readings of the clock (one poll in poll_stride reads it) — and its stopped
    // seed is finalised before its text is built: the time from the walk-until to the walk's
    // end assumed until the server measured a longer one (ms; the server sets chunk_target_ms +
    // 950: a walk stopped 124 ms after its walk-until with none assumed would leave its
    // delivery that much short of the reserve, 503). Calibrated: SRA attempts stopped 352 ms
    // after their walk-until on a quiet fresh server, 1,001 ms on a loaded one (1,699 ms once
    // under load, on a 400 MB result): the middle one is assumed, a longer one measured
    // replaces it; a finalisation that grows with the result is covered by the reserve's margin
    // as long as it runs faster than 4 x build_mbps (about 235 MB/s for that 400 MB result,
    // against 40 MB/s needed)
    double delivery_stop_ms = 1000;
};

// The delivery model's rates and ratios vary between responses: the reserve keeps this much
// more than the model's time (a reserve with no margin would deliver a response the model
// fits exactly only about half the time)
constexpr double kReserveMargin = 1.25;

// The most text one byte of the record coordinates' share of a seed's account (DeliveryCosts:
// coordinate_fixed, coordinate_run, coordinate_seed, occurrence) can write, inverted: every
// byte of coordinate text costs at least this many account bytes, whatever the digits — an
// occurrence's widest text (two 20-digit numbers, 44 bytes with its brackets and comma) is
// priced 872, a run's entry (221 bytes at its widest) 3,680, a seed label's (79) 1,728, the
// block's skeleton with a cut list's limitation, K record and drop_coordinates (960 bytes in a
// graphlet) 15,151, the null form (106) 1,584: 14.9 at the least
// (GraphletCoordinates.CoordinateAccountBoundsItsText). The
// reserve estimates the coordinate share's text with it, not with the measured ratio of the
// rest of the output (115-129 for a tree, against 20-73 for an occurrence: the measured
// ratio would understate a coordinate-heavy seed's text up to 4 times, DESIGN §18.8, M2)
constexpr uint64_t kCoordinateAccountPerTextByte = 12;

// a seed's text (or a response) of at least this many bytes measures a delivery rate (smaller
// ones are dominated by fixed costs)
constexpr uint64_t kMeasuredTextBytes = 1 << 20;
// the delivery rates the server measured on its own responses: the slowest of the last this
// many measurements replaces the configured rate (a machine slower than assumed is seen at
// once, a temporarily slow response is remembered for a while)
constexpr size_t kRateWindow = 16;

// What a server measured of its own deliveries (0: not measured yet): the slowest
// rates (MB/s) at which it built a seed's result text and compressed a response, and per
// detail the smallest ratio of a seed's modelled account to its text — each over its last
// kRateWindow measurements of at least kMeasuredTextBytes
struct DeliveryMeasurements {
    double build_mbps = 0;
    double compress_mbps = 0;
    std::map<std::string, double> account_per_text_byte;    // by output.detail
    // the longest time from the walk-until to the walk's end (ms) of its last kRateWindow
    // attempts that walked past their walk-until
    double stop_ms = 0;
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
    // |soft_excess| what its memory_bound_soft states (the echo included); |memory_peak| is at
    // most the budget (what the budget admitted)
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
    // bound - max(allowance / 2, the delivery reserve), leaving the rest to deliver what was
    // walked.
    // Set once the request is parsed; before that the hard cap alone bounds the attempt.
    // |memory_budget|: the request's bounds.max_memory_mb in bytes, 0 for none (the soft
    // excess is observed under one, and what was held is bounded by it)
    void set_bound(size_t seeds, double time_budget_ms, uint64_t memory_budget = 0);
    double bound_ms() const { return bound_ms_; }
    double walk_until_ms() const { return walk_until_ms_; }
    bool enforced() const { return enforced_ && managed(); }
    // the time left before the seeds stop being walked (ms; infinity when the bound is not
    // enforced): what sizes the chunks of the walk's annotation reads (the handler's thread)
    double ms_left() const;

    // ---- the delivery reserve: the seeds stop being walked at
    // bound - max(allowance / 2, reserve), the reserve being kReserveMargin times the time to
    // build and compress what the response will hold — the text of the seeds finished so far
    // (exact, written as each was built) and of the seed being walked (its modelled account
    // less its coordinate share / the account per text byte of the requested detail, plus the
    // coordinate share / kCoordinateAccountPerTextByte) — at the stated rates, or the build
    // rate measured on this attempt's own seeds, plus the time from the walk-until to the
    // walk's end (delivery_stop_ms, or this server's measured stop latency). The handler's
    // thread: each moves the walk-until.
    // |detail|: the requested output.detail (its account per text byte: measured on this
    // server, else the configured one for JSON or for a graphlet)
    void set_delivery_detail(const std::string &detail);
    const std::string& delivery_detail() const { return detail_; }
    // the seed being walked holds |account| modelled bytes, |coordinates| of them its record
    // coordinates' share (ResourceAccount::coordinates; 0 without output.coordinates), at
    // every level's end
    void progress(uint64_t account, uint64_t coordinates = 0);
    // a seed's result was written as |text_bytes| of text, built in |build_seconds|, from a
    // walk whose modelled account ended at |account| bytes (0: no walk, a failed seed), of
    // which |coordinate_account| is the coordinate share and |coordinate_text| the exact text
    // it wrote (coordinates_text_bytes): the ratio of the rest is measured without both, so
    // that a request with coordinates measures what the same walk without them would (their
    // ratio would lower the estimate of every later attempt without them)
    void note_delivered(uint64_t text_bytes, double build_seconds, uint64_t account = 0,
                        uint64_t coordinate_account = 0, uint64_t coordinate_text = 0);
    // what this server measured before the attempt started: it replaces the configured rates
    // and ratios (a measurement of the attempt's own seeds replaces it when more conservative)
    void set_measured(const DeliveryMeasurements &measured);
    // the slowest build rate and the smallest account-to-text ratio measured on this
    // attempt's own seeds (0: none; the ratio without the coordinate share), and the longest
    // time a seed's walk ended after the walk-until (0: none did), for the server's
    // measurements. The ratio is pooled whether or not the output carried coordinates: their
    // account and text are left out of it exactly (note_delivered), so a request with them
    // adds what the same request without them would
    double own_build_mbps() const;
    double own_account_per_text_byte() const;
    double own_stop_ms() const;
    // the reserve now (ms)
    double reserve_ms() const;

    // ---- stops
    // asks the walk to stop; false when another reason was first (the first one wins)
    bool request_stop(graph::traversal::ExternalStop reason);
    // POST /traverse/cancel: records when it was asked, then request_stop(CANCELLED)
    bool cancel();
    // the handler's poll, at the walker's checkpoints and between seeds: the stop flag, and at
    // every poll_stride-th poll (or with |force|: between seeds, before a paced read's chunk,
    // in the lookahead) the walk-until (enforced only) and the client's connection (at most
    // every client_check_ms). Records the first poll that
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
    // The longest uninterruptible piece of the request's walks so far (ms; an annotation read
    // or chunk, or a head piece — the walk between two readings of the clock for a stop,
    // DecodePacer::max_uninterruptible_ms),
    // stated in usage as observed_max_uninterruptible_ms
    void note_max_uninterruptible_ms(double ms);
    double max_uninterruptible_ms() const;
    // The longest time between two delivery checks while a seed's text or the response was
    // written (json_text's max_gap_ms: e.g. the preparation of one large token), which the
    // server adds to deadline_check.observed_max_uninterruptible_ms (not to usage)
    void note_delivery_gap_ms(double ms);
    double max_delivery_gap_ms() const;

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
    friend class AttemptRegistry;
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
    double max_uninterruptible_ms_ = 0;               // note_max_uninterruptible_ms
    double max_delivery_gap_ms_ = 0;                  // note_delivery_gap_ms
    // the delivery reserve's state (the handler's thread; written under |mutex_|)
    std::string detail_;
    double configured_ratio_ = 20;                    // the configured account per text byte
    double server_ratio_ = 0;                         // measured on this server (0: none)
    double own_ratio_ = 0;                            // measured on this attempt (0: none)
    uint64_t delivered_bytes_ = 0;                    // the text of the seeds finished
    uint64_t walking_account_ = 0;                    // the seed being walked, its account
    uint64_t walking_coordinates_ = 0;                // ... and its coordinate share
    double measured_build_mbps_ = 0;                  // this attempt's own (0: none yet)
    DeliveryMeasurements server_;                     // set_measured
    // the lowest walk-until a poll checked (the handler's thread writes it; usage reads it) and
    // the walk-until in force when it stopped the walk (usage.bound.walk_until_ms states the
    // latter, else the former)
    std::atomic<double> checked_walk_until_ms_ { std::numeric_limits<double>::infinity() };
    std::optional<double> tripped_walk_until_ms_;
    double own_stop_ms_ = 0;                          // own_stop_ms()
    // the account per text byte in use, and the walked seed's estimated text
    double ratio_locked() const;
    // bound - max(allowance / 2, reserve), under |mutex_|
    void update_walk_until_locked();
    std::string reason_;                              // set by finish()
    int status_ = 0;
    std::optional<size_t> bytes_;
    Json::Value final_usage_;                         // the usage without per_seed, at finish
    bool tombstone_ = false;
    // The hold of a tombstone, or of a finished attempt's id (the registry's, under its
    // mutex): it is kept while EITHER clock says it is live — the steady clock until
    // |tomb_steady_until_| (exclusive: the duration computed when it was set or extended, plus
    // the one ms the wall clock's inclusive expiry adds; a forward step of the wall clock does
    // not shorten it), and the wall clock while it reads at most |tomb_wall_until_ms_|
    // (inclusive: the strict not_after check still admits a copy at not_after_ms itself, so a
    // hold to not_after_ms + skew must include that instant) — so a covered not_after_ms
    // has passed on that clock whenever the hold is gone. Never shortened
    Clock::time_point tomb_steady_until_;
    uint64_t tomb_wall_until_ms_ = 0;
    // the largest not_after_ms a cancel or a refused copy named for the id (none: none did)
    std::optional<uint64_t> tomb_not_after_ms_;
    // a finished attempt whose retention ended while its hold had not: kept among the
    // tombstones (in the registry's |tomb_expiry_|) until the hold ends
    bool held_ = false;
};

/**
 * The attempts of one server process: the running ones by attempt_id, and the finished ones
 * (and tombstones) for the retention period. An attempt_id is unique per server process while
 * it is running or retained: a second request with it is refused (409), never run twice. A
 * finished attempt sent with not_after_ms is retained at least until not_after_ms +
 * clock_skew_ms (within the cap from its finish), so that a replay of the request is refused
 * while it could still be admitted (a ledger releases on a finished state, and a replay
 * arriving after retention_s, or after retention_count later finishes, would otherwise run
 * again).
 */
class AttemptRegistry {
  public:
    explicit AttemptRegistry(AttemptSettings settings);

    const AttemptSettings& settings() const { return settings_; }
    // random per process (16 hex digits): a ledger tells a restarted backend (whose attempts
    // are gone) from one that forgot an attempt after the retention period
    const std::string& server_instance() const { return instance_; }

    // Why start() registered nothing (a 409): the request names another server_instance
    // (|instance_mismatch|, |body| the whole 409 body, instance_mismatch_json), the id is
    // running, retained or tombstoned (|body| the other attempt's state; a tombstone's with
    // its suppression, judged against the refused request's not_after_ms), or the request's
    // not_after_ms has passed (|expired|, |body| the whole 409 body, expired_json)
    struct StartRefusal {
        bool expired = false;
        Json::Value body;
        bool instance_mismatch = false;
    };
    // Registers a managed attempt, unless it names another server_instance (checked first:
    // the ids of another process mean nothing here), its id is running or retained (or
    // tombstoned) — checked before the not_after check, so that an expired refusal always
    // means that no attempt with the id exists on this server_instance — or its not_after_ms
    // has passed on the wall clock: then nothing is registered (an expired request is kept
    // nowhere: any later copy of it is expired too). A tombstoned id's refusal extends the
    // tombstone to the request's not_after_ms + clock_skew_ms (within the cap, never
    // shortened), so that every later copy of it is refused while it could still be admitted;
    // a running or finished attempt's refusal extends its hold likewise (a running one's from
    // its finish)
    std::optional<StartRefusal> start(const std::shared_ptr<Attempt> &attempt);
    // the wall clock not_after_ms is compared with (Unix epoch ms)
    uint64_t now_ms() const;
    // the `attempts` block of both capabilities routes: the fields, routes, retention, the
    // bound's parts with their number types (integers, ms), the clock skew a ledger adds to
    // not_after_ms, and the rules as references to the SPEC sections that state them
    Json::Value capabilities_json() const;
    // POST /traverse/cancel {attempt_id, wait_ms, not_after_ms}: (HTTP status, body). An
    // unknown id is tombstoned (404, tombstone: true) for at least retention_s and, when the
    // cancel names the request's not_after_ms, until not_after_ms + clock_skew_ms (at most
    // max(tombstone_max_s, retention_s) from now); every tombstone answer states
    // suppressed_until_ms and covers_admission. A repeated cancel never shortens the
    // tombstone, extends it to its own not_after_ms, and keeps it at least retention_s from
    // now. When retention_count tombstones are held already, or retention_s is 0, the id is
    // not tombstoned and the cancel is refused (429, tombstone: false, reason): nothing is
    // promised for it
    std::pair<int, Json::Value> cancel(const std::string &id, uint64_t wait_ms,
                                       std::optional<uint64_t> not_after_ms = std::nullopt);
    // GET /traverse/attempt/{id}: (HTTP status, body)
    std::pair<int, Json::Value> state(const std::string &id);
    // the handler returned: the attempt is finished, and kept for the retention period
    void finish(const std::shared_ptr<Attempt> &attempt, const std::string &reason, int status,
                std::optional<size_t> bytes);
    // the retention rule, stated by 404s and GET /traverse/capabilities
    std::string retention_text() const;
    // The longest single annotation read — one uninterruptible piece of a deadline — that a
    // /traverse of this process made (ms; with or without attempt_id): stated as
    // deadline_check.observed_max_uninterruptible_ms, which is an observation, not a bound
    void note_uninterruptible(double ms);
    uint64_t observed_max_uninterruptible_ms() const;
    // The delivery rates measured on this server's /traverse responses of at least
    // kMeasuredTextBytes (MB/s): a seed's result text built, a response compressed; and per
    // output.detail the ratio of a seed's modelled account to its text. The slowest rate and
    // the smallest ratio of the last kRateWindow of each replace the configured ones in every
    // attempt started afterwards (0: none measured yet)
    void note_build_rate(double mbps);
    void note_compress_rate(double mbps);
    void note_account_per_text_byte(const std::string &detail, double ratio);
    // the time an attempt's walk ended after its walk-until (ms; 0: it did not): the longest
    // of the last kRateWindow replaces delivery_stop_ms when longer
    void note_stop_latency(double ms);
    DeliveryMeasurements measured() const;

  private:
    void expire_locked();
    Attempt::Clock::time_point now() const;
    Json::Value unknown_json(const std::string &id) const;
    // extends |held|'s hold — a tombstone's, or a finished (or running) attempt's — to at
    // least retention_s from now (|at_least_retention|: tombstones) and, given |not_after_ms|,
    // to not_after_ms + clock_skew_ms within the cap from now: never shortened, re-keyed in
    // |tomb_expiry_| when it is there (under |mutex_|)
    void hold_locked(Attempt &held, std::optional<uint64_t> not_after_ms,
                     bool at_least_retention = true);
    // a tombstone's suppression, judged against |against| — and, with |fallback|, without it
    // against the largest not_after_ms named for the id (a cancel's and GET's answers; never a
    // refused request's, whose own not_after_ms alone says whether a copy of IT is covered):
    // suppressed_until_ms, not_after_ms, covers_admission and, when that is false,
    // covers_admission_reason (under |mutex_|)
    void add_suppression_locked(Json::Value *j, const Attempt &tomb,
                                std::optional<uint64_t> against, bool fallback) const;
    // whether |a|'s hold is live now (either clock; under |mutex_|)
    bool hold_live_locked(const Attempt &a, Attempt::Clock::time_point t, uint64_t wall) const;
    // the cap of a tombstone's hold (ms): max(tombstone_max_s, retention_s)
    uint64_t tombstone_cap_ms() const;

    AttemptSettings settings_;
    std::string instance_;
    mutable std::mutex mutex_;
    std::map<std::string, std::shared_ptr<Attempt>> attempts_;
    // finished attempts, oldest first: kept retention_s, at most retention_count of them
    std::deque<std::shared_ptr<Attempt>> retained_;
    // the tombstones, and the finished attempts held past their retention, by their
    // steady-clock expiry (and id), apart from the retained attempts: each is kept its whole
    // hold, whatever finishes after it (sharing retention_count, later finishes would evict a
    // tombstone early and the cancelled id would run); one whose steady expiry passed while its
    // wall-clock one has not (the wall clock stepped back)
    // is re-keyed to the time its wall clock still needs. A held finished attempt is never
    // dropped early, so the table can exceed retention_count by them (a cancel of an unknown id
    // is refused, 429, while it holds that many): what it holds is bounded by the attempts
    // that finish within the cap
    std::set<std::pair<Attempt::Clock::time_point, std::string>> tomb_expiry_;
    // microseconds, rounded up
    std::atomic<uint64_t> max_uninterruptible_us_ { 0 };
    mutable std::mutex rates_mutex_;
    std::deque<double> build_rates_;
    std::deque<double> compress_rates_;
    std::deque<double> stop_latencies_;
    std::map<std::string, std::deque<double>> ratios_;
};

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_TRAVERSE_ATTEMPTS_HPP__
