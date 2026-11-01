#include "gtest/gtest.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <json/json.h>

#include "cli/server_checks.hpp"
#include "cli/server_envelope.hpp"
#include "cli/server_fanout.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"


// The multi-graph fan-out of /resolve and /pattern (server_fanout.hpp; SPEC-multigraph-envelope.md
// sections 4 and 5): its deadlines, the work cutoff with a completed pair beside a blocked or a
// cut sibling, the lifecycle of workers that outlive the answer, and the registry by
// request_id. Each test waits until every worker of its fan-outs ended, so that a sanitizer
// build (CMAKE_BUILD_TYPE Address or Threads) checks the workers' whole lives

namespace {

using namespace mtg::cli;
using namespace std::chrono_literals;

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errors;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errors)) << errors;
    return v;
}

// A reference sample of SPEC-multigraph-envelope.md, {route, status, note, request?, body}
Json::Value sample(const std::string &name) {
    // run from the build directory, like every test reading ../tests/data
    const std::string path = "../tests/data/multigraph_contract/" + name;
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.good()) << path << " not found (run from the build directory)";
    const std::string text { std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>() };
    return parse(text);
}

std::set<std::string> members(const Json::Value &v) {
    const auto names = v.getMemberNames();
    return std::set<std::string>(names.begin(), names.end());
}

// every worker of the test's fan-outs ended (their shared state freed)
void wait_for_workers() {
    for (int i = 0; i < 2000 && fanouts_alive(); ++i) {
        std::this_thread::sleep_for(5ms);
    }
    ASSERT_EQ(0u, fanouts_alive());
}

// Runs each task on a thread of its own (the graphs' pool in the server); joined by the
// destructor, after the workers ended
class Threads {
  public:
    Launch launch() {
        return [this](std::function<void()> task) {
            std::lock_guard<std::mutex> lock(mutex_);
            threads_.emplace_back(std::move(task));
        };
    }
    ~Threads() {
        for (auto &t : threads_) {
            t.join();
        }
    }

  private:
    std::mutex mutex_;
    std::vector<std::thread> threads_;
};

// A gate a test job waits at (ignoring the request's stop: an uninterruptible piece of work, a
// load) until the test opens it
class Gate {
  public:
    void open() {
        std::lock_guard<std::mutex> lock(mutex_);
        open_ = true;
        cv_.notify_all();
    }
    void wait() {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [&]() { return open_; });
    }

  private:
    std::mutex mutex_;
    std::condition_variable cv_;
    bool open_ = false;
};

const std::vector<SelectedPair> kPairs {
    { "a", { "a.dbg", "a.anno" } },
    { "b", { "b.dbg", "b.anno" } },
    { "c", { "c.dbg", "c.anno" } },
};

std::vector<Json::Value> refs_of(const std::vector<SelectedPair> &pairs) {
    std::vector<Json::Value> refs;
    for (const auto &pair : pairs) {
        refs.push_back(envelope_pair_ref(pair));
    }
    return refs;
}

// The deadlines of a request arriving now with |budget_ms| and |handles| (the skew allowance
// 2000 ms, the reserve 50 ms, the transport's cap 899 s)
RequestDeadlines deadlines_now(double budget_ms, const RequestHandles &handles = {}) {
    const auto wall = std::chrono::system_clock::now();
    return request_deadlines(wall, wall, std::chrono::steady_clock::now(), budget_ms, handles,
                             2000, 50, 899'000);
}

uint64_t now_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count();
}

FanoutSettings settings_of(const std::vector<SelectedPair> &pairs,
                           const RequestDeadlines &deadlines,
                           const std::string &route = "/resolve") {
    FanoutSettings s;
    s.route = route;
    s.pairs = pairs;
    s.index_fps.assign(pairs.size(), "");
    s.index_fps[0] = "fp-a";
    s.deadlines = deadlines;
    s.poll_ms = 10;
    return s;
}

std::shared_ptr<RequestEntry> entry_of(const std::vector<SelectedPair> &pairs,
                                       std::optional<std::string> id = std::nullopt) {
    return std::make_shared<RequestEntry>(std::move(id), "/resolve", "inst",
                                          std::chrono::system_clock::now(), refs_of(pairs));
}

// a pair's answer: its single-graph answer with the pair's tags
PairResult answered(const PairControl &control, uint64_t memory_bytes = 7) {
    PairResult r;
    Json::Value answer;
    answer["work"]["memory_bytes"] = Json::UInt64(memory_bytes);
    answer["candidates"] = Json::Value(Json::arrayValue);
    r.entry = envelope_entry_answered(control.pair(), "", answer);
    r.answered = true;
    r.memory_bytes = memory_bytes;
    return r;
}

const std::set<std::string> kNotReachedMembers { "code", "reason", "error", "http_status" };
const std::set<std::string> kUnfinishedMembers { "code", "reason", "error", "usage",
                                                 "http_status" };

} // namespace


// ---------------------------------------------------------------- the request's handles

TEST(ServerFanout, RequestHandles) {
    const RequestHandles none = read_request_handles(parse(R"({"sequence": "ACGT"})"));
    EXPECT_FALSE(none.request_id);
    EXPECT_FALSE(none.not_after_ms);
    EXPECT_FALSE(none.finish_by_ms);
    const RequestHandles all = read_request_handles(parse(
            R"({"request_id": "svc-res-1:a_b.c", "not_after_ms": 5, "finish_by_ms": 9})"));
    EXPECT_EQ("svc-res-1:a_b.c", *all.request_id);
    EXPECT_EQ(5u, *all.not_after_ms);
    EXPECT_EQ(9u, *all.finish_by_ms);
    // attempt_id's token rule, not_after_ms's integer rule for both instants
    for (const char *bad : { R"({"request_id": ""})", R"({"request_id": "a b"})",
                             R"({"request_id": 7})", R"({"finish_by_ms": -1})",
                             R"({"finish_by_ms": 1.5})", R"({"finish_by_ms": "9"})",
                             R"({"not_after_ms": 1e20})" }) {
        EXPECT_THROW(read_request_handles(parse(bad)), InvalidRequest) << bad;
    }
    // the strict parser of /resolve accepts them, and refuses what the reader refuses
    EXPECT_NO_THROW(parse_resolve_request(parse(
            R"({"sequence": "ACGT", "labels": ["x"], "request_id": "r1", "not_after_ms": 1,
                "finish_by_ms": 2})")));
    EXPECT_THROW(parse_resolve_request(parse(
            R"({"sequence": "ACGT", "labels": ["x"], "request_id": "r 1"})")), InvalidRequest);
    // /traverse reads finish_by_ms beside its attempt fields
    EXPECT_EQ(42u, *read_finish_by_ms(parse(R"({"finish_by_ms": 42})"), "request"));
    EXPECT_THROW(read_finish_by_ms(parse(R"({"finish_by_ms": -42})"), "request"),
                 InvalidRequest);
}

// ---------------------------------------------------------------- the deadlines

TEST(ServerFanout, EffectiveDeadlineAndCutoff) {
    // (the samples' instants are illustrative: their ms and ISO forms are not one instant)
    const auto wall = std::chrono::system_clock::time_point(std::chrono::milliseconds(
            int64_t(1793520000012)));
    const auto now = std::chrono::steady_clock::now();
    // arrival 08:00:00.012, budget 60 s, finish_by_ms 08:01:00.000 less the 2 s allowance:
    // the effective deadline is finish_by's, the cutoff 250 ms before it
    RequestHandles h;
    h.not_after_ms = 1793520010000;
    h.finish_by_ms = 1793520060000;
    RequestDeadlines d = request_deadlines(wall, wall, now, 60'000, h, 2000, 250, 899'000);
    Json::Value j = d.json();
    EXPECT_EQ("2026-11-01T08:00:58.000Z", j["effective_deadline_at"].asString());
    EXPECT_EQ("2026-11-01T08:00:57.750Z", j["cutoff_at"].asString());
    EXPECT_EQ(1793520010000u, j["not_after_ms"].asUInt64());
    EXPECT_EQ(1793520060000u, j["finish_by_ms"].asUInt64());
    EXPECT_FALSE(d.transport_bound);
    ASSERT_TRUE(d.finish);
    EXPECT_EQ(*d.finish, d.effective);
    // a short route budget with a later finish_by_ms: the effective deadline is the arrival
    // plus the budget, never finish_by's alone
    d = request_deadlines(wall, wall, now, 2000, h, 2000, 250, 899'000);
    EXPECT_EQ("2026-11-01T08:00:02.012Z", d.json()["effective_deadline_at"].asString());
    EXPECT_EQ("2026-11-01T08:00:01.762Z", d.json()["cutoff_at"].asString());
    EXPECT_LT(d.effective, *d.finish);
    // nothing given: the instants are not stated (null), though the cutoff applies; a
    // budget at the transport's cap is bound by the transport
    d = request_deadlines(wall, wall, now, 899'000, {}, 2000, 250, 899'000);
    EXPECT_EQ(deadline_json(std::nullopt, std::nullopt, std::nullopt, std::nullopt), d.json());
    EXPECT_TRUE(d.transport_bound);
    EXPECT_EQ(now + 899'000ms - 250ms, d.cutoff);
    // the arrival read later by the handler: its steady instant lies before now
    d = request_deadlines(wall - 1500ms, wall, now, 10'000, {}, 0, 0, 899'000);
    EXPECT_EQ(now - 1500ms, d.arrival);
    EXPECT_EQ(now + 8500ms, d.effective);
    // a finish_by_ms below the skew allowance: already passed
    h.finish_by_ms = 1;
    d = request_deadlines(wall, wall, now, 60'000, h, 2000, 250, 899'000);
    EXPECT_LT(d.effective, now);
}

// ---------------------------------------------------------------- the fan-out

// fanout F4 of the release-10 review: a finish_by_ms far in the future (up to 2^53 - 1, as
// read_finish_by_ms accepts) binds nothing and overflows no clock: the effective deadline is
// the budget's (or the transport's), never an instant in the past
TEST(ServerFanout, FarFinishByMsBindsNothing) {
    for (uint64_t finish_by : { uint64_t(10'000'000'000'000), uint64_t(12'000'000'000'000),
                                (uint64_t(1) << 53) - 1 }) {
        RequestHandles h;
        h.finish_by_ms = finish_by;
        const RequestDeadlines d = deadlines_now(60'000, h);
        const RequestDeadlines without = deadlines_now(60'000);
        const double effective_ms = std::chrono::duration<double, std::milli>(
                d.effective - d.arrival).count();
        EXPECT_NEAR(60000.0, effective_ms, 1.0) << finish_by;
        EXPECT_FALSE(d.transport_bound) << finish_by;
        EXPECT_GT(d.cutoff, d.arrival) << finish_by;
        ASSERT_TRUE(d.finish) << finish_by;
        EXPECT_GT(*d.finish, d.effective) << finish_by;
        // the instants written: the budget's, after the arrival, in this millennium
        EXPECT_EQ(iso_utc(without.effective_wall).substr(0, 4),
                  iso_utc(d.effective_wall).substr(0, 4)) << finish_by;
        EXPECT_EQ(finish_by, d.json()["finish_by_ms"].asUInt64()) << finish_by;

        // beyond the transport's cap the transport binds
        const RequestDeadlines t = deadlines_now(2e6, h);
        EXPECT_TRUE(t.transport_bound) << finish_by;
        const double transport_ms = std::chrono::duration<double, std::milli>(
                t.effective - t.arrival).count();
        EXPECT_NEAR(899000.0, transport_ms, 1.0) << finish_by;
    }
}

// One pair runs on the calling thread, as a single-graph request (no pool, no cutoff wait)
TEST(ServerFanout, OnePairRunsInline) {
    const std::thread::id handler = std::this_thread::get_id();
    const std::vector<SelectedPair> pairs { kPairs[0] };
    auto entry = entry_of(pairs);
    bool launched = false;
    FanoutOutcome out = run_fanout(settings_of(pairs, deadlines_now(60'000)), entry,
        [&](PairControl &control) {
            EXPECT_EQ(handler, std::this_thread::get_id());
            EXPECT_TRUE(control.inline_run());
            // inline: only finish_by_ms caps the pair's own deadline, never the cutoff (the
            // pair is the request, answered as a single-graph one)
            EXPECT_FALSE(control.deadline_cap());
            EXPECT_FALSE(control.cutoff_cap());
            EXPECT_TRUE(control.start_running());
            return answered(control);
        },
        [&](std::function<void()>) { launched = true; });
    EXPECT_FALSE(launched);
    EXPECT_FALSE(out.client_gone);
    ASSERT_EQ(1u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    // the job tags its own entry (here without a manifest)
    EXPECT_TRUE(out.entries[0]["index_fp"].isNull());
    entry->response(200, 10);
    EXPECT_TRUE(entry->ended());
    const Json::Value j = entry->json(3600);
    EXPECT_EQ("ended", j["state"].asString());
    EXPECT_EQ("answered", j["reason"].asString());
    EXPECT_EQ("answered", j["pairs"][0]["outcome"].asString());
    EXPECT_EQ(7u, j["pairs"][0]["usage"]["memory_bytes"].asUInt64());
    wait_for_workers();
}

// fanout F6 of the release-10 review: a pooled pair whose memory wait gave up is not reached
// for the same reason whichever thread decides it, the worker (here, before the handler's
// cutoff) or the handler at the cutoff: on a transport-bound request "transport"
TEST(ServerFanout, NotReachedReasonIsTheRequestsWhoeverDecides) {
    for (const double budget : { 2e6, 60'000.0 }) {
        const RequestDeadlines d = deadlines_now(budget);
        const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
        auto entry = entry_of(pairs);
        FanoutOutcome out;
        {
            Threads threads;
            out = run_fanout(settings_of(pairs, d), entry, [&](PairControl &control) {
                if (control.index() == 0) {
                    EXPECT_TRUE(control.start_running());
                    return answered(control);
                }
                // the wait gave up (as at its limit): not started
                control.waiting_memory();
                PairResult r;
                r.not_started = true;
                return r;
            }, threads.launch());
        }
        ASSERT_EQ(2u, out.entries.size());
        const Json::Value &refusal = out.entries[1]["refusal"];
        EXPECT_EQ("not_reached", refusal["code"].asString()) << budget;
        EXPECT_EQ(d.transport_bound ? "transport" : "cutoff", refusal["reason"].asString())
            << budget;
        EXPECT_EQ(budget > 899'000, d.transport_bound);
    }
    wait_for_workers();
}

// fanout F2 of the release-10 review: a pair does not start loading its index once its wait
// limit passed (inline: finish_by_ms less the skew allowance) or its client left, also when
// the memory is free: not reached (its client gone: nothing answered), never loaded
TEST(ServerFanout, NoLoadPastTheWaitLimit) {
    const std::vector<SelectedPair> pairs { kPairs[0] };
    RequestHandles h;
    h.finish_by_ms = now_ms() - 60'000;
    auto entry = entry_of(pairs);
    bool loaded = false;
    FanoutOutcome out = run_fanout(settings_of(pairs, deadlines_now(60'000, h)), entry,
        [&](PairControl &control) {
            PairResult r;
            control.waiting_memory();
            // the memory is free (the reservation is granted at once): the predicate is the
            // only check before the load
            EXPECT_TRUE(control.wait_predicate()());
            loaded = control.start_loading(100);
            r.not_started = !loaded;
            return r;
        },
        [&](std::function<void()>) { ADD_FAILURE() << "one pair runs inline"; });
    EXPECT_FALSE(loaded);
    EXPECT_FALSE(out.client_gone);
    ASSERT_EQ(1u, out.entries.size());
    EXPECT_EQ("refused", out.entries[0]["outcome"].asString());
    const Json::Value &refusal = out.entries[0]["refusal"];
    EXPECT_EQ(kNotReachedMembers, members(refusal));
    EXPECT_EQ(503, refusal["http_status"].asInt());
    EXPECT_EQ("not_reached", refusal["code"].asString());
    EXPECT_EQ("cutoff", refusal["reason"].asString());
    EXPECT_EQ("never_started", entry->json(3600)["pairs"][0]["state"].asString());

    // the client gone before the load: nothing loaded, nothing answered
    auto gone_entry = entry_of(pairs);
    FanoutSettings gone_settings = settings_of(pairs, deadlines_now(60'000));
    gone_settings.gone = []() { return true; };
    loaded = false;
    out = run_fanout(gone_settings, gone_entry, [&](PairControl &control) {
        PairResult r;
        control.waiting_memory();
        loaded = control.start_loading(100);
        r.not_started = !loaded;
        return r;
    }, [&](std::function<void()>) { ADD_FAILURE() << "one pair runs inline"; });
    EXPECT_FALSE(loaded);
    EXPECT_TRUE(out.client_gone);

    // before the limit, with the client there: the load starts
    auto later = entry_of(pairs);
    loaded = false;
    out = run_fanout(settings_of(pairs, deadlines_now(60'000)), later, [&](PairControl &control) {
        control.waiting_memory();
        loaded = control.start_loading(100);
        EXPECT_TRUE(control.start_running());
        return answered(control);
    }, [&](std::function<void()>) { ADD_FAILURE() << "one pair runs inline"; });
    EXPECT_TRUE(loaded);
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    wait_for_workers();
}

// E4's test: one completed pair beside a sibling blocked waiting for the memory to load its
// index, under a short route budget and a later finish_by_ms: at the cutoff (from the
// effective deadline, not finish_by's) the completed answer survives and the sibling is
// not_reached; the blocked worker gives up its wait and ends
TEST(ServerFanout, CompletedPairSurvivesABlockedSibling) {
    RequestHandles h;
    h.finish_by_ms = now_ms() + 60'000;
    const RequestDeadlines d = deadlines_now(300, h);
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto entry = entry_of(pairs, "req-1");
    // the memory is held elsewhere: the sibling's wait cannot end before the cutoff
    LoadReservations reservations(1);
    ASSERT_TRUE(reservations.reserve(1));
    const auto started = std::chrono::steady_clock::now();
    FanoutOutcome out;
    {
        Threads threads;
        out = run_fanout(settings_of(pairs, d), entry, [&](PairControl &control) {
            PairResult r;
            if (control.index() == 0) {
                EXPECT_TRUE(control.start_running());
                return answered(control);
            }
            control.waiting_memory();
            if (!reservations.reserve(1, control.wait_predicate(), 5)) {
                r.not_started = true;
                return r;
            }
            ADD_FAILURE() << "the memory was never free";
            reservations.release(1);
            return r;
        }, threads.launch());
        const double waited = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - started).count();
        // at the cutoff: 300 ms less the 50 ms reserve, long before finish_by_ms
        EXPECT_GE(waited, 240);
        EXPECT_LT(waited, 5000);
    }
    ASSERT_EQ(2u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    const Json::Value &refusal = out.entries[1]["refusal"];
    EXPECT_EQ("refused", out.entries[1]["outcome"].asString());
    EXPECT_EQ(kNotReachedMembers, members(refusal));
    EXPECT_EQ(503, refusal["http_status"].asInt());
    EXPECT_EQ("not_reached", refusal["code"].asString());
    EXPECT_EQ("cutoff", refusal["reason"].asString());
    EXPECT_NE(std::string::npos, refusal["error"].asString().find("memory")) << refusal;
    // the envelope lists it (sample 03's not_reached shape)
    EnvelopeContext context;
    context.request_id = "req-1";
    context.server_instance = "inst";
    context.deadline = d.json();
    const Json::Value envelope = post_envelope("/resolve", { "a", "b" }, out.entries, context, 1);
    EXPECT_EQ(1u, envelope["answered"].asUInt());
    EXPECT_EQ(1u, envelope["refused"].asUInt());
    ASSERT_EQ(1u, envelope["not_reached"].size());
    EXPECT_EQ(envelope_pair_ref(kPairs[1]), envelope["not_reached"][0]);
    EXPECT_EQ(7u, envelope["work"]["memory_bytes"].asUInt64());
    const Json::Value sample_refusal
            = sample("03_resolve_envelope_cutoff.json")["body"]["answers"][1]["refusal"];
    EXPECT_EQ(members(sample_refusal), members(refusal));
    wait_for_workers();
    entry->response(200, 99);
    const Json::Value j = entry->json(3600);
    EXPECT_EQ("ended", j["state"].asString());
    EXPECT_EQ("never_started", j["pairs"][1]["state"].asString());
    EXPECT_EQ("not_reached", j["pairs"][1]["outcome"].asString());
    EXPECT_TRUE(j["pairs"][1]["usage"].isNull());
    reservations.release(1);
}

// A pair started and not finished at the cutoff (an uninterruptible piece of work): unfinished
// with its usage so far and told to stop; its worker ends after the answer, and the registry
// states the end with its final usage — the entry stays running until then
TEST(ServerFanout, StartedPairCutAtTheCutoffEndsLater) {
    const RequestDeadlines d = deadlines_now(300);
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1], kPairs[2] };
    RequestSettings rs;
    rs.server_instance = "inst";
    RequestRegistry registry(rs);
    auto entry = std::make_shared<RequestEntry>(std::string("req-cut"), "/resolve", "inst",
                                                std::chrono::system_clock::now(),
                                                refs_of(pairs));
    ASSERT_FALSE(registry.start(entry));
    Gate gate;
    std::atomic<bool> told_to_stop { false };
    Threads threads;
    FanoutOutcome out = run_fanout(settings_of(pairs, d), entry, [&](PairControl &control) {
        PairResult r;
        if (control.index() == 0) {
            EXPECT_TRUE(control.start_running());
            return answered(control, 11);
        }
        if (control.index() == 1) {
            // a pair that loads (unfinished, reason loading)
            EXPECT_TRUE(control.start_loading(4096));
            gate.wait();
            control.start_running();
            r.aborted = control.stopped();
            return r;
        }
        EXPECT_TRUE(control.start_running());
        // a pooled pair's own budget is capped at the cutoff (the owner's decision of
        // 2026-11-01); finish_by_ms caps nothing here (none given)
        EXPECT_FALSE(control.deadline_cap());
        EXPECT_EQ(d.cutoff, control.cutoff_cap());
        // some CPU time to state, then uninterruptible until the gate opens
        const auto until = std::chrono::steady_clock::now() + 30ms;
        volatile uint64_t x = 0;
        while (std::chrono::steady_clock::now() < until) {
            x = x + 1;
            control.stopped();
        }
        gate.wait();
        told_to_stop = control.stopped();
        r.aborted = true;
        r.memory_bytes = 123;
        return r;
    }, threads.launch());
    ASSERT_EQ(3u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    const Json::Value &loading = out.entries[1]["refusal"];
    EXPECT_EQ(kUnfinishedMembers, members(loading));
    EXPECT_EQ("unfinished", loading["code"].asString());
    EXPECT_EQ("loading", loading["reason"].asString());
    EXPECT_EQ("loading", loading["usage"]["state"].asString());
    EXPECT_EQ(4096u, loading["usage"]["memory_bytes"].asUInt64());
    const Json::Value &cut = out.entries[2]["refusal"];
    EXPECT_EQ(kUnfinishedMembers, members(cut));
    EXPECT_EQ(503, cut["http_status"].asInt());
    EXPECT_EQ("cutoff", cut["reason"].asString());
    const Json::Value sample_usage = sample("03_resolve_envelope_cutoff.json")
            ["body"]["answers"][2]["refusal"]["usage"];
    EXPECT_EQ(members(sample_usage), members(cut["usage"]));
    EXPECT_EQ("stopping", cut["usage"]["state"].asString());
    EXPECT_GE(cut["usage"]["elapsed_ms"].asDouble(), 200);
    EXPECT_TRUE(cut["usage"]["cpu_ms"].isNumeric());
    EXPECT_NE(std::string::npos, cut["error"].asString().find("GET /requests/req-cut"));
    // the answer is written while two workers still run: the request is running
    entry->response(200, 4821);
    Json::Value j = registry.state("req-cut").second;
    EXPECT_EQ("running", j["state"].asString());
    EXPECT_TRUE(j["finished_at"].isNull());
    EXPECT_TRUE(j["reason"].isNull());
    EXPECT_EQ(4821u, j["response"]["bytes"].asUInt64());
    EXPECT_TRUE(j["response"]["written"].asBool());
    EXPECT_EQ("running", j["pairs"][2]["state"].asString());
    EXPECT_EQ("unfinished", j["pairs"][2]["outcome"].asString());
    EXPECT_TRUE(j["pairs"][2]["usage"]["finished_at"].isNull());
    EXPECT_EQ(members(sample("17_requests_running.json")["body"]), members(j));
    // the workers end: the request ended, with each pair's final usage
    gate.open();
    wait_for_workers();
    EXPECT_TRUE(told_to_stop);
    j = registry.state("req-cut").second;
    EXPECT_EQ("ended", j["state"].asString());
    EXPECT_EQ("answered", j["reason"].asString());
    EXPECT_FALSE(j["retained_until"].isNull());
    EXPECT_EQ("ended", j["pairs"][2]["state"].asString());
    EXPECT_EQ("unfinished", j["pairs"][2]["outcome"].asString());
    EXPECT_EQ(123u, j["pairs"][2]["usage"]["memory_bytes"].asUInt64());
    EXPECT_FALSE(j["pairs"][2]["usage"]["finished_at"].isNull());
    EXPECT_EQ(11u, j["pairs"][0]["usage"]["memory_bytes"].asUInt64());
    const Json::Value ended = sample("18_requests_ended.json")["body"];
    EXPECT_EQ(members(ended), members(j));
    EXPECT_EQ(members(ended["pairs"][2]), members(j["pairs"][2]));
    EXPECT_EQ(members(ended["pairs"][2]["usage"]), members(j["pairs"][2]["usage"]));
}

// A pair a worker takes after not_after_ms is not_reached (reason not_after_ms), never started;
// pairs still queued then are decided at once rather than waited for
TEST(ServerFanout, NotAfterMsLeavesQueuedPairsUnreached) {
    RequestHandles h;
    h.not_after_ms = now_ms() + 150;
    const RequestDeadlines d = deadlines_now(60'000, h);
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto entry = entry_of(pairs);
    // the pool runs one task now, the other once the test lets it (a pool thread busy with
    // another request)
    std::vector<std::function<void()>> held;
    std::mutex held_mutex;
    std::vector<std::thread> threads;
    int launched = 0;
    FanoutOutcome out = run_fanout(settings_of(pairs, d, "/pattern"), entry,
        [&](PairControl &control) {
            EXPECT_TRUE(control.start_running());
            return answered(control);
        },
        [&](std::function<void()> task) {
            if (launched++ == 0) {
                threads.emplace_back(std::move(task));
            } else {
                std::lock_guard<std::mutex> lock(held_mutex);
                held.push_back(std::move(task));
            }
        });
    // decided once not_after_ms passed, not at the cutoff 60 s later
    ASSERT_EQ(2u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    const Json::Value &refusal = out.entries[1]["refusal"];
    EXPECT_EQ("not_reached", refusal["code"].asString());
    EXPECT_EQ("not_after_ms", refusal["reason"].asString());
    EXPECT_EQ(0u, refusal["error"].asString().find("pattern: pair (b, b.dbg, b.anno) was not "
                                                   "started: not_after_ms"));
    // the held task runs later: it starts nothing
    for (auto &task : held) {
        threads.emplace_back(std::move(task));
    }
    for (auto &t : threads) {
        t.join();
    }
    wait_for_workers();
    entry->response(200, 1);
    EXPECT_EQ("never_started", entry->json(0)["pairs"][1]["state"].asString());
    EXPECT_TRUE(entry->ended());
}

// The client leaves while the pairs run: nothing is answered, the workers are told to stop,
// and the registry records the end (client_gone) once they ended
TEST(ServerFanout, ClientGoneStopsTheWorkers) {
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    RequestSettings rs;
    rs.server_instance = "inst";
    RequestRegistry registry(rs);
    auto entry = std::make_shared<RequestEntry>(std::string("req-gone"), "/pattern", "inst",
                                                std::chrono::system_clock::now(),
                                                refs_of(pairs));
    ASSERT_FALSE(registry.start(entry));
    std::atomic<bool> gone { false };
    FanoutSettings s = settings_of(pairs, deadlines_now(60'000));
    s.gone = [&]() { return gone.load(); };
    std::atomic<int> running { 0 };
    Threads threads;
    std::thread leaver([&]() {
        while (running < 2) {
            std::this_thread::sleep_for(1ms);
        }
        gone = true;
    });
    FanoutOutcome out = run_fanout(s, entry, [&](PairControl &control) {
        EXPECT_TRUE(control.start_running());
        running++;
        PairResult r;
        while (!control.stopped()) {
            std::this_thread::sleep_for(1ms);
        }
        r.aborted = true;
        return r;
    }, threads.launch());
    leaver.join();
    EXPECT_TRUE(out.client_gone);
    entry->response(0, std::nullopt);
    wait_for_workers();
    const Json::Value j = registry.state("req-gone").second;
    EXPECT_EQ("ended", j["state"].asString());
    EXPECT_EQ("client_gone", j["reason"].asString());
    EXPECT_FALSE(j["response"]["written"].asBool());
    EXPECT_EQ("ended", j["pairs"][0]["state"].asString());
    EXPECT_EQ("unfinished", j["pairs"][0]["outcome"].asString());
    EXPECT_FALSE(j["pairs"][0]["usage"].isNull());
}

// ---------------------------------------------------------------- the registry

// A repeated request_id that is running or retained starts no work: 409, its entry with an
// error; no second worker. The lost response: written (status, bytes) while the client
// dropped it, the request ended once its worker ended
TEST(ServerFanout, DuplicateRequestIdStartsNoWork) {
    RequestSettings rs;
    rs.server_instance = "inst";
    rs.retention_s = 3600;
    rs.retention_count = 10;
    RequestRegistry registry(rs);
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto first = entry_of(pairs, "svc-res-1");
    ASSERT_FALSE(registry.start(first));
    Gate gate;
    std::atomic<int> jobs { 0 };
    Threads threads;
    PairJob job = [&](PairControl &control) {
        jobs++;
        EXPECT_TRUE(control.start_running());
        if (control.index() == 1)
            gate.wait();
        return answered(control);
    };
    // the first request is cut at its cutoff while pair b still runs, and answered
    FanoutOutcome out = run_fanout(settings_of(pairs, deadlines_now(200)), first, job,
                                   threads.launch());
    EXPECT_EQ("unfinished", out.entries[1]["refusal"]["code"].asString());
    first->response(200, 321);
    // a repeat while it runs: refused with the entry, and nothing of it runs
    auto second = entry_of(pairs, "svc-res-1");
    bool admitted = false;
    const auto refused = registry.start(second, [&]() -> std::optional<Json::Value> {
        admitted = true;
        return std::nullopt;
    });
    ASSERT_TRUE(refused);
    EXPECT_FALSE(admitted);
    EXPECT_EQ(2, jobs.load());
    const Json::Value sample19 = sample("19_requests_409_duplicate.json")["body"];
    EXPECT_EQ(members(sample19), members(*refused));
    EXPECT_EQ("running", (*refused)["state"].asString());
    EXPECT_EQ(0u, (*refused)["error"].asString().find("request_id 'svc-res-1' is running, or "
                                                     "retained after it finished"));
    EXPECT_NE(std::string::npos,
              (*refused)["error"].asString().find("requests are kept 3600 s after they "
                                                  "finish, at most 10"));
    gate.open();
    wait_for_workers();
    // ended and retained: still refused, the response that was written stated
    const auto again = registry.start(entry_of(pairs, "svc-res-1"));
    ASSERT_TRUE(again);
    EXPECT_EQ("ended", (*again)["state"].asString());
    EXPECT_EQ(200, (*again)["response"]["status"].asInt());
    EXPECT_EQ(321u, (*again)["response"]["bytes"].asUInt64());
    EXPECT_TRUE((*again)["response"]["written"].asBool());
    // an admission refusal (an expired not_after_ms) registers nothing
    const auto expired = registry.start(entry_of(pairs, "svc-res-2"),
                                        []() -> std::optional<Json::Value> {
        Json::Value body;
        body["state"] = "expired";
        return body;
    });
    ASSERT_TRUE(expired);
    EXPECT_EQ("expired", (*expired)["state"].asString());
    EXPECT_EQ(404, registry.state("svc-res-2").first);
}

// Retention by period and count, then 404 {error, request_id, server_instance}
TEST(ServerFanout, RetentionThen404) {
    auto t = std::make_shared<std::atomic<int64_t>>(0);
    RequestSettings rs;
    rs.server_instance = "inst";
    rs.retention_s = 60;
    rs.retention_count = 2;
    rs.clocks.steady = [t]() {
        return std::chrono::steady_clock::time_point(std::chrono::milliseconds(t->load()));
    };
    RequestRegistry registry(rs);
    auto end = [&](const std::string &id) {
        auto e = std::make_shared<RequestEntry>(id, "/pattern", "inst",
                                                std::chrono::system_clock::now(),
                                                refs_of({ kPairs[0] }), rs.clocks);
        EXPECT_FALSE(registry.start(e));
        e->pair_running(0);
        e->pair_outcome(0, "answered");
        e->pair_ended(0, 1, 0.5);
        e->response(200, 10);
        EXPECT_TRUE(e->ended());
    };
    end("r1");
    *t += 1000;
    end("r2");
    EXPECT_EQ(200, registry.state("r1").first);
    // by count: a third ended request drops the oldest
    *t += 1000;
    end("r3");
    auto [status, body] = registry.state("r1");
    EXPECT_EQ(404, status);
    const Json::Value sample20 = sample("20_requests_404.json")["body"];
    EXPECT_EQ(members(sample20), members(body));
    EXPECT_EQ("r1", body["request_id"].asString());
    EXPECT_EQ("inst", body["server_instance"].asString());
    EXPECT_EQ(0u, body["error"].asString().find("request_id 'r1' is not known on this server"));
    EXPECT_EQ(200, registry.state("r2").first);
    // by period: 60 s after it ended
    *t += 58'999;
    EXPECT_EQ(200, registry.state("r2").first);
    *t += 1;
    EXPECT_EQ(404, registry.state("r2").first);
    EXPECT_EQ(200, registry.state("r3").first);
    // a running request is never dropped, and its id runs again once dropped
    auto running = std::make_shared<RequestEntry>(std::string("r4"), "/pattern", "inst",
                                                  std::chrono::system_clock::now(),
                                                  refs_of({ kPairs[0] }), rs.clocks);
    EXPECT_FALSE(registry.start(running));
    *t += 3'600'000;
    EXPECT_EQ(200, registry.state("r4").first);
    EXPECT_EQ(404, registry.state("r3").first);
    EXPECT_FALSE(registry.start(std::make_shared<RequestEntry>(
            std::string("r3"), "/pattern", "inst", std::chrono::system_clock::now(),
            refs_of({ kPairs[0] }), rs.clocks)));
    // retention 0: nothing kept once ended
    RequestSettings none = rs;
    none.retention_s = 0;
    RequestRegistry forgetful(none);
    auto e = std::make_shared<RequestEntry>(std::string("x"), "/resolve", "inst",
                                            std::chrono::system_clock::now(),
                                            refs_of({ kPairs[0] }), rs.clocks);
    EXPECT_FALSE(forgetful.start(e));
    e->pair_never_started(0);
    e->response(400, 20);
    EXPECT_EQ(404, forgetful.state("x").first);
    EXPECT_NE(std::string::npos, forgetful.retention_text().find("not kept"));
}

// GET /requests/{attempt_id} of a /traverse attempt: the request object built from the
// attempt's state (sample 21)
TEST(ServerFanout, TraverseRequestFromTheAttempt) {
    const Json::Value s21 = sample("21_requests_traverse.json")["body"];
    const Json::Value pair = parse(R"({"graph": "bact_other",
        "graph_path": "/scratch/nvme8/refseq_r237/bact_other/graph.dbg",
        "annotation_path": "/scratch/nvme8/refseq_r237/bact_other/annotation.row_diff_brwt_coord.annodbg"})");
    const Json::Value j = traverse_request_json(s21["attempt"], pair, 3600);
    EXPECT_EQ(members(s21), members(j));
    EXPECT_EQ("svc-trv-91aa07-1", j["request_id"].asString());
    EXPECT_EQ("/traverse", j["route"].asString());
    EXPECT_EQ("ended", j["state"].asString());
    EXPECT_EQ("completed", j["reason"].asString());
    EXPECT_EQ("2026-10-21T10:14:29.385Z", j["retained_until"].asString());
    EXPECT_EQ(s21["pairs"], j["pairs"]);
    EXPECT_EQ(s21["response"], j["response"]);
    // a running attempt on a single-graph server: null tags, no response yet
    Json::Value running = s21["attempt"];
    running["state"] = "running";
    const Json::Value r = traverse_request_json(running, Json::Value(), 3600);
    EXPECT_EQ("running", r["state"].asString());
    EXPECT_TRUE(r["response"].isNull());
    EXPECT_TRUE(r["retained_until"].isNull());
    EXPECT_TRUE(r["pairs"][0]["graph"].isNull());
    EXPECT_TRUE(r["pairs"][0]["outcome"].isNull());
}

// finish_by_ms bounds a /traverse attempt's bound (and walk-until), enforced without
// attempt_id too
TEST(ServerFanout, FinishByCapsTheAttemptBound) {
    AttemptSettings settings;
    settings.clock_skew_ms = 2000;
    Attempt attempt(1, std::chrono::system_clock::now(), settings, "inst", true);
    attempt.set_ids(AttemptIds());
    attempt.set_finish_by(now_ms() + 2000 + 5000);
    EXPECT_FALSE(attempt.managed());
    EXPECT_TRUE(attempt.enforced());
    EXPECT_LE(attempt.bound_ms(), 5001);
    EXPECT_GE(attempt.bound_ms(), 4000);
    attempt.set_bound(1, 60'000);
    EXPECT_LE(attempt.bound_ms(), 5001);
    EXPECT_EQ("finish_by_ms", attempt.usage_json("completed")["bound"]["capped_by"].asString());
    // without finish_by_ms an unmanaged attempt is not enforced, as before
    Attempt plain(2, std::chrono::system_clock::now(), settings, "inst", true);
    plain.set_ids(AttemptIds());
    EXPECT_FALSE(plain.enforced());
}
