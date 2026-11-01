#include "gtest/gtest.h"

#include <atomic>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <json/json.h>

#include "cli/server_checks.hpp"
#include "cli/server_envelope.hpp"
#include "cli/server_fanout.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "common/threads/threading.hpp"


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
        return [this](std::function<void()> task, SteadyClock::time_point) {
            std::lock_guard<std::mutex> lock(mutex_);
            threads_.emplace_back(std::move(task));
            return true;
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

// A fan-out of one pair runs it inline: nothing is handed over
const Launch not_launched = [](std::function<void()>, SteadyClock::time_point) {
    ADD_FAILURE() << "one pair runs inline";
    return true;
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
// GPT A7 on every clock: no accepted finish_by_ms (up to 2^53 - 1) nor budget (the operator's
// cap of /pattern's, up to 2^53 - 1 ms) reaches a clock past what it holds (the steady clock
// counts nanoseconds, as Linux's system clock does): either binds nothing beyond the transport
TEST(ServerFanout, NoFarInstantOverflowsTheClocks) {
    const double far = static_cast<double>((uint64_t(1) << 53) - 1);
    for (const double budget : { far, 1e300, 1e13, 899'001.0 }) {
        for (const bool with_finish : { false, true }) {
            RequestHandles h;
            if (with_finish)
                h.finish_by_ms = (uint64_t(1) << 53) - 1;
            const RequestDeadlines d = deadlines_now(budget, h);
            const double effective_ms = std::chrono::duration<double, std::milli>(
                    d.effective - d.arrival).count();
            EXPECT_NEAR(899'000.0, effective_ms, 1.0) << budget << " " << with_finish;
            EXPECT_TRUE(d.transport_bound) << budget << " " << with_finish;
            EXPECT_LT(d.cutoff, d.effective);
            EXPECT_GT(d.cutoff, d.arrival);
            EXPECT_EQ(iso_utc(d.wall_of(d.arrival + std::chrono::milliseconds(899'000))),
                      iso_utc(d.effective_wall));
        }
    }
}

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
        [&](std::function<void()>, SteadyClock::time_point) { launched = true; return true; });
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

// Samples 03 and 16 name a not_reached or unfinished pair as the code does (contract F6 of
// the release-10 review): "(name, graph_path, annotation_path)", the entry's own tags, since a
// name may hold several pairs
TEST(ServerFanout, SamplesNameThePairAsTheFanOutDoes) {
    size_t checked = 0;
    for (const char *name : { "03_resolve_envelope_cutoff.json",
                              "16_pattern_envelope_not_reached.json" }) {
        const Json::Value body = sample(name)["body"];
        const std::string route = body["route"].asString().substr(1);
        for (const Json::Value &entry : body["answers"]) {
            if (entry["outcome"].asString() != "refused")
                continue;
            const std::string code = entry["refusal"]["code"].asString();
            if (code != "not_reached" && code != "unfinished")
                continue;
            const std::string label = "(" + entry["graph"].asString() + ", "
                    + entry["graph_path"].asString() + ", "
                    + entry["annotation_path"].asString() + ")";
            EXPECT_EQ(0u, entry["refusal"]["error"].asString().find(route + ": pair " + label + " "))
                << name << ": " << entry["refusal"]["error"].asString();
            ++checked;
        }
    }
    EXPECT_EQ(3u, checked);
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
        not_launched);
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
    }, not_launched);
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
    }, not_launched);
    EXPECT_TRUE(loaded);
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    wait_for_workers();
}

// A worker that reaches a pair's start after the cutoff, before the handler closed the request
// (GPT A2: an immediately granted reservation asks no wait predicate, and the load then cannot
// be interrupted): the cutoff is read with the transition, so nothing starts, a load or the
// work, and the pair is not_reached. The worker runs inside the handler's client check, where
// the handler has read the clock and holds no lock: the order that lost the race
TEST(ServerFanout, NoPooledPairStartsPastTheCutoff) {
    for (const bool load : { true, false }) {
        SCOPED_TRACE(load);
        const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
        auto entry = entry_of(pairs);
        const RequestDeadlines d = deadlines_now(60'000);
        auto time = std::make_shared<std::atomic<SteadyClock::rep>>(
                d.arrival.time_since_epoch().count());
        std::vector<std::function<void()>> tasks;
        size_t asked = 0;
        size_t started = 0;
        FanoutSettings settings = settings_of(pairs, d, "/pattern");
        settings.poll_ms = 1;
        settings.clock = [time]() {
            return SteadyClock::time_point(SteadyClock::duration(time->load()));
        };
        settings.gone = [&]() {
            if (++asked == 2) {
                EXPECT_EQ(2u, tasks.size());
                time->store((d.cutoff + 1ms).time_since_epoch().count());
                for (auto &task : tasks) {
                    task();
                }
            }
            return false;
        };
        FanoutOutcome out = run_fanout(settings, entry, [&](PairControl &control) {
            PairResult r;
            if (load) {
                // (the reservation was granted at once: no wait polled the predicate)
                control.waiting_memory();
            }
            r.not_started = !(load ? control.start_loading(100) : control.start_running());
            started += !r.not_started;
            return r;
        }, [&](std::function<void()> task, SteadyClock::time_point) {
            tasks.push_back(std::move(task));
            return true;
        });
        EXPECT_EQ(0u, started);
        ASSERT_EQ(2u, out.entries.size());
        for (size_t i = 0; i < 2; ++i) {
            const Json::Value &refusal = out.entries[i]["refusal"];
            EXPECT_EQ("not_reached", refusal["code"].asString());
            EXPECT_EQ("cutoff", refusal["reason"].asString());
            EXPECT_NE(std::string::npos, refusal["error"].asString().find(
                    load ? "it waited for the memory to load its index"
                         : "it waited for a pool thread")) << refusal["error"].asString();
            EXPECT_EQ("never_started", entry->json(0)["pairs"][int(i)]["state"].asString());
        }
    }
    wait_for_workers();
}

// The cutoff is no limit of an inline pair (section 4): it starts and answers as a single-graph
// request does, past the cutoff too
TEST(ServerFanout, AnInlinePairStartsPastTheCutoff) {
    const std::vector<SelectedPair> pairs { kPairs[0] };
    auto entry = entry_of(pairs);
    // (the budget is 0: the cutoff passed before the handler started)
    FanoutOutcome out = run_fanout(settings_of(pairs, deadlines_now(0)), entry,
        [&](PairControl &control) {
            EXPECT_TRUE(control.start_running());
            return answered(control);
        }, not_launched);
    ASSERT_EQ(1u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    wait_for_workers();
}

// The spec states what an inline pair does (GPT A-Q1, A-Q2), and the fan-out does it: started,
// it is answered past the cutoff (its own deadline ends it, never the cutoff), and its entry is
// the request's as soon as it ends; before it starts, an in_ram load past its wait limit is not
// started and the pair is refused 503 not_reached (NoLoadPastTheWaitLimit)
TEST(ServerFanout, TheSpecStatesTheInlinePairRule) {
    const std::vector<SelectedPair> pairs { kPairs[0] };
    auto entry = entry_of(pairs);
    // the cutoff 10 ms after the arrival; the pair works 60 ms past it
    const RequestDeadlines d = deadlines_now(60);
    FanoutOutcome out = run_fanout(settings_of(pairs, d), entry, [&](PairControl &control) {
        EXPECT_TRUE(control.start_running());
        std::this_thread::sleep_for(70ms);
        EXPECT_FALSE(control.stopped());
        return answered(control);
    }, not_launched);
    ASSERT_EQ(1u, out.entries.size());
    EXPECT_EQ("answered", out.entries[0]["outcome"].asString());
    EXPECT_GT(SteadyClock::now(), d.cutoff + 50ms);
    wait_for_workers();

    std::ifstream in("../docs/SPEC-multigraph-envelope.md", std::ios::binary);
    ASSERT_TRUE(in.good()) << "the spec not found (run from the build directory)";
    std::string spec;
    bool space = false;
    for (char c; in.get(c); ) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            space = true;
            continue;
        }
        if (space && spec.size())
            spec += ' ';
        space = false;
        spec += c;
    }
    for (const char *clause : {
            "A one-pair request (the selection holds one pair) runs inline, as its single-graph "
            "route does.",
            "is capped by `finish_by_ms` - clock_skew_allowance_ms, never by the work cutoff",
            "a load in progress is not interrupted; once the pair has started it is never "
            "`not_reached` or `unfinished` by the cutoff",
            "its envelope, answered or refused, is written as soon as the pair ends",
            "nothing is loaded and its entry is refused 503 `not_reached`",
            "`metagraph query --query-mode labels`",
            "so it is not reachable over HTTP" }) {
        EXPECT_NE(std::string::npos, spec.find(clause)) << clause;
    }
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
        [&](std::function<void()> task, SteadyClock::time_point) {
            if (launched++ == 0) {
                threads.emplace_back(std::move(task));
            } else {
                std::lock_guard<std::mutex> lock(held_mutex);
                held.push_back(std::move(task));
            }
            return true;
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

// The graphs' pool with its bounded queue full behind a blocked worker: what a fan-out hands
// over is |launch|'s in the server (ThreadPool::try_enqueue_until), whose wait must leave the
// handler its cutoff and its client
class FullPool {
  public:
    explicit FullPool(size_t room_left = 0) : pool_(1, kQueue) {
        std::promise<void> running;
        pool_.enqueue([this, &running]() {
            running.set_value();
            gate_.wait();
        });
        running.get_future().wait();
        for (size_t i = room_left; i < kQueue; ++i) {
            pool_.enqueue([]() {});
        }
    }
    ~FullPool() {
        gate_.open();
        pool_.join();
    }
    void release() { gate_.open(); }
    Launch launch() {
        return [this](std::function<void()> task, SteadyClock::time_point until) {
            return pool_.try_enqueue_until(until, std::move(task));
        };
    }

  private:
    static constexpr size_t kQueue = 20;
    Gate gate_;
    ThreadPool pool_;
};

// A cutoff at a full queue: the handler answers at the cutoff, every pair not_reached by the
// reason of a queued one, nothing handed over, and its client is asked while it waits (GPT A1:
// run_fanout blocked in the pool's enqueue past the cutoff, until a worker freed a slot)
TEST(ServerFanout, AFullPoolLeavesTheHandlerItsCutoff) {
    FullPool pool;
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto entry = entry_of(pairs);
    std::atomic<size_t> jobs { 0 };
    std::atomic<size_t> asked { 0 };
    // (the deadline less the reserve of 50 ms)
    FanoutSettings settings = settings_of(pairs, deadlines_now(130), "/pattern");
    settings.poll_ms = 5;
    settings.gone = [&]() { ++asked; return false; };
    const auto start = SteadyClock::now();
    std::promise<FanoutOutcome> answer;
    std::thread handler([&]() {
        answer.set_value(run_fanout(settings, entry, [&](PairControl &control) {
            ++jobs;
            return answered(control);
        }, pool.launch()));
    });
    auto future = answer.get_future();
    const bool returned = future.wait_for(30s) == std::future_status::ready;
    const auto elapsed = SteadyClock::now() - start;
    pool.release();
    handler.join();
    ASSERT_TRUE(returned);
    // the cutoff is 80 ms from the start (a scheduler's delay aside, well before the worker
    // was released: nothing else could free the queue)
    EXPECT_GE(elapsed, 70ms);
    EXPECT_LT(elapsed, 20s);
    EXPECT_GT(asked.load(), 1u);
    FanoutOutcome out = future.get();
    ASSERT_EQ(2u, out.entries.size());
    for (const Json::Value &e : out.entries) {
        EXPECT_EQ("refused", e["outcome"].asString());
        EXPECT_EQ("not_reached", e["refusal"]["code"].asString());
        EXPECT_EQ("cutoff", e["refusal"]["reason"].asString());
        EXPECT_NE(std::string::npos,
                  e["refusal"]["error"].asString().find("waited for a pool thread"));
    }
    wait_for_workers();
    EXPECT_EQ(0u, jobs.load());
    EXPECT_EQ("never_started", entry->json(0)["pairs"][0]["state"].asString());
    EXPECT_EQ("never_started", entry->json(0)["pairs"][1]["state"].asString());
    entry->response(200, 1);
    EXPECT_TRUE(entry->ended());
}

// A client that leaves while the queue is full stops the fan-out as it does otherwise
TEST(ServerFanout, AFullPoolLeavesTheHandlerItsClient) {
    FullPool pool;
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto entry = entry_of(pairs);
    std::atomic<size_t> asked { 0 };
    FanoutSettings settings = settings_of(pairs, deadlines_now(600'000), "/pattern");
    settings.poll_ms = 5;
    settings.gone = [&]() { return ++asked > 3; };
    std::promise<FanoutOutcome> answer;
    std::thread handler([&]() {
        answer.set_value(run_fanout(settings, entry, [&](PairControl &control) {
            return answered(control);
        }, pool.launch()));
    });
    auto future = answer.get_future();
    const bool returned = future.wait_for(30s) == std::future_status::ready;
    pool.release();
    handler.join();
    ASSERT_TRUE(returned);
    EXPECT_TRUE(future.get().client_gone);
    wait_for_workers();
    entry->response(0, std::nullopt);
    EXPECT_TRUE(entry->ended());
    EXPECT_EQ("client_gone", entry->json(0)["reason"].asString());
}

// A queue that frees up before the cutoff delays the pairs, and loses none: those that did not
// fit at first are handed over as room appears, in their order
TEST(ServerFanout, PairsWaitingForRoomAreHandedOver) {
    FullPool pool(/* room_left */ 1);
    const std::vector<SelectedPair> pairs = kPairs;
    auto entry = entry_of(pairs);
    FanoutSettings settings = settings_of(pairs, deadlines_now(600'000), "/pattern");
    settings.poll_ms = 5;
    std::promise<FanoutOutcome> answer;
    std::thread handler([&]() {
        answer.set_value(run_fanout(settings, entry, [&](PairControl &control) {
            EXPECT_TRUE(control.start_running());
            return answered(control);
        }, pool.launch()));
    });
    auto future = answer.get_future();
    // the handler keeps waiting for a full queue
    EXPECT_EQ(std::future_status::timeout, future.wait_for(30ms));
    pool.release();
    handler.join();
    FanoutOutcome out = future.get();
    ASSERT_EQ(3u, out.entries.size());
    for (const Json::Value &e : out.entries) {
        EXPECT_EQ("answered", e["outcome"].asString());
    }
    wait_for_workers();
}

// Pairs that end while the handler hands them over (its lock released, their notices sent to
// nobody) are seen at once, not after a whole poll (recheck4: the handler slept poll_ms with
// every pair resolved)
TEST(ServerFanout, PairsEndedDuringTheHandOverAreSeenAtOnce) {
    const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
    auto entry = entry_of(pairs);
    FanoutSettings settings = settings_of(pairs, deadlines_now(600'000), "/pattern");
    settings.poll_ms = 10'000;
    const auto start = SteadyClock::now();
    FanoutOutcome out = run_fanout(settings, entry, [&](PairControl &control) {
        EXPECT_TRUE(control.start_running());
        return answered(control);
    }, [](std::function<void()> task, SteadyClock::time_point) {
        // the pair runs to its end before the launch returns
        std::thread(std::move(task)).join();
        return true;
    });
    const auto elapsed = SteadyClock::now() - start;
    EXPECT_LT(elapsed, 5s) << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count()
                           << " ms";
    ASSERT_EQ(2u, out.entries.size());
    for (const Json::Value &e : out.entries) {
        EXPECT_EQ("answered", e["outcome"].asString());
    }
    wait_for_workers();
}

// The same on the graphs' pool as the server uses it: no fan-out of trivial pairs waits a
// whole poll after its pairs ended (recheck4: 1 in 300 did)
TEST(ServerFanout, NoFanOutOnARealPoolWaitsAPollAfterItsPairs) {
    ThreadPool pool(4, 1000);
    const Launch launch = [&pool](std::function<void()> task, SteadyClock::time_point until) {
        return pool.try_enqueue_until(until, std::move(task));
    };
    size_t slow = 0;
    for (size_t r = 0; r < 300; ++r) {
        const std::vector<SelectedPair> pairs { kPairs[0], kPairs[1] };
        auto entry = entry_of(pairs);
        FanoutSettings settings = settings_of(pairs, deadlines_now(600'000), "/pattern");
        settings.poll_ms = 2'000;
        const auto start = SteadyClock::now();
        run_fanout(settings, entry, [&](PairControl &control) {
            control.start_running();
            return answered(control);
        }, launch);
        if (SteadyClock::now() - start >= 1s)
            ++slow;
    }
    pool.join();
    wait_for_workers();
    EXPECT_EQ(0u, slow);
}

// An exception of the submission (GPT A4: a task's allocation, the second of two): the first
// pair is told to stop and the other ends at once, so that the entry ends with the response,
// its id free again after the retention (none here), instead of "running" for good. The
// exception reaches the caller, which answers its error
TEST(ServerFanout, ASubmissionThatThrowsEndsTheRequest) {
    for (const bool first_started : { true, false }) {
        SCOPED_TRACE(first_started);
        const std::vector<SelectedPair> pairs = kPairs;
        RequestSettings rs;
        rs.server_instance = "inst";
        rs.retention_s = 0;
        rs.retention_count = 0;
        RequestRegistry registry(rs);
        auto entry = std::make_shared<RequestEntry>(std::string("r"), "/pattern", "inst",
                                                    std::chrono::system_clock::now(),
                                                    refs_of(pairs));
        ASSERT_FALSE(registry.start(entry));
        Gate running;
        std::atomic<bool> told { false };
        std::vector<std::function<void()>> queued;
        std::thread worker;
        size_t offered = 0;
        EXPECT_THROW(run_fanout(settings_of(pairs, deadlines_now(600'000), "/pattern"), entry,
            [&](PairControl &control) {
                EXPECT_TRUE(control.start_running());
                running.open();
                while (!control.stopped()) {
                    std::this_thread::sleep_for(1ms);
                }
                told = true;
                PairResult r;
                r.aborted = true;
                return r;
            },
            [&](std::function<void()> task, SteadyClock::time_point) -> bool {
                if (offered++ == 1)
                    throw std::bad_alloc();
                if (first_started) {
                    worker = std::thread(std::move(task));
                    running.wait();
                } else {
                    queued.push_back(std::move(task));
                }
                return true;
            }), std::bad_alloc);
        EXPECT_EQ(2u, offered);
        if (worker.joinable())
            worker.join();
        EXPECT_EQ(first_started, told.load());
        // the pairs never handed over are ended; one still queued ends with its task
        Json::Value state = registry.state("r").second;
        EXPECT_EQ("running", state["state"].asString());
        EXPECT_EQ(first_started ? "ended" : "never_started",
                  state["pairs"][0]["state"].asString());
        EXPECT_EQ(first_started ? "unfinished" : "not_reached",
                  state["pairs"][0]["outcome"].asString());
        for (int i : { 1, 2 }) {
            EXPECT_EQ("never_started", state["pairs"][i]["state"].asString());
            EXPECT_EQ("not_reached", state["pairs"][i]["outcome"].asString());
        }
        for (auto &task : queued) {
            task();
        }
        queued.clear();
        // what the handler answers: an error, and the request has ended
        entry->response(500, 10);
        EXPECT_TRUE(entry->ended());
        EXPECT_EQ("error", entry->json(0)["reason"].asString());
        EXPECT_EQ(404, registry.state("r").first);
        EXPECT_FALSE(registry.start(std::make_shared<RequestEntry>(
                std::string("r"), "/pattern", "inst", std::chrono::system_clock::now(),
                refs_of(pairs))));
        wait_for_workers();
    }
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

// A request refused 400 after its admission (only its work could find the fault) frees its id:
// nothing is retained, a retry runs, and a late worker of the freed entry touches no newer
// request of the id. Any other status keeps the entry
TEST(ServerFanout, A400FreesTheRequestId) {
    RequestSettings rs;
    rs.server_instance = "inst";
    rs.retention_s = 3600;
    RequestRegistry registry(rs);
    auto entry = [&](const std::string &id) {
        return std::make_shared<RequestEntry>(id, "/resolve", "inst",
                                              std::chrono::system_clock::now(),
                                              refs_of({ kPairs[0] }));
    };
    auto first = entry("r");
    ASSERT_FALSE(registry.start(first));
    first->pair_running(0);
    EXPECT_EQ(200, registry.state("r").first);
    EXPECT_TRUE(registry.start(entry("r")));
    // (the route's order: the id, then the response)
    registry.forget(first);
    first->response(400, 20);
    EXPECT_EQ(404, registry.state("r").first);
    EXPECT_EQ(0u, registry.size());

    auto retry = entry("r");
    EXPECT_FALSE(registry.start(retry));
    // the worker of the freed request ends now: nothing of it reaches the registry
    first->pair_outcome(0, "refused");
    first->pair_ended(0);
    EXPECT_TRUE(first->ended());
    EXPECT_EQ("running", registry.state("r").second["state"].asString());
    // forgetting the old entry again leaves the new one
    registry.forget(first);
    EXPECT_EQ(200, registry.state("r").first);
    retry->pair_running(0);
    retry->pair_ended(0);
    retry->response(200, 10);
    EXPECT_EQ("answered", registry.state("r").second["reason"].asString());

    // every other status stays retained (reason error, deadline)
    for (const auto &[id, status, reason] : { std::make_tuple("e500", 500, "error"),
                                              std::make_tuple("e404", 404, "error"),
                                              std::make_tuple("e503", 503, "deadline") }) {
        auto e = entry(id);
        ASSERT_FALSE(registry.start(e));
        e->pair_never_started(0);
        e->response(status, 10);
        EXPECT_EQ(reason, registry.state(id).second["reason"].asString());
        EXPECT_TRUE(registry.start(entry(id)));
    }

    // a request without an id is held by nobody
    auto anonymous = std::make_shared<RequestEntry>(std::nullopt, "/resolve", "inst",
                                                    std::chrono::system_clock::now(),
                                                    refs_of({ kPairs[0] }));
    registry.forget(anonymous);
}

// request_id and attempt_id are one namespace (GET /requests/{id} looks up both): an id held by
// either registry refuses the other route, running or retained, in that route's 409, until its
// retention ended; a repeat on one route answers as before (GPT A3: a retained /resolve "x" hid
// a running /traverse attempt "x")
class IdNamespace : public ::testing::Test {
  protected:
    IdNamespace() : requests_(request_settings()), attempts_(attempt_settings()),
                    admission_(requests_, attempts_) {}

    RequestSettings request_settings() {
        RequestSettings rs;
        rs.server_instance = "inst";
        rs.retention_s = 60;
        rs.clocks.steady = [this]() { return now(); };
        return rs;
    }
    AttemptSettings attempt_settings() {
        AttemptSettings as;
        as.retention_s = 60;
        as.retention_count = 100;
        as.clock = [this]() { return now(); };
        return as;
    }
    SteadyClock::time_point now() const { return base_ + std::chrono::seconds(seconds_.load()); }
    void pass(int64_t seconds) { seconds_ += seconds; }

    std::shared_ptr<RequestEntry> entry(const std::string &id) {
        RequestEntry::Clocks clocks;
        clocks.steady = [this]() { return now(); };
        return std::make_shared<RequestEntry>(id, "/resolve", "inst",
                                              std::chrono::system_clock::now(),
                                              refs_of({ kPairs[0] }), clocks);
    }
    std::shared_ptr<Attempt> attempt(const std::string &id) {
        auto a = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(),
                                           attempts_.settings(), attempts_.server_instance(),
                                           true, nullptr);
        AttemptIds ids;
        ids.attempt_id = id;
        a->set_ids(ids);
        return a;
    }
    void end(const std::shared_ptr<RequestEntry> &e) {
        e->pair_running(0);
        e->pair_outcome(0, "answered");
        e->pair_ended(0, 1, 0.5);
        e->response(200, 10);
    }

    SteadyClock::time_point base_ = SteadyClock::now();
    std::atomic<int64_t> seconds_ { 0 };
    RequestRegistry requests_;
    AttemptRegistry attempts_;
    IdAdmission admission_;
};

TEST_F(IdNamespace, ARequestHoldsItsIdAgainstAnAttempt) {
    auto e = entry("x");
    ASSERT_FALSE(admission_.start_request(e));
    // running
    auto refused = admission_.start_attempt(attempt("x"));
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_FALSE(refused->instance_mismatch);
    EXPECT_EQ("x", refused->body["attempt_id"].asString());
    EXPECT_EQ("running", refused->body["state"].asString());
    EXPECT_EQ("/resolve", refused->body["route"].asString());
    EXPECT_EQ(404, attempts_.state("x").first);
    // retained
    end(e);
    refused = admission_.start_attempt(attempt("x"));
    ASSERT_TRUE(refused);
    EXPECT_EQ("finished", refused->body["state"].asString());
    EXPECT_EQ("answered", refused->body["reason"].asString());
    EXPECT_EQ(404, attempts_.state("x").first);
    // past its retention: the attempt runs
    pass(60);
    auto a = attempt("x");
    EXPECT_FALSE(admission_.start_attempt(a));
    EXPECT_EQ(200, attempts_.state("x").first);
    EXPECT_EQ(404, requests_.state("x").first);
    // and holds the id in its turn
    EXPECT_TRUE(admission_.start_request(entry("x")));
}

TEST_F(IdNamespace, AnAttemptHoldsItsIdAgainstARequest) {
    auto a = attempt("y");
    ASSERT_FALSE(admission_.start_attempt(a));
    // running: the attempt as GET /requests/{id} answers it, with the route's error
    auto e = entry("y");
    bool asked = false;
    auto refused = admission_.start_request(e, [&]() -> std::optional<Json::Value> {
        asked = true;
        return std::nullopt;
    });
    ASSERT_TRUE(refused);
    EXPECT_FALSE(asked);
    EXPECT_EQ("y", (*refused)["request_id"].asString());
    EXPECT_EQ("/traverse", (*refused)["route"].asString());
    EXPECT_EQ("running", (*refused)["state"].asString());
    EXPECT_EQ("y", (*refused)["attempt"]["attempt_id"].asString());
    EXPECT_EQ(requests_.repeated_text("y"), (*refused)["error"].asString());
    EXPECT_EQ(404, requests_.state("y").first);
    // retained
    attempts_.finish(a, "completed", 200, 10);
    refused = admission_.start_request(entry("y"));
    ASSERT_TRUE(refused);
    EXPECT_EQ("ended", (*refused)["state"].asString());
    EXPECT_EQ(404, requests_.state("y").first);
    // past its retention: the request runs
    pass(60);
    EXPECT_FALSE(admission_.start_request(e));
    EXPECT_EQ(200, requests_.state("y").first);
    EXPECT_EQ(404, attempts_.state("y").first);
    EXPECT_TRUE(admission_.start_attempt(attempt("y")));
}

// a repeat on one route answers as before: its own registry's body; a tombstone of a cancel holds
// an attempt_id only (no request ran with it, GET /requests does not answer it)
TEST_F(IdNamespace, ARepeatOnOneRouteAnswersAsBefore) {
    auto e = entry("r");
    ASSERT_FALSE(admission_.start_request(e));
    auto repeated = admission_.start_request(entry("r"));
    auto direct = requests_.start(entry("r"));
    ASSERT_TRUE(repeated);
    ASSERT_TRUE(direct);
    EXPECT_EQ(*direct, *repeated);

    auto a = attempt("t");
    ASSERT_FALSE(admission_.start_attempt(a));
    auto again = admission_.start_attempt(attempt("t"));
    ASSERT_TRUE(again);
    EXPECT_EQ("t", again->body["attempt_id"].asString());
    EXPECT_FALSE(again->body.isMember("route"));

    attempts_.cancel("cancelled", 0);
    EXPECT_FALSE(admission_.start_request(entry("cancelled")));
    ASSERT_TRUE(admission_.start_attempt(attempt("cancelled")));
}

// an admission's check of the other registry, whatever the order: never both
TEST_F(IdNamespace, NoWindowAdmitsBoth) {
    const size_t rounds = 2000;
    std::atomic<size_t> go { 0 };
    std::vector<std::shared_ptr<RequestEntry>> entries;
    std::vector<std::shared_ptr<Attempt>> tries;
    for (size_t i = 0; i < rounds; ++i) {
        entries.push_back(entry("id" + std::to_string(i)));
        tries.push_back(attempt("id" + std::to_string(i)));
    }
    std::vector<int> request_won(rounds, 0);
    std::vector<int> attempt_won(rounds, 0);
    auto wait_turn = [&](size_t i) {
        ++go;
        while (go.load() < 2 * (i + 1)) {
            std::this_thread::yield();
        }
    };
    std::thread a([&]() {
        for (size_t i = 0; i < rounds; ++i) {
            wait_turn(i);
            request_won[i] = !admission_.start_request(entries[i]);
        }
    });
    std::thread b([&]() {
        for (size_t i = 0; i < rounds; ++i) {
            wait_turn(i);
            attempt_won[i] = !admission_.start_attempt(tries[i]);
        }
    });
    a.join();
    b.join();
    for (size_t i = 0; i < rounds; ++i) {
        EXPECT_EQ(1, request_won[i] + attempt_won[i]) << i;
    }
}

// One server answers each sample and one execution is described in it: every server_instance in
// a sample (the envelope's, an attempt's, its usage's, a capabilities block's) is the same, and
// so is every attempt_id and request_id (the request's, its attempt's, the attempt's usage). A
// reader that checks identity must accept every sample (GPT review 3, finding 10).
TEST(ServerFanout, EachSampleNamesOneServerAndOneExecution) {
    const std::string dir = "../tests/data/multigraph_contract";
    std::vector<std::string> names;
    for (const auto &file : std::filesystem::directory_iterator(dir)) {
        if (file.path().extension() == ".json")
            names.push_back(file.path().filename().string());
    }
    ASSERT_LE(22u, names.size()) << dir << " (run from the build directory)";
    for (const std::string &name : names) {
        std::set<std::string> instances;
        std::set<std::string> ids;
        std::function<void(const Json::Value &, const std::string &)> walk
                = [&](const Json::Value &v, const std::string &path) {
            if (v.isObject()) {
                for (const std::string &key : v.getMemberNames()) {
                    const Json::Value &child = v[key];
                    if (child.isString() && key == "server_instance")
                        instances.insert(child.asString());
                    if (child.isString() && (key == "attempt_id" || key == "request_id"))
                        ids.insert(child.asString());
                    walk(child, path + "." + key);
                }
            } else if (v.isArray()) {
                for (Json::ArrayIndex i = 0; i < v.size(); ++i) {
                    walk(v[i], path + "[" + std::to_string(i) + "]");
                }
            }
        };
        walk(sample(name), "");
        EXPECT_GE(1u, instances.size()) << name << ": " << instances.size() << " server instances";
        EXPECT_GE(1u, ids.size()) << name << ": " << ids.size() << " execution ids";
    }
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
