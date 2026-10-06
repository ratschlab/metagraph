#include "gtest/gtest.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <arpa/inet.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <json/json.h>
#include <zlib.h>

#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "cli/load/load_annotation.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace {

using namespace mtg::cli;
using mtg::graph::traversal::AttemptAborted;
using mtg::graph::traversal::ExternalStop;
using Clock = Attempt::Clock;

// A clock the test moves: the registry's retention and an attempt's bound read it. The wall
// clock (not_after_ms, a tombstone's wall-clock hold) moves with it, |wall_step| ms apart (a
// test steps the wall clock alone by changing that)
struct FakeClock {
    Clock::time_point base = Clock::now();
    std::atomic<int64_t> ms { 0 };
    std::atomic<int64_t> wall_step { 0 };
    static constexpr int64_t kWallBase = 1791137002417;   // 2026-10-04T18:03:22.417Z
    std::function<Clock::time_point()> fn() {
        return [this]() { return base + std::chrono::milliseconds(ms.load()); };
    }
    std::function<std::chrono::system_clock::time_point()> wall_fn() {
        return [this]() {
            return std::chrono::system_clock::time_point(
                    std::chrono::milliseconds(kWallBase + ms.load() + wall_step.load()));
        };
    }
    uint64_t wall_ms() const { return kWallBase + ms.load() + wall_step.load(); }
};

AttemptSettings settings_with(FakeClock *clock, uint64_t retention_s = 60,
                              size_t retention_count = 100) {
    AttemptSettings s;
    s.allowance_ms = 1000;
    s.hard_cap_ms = 899'000;
    s.retention_s = retention_s;
    s.retention_count = retention_count;
    s.client_check_ms = 0;
    s.poll_stride = 1;
    // a stop time below the floor (allowance / 2), so that the bounds these tests compute are
    // the floor's unless a test sets the reserve's parts (the default, 1000 since the
    // efficiency pass's calibration, would exceed this allowance's floor)
    s.delivery_stop_ms = 250;
    if (clock) {
        s.clock = clock->fn();
        s.wall_clock = clock->wall_fn();
    }
    return s;
}

std::shared_ptr<Attempt> attempt_of(const AttemptRegistry &registry, const std::string &id,
                                    std::function<bool()> peer_gone = nullptr) {
    // the epoch: the header was read now, on the registry's clock
    auto a = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(),
                                       registry.settings(), registry.server_instance(), true,
                                       std::move(peer_gone));
    AttemptIds ids;
    ids.attempt_id = id;
    a->set_ids(ids);
    return a;
}

} // namespace


TEST(GraphletAttempt, IdsAreValidated) {
    EXPECT_TRUE(valid_attempt_id("a"));
    EXPECT_TRUE(valid_attempt_id("run-7.locus_3:try-2"));
    EXPECT_TRUE(valid_attempt_id(std::string(128, 'x')));
    EXPECT_FALSE(valid_attempt_id(""));
    EXPECT_FALSE(valid_attempt_id(std::string(129, 'x')));
    EXPECT_FALSE(valid_attempt_id("a b"));
    EXPECT_FALSE(valid_attempt_id("a/b"));
    EXPECT_FALSE(valid_attempt_id("\xc3\xa9"));

    Json::Value r;
    EXPECT_TRUE(attempt_ids(r).attempt_id.empty());
    r["attempt_id"] = "a1";
    r["budget_id"] = "b1";
    r["locus_id"] = "l1";
    const AttemptIds ids = attempt_ids(r);
    EXPECT_EQ("a1", ids.attempt_id);
    EXPECT_EQ("b1", ids.budget_id);
    EXPECT_EQ("l1", ids.locus_id);
    r["attempt_id"] = 7;
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    r["attempt_id"] = "a 1";
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    // an echo without the attempt it belongs to
    r.removeMember("attempt_id");
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    r.removeMember("budget_id");
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    r.removeMember("locus_id");
    EXPECT_NO_THROW(attempt_ids(r));
    // expect_server_instance (review of pass 5, finding 1): a token, with attempt_id only
    r["expect_server_instance"] = "0123456789abcdef";
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    r["attempt_id"] = "a1";
    EXPECT_EQ("0123456789abcdef", attempt_ids(r).expect_server_instance);
    r["expect_server_instance"] = 7;
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    r["expect_server_instance"] = "";
    EXPECT_THROW(attempt_ids(r), InvalidRequest);
    // declared by strict parsing
    Json::Value req;
    req["seeds"][0]["sequence"] = std::string(40, 'A');
    req["strategy"] = Json::Value(Json::objectValue);
    req["attempt_id"] = "a1";
    req["expect_server_instance"] = "0123456789abcdef";
    EXPECT_NO_THROW(parse_traverse_request(req));
}

// The bound is n_seeds x T + allowance on the attempt's clock, capped by the HTTP server's;
// the seeds stop being walked at bound - allowance / 2, and what is built or written stops at
// the bound itself
TEST(GraphletAttempt, BoundIsEnforcedOnTheAttemptsClock) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    auto a = attempt_of(registry, "bound");
    a->set_bound(4, 1000);
    EXPECT_EQ(5000, a->bound_ms());
    EXPECT_EQ(4500, a->walk_until_ms());
    EXPECT_TRUE(a->enforced());
    clock.ms = 4499;
    EXPECT_EQ(ExternalStop::NONE, a->poll());
    EXPECT_NO_THROW(a->check_delivery());
    clock.ms = 4500;
    EXPECT_EQ(ExternalStop::ATTEMPT_DEADLINE, a->poll());
    EXPECT_STREQ("deadline", a->walk_reason());
    // the response may still be written until the bound
    EXPECT_NO_THROW(a->check_delivery());
    clock.ms = 5000;
    EXPECT_THROW(a->check_delivery(), AttemptAtBound);
    const Json::Value usage = a->usage_json("deadline");
    EXPECT_EQ(5000u, usage["bound_ms"].asUInt64());
    EXPECT_EQ(4u, usage["bound"]["seeds"].asUInt64());
    EXPECT_EQ(1000.0, usage["bound"]["time_budget_ms"].asDouble());
    EXPECT_EQ(1000u, usage["bound"]["allowance_ms"].asUInt64());
    EXPECT_EQ(4500u, usage["bound"]["walk_until_ms"].asUInt64());
    EXPECT_TRUE(usage["bound"]["capped_by"].isNull());
    EXPECT_TRUE(usage["bound"]["enforced"].asBool());

    // capped by the content timeout
    auto big = attempt_of(registry, "big");
    big->set_bound(10'000, 1000);
    EXPECT_EQ(899'000, big->bound_ms());
    EXPECT_EQ("content_timeout", big->usage_json("completed")["bound"]["capped_by"].asString());
    // a non-positive time budget walks nothing: the allowance alone
    auto zero = attempt_of(registry, "zero");
    zero->set_bound(3, 0);
    EXPECT_EQ(1000, zero->bound_ms());

    // without an attempt_id the bound is not enforced (the request keeps today's form)
    auto plain = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(),
                                           registry.settings(), registry.server_instance(), true);
    plain->set_bound(1, 1);
    clock.ms = 1'000'000;
    EXPECT_FALSE(plain->enforced());
    EXPECT_EQ(ExternalStop::NONE, plain->poll());
    EXPECT_NO_THROW(plain->check_delivery());

    // the CLI: no cap, nothing enforced, an infinite bound is stated as null
    AttemptSettings cli;
    cli.hard_cap_ms = 0;
    Attempt local(0, std::chrono::system_clock::now(), cli, "cli", false);
    AttemptIds ids;
    ids.attempt_id = "cli";
    local.set_ids(ids);
    local.set_bound(2, std::numeric_limits<double>::infinity());
    EXPECT_FALSE(local.enforced());
    EXPECT_TRUE(local.usage_json("completed")["bound_ms"].isNull());
    EXPECT_FALSE(local.usage_json("completed")["bound"]["enforced"].asBool());
}

// A client that is gone stops the walk (CLIENT_GONE) and anything built after it
// (AttemptAborted), whatever stopped the walk first
TEST(GraphletAttempt, GoneClientAbandonsTheAttempt) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    bool gone = false;
    auto a = attempt_of(registry, "gone", [&]() { return gone; });
    a->set_bound(1, 1000);
    EXPECT_EQ(ExternalStop::NONE, a->poll());
    gone = true;
    EXPECT_EQ(ExternalStop::CLIENT_GONE, a->poll());
    EXPECT_THROW(a->check_delivery(), AttemptAborted);

    // after a cancel the walk stops as cancelled; a client that then goes is still gone
    bool gone2 = false;
    auto b = attempt_of(registry, "gone2", [&]() { return gone2; });
    b->set_bound(1, 1000);
    ASSERT_FALSE(registry.start(b));
    EXPECT_EQ(200, registry.cancel("gone2", 0).first);
    EXPECT_EQ(ExternalStop::CANCELLED, b->poll());
    EXPECT_NO_THROW(b->check_delivery());
    gone2 = true;
    EXPECT_THROW(b->check_delivery(), AttemptAborted);
}

// The review of the stage-4 backend, F5: a seed whose walk the client's departure cut is
// abandoned, not finished — in the usage and in the state kept after the attempt finished; a
// seed whose walk ended before its result was abandoned stays finished
TEST(GraphletAttempt, AbandonedWalksAreNotFinished) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    auto a = attempt_of(registry, "ab");
    a->set_bound(3, 1000);
    ASSERT_FALSE(registry.start(a));
    SeedUsage walked;
    walked.outcome = "complete";
    a->seed_started(0);
    a->seed_walked(0, walked);
    // seed 1: walked, then abandoned while its result was built
    a->seed_started(1);
    a->seed_walked(1, walked);
    SeedUsage cut;
    cut.outcome = "abandoned";
    cut.stopped_by = "client_gone";
    a->seed_walked(1, cut);
    // seed 2: abandoned in its walk
    a->seed_started(2);
    a->seed_walked(2, cut);
    const Json::Value seeds = a->usage_json("client_gone")["seeds"];
    EXPECT_EQ(3u, seeds["started"].asUInt64());
    EXPECT_EQ(2u, seeds["finished"].asUInt64());
    EXPECT_EQ(1u, seeds["abandoned"].asUInt64());
    registry.finish(a, "client_gone", 0, std::nullopt);
    const Json::Value state = registry.state("ab").second;
    EXPECT_EQ(seeds, state["seeds"]);
    EXPECT_EQ(seeds, state["usage"]["seeds"]);
    EXPECT_FALSE(state["usage"].isMember("per_seed"));
    EXPECT_NE(std::string::npos, a->stop_summary().find("1 abandoned"));
}

// The stop latency is measured where a seed's walk ends: its first seed_walked. A second call
// for the seed comes while its result is built (the client left, or a writer refused it), and
// the time since the walk-until then includes building, which the reserve prices apart. Taken
// as a stop, it became the server's measured_stop_ms, and later attempts whose bound was no
// longer walked nothing (review of 2026-10-06, U11-01: 10.5 s recorded after a 400 ms stop)
TEST(GraphletAttempt, StopLatencyIsMeasuredWhereTheWalkEnds) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    // the walk-until trips; the walk ends 400 ms later; the client leaves 10 s after that
    for (const char *second : { "abandoned", "failed" }) {
        clock.ms = 0;
        auto a = attempt_of(registry, std::string("late-") + second);
        a->set_bound(1, 10'000);
        ASSERT_FALSE(registry.start(a));
        const double until = a->walk_until_ms();
        a->seed_started(0);
        clock.ms = static_cast<int64_t>(std::ceil(until)) + 1;
        ASSERT_EQ(ExternalStop::ATTEMPT_DEADLINE, a->poll(true));
        clock.ms = static_cast<int64_t>(std::ceil(until)) + 400;
        SeedUsage walked;
        walked.outcome = "partial";
        walked.stopped_by = "attempt_deadline";
        a->seed_walked(0, walked);
        const double stop = clock.ms - until;
        EXPECT_NEAR(stop, a->own_stop_ms(), 1e-6) << second;
        clock.ms += 10'000;
        SeedUsage cut;
        cut.outcome = second;
        cut.stopped_by = std::string(second) == "abandoned" ? "client_gone" : "";
        a->seed_walked(0, cut);
        EXPECT_NEAR(stop, a->own_stop_ms(), 1e-6) << second;
        registry.note_stop_latency(a->own_stop_ms());
        registry.finish(a, "client_gone", 0, std::nullopt);
        EXPECT_NEAR(stop, registry.measured().stop_ms, 1e-6) << second;
    }
    // a walk that ended before its walk-until measures nothing, however late a second call
    // comes (the reviewer's variant B: 6,000 ms recorded for a walk that never passed it)
    AttemptRegistry other(settings_with(&clock));
    clock.ms = 0;
    auto b = attempt_of(other, "early");
    b->set_bound(1, 10'000);
    ASSERT_FALSE(other.start(b));
    const double until = b->walk_until_ms();
    b->seed_started(0);
    clock.ms = static_cast<int64_t>(until) - 2000;
    EXPECT_EQ(ExternalStop::NONE, b->poll(true));
    SeedUsage complete;
    complete.outcome = "complete";
    b->seed_walked(0, complete);
    EXPECT_EQ(0.0, b->own_stop_ms());
    clock.ms = static_cast<int64_t>(until) + 6000;
    SeedUsage cut;
    cut.outcome = "abandoned";
    cut.stopped_by = "client_gone";
    b->seed_walked(0, cut);
    EXPECT_EQ(0.0, b->own_stop_ms());
    other.note_stop_latency(b->own_stop_ms());
    EXPECT_EQ(0.0, other.measured().stop_ms);
    EXPECT_TRUE(other.capabilities_json()["delivery_reserve"]["measured_stop_ms"].isNull());
    // the counters are once per seed, as before
    EXPECT_EQ(1u, b->usage_json("client_gone")["seeds"]["finished"].asUInt64());
}

// not_after_ms (pass 5, W1): an integer in [0, 2^53 - 1], with or without attempt_id; a
// fraction, a sign, a string or a value no JSON reader keeps exactly is a 400 naming the field
TEST(GraphletAttempt, NotAfterMsIsParsedStrictly) {
    Json::Value r;
    r["not_after_ms"] = Json::UInt64(1759601000000ull);
    EXPECT_EQ(1759601000000ull, attempt_ids(r).not_after_ms.value());
    EXPECT_TRUE(attempt_ids(r).attempt_id.empty());
    r["not_after_ms"] = Json::UInt64(kMaxNotAfterMs);
    EXPECT_EQ(kMaxNotAfterMs, attempt_ids(r).not_after_ms.value());
    r["not_after_ms"] = 0;
    EXPECT_EQ(0u, attempt_ids(r).not_after_ms.value());
    // a JSON number written with an exponent but without a fraction is the integer it denotes
    r["not_after_ms"] = 1.759601e12;
    EXPECT_EQ(1759601000000ull, attempt_ids(r).not_after_ms.value());
    for (const Json::Value &bad : { Json::Value(1759601000000.5), Json::Value(-1),
                                    Json::Value(Json::UInt64(kMaxNotAfterMs + 1)),
                                    Json::Value("1759601000000"), Json::Value(true),
                                    Json::Value(Json::nullValue) }) {
        r["not_after_ms"] = bad;
        try {
            attempt_ids(r);
            ADD_FAILURE() << "accepted " << bad.toStyledString();
        } catch (const InvalidRequest &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find("request.not_after_ms"));
        }
    }
    r.removeMember("not_after_ms");
    EXPECT_FALSE(attempt_ids(r).not_after_ms);

    // strict parsing of the whole request declares it (no unknown-field 400)
    Json::Value req;
    req["seeds"][0]["sequence"] = std::string(40, 'A');
    req["strategy"] = Json::Value(Json::objectValue);
    req["not_after_ms"] = Json::UInt64(1759601000000ull);
    EXPECT_EQ(1759601000000ull, parse_traverse_request(req).not_after_ms.value());
    req["not_after_ms"] = -5;
    EXPECT_THROW(parse_traverse_request(req), InvalidRequest);
}

// Refused at handler start once passed, strictly (no allowance added), with the exact body;
// nothing is registered, so the state route answers 404; an id that is running, retained or
// tombstoned is answered first (an "expired" 409 always means no attempt with the id exists)
TEST(GraphletAttemptRegistry, NotAfterMsRefusesAtStart) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock);
    std::atomic<int64_t> wall_ms { 1791137002417 };
    settings.wall_clock = [&]() {
        return std::chrono::system_clock::time_point(std::chrono::milliseconds(wall_ms.load()));
    };
    AttemptRegistry registry(settings);
    EXPECT_EQ(1791137002417u, registry.now_ms());

    // not passed: registered as usual, and echoed in usage and in the state
    auto a = attempt_of(registry, "a-17");
    AttemptIds ids = a->ids();
    ids.not_after_ms = 1791137002417;   // equal: not passed
    a->set_ids(ids);
    EXPECT_FALSE(registry.start(a));
    a->set_bound(1, 1000);
    EXPECT_EQ(1791137002417u, a->usage_json("completed")["not_after_ms"].asUInt64());
    EXPECT_EQ(1791137002417u, registry.state("a-17").second["not_after_ms"].asUInt64());

    // passed by one ms: refused, nothing registered
    auto b = attempt_of(registry, "a-18");
    ids.attempt_id = "a-18";
    ids.budget_id = "b-3";
    ids.locus_id = "l-9";
    ids.not_after_ms = 1791137002416;
    b->set_ids(ids);
    auto refused = registry.start(b);
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->expired);
    const Json::Value &body = refused->body;
    std::vector<std::string> keys = body.getMemberNames();
    EXPECT_EQ((std::vector<std::string>{ "attempt_id", "budget_id", "error", "locus_id",
                                         "not_after_ms", "server_instance", "server_time_ms",
                                         "state" }), keys);
    EXPECT_EQ("expired", body["state"].asString());
    EXPECT_EQ(1791137002416u, body["not_after_ms"].asUInt64());
    EXPECT_EQ(Json::uintValue, body["not_after_ms"].type());
    EXPECT_EQ(Json::uintValue, body["server_time_ms"].type());
    EXPECT_EQ(1791137002417u, body["server_time_ms"].asUInt64());
    EXPECT_EQ(registry.server_instance(), body["server_instance"].asString());
    EXPECT_NE(std::string::npos, body["error"].asString().find("2026-10-04T18:03:22.416Z"));
    EXPECT_FALSE(body.isMember("usage"));
    EXPECT_EQ(404, registry.state("a-18").first);
    // and a later copy is expired too: still nothing registered
    EXPECT_TRUE(registry.start(b)->expired);
    EXPECT_EQ(404, registry.state("a-18").first);

    // without attempt_id the same rule and body, without the ids
    AttemptIds plain;
    plain.not_after_ms = 1791137002416;
    EXPECT_TRUE(not_after_passed(plain, registry.now_ms()));
    EXPECT_FALSE(not_after_passed(plain, 1791137002416));
    const Json::Value pb = expired_json(plain, registry.now_ms(), registry.server_instance());
    keys = pb.getMemberNames();
    EXPECT_EQ((std::vector<std::string>{ "error", "not_after_ms", "server_instance",
                                         "server_time_ms", "state" }), keys);
    EXPECT_FALSE(not_after_passed(AttemptIds(), registry.now_ms()));

    // a running id is refused as a duplicate even when the new request has also expired
    auto dup = attempt_of(registry, "a-17");
    ids.attempt_id = "a-17";
    dup->set_ids(ids);
    refused = registry.start(dup);
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_EQ("running", refused->body["state"].asString());
    // a tombstoned id too
    EXPECT_EQ(404, registry.cancel("t-1", 0).first);
    auto tomb = attempt_of(registry, "t-1");
    ids.attempt_id = "t-1";
    tomb->set_ids(ids);
    refused = registry.start(tomb);
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_TRUE(refused->body["tombstone"].asBool());

    // without not_after_ms the usage and the state keep their keys
    auto c = attempt_of(registry, "c");
    EXPECT_FALSE(registry.start(c));
    c->set_bound(1, 1000);
    EXPECT_FALSE(c->usage_json("completed").isMember("not_after_ms"));
    EXPECT_FALSE(registry.state("c").second.isMember("not_after_ms"));
}

// The attempts block of both capabilities routes (W1, W4): integers where usage.bound states
// integers, the hard cap, the clock skew a ledger adds, the not_after rule
TEST(GraphletAttemptRegistry, CapabilitiesStateIntegers) {
    AttemptSettings s;
    s.clock_skew_ms = 2500;
    AttemptRegistry registry(s);
    const Json::Value att = registry.capabilities_json();
    for (const char *f : { "allowance_ms", "hard_cap_ms", "clock_skew_allowance_ms",
                           "retention_s", "retention_count", "content_timeout_s",
                           "client_check_ms" }) {
        EXPECT_EQ(Json::uintValue, att[f].type()) << f;
    }
    EXPECT_EQ(10000u, att["allowance_ms"].asUInt64());
    EXPECT_EQ(899000u, att["hard_cap_ms"].asUInt64());
    EXPECT_EQ(2500u, att["clock_skew_allowance_ms"].asUInt64());
    EXPECT_EQ(900u, att["content_timeout_s"].asUInt64());
    // written as integers (jsoncpp writes an integral double as 10000.0)
    EXPECT_NE(std::string::npos, json_text(att, true).find("\"allowance_ms\":10000,"));
    Json::Value fields(Json::arrayValue);
    for (const char *f : { "attempt_id", "budget_id", "locus_id", "not_after_ms",
                           "expect_server_instance" }) {
        fields.append(f);
    }
    EXPECT_EQ(fields, att["fields"]);
    EXPECT_NE(std::string::npos, att["not_after"].asString().find("clock_skew_allowance_ms"));
    // the review's sentence: an unanswered request may already be running
    EXPECT_NE(std::string::npos, att["not_after"].asString().find("cannot start subsequently"));
    EXPECT_EQ(std::string::npos, att["not_after"].asString().find("never started"));
    EXPECT_EQ(registry.server_instance(), att["server_instance"].asString());
    // the tombstone's cap as applied, the cancel's fields, and the normative release rule
    EXPECT_EQ(Json::uintValue, att["tombstone_max_s"].type());
    EXPECT_EQ(86400u, att["tombstone_max_s"].asUInt64());
    Json::Value cancel_fields(Json::arrayValue);
    for (const char *f : { "attempt_id", "wait_ms", "not_after_ms" }) {
        cancel_fields.append(f);
    }
    EXPECT_EQ(cancel_fields, att["cancel_fields"]);
    for (const char *phrase : { "covers_admission: true", "same server_instance",
                                "exactly that not_after_ms", "expect_server_instance",
                                "never released early",
                                "not_after_ms + clock_skew_allowance_ms + bound_ms",
                                // the review of the pass-5 fixes, finding 1: what a finished
                                // state promises, and what it assumes
                                "A finished attempt's id stays refused (409)",
                                "never dropped early",
                                "assumes that no copy of the request arrives after the attempt "
                                "left retention" }) {
        EXPECT_NE(std::string::npos, att["release_rule"].asString().find(phrase)) << phrase;
    }
    // findings 2 and 3: the hold is live through suppressed_until_ms inclusive, and a refused
    // copy is judged by its own not_after_ms alone
    for (const char *phrase : { "suppressed_until_ms",
                                "the wall clock while it reads at most suppressed_until_ms "
                                "(inclusive)",
                                "for the refused request's 409 its own, absent when it has none",
                                "finished attempts held past their retention" }) {
        EXPECT_NE(std::string::npos, att["suppression"].asString().find(phrase)) << phrase;
    }
    EXPECT_NE(std::string::npos, att["instance"].asString().find("instance_mismatch"));
    // a retention below the cap's setting: the cap is never below it
    s.retention_s = 100'000;
    s.tombstone_max_s = 10;
    EXPECT_EQ(100'000u, AttemptRegistry(s).capabilities_json()["tombstone_max_s"].asUInt64());
    // retention 0: no tombstones, said so
    s.retention_s = 0;
    EXPECT_NE(std::string::npos, AttemptRegistry(s).capabilities_json()["suppression"].asString()
                                         .find("no_suppression"));
}

// Pass 5, W6: the seeds stop being walked at bound - max(allowance / 2, reserve), the reserve
// 1.25 times the time to compress the text written so far and the walked seed's estimated text
// (its account / the account per text byte), and to build the latter — at the configured rates
// and ratios until the server or the attempt measured its own — plus the time from the
// walk-until to the walk's end (review of pass 5, F3: without the margin and that time, a
// response the model fitted exactly was a 503 about half the time)
TEST(GraphletAttempt, DeliveryReserveMovesTheWalkUntil) {
    FakeClock clock;
    AttemptSettings s = settings_with(&clock);
    s.allowance_ms = 10'000;
    s.delivery_compress_mbps = 50;     // 50,000 bytes per ms
    s.delivery_build_mbps = 5;         // 5,000 bytes per ms
    s.delivery_stop_ms = 100;
    // the ratios this arithmetic is written for (the defaults are 30 and 50 since the
    // efficiency pass's calibration)
    s.account_per_text_byte_json = 20;
    s.account_per_text_byte_graphlet = 40;
    auto reserve = [](double model, double stop = 100) { return 1.25 * model + stop; };
    AttemptRegistry registry(s);
    auto a = attempt_of(registry, "reserve");
    a->set_delivery_detail("full");    // 20 bytes of account per text byte
    a->set_bound(2, 30'000);           // 70,000 ms
    EXPECT_EQ(65'000, a->walk_until_ms());
    EXPECT_NEAR(100, a->reserve_ms(), 1e-9);    // nothing to deliver yet: the stop alone
    // small outputs: the floor (allowance / 2)
    a->progress(20 * 1'000'000);       // 1 MB of text estimated
    EXPECT_NEAR(reserve(1e6 / 5e4 + 1e6 / 5e3), a->reserve_ms(), 1e-9);
    EXPECT_EQ(65'000, a->walk_until_ms());
    // a large walked seed: 40 MB of text estimated -> 1.25 x (800 + 8000) + 100 ms
    a->progress(20 * 40'000'000ull);
    EXPECT_NEAR(reserve(800 + 8000), a->reserve_ms(), 1e-9);
    EXPECT_NEAR(70'000 - reserve(8800), a->walk_until_ms(), 1e-9);
    // a graphlet's account counts 40 per text byte
    a->set_delivery_detail("graphlet");
    a->progress(40 * 40'000'000ull);
    EXPECT_NEAR(70'000 - reserve(8800), a->walk_until_ms(), 1e-9);
    // the seed was written: its exact bytes replace the estimate; built at 40 MB/s, from an
    // account of 80 bytes per text byte (both measured: a text of 1 MiB or more)
    a->note_delivered(40'000'000, 1.0, 80 * 40'000'000ull);
    EXPECT_NEAR(reserve(40e6 / 5e4), a->reserve_ms(), 1e-9);          // compression only
    EXPECT_EQ(65'000, a->walk_until_ms());
    EXPECT_EQ(40.0, a->own_build_mbps());
    EXPECT_EQ(80.0, a->own_account_per_text_byte());
    // the next seed's estimate uses the measured ratio (80) and build rate (40 MB/s)
    a->progress(80 * 100'000'000ull);
    EXPECT_NEAR(reserve((40e6 + 100e6) / 5e4 + 100e6 / 4e4), a->reserve_ms(), 1e-6);
    const double until = 70'000 - a->reserve_ms();
    EXPECT_NEAR(until, a->walk_until_ms(), 1e-9);
    // the walk stops at the moved instant; usage states the walk-until in force when it stopped
    // the walk (not the lowest computed, 70,000 - reserve(8800), before the first seed was
    // delivered: review of pass 5, F4), also once the seed is delivered and the reserve shrinks
    ASSERT_FALSE(registry.start(a));
    clock.ms = static_cast<int64_t>(until) - 1;
    EXPECT_EQ(ExternalStop::NONE, a->poll(true));
    clock.ms = static_cast<int64_t>(until) + 1;
    EXPECT_EQ(ExternalStop::ATTEMPT_DEADLINE, a->poll(true));
    // the walk ended 120 ms after its walk-until: what the reserve's stop is measured on
    clock.ms = static_cast<int64_t>(until) + 120;
    a->seed_walked(1, SeedUsage());
    EXPECT_NEAR(120 - (until - std::floor(until)), a->own_stop_ms(), 1e-6);
    a->note_delivered(1000, 0.001);
    EXPECT_EQ(static_cast<uint64_t>(std::ceil(until)),
              a->usage_json("deadline")["bound"]["walk_until_ms"].asUInt64());
    // a slower measured seed lowers the rate used (the slowest wins); a small one measures
    // nothing
    auto b = attempt_of(registry, "reserve2");
    b->set_delivery_detail("full");
    b->set_bound(1, 30'000);
    b->note_delivered(1000, 10.0, 1'000'000);
    EXPECT_EQ(0.0, b->own_build_mbps());
    EXPECT_EQ(0.0, b->own_account_per_text_byte());
    b->progress(20 * 10'000'000ull);
    EXPECT_NEAR(reserve(10e6 / 5e3 + (1000 + 10e6) / 5e4), b->reserve_ms(), 1e-6);
    b->note_delivered(2'000'000, 1.0);         // 2 MB/s, no walk: no ratio
    b->progress(20 * 10'000'000ull);
    EXPECT_NEAR(reserve(10e6 / 2e3 + (2'001'000 + 10e6) / 5e4), b->reserve_ms(), 1e-6);
    // the floor holds in usage while the reserve stays below it (no walk-until stopped it)
    auto f = attempt_of(registry, "reserve-floor");
    f->set_delivery_detail("full");
    f->set_bound(1, 30'000);
    f->progress(1000);
    EXPECT_EQ(35'000u, f->usage_json("completed")["bound"]["walk_until_ms"].asUInt64());
    // ... also when the walk-until falls only once the last seed's text is written (its exact
    // bytes, no walk left): no poll read that one, so it bounded no walk (review of pass 5:
    // usage read 20174 for a walk its own 30 s budget ended)
    clock.ms = 0;
    ASSERT_FALSE(registry.start(f));
    EXPECT_EQ(ExternalStop::NONE, f->poll(true));          // before the seed (between seeds)
    f->seed_started(0);
    f->progress(20 * 1'000'000);
    EXPECT_EQ(ExternalStop::NONE, f->poll(true));
    f->seed_walked(0, SeedUsage());
    f->note_delivered(400'000'000, 10.0, 20 * 400'000'000ull);
    EXPECT_NEAR(40'000 - reserve(400e6 / 5e4), f->walk_until_ms(), 1e-6);
    EXPECT_EQ(35'000u, f->usage_json("completed")["bound"]["walk_until_ms"].asUInt64());
    // the lowest walk-until a poll read counts, the reserve having moved it, though it rose
    // again before the walk ended
    auto w = attempt_of(registry, "reserve-walked");
    w->set_delivery_detail("full");
    w->set_bound(1, 30'000);
    ASSERT_FALSE(registry.start(w));
    w->seed_started(0);
    w->progress(20 * 40'000'000ull);
    EXPECT_EQ(ExternalStop::NONE, w->poll(true));
    w->progress(20 * 1'000ull);
    EXPECT_EQ(ExternalStop::NONE, w->poll(true));
    w->seed_walked(0, SeedUsage());
    EXPECT_EQ(static_cast<uint64_t>(std::ceil(40'000 - reserve(8800))),
              w->usage_json("completed")["bound"]["walk_until_ms"].asUInt64());
    // an attempt whose bound is not enforced (the CLI's) keeps the floor: no delivery is bounded
    auto local = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(),
                                           registry.settings(), registry.server_instance(),
                                           /* enforced */ false);
    AttemptIds local_ids;
    local_ids.attempt_id = "local-1";
    local->set_ids(local_ids);
    local->set_delivery_detail("full");
    local->set_bound(1, 30'000);
    local->progress(20 * 40'000'000ull);
    EXPECT_EQ(35'000, local->walk_until_ms());
    // the capabilities state the rule and its parts
    Json::Value r = registry.capabilities_json()["delivery_reserve"];
    EXPECT_EQ(50, r["compress_mbps"].asDouble());
    EXPECT_EQ(5, r["build_mbps"].asDouble());
    EXPECT_EQ(20, r["account_per_text_byte"]["json"].asDouble());
    EXPECT_EQ(40, r["account_per_text_byte"]["graphlet"].asDouble());
    EXPECT_NE(std::string::npos, r["calibration"].asString().find("starting estimates"));
    // the calibrated starting estimates (the efficiency pass; 20, 40 and 250 before)
    const AttemptSettings defaults;
    EXPECT_EQ(30, defaults.account_per_text_byte_json);
    EXPECT_EQ(50, defaults.account_per_text_byte_graphlet);
    EXPECT_EQ(1000, defaults.delivery_stop_ms);
    EXPECT_EQ(1.25, r["margin"].asDouble());
    EXPECT_TRUE(r["stop_ms"].isIntegral());
    EXPECT_EQ(100u, r["stop_ms"].asUInt64());
    EXPECT_TRUE(r["measured_stop_ms"].isNull());
    EXPECT_TRUE(r["measured_build_mbps"].isNull());
    EXPECT_TRUE(r["measured_compress_mbps"].isNull());
    EXPECT_TRUE(r["measured_account_per_text_byte"]["full"].isNull());
    EXPECT_NE(std::string::npos, registry.capabilities_json()["bound"].asString()
                                         .find("the delivery reserve"));

    // the server's measurements replace the configured values: the slowest rate and the
    // smallest ratio of the last 16, the longest stop latency
    for (int i = 0; i < 20; ++i) {
        registry.note_build_rate(i == 2 ? 1.0 : 100.0 + i);     // the slow one ages out
        registry.note_compress_rate(400.0 + i);
        registry.note_account_per_text_byte("full", 60.0 + i);
        registry.note_stop_latency(i == 1 ? 5000.0 : 40.0 + i); // the long one ages out
    }
    registry.note_build_rate(0);                                 // nothing measured: ignored
    registry.note_account_per_text_byte("graphlet", 0);
    registry.note_stop_latency(0);
    const DeliveryMeasurements m = registry.measured();
    EXPECT_EQ(104.0, m.build_mbps);
    EXPECT_EQ(404.0, m.compress_mbps);
    EXPECT_EQ(64.0, m.account_per_text_byte.at("full"));
    EXPECT_EQ(0u, m.account_per_text_byte.count("graphlet"));
    EXPECT_EQ(59.0, m.stop_ms);
    r = registry.capabilities_json()["delivery_reserve"];
    EXPECT_EQ(104.0, r["measured_build_mbps"].asDouble());
    EXPECT_EQ(404.0, r["measured_compress_mbps"].asDouble());
    EXPECT_EQ(64.0, r["measured_account_per_text_byte"]["full"].asDouble());
    EXPECT_TRUE(r["measured_account_per_text_byte"]["graphlet"].isNull());
    EXPECT_EQ(59u, r["measured_stop_ms"].asUInt64());
    auto c = attempt_of(registry, "reserve3");
    c->set_measured(m);
    c->set_delivery_detail("full");
    c->set_bound(1, 30'000);
    c->progress(64 * 40'000'000ull);                             // 40 MB of text estimated
    // a measured stop shorter than the configured one does not lower it
    EXPECT_NEAR(reserve(40e6 / 404e3 + 40e6 / 104e3), c->reserve_ms(), 1e-6);
    EXPECT_EQ(35'000, c->walk_until_ms());                       // within the floor
    // the attempt's own slower seed, and smaller ratio, win over the server's
    c->note_delivered(2'000'000, 1.0, 30 * 2'000'000ull);        // 2 MB/s, 30 per byte
    c->progress(30 * 40'000'000ull);
    EXPECT_NEAR(reserve((2e6 + 40e6) / 404e3 + 40e6 / 2e3), c->reserve_ms(), 1e-6);
    // a detail the server has not measured keeps its configured ratio; a longer measured stop
    // replaces the configured one
    DeliveryMeasurements slow_stop = m;
    slow_stop.stop_ms = 300;
    auto g = attempt_of(registry, "reserve4");
    g->set_measured(slow_stop);
    g->set_delivery_detail("graphlet");
    g->set_bound(1, 30'000);
    g->progress(40 * 1'000'000ull);
    EXPECT_NEAR(reserve(1e6 / 404e3 + 1e6 / 104e3, 300), g->reserve_ms(), 1e-6);
}

// usage.bound.walk_until_ms, when no walk-until stopped the walk, is the lowest walk-until that
// a poll READING THE CLOCK compared with: every poll_stride-th poll of a walk, every forced poll
// (between seeds, before a paced read's chunk, in the lookahead). A lower one that the reserve
// set at a level's end and that only polls not reading the clock saw is not stated — the rule
// "a value no poll read bounded no walk" (review of pass 5, F4). The review of 2026-10-06 (C20,
// U11-02) found the texts claiming "the lowest in force while a seed was walked"; this pins the
// code's rule with the server's stride (8) as the corrected texts state it. The reviewer's
// driver: one head per level and one non-forced poll before each; at the end of level 8 the
// seed's estimated text jumps from 1 MB to 40 MB (walk-until 65,000 -> 58,900); the clock-reading
// poll was the 8th (level 7), so levels 9-11 are walked under 58,900 with polls that do not read
// the clock, and the delivered seed raises the walk-until again: stride 8 states 65,000, stride 1
// (every poll reads the clock) 58,900
TEST(GraphletAttempt, WalkUntilStatedIsTheLowestAClockReadingPollSaw) {
    for (uint32_t stride : { 8u, 1u }) {
        FakeClock clock;
        AttemptSettings s = settings_with(&clock);
        s.allowance_ms = 10'000;
        s.poll_stride = stride;
        s.delivery_compress_mbps = 50;
        s.delivery_build_mbps = 5;
        s.delivery_stop_ms = 100;
        s.account_per_text_byte_json = 20;
        s.account_per_text_byte_graphlet = 40;
        AttemptRegistry registry(s);
        auto a = attempt_of(registry, "u11-02-" + std::to_string(stride));
        a->set_delivery_detail("full");
        a->set_bound(2, 30'000);           // bound 70,000 ms, floor 65,000
        ASSERT_FALSE(registry.start(a));
        double lowest_in_force = std::numeric_limits<double>::infinity();
        // seed 0: one head per level, its poll before it
        EXPECT_EQ(ExternalStop::NONE, a->poll(true));
        a->seed_started(0);
        for (int level = 0; level < 12; ++level) {
            clock.ms += 1000;
            EXPECT_EQ(ExternalStop::NONE, a->poll());
            lowest_in_force = std::min(lowest_in_force, a->walk_until_ms());
            a->progress(level >= 8 ? 20 * 40'000'000ull : 20 * 1'000'000ull);
            if (level + 1 < 12)
                lowest_in_force = std::min(lowest_in_force, a->walk_until_ms());
        }
        EXPECT_NEAR(58'900, lowest_in_force, 1e-6);
        a->seed_walked(0, SeedUsage());
        a->note_delivered(40'000'000, 1.0, 20 * 40'000'000ull);
        EXPECT_EQ(65'000, a->walk_until_ms());
        // seed 1: a small walk
        clock.ms += 10;
        EXPECT_EQ(ExternalStop::NONE, a->poll(true));
        a->seed_started(1);
        for (int level = 0; level < 3; ++level) {
            clock.ms += 1000;
            EXPECT_EQ(ExternalStop::NONE, a->poll());
            a->progress(20 * 100'000ull);
        }
        a->seed_walked(1, SeedUsage());
        a->note_delivered(100'000, 0.01);
        a->walk_ended();
        EXPECT_EQ(stride == 8 ? 65'000u : 58'900u,
                  a->usage_json("completed")["bound"]["walk_until_ms"].asUInt64())
            << "poll_stride " << stride;
    }
    // the texts state the clock-reading polls, not "every poll" (C20)
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    const std::string bound = registry.capabilities_json()["bound"].asString();
    EXPECT_NE(std::string::npos, bound.find("its first poll that reads the clock after it"))
        << bound;
    EXPECT_NE(std::string::npos, bound.find("not a lower one in force only between two such "
                                            "polls")) << bound;
}

// The delivery reserve's coordinate share (C3, plan revision 3): the walked seed's record
// coordinates (their part of its account) are estimated at kCoordinateAccountPerTextByte, the
// rest at the ratio in use; without coordinates the estimate is the one before the split
TEST(GraphletAttempt, ReserveCountsCoordinateText) {
    FakeClock clock;
    AttemptSettings s = settings_with(&clock);
    s.allowance_ms = 10'000;
    s.delivery_compress_mbps = 50;     // 50,000 bytes per ms
    s.delivery_build_mbps = 5;         // 5,000 bytes per ms
    s.delivery_stop_ms = 100;
    s.account_per_text_byte_json = 20;
    auto reserve = [](double text) { return 1.25 * (text / 5e4 + text / 5e3) + 100; };
    AttemptRegistry registry(s);
    auto a = attempt_of(registry, "coordinates");
    a->set_delivery_detail("full");
    a->set_bound(2, 30'000);
    // 20 account bytes a text byte for the rest, 12 for the coordinates: 1e6 + 1e5 bytes of text
    a->progress(20 * 1'000'000ull + 12 * 100'000ull, 12 * 100'000ull);
    EXPECT_NEAR(reserve(1'100'000), a->reserve_ms(), 1e-9);
    // the same account without a coordinate share: all of it at 20 (the formula before C3)
    a->progress(20 * 1'000'000ull + 12 * 100'000ull);
    EXPECT_NEAR(reserve(1'060'000), a->reserve_ms(), 1e-9);
    // both rounded up apart (a byte of text never priced below its account)
    a->progress(41, 13);
    EXPECT_NEAR(reserve(2 + 2), a->reserve_ms(), 1e-9);
    // a share larger than the account (never reported) is taken as the whole account
    a->progress(24, 1000);
    EXPECT_NEAR(reserve(2), a->reserve_ms(), 1e-9);
    // the ratio sample leaves the coordinate share's account and its exact text out: 2 MB of
    // text of which 0.5 MB coordinates, from 80 x 1.5 MB + 6 MB
    a->note_delivered(2'000'000, 1.0, 80 * 1'500'000ull + 6'000'000, 6'000'000, 500'000);
    EXPECT_EQ(80.0, a->own_account_per_text_byte());
    EXPECT_EQ(2.0, a->own_build_mbps());               // the build rate over the whole text
    // a rest below measured_text_bytes measures no ratio, even in a large text
    auto b = attempt_of(registry, "coordinates-heavy");
    b->set_delivery_detail("full");
    b->set_bound(1, 30'000);
    b->note_delivered(2'000'000, 1.0, 30 * 1'000'000ull + 12 * 1'500'000ull, 12 * 1'500'000ull,
                      1'500'000);
    EXPECT_EQ(0.0, b->own_account_per_text_byte());
    EXPECT_EQ(2.0, b->own_build_mbps());
    // stated in the capabilities, a number a ledger can compute with
    const Json::Value r = registry.capabilities_json()["delivery_reserve"];
    ASSERT_TRUE(r["coordinate_account_per_text_byte"].isUInt64());
    EXPECT_EQ(12u, r["coordinate_account_per_text_byte"].asUInt64());
    EXPECT_EQ(kCoordinateAccountPerTextByte, r["coordinate_account_per_text_byte"].asUInt64());
    EXPECT_NE(std::string::npos, r["rule"].asString().find("ceil(C / coordinate_account_per_text_byte)"));
}

// Plan revision 3 (review of W1, finding 3): an attempt with record coordinates feeds the
// server's measured ratio with the sample the same walk without them gives (their account and
// their exact text left out), so a later attempt without coordinates reads the same
// measured_account_per_text_byte, and walks to the same walk-until, whether or not one with
// coordinates came first. The exactness of the two parts on real responses is
// MiniRefSeq.CoordinateShareIsExact's
TEST(GraphletAttempt, CoordinatesLeaveTheServersRatioUnchanged) {
    struct Seen {
        double measured;
        double walk_until;
    };
    auto run = [](bool with_coordinates) -> Seen {
        FakeClock clock;
        AttemptRegistry registry(settings_with(&clock));
        auto deliver = [&](const std::string &id, uint64_t coordinate_account,
                           uint64_t coordinate_text, double ratio) {
            auto a = attempt_of(registry, id);
            a->set_measured(registry.measured());
            a->set_delivery_detail("full");
            a->set_bound(1, 30'000);
            const uint64_t text = 2'000'000;
            a->note_delivered(text + coordinate_text, 1.0,
                              static_cast<uint64_t>(ratio * text) + coordinate_account,
                              coordinate_account, coordinate_text);
            EXPECT_EQ(ratio, a->own_account_per_text_byte()) << id;
            // as the server does once the response is written (server.cpp, on_written)
            registry.note_account_per_text_byte(a->delivery_detail(),
                                                a->own_account_per_text_byte());
        };
        deliver("first", 0, 0, 100);
        // the same walk with coordinates (about 20 account bytes a byte of their text), or
        // without them
        if (with_coordinates) {
            deliver("second", 20 * 300'000, 300'000, 100);
        } else {
            deliver("second", 0, 0, 100);
        }
        const DeliveryMeasurements m = registry.measured();
        auto later = attempt_of(registry, "later");
        later->set_measured(m);
        later->set_delivery_detail("full");
        later->set_bound(1, 30'000);
        later->progress(100 * 40'000'000ull);
        return { m.account_per_text_byte.at("full"), later->walk_until_ms() };
    };
    const Seen without = run(false), with = run(true);
    EXPECT_EQ(100.0, without.measured);
    EXPECT_EQ(without.measured, with.measured);
    EXPECT_EQ(without.walk_until, with.walk_until);
}

TEST(GraphletAttemptRegistry, AnIdRunsOnce) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock, 60));
    EXPECT_EQ(16u, registry.server_instance().size());
    auto a = attempt_of(registry, "x");
    EXPECT_FALSE(registry.start(a));
    // running: refused, with the running attempt's state
    auto b = attempt_of(registry, "x");
    auto conflict = registry.start(b);
    ASSERT_TRUE(conflict);
    EXPECT_EQ("running", conflict->body["state"].asString());
    registry.finish(a, "completed", 200, 123);
    // retained: refused, the finished state
    conflict = registry.start(b);
    ASSERT_TRUE(conflict);
    EXPECT_EQ("finished", conflict->body["state"].asString());
    EXPECT_EQ("completed", conflict->body["reason"].asString());
    // expired: the id is free again
    clock.ms = 60'000;
    EXPECT_FALSE(registry.start(b));
    EXPECT_EQ("running", registry.state("x").second["state"].asString());
}

// A cancel is acknowledged (200, stopping) and repeatable; the walk sees it at its next poll;
// once the attempt finished a cancel finds nothing to stop (404) and GET states when it stopped
TEST(GraphletAttemptRegistry, CancelIsIdempotentAndAcknowledged) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    auto a = attempt_of(registry, "c");
    a->set_bound(2, 1000);
    ASSERT_FALSE(registry.start(a));
    a->seed_started(0);
    clock.ms = 10;
    auto [status, body] = registry.cancel("c", 0);
    EXPECT_EQ(200, status);
    EXPECT_TRUE(body["cancelled"].asBool());
    EXPECT_EQ("stopping", body["state"].asString());
    EXPECT_EQ("cancel", body["attempt"]["stop_requested_by"].asString());
    EXPECT_FALSE(body["attempt"]["cancel_requested_at"].isNull());
    // the walk has not seen it yet
    EXPECT_TRUE(body["attempt"]["stopped_at"].isNull());
    EXPECT_EQ(200, registry.cancel("c", 0).first);
    clock.ms = 25;
    EXPECT_EQ(ExternalStop::CANCELLED, a->poll());
    EXPECT_STREQ("cancelled", a->walk_reason());
    a->walk_ended();
    clock.ms = 40;
    registry.finish(a, a->walk_reason(), 200, 4567);
    auto [after, finished] = registry.cancel("c", 0);
    EXPECT_EQ(404, after);
    EXPECT_FALSE(finished["cancelled"].asBool());
    EXPECT_EQ("finished", finished["state"].asString());
    const auto [ok, state] = registry.state("c");
    EXPECT_EQ(200, ok);
    EXPECT_EQ("finished", state["state"].asString());
    EXPECT_EQ("cancelled", state["reason"].asString());
    EXPECT_TRUE(state["response"]["written"].asBool());
    EXPECT_EQ(4567u, state["response"]["bytes"].asUInt64());
    EXPECT_EQ(40u, state["elapsed_ms"].asUInt64());
    // the walk stopped at the poll that saw the cancel: 15 ms after it was asked
    EXPECT_EQ("cancelled", state["usage"]["reason"].asString());
    EXPECT_FALSE(state["usage"].isMember("per_seed"));
    EXPECT_NE(state["stopped_at"].asString(), state["cancel_requested_at"].asString());
}

// A cancel that overtakes its request tombstones the id: 404 (nothing ran), and a request that
// arrives later with it is refused, so the 404 is an acknowledgement a ledger can release on
TEST(GraphletAttemptRegistry, CancelOfAnUnknownIdRefusesItLater) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock, 30));
    auto [status, body] = registry.cancel("early", 0);
    EXPECT_EQ(404, status);
    EXPECT_EQ("unknown", body["state"].asString());
    EXPECT_TRUE(body["tombstone"].asBool());
    EXPECT_EQ(registry.server_instance(), body["server_instance"].asString());
    EXPECT_EQ(404, registry.cancel("early", 0).first);
    EXPECT_EQ(404, registry.state("early").first);
    auto late = attempt_of(registry, "early");
    auto conflict = registry.start(late);
    ASSERT_TRUE(conflict);
    EXPECT_TRUE(conflict->body["tombstone"].asBool());
    // an id never named is simply unknown
    auto [unknown, j] = registry.state("never");
    EXPECT_EQ(404, unknown);
    EXPECT_FALSE(j.isMember("tombstone"));
    // held through retention_s inclusive (the wall clock's expiry; GET does not extend it, a
    // refused copy would), gone after it
    clock.ms = 30'000;
    EXPECT_TRUE(registry.state("early").second.isMember("tombstone"));
    clock.ms = 30'001;
    EXPECT_FALSE(registry.start(late));
}

// The review of the stage-4 backend, F6: a tombstone is kept its whole retention period,
// whatever finishes after it (they shared retention_count, so later finishes evicted it and the
// cancelled id ran); at most retention_count tombstones are held, and a cancel of an unknown id
// beyond that is refused (429, tombstone: false) rather than promised for less
TEST(GraphletAttemptRegistry, TombstonesOutliveLaterFinishes) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock, 30, 3));
    auto [status, body] = registry.cancel("Y", 0);
    ASSERT_EQ(404, status);
    EXPECT_TRUE(body["tombstone"].asBool());
    for (int i = 0; i < 10; ++i) {
        auto a = attempt_of(registry, "done" + std::to_string(i));
        ASSERT_FALSE(registry.start(a));
        clock.ms += 100;
        registry.finish(a, "completed", 200, 1);
    }
    // ten finishes later, within retention_s: Y is still refused
    auto late = attempt_of(registry, "Y");
    auto conflict = registry.start(late);
    ASSERT_TRUE(conflict);
    EXPECT_TRUE(conflict->body["tombstone"].asBool());
    // the finished attempts are kept by count as before
    EXPECT_EQ(404, registry.state("done0").first);
    EXPECT_EQ(200, registry.state("done9").first);
    // two more tombstones fill the room; a fourth is refused, nothing promised
    EXPECT_EQ(404, registry.cancel("Z1", 0).first);
    EXPECT_EQ(404, registry.cancel("Z2", 0).first);
    auto [full, refused] = registry.cancel("Z3", 0);
    EXPECT_EQ(429, full);
    EXPECT_FALSE(refused["tombstone"].asBool());
    EXPECT_FALSE(refused["cancelled"].asBool());
    EXPECT_EQ("unknown", refused["state"].asString());
    EXPECT_NE(std::string::npos, refused["error"].asString().find("NOT tombstoned"));
    // so a request with that id runs
    auto z3 = attempt_of(registry, "Z3");
    EXPECT_FALSE(registry.start(z3));
    registry.finish(z3, "completed", 200, 1);
    // a known tombstone answers as before, without taking more room
    EXPECT_EQ(404, registry.cancel("Z1", 0).first);
    // by age only: after retention_s (inclusive) they go, and the ids are free
    clock.ms += 30'001;
    EXPECT_FALSE(registry.start(late));
    EXPECT_EQ(404, registry.cancel("Z4", 0).first);
}

// Review of pass 5, finding 1: a tombstone held for retention_s alone let a half-uploaded
// request whose not_after_ms lay beyond it run once it expired. A cancel naming the request's
// not_after_ms holds the tombstone until not_after_ms + clock_skew_ms (within the cap), every
// tombstone answer states suppressed_until_ms and covers_admission, and the timeline the review
// walked ends with the late copy refused: by the tombstone while it is live, by the strict
// not_after check once it is gone
TEST(GraphletAttemptRegistry, CancelNotAfterMsHoldsTheTombstoneThroughAdmission) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 1);   // retention 1 s, as the review's
    settings.clock_skew_ms = 2000;
    AttemptRegistry registry(settings);
    const uint64_t t0 = clock.wall_ms();
    const uint64_t not_after = t0 + 60'000;

    auto [status, body] = registry.cancel("delayed-original", 0, not_after);
    ASSERT_EQ(404, status);
    EXPECT_TRUE(body["tombstone"].asBool());
    EXPECT_FALSE(body["cancelled"].asBool());
    EXPECT_EQ(Json::uintValue, body["suppressed_until_ms"].type());
    EXPECT_EQ(not_after + 2000, body["suppressed_until_ms"].asUInt64());
    EXPECT_EQ(not_after, body["not_after_ms"].asUInt64());
    EXPECT_TRUE(body["covers_admission"].asBool());
    EXPECT_FALSE(body.isMember("covers_admission_reason"));
    EXPECT_EQ(registry.server_instance(), body["server_instance"].asString());

    auto copy = [&](const std::string &id, std::optional<uint64_t> na) {
        auto a = attempt_of(registry, id);
        AttemptIds ids = a->ids();
        ids.not_after_ms = na;
        a->set_ids(ids);
        return a;
    };
    // the repeated cancel at 0.763 s states the same expiry
    clock.ms = 763;
    auto [again, repeated] = registry.cancel("delayed-original", 0);
    EXPECT_EQ(404, again);
    EXPECT_EQ(not_after + 2000, repeated["suppressed_until_ms"].asUInt64());
    EXPECT_EQ(not_after, repeated["not_after_ms"].asUInt64());
    EXPECT_TRUE(repeated["covers_admission"].asBool());
    // the original upload completes at 1.177 s, past retention_s: refused, nothing registered
    clock.ms = 1177;
    auto refused = registry.start(copy("delayed-original", not_after));
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_TRUE(refused->body["tombstone"].asBool());
    EXPECT_TRUE(refused->body["covers_admission"].asBool());
    EXPECT_EQ(not_after + 2000, refused->body["suppressed_until_ms"].asUInt64());
    // GET states it too (read only)
    auto [get_status, got] = registry.state("delayed-original");
    EXPECT_EQ(404, get_status);
    EXPECT_TRUE(got["tombstone"].asBool());
    EXPECT_TRUE(got["covers_admission"].asBool());
    EXPECT_EQ(not_after + 2000, got["suppressed_until_ms"].asUInt64());
    // through the expiry still held (GET, which does not extend it: the expiry is inclusive);
    // after it the tombstone is gone, and a copy is refused by its own not_after_ms (409
    // expired), which passed 2 s ago
    auto tombstoned = [&](const std::string &id) {
        return registry.state(id).second.isMember("tombstone");
    };
    clock.ms = static_cast<int64_t>(not_after + 2000 - t0) - 1;
    EXPECT_TRUE(tombstoned("delayed-original"));
    clock.ms = static_cast<int64_t>(not_after + 2000 - t0);
    EXPECT_TRUE(tombstoned("delayed-original"));
    clock.ms = static_cast<int64_t>(not_after + 2000 - t0) + 1;
    EXPECT_FALSE(tombstoned("delayed-original"));
    refused = registry.start(copy("delayed-original", not_after));
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->expired);
    EXPECT_EQ(404, registry.state("delayed-original").first);
    EXPECT_FALSE(registry.state("delayed-original").second.isMember("tombstone"));

    // Without not_after_ms the tombstone is held retention_s and covers nothing: a ledger may
    // not release on it, and the review's late upload runs once it expired
    const uint64_t t1 = clock.wall_ms();
    auto [plain_status, plain] = registry.cancel("no-not-after", 0);
    EXPECT_EQ(404, plain_status);
    EXPECT_EQ(t1 + 1000, plain["suppressed_until_ms"].asUInt64());
    EXPECT_FALSE(plain["covers_admission"].asBool());
    EXPECT_EQ("no_not_after_ms", plain["covers_admission_reason"].asString());
    EXPECT_FALSE(plain.isMember("not_after_ms"));
    clock.ms += 1001;
    EXPECT_FALSE(registry.start(copy("no-not-after", std::nullopt)));
}

// A repeated cancel never shortens a tombstone, extends it with a later not_after_ms, and keeps
// it at least retention_s from then; the cap bounds the extension (and covers_admission says
// when the cap left the not_after_ms uncovered)
TEST(GraphletAttemptRegistry, TombstonesAreNeverShortenedAndCapped) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 10);
    settings.clock_skew_ms = 500;
    settings.tombstone_max_s = 100;
    AttemptRegistry registry(settings);
    const uint64_t t0 = clock.wall_ms();
    auto [s1, b1] = registry.cancel("x", 0, t0 + 50'000);
    ASSERT_EQ(404, s1);
    EXPECT_EQ(t0 + 50'500, b1["suppressed_until_ms"].asUInt64());
    // an earlier not_after_ms: unchanged (never shortened), and judged against the cancel's own
    auto [s2, b2] = registry.cancel("x", 0, t0 + 1'000);
    EXPECT_EQ(t0 + 50'500, b2["suppressed_until_ms"].asUInt64());
    EXPECT_EQ(t0 + 1'000, b2["not_after_ms"].asUInt64());
    EXPECT_TRUE(b2["covers_admission"].asBool());
    // a later one extends it
    auto [s3, b3] = registry.cancel("x", 0, t0 + 70'000);
    EXPECT_EQ(t0 + 70'500, b3["suppressed_until_ms"].asUInt64());
    // without one: judged against the largest named, and kept at least retention_s from now
    clock.ms = 65'000;
    auto [s4, b4] = registry.cancel("x", 0);
    EXPECT_EQ(t0 + 75'000, b4["suppressed_until_ms"].asUInt64());
    EXPECT_EQ(t0 + 70'000, b4["not_after_ms"].asUInt64());
    EXPECT_TRUE(b4["covers_admission"].asBool());
    // beyond the cap (100 s from now): held to the cap, not covered
    auto [s5, b5] = registry.cancel("y", 0, clock.wall_ms() + 500'000);
    EXPECT_EQ(404, s5);
    EXPECT_EQ(clock.wall_ms() + 100'000, b5["suppressed_until_ms"].asUInt64());
    EXPECT_FALSE(b5["covers_admission"].asBool());
    EXPECT_EQ("beyond_tombstone_max", b5["covers_admission_reason"].asString());
    // a not_after_ms already passed is covered by retention_s alone
    auto [s6, b6] = registry.cancel("z", 0, clock.wall_ms() - 5'000);
    EXPECT_EQ(clock.wall_ms() + 10'000, b6["suppressed_until_ms"].asUInt64());
    EXPECT_TRUE(b6["covers_admission"].asBool());
}

// A copy of a cancelled request refused at dispatch (409) extends the tombstone to its own
// not_after_ms + clock_skew_ms and states the suppression against it, so that every later copy
// (a proxy's replay carries the same not_after_ms) is refused while it could still be admitted
TEST(GraphletAttemptRegistry, ARefusedCopyExtendsTheTombstone) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 1);
    settings.clock_skew_ms = 2000;
    AttemptRegistry registry(settings);
    ASSERT_EQ(404, registry.cancel("r", 0).first);   // no not_after_ms: 1 s
    const uint64_t not_after = clock.wall_ms() + 30'000;
    auto copy = [&]() {
        auto a = attempt_of(registry, "r");
        AttemptIds ids = a->ids();
        ids.not_after_ms = not_after;
        a->set_ids(ids);
        return a;
    };
    clock.ms = 500;
    auto refused = registry.start(copy());
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_EQ("unknown", refused->body["state"].asString());
    EXPECT_TRUE(refused->body["tombstone"].asBool());
    EXPECT_EQ(not_after + 2000, refused->body["suppressed_until_ms"].asUInt64());
    EXPECT_EQ(not_after, refused->body["not_after_ms"].asUInt64());
    EXPECT_TRUE(refused->body["covers_admission"].asBool());
    // past the 1 s it was cancelled for: a later copy is still refused
    clock.ms = 20'000;
    refused = registry.start(copy());
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_TRUE(refused->body["covers_admission"].asBool());
    // and once the tombstone is gone, by its own not_after_ms (each refusal also kept the
    // tombstone at least retention_s from then, which the last one, at 20 s, did not lengthen)
    clock.ms = 32'000;
    EXPECT_TRUE(registry.state("r").second.isMember("tombstone"));
    clock.ms = 32'001;
    EXPECT_FALSE(registry.state("r").second.isMember("tombstone"));
    EXPECT_TRUE(registry.start(copy())->expired);
}

// The hold is kept while either clock says it is live: a forward step of the wall clock does
// not shorten it (the steady clock holds it), a backward one lengthens it (the wall clock must
// read suppressed_until_ms too)
TEST(GraphletAttemptRegistry, TombstonesSurviveWallClockSteps) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 10);
    AttemptRegistry registry(settings);
    // GET reads the hold without extending it (a refused copy would extend it)
    auto tombstoned = [&](const std::string &id) {
        return registry.state(id).second.isMember("tombstone");
    };
    ASSERT_EQ(404, registry.cancel("fwd", 0).first);
    // the wall clock jumps an hour ahead: the steady clock still holds it for 10 s (through
    // its last ms, as the wall clock would)
    clock.wall_step = 3'600'000;
    clock.ms = 10'000;
    EXPECT_TRUE(tombstoned("fwd"));
    clock.ms = 10'001;
    EXPECT_FALSE(tombstoned("fwd"));
    // a minute back: the steady hold of a new tombstone ends 10 s later, the wall clock reads
    // past its expiry only 60 s after that, and until then the tombstone stays
    clock.wall_step = 0;
    ASSERT_EQ(404, registry.cancel("back", 0).first);   // at 10'001: held through 20'001
    clock.wall_step = -60'000;
    clock.ms = 20'002;
    EXPECT_TRUE(tombstoned("back"));
    clock.ms = 80'001;
    EXPECT_TRUE(tombstoned("back"));
    clock.ms = 80'002;
    EXPECT_FALSE(tombstoned("back"));
    EXPECT_FALSE(registry.start(attempt_of(registry, "back")));
}

// Review of pass 5, finding 4: retention 0 means no suppression — a cancel of an unknown id is
// refused (429, tombstone: false, reason no_suppression), never promised and expired at once
TEST(GraphletAttemptRegistry, RetentionZeroKeepsNoTombstones) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock, 0));
    auto [status, body] = registry.cancel("c0", 0, clock.wall_ms() + 60'000);
    EXPECT_EQ(429, status);
    EXPECT_FALSE(body["tombstone"].asBool());
    EXPECT_FALSE(body["cancelled"].asBool());
    EXPECT_EQ("no_suppression", body["reason"].asString());
    EXPECT_FALSE(body.isMember("suppressed_until_ms"));
    EXPECT_NE(std::string::npos, body["error"].asString().find("NOT tombstoned"));
    // so a request with the id runs
    EXPECT_FALSE(registry.start(attempt_of(registry, "c0")));
    EXPECT_NE(std::string::npos, registry.retention_text().find("not kept"));
}

// The restart hole of finding 1: a request naming another server_instance is refused before
// anything (409 instance_mismatch, nothing registered, ahead of a duplicate or an expired
// not_after_ms); naming this one, it runs
TEST(GraphletAttemptRegistry, ExpectServerInstanceRefusesAnotherProcess) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock));
    auto with = [&](const std::string &id, const std::string &instance,
                    std::optional<uint64_t> not_after = std::nullopt) {
        auto a = attempt_of(registry, id);
        AttemptIds ids = a->ids();
        ids.budget_id = "b";
        ids.expect_server_instance = instance;
        ids.not_after_ms = not_after;
        a->set_ids(ids);
        return a;
    };
    const std::string other = registry.server_instance() == "0123456789abcdef"
                            ? "fedcba9876543210" : "0123456789abcdef";
    auto refused = registry.start(with("i1", other, clock.wall_ms() - 1));
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->instance_mismatch);
    EXPECT_FALSE(refused->expired);
    const Json::Value &body = refused->body;
    EXPECT_EQ((std::vector<std::string>{ "attempt_id", "budget_id", "error",
                                         "expect_server_instance", "not_after_ms",
                                         "server_instance", "state" }), body.getMemberNames());
    EXPECT_EQ("instance_mismatch", body["state"].asString());
    EXPECT_EQ(other, body["expect_server_instance"].asString());
    EXPECT_EQ(registry.server_instance(), body["server_instance"].asString());
    EXPECT_EQ(404, registry.state("i1").first);
    // ahead of a tombstone: another process's id means nothing here
    ASSERT_EQ(404, registry.cancel("i2", 0).first);
    EXPECT_TRUE(registry.start(with("i2", other))->instance_mismatch);
    // this process's instance: as without the field
    EXPECT_FALSE(registry.start(with("i3", registry.server_instance())));
    EXPECT_EQ("running", registry.state("i3").second["state"].asString());
    EXPECT_TRUE(instance_mismatch(with("x", other)->ids(), registry.server_instance()));
    EXPECT_FALSE(instance_mismatch(attempt_of(registry, "y")->ids(), registry.server_instance()));
}

// Review of the pass-5 fixes, finding 1: a ledger releases on a finished state, so a finished
// attempt's id must stay refused while a replay of the request could still be admitted — until
// its not_after_ms + clock_skew_ms (within the cap from its finish), past retention_s and past
// retention_count, among the tombstones (a cancel of an unknown id is refused while they fill
// the table). Without not_after_ms nothing bounds a replay: the id goes with its retention, as
// stated, and with retention 0 nothing is kept at all
TEST(GraphletAttemptRegistry, FinishedAttemptsAreHeldThroughTheirNotAfterMs) {
    auto sent = [](AttemptRegistry &registry, const std::string &id,
                   std::optional<uint64_t> not_after) {
        auto a = attempt_of(registry, id);
        AttemptIds ids = a->ids();
        ids.not_after_ms = not_after;
        ids.expect_server_instance = registry.server_instance();
        a->set_ids(ids);
        return a;
    };
    {
        // by age (the reviewer's p1: retention 1 s, the replay at 1.3 s ran)
        FakeClock clock;
        AttemptSettings settings = settings_with(&clock, 1);
        settings.clock_skew_ms = 2000;
        AttemptRegistry registry(settings);
        const uint64_t not_after = clock.wall_ms() + 60'000;
        auto original = sent(registry, "fin-1", not_after);
        ASSERT_FALSE(registry.start(original));
        registry.finish(original, "completed", 200, 1);
        auto [cs, cancelled] = registry.cancel("fin-1", 0, not_after);
        EXPECT_EQ(404, cs);
        EXPECT_EQ("finished", cancelled["state"].asString());
        EXPECT_EQ("finished", registry.state("fin-1").second["state"].asString());
        clock.ms = 1'300;
        auto replay = registry.start(sent(registry, "fin-1", not_after));
        ASSERT_TRUE(replay);
        EXPECT_FALSE(replay->expired);
        EXPECT_EQ("finished", replay->body["state"].asString());
        EXPECT_EQ(200, registry.state("fin-1").first);
        // through not_after_ms + the skew, inclusive; after it a replay is expired
        clock.ms = 62'000;
        replay = registry.start(sent(registry, "fin-1", not_after));
        ASSERT_TRUE(replay);
        EXPECT_FALSE(replay->expired);
        clock.ms = 62'001;
        EXPECT_EQ(404, registry.state("fin-1").first);
        replay = registry.start(sent(registry, "fin-1", not_after));
        ASSERT_TRUE(replay);
        EXPECT_TRUE(replay->expired);
    }
    {
        // by count (the reviewer's p1b: retention_count 1, one later finish evicted it)
        FakeClock clock;
        AttemptRegistry registry(settings_with(&clock, 3600, 1));
        const uint64_t not_after = clock.wall_ms() + 60'000;
        auto original = sent(registry, "ev-1", not_after);
        ASSERT_FALSE(registry.start(original));
        registry.finish(original, "completed", 200, 1);
        auto other = sent(registry, "ev-2", std::nullopt);
        ASSERT_FALSE(registry.start(other));
        registry.finish(other, "completed", 200, 1);
        auto replay = registry.start(sent(registry, "ev-1", not_after));
        ASSERT_TRUE(replay);
        EXPECT_FALSE(replay->expired);
        EXPECT_EQ("finished", replay->body["state"].asString());
        // held among the tombstones: the table (retention_count 1) is full, so a cancel of an
        // unknown id promises nothing
        auto [full, refused] = registry.cancel("unknown-1", 0);
        EXPECT_EQ(429, full);
        EXPECT_EQ("tombstones_full", refused["reason"].asString());
        // one without not_after_ms is dropped by count as before: a replay of it runs
        auto bare = sent(registry, "ev-3", std::nullopt);
        ASSERT_FALSE(registry.start(bare));
        registry.finish(bare, "completed", 200, 1);
        auto last = sent(registry, "ev-4", std::nullopt);
        ASSERT_FALSE(registry.start(last));
        registry.finish(last, "completed", 200, 1);
        EXPECT_FALSE(registry.start(sent(registry, "ev-3", std::nullopt)));
        // the held one stays until its not_after_ms + skew (2 s), then frees the table
        clock.ms = 62'000;
        EXPECT_TRUE(registry.start(sent(registry, "ev-1", not_after)));
        clock.ms = 62'001;
        EXPECT_TRUE(registry.start(sent(registry, "ev-1", not_after))->expired);
        EXPECT_EQ(404, registry.cancel("unknown-2", 0).first);
    }
    {
        // the cap: held at most tombstone_max_s after the finish (stated: a not_after_ms
        // further out is not covered); a refused copy with a later not_after_ms extends the
        // hold of a running attempt, applied from its finish
        FakeClock clock;
        AttemptSettings settings = settings_with(&clock, 1);
        settings.clock_skew_ms = 0;
        settings.tombstone_max_s = 100;
        AttemptRegistry registry(settings);
        auto far = sent(registry, "far", clock.wall_ms() + 500'000);
        ASSERT_FALSE(registry.start(far));
        registry.finish(far, "completed", 200, 1);
        clock.ms = 100'000;
        EXPECT_TRUE(registry.start(sent(registry, "far", clock.wall_ms() + 400'000)));
        // that refusal extended it (the copy's not_after_ms, capped 100 s from then)
        clock.ms = 200'001;
        EXPECT_FALSE(registry.start(sent(registry, "far", clock.wall_ms() + 300'000)));

        const uint64_t own = clock.wall_ms() + 10'000;
        auto running = sent(registry, "run", own);
        ASSERT_FALSE(registry.start(running));
        EXPECT_TRUE(registry.start(sent(registry, "run", own + 20'000)));
        registry.finish(running, "completed", 200, 1);
        clock.ms += 25'000;   // past its own not_after_ms, before the copy's
        EXPECT_TRUE(registry.start(sent(registry, "run", own + 20'000)));
    }
    {
        // retention 0: nothing kept, a finished state promises nothing against a replay
        FakeClock clock;
        AttemptRegistry registry(settings_with(&clock, 0));
        const uint64_t not_after = clock.wall_ms() + 60'000;
        auto original = sent(registry, "z", not_after);
        ASSERT_FALSE(registry.start(original));
        registry.finish(original, "completed", 200, 1);
        EXPECT_FALSE(registry.start(sent(registry, "z", not_after)));
        EXPECT_NE(std::string::npos, registry.capabilities_json()["release_rule"].asString()
                                             .find("retention_s is 0"));
    }
}

// Review of levels 4-5, finding 6 (the reviewer's finish-restart-replay probe; the restarted
// process is a second registry, with its own server_instance): a finished attempt's hold is its
// process's. A request finished with not_after_ms 60 s ahead and replayed within its hold is
// refused by that process whether pinned or not; after a restart the unpinned replay runs again
// (on the reviewer's server: 200, 170 work units), which the release rule now states — a finished
// state is replay-safe only for an attempt sent with expect_server_instance — and the pinned one
// is refused there, 409 instance_mismatch
TEST(GraphletAttemptRegistry, AFinishedStateIsReplaySafeOnlyWhenPinned) {
    FakeClock clock;
    AttemptRegistry first(settings_with(&clock, 1, 1));
    const uint64_t not_after = clock.wall_ms() + 60'000;
    auto sent = [&](AttemptRegistry &registry, const std::string &id,
                    const std::string &instance) {
        auto a = attempt_of(registry, id);
        AttemptIds ids = a->ids();
        ids.not_after_ms = not_after;
        ids.expect_server_instance = instance;
        a->set_ids(ids);
        return a;
    };
    const std::string old_instance = first.server_instance();
    auto unpinned = sent(first, "finished-unpinned", "");
    ASSERT_FALSE(first.start(unpinned));
    first.finish(unpinned, "completed", 200, 1);
    auto pinned = sent(first, "finished-pinned", old_instance);
    ASSERT_FALSE(first.start(pinned));
    first.finish(pinned, "completed", 200, 1);
    // past retention_s and retention_count: both held by this process through their hold
    clock.ms = 1'100;
    for (const auto &[id, instance] : { std::make_pair("finished-unpinned", std::string()),
                                        std::make_pair("finished-pinned", old_instance) }) {
        auto replay = first.start(sent(first, id, instance));
        ASSERT_TRUE(replay) << id;
        EXPECT_EQ("finished", replay->body["state"].asString()) << id;
    }
    // the restart: a new server_instance, no holds
    AttemptRegistry restarted(settings_with(&clock, 1, 1));
    ASSERT_NE(old_instance, restarted.server_instance());
    auto again = sent(restarted, "finished-unpinned", "");
    EXPECT_FALSE(restarted.start(again));   // it runs again
    restarted.finish(again, "completed", 200, 1);
    auto refused = restarted.start(sent(restarted, "finished-pinned", old_instance));
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->instance_mismatch);
    EXPECT_EQ("instance_mismatch", refused->body["state"].asString());
    EXPECT_EQ(404, restarted.state("finished-pinned").first);
    // what the contract states: the hold is the process's, a finished state is replay-safe only
    // when pinned, and what an unpinned one assumes
    const Json::Value caps = restarted.capabilities_json();
    const std::string rule = caps["release_rule"].asString();
    for (const char *phrase : { "The hold is this process's, in memory: a restarted process (a "
                                "new server_instance) holds none",
                                "replay-safe only as stated next",
                                "So a finished state is replay-safe",
                                "only for an attempt sent with expect_server_instance equal to "
                                "this server_instance and with not_after_ms",
                                "pins the instance, as for a tombstone",
                                "A finished state of an attempt sent without "
                                "expect_server_instance assumes that no copy of the request "
                                "reaches a restarted process",
                                // kept from the review of the pass-5 fixes
                                "A finished attempt's id stays refused (409)",
                                "assumes that no copy of the request arrives after the attempt "
                                "left retention" }) {
        EXPECT_NE(std::string::npos, rule.find(phrase)) << phrase;
    }
    EXPECT_NE(std::string::npos, caps["instance"].asString().find(
            "a delayed copy of a cancelled or finished request would otherwise run there"));
    AttemptRegistry none(settings_with(&clock, 0));
    EXPECT_NE(std::string::npos, none.capabilities_json()["release_rule"].asString().find(
            "a copy sent with expect_server_instance is still refused by a restarted process"));
}

// The review of 2026-10-06 (X2, the search service's release parity LRG-R1/R2, and C16, C20,
// C24, C30, X4 on the same strings), text only: the clock release takes the attempt as stopped
// at its bound, which it is apart from what it runs past it until its next delivery check, of
// no stated length (the review of the P2 fixes: not "one uninterruptible step", since only the
// delivery checks compare the bound), and an answer of running or stopping past that instant
// shows that run; the 409 refusing a copy
// carries the id's state as GET answers it, so its finished state is a finished state; an
// expired 409 releases, because the id's registration is checked before the expiry (the order
// the not_after text states), and settles nothing; refused copies extend a finished attempt's
// hold; expect_server_instance has the ids' pattern and a 400; the bound names its cap and the
// walk-until it states. With retention and without
TEST(GraphletAttemptRegistry, ReleaseTextsStateTheOverrunAndTheGrounds) {
    FakeClock clock;
    for (uint64_t retention_s : { uint64_t(60), uint64_t(0) }) {
        AttemptRegistry registry(settings_with(&clock, retention_s));
        const Json::Value caps = registry.capabilities_json();
        const std::string what = "retention_s " + std::to_string(retention_s);
        auto has = [&](const char *field, const char *phrase) {
            EXPECT_NE(std::string::npos, caps[field].asString().find(phrase))
                << what << ", " << field << ": " << phrase;
        };
        for (const char *phrase : { "the cap whatever the seeds' budgets (usage.bound.capped_by: "
                                    "content_timeout when it applied)",
                                    "the walk stops at its first poll that reads the clock "
                                    "after it",
                                    "else the lowest walk-until seen: the lowest that such a "
                                    "poll compared with",
                                    "not a lower one in force only between two such polls",
                                    "The bound itself is compared only at the delivery checks",
                                    "past it the attempt runs on until its next delivery check "
                                    "(then the 503) or its handler's return",
                                    "when its walk had not stopped by then, the walk up to its "
                                    "next poll that reads the clock, the stopped seed's "
                                    "finalisation and the building of its result up to the "
                                    "first delivery check",
                                    "a run of no stated length" }) {
            has("bound", phrase);
        }
        // the overrun is not one step (the review of the P2 fixes)
        for (const char *field : { "bound", "not_after", "release_rule" }) {
            EXPECT_EQ(std::string::npos, caps[field].asString().find("uninterruptible step"))
                << what << ", " << field;
            EXPECT_EQ(std::string::npos, caps[field].asString().find("first poll after"))
                << what << ", " << field;
        }
        EXPECT_EQ(std::string::npos, caps["delivery_reserve"]["rule"].asString()
                                             .find("first poll after")) << what;
        for (const char *phrase : { "It is the last of the refusals, all made under one lock",
                                    "a request naming another expect_server_instance is refused "
                                    "first (409 instance_mismatch)",
                                    "whatever its not_after_ms (the duplicate's 409",
                                    "so an expired 409 means that no attempt with that id existed "
                                    "on this server_instance when it was judged",
                                    "cannot start subsequently",
                                    "as past its bound once its clock passes that + bound_ms — "
                                    "apart from what it runs past its bound, up to its next "
                                    "delivery check, a run of no stated length" }) {
            has("not_after", phrase);
        }
        EXPECT_EQ(std::string::npos, caps["not_after"].asString().find("as stopped once")) << what;
        for (const char *phrase : { "a string matching id_pattern",
                                    "checked first (ahead of a duplicate id and of not_after_ms)",
                                    "and the field without attempt_id are a 400 naming the field, "
                                    "without usage" }) {
            has("instance", phrase);
        }
        for (const char *phrase : { "the 409 refusing a copy of the request whose attempt — the "
                                    "id's state, as GET answers it, with its attempt_id and "
                                    "server_instance — says state finished",
                                    "on an expired 409 (below)",
                                    "An expired 409 (state: expired) for an attempt sent with "
                                    "exactly that not_after_ms",
                                    "checks the id's registration (running, retained, held, "
                                    "tombstoned) before the expiry, under one lock",
                                    "It settles nothing",
                                    "the 409's server_time_ms - not_after_ms is the step it "
                                    "survives",
                                    "The clock release takes the attempt as stopped",
                                    "past its bound apart from what it runs past it until its "
                                    "next delivery check",
                                    "a run of no stated length",
                                    "still answers running or stopping after that instant shows "
                                    "that run still going",
                                    "a 409 whose attempt says running or stopping or is a "
                                    "tombstone, and a 409 instance_mismatch release nothing",
                                    "no copy of the request reaches another server that serves "
                                    "the same ledger",
                                    "for the clock release, that the attempt's run past its "
                                    "bound has ended by then",
                                    // kept
                                    "replay-safe only as stated next",
                                    "not_after_ms + clock_skew_allowance_ms + bound_ms" }) {
            has("release_rule", phrase);
        }
        if (retention_s) {
            for (const char *phrase : { "at most tombstone_max_s after it finished or after its "
                                        "latest refused copy",
                                        "their number is at most the finishes and refused copies "
                                        "within tombstone_max_s" }) {
                has("release_rule", phrase);
            }
            EXPECT_EQ(std::string::npos, caps["release_rule"].asString().find(
                    "tombstone_max_s bounds how many are held")) << what;
            EXPECT_NE(std::string::npos,
                      registry.retention_text().find("after it finished or after its latest "
                                                     "refused copy")) << what;
        }
    }
}

// Review of the pass-5 fixes, finding 2: with clock_skew_ms 0 a cancel's covers_admission
// (suppressed_until_ms == not_after_ms) must hold at the instant the wall clock reads
// not_after_ms, which the strict not_after check still admits — the tombstone is live through
// suppressed_until_ms inclusive, also after its steady hold passed
TEST(GraphletAttemptRegistry, TombstonesHoldThroughSuppressedUntilInclusive) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 1);
    settings.clock_skew_ms = 0;
    AttemptRegistry registry(settings);
    const uint64_t not_after = clock.wall_ms() + 1200;
    auto [status, body] = registry.cancel("sk-0", 0, not_after);
    ASSERT_EQ(404, status);
    EXPECT_TRUE(body["covers_admission"].asBool());
    EXPECT_EQ(not_after, body["suppressed_until_ms"].asUInt64());
    auto copy = [&]() {
        auto a = attempt_of(registry, "sk-0");
        AttemptIds ids = a->ids();
        ids.not_after_ms = not_after;
        a->set_ids(ids);
        return a;
    };
    // the steady hold over, the wall clock (one ms behind) reading not_after_ms exactly
    clock.ms = 1'201;
    clock.wall_step = -1;
    ASSERT_EQ(not_after, clock.wall_ms());
    EXPECT_FALSE(not_after_passed(copy()->ids(), clock.wall_ms()));
    auto refused = registry.start(copy());
    ASSERT_TRUE(refused);
    EXPECT_FALSE(refused->expired);
    EXPECT_TRUE(refused->body["tombstone"].asBool());
    // the refusal held it again through not_after_ms (and its retention): once the clock read
    // later, the copy is refused by its own not_after_ms
    clock.ms = 2'202;
    clock.wall_step = 0;
    refused = registry.start(copy());
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->expired);
}

// Review of the pass-5 fixes, finding 3: a refused copy's suppression is judged against its own
// not_after_ms alone — a copy without one is not covered (no_not_after_ms), whatever not_after_ms
// a cancel named, and indeed runs once re-sent after the tombstone
TEST(GraphletAttemptRegistry, ARefusedCopyIsJudgedByItsOwnNotAfterMs) {
    FakeClock clock;
    AttemptSettings settings = settings_with(&clock, 1);
    settings.clock_skew_ms = 100;
    AttemptRegistry registry(settings);
    const uint64_t not_after = clock.wall_ms() + 1500;
    auto [status, body] = registry.cancel("nna-1", 0, not_after);
    ASSERT_EQ(404, status);
    EXPECT_TRUE(body["covers_admission"].asBool());
    auto refused = registry.start(attempt_of(registry, "nna-1"));   // no not_after_ms
    ASSERT_TRUE(refused);
    EXPECT_TRUE(refused->body["tombstone"].asBool());
    EXPECT_FALSE(refused->body.isMember("not_after_ms"));
    EXPECT_FALSE(refused->body["covers_admission"].asBool());
    EXPECT_EQ("no_not_after_ms", refused->body["covers_admission_reason"].asString());
    // GET and a repeated cancel still answer for the largest not_after_ms named
    EXPECT_TRUE(registry.state("nna-1").second["covers_admission"].asBool());
    const uint64_t until = refused->body["suppressed_until_ms"].asUInt64();
    clock.ms = static_cast<int64_t>(until - FakeClock::kWallBase) + 1;
    EXPECT_FALSE(registry.start(attempt_of(registry, "nna-1")));
}

TEST(GraphletAttemptRegistry, RetentionByCountAndAge) {
    FakeClock clock;
    AttemptRegistry registry(settings_with(&clock, 100, 3));
    std::vector<std::shared_ptr<Attempt>> done;
    for (int i = 0; i < 5; ++i) {
        auto a = attempt_of(registry, "r" + std::to_string(i));
        ASSERT_FALSE(registry.start(a));
        clock.ms += 10;
        registry.finish(a, "completed", 200, 1);
        done.push_back(a);
    }
    // the last three finished attempts are kept
    EXPECT_EQ(404, registry.state("r0").first);
    EXPECT_EQ(404, registry.state("r1").first);
    for (int i = 2; i < 5; ++i) {
        EXPECT_EQ(200, registry.state("r" + std::to_string(i)).first) << i;
    }
    // a running attempt is never dropped, however old
    auto running = attempt_of(registry, "run");
    ASSERT_FALSE(registry.start(running));
    clock.ms += 100'000;
    for (int i = 2; i < 5; ++i) {
        EXPECT_EQ(404, registry.state("r" + std::to_string(i)).first) << i;
    }
    EXPECT_EQ(200, registry.state("run").first);
}

// wait_ms: the cancel answers once the attempt finished, within the wait
TEST(GraphletAttemptRegistry, CancelWaitsForTheEnd) {
    AttemptRegistry registry(settings_with(nullptr));
    auto a = attempt_of(registry, "w");
    a->set_bound(1, 60'000);
    ASSERT_FALSE(registry.start(a));
    std::thread walker([&]() {
        // the walk polls until it sees the stop, then the handler finishes
        while (a->poll() == ExternalStop::NONE) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(30));
        registry.finish(a, a->walk_reason(), 200, 10);
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    auto [status, body] = registry.cancel("w", 5000);
    walker.join();
    EXPECT_EQ(200, status);
    EXPECT_TRUE(body["cancelled"].asBool());
    EXPECT_EQ("finished", body["state"].asString());
    EXPECT_EQ("cancelled", body["attempt"]["reason"].asString());
}

// Many threads starting, cancelling, reading and finishing attempts: an id runs once, every
// started attempt finishes once with its reason, its state only moves forward, and a cancel is
// 200 exactly while the attempt has not finished
TEST(GraphletAttemptRegistry, ConcurrentAttempts) {
    AttemptRegistry registry(settings_with(nullptr, 3600, 100'000));
    constexpr int kThreads = 8;
    constexpr int kIds = 64;
    std::atomic<int> started { 0 }, conflicts { 0 }, cancelled_ok { 0 }, bad { 0 };
    std::vector<std::thread> threads;
    std::mutex m;
    std::map<std::string, int> runs;
    for (int t = 0; t < kThreads; ++t) {
        threads.emplace_back([&, t]() {
            std::mt19937 rng(t);
            for (int i = 0; i < 300; ++i) {
                const std::string id = "id" + std::to_string(rng() % kIds) + "-"
                                     + std::to_string(i % 5);
                switch (rng() % 3) {
                    case 0: {
                        auto a = attempt_of(registry, id);
                        a->set_bound(1, 1000);
                        if (registry.start(a)) {
                            conflicts++;
                            break;
                        }
                        {
                            std::lock_guard<std::mutex> lock(m);
                            runs[id]++;
                        }
                        started++;
                        std::string last = "running";
                        for (int p = 0; p < 3; ++p) {
                            a->poll();
                            const std::string s = registry.state(id).second["state"].asString();
                            // running -> stopping -> finished, never back
                            if ((last == "stopping" && s == "running") || s == "finished")
                                bad++;
                            last = s;
                        }
                        registry.finish(a, a->walk_reason(), 200, 1);
                        if (registry.state(id).second["state"].asString() != "finished")
                            bad++;
                        break;
                    }
                    case 1: {
                        const int status = registry.cancel(id, rng() % 2 ? 0 : 5).first;
                        if (status == 200) {
                            cancelled_ok++;
                        } else if (status != 404) {
                            bad++;
                        }
                        break;
                    }
                    default: {
                        const int status = registry.state(id).first;
                        if (status != 200 && status != 404)
                            bad++;
                    }
                }
            }
        });
    }
    for (auto &t : threads) {
        t.join();
    }
    EXPECT_EQ(0, bad.load());
    EXPECT_GT(started.load(), 50);
    EXPECT_GT(conflicts.load(), 0);
    for (const auto &[id, n] : runs) {
        EXPECT_EQ(1, n) << id;
    }
}

namespace {

// a connected TCP pair on the loopback: (client, server side)
std::pair<int, int> tcp_pair() {
    const int listener = ::socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in addr {};
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    addr.sin_port = 0;
    EXPECT_EQ(0, ::bind(listener, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)));
    EXPECT_EQ(0, ::listen(listener, 1));
    socklen_t len = sizeof(addr);
    EXPECT_EQ(0, ::getsockname(listener, reinterpret_cast<sockaddr*>(&addr), &len));
    const int client = ::socket(AF_INET, SOCK_STREAM, 0);
    EXPECT_EQ(0, ::connect(client, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)));
    const int server = ::accept(listener, nullptr, nullptr);
    ::close(listener);
    return { client, server };
}

// the peer's close takes a moment to arrive
bool becomes_closed(int fd) {
    for (int i = 0; i < 200; ++i) {
        if (peer_closed(fd))
            return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return false;
}

} // namespace

// B1's check: a peek that consumes nothing tells a connected client (idle, or with a request
// waiting) from one that closed, half-closed or reset its connection
TEST(GraphletServer, PeerClosedTellsAGoneClient) {
    {
        auto [client, server] = tcp_pair();
        EXPECT_FALSE(peer_closed(server));
        const char data[] = "GET / HTTP/1.1\r\n\r\n";
        ASSERT_EQ(static_cast<ssize_t>(sizeof(data) - 1), ::send(client, data, sizeof(data) - 1, 0));
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        EXPECT_FALSE(peer_closed(server)) << "a pipelined request waiting";
        // the peek consumed nothing
        char buf[64];
        EXPECT_EQ(static_cast<ssize_t>(sizeof(data) - 1), ::recv(server, buf, sizeof(buf), 0));
        EXPECT_FALSE(peer_closed(server));
        ::close(client);
        EXPECT_TRUE(becomes_closed(server)) << "closed";
        ::close(server);
    }
    {
        auto [client, server] = tcp_pair();
        ::shutdown(client, SHUT_WR);
        EXPECT_TRUE(becomes_closed(server)) << "half-closed";
        ::close(client);
        ::close(server);
    }
    {
        auto [client, server] = tcp_pair();
        linger l { 1, 0 };
        ::setsockopt(client, SOL_SOCKET, SO_LINGER, &l, sizeof(l));
        ::close(client);
        EXPECT_TRUE(becomes_closed(server)) << "reset";
        ::close(server);
    }
    EXPECT_TRUE(peer_closed(-1));
}

// The review of the stage-4 backend, F3: bytes waiting past the request — a trailing CRLF (RFC
// 9112 §2.2), a pipelined request — hide a close behind them from a peek, which returns the
// bytes: the connection's TCP state tells. A client that sent them and stays connected is not
// gone; nothing is consumed either way
TEST(GraphletServer, PeerClosedSeesACloseBehindWaitingBytes) {
#if defined(__linux__) || defined(__APPLE__)
    for (const std::string &waiting : { std::string("\r\n"),
                                        std::string("GET /stats HTTP/1.1\r\nHost: x\r\n\r\n") }) {
        for (const char *how : { "close", "half-close", "reset" }) {
            auto [client, server] = tcp_pair();
            ASSERT_EQ(static_cast<ssize_t>(waiting.size()),
                      ::send(client, waiting.data(), waiting.size(), 0));
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
            EXPECT_FALSE(peer_closed(server)) << how << ": connected, bytes waiting";
            if (std::string(how) == "reset") {
                linger l { 1, 0 };
                ::setsockopt(client, SOL_SOCKET, SO_LINGER, &l, sizeof(l));
                ::close(client);
            } else if (std::string(how) == "half-close") {
                ::shutdown(client, SHUT_WR);
            } else {
                ::close(client);
            }
            EXPECT_TRUE(becomes_closed(server)) << how << " behind " << waiting.size()
                                                << " waiting byte(s)";
            if (std::string(how) != "reset") {
                // the peek consumed nothing: the bytes are still there to read
                char buf[64];
                EXPECT_EQ(static_cast<ssize_t>(waiting.size()),
                          ::recv(server, buf, sizeof(buf), MSG_DONTWAIT)) << how;
            }
            if (std::string(how) == "half-close")
                ::close(client);
            ::close(server);
        }
    }
#endif
}

// The server's own shutdown (the HTTP server's content timeout) behind waiting bytes: Linux
// keeps the bytes, so the peek returns them and the state is FIN_WAIT2, which read "connected"
// while the client kept its end open — the walk computed on past the timeout (review of
// 2026-10-06, U13-02). macOS discards them on SHUT_RD: the peek reads the end. Gone on both.
// A descriptor that holds no connection — not a socket, or closed — is gone too, as the
// header says (ENOTSOCK and EBADF read "connected" before)
TEST(GraphletServer, PeerClosedSeesTheServersOwnShutdown) {
#if defined(__linux__) || defined(__APPLE__)
    for (const std::string &waiting : { std::string("\r\n"),
                                        std::string("GET /stats HTTP/1.1\r\nHost: x\r\n\r\n") }) {
        auto [client, server] = tcp_pair();
        ASSERT_EQ(static_cast<ssize_t>(waiting.size()),
                  ::send(client, waiting.data(), waiting.size(), 0));
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        EXPECT_FALSE(peer_closed(server)) << "connected, " << waiting.size() << " byte(s) waiting";
        ::shutdown(server, SHUT_RDWR);
        // the client keeps its end open: only the server's state tells
        EXPECT_TRUE(becomes_closed(server)) << "own shutdown behind " << waiting.size()
                                            << " waiting byte(s)";
        ::close(client);
        ::close(server);
    }
    {
        // nothing waiting: the peek reads the end on both platforms (as before)
        auto [client, server] = tcp_pair();
        ::shutdown(server, SHUT_RDWR);
        EXPECT_TRUE(becomes_closed(server)) << "own shutdown, nothing waiting";
        ::close(client);
        ::close(server);
    }
#endif
    int pipe_fds[2];
    ASSERT_EQ(0, ::pipe(pipe_fds));
    EXPECT_TRUE(peer_closed(pipe_fds[0])) << "not a socket";
    ::close(pipe_fds[0]);
    ::close(pipe_fds[1]);
    // the descriptor just closed: nothing in this test opens another one before the call
    EXPECT_TRUE(peer_closed(pipe_fds[0])) << "a closed descriptor";
}

// The response text is written under the attempt's check, byte for byte what the writer wrote
// before, and the check's exception stops it
TEST(GraphletServer, CheckedWriterIsByteIdentical) {
    Json::Value v;
    std::mt19937 rng(7);
    for (int i = 0; i < 3000; ++i) {
        Json::Value e;
        e["id"] = i;
        e["x"] = static_cast<double>(rng()) / 7.0;
        e["s"] = std::string(rng() % 200, static_cast<char>('a' + rng() % 26)) + "\"\n\x01\xc3\xa9";
        e["a"].append(Json::Value::null);
        e["a"].append(true);
        v["results"].append(e);
    }
    for (bool compact : { true, false }) {
        Json::StreamWriterBuilder builder;
        if (compact)
            builder["indentation"] = "";
        const std::string expected = Json::writeString(builder, v);
        size_t checks = 0;
        EXPECT_EQ(expected, json_text(v, compact, [&]() { checks++; }));
        EXPECT_EQ(expected, json_text(v, compact));
        EXPECT_GE(checks, expected.size() / 65536 - 1);
        EXPECT_GT(checks, 0u);
        struct Stop : std::runtime_error {
            using std::runtime_error::runtime_error;
        };
        size_t calls = 0;
        EXPECT_THROW(json_text(v, compact, [&]() {
                         if (++calls == 2)
                             throw Stop("stop");
                     }),
                     Stop);
        EXPECT_EQ(2u, calls);
    }
}

// Review of pass 5, finding 6 (the reviewer's pass5_delivery_probe.cpp): one JSON string value
// of 16 MiB was appended whole by the checked writer, and a seed's 16 MiB text whole by the
// response's assembly — one check each. Both now copy in pieces up to the next check: at
// least size / 64 KiB checks, the bytes unchanged, and the longest stretch between two checks
// (the value's escaping, before it is copied) measured
TEST(GraphletServer, DeliveryChecksEvery64KiBOfALargeToken) {
    Json::Value v;
    v["text"] = std::string(16 * 1024 * 1024, 'A');
    size_t checks = 0;
    double gap = 0;
    const std::string text = json_text(v, true, [&]() { checks++; }, &gap);
    EXPECT_EQ(json_text(v, true), text);
    EXPECT_GE(checks, text.size() / kDeliveryCheckBytes);
    EXPECT_GT(gap, 0);
    Json::Value envelope(Json::objectValue);
    envelope["algorithm_version"] = "x";
    envelope["timing"]["elapsed_ms"] = 1;
    Json::Value whole = envelope;
    whole["results"].append(v);
    whole["results"].append(Json::Value(Json::objectValue));
    checks = 0;
    gap = 0;
    const std::string response = assemble_traverse_response(envelope, { text, "{}" },
                                                            [&]() { checks++; }, &gap);
    EXPECT_EQ(json_text(whole, true), response);
    EXPECT_EQ(assemble_traverse_response(envelope, { text, "{}" }), response);
    EXPECT_GE(checks, text.size() / kDeliveryCheckBytes);
    EXPECT_GT(gap, 0);
    // a check's exception stops the copying between two pieces
    struct Stop : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    size_t calls = 0;
    EXPECT_THROW(assemble_traverse_response(envelope, { text }, [&]() {
                     if (++calls == 3)
                         throw Stop("stop");
                 }),
                 Stop);
    EXPECT_EQ(3u, calls);
}

// The review of 2026-10-06, C9 (U12-02, U13-03, X-EFFICIENCY-02): the server held 2.5x (gzip) to
// 3.5x (identity) a /traverse response's text while delivering it — the per-seed texts lived
// until the handler returned, beside the assembled copy, which grew by doubling. The texts are
// now moved into the assembly, each freed once copied, and the response is reserved at its exact
// size: its capacity is its size (no doubling step held an old and a new buffer), and its bytes
// are unchanged, with and without checks, whichever side of "results" the envelope's members
// are on
TEST(GraphletServer, AssemblyReservesTheResponseOnce) {
    std::mt19937 rng(17);
    std::vector<std::string> texts;
    for (size_t i = 0; i < 5; ++i) {
        Json::Value r;
        r["seed"]["seed_id"] = "s" + std::to_string(i);
        r["text"] = std::string(200'000 + rng() % 100'000, static_cast<char>('a' + i));
        texts.push_back(json_text(r, true));
    }
    Json::Value both(Json::objectValue), before(Json::objectValue), after(Json::objectValue);
    both["algorithm_version"] = "x";
    both["usage"]["seeds"]["started"] = 5;
    before["algorithm_version"] = "x";
    after["timing"]["elapsed_ms"] = 1.5;
    auto parse = [](const std::string &t) {
        Json::Value v;
        std::unique_ptr<Json::CharReader> reader(Json::CharReaderBuilder().newCharReader());
        std::string errors;
        EXPECT_TRUE(reader->parse(t.data(), t.data() + t.size(), &v, &errors)) << errors;
        return v;
    };
    for (const Json::Value *envelope : { &both, &before, &after }) {
        for (size_t n : { size_t(0), size_t(1), texts.size() }) {
            const std::vector<std::string> some(texts.begin(), texts.begin() + n);
            Json::Value whole = *envelope;
            whole["results"] = Json::Value(Json::arrayValue);
            for (const std::string &t : some) {
                whole["results"].append(parse(t));
            }
            const std::string expected = json_text(whole, true);
            for (bool checked : { false, true }) {
                std::vector<std::string> moved = some;
                size_t checks = 0;
                const std::string response = checked
                    ? assemble_traverse_response(*envelope, std::move(moved), [&]() { checks++; })
                    : assemble_traverse_response(*envelope, std::move(moved));
                EXPECT_EQ(expected, response) << n << " " << checked;
                // reserved once at the size written (a string's capacity rounds up by less
                // than 16 bytes; one doubling step would have left up to the size again)
                if (n) {
                    EXPECT_LE(response.capacity(), response.size() + 16) << n << " " << checked;
                }
                if (checked && n == texts.size()) {
                    EXPECT_GE(checks, expected.size() / kDeliveryCheckBytes - 1);
                }
            }
        }
    }
}

// zlib counts its input in 32 bits (uInt avail_in): a text of 4 GiB or more was cut to its size
// modulo 2^32 and the server answered 200 with a well-formed stream of that prefix (review of
// 2026-10-06, U13-04). compress_string hands the text over in pieces of at most 2^32 - 1 bytes;
// tested with small pieces, which take the same path: the stream inflates to the whole text
// across every boundary, in both containers, under the check; a text of one piece is compressed
// by the same calls as before, so its bytes are those of the whole-text call
TEST(GraphletServer, CompressionTakesTheTextInPieces) {
    auto inflate_all = [](const std::string &compressed) {
        z_stream zs;
        memset(&zs, 0, sizeof(zs));
        // 15 + 32: a zlib or a gzip header, whichever the stream has
        EXPECT_EQ(Z_OK, inflateInit2(&zs, 15 + 32));
        zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(compressed.data()));
        zs.avail_in = static_cast<uInt>(compressed.size());
        std::string out;
        char buffer[65536];
        int ret;
        do {
            zs.next_out = reinterpret_cast<Bytef *>(buffer);
            zs.avail_out = sizeof(buffer);
            ret = inflate(&zs, Z_NO_FLUSH);
            out.append(buffer, sizeof(buffer) - zs.avail_out);
        } while (ret == Z_OK);
        EXPECT_EQ(Z_STREAM_END, ret);
        inflateEnd(&zs);
        return out;
    };
    std::mt19937 rng(5);
    std::string text;
    for (size_t i = 0; i < 300'000; ++i) {
        text.push_back("ACGT{}\":,0123456789"[rng() % 19]);
    }
    for (bool gzip : { false, true }) {
        const std::string whole = compress_string(text, 1, gzip);
        EXPECT_EQ(text, inflate_all(whole)) << gzip;
        EXPECT_EQ(whole, compress_string(text, 1, gzip, nullptr, text.size())) << gzip;
        for (size_t piece : { size_t(1), size_t(7), size_t(32768), size_t(65537),
                              text.size() - 1 }) {
            size_t checks = 0;
            const std::string compressed
                = compress_string(text, 1, gzip, [&]() { checks++; }, piece);
            EXPECT_EQ(text, inflate_all(compressed)) << gzip << " " << piece;
            EXPECT_GT(checks, 0u);
        }
        EXPECT_EQ("", inflate_all(compress_string("", 1, gzip, nullptr, 1))) << gzip;
        EXPECT_EQ("x", inflate_all(compress_string("x", 9, gzip, nullptr, 1))) << gzip;
    }
    // a check's exception leaves the call (the stream released) between two blocks
    struct Stop : std::runtime_error {
        using std::runtime_error::runtime_error;
    };
    size_t calls = 0;
    EXPECT_THROW(compress_string(text, 1, true, [&]() {
                     if (++calls == 3)
                         throw Stop("stop");
                 }, 1000),
                 Stop);
    EXPECT_EQ(3u, calls);
}

// The multi-graph list (pass 5, W2): three columns read as they always were, two optional
// ones (manifest_path, index_ns), empty meaning none; more than five columns, fewer than
// three, or an index_ns that is no token refuse the line
TEST(GraphletServer, GraphListLinesAreParsedStrictly) {
    GraphListEntry e = parse_graph_list_line("uhgg,/d/g.dbg,/d/a.annodbg", 3);
    EXPECT_EQ(3u, e.line);
    EXPECT_EQ("uhgg", e.name);
    EXPECT_EQ("/d/g.dbg", e.graph_path);
    EXPECT_EQ("/d/a.annodbg", e.annotation_path);
    EXPECT_EQ("", e.manifest_path);
    EXPECT_EQ("", e.index_ns);
    // as before: an empty annotation column is read as given (the loader names it)
    EXPECT_EQ("", parse_graph_list_line("n,g,", 1).annotation_path);
    e = parse_graph_list_line("n,g,a,/m/a.manifest.json", 1);
    EXPECT_EQ("/m/a.manifest.json", e.manifest_path);
    EXPECT_EQ("", e.index_ns);
    e = parse_graph_list_line("n,g,a,/m/a.manifest.json,uhgg_k31", 1);
    EXPECT_EQ("uhgg_k31", e.index_ns);
    // a name without a manifest, and empty optional columns
    e = parse_graph_list_line("n,g,a,,uhgg.v2", 1);
    EXPECT_EQ("", e.manifest_path);
    EXPECT_EQ("uhgg.v2", e.index_ns);
    e = parse_graph_list_line("n,g,a,,", 1);
    EXPECT_EQ("", e.manifest_path);
    EXPECT_EQ("", e.index_ns);
    for (const char *bad : { "n", "n,g", "n,g,a,m,ns,extra", "n,g,a,m,bad ns", "n,g,a,m,a/b" }) {
        try {
            parse_graph_list_line(bad, 7);
            ADD_FAILURE() << "accepted " << bad;
        } catch (const std::invalid_argument &ex) {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find("line 7")) << ex.what();
        }
    }
}

// One index, one identity: lines naming the same pair must agree (an empty column states
// nothing), each (pair, manifest) is checked once, and a conflict names both lines
TEST(GraphletServer, GraphListIdentitiesAgreePerPair) {
    std::vector<GraphListEntry> entries {
        parse_graph_list_line("A,g1,a1,m1,ns1", 1),
        parse_graph_list_line("B,g1,a1", 2),            // the same pair under another name
        parse_graph_list_line("B,g2,a2,m2", 3),
        parse_graph_list_line("C,g3,a3", 4),
        parse_graph_list_line("D,g1,a1,m1", 5),          // the same manifest again
    };
    std::map<std::string, int> calls;
    auto fp = [&](const GraphListEntry &e) {
        calls[e.manifest_path]++;
        return "fp-" + e.manifest_path;
    };
    auto ids = graph_list_identities(entries, fp);
    EXPECT_EQ(3u, ids.size());
    EXPECT_EQ(std::make_pair(std::string("ns1"), std::string("fp-m1")), ids.at({ "g1", "a1" }));
    EXPECT_EQ(std::make_pair(std::string(""), std::string("fp-m2")), ids.at({ "g2", "a2" }));
    EXPECT_EQ(std::make_pair(std::string(""), std::string("")), ids.at({ "g3", "a3" }));
    EXPECT_EQ(1, calls["m1"]);
    EXPECT_EQ(1, calls["m2"]);

    // a different name for the same pair
    entries.push_back(parse_graph_list_line("E,g1,a1,,ns2", 6));
    try {
        graph_list_identities(entries, fp);
        ADD_FAILURE() << "accepted two names for one index";
    } catch (const std::invalid_argument &ex) {
        EXPECT_NE(std::string::npos, std::string(ex.what()).find("lines 1 and 6")) << ex.what();
    }
    entries.pop_back();
    // a different manifest digest for the same pair
    entries.push_back(parse_graph_list_line("E,g1,a1,m9", 6));
    EXPECT_THROW(graph_list_identities(entries, fp), std::invalid_argument);
    entries.pop_back();
    // a manifest that does not describe the pair's files: the fingerprint's error
    entries.push_back(parse_graph_list_line("F,g4,a4,bad", 7));
    EXPECT_THROW(graph_list_identities(entries, [&](const GraphListEntry &e) -> std::string {
        if (e.manifest_path == "bad")
            throw std::invalid_argument("the loaded file g4 is not listed");
        return fp(e);
    }), std::invalid_argument);
    entries.pop_back();
    // Review of pass 5: one index is one pair of files however its paths are spelled — a
    // second spelling takes the pair's identity, and may not state another one
    entries.push_back(parse_graph_list_line("G,./g1,./a1", 8));
    ids = graph_list_identities(entries, fp);
    EXPECT_EQ(std::make_pair(std::string("ns1"), std::string("fp-m1")), ids.at({ "./g1", "./a1" }));
    entries.push_back(parse_graph_list_line("H,./g1,a1,m9", 9));
    EXPECT_THROW(graph_list_identities(entries, fp), std::invalid_argument);
    entries.pop_back();
    // ... and two different pairs never state one index_fp (one manifest written for both, as
    // for two annotations of one graph with swapped memberships, which sizes cannot tell apart)
    entries.push_back(parse_graph_list_line("I,g1,a5,m1", 10));
    try {
        graph_list_identities(entries, fp);
        ADD_FAILURE() << "accepted one index_fp for two indexes";
    } catch (const std::invalid_argument &ex) {
        EXPECT_NE(std::string::npos, std::string(ex.what()).find("one index_fp")) << ex.what();
        EXPECT_NE(std::string::npos, std::string(ex.what()).find("lines 1 and 10")) << ex.what();
    }
}


// Review of pass 5, findings 2 and 3: one inventory of what the loaders open. Its table of
// annotation types is checked against the loader itself — for every annotation type the CLI
// knows, the extension of the object initialize_annotation constructs is in the table, maps
// back to that type, and the table says it reads the row-diff anchors beside the graph exactly
// when its matrix is RowDiff<ColumnMajor>, and the sequence headers exactly when it is a
// MultiIntMatrix (build_annotated_dbg's and load_coord_to_header's own tests)
TEST(GraphletServer, InventoryTableMatchesTheLoadersTypes) {
    using namespace mtg::annot;
    size_t seen = 0;
    for (int t = Config::ColumnCompressed; t <= Config::RowDiffDiskCoord; ++t) {
        const auto type = static_cast<Config::AnnotationType>(t);
        auto annotation = initialize_annotation(type);
        ASSERT_TRUE(annotation) << t;
        const std::string ext = annotation->file_extension();
        EXPECT_EQ(type, parse_annotation_type("x" + ext)) << ext;
        const auto &matrix = annotation->get_matrix();
        const bool coordinates = dynamic_cast<const matrix::MultiIntMatrix *>(&matrix);
        const bool anchors = dynamic_cast<const matrix::RowDiff<matrix::ColumnMajor> *>(&matrix);
        const IndexAnnotationKind *kind = nullptr;
        for (const IndexAnnotationKind &k : index_annotation_kinds()) {
            // the first suffix that matches decides, as in parse_annotation_type
            if (ext.size() >= k.extension.size()
                    && ext.compare(ext.size() - k.extension.size(), k.extension.size(),
                                   k.extension) == 0) {
                kind = &k;
                break;
            }
        }
        ASSERT_TRUE(kind) << ext;
        EXPECT_EQ(ext, kind->extension);
        EXPECT_EQ(coordinates, kind->coordinates) << ext;
        EXPECT_EQ(anchors, kind->row_diff_anchors) << ext;
        ++seen;
    }
    EXPECT_EQ(index_annotation_kinds().size(), seen);
}

// Finding 2 itself: two lines whose graph and annotation are symlinks to the same files, but
// whose .seqs beside the symlinks differ, are two indexes — each validated against the manifest
// it names, so the shared manifest (A's) refuses B, naming B's own .seqs; with B's own
// manifest they state different index_fp
TEST(GraphletServer, SymlinkedMainFilesDoNotHideTheirSidecars) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::absolute("temp_inventory_symlinks_" + std::to_string(getpid()));
    fs::remove_all(dir);
    fs::create_directories(dir / "A");
    fs::create_directories(dir / "B");
    auto write = [&](const fs::path &p, const std::string &content) {
        std::ofstream(p, std::ios::binary) << content;
        return p.string();
    };
    write(dir / "shared.dbg", "graph bytes");
    write(dir / "sharedanno.column_coord.annodbg", "annotation bytes");
    for (const char *side : { "A", "B" }) {
        fs::create_symlink(dir / "shared.dbg", dir / side / "graph.dbg");
        fs::create_symlink(dir / "sharedanno.column_coord.annodbg",
                           dir / side / "annotation.column_coord.annodbg");
    }
    write(dir / "A" / "annotation.seqs", "sampleA");
    write(dir / "B" / "annotation.seqs", "sampleBBBB");
    const std::string ga = (dir / "A" / "graph.dbg").string(),
                      aa = (dir / "A" / "annotation.column_coord.annodbg").string(),
                      gb = (dir / "B" / "graph.dbg").string(),
                      ab = (dir / "B" / "annotation.column_coord.annodbg").string();
    // the sidecars are looked up beside the listed spelling, not beside the symlinks' target
    EXPECT_EQ((std::vector<std::string> { ga, aa, (dir / "A" / "annotation.seqs").string() }),
              index_bundle_files(ga, aa));
    EXPECT_EQ((std::vector<std::string> { gb, ab, (dir / "B" / "annotation.seqs").string() }),
              index_bundle_files(gb, ab));
    const Json::Value inventory = index_inventory_json(gb, ab);
    EXPECT_EQ("coord_to_header", inventory["files"][2]["role"].asString());
    EXPECT_TRUE(inventory["files"][2]["exists"].asBool());
    EXPECT_FALSE(inventory["files"][2]["required"].asBool());
    auto manifest = [&](const std::string &name, const std::string &seqs) {
        Json::Value m;
        for (const auto &[path, content] : std::vector<std::pair<std::string, std::string>> {
                 { "graph.dbg", "graph bytes" },
                 { "annotation.column_coord.annodbg", "annotation bytes" },
                 { "annotation.seqs", seqs } }) {
            Json::Value e;
            e["path"] = path;
            e["size"] = Json::UInt64(content.size());
            e["sha256"] = sha256_hex(content);
            m["files"].append(e);
        }
        return write(dir / name, json_text(m, true));
    };
    const std::string shared = manifest("shared-manifest.json", "sampleA");
    std::vector<GraphListEntry> entries {
        parse_graph_list_line("A," + ga + "," + aa + "," + shared + ",bundle", 1),
        parse_graph_list_line("B," + gb + "," + ab + "," + shared + ",bundle", 2),
    };
    auto fp = [&](const GraphListEntry &e) {
        return index_manifest_fingerprint(e.manifest_path,
                                          index_bundle_files(e.graph_path, e.annotation_path));
    };
    try {
        graph_list_identities(entries, fp);
        ADD_FAILURE() << "B's .seqs was hidden behind A's";
    } catch (const std::exception &ex) {
        EXPECT_NE(std::string::npos,
                  std::string(ex.what()).find((dir / "B" / "annotation.seqs").string()))
                << ex.what();
    }
    // each with its own manifest: two indexes, two identities
    entries[1] = parse_graph_list_line("B," + gb + "," + ab + ","
                                       + manifest("b-manifest.json", "sampleBBBB") + ",bundle", 2);
    const auto ids = graph_list_identities(entries, fp);
    EXPECT_NE(ids.at({ ga, aa }).second, ids.at({ gb, ab }).second);
    EXPECT_FALSE(ids.at({ ga, aa }).second.empty());
    // the same pair spelled twice (here through the same symlinks) is still one index
    entries.push_back(parse_graph_list_line("A2," + ga + "," + aa, 3));
    EXPECT_EQ(ids.at({ ga, aa }), graph_list_identities(entries, fp).at({ ga, aa }));
    // without the .seqs beside B, B loads no headers: a third bundle, refused by A's manifest
    // only if it lists... (it does: the manifest lists a file B does not load, which is
    // allowed; the loaded files are covered) — so B states A's identity only when it loads A's
    // files: here its sidecar set differs, and one index_fp for two indexes is refused
    fs::remove(dir / "B" / "annotation.seqs");
    entries.pop_back();
    entries[1] = parse_graph_list_line("B," + gb + "," + ab + "," + shared + ",bundle", 2);
    try {
        graph_list_identities(entries, fp);
        ADD_FAILURE() << "two indexes stated one index_fp";
    } catch (const std::invalid_argument &ex) {
        EXPECT_NE(std::string::npos, std::string(ex.what()).find("one index_fp")) << ex.what();
    }
    fs::remove_all(dir);
}
