#include "gtest/gtest.h"

#include <atomic>
#include <chrono>
#include <map>
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

#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"


namespace {

using namespace mtg::cli;
using mtg::graph::traversal::AttemptAborted;
using mtg::graph::traversal::ExternalStop;
using Clock = Attempt::Clock;

// A clock the test moves: the registry's retention and an attempt's bound read it
struct FakeClock {
    Clock::time_point base = Clock::now();
    std::atomic<int64_t> ms { 0 };
    std::function<Clock::time_point()> fn() {
        return [this]() { return base + std::chrono::milliseconds(ms.load()); };
    }
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
    if (clock)
        s.clock = clock->fn();
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
    for (const char *f : { "attempt_id", "budget_id", "locus_id", "not_after_ms" }) {
        fields.append(f);
    }
    EXPECT_EQ(fields, att["fields"]);
    EXPECT_NE(std::string::npos, att["not_after"].asString().find("clock_skew_allowance_ms"));
    EXPECT_EQ(registry.server_instance(), att["server_instance"].asString());
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
    clock.ms = 30'000;
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
    // by age only: after retention_s they go, and the ids are free
    clock.ms += 30'000;
    EXPECT_FALSE(registry.start(late));
    EXPECT_EQ(404, registry.cancel("Z4", 0).first);
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
