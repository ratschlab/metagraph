#include "gtest/gtest.h"

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <json/json.h>

#include "cli/config/config.hpp"
#include "cli/server_checks.hpp"
#include "cli/server_envelope.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/file_utils.hpp"
#include "tests/source_path.hpp"


// Handler capacity (E19, SPEC-labeled-traversal-core.md §10.3, "Handler capacity";
// SPEC-multigraph-envelope.md section 1.1): the caps on the requests of /traverse, /resolve
// and /pattern whose handlers run at once (RouteSlots), the busy refusal they answer above
// them (ServerBusy: 503, Retry-After, {code: "busy", error}), the flags and their start-up
// rule, and the capabilities that state them (max_concurrent) against the reference samples
// 10, 12 and 24

namespace {

using namespace mtg::cli;

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errors;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errors)) << errors;
    return v;
}

std::string text(const Json::Value &value) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    return Json::writeString(builder, value);
}

// A reference sample, {route, status, note, request?, body}
Json::Value sample(const std::string &name) {
    // by its absolute path in the source tree: read from whatever directory the test runs in (a
    // shard directory of unit_tests_parallel too)
    const std::string path = source_path("tests/data/multigraph_contract/" + name);
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.good()) << path << " not found";
    const std::string body { std::istreambuf_iterator<char>(in),
                             std::istreambuf_iterator<char>() };
    return parse(body);
}

// ---------------------------------------------------------------- RouteSlots

// Without a cap every request takes a slot, and the running ones are counted
TEST(ServerCapacity, NoCapAdmitsEveryRequest) {
    RouteSlots slots("/traverse", "--traverse-max-concurrent", 0);
    std::vector<RouteSlots::Slot> held;
    for (size_t i = 0; i < 100; ++i) {
        held.push_back(slots.try_acquire());
        EXPECT_TRUE(held.back()) << i;
    }
    EXPECT_EQ(100u, slots.running());
    EXPECT_EQ(0u, slots.max());
    held.clear();
    EXPECT_EQ(0u, slots.running());
}

// At its cap a route admits no further request until a slot is given back, by the slot's
// destructor (the handler returning) or by moving another slot over it
TEST(ServerCapacity, CapRefusesAboveItUntilASlotIsFree) {
    RouteSlots slots("/resolve", "--resolve-max-concurrent", 2);
    RouteSlots::Slot a = slots.try_acquire();
    RouteSlots::Slot b = slots.try_acquire();
    ASSERT_TRUE(a);
    ASSERT_TRUE(b);
    EXPECT_EQ(2u, slots.running());
    RouteSlots::Slot c = slots.try_acquire();
    EXPECT_FALSE(c);
    // a refusal holds nothing
    EXPECT_EQ(2u, slots.running());
    {
        // a slot moved keeps its place: one release, by its last owner
        RouteSlots::Slot moved = std::move(a);
        EXPECT_FALSE(a);
        EXPECT_TRUE(moved);
        EXPECT_EQ(2u, slots.running());
        EXPECT_FALSE(slots.try_acquire());
    }
    EXPECT_EQ(1u, slots.running());
    c = slots.try_acquire();
    EXPECT_TRUE(c);
    EXPECT_EQ(2u, slots.running());
    // assigned over, a slot gives its place back first
    b = RouteSlots::Slot();
    EXPECT_EQ(1u, slots.running());
    b = slots.try_acquire();
    EXPECT_TRUE(b);
    EXPECT_FALSE(slots.try_acquire());
    // a slot moved onto itself keeps it
    RouteSlots::Slot &same = b;
    b = std::move(same);
    EXPECT_TRUE(b);
    EXPECT_EQ(2u, slots.running());
}

// Handlers on many threads never hold more slots than the cap, and every slot taken is given
// back: the count ends at 0
TEST(ServerCapacity, ConcurrentHandlersNeverExceedTheCap) {
    const uint64_t cap = 3;
    RouteSlots slots("/pattern", "--pattern-max-concurrent", cap);
    std::atomic<uint64_t> holding { 0 };
    std::atomic<uint64_t> most { 0 };
    std::atomic<uint64_t> admitted { 0 };
    std::atomic<uint64_t> refused { 0 };
    std::vector<std::thread> threads;
    for (size_t t = 0; t < 8; ++t) {
        threads.emplace_back([&]() {
            for (size_t i = 0; i < 2000; ++i) {
                RouteSlots::Slot slot = slots.try_acquire();
                if (!slot) {
                    ++refused;
                    continue;
                }
                ++admitted;
                const uint64_t now = ++holding;
                uint64_t seen = most.load();
                while (now > seen && !most.compare_exchange_weak(seen, now)) {}
                std::this_thread::yield();
                --holding;
            }
        });
    }
    for (std::thread &t : threads) {
        t.join();
    }
    EXPECT_LE(most.load(), cap);
    EXPECT_EQ(8u * 2000, admitted.load() + refused.load());
    EXPECT_GT(admitted.load(), 0u);
    EXPECT_EQ(0u, slots.running());
}

// ---------------------------------------------------------------- the busy refusal

// The message names the route, its cap and its flag, and says that nothing of the request was
// run, parsed or registered
TEST(ServerCapacity, BusyMessage) {
    const RouteSlots slots("/traverse", "--traverse-max-concurrent", 4);
    const std::string message = slots.busy_message();
    EXPECT_EQ("/traverse: this server runs at most 4 /traverse request(s) at once "
              "(--traverse-max-concurrent) and runs that many now: this one was not run, its body "
              "not parsed and nothing of it registered; send it again after Retry-After (1 s)",
              message);
    EXPECT_EQ(1, kBusyRetryAfterS);
}

// answer_request answers a ServerBusy with 503, Retry-After and {code: "busy", error,
// server_instance} (busy_answer, which the handlers write themselves when they refuse at their
// start): never compressed (as every error), and not written to a client that is gone (a route's
// |gone|)
TEST(ServerCapacity, BusyIsAnswered503WithRetryAfter) {
    const RouteSlots slots("/resolve", "--resolve-max-concurrent", 1);
    const auto busy = [&](const std::string &) -> Json::Value {
        throw ServerBusy(slots.busy_message(), "7c1e0f3a9b24d6e8");
    };
    const RequestAnswer direct = busy_answer(ServerBusy(slots.busy_message(), "7c1e0f3a9b24d6e8"));
    EXPECT_EQ(503, direct.status);
    EXPECT_EQ(R"({"code":"busy","error":")" + slots.busy_message()
                      + R"(","server_instance":"7c1e0f3a9b24d6e8"})", direct.body);
    ASSERT_EQ(1u, direct.header.size());
    EXPECT_EQ(std::make_pair(std::string("Retry-After"), std::string("1")), direct.header[0]);
    for (const char *encoding : { "", "gzip", "deflate" }) {
        for (bool compact : { true, false }) {
            const RequestAnswer answer = answer_request("{}", encoding, 7, busy, compact,
                                                        nullptr);
            EXPECT_EQ(503, answer.status) << encoding;
            ASSERT_EQ(1u, answer.header.size()) << encoding;
            EXPECT_EQ("Retry-After", answer.header[0].first);
            EXPECT_EQ("1", answer.header[0].second);
            const Json::Value body = parse(answer.body);
            EXPECT_EQ(3u, body.size()) << answer.body;
            EXPECT_EQ("busy", body["code"].asString());
            EXPECT_EQ(slots.busy_message(), body["error"].asString());
            EXPECT_EQ("7c1e0f3a9b24d6e8", body["server_instance"].asString());
            EXPECT_FALSE(body.isMember("usage"));
            // the members in byte order, as every answer of the server; compact on the
            // traversal routes
            EXPECT_LT(answer.body.find("\"code\""), answer.body.find("\"error\""));
            if (compact) {
                EXPECT_EQ(direct.body, answer.body);
            } else {
                EXPECT_EQ(parse(direct.body), body);
            }
        }
    }
    // a route whose client is gone (or a server that stops) writes nothing
    ResponseControl control;
    control.gone = []() { return true; };
    EXPECT_EQ(0, answer_request("{}", "gzip", 7, busy, true, &control).status);
    control.gone = []() { return false; };
    EXPECT_EQ(503, answer_request("{}", "gzip", 7, busy, true, &control).status);
}

// The loading 503 keeps its shape (no code, Retry-After 60): the busy one is told apart from it
// by its code and Retry-After, from the 503s with usage or code "deadline" by having neither
// usage nor that code
TEST(ServerCapacity, BusyIsNotTheLoading503) {
    const RequestAnswer loading = answer_request("{}", "", 1, [](const std::string &) -> Json::Value {
        throw CurrentlyInitializingError();
    }, true, nullptr);
    EXPECT_EQ(503, loading.status);
    ASSERT_EQ(1u, loading.header.size());
    EXPECT_EQ("60", loading.header[0].second);
    EXPECT_FALSE(parse(loading.body).isMember("code"));
}

// ---------------------------------------------------------------- the capabilities

// max_concurrent: every route's cap, an integer, 0 without one
TEST(ServerCapacity, MaxConcurrentJson) {
    EXPECT_EQ(R"({"/pattern":0,"/resolve":0,"/traverse":0})", text(max_concurrent_json({})));
    RouteCaps caps;
    caps.traverse = 4;
    caps.resolve = 2;
    caps.pattern = 1;
    const Json::Value m = max_concurrent_json(caps);
    EXPECT_EQ(R"({"/pattern":1,"/resolve":2,"/traverse":4})", text(m));
    for (const std::string &route : m.getMemberNames()) {
        EXPECT_TRUE(m[route].isUInt64()) << route;
    }
    // the multigraph block states the same object
    EXPECT_EQ(text(m),
              text(multigraph_capabilities_json(0, caps, Json::Value())["max_concurrent"]));
    EXPECT_EQ(text(max_concurrent_json({})),
              text(multigraph_capabilities_json(0, {}, Json::Value())["max_concurrent"]));
}

// Sample 12's multigraph block and every /traverse probe of sample 10 state max_concurrent of a
// server without caps (as the recorded servers ran)
TEST(ServerCapacity, Samples10And12StateTheCaps) {
    const Json::Value s12 = sample("12_capabilities_multigraph.json")["body"];
    EXPECT_EQ(text(max_concurrent_json({})), text(s12["multigraph"]["max_concurrent"]));
    const Json::Value s10 = sample("10_traverse_capabilities_envelope.json")["body"];
    ASSERT_EQ(2u, s10["answers"].size());
    for (const Json::Value &entry : s10["answers"]) {
        EXPECT_EQ(text(max_concurrent_json({})), text(entry["max_concurrent"]))
            << entry["graph"].asString();
    }
}

// Sample 24: the busy refusal of a /traverse at its cap, unwrapped: 503, {code, error,
// server_instance} and nothing else (no usage, no envelope members), its message the server's for
// the cap it ran with; its request is a managed one (attempt_id, not_after_ms), which nothing
// registered
TEST(ServerCapacity, Sample24IsTheBusyRefusal) {
    const Json::Value s = sample("24_traverse_503_busy.json");
    EXPECT_EQ("POST /traverse", s["route"].asString());
    EXPECT_EQ(503, s["status"].asInt());
    const Json::Value &body = s["body"];
    EXPECT_EQ(parse(R"(["code", "error", "server_instance"])"), Json::Value([&]() {
        Json::Value names(Json::arrayValue);
        for (const std::string &n : body.getMemberNames()) {
            names.append(n);
        }
        return names;
    }()));
    EXPECT_EQ("busy", body["code"].asString());
    EXPECT_EQ(RouteSlots("/traverse", "--traverse-max-concurrent", 1).busy_message(),
              body["error"].asString());
    // the sample's server, as every sample's (SERVER_INSTANCE of the generator)
    EXPECT_EQ(sample("12_capabilities_multigraph.json")["body"]["server_instance"].asString(),
              body["server_instance"].asString());
    EXPECT_TRUE(s["request"]["attempt_id"].isString());
    EXPECT_TRUE(s["request"]["not_after_ms"].isIntegral());
    // the note states the header a client reads
    EXPECT_NE(std::string::npos, s["note"].asString().find("Retry-After: 1"));
}

// ---------------------------------------------------------------- the flags

std::unique_ptr<Config> server_config(const std::vector<std::string> &flags) {
    std::vector<std::string> args = { "metagraph", "server_query" };
    args.insert(args.end(), flags.begin(), flags.end());
    for (const char *a : { "-i", "graph.dbg", "-a", "annotation.column.annodbg" }) {
        args.push_back(a);
    }
    std::vector<char*> argv;
    for (std::string &a : args) {
        argv.push_back(a.data());
    }
    argv.push_back(nullptr);
    // the Config sets the swap path and -p the threads: those the other tests use are restored
    const std::filesystem::path swap = utils::get_swap_path();
    const unsigned int threads = get_num_threads();
    auto config = std::make_unique<Config>(static_cast<int>(args.size()), argv.data());
    utils::set_swap_path(swap);
    set_num_threads(threads);
    return config;
}

// --http-threads and the three caps: integers in [0, 4096], 0 by default (today's server: -p
// handler threads, no cap)
TEST(ServerCapacity, FlagsAndTheirRanges) {
    const auto plain = server_config({});
    EXPECT_EQ(0u, plain->http_threads);
    EXPECT_EQ(0u, plain->traverse_max_concurrent);
    EXPECT_EQ(0u, plain->resolve_max_concurrent);
    EXPECT_EQ(0u, plain->pattern_max_concurrent);
    const auto set = server_config({ "--http-threads", "16", "--traverse-max-concurrent", "4",
                                     "--resolve-max-concurrent", "2",
                                     "--pattern-max-concurrent", "3" });
    EXPECT_EQ(16u, set->http_threads);
    EXPECT_EQ(4u, set->traverse_max_concurrent);
    EXPECT_EQ(2u, set->resolve_max_concurrent);
    EXPECT_EQ(3u, set->pattern_max_concurrent);
    EXPECT_EQ(4096u, server_config({ "--http-threads", "4096" })->http_threads);
    for (const char *flag : { "--http-threads", "--traverse-max-concurrent",
                              "--resolve-max-concurrent", "--pattern-max-concurrent" }) {
        for (const char *value : { "-1", "4097", "1.5", "x", "" }) {
            EXPECT_DEATH(server_config({ "--http-threads", "4096", flag, value }),
                         std::string(flag) + " must be an integer in \\[0, 4096\\]")
                << flag << " " << value;
        }
    }
}

// The caps that are set sum to less than the handler threads (--http-threads, else -p's
// value), else the server refuses to start, naming the flags and the threads it needs, and that
// the uncapped long calls need more (the rule is necessary, not sufficient)
TEST(ServerCapacity, CapsMustLeaveAThread) {
    // one thread left: started
    EXPECT_EQ(3u, server_config({ "--http-threads", "4", "--traverse-max-concurrent", "3" })
                          ->traverse_max_concurrent);
    EXPECT_EQ(1u, server_config({ "--http-threads", "4", "--traverse-max-concurrent", "1",
                                  "--resolve-max-concurrent", "1",
                                  "--pattern-max-concurrent", "1" })->pattern_max_concurrent);
    // the threads of -p when --http-threads is not given
    EXPECT_EQ(2u, server_config({ "-p", "3", "--resolve-max-concurrent", "2" })
                          ->resolve_max_concurrent);
    // --http-threads without caps: nothing to check
    EXPECT_EQ(1u, server_config({ "--http-threads", "1" })->http_threads);
    EXPECT_DEATH(server_config({ "--http-threads", "4", "--traverse-max-concurrent", "4" }),
                 "--traverse-max-concurrent 4, --resolve-max-concurrent 0 and "
                 "--pattern-max-concurrent 0 sum to 4, which needs more handler threads: "
                 "--http-threads must be at least 5, got 4; and above the capped sum by the "
                 "uncapped long calls that run at once \\(/search, /align, a route without a "
                 "cap, POST /traverse/cancel with wait_ms\\)");
    EXPECT_DEATH(server_config({ "--http-threads", "4", "--traverse-max-concurrent", "2",
                                 "--resolve-max-concurrent", "1",
                                 "--pattern-max-concurrent", "2" }),
                 "sum to 5, which needs more handler threads: --http-threads must be at least "
                 "6");
    EXPECT_DEATH(server_config({ "-p", "2", "--pattern-max-concurrent", "2" }),
                 "--http-threads must be at least 3, got 2 \\(-p's value: --http-threads not "
                 "given\\)");
    // a single-thread server (the default -p 1) takes no cap
    EXPECT_DEATH(server_config({ "-p", "1", "--traverse-max-concurrent", "1" }),
                 "--http-threads must be at least 2");
}

} // namespace
