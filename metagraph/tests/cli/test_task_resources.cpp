#include <atomic>
#include <chrono>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/json_helpers.hpp"
#include "cli/pattern.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "common/threads/thread_usage.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

// The per-task resource figures every /traverse, /resolve and /pattern answer states of itself
// (timing.cpu_ms, timing.page_faults; SPEC-labeled-traversal-core.md §6.8,
// SPEC-pattern-search.md §8.4): read from per-thread clocks, so another thread's work never
// counts; present on every route, per seed and per pattern entry; page faults on Linux only.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;

const size_t kK = 15;

std::string random_dna(size_t length, uint32_t seed) {
    std::mt19937 rng(seed);
    std::string s(length, 'A');
    for (char &c : s) {
        c = "ACGT"[rng() % 4];
    }
    return s;
}

// three records sharing a core: the walks branch where they part
struct Fixture {
    std::string core = random_dna(400, 1);
    std::vector<std::string> records = { random_dna(300, 2) + core + random_dna(300, 3),
                                         random_dna(300, 4) + core + random_dna(300, 5),
                                         random_dna(300, 6) + core };
    std::unique_ptr<AnnotatedDBG> anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, records, { "A", "B", "C" }, DeBruijnGraph::BASIC);
};

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errs;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errs)) << errs;
    return v;
}

// spins until this thread used at least |ms| of CPU
void burn_cpu_ms(double ms) {
    const common::ThreadUsageMeter meter;
    volatile uint64_t x = 0;
    while (meter.elapsed().cpu_ms() < ms) {
        for (int i = 0; i < 10000; ++i) {
            x = x + i;
        }
    }
}

// the volatile figures of one timing block: cpu_ms always, page_faults exactly on Linux
void expect_figures(const Json::Value &timing, const std::string &where) {
    ASSERT_TRUE(timing.isObject()) << where;
    ASSERT_TRUE(timing.isMember("cpu_ms")) << where;
    EXPECT_TRUE(timing["cpu_ms"].isDouble()) << where;
    EXPECT_GT(timing["cpu_ms"].asDouble(), 0) << where;
#ifdef __linux__
    ASSERT_TRUE(timing.isMember("page_faults")) << where;
    const Json::Value &faults = timing["page_faults"];
    ASSERT_TRUE(faults.isObject()) << where;
    EXPECT_EQ(2u, faults.size()) << where;
    EXPECT_TRUE(faults["minor"].isIntegral()) << where;
    EXPECT_TRUE(faults["major"].isIntegral()) << where;
#else
    // no per-thread rusage: absent, never 0
    EXPECT_FALSE(timing.isMember("page_faults")) << where;
#endif
}


TEST(ThreadUsage, CountsTheCallingThreadOnly) {
    const common::ThreadUsage zero = common::ThreadUsage::zero();
    EXPECT_TRUE(zero.has_cpu);
    EXPECT_EQ(0u, zero.cpu_ns);
#ifdef __linux__
    EXPECT_TRUE(zero.has_faults);
#else
    EXPECT_FALSE(zero.has_faults);
#endif
    EXPECT_EQ(0u, zero.minor_faults);
    EXPECT_EQ(0u, zero.major_faults);

    // this thread's own work counts
    common::ThreadUsageMeter meter;
    burn_cpu_ms(30);
    EXPECT_GE(meter.elapsed().cpu_ms(), 30);

    // another thread's does not: 200 ms burnt there while this one waits for it
    meter.reset();
    std::atomic<double> other_ms { 0 };
    std::thread other([&]() {
        const common::ThreadUsageMeter own;
        burn_cpu_ms(200);
        other_ms = own.elapsed().cpu_ms();
    });
    other.join();
    EXPECT_GE(other_ms.load(), 200);
    EXPECT_LT(meter.elapsed().cpu_ms(), 50);

    // parts add up; a figure one part lacks is dropped
    common::ThreadUsage a = common::ThreadUsage::zero();
    a.cpu_ns = 5;
    a.minor_faults = a.has_faults ? 2 : 0;
    common::ThreadUsage b = a;
    a.add(b);
    EXPECT_EQ(10u, a.cpu_ns);
    EXPECT_EQ(b.has_faults ? 4u : 0u, a.minor_faults);
    common::ThreadUsage none;
    a.add(none);
    EXPECT_FALSE(a.has_cpu);
    EXPECT_FALSE(a.has_faults);

    // written into a timing block and read back
    common::ThreadUsage c = common::ThreadUsage::zero();
    c.cpu_ns = 1'500'000;
    Json::Value timing(Json::objectValue);
    put_thread_usage(&timing, c);
    EXPECT_DOUBLE_EQ(1.5, timing["cpu_ms"].asDouble());
    EXPECT_EQ(c.has_faults, timing.isMember("page_faults"));
    EXPECT_EQ(c.cpu_ns, thread_usage_of(timing).cpu_ns);
    // nothing measured: nothing written
    Json::Value empty(Json::objectValue);
    put_thread_usage(&empty, common::ThreadUsage());
    EXPECT_TRUE(empty.empty());
}

TEST(TaskResources, TraverseStatesCpuPerRequestAndPerSeed) {
    Fixture f;
    for (const char *detail : { "full", "summary", "graphlet" }) {
        Json::Value request = parse(R"({"strategy": {"bounds": {"max_extension_bp": 200}}})");
        request["strategy"]["output"]["detail"] = detail;
        for (size_t i = 0; i < 3; ++i) {
            Json::Value seed;
            seed["sequence"] = f.core.substr(100 + 50 * i, 40);
            request["seeds"].append(seed);
        }
        const Json::Value out = process_traverse_request(request, *f.anno, "");
        expect_figures(out["timing"], std::string(detail) + ": top level");
        ASSERT_EQ(3u, out["results"].size());
        double seeds_ms = 0;
        for (const Json::Value &r : out["results"]) {
            ASSERT_FALSE(r.isMember("error")) << r["error"];
            expect_figures(r["timing"], std::string(detail) + ": a seed");
            seeds_ms += r["timing"]["cpu_ms"].asDouble();
        }
        // the request's span holds its seeds' spans
        EXPECT_GE(out["timing"]["cpu_ms"].asDouble(), seeds_ms) << detail;
        if (std::string(detail) == "graphlet") {
            // nothing of it in the body: the figures are JSON beside it
            for (const Json::Value &r : out["results"]) {
                const std::string body = r["graphlet"].asString();
                EXPECT_EQ(std::string::npos, body.find("cpu_ms"));
                EXPECT_EQ(std::string::npos, body.find("page_faults"));
            }
        }
    }
    // without timing: no figure anywhere (they live only inside the timing blocks)
    Json::Value request = parse(R"({"strategy": {"output": {"timing": false}}})");
    request["seeds"][0]["sequence"] = f.core.substr(100, 40);
    const Json::Value out = process_traverse_request(request, *f.anno, "");
    EXPECT_FALSE(out.isMember("timing"));
    EXPECT_FALSE(out["results"][0].isMember("timing"));
}

TEST(TaskResources, AttemptUsageStatesTheHandlersCpu) {
    Fixture f;
    Json::Value request = parse(R"({"attempt_id": "a-1", "strategy": {}})");
    request["seeds"][0]["sequence"] = f.core.substr(100, 40);
    const Json::Value out = process_traverse_request(request, *f.anno, "");
    ASSERT_TRUE(out.isMember("usage")) << out;
    // the handler's figures since the attempt was created (here, the CLI's local attempt, in
    // the request's processing; on the server at the handler's start)
    expect_figures(out["usage"]["timing"], "usage");

    // built on another thread than the attempt's handler (never done by the server): none
    AttemptSettings settings;
    auto attempt = std::make_shared<Attempt>(1, std::chrono::system_clock::now(), settings,
                                             "0123456789abcdef", true);
    attempt->set_ids({ "a-2", "", "", std::nullopt, {} });
    EXPECT_TRUE(attempt->usage_json("completed").isMember("timing"));
    Json::Value elsewhere;
    std::thread([&]() { elsewhere = attempt->usage_json("completed"); }).join();
    EXPECT_FALSE(elsewhere.isMember("timing"));
}

TEST(TaskResources, ResolveStatesCpu) {
    Fixture f;
    for (const std::string body : { R"({"labels": ["A", "B"]})",
                                    R"({"discover": {"max_labels": 10}, "select": {"max_seeds": 3}})" }) {
        Json::Value request = parse(body);
        request["sequence"] = f.records[0].substr(200, 500);
        const Json::Value out = process_resolve_request(request, *f.anno, "");
        expect_figures(out["timing"], body);
    }
}

TEST(TaskResources, PatternStatesCpuPerRequestAndPerEntry) {
    Fixture f;
    PatternLimits caps;
    caps.min_information_bits = 4;
    for (const std::string mode : { "count", "partial", "all_or_count" }) {
        Json::Value request = parse(R"({"patterns": [{"dna": "ACGTAC"}, {"dna": "GGATCC"}]})");
        request["patterns"].append(Json::Value(Json::objectValue));
        request["patterns"][2]["dna"] = f.core.substr(10, kK + 5);
        request["mode"] = mode;
        if (mode != "count") {
            request["output"]["labels"] = "all";
            request["allow_unbudgeted_annotation"] = true;
        }
        const Json::Value out = process_pattern_request(request, *f.anno, caps, "rel");
        expect_figures(out["timing"], mode + ": top level");
        ASSERT_EQ(3u, out["patterns"].size());
        for (const Json::Value &entry : out["patterns"]) {
            if (entry.isMember("error"))
                continue;
            expect_figures(entry["timing"], mode + ": an entry");
        }
    }
}

} // namespace
