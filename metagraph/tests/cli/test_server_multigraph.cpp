#include "gtest/gtest.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <json/json.h>
#include <sdsl/int_vector.hpp>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "cli/config/config.hpp"
#include "cli/pattern.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"


// The multi-graph server's parts without the HTTP library (server_checks.hpp): the graph a
// traversal request selects (graph, graphs), the names a /pattern request selects, the
// reservations and the plan of `in_ram`, the columns' overlap; the sampled check of a mask at
// load (load_annotated_graph.hpp); a seed a graph does not hold as a per-seed result
// (process_traverse_request); an attempt's bound after a load (traverse_attempts.hpp); and
// the record mapping an `in_ram` load's copy shares with the resident pair
// (initialize_annotated_dbg with AnnotatedDBG::share_coord_to_header)

namespace {

using namespace mtg;
using namespace mtg::cli;
using mtg::graph::DBGSuccinct;
using mtg::graph::DeBruijnGraph;
using mtg::graph::AnnotatedDBG;
using mtg::graph::boss::BOSS;
namespace fs = std::filesystem;

Json::Value parse(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errors;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errors)) << errors;
    return v;
}

// the message of |f|'s std::invalid_argument, "" when it throws none
template <class F>
std::string refusal_of(const F &f) {
    try {
        f();
    } catch (const std::invalid_argument &e) {
        return e.what();
    }
    return "";
}

std::string random_seq(size_t length, std::mt19937 *rng) {
    std::string s(length, 'A');
    for (char &c : s) {
        c = "ACGT"[(*rng)() % 4];
    }
    return s;
}

} // namespace


// ---------------------------------------------------------------- the graph a request selects

// SPEC-T §10.3, routing: on a single-graph server graph / graph_path and graphs are refused
// (the first with the message it always had); on a multi-graph server graph or graphs ([name],
// /search's form) names the graph, never both, and both spell one selection
TEST(MultiGraphSelection, GraphOrGraphsNamesOneGraph) {
    const std::string single_graph = "Bad request: this server hosts a single graph; remove "
                                     "the 'graph' / 'graph_path' field";
    EXPECT_EQ(single_graph, refusal_of([] {
        traverse_graph_selection(parse(R"({"graph": "a"})"), false); }));
    EXPECT_EQ(single_graph, refusal_of([] {
        traverse_graph_selection(parse(R"({"graph_path": "x", "graphs": ["a"]})"), false); }));
    EXPECT_EQ("Bad request: this server hosts a single graph; remove the 'graphs' field",
              refusal_of([] { traverse_graph_selection(parse(R"({"graphs": ["a"]})"), false); }));
    EXPECT_EQ("", refusal_of([] {
        traverse_graph_selection(parse(R"({"seeds": [], "in_ram": true})"), false); }));
    EXPECT_TRUE(traverse_graph_selection(parse("{}"), false).name.empty());

    // multi-graph: graph as it always was
    GraphSelection s = traverse_graph_selection(parse(R"({"graph": "a"})"), true);
    EXPECT_EQ("a", s.name);
    EXPECT_FALSE(s.graph_path);
    s = traverse_graph_selection(parse(R"({"graph": "a", "graph_path": "p.dbg"})"), true);
    EXPECT_EQ("p.dbg", s.graph_path.value_or(""));
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([] { traverse_graph_selection(parse("{}"), true); }));
    EXPECT_EQ("Bad request: 'graph' (index name) is required in multi-graph mode",
              refusal_of([] { traverse_graph_selection(parse(R"({"graph": 1})"), true); }));
    // a graph_path that is no string is left to the request's parser (it names the field)
    EXPECT_FALSE(traverse_graph_selection(parse(R"({"graph": "a", "graph_path": 1})"), true)
                         .graph_path);

    // graphs: one name, as /search takes it — the selection graph makes, under another spelling
    s = traverse_graph_selection(parse(R"({"graphs": ["b"], "graph_path": "q.dbg"})"), true);
    EXPECT_EQ("b", s.name);
    EXPECT_EQ("q.dbg", s.graph_path.value_or(""));
    const GraphSelection spelled_graph
            = traverse_graph_selection(parse(R"({"graph": "b", "graph_path": "q.dbg"})"), true);
    EXPECT_EQ(spelled_graph.name, s.name);
    EXPECT_EQ(spelled_graph.graph_path, s.graph_path);
    const std::string one = "Bad request: 'graphs' names the one graph a traversal reads: "
                            "expected [name]";
    for (const char *body : { R"({"graphs": []})", R"({"graphs": ["a", "b"]})",
                              R"({"graphs": "a"})", R"({"graphs": [1]})",
                              R"({"graphs": null})" }) {
        EXPECT_EQ(one, refusal_of([&] { traverse_graph_selection(parse(body), true); })) << body;
    }
    EXPECT_EQ("Bad request: give the graph as 'graph' or as 'graphs', not both",
              refusal_of([] {
                  traverse_graph_selection(parse(R"({"graph": "a", "graphs": ["a"]})"), true);
              }));
}

// The parsers accept the server's fields (strictly typed), so that the unknown-field check
// passes them
TEST(MultiGraphSelection, ParsersDeclareGraphsAndInRam) {
    Json::Value t = parse(R"({"seeds": [{"sequence": "ACGTACGTACGT"}], "strategy": {},
                              "graphs": ["a"], "in_ram": true})");
    EXPECT_NO_THROW(parse_traverse_request(t));
    t["in_ram"] = "yes";
    EXPECT_EQ("request.in_ram: expected a boolean",
              refusal_of([&] { parse_traverse_request(t); }));
    t["in_ram"] = false;
    t["graphs"] = 5;
    EXPECT_EQ("request.graphs: expected an array of strings",
              refusal_of([&] { parse_traverse_request(t); }));

    Json::Value r = parse(R"({"sequence": "ACGTACGTACGT", "labels": ["x"], "graphs": ["a"],
                              "in_ram": false})");
    EXPECT_NO_THROW(parse_resolve_request(r));
    r["in_ram"] = 1;
    EXPECT_EQ("request.in_ram: expected a boolean",
              refusal_of([&] { parse_resolve_request(r); }));

    EXPECT_FALSE(in_ram_field(parse("{}")));
    EXPECT_TRUE(in_ram_field(parse(R"({"in_ram": true})")).value());
    EXPECT_FALSE(in_ram_field(parse(R"({"in_ram": false})")).value());
    EXPECT_EQ("request.in_ram: expected a boolean",
              refusal_of([] { in_ram_field(parse(R"({"in_ram": "true"})")); }));
}

// POST /pattern on a multi-graph server: /search's selection (deduplicated, sorted; every name
// of a small server without the field), every name known
TEST(MultiGraphSelection, PatternGraphNames) {
    const std::vector<std::string> known = { "c", "a", "b" };
    EXPECT_EQ(std::vector<std::string>({ "a", "c" }),
              pattern_graph_names(parse(R"({"graphs": ["c", "a", "c"]})"), known));
    EXPECT_EQ(std::vector<std::string>({ "a", "b", "c" }), pattern_graph_names(parse("{}"), known));
    EXPECT_EQ("request.graphs: unknown graph 'd' (GET /capabilities lists the graphs)",
              refusal_of([&] { pattern_graph_names(parse(R"({"graphs": ["a", "d"]})"), known); }));
    for (const char *body : { R"({"graphs": []})", R"({"graphs": "a"})", R"({"graphs": [1]})",
                              R"({"graphs": null})" }) {
        EXPECT_EQ("request.graphs: expected a non-empty array of graph names (GET "
                  "/capabilities lists them)",
                  refusal_of([&] { pattern_graph_names(parse(body), known); })) << body;
    }
    // more names than a request may leave out (/search's 10)
    std::vector<std::string> many;
    for (int i = 0; i < 11; ++i) {
        many.push_back("g" + std::to_string(i));
    }
    EXPECT_EQ(std::vector<std::string>({ "g3" }),
              pattern_graph_names(parse(R"({"graphs": ["g3"]})"), many));
    EXPECT_EQ("request.graphs: required on this server, which hosts 11 graph names (more than "
              "10; GET /capabilities lists them)",
              refusal_of([&] { pattern_graph_names(parse("{}"), many); }));
    EXPECT_EQ(many.size(), pattern_graph_names(parse("{}"), many, 11).size());
}


// ---------------------------------------------------------------- in_ram

// /search's rule: a load only for in_ram, on a server on mmap, of a pair within the capacity
TEST(MultiGraphInRam, PlanIsTheRuleOfSearch) {
    EXPECT_EQ(InRamPlan::RESIDENT, in_ram_plan(false, true, 10, 100));
    EXPECT_EQ(InRamPlan::RESIDENT, in_ram_plan(false, false, 10, 100));
    EXPECT_EQ(InRamPlan::RESIDENT_IN_RAM, in_ram_plan(true, false, 10, 100));
    EXPECT_EQ(InRamPlan::RESIDENT_TOO_LARGE, in_ram_plan(true, true, 101, 100));
    EXPECT_EQ(InRamPlan::RESIDENT_TOO_LARGE, in_ram_plan(true, true, 1, 0));
    EXPECT_EQ(InRamPlan::LOAD, in_ram_plan(true, true, 100, 100));
    EXPECT_EQ(InRamPlan::LOAD, in_ram_plan(true, true, 0, 0));
}

// A reservation waits until its bytes are free (as /search's wait), and a waiter whose client
// left gives up without taking anything
TEST(MultiGraphInRam, ReservationsWaitForTheMemory) {
    LoadReservations r(100);
    ASSERT_TRUE(r.reserve(60));
    EXPECT_EQ(40u, r.left());
    std::atomic<bool> second { false };
    std::thread waiter([&]() {
        EXPECT_TRUE(r.reserve(70));
        second = true;
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    EXPECT_FALSE(second.load());
    r.release(60);
    waiter.join();
    EXPECT_TRUE(second.load());
    EXPECT_EQ(30u, r.left());

    // a waiter that is gone: nothing taken, the memory stays as it was
    std::atomic<int> asked { 0 };
    EXPECT_FALSE(r.reserve(50, [&]() { return ++asked >= 3; }, 1));
    EXPECT_GE(asked.load(), 3);
    EXPECT_EQ(30u, r.left());
    // one that fits is not asked
    asked = 0;
    EXPECT_TRUE(r.reserve(30, [&]() { ++asked; return true; }));
    EXPECT_EQ(0, asked.load());
    EXPECT_EQ(0u, r.left());
    r.release(100);
    EXPECT_EQ(100u, r.left());
}


// ---------------------------------------------------------------- columns_disjoint

namespace {

// the oracle: every name with the set of pairs holding it
ColumnOverlap naive_overlap(const std::vector<std::vector<std::string>> &columns) {
    std::map<std::string, std::set<size_t>> pairs_of;
    ColumnOverlap result;
    for (size_t p = 0; p < columns.size(); ++p) {
        for (const std::string &name : columns[p]) {
            pairs_of[name].insert(p);
            result.columns++;
        }
    }
    for (const auto &[name, pairs] : pairs_of) {
        if (pairs.size() > 1) {
            if (!result.shared)
                result.example = name;
            result.shared++;
        }
    }
    result.disjoint = !result.shared;
    return result;
}

} // namespace

// Disjoint exactly when no name is a column of two pairs; the shared names counted once each,
// the smallest named; with a weak hash (every name of one length in one group) the names are
// still compared themselves
TEST(MultiGraphColumns, OverlapIsTheNaiveOne) {
    std::mt19937 rng(20261009);
    auto weak = [](const std::string &s) { return s.size(); };
    for (int round = 0; round < 300; ++round) {
        std::vector<std::vector<std::string>> columns(1 + rng() % 5);
        const size_t alphabet = 2 + rng() % 30;
        for (auto &names : columns) {
            std::set<std::string> distinct;
            const size_t n = rng() % 8;
            for (size_t i = 0; i < n; ++i) {
                distinct.insert("c" + std::to_string(rng() % alphabet));
            }
            names.assign(distinct.begin(), distinct.end());
            std::shuffle(names.begin(), names.end(), rng);
        }
        const ColumnOverlap expected = naive_overlap(columns);
        for (const ColumnOverlap &got : { column_overlap(columns), column_overlap(columns, weak) }) {
            EXPECT_EQ(expected.disjoint, got.disjoint) << round;
            EXPECT_EQ(expected.shared, got.shared) << round;
            EXPECT_EQ(expected.columns, got.columns) << round;
            EXPECT_EQ(expected.example, got.example) << round;
        }
    }
    // the chunks of an index partition its samples: disjoint
    EXPECT_TRUE(column_overlap({ { "s1", "s2" }, { "s3" }, {} }).disjoint);
    // one name in three pairs is one shared name
    const ColumnOverlap three = column_overlap({ { "x", "a" }, { "x" }, { "x", "b" } });
    EXPECT_FALSE(three.disjoint);
    EXPECT_EQ(1u, three.shared);
    EXPECT_EQ("x", three.example);
}


// ---------------------------------------------------------------- the mask's sampled check

namespace {

// A sparse succinct graph (k = 12) with many sink dummies (one per record end), written with its
// mask to |base|.dbg and .edgemask
std::string masked_graph(const std::string &base) {
    std::mt19937 rng(7);
    DBGSuccinct graph(12);
    for (int r = 0; r < 240; ++r) {
        graph.add_sequence(random_seq(26, &rng));
    }
    graph.mask_dummy_kmers(1, false);
    // as servers load graphs: the succinct state, its mask a bit_vector_small
    graph.switch_state(BOSS::State::STAT);
    graph.serialize(base);
    return base;
}

std::shared_ptr<DBGSuccinct> load_graph(const std::string &base) {
    auto graph = std::make_shared<DBGSuccinct>(2);
    EXPECT_TRUE(graph->load(base + ".dbg"));
    return graph;
}

// the edges with W = $ by a scan of W (the oracle: no rank, no select), plain then marked
std::pair<std::vector<uint64_t>, std::vector<uint64_t>> sentinel_edges(const DBGSuccinct &g) {
    const BOSS &boss = g.get_boss();
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> edges;
    for (uint64_t i = 1; i <= boss.num_edges(); ++i) {
        const auto w = boss.get_W(i);
        if (w == BOSS::kSentinelCode)
            edges.first.push_back(i);
        else if (w == BOSS::kSentinelCode + boss.alph_size)
            edges.second.push_back(i);
    }
    return edges;
}

// |base|'s mask with every edge of |valid| marked valid, written to |out| (graph and mask)
std::string with_valid(const std::string &base, const std::vector<uint64_t> &valid,
                       const std::string &out) {
    fs::copy_file(base + ".dbg", out + ".dbg", fs::copy_options::overwrite_existing);
    auto graph = load_graph(base);
    const bit_vector &mask = *graph->get_mask();
    sdsl::bit_vector bits(mask.size(), false);
    for (uint64_t i = 0; i < mask.size(); ++i) {
        bits[i] = mask[i];
    }
    for (uint64_t i : valid) {
        bits[i] = true;
    }
    std::ofstream file(out + ".edgemask", std::ios::binary);
    bit_vector_small(std::move(bits)).serialize(file);
    return out;
}

} // namespace

// The sample looks up the main dummy edge and |samples| edges with W = $ drawn from both
// forms, the same ones every time; all of them when there are no more (the full check); a
// mask_dummy_kmers mask passes, one that marks such edges valid is found
TEST(MaskCheck, SampleOfTheSentinelEdges) {
    const std::string dir = test_dump_dir() + "/multigraph_mask_check";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string base = masked_graph(dir + "/good");
    auto good = load_graph(base);
    const auto [plain, marked] = sentinel_edges(*good);
    std::vector<uint64_t> sentinels = plain;
    sentinels.insert(sentinels.end(), marked.begin(), marked.end());
    const uint64_t total = sentinels.size();
    ASSERT_GT(total, 200u);
    ASSERT_EQ(static_cast<uint64_t>(BOSS::kSentinelCode),
              static_cast<uint64_t>(good->get_boss().get_W(1))) << "the main dummy edge";

    // the good mask: nothing valid, by the sample, the full check and the oracle
    MaskSentinelSample s = sample_mask_sentinels(*good, 16);
    EXPECT_EQ(total, s.sentinel_edges);
    EXPECT_EQ(17u, s.looked_up);
    EXPECT_EQ(0u, s.valid);
    EXPECT_EQ(0u, good->count_valid_sentinel_edges());
    s = sample_mask_sentinels(*good, total);
    EXPECT_EQ(total, s.looked_up);
    EXPECT_EQ(0u, s.valid);
    // the default: a small sample of a graph with more
    s = sample_mask_sentinels(*good);
    EXPECT_EQ(std::min(total, kMaskCheckSamples + 1), s.looked_up);
    EXPECT_EQ(0u, s.valid);

    // every such edge valid; only the main dummy edge (always looked up); only the last
    // quarter of them (drawn from the whole range): each found by a small sample, the same
    // edges every time, and all of them by a sample as large as the edges (the full count)
    const std::vector<uint64_t> last_quarter(sentinels.end() - total / 4, sentinels.end());
    for (const auto &[name, valid] : { std::make_pair(std::string("all"), sentinels),
                                       std::make_pair(std::string("main"),
                                                      std::vector<uint64_t>{ 1 }),
                                       std::make_pair(std::string("last_quarter"),
                                                      last_quarter) }) {
        auto bad = load_graph(with_valid(base, valid, dir + "/" + name));
        EXPECT_EQ(valid.size(), bad->count_valid_sentinel_edges()) << name;
        const MaskSentinelSample a = sample_mask_sentinels(*bad, 32);
        EXPECT_EQ(33u, a.looked_up) << name;
        EXPECT_GT(a.valid, 0u) << name;
        EXPECT_LE(a.valid, a.looked_up) << name;
        if (name == "all") {
            EXPECT_EQ(a.looked_up, a.valid);
        }
        const MaskSentinelSample b = sample_mask_sentinels(*load_graph(dir + "/" + name), 32);
        EXPECT_EQ(a.valid, b.valid) << name;
        EXPECT_EQ(valid.size(), sample_mask_sentinels(*bad, total).valid) << name;
    }

    // a random half marked valid: the full sample counts exactly them (the oracle), a sample
    // finds about half of what it looks up
    std::vector<uint64_t> half;
    std::mt19937 rng(11);
    for (uint64_t e : sentinels) {
        if (rng() % 2)
            half.push_back(e);
    }
    auto bad = load_graph(with_valid(base, half, dir + "/half"));
    EXPECT_EQ(half.size(), sample_mask_sentinels(*bad, total).valid);
    const MaskSentinelSample sampled = sample_mask_sentinels(*bad, 100);
    EXPECT_EQ(101u, sampled.looked_up);
    EXPECT_GT(sampled.valid, 25u);
    EXPECT_LT(sampled.valid, 76u);

    // at load: found and recorded (the pattern search then answers mask_invalid), or not
    std::shared_ptr<DeBruijnGraph> served = load_graph(dir + "/half");
    EXPECT_GT(check_mask_at_load(served), 0u);
    EXPECT_TRUE(mask_invalid_at_load(*served));
    std::shared_ptr<DeBruijnGraph> fine = load_graph(base);
    EXPECT_EQ(0u, check_mask_at_load(fine));
    EXPECT_FALSE(mask_invalid_at_load(*fine));
    // a graph without a mask is not looked up
    auto bare = std::make_shared<DBGSuccinct>(2);
    ASSERT_TRUE(bare->load_without_mask(base + ".dbg"));
    EXPECT_EQ(0u, sample_mask_sentinels(*bare).looked_up);
}


// ---------------------------------------------------------------- not_in_graph

namespace {

// the k-mers of |records| (a BASIC graph: as written)
std::set<std::string> kmers_of(const std::vector<std::string> &records, size_t k) {
    std::set<std::string> kmers;
    for (const std::string &r : records) {
        for (size_t i = 0; i + k <= r.size(); ++i) {
            kmers.insert(r.substr(i, k));
        }
    }
    return kmers;
}

uint64_t present_in(const std::string &seed, const std::set<std::string> &kmers, size_t k) {
    uint64_t n = 0;
    for (size_t i = 0; i + k <= seed.size(); ++i) {
        n += kmers.count(seed.substr(i, k));
    }
    return n;
}

} // namespace

// A request to a graph that need not hold its seeds (a multi-graph server's chunk, whichever
// spelling named it; TraverseLimits::not_in_graph_per_seed): a seed with a k-mer the graph does
// not have is answered per seed — outcome.walks not_in_graph, its k-mers and those the graph
// has (counted here on the records) — and the other seeds are walked as alone; on the only
// graph of a server (the limits' default) the request fails as it always did
TEST(TraverseNotInGraph, AbsentSeedsAreResults) {
    const size_t k = 11;
    std::mt19937 rng(5);
    const std::vector<std::string> records = { random_seq(80, &rng), random_seq(80, &rng) };
    auto index = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, records, { "A", "B" }, DeBruijnGraph::BASIC);
    const std::set<std::string> kmers = kmers_of(records, k);
    const std::string present = records[0].substr(10, 40);
    const std::string absent = random_seq(30, &rng);
    const std::string partial = records[1].substr(0, 25) + random_seq(25, &rng);
    ASSERT_EQ(0u, present_in(absent, kmers, k));
    const uint64_t partly = present_in(partial, kmers, k);
    ASSERT_GT(partly, 0u);
    ASSERT_LT(partly, partial.size() - k + 1);

    auto request = [&](const std::vector<std::string> &seeds) {
        Json::Value r = parse(R"({"strategy": {"direction": "right", "bounds":
                                  {"max_extension_bp": 20}, "output": {"detail": "full",
                                  "timing": false}}})");
        for (size_t i = 0; i < seeds.size(); ++i) {
            Json::Value s;
            s["sequence"] = seeds[i];
            s["seed_id"] = "s" + std::to_string(i);
            r["seeds"].append(s);
        }
        return r;
    };
    TraverseLimits fan_out;
    fan_out.not_in_graph_per_seed = true;

    // the only graph of a server: the request fails, naming the seed and its graph runs
    try {
        process_traverse_request(request({ present, absent }), *index, "");
        ADD_FAILURE() << "an absent seed was answered";
    } catch (const InvalidRequest &e) {
        EXPECT_EQ("seed 's1': Seed is not fully present in the graph; graph runs: o20",
                  std::string(e.what()));
    }

    const Json::Value out = process_traverse_request(request({ present, absent, partial }),
                                                     *index, "", fan_out);
    ASSERT_EQ(3u, out["results"].size());
    for (size_t i : { 1u, 2u }) {
        const Json::Value &r = out["results"][static_cast<Json::ArrayIndex>(i)];
        const std::string &seed = i == 1 ? absent : partial;
        EXPECT_EQ("not_in_graph", r["outcome"]["walks"].asString()) << i;
        EXPECT_EQ(seed.size() - k + 1, r["not_in_graph"]["kmers"].asUInt64()) << i;
        EXPECT_EQ(present_in(seed, kmers, k), r["not_in_graph"]["kmers_present"].asUInt64()) << i;
        EXPECT_FALSE(r.isMember("arms")) << i;
        EXPECT_EQ("s" + std::to_string(i), r["seed"]["seed_id"].asString());
        EXPECT_EQ(seed.size(), r["seed"]["length_bp"].asUInt64());
        EXPECT_NE(std::string::npos,
                  r["error"].asString().find("Seed is not fully present in the graph")) << i;
        EXPECT_EQ(0u, r["limitations"].size()) << i;
    }
    // the present seed is walked as in a request of its own
    Json::Value alone = process_traverse_request(request({ present }), *index, "", fan_out);
    EXPECT_EQ(alone["results"][0], out["results"][0]);
    EXPECT_EQ("complete", out["results"][0]["outcome"]["walks"].asString());

    // under a memory budget it states memory_bound_soft, as every result there; asked for
    // coordinates, it states it has none, and why
    Json::Value budget = request({ absent });
    budget["strategy"]["bounds"]["max_memory_mb"] = 64;
    const Json::Value budgeted = process_traverse_request(budget, *index, "", fan_out);
    const Json::Value &lims = budgeted["results"][0]["limitations"];
    ASSERT_EQ(1u, lims.size()) << lims;
    EXPECT_EQ("memory_bound_soft", lims[0]["kind"].asString());
    Json::Value coordinates = request({ absent });
    coordinates["strategy"]["output"]["coordinates"] = true;
    const Json::Value with_coordinates = process_traverse_request(coordinates, *index, "",
                                                                  fan_out);
    EXPECT_TRUE(with_coordinates["results"][0]["coordinates"].isNull());
    EXPECT_TRUE(with_coordinates["results"][0].isMember("coordinates_reason"));

    // a graphlet is written for walked seeds only
    Json::Value graphlets = request({ present, absent });
    graphlets["strategy"]["output"]["detail"] = "graphlet";
    const Json::Value g = process_traverse_request(graphlets, *index, "", fan_out);
    EXPECT_TRUE(g["results"][0].isMember("graphlet"));
    EXPECT_FALSE(g["results"][1].isMember("graphlet"));

    // a ledger-managed request: its usage states the seed's outcome
    Json::Value managed = request({ absent, present });
    managed["attempt_id"] = "fan-out-1";
    const Json::Value usage = process_traverse_request(managed, *index, "", fan_out)["usage"];
    EXPECT_EQ("not_in_graph", usage["per_seed"][0]["outcome"].asString()) << usage;
    EXPECT_EQ("complete", usage["per_seed"][1]["outcome"].asString()) << usage;
}

// in_ram: the time before the work (the load) is stated in timing, and the attempt's bound
// starts after it
TEST(TraverseNotInGraph, LoadIsStatedInTiming) {
    std::mt19937 rng(9);
    const std::vector<std::string> records = { random_seq(60, &rng) };
    auto index = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, records, { "A" }, DeBruijnGraph::BASIC);
    Json::Value r = parse(R"({"strategy": {"direction": "right", "bounds":
                              {"max_extension_bp": 10}}})");
    r["seeds"][0]["sequence"] = records[0].substr(0, 30);
    TraverseLimits limits;
    Json::Value out = process_traverse_request(r, *index, "", limits);
    EXPECT_FALSE(out["timing"].isMember("load_ms"));
    limits.load_ms = 1234.5;
    out = process_traverse_request(r, *index, "", limits);
    EXPECT_EQ(1234.5, out["timing"]["load_ms"].asDouble());
    r["attempt_id"] = "loaded-1";
    out = process_traverse_request(r, *index, "", limits);
    EXPECT_EQ(1235u, out["usage"]["bound"]["load_ms"].asUInt64());
}

TEST(GraphletAttempt, LoadComesBeforeTheBound) {
    AttemptSettings settings;
    settings.allowance_ms = 1000;
    settings.hard_cap_ms = 899'000;
    AttemptRegistry registry(settings);
    auto attempt = [&](const std::string &id) {
        auto a = std::make_shared<Attempt>(0, std::chrono::system_clock::now(),
                                           registry.settings(), registry.server_instance(),
                                           true);
        AttemptIds ids;
        ids.attempt_id = id;
        a->set_ids(ids);
        return a;
    };
    auto a = attempt("a");
    a->set_bound(2, 1000, 0, 5000.0);
    EXPECT_EQ(8000, a->bound_ms());
    EXPECT_EQ(5000u, a->usage_json("completed")["bound"]["load_ms"].asUInt64());
    // none without a load: the bound and its statement as they always were
    auto b = attempt("b");
    b->set_bound(2, 1000);
    EXPECT_EQ(3000, b->bound_ms());
    EXPECT_FALSE(b->usage_json("completed")["bound"].isMember("load_ms"));
    // a load of 0 is stated (the request asked for in_ram)
    auto c = attempt("c");
    c->set_bound(2, 1000, 0, 0.0);
    EXPECT_EQ(3000, c->bound_ms());
    EXPECT_EQ(0u, c->usage_json("completed")["bound"]["load_ms"].asUInt64());
    // the transport's cap counts the load too
    auto d = attempt("d");
    d->set_bound(1, 1000, 0, 898'500.0);
    EXPECT_EQ(899'000, d->bound_ms());
    EXPECT_EQ("content_timeout", d->usage_json("completed")["bound"]["capped_by"].asString());
}


// ---------------------------------------------------------------- the mapping of an in_ram load

namespace {

// |value| without its "timing" members at every level: the answers compared apart from timing
Json::Value untimed(const Json::Value &value) {
    if (value.isObject()) {
        Json::Value out(Json::objectValue);
        for (const std::string &key : value.getMemberNames()) {
            if (key != "timing")
                out[key] = untimed(value[key]);
        }
        return out;
    }
    if (value.isArray()) {
        Json::Value out(Json::arrayValue);
        for (const Json::Value &v : value) {
            out.append(untimed(v));
        }
        return out;
    }
    return value;
}

// A pair on disk as a graph list names it: a BASIC succinct graph of two records in one
// column (`annotate --anno-filename`) with the records' coordinates, and its record mapping
// (the .seqs: the headers acc1 and acc2), in a fresh directory
struct PairOnDisk {
    std::string dir;
    std::string graph;
    std::string annotation;
    std::vector<std::string> records;
};

PairOnDisk write_pair(const std::string &name) {
    const size_t k = 11;
    std::mt19937 rng(17);
    PairOnDisk pair;
    pair.dir = test_dump_dir() + "/multigraph_in_ram_" + name;
    fs::remove_all(pair.dir);
    fs::create_directories(pair.dir);
    pair.records = { random_seq(120, &rng), random_seq(90, &rng) };
    const uint64_t n1 = pair.records[0].size() - k + 1;
    const uint64_t n2 = pair.records[1].size() - k + 1;
    auto index = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, pair.records, { "F", "F" }, DeBruijnGraph::BASIC, true, { 0, n1 });
    pair.graph = pair.dir + "/graph" + DBGSuccinct::kExtension;
    dynamic_cast<const DBGSuccinct&>(index->get_graph()).serialize(pair.graph);
    const std::string anno_base = pair.dir + "/anno";
    index->get_annotator().serialize(anno_base);
    pair.annotation = anno_base + annot::ColumnCoordAnnotator::kExtension;
    annot::CoordToHeader cth({ { "acc1", "acc2" } }, { { n1, n2 } });
    cth.serialize(anno_base);
    return pair;
}

// The Config a server loads |pair| with (-i, -a); made from a command line as main() makes
// it (the Config sets the swap path; the one the other tests use is restored)
std::unique_ptr<Config> config_of(const PairOnDisk &pair) {
    const fs::path swap = utils::get_swap_path();
    std::ofstream(pair.dir + "/request.json") << "{}";
    std::vector<std::string> args = { "metagraph", "pattern", "-i", pair.graph, "-a",
                                      pair.annotation, pair.dir + "/request.json" };
    std::vector<char*> argv;
    for (std::string &a : args) {
        argv.push_back(a.data());
    }
    argv.push_back(nullptr);
    auto config = std::make_unique<Config>(static_cast<int>(args.size()), argv.data());
    utils::set_swap_path(swap);
    return config;
}

} // namespace

// The copy an `in_ram` load makes is built with the resident pair's record mapping: the same
// object, whose reverse header index was built once at start-up and is not built again, for
// the traversal routes' copy and for /pattern's (its mask checked, its dummy fraction
// sampled); header-label requests on the copy answer as on the resident. Without the
// resident's mapping the .seqs is loaded, as every load does; a mapping of another
// annotation is not shared
TEST(MultiGraphInRam, LoadSharesTheResidentRecordMapping) {
    const PairOnDisk pair = write_pair("shared");
    ASSERT_TRUE(fs::exists(pair.graph));
    ASSERT_TRUE(fs::exists(pair.annotation));
    ASSERT_TRUE(fs::exists(pair.dir + "/anno" + annot::CoordToHeader::kExtension));
    auto config = config_of(pair);
    PatternPreparation none;
    none.progress = false;
    PatternPreparation for_pattern;
    for_pattern.check_mask = true;
    for_pattern.sample_fraction = true;
    for_pattern.progress = false;

    // the resident pair: its header index built as the server builds it at start-up
    auto resident = initialize_annotated_dbg(*config, none);
    ASSERT_TRUE(resident->get_coord_to_header());
    EXPECT_EQ(0u, resident->get_coord_to_header()->num_header_index_builds());
    EXPECT_EQ(2u, resident->get_coord_to_header()->build_header_index());
    EXPECT_EQ(1u, resident->get_coord_to_header()->num_header_index_builds());

    // the copies: the resident's object, nothing loaded or built for them
    auto copy = initialize_annotated_dbg(*config, none, resident->share_coord_to_header());
    auto pattern_copy = initialize_annotated_dbg(*config, for_pattern,
                                                 resident->share_coord_to_header());
    EXPECT_EQ(resident->get_coord_to_header(), copy->get_coord_to_header());
    EXPECT_EQ(resident->get_coord_to_header(), pattern_copy->get_coord_to_header());
    EXPECT_EQ(resident->share_coord_to_header(), copy->share_coord_to_header());
    EXPECT_EQ(1u, resident->get_coord_to_header()->num_header_index_builds());

    // header-label requests: the copies answer as the resident does
    const std::string &acc1 = pair.records[0];
    Json::Value resolve = parse(R"({"labels": ["acc1", "acc2", "F"]})");
    resolve["sequence"] = acc1.substr(20, 60);
    Json::Value discover = parse(R"({"discover": {"max_labels": 10, "kind": "header"}})");
    discover["sequence"] = acc1.substr(20, 60);
    Json::Value traverse = parse(R"({"seeds": [{"labels": ["acc1"]}], "strategy":
                                    {"direction": "right", "bounds": {"max_extension_bp": 30},
                                     "output": {"detail": "full", "timing": false}}})");
    traverse["seeds"][0]["sequence"] = acc1.substr(0, 40);
    for (const AnnotatedDBG *index : { copy.get(), pattern_copy.get() }) {
        EXPECT_EQ(untimed(process_resolve_request(resolve, *resident, "")),
                  untimed(process_resolve_request(resolve, *index, "")));
        EXPECT_EQ(untimed(process_resolve_request(discover, *resident, "")),
                  untimed(process_resolve_request(discover, *index, "")));
        EXPECT_EQ(untimed(process_traverse_request(traverse, *resident, "")),
                  untimed(process_traverse_request(traverse, *index, "")));
    }
    // (the answers name the headers: the mapping was read)
    const Json::Value answered = process_resolve_request(resolve, *copy, "");
    ASSERT_EQ(3u, answered["labels"].size()) << answered;
    EXPECT_EQ("acc1", answered["labels"][0]["label"].asString());
    EXPECT_EQ("header", answered["labels"][0]["kind"].asString());
    const Json::Value discovered = process_resolve_request(discover, *copy, "");
    ASSERT_EQ(1u, discovered["labels"].size()) << discovered;
    EXPECT_EQ("acc1", discovered["labels"][0]["label"].asString());
    // /pattern on its copy: the record labels of its paths as the resident's (a column
    // annotation's reads are unbudgeted: allowed, as a request on such an index asks)
    const Json::Value pattern = parse_pattern_body(
            "{\"patterns\": [{\"dna\": \"" + acc1.substr(30, 40) + "\"}], "
            "\"long_search\": \"paths\", \"output\": {\"labels\": \"all\"}, "
            "\"allow_unbudgeted_annotation\": true}");
    PatternLimits caps;
    caps.min_information_bits = 0;
    const Json::Value on_resident = process_pattern_request(pattern, *resident, caps, "");
    EXPECT_EQ(untimed(on_resident),
              untimed(process_pattern_request(pattern, *pattern_copy, caps, "")));
    ASSERT_EQ(1u, on_resident["patterns"][0]["results"].size()) << on_resident;
    EXPECT_EQ("record_verified",
              on_resident["patterns"][0]["results"][0]["support"].asString()) << on_resident;
    // one index, the resident's, for all of it
    EXPECT_EQ(1u, resident->get_coord_to_header()->num_header_index_builds());

    // without the resident's mapping the .seqs is loaded: another object with an index of
    // its own, built by its first header lookup
    auto loaded = initialize_annotated_dbg(*config, none);
    ASSERT_TRUE(loaded->get_coord_to_header());
    EXPECT_NE(resident->get_coord_to_header(), loaded->get_coord_to_header());
    EXPECT_EQ(0u, loaded->get_coord_to_header()->num_header_index_builds());
    EXPECT_EQ(untimed(process_resolve_request(resolve, *resident, "")),
              untimed(process_resolve_request(resolve, *loaded, "")));
    EXPECT_EQ(1u, loaded->get_coord_to_header()->num_header_index_builds());
    EXPECT_EQ(1u, resident->get_coord_to_header()->num_header_index_builds());

    // a mapping of another annotation (two columns to this one's label) is not shared
    auto other = std::make_shared<annot::CoordToHeader>(
            std::vector<std::vector<std::string>>{ { "x" }, { "y" } },
            std::vector<std::vector<uint64_t>>{ { 3 }, { 4 } });
    auto refused = initialize_annotated_dbg(*config, none, other);
    ASSERT_TRUE(refused->get_coord_to_header());
    EXPECT_NE(other.get(), refused->get_coord_to_header());
    EXPECT_EQ(2u, refused->get_coord_to_header()->build_header_index());
}
