#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// The JSON of POST /pattern and `metagraph pattern` (src/cli/pattern.cpp) on tiny graphs:
// refusals, the answer's shape, rows, the capabilities block, and the deadline paths the
// integration tests cannot reach deterministically (a time stop answered within the
// finalisation reserve; 503 when even the reserve is overrun), driven by an injected clock.
// Expected contexts come from a brute force over the records' k-mers, never from the engine.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;
using Clock = pattern::Deadline::Clock;

const size_t kK = 7;
const std::vector<std::string> kRecords = { "ACGTTGCAACGTAAGGCTTACGATCCA",
                                            "TTGGCCAACGTACGTTTGCA" };

std::unique_ptr<AnnotatedDBG> tiny(DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, kRecords, { "r1", "r2" }, mode);
}

// a floor low enough for the 4-mers of the tiny graphs
PatternLimits limits() {
    PatternLimits l;
    l.min_information_bits = 4;
    return l;
}

std::string rc(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

// {(strand, k-mer, offset)}: the graph contexts of an exact |p| on a BASIC graph of the
// records (every k-mer as deposited)
std::set<std::tuple<std::string, std::string, uint64_t>> oracle(const std::string &p) {
    std::set<std::string> kmers;
    for (const std::string &r : kRecords) {
        for (size_t i = 0; i + kK <= r.size(); ++i) {
            kmers.insert(r.substr(i, kK));
        }
    }
    std::vector<std::pair<std::string, std::string>> oriented;
    if (rc(p) == p) {
        oriented.emplace_back("=", p);
    } else {
        oriented.emplace_back("+", p);
        oriented.emplace_back("-", rc(p));
    }
    std::set<std::tuple<std::string, std::string, uint64_t>> out;
    for (const auto &[strand, q] : oriented) {
        for (const std::string &kmer : kmers) {
            for (size_t o = 0; o + q.size() <= kK; ++o) {
                if (kmer.compare(o, q.size(), q) == 0)
                    out.emplace(strand, kmer, o);
            }
        }
    }
    return out;
}

Json::Value run(const AnnotatedDBG &anno_graph, const std::string &body,
                PatternDelivery *delivery = nullptr,
                const std::function<Clock::time_point()> &clock = nullptr) {
    return process_pattern_request(parse_pattern_body(body), anno_graph, limits(), "rel",
                                   nullptr, delivery, clock);
}

// the refusal of |body|: {status, code}
std::pair<int, std::string> refusal(const AnnotatedDBG &anno_graph, const std::string &body) {
    try {
        run(anno_graph, body);
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(e.what(), e.body()["error"].asString());
        EXPECT_EQ(e.code(), e.body()["code"].asString());
        EXPECT_EQ(2u, e.body().size());
        return { e.status(), e.code() };
    }
    return { 200, "" };
}

// a clock that reads |start| once (the request's start) and |start + later_ms| after
std::function<Clock::time_point()> jumping_clock(double later_ms) {
    const Clock::time_point start = Clock::now();
    auto calls = std::make_shared<size_t>(0);
    return [=]() {
        return (*calls)++ ? start + std::chrono::duration_cast<Clock::duration>(
                                   std::chrono::duration<double, std::milli>(later_ms))
                          : start;
    };
}


TEST(PatternRoute, Refusals) {
    auto g = tiny();
    const std::string p = "\"patterns\": [{\"dna\": \"AACG\"}]";
    const std::vector<std::tuple<std::string, int, std::string>> cases = {
        { "[1]", 400, "invalid_request" },
        { "{\"patterns\": [", 400, "invalid_request" },
        { "{}", 400, "invalid_request" },
        { "{\"patterns\": []}", 400, "invalid_request" },
        { "{" + p + ", \"bogus\": 1}", 400, "invalid_request" },
        { "{\"patterns\": [{\"dna\": \"AACG\", \"x\": 1}]}", 400, "invalid_request" },
        { "{\"patterns\": [{\"dna\": \"AACG\", \"iupac\": \"AACG\"}]}", 400, "invalid_request" },
        { "{" + p + ", \"mode\": \"labels\"}", 400, "invalid_request" },
        { "{" + p + ", \"max_steps\": 0}", 400, "invalid_request" },
        { "{" + p + ", \"time_budget_ms\": 250}", 400, "invalid_request" },
        { "{" + p + ", \"output\": {\"labels\": \"all\"}}", 400, "later_increment" },
        { "{" + p + ", \"mode\": \"count\", \"output\": {\"labels\": \"predicate_only\"}}",
          400, "later_increment" },
        { "{" + p + ", \"output\": {\"occurrences\": true}}", 400, "later_increment" },
        { "{" + p + ", \"predicate\": null}", 400, "later_increment" },
        { "{\"patterns\": [{\"protein\": \"MK\"}]}", 400, "later_increment" },
        { "{" + p + ", \"in_ram\": true}", 400, "resident_only" },
        { "{" + p + ", \"output\": {\"labels\": \"none\", \"paths\": false}}", 200, "" },
    };
    for (const auto &[body, status, code] : cases) {
        EXPECT_EQ(std::make_pair(status, code), refusal(*g, body)) << body;
    }
}

TEST(PatternRoute, CountAnswer) {
    auto g = tiny();
    Json::Value out = run(*g, "{\"patterns\": [{\"id\": \"pal\", \"dna\": \"acgt\"}, "
                              "{\"dna\": \"ACGU\"}], \"mode\": \"count\"}");
    EXPECT_EQ(1, out["pattern_contract_version"].asInt());
    EXPECT_EQ("count", out["mode"].asString());
    EXPECT_TRUE(out["output"].isNull());
    EXPECT_EQ("rel", out["index"]["release"].asString());
    EXPECT_EQ(kK, out["index"]["k"].asUInt64());
    EXPECT_EQ("basic", out["index"]["graph_mode"].asString());
    EXPECT_EQ(4.0, out["limits"]["min_information_bits"].asDouble());
    EXPECT_EQ(0u, out["limits"]["clamped"].size());
    ASSERT_EQ(2u, out["patterns"].size());

    const Json::Value &pal = out["patterns"][0];
    EXPECT_EQ("pal", pal["id"].asString());
    EXPECT_EQ("ACGT", pal["pattern"].asString());
    EXPECT_TRUE(pal["palindromic"].asBool());
    ASSERT_EQ(1u, pal["strands"].size());
    EXPECT_EQ("=", pal["strands"][0].asString());
    const Json::Value &c = pal["counts"]["contexts"];
    EXPECT_EQ(oracle("ACGT").size(), c["value"].asUInt64());
    EXPECT_GT(c["value"].asUInt64(), 0u);
    EXPECT_EQ("exact", c["relation"].asString());
    EXPECT_EQ("graph_contexts", c["unit"].asString());
    EXPECT_EQ(std::vector<std::string>({ "both" }), c["by_strand"].getMemberNames());
    EXPECT_EQ(kK - 4 + 1, c["by_offset"].size());
    EXPECT_TRUE(pal["counts"]["labels"]["value"].isNull());
    EXPECT_EQ("unknown", pal["counts"]["labels"]["relation"].asString());
    EXPECT_EQ("unknown", pal["counts"]["occurrences"]["relation"].asString());
    EXPECT_FALSE(pal["retrieval_complete"].asBool());
    EXPECT_FALSE(pal.isMember("results"));
    EXPECT_FALSE(pal.isMember("withheld"));
    EXPECT_EQ("full", pal["determinism"].asString());
    EXPECT_TRUE(pal["id"].isString());

    const Json::Value &bad = out["patterns"][1];
    EXPECT_TRUE(bad["id"].isNull());
    EXPECT_EQ("bad_alphabet", bad["error"]["code"].asString());
    EXPECT_FALSE(bad.isMember("pattern"));
    EXPECT_FALSE(bad.isMember("counts"));
}

TEST(PatternRoute, RetrievalRows) {
    auto g = tiny();
    // the projection omitted: "none" in this increment, stated in the answer
    Json::Value out = run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}]}");
    EXPECT_EQ("all_or_count", out["mode"].asString());
    EXPECT_EQ("none", out["output"]["labels"].asString());
    const Json::Value &e = out["patterns"][0];
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    EXPECT_TRUE(e["withheld"].isNull());
    EXPECT_TRUE(e["cut"].isNull());
    const auto expected = oracle("AACG");
    ASSERT_EQ(expected.size(), e["results"].size());
    EXPECT_EQ(expected.size(), e["returned"].asUInt64());
    std::set<std::tuple<std::string, std::string, uint64_t>> got;
    std::tuple<uint64_t, uint64_t> previous { 0, 0 };
    for (const Json::Value &r : e["results"]) {
        const std::string kmer = r["kmer"].asString();
        const uint64_t offset = r["offset"].asUInt64();
        EXPECT_EQ(kmer.substr(offset, 4), r["instance"].asString());
        EXPECT_EQ(r["strand"].asString() == "-" ? rc("AACG") : "AACG", r["instance"].asString());
        // BASIC: the node is the stored k-mer's, its row the annotation's row of it
        EXPECT_EQ(r["node"].asUInt64() - 1, r["row"].asUInt64());
        EXPECT_LT(r["row"].asUInt64(), g->get_annotator().num_objects());
        EXPECT_EQ(kmer, g->get_graph().get_node_sequence(r["node"].asUInt64()));
        std::tuple<uint64_t, uint64_t> key { r["node"].asUInt64(), offset };
        EXPECT_LT(previous, key);
        previous = key;
        got.emplace(r["strand"].asString(), kmer, offset);
    }
    EXPECT_EQ(expected, got);
}

TEST(PatternRoute, CanonicalRowsAreTheAnnotationKeys) {
    // both orientations of a k-mer are stored; the annotation is keyed on the canonical one,
    // so a context's row is its k-mer's canonical node, shared with its reverse complement
    auto g = tiny(DeBruijnGraph::CANONICAL);
    Json::Value out = run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}]}");
    EXPECT_EQ("canonical", out["index"]["graph_mode"].asString());
    const Json::Value &results = out["patterns"][0]["results"];
    ASSERT_GT(results.size(), 0u);
    std::map<std::string, std::pair<uint64_t, uint64_t>> node_row;
    for (const Json::Value &r : results) {
        EXPECT_TRUE(r.isMember("orientation"));
        EXPECT_FALSE(r.isMember("strand"));
        node_row[r["kmer"].asString()] = { r["node"].asUInt64(), r["row"].asUInt64() };
    }
    for (const auto &[kmer, nr] : node_row) {
        auto it = node_row.find(rc(kmer));
        ASSERT_TRUE(it != node_row.end()) << kmer;
        EXPECT_EQ(nr.second, it->second.second) << kmer;
        EXPECT_EQ(std::min(nr.first, it->second.first) - 1, nr.second) << kmer;
    }
}

TEST(PatternRoute, TimeStopAnsweredWithinTheReserve) {
    auto g = tiny();
    PatternDelivery delivery;
    // past the work time (1000 - 250 ms), before the deadline: the counts stop, the answer
    // is written
    Json::Value out = run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}], \"time_budget_ms\": 1000}",
                          &delivery, jumping_clock(900));
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("time", e["stop"]["reason"].asString());
    EXPECT_EQ("time_limited", e["determinism"].asString());
    EXPECT_NE("exact", e["counts"]["contexts"]["relation"].asString());
    EXPECT_EQ("deadline", e["withheld"]["reason"].asString());
    EXPECT_EQ(0u, e["results"].size());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
}

TEST(PatternRoute, DeadlineOverrunIs503) {
    auto g = tiny();
    PatternDelivery delivery;
    try {
        run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}], \"time_budget_ms\": 1000}", &delivery,
            jumping_clock(1001));
        FAIL() << "answered past the deadline";
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(503, e.status());
        EXPECT_EQ("deadline", e.code());
    }
}

TEST(PatternRoute, DeliveryCheck) {
    PatternDelivery delivery;
    // before the request was parsed there is no deadline to miss
    EXPECT_NO_THROW(delivery.check());
    const Clock::time_point start = Clock::now();
    delivery.set_deadline(pattern::Deadline(start, 1000, 250, [start]() {
        return start + std::chrono::milliseconds(999);
    }));
    EXPECT_NO_THROW(delivery.check());
    delivery.set_deadline(pattern::Deadline(start, 1000, 250, [start]() {
        return start + std::chrono::milliseconds(1000);
    }));
    EXPECT_THROW(delivery.check(), PatternRefusal);
}

TEST(PatternRoute, MaskRequired) {
    auto graph = std::make_shared<DBGSuccinct>(kK);
    for (const std::string &r : kRecords) {
        graph->add_sequence(r);
    }
    ASSERT_EQ(nullptr, graph->get_mask());
    AnnotatedDBG anno_graph(graph,
                            std::make_unique<annot::ColumnCompressed<>>(graph->max_index()));
    EXPECT_EQ(std::make_pair(400, std::string("mask_required")),
              refusal(anno_graph, "{\"patterns\": [{\"dna\": \"AACG\"}]}"));
    Json::Value caps = pattern_capabilities_json(&anno_graph, limits(), false);
    EXPECT_FALSE(caps["available"].asBool());
    EXPECT_EQ("mask_required", caps["unavailable_reason"].asString());
    EXPECT_EQ("absent", caps["mask"].asString());
}

TEST(PatternRoute, Capabilities) {
    auto g = tiny();
    Json::Value caps = pattern_capabilities_json(g.get(), limits(), false);
    EXPECT_TRUE(caps["available"].asBool());
    EXPECT_TRUE(caps["unavailable_reason"].isNull());
    EXPECT_EQ("file", caps["mask"].asString());
    EXPECT_EQ("basic", caps["graph_mode"].asString());
    EXPECT_EQ(kK, caps["k"].asUInt64());
    EXPECT_EQ("none", caps["default_projection"].asString());
    EXPECT_EQ(4.0, caps["caps"]["min_information_bits"].asDouble());
    EXPECT_EQ("none", caps["placement"].asString());

    // loading: nothing about the graph is known yet
    caps = pattern_capabilities_json(nullptr, limits(), false);
    EXPECT_TRUE(caps["available"].isNull());
    EXPECT_TRUE(caps["graph_mode"].isNull());
    EXPECT_TRUE(caps["k"].isNull());

    caps = pattern_capabilities_json(g.get(), limits(), true);
    EXPECT_EQ(3u, caps.size());
    EXPECT_FALSE(caps["available"].asBool());
    EXPECT_EQ("multi_graph_later_increment", caps["unavailable_reason"].asString());
}

} // namespace
