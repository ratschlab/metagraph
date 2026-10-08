#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
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
                const std::function<Clock::time_point()> &clock = nullptr,
                const PatternLimits &caps = limits(), const RetrievalHooks *hooks = nullptr) {
    return process_pattern_request(parse_pattern_body(body), anno_graph, caps, "rel",
                                   nullptr, delivery, clock, hooks);
}

// the refusal of |body|: {status, code}, and its message in |error|
std::pair<int, std::string> refusal(const AnnotatedDBG &anno_graph, const std::string &body,
                                    std::string *error = nullptr) {
    try {
        run(anno_graph, body);
    } catch (const PatternRefusal &e) {
        EXPECT_EQ(e.what(), e.body()["error"].asString());
        EXPECT_EQ(e.code(), e.body()["code"].asString());
        EXPECT_EQ(2u, e.body().size());
        if (error)
            *error = e.what();
        return { e.status(), e.code() };
    }
    return { 200, "" };
}

// a clock under the test's control: |ms| past its start
struct VirtualClock {
    const Clock::time_point start = Clock::now();
    std::shared_ptr<double> ms = std::make_shared<double>(0);
    std::shared_ptr<uint64_t> reads = std::make_shared<uint64_t>(0);
    std::function<Clock::time_point()> fn() const {
        return [start = start, ms = ms, reads = reads]() {
            ++*reads;
            return start + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double, std::milli>(*ms));
        };
    }
};

std::string repeat(const std::string &s, size_t n) {
    std::string out;
    for (size_t i = 0; i < n; ++i) {
        out += s;
    }
    return out;
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
        // labels "all" (increment 3) on a column annotation: no budget-aware decode
        { "{" + p + ", \"output\": {\"labels\": \"all\"}}", 400, "annotation_unbudgeted" },
        { "{" + p + ", \"output\": {\"labels\": \"all\"}, \"allow_unbudgeted_annotation\": "
          "true}", 200, "" },
        // mode count reads no annotation, whatever the projection
        { "{" + p + ", \"mode\": \"count\", \"output\": {\"labels\": \"all\"}}", 200, "" },
        { "{" + p + ", \"mode\": \"count\", \"output\": {\"labels\": \"predicate_only\"}}",
          400, "later_increment" },
        // occurrences are placed per label: they need labels "all"
        { "{" + p + ", \"output\": {\"occurrences\": true}}", 400, "invalid_request" },
        { "{" + p + ", \"max_labels_per_anchor\": 0}", 400, "invalid_request" },
        { "{" + p + ", \"max_memory_mb\": 0}", 400, "invalid_request" },
        { "{" + p + ", \"max_annotation_work\": 0}", 400, "invalid_request" },
        { "{" + p + ", \"max_labels\": -1}", 400, "invalid_request" },
        { "{" + p + ", \"allow_unbudgeted_annotation\": 1}", 400, "invalid_request" },
        { "{" + p + ", \"max_labels\": 0, \"max_occurrences_per_label\": 0}", 200, "" },
        { "{" + p + ", \"predicate\": null}", 400, "later_increment" },
        // owner decision #13: long_search, reserved until increment 4, is served (paths
        // opt-in); its values are "anchors" (the default) and "paths", nothing else
        { "{" + p + ", \"long_search\": \"paths\"}", 200, "" },
        { "{" + p + ", \"long_search\": \"anchors\"}", 200, "" },
        { "{" + p + ", \"long_search\": null}", 400, "invalid_request" },
        { "{" + p + ", \"long_search\": \"Paths\"}", 400, "invalid_request" },
        { "{" + p + ", \"long_search\": true}", 400, "invalid_request" },
        // max_paths: an integer >= 0 like max_anchors, accepted with any request
        { "{" + p + ", \"max_paths\": 0}", 200, "" },
        { "{" + p + ", \"max_paths\": -1}", 400, "invalid_request" },
        { "{" + p + ", \"max_paths\": 1.5}", 400, "invalid_request" },
        { "{" + p + ", \"max_paths\": null}", 400, "invalid_request" },
        // require_support (owner decision #14): label_intersection or record_verified
        { "{" + p + ", \"require_support\": \"label_intersection\"}", 200, "" },
        { "{" + p + ", \"require_support\": \"record_verified\"}", 200, "" },
        { "{" + p + ", \"require_support\": \"kmer\"}", 400, "invalid_request" },
        { "{" + p + ", \"require_support\": null}", 400, "invalid_request" },
        // the verification reads the coordinates, which occurrences false declines
        { "{" + p + ", \"require_support\": \"record_verified\", \"output\": {\"labels\": "
          "\"all\", \"occurrences\": false}, \"allow_unbudgeted_annotation\": true}", 400,
          "invalid_request" },
        // record_verified on an index that cannot verify (no coordinates): refused where it
        // acts (labels read, long_search paths), not in a request that reads no labels
        { "{" + p + ", \"require_support\": \"record_verified\", \"long_search\": \"paths\", "
          "\"output\": {\"labels\": \"all\"}, \"allow_unbudgeted_annotation\": true}", 400,
          "support_unavailable" },
        { "{" + p + ", \"require_support\": \"record_verified\", \"long_search\": \"paths\", "
          "\"mode\": \"count\", \"output\": {\"labels\": \"all\"}}", 200, "" },
        { "{" + p + ", \"require_support\": \"label_intersection\", \"long_search\": "
          "\"paths\", \"output\": {\"labels\": \"all\"}, \"allow_unbudgeted_annotation\": "
          "true}", 200, "" },
        // increment 5 (owner decision #15): protein is the third kind, genetic_code an NCBI
        // translation table id (default 1); an unknown id has its own code
        { "{\"patterns\": [{\"protein\": \"MK\"}]}", 200, "" },
        { "{\"patterns\": [{\"protein\": \"MK\"}], \"genetic_code\": 2}", 200, "" },
        { "{\"patterns\": [{\"protein\": \"MK\"}], \"genetic_code\": 33.0}", 200, "" },
        { "{" + p + ", \"genetic_code\": 11}", 200, "" },
        { "{\"patterns\": [{\"protein\": 5}]}", 400, "invalid_request" },
        { "{\"patterns\": [{\"protein\": \"MK\", \"dna\": \"ACG\"}]}", 400,
          "invalid_request" },
        { "{\"patterns\": [{\"protein\": \"MK\", \"iupac\": \"ACG\"}]}", 400,
          "invalid_request" },
        { "{\"patterns\": [{\"protein\": null}]}", 400, "invalid_request" },
        { "{\"patterns\": [{\"protein\": \"MK\"}], \"genetic_code\": 7}", 400,
          "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": 0}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": -1}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": 34}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": 17}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": 18446744073709551615}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": 4294967297}", 400, "genetic_code_unknown" },
        { "{" + p + ", \"genetic_code\": \"1\"}", 400, "invalid_request" },
        { "{" + p + ", \"genetic_code\": 1.5}", 400, "invalid_request" },
        { "{" + p + ", \"genetic_code\": null}", 400, "invalid_request" },
        { "{" + p + ", \"genetic_code\": true}", 400, "invalid_request" },
        { "{" + p + ", \"in_ram\": true}", 400, "resident_only" },
        { "{" + p + ", \"output\": {\"labels\": \"none\", \"paths\": false}}", 200, "" },
        // increment 4: output.paths is accepted with either value and changes nothing (a path
        // result always carries its node path)
        { "{" + p + ", \"output\": {\"labels\": \"none\", \"paths\": true}}", 200, "" },
        { "{" + p + ", \"output\": {\"paths\": 1}}", 400, "invalid_request" },
        // review of 2026-10-07, R1-02: one RFC 8259 JSON text with unique member names, nothing
        // else (each was answered 200 as if it were the leading object, or with the last of a
        // duplicated member's values)
        { "{" + p + "} GARBAGE", 400, "invalid_request" },
        { "{" + p + ",}", 400, "invalid_request" },
        { "{\"patterns\": [{\"dna\": \"AACG\"},]}", 400, "invalid_request" },
        { "{" + p + ", /* x */ \"mode\": \"count\"}", 400, "invalid_request" },
        { "{" + p + ", // x\n \"mode\": \"count\"}", 400, "invalid_request" },
        { "{" + p + "} {" + p + "}", 400, "invalid_request" },
        { "{" + p + "}]", 400, "invalid_request" },
        { "{" + p + ", \"mode\": \"count\", \"mode\": \"partial\"}", 400, "invalid_request" },
        { "{" + p + ", " + p + "}", 400, "invalid_request" },
        { "{" + p + ", \"max_steps\": 1, \"max_steps\": 100000}", 400, "invalid_request" },
        { "{\"patterns\": [{\"dna\": \"AACG\", \"dna\": \"ACGT\"}]}", 400, "invalid_request" },
        // trailing white space is not extra
        { "{" + p + "} \n\t ", 200, "" },
        // review of 2026-10-07, R1-03, R2-01: past jsoncpp's nesting limit (it throws, and
        // the route answered 400 without a code), at the top level and inside a field
        { repeat("[", 1001) + repeat("]", 1001), 400, "invalid_request" },
        { "{" + p + ", \"x\": " + repeat("[", 1500) + repeat("]", 1500) + "}", 400,
          "invalid_request" },
        { repeat("{\"a\": ", 1001) + "1" + repeat("}", 1001), 400, "invalid_request" },
        // at the limit: parsed, then refused for what it is
        { repeat("[", 1000) + repeat("]", 1000), 400, "invalid_request" },
    };
    for (const auto &[body, status, code] : cases) {
        EXPECT_EQ(std::make_pair(status, code), refusal(*g, body)) << body.substr(0, 200);
    }
    std::string error;
    refusal(*g, repeat("[", 1001) + repeat("]", 1001), &error);
    EXPECT_NE(std::string::npos, error.find("request: not JSON")) << error;
    refusal(*g, repeat("[", 1000) + repeat("]", 1000), &error);
    EXPECT_EQ("request: expected an object", error);
    refusal(*g, "{" + p + ", \"mode\": \"count\", \"mode\": \"partial\"}", &error);
    EXPECT_NE(std::string::npos, error.find("Duplicate key: 'mode'")) << error;
}

// SPEC §5: a request is refused by the first check it fails, in order (review of 2026-10-07,
// R2-03: every other refusal test sends one fault). Each body has two faults; the code and the
// field named first in the message are the earlier check's
TEST(PatternRoute, RefusalOrder) {
    auto g = tiny();
    const std::string p = "\"patterns\": [{\"dna\": \"AACG\"}]";
    const std::vector<std::tuple<std::string, std::string, std::string>> cases = {
        // 6 before 7: a later-increment field before the patterns
        { "{\"patterns\": \"x\", \"predicate\": 1}", "later_increment", "request.predicate" },
        // 6, alphabetical: graphs < in_ram; in_ram < max_paths
        { "{" + p + ", \"in_ram\": true, \"graphs\": []}", "later_increment", "request.graphs" },
        { "{" + p + ", \"in_ram\": true, \"max_paths\": 1}", "resident_only", "request.in_ram" },
        // (long_search is served since increment 4: checked after the patterns, step 8)
        { "{" + p + ", \"in_ram\": true, \"long_search\": \"paths\"}", "resident_only",
          "request.in_ram" },
        { "{\"patterns\": \"x\", \"long_search\": \"x\"}", "invalid_request",
          "request.patterns" },
        // 7 before 8
        { "{\"patterns\": [], \"mode\": \"x\"}", "invalid_request", "request.patterns" },
        // within 7: id, exactly one of dna / iupac / protein, its type, an unknown field
        { "{\"patterns\": [{\"protein\": \"MK\", \"id\": 1}]}", "invalid_request",
          "request.patterns[0].id" },
        { "{\"patterns\": [{\"id\": 1, \"dna\": \"A\", \"iupac\": \"A\"}]}",
          "invalid_request", "request.patterns[0].id" },
        { "{\"patterns\": [{\"protein\": \"M\", \"dna\": 1}], \"mode\": \"x\"}",
          "invalid_request", "request.patterns[0]: expected exactly one of 'dna', 'iupac', "
          "'protein'" },
        { "{\"patterns\": [{\"protein\": 1, \"x\": 1}]}", "invalid_request",
          "request.patterns[0].protein" },
        { "{\"patterns\": [{\"dna\": 1, \"x\": 1}]}", "invalid_request",
          "request.patterns[0].dna" },
        { "{\"patterns\": [{\"dna\": \"A\", \"x\": 1}], \"mode\": \"x\"}",
          "invalid_request", "request.patterns[0]: unknown field" },
        // within 8: mode, output, scope, strands, stop_at_threshold, the caps, the time budget,
        // then increment 3's
        { "{" + p + ", \"mode\": \"x\", \"output\": {\"labels\": \"predicate_only\"}}",
          "invalid_request", "request.mode" },
        { "{" + p + ", \"output\": {\"labels\": \"x\"}, \"scope\": \"x\"}", "invalid_request",
          "request.output.labels" },
        { "{" + p + ", \"output\": {\"paths\": 1, \"x\": 1}}", "invalid_request",
          "request.output.paths" },
        { "{" + p + ", \"output\": {\"paths\": true, \"x\": 1}}", "invalid_request",
          "request.output: unknown field 'x'" },
        { "{" + p + ", \"scope\": \"x\", \"strands\": \"x\"}", "invalid_request",
          "request.scope" },
        { "{" + p + ", \"strands\": \"x\", \"stop_at_threshold\": 1}", "invalid_request",
          "request.strands" },
        { "{" + p + ", \"stop_at_threshold\": 1, \"max_contexts\": -1}", "invalid_request",
          "request.stop_at_threshold" },
        { "{" + p + ", \"max_contexts\": -1, \"max_anchors\": -1}", "invalid_request",
          "request.max_contexts" },
        { "{" + p + ", \"max_anchors\": -1, \"max_steps\": 0}", "invalid_request",
          "request.max_anchors" },
        { "{" + p + ", \"max_steps\": 0, \"time_budget_ms\": 1}", "invalid_request",
          "request.max_steps" },
        { "{" + p + ", \"time_budget_ms\": 1, \"max_labels_per_anchor\": 0}", "invalid_request",
          "request.time_budget_ms" },
        { "{" + p + ", \"max_labels_per_anchor\": 0, \"max_annotation_work\": 0}",
          "invalid_request", "request.max_labels_per_anchor" },
        { "{" + p + ", \"max_occurrences_per_label\": -1, "
          "\"allow_unbudgeted_annotation\": 1}", "invalid_request",
          "request.max_occurrences_per_label" },
        // then increment 4's: long_search, max_paths, require_support (and its contradiction
        // with output.occurrences false)
        { "{" + p + ", \"allow_unbudgeted_annotation\": 1, \"long_search\": \"x\"}",
          "invalid_request", "request.allow_unbudgeted_annotation" },
        { "{" + p + ", \"long_search\": \"x\", \"max_paths\": -1}", "invalid_request",
          "request.long_search" },
        { "{" + p + ", \"max_paths\": -1, \"require_support\": \"x\"}", "invalid_request",
          "request.max_paths" },
        { "{" + p + ", \"require_support\": \"x\", \"bogus\": 1}", "invalid_request",
          "request.require_support" },
        { "{" + p + ", \"require_support\": \"record_verified\", \"output\": {\"labels\": "
          "\"all\", \"occurrences\": false}, \"genetic_code\": 7}", "invalid_request",
          "request.require_support" },
        // then increment 5's genetic_code (its type, then the table)
        { "{" + p + ", \"genetic_code\": \"x\", \"bogus\": 1}", "invalid_request",
          "request.genetic_code: expected an integer" },
        { "{" + p + ", \"genetic_code\": 7, \"bogus\": 1}", "genetic_code_unknown",
          "request.genetic_code: 7 is not an NCBI translation table id" },
        // a peptide's slot error is no refusal: the genetic code's comes first
        { "{\"patterns\": [{\"protein\": \"M*U\"}], \"genetic_code\": 8}",
          "genetic_code_unknown", "request.genetic_code" },
        // 8 before 9
        { "{" + p + ", \"max_steps\": 0, \"bogus\": 1}", "invalid_request", "request.max_steps" },
        // 9 before 10: on this column annotation labels "all" would be annotation_unbudgeted
        { "{" + p + ", \"output\": {\"labels\": \"all\"}, \"bogus\": 1}", "invalid_request",
          "request: unknown field 'bogus'" },
        // 10 before 11: the unbudgeted annotation before the support it cannot give
        { "{" + p + ", \"output\": {\"labels\": \"all\"}, \"long_search\": \"paths\", "
          "\"require_support\": \"record_verified\"}", "annotation_unbudgeted",
          "pattern: output.labels" },
    };
    for (const auto &[body, code, first] : cases) {
        std::string error;
        EXPECT_EQ(std::make_pair(400, code), refusal(*g, body, &error)) << body;
        EXPECT_EQ(0u, error.find(first)) << body << ": " << error;
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

TEST(PatternRoute, AnchorInformationBitsKeepTheirMeaning) {
    // review of 2026-10-07, X-GUARANTEES-01 and the owner's decision of the same day:
    // anchor_information_bits stays the bits of P[0, k) (contract version 1);
    // min_anchor_information_bits, an addition, is the least searched anchor window's, the
    // floor's operand for L > k
    auto g = tiny();
    Json::Value out = run(*g, "{\"patterns\": [{\"iupac\": \"ACGTTGCAAN\"}, "
                              "{\"iupac\": \"ACGTTGCNNNNNN\"}, {\"dna\": \"ACGT\"}], "
                              "\"mode\": \"count\"}");
    ASSERT_EQ(3u, out["patterns"].size());
    // both windows searched: P[0, 7) = ACGTTGC 14 bits, P[3, 10) = TTGCAAN 12 bits
    const Json::Value &both = out["patterns"][0];
    EXPECT_FALSE(both.isMember("error"));
    EXPECT_EQ(14.0, both["anchor_information_bits"].asDouble());
    EXPECT_EQ(12.0, both["min_anchor_information_bits"].asDouble());
    // the reverse window CNNNNNN (2 bits) is below the floor of 4: refused, both stated
    const Json::Value &refused = out["patterns"][1];
    EXPECT_EQ("information_below_floor", refused["error"]["code"].asString());
    EXPECT_EQ(14.0, refused["anchor_information_bits"].asDouble());
    EXPECT_EQ(2.0, refused["min_anchor_information_bits"].asDouble());
    // L <= k: neither
    const Json::Value &shorter = out["patterns"][2];
    EXPECT_TRUE(shorter.isMember("anchor_information_bits"));
    EXPECT_TRUE(shorter["anchor_information_bits"].isNull());
    EXPECT_TRUE(shorter.isMember("min_anchor_information_bits"));
    EXPECT_TRUE(shorter["min_anchor_information_bits"].isNull());

    // forward only: the forward window alone, 14 bits, admitted
    out = run(*g, "{\"patterns\": [{\"iupac\": \"ACGTTGCNNNNNN\"}], \"mode\": \"count\", "
                  "\"strands\": \"forward\"}");
    const Json::Value &forward = out["patterns"][0];
    EXPECT_FALSE(forward.isMember("error"));
    EXPECT_EQ(14.0, forward["anchor_information_bits"].asDouble());
    EXPECT_EQ(14.0, forward["min_anchor_information_bits"].asDouble());
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

// X-TESTS-02 (review of 2026-10-07): the work done, its time passing in the finalisation
// reserve: the answer is the work's (exact counts, determinism full, no stop), written in the
// reserve; past the deadline, 503
TEST(PatternRoute, WorkTimePassingAfterTheWork) {
    auto g = tiny();
    for (double after : { 900.0, 1001.0 }) {
        VirtualClock clock;
        RetrievalHooks hooks;
        hooks.work_done_hook = [&]() { *clock.ms = after; };
        PatternDelivery delivery;
        try {
            Json::Value out = run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}], "
                                      "\"time_budget_ms\": 1000}", &delivery, clock.fn(),
                                  limits(), &hooks);
            ASSERT_EQ(900.0, after) << "answered past the deadline";
            const Json::Value &e = out["patterns"][0];
            EXPECT_EQ("exact", e["counts"]["contexts"]["relation"].asString());
            EXPECT_EQ(oracle("AACG").size(), e["counts"]["contexts"]["value"].asUInt64());
            EXPECT_EQ("full", e["determinism"].asString());
            EXPECT_TRUE(e["stop"].isNull());
            EXPECT_TRUE(e["retrieval_complete"].asBool());
            EXPECT_EQ(oracle("AACG").size(), e["results"].size());
            EXPECT_EQ(900.0, out["timing"]["elapsed_ms"].asDouble());
            EXPECT_NO_THROW(delivery.check());
        } catch (const PatternRefusal &e) {
            EXPECT_EQ(1001.0, after) << e.what();
            EXPECT_EQ(503, e.status());
            EXPECT_EQ("deadline", e.code());
        }
    }
}

// random records whose graph has more than kDeliveryStride contexts of a 2-mer
std::unique_ptr<AnnotatedDBG> stride_graph() {
    std::mt19937 rng(7);
    const char bases[] = "ACGT";
    std::vector<std::string> records;
    for (size_t r = 0; r < 20; ++r) {
        std::string seq;
        for (size_t i = 0; i < 1000; ++i) {
            seq += bases[rng() % 4];
        }
        records.push_back(seq);
    }
    std::vector<std::string> labels(records.size(), "r");
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(15, records, labels);
}

// X-TESTS-02 (review of 2026-10-07): the answer's assembly reads the deadline every
// kDeliveryStride objects (two patterns of 10,000 results: once after each), and a deadline
// passing there is a 503 before the assembly ends — the final check would read a clock not
// yet past it
TEST(PatternRoute, DeadlineReadWhileTheAnswerIsAssembled) {
    auto g = stride_graph();
    const std::string body = "{\"patterns\": [{\"dna\": \"AC\"}, {\"dna\": \"GG\"}], "
                             "\"mode\": \"partial\", \"time_budget_ms\": 1000}";
    for (uint64_t pass_at : { 3u, 4u, 0u }) {
        VirtualClock clock;
        RetrievalHooks hooks;
        // from the work's end, the clock reads 900 ms, and 1001 ms from the |pass_at|-th
        // reading on (the first is timing.elapsed_ms, then one per stride, then the final)
        auto after = std::make_shared<uint64_t>(0);
        bool done = false;
        hooks.work_done_hook = [&]() { done = true; *clock.ms = 900; };
        auto inner = clock.fn();
        auto fn = [&, inner, after]() {
            if (done && ++*after >= pass_at && pass_at)
                *clock.ms = 1001;
            return inner();
        };
        PatternDelivery delivery;
        try {
            Json::Value out = run(*g, body, &delivery, fn, limits(), &hooks);
            ASSERT_EQ(0u, pass_at) << "answered past the deadline";
            ASSERT_EQ(2u, out["patterns"].size());
            EXPECT_EQ(10000u, out["patterns"][0]["returned"].asUInt64());
            EXPECT_EQ(10000u, out["patterns"][1]["returned"].asUInt64());
            // timing, a check after each pattern's 10,001 objects, the final check
            EXPECT_EQ(4u, *after);
        } catch (const PatternRefusal &e) {
            EXPECT_EQ(503, e.status());
            EXPECT_NE(0u, pass_at);
            // thrown at the first check that read the clock past the deadline: from the 3rd
            // reading on, the second stride check, before the assembly ended (without the
            // stride checks the final one, the 2nd reading, would answer 200); from the 4th,
            // the final check
            EXPECT_EQ(pass_at, *after);
        }
    }
}

// X-EFFICIENCY-04 (review of 2026-10-07): the work stops earlier by the estimated time to
// write what the answer holds, so that a request whose patterns buffered many results is
// answered (a time stop, its counts kept) rather than 503. A clock that never moves: only the
// answer's volume can stop the work. At a build rate of 1 byte per second the first pattern's
// results alone take longer to write than the default budget: the second pattern is stopped
// by time where it starts, the answer 200
TEST(PatternRoute, TheAnswersVolumeStopsTheWork) {
    auto g = tiny();
    const std::string body = "{\"patterns\": [{\"dna\": \"AACG\"}, {\"dna\": \"ACGT\"}], "
                             "\"mode\": \"partial\"}";
    VirtualClock clock;
    PatternLimits slow = limits();
    slow.delivery_build_mbps = 1e-6;
    PatternDelivery delivery;
    Json::Value out = run(*g, body, &delivery, clock.fn(), slow);
    const Json::Value &first = out["patterns"][0];
    EXPECT_TRUE(first["retrieval_complete"].asBool());
    EXPECT_EQ(oracle("AACG").size(), first["results"].size());
    EXPECT_EQ("full", first["determinism"].asString());
    const Json::Value &second = out["patterns"][1];
    EXPECT_EQ("{\"phase\":\"discovery\",\"reason\":\"time\"}", [&]() {
        Json::StreamWriterBuilder b;
        b["indentation"] = "";
        return Json::writeString(b, second["stop"]);
    }());
    EXPECT_EQ("unknown", second["counts"]["contexts"]["relation"].asString());
    EXPECT_EQ("time", second["cut"]["reason"].asString());
    EXPECT_EQ(0u, second["returned"].asUInt64());
    EXPECT_EQ("time_limited", second["determinism"].asString());
    // the deadline itself never passed: the answer is written
    EXPECT_NO_THROW(delivery.check());

    // the same request without the model: both patterns complete
    PatternLimits none = limits();
    none.delivery_build_mbps = 0;
    none.delivery_compress_mbps = 0;
    VirtualClock still;
    out = run(*g, body, nullptr, still.fn(), none);
    EXPECT_TRUE(out["patterns"][1]["retrieval_complete"].asBool());
    EXPECT_EQ(oracle("ACGT").size(), out["patterns"][1]["results"].size());
    // and with the default rates: the tiny answer leaves the default budget's work time
    VirtualClock defaults;
    out = run(*g, body, nullptr, defaults.fn());
    EXPECT_TRUE(out["patterns"][1]["retrieval_complete"].asBool());
}

TEST(PatternRoute, AnswerVolumeModel) {
    // 1 MB at 10 MB/s (100 ms) and at 50 MB/s (20 ms), with the margin
    AnswerVolume v(10, 50);
    EXPECT_DOUBLE_EQ(0, v.finalize_ms());
    v.add(1'000'000);
    EXPECT_DOUBLE_EQ(kAnswerVolumeMargin * 120, v.finalize_ms());
    // what is still to be built is built at the build rate as well
    v.add_pending(1'000'000);
    EXPECT_DOUBLE_EQ(kAnswerVolumeMargin * (240 + 100), v.finalize_ms());
    v.settle(1'000'000);
    EXPECT_EQ(0u, v.pending_bytes());
    EXPECT_EQ(2'000'000u, v.text_bytes());
    EXPECT_DOUBLE_EQ(kAnswerVolumeMargin * 240, v.finalize_ms());
    // the CLI's indented text
    AnswerVolume indented(10, 50, 2);
    indented.add(1'000'000);
    EXPECT_DOUBLE_EQ(kAnswerVolumeMargin * 240, indented.finalize_ms());
    // no model
    AnswerVolume off(0, std::numeric_limits<double>::infinity());
    off.add(1'000'000'000);
    EXPECT_DOUBLE_EQ(0, off.finalize_ms());
}

// compact_json_bytes is the compact text's length: exact for ASCII without doubles, from above
// otherwise (escapes, non-ASCII, doubles)
TEST(PatternRoute, CompactJsonBytes) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    auto text = [&](const Json::Value &v) { return Json::writeString(b, v).size(); };
    auto g = tiny();
    const Json::Value out = run(*g, "{\"patterns\": [{\"dna\": \"AACG\"}], \"mode\": \"partial\"}");
    ASSERT_GT(out["patterns"][0]["results"].size(), 0u);
    for (const Json::Value &r : out["patterns"][0]["results"]) {
        EXPECT_EQ(text(r), compact_json_bytes(r));
    }
    Json::Value v;
    v["a"] = Json::Value(Json::arrayValue);
    v["b"] = Json::Value();
    v["c"] = true;
    v["d"] = false;
    v["e"] = -12345;
    v["f"] = Json::UInt64(18446744073709551615ull);
    v["g"] = "plain";
    v["h"]["nested"].append(1);
    v["h"]["nested"].append("x");
    EXPECT_EQ(text(v), compact_json_bytes(v));
    for (const Json::Value &x : { Json::Value("a\"b\\c\n\x01"), Json::Value("caf\xc3\xa9"),
                                  Json::Value(0.1), Json::Value(-1.5e300), Json::Value(1.0 / 3) }) {
        EXPECT_GE(compact_json_bytes(x), text(x)) << text(x);
    }
}

// T3-07 (review of 2026-10-07): the 503 states a fractional budget as given (a cast floored
// 1000.5 to 1000)
TEST(PatternRoute, DeadlineMessageStatesTheBudget) {
    PatternDelivery delivery;
    const Clock::time_point start = Clock::now();
    delivery.set_deadline(pattern::Deadline(start, 1000.5, 250, [start]() {
        return start + std::chrono::milliseconds(2000);
    }));
    try {
        delivery.check();
        FAIL() << "no 503";
    } catch (const PatternRefusal &e) {
        EXPECT_EQ("pattern: the answer could not be written within time_budget_ms (1000.5 ms, "
                  "the finalisation reserve of 250 ms included): nothing partial is sent",
                  std::string(e.what()));
    }
}

// R1-03, R2-01 (review of 2026-10-07): `metagraph pattern` answers every request file as the
// server would — a refusal's body, or for any other failure the server's 400 body without a
// code — and goes on with the next (an exception escaped and aborted the run)
TEST(PatternRoute, TheCliAnswersEveryRequest) {
    auto g = tiny();
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    auto answer = [&](const std::string &body, const RetrievalHooks *hooks, bool *ok) {
        std::ostringstream out;
        *ok = write_pattern_answer(body, *g, limits(), "rel", nullptr, b, out, "test", hooks);
        Json::Value v;
        std::istringstream in(out.str());
        in >> v;
        EXPECT_EQ('\n', out.str().back());
        return v;
    };
    bool ok = false;
    Json::Value v = answer("{\"patterns\": [{\"dna\": \"AACG\"}]}", nullptr, &ok);
    EXPECT_TRUE(ok);
    EXPECT_EQ(1, v["pattern_contract_version"].asInt());
    v = answer(repeat("[", 1001) + repeat("]", 1001), nullptr, &ok);
    EXPECT_FALSE(ok);
    EXPECT_EQ("invalid_request", v["code"].asString());
    // an unexpected failure (here a read that throws): {"error"} without a code, as the server
    RetrievalHooks failing;
    failing.read_hook = [](size_t) { throw std::runtime_error("a read failed"); };
    v = answer("{\"patterns\": [{\"dna\": \"AACG\"}], \"output\": {\"labels\": \"all\"}, "
               "\"allow_unbudgeted_annotation\": true}", &failing, &ok);
    EXPECT_FALSE(ok);
    EXPECT_EQ("{\"error\":\"a read failed\"}", Json::writeString(b, v));
}

// X-CONCURRENCY-01, R2-02 (review of 2026-10-07): a caller that left (or a server that stops)
// is not answered — the work ends at its first clock reading (Budget::set_abort), the writing
// at its next check — where the route ran to its deadline and wrote into a closed socket
TEST(PatternRoute, AnAbortedRequestIsNotAnswered) {
    auto g = tiny();
    const std::string body = "{\"patterns\": [{\"dna\": \"AACG\"}, {\"dna\": \"ACGT\"}]}";
    size_t asked = 0;
    PatternDelivery delivery;
    delivery.set_abort([&asked]() { ++asked; return true; });
    EXPECT_THROW(run(*g, body, &delivery), pattern::Aborted);
    // asked at the work's readings of the clock, the first of which ended it
    EXPECT_GE(asked, 1u);

    bool gone = false;
    PatternDelivery live;
    live.set_abort([&gone]() { return gone; });
    EXPECT_NO_THROW(run(*g, body, &live));
    EXPECT_NO_THROW(live.check());
    // the caller left after the work: its writing stops at the next check
    gone = true;
    EXPECT_THROW(live.check(), pattern::Aborted);
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

// SPEC §5: the body not JSON (3), then the graph (4), then the body not an object (5)
// (review of 2026-10-07, R2-03), on a graph the engine does not serve (a graph without its
// mask, which mask_required refused before, is served since owner decision #16)
TEST(PatternRoute, TheGraphIsCheckedBetweenTheJsonAndTheObject) {
    auto hash = std::make_shared<DBGHashFast>(kK);
    for (const std::string &r : kRecords) {
        hash->add_sequence(r);
    }
    AnnotatedDBG anno_graph(hash, std::make_unique<annot::ColumnCompressed<>>(hash->max_index()));
    const std::string code = "representation_unsupported";
    EXPECT_EQ(std::make_pair(400, code),
              refusal(anno_graph, "{\"patterns\": [{\"dna\": \"AACG\"}]}"));
    EXPECT_EQ(std::make_pair(400, code), refusal(anno_graph, "[1]"));
    EXPECT_EQ(std::make_pair(400, code),
              refusal(anno_graph, "{\"patterns\": 1, \"predicate\": 1}"));
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(anno_graph, "{\"patterns\": ["));
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(anno_graph, "{\"patterns\": [{\"dna\": \"AACG\"}]} GARBAGE"));
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(anno_graph, repeat("[", 1001) + repeat("]", 1001)));
}

// C2-01, C1-05 (review of 2026-10-07): the two graph reasons of a graph the engine does not
// recognise — a fixture shows representation_unsupported, none can show primary_unwrapped
// (server_query and the CLI always wrap a PRIMARY graph) — each a 400 with its code and
// message, and the capabilities block stating only k
TEST(PatternRoute, GraphsTheEngineDoesNotServe) {
    auto hash = std::make_shared<DBGHashFast>(kK);
    for (const std::string &r : kRecords) {
        hash->add_sequence(r);
    }
    auto primary = std::make_shared<DBGSuccinct>(kK, DeBruijnGraph::PRIMARY);
    const std::vector<std::tuple<std::shared_ptr<DeBruijnGraph>, std::string, std::string>> cases = {
        { hash, "representation_unsupported",
          "pattern: the graph is not a succinct graph (a DBGSuccinct, or a PRIMARY one wrapped "
          "in CanonicalDBG): the pattern lookup narrows BOSS ranges" },
        { primary, "primary_unwrapped",
          "pattern: a PRIMARY graph is served only wrapped in CanonicalDBG" },
    };
    for (const auto &[graph, code, message] : cases) {
        AnnotatedDBG anno_graph(graph,
                                std::make_unique<annot::ColumnCompressed<>>(graph->max_index()));
        std::string error;
        EXPECT_EQ(std::make_pair(400, code),
                  refusal(anno_graph, "{\"patterns\": [{\"dna\": \"AACG\"}]}", &error));
        EXPECT_EQ(message, error);
        const Json::Value caps = pattern_capabilities_json(&anno_graph, limits(), false);
        EXPECT_FALSE(caps["available"].asBool()) << code;
        EXPECT_EQ(code, caps["unavailable_reason"].asString());
        EXPECT_EQ(kK, caps["k"].asUInt64()) << code;
        for (const char *f : { "graph_mode", "alphabet", "strand_stated", "mask", "counting",
                                "dummy_fraction", "scopes", "placement", "support",
                                "annotation" }) {
            EXPECT_TRUE(caps[f].isNull()) << code << " " << f;
        }
    }
}

// Owner decision #4 of 2026-10-07 (review I26): /pattern and `metagraph pattern` serve $ACGT
// graphs only; $ACGTN is refused as alphabet_untested (400, and the capabilities' reason) until
// a DNA5 build passes the pattern tests, while the engine keeps its DNA5 paths. A DNA4 build
// cannot load a DNA5 graph, so the decision is pinned as the pure function of the alphabet it
// is, and through the route on this build's own alphabet (refused on a DNA5 build)
TEST(PatternRoute, AlphabetRefusal) {
    EXPECT_EQ("", alphabet_refusal("$ACGT"));
    EXPECT_EQ("alphabet_untested", alphabet_refusal("$ACGTN"));
    for (const char *other : { "", "ACGT", "$ACGTNacgt", "$ACGTX", "$ACDEFGHIKLMNPQRSTVWYX" }) {
        EXPECT_EQ("alphabet_unsupported", alphabet_refusal(other)) << other;
    }

    auto g = tiny();
    const pattern::GraphSupport engine = pattern::PatternSearch::support(g->get_graph());
    const pattern::GraphSupport route = route_support(g->get_graph());
    // the engine serves both alphabets; the route only $ACGT
    ASSERT_TRUE(engine.supported) << engine.reason;
    const std::string expected = alphabet_refusal(engine.alphabet);
    EXPECT_EQ(expected.empty(), route.supported);
    EXPECT_EQ(expected, route.reason);
    const Json::Value caps = pattern_capabilities_json(g.get(), limits(), false);
    if (expected.empty()) {
        EXPECT_EQ("$ACGT", engine.alphabet);
        EXPECT_TRUE(caps["available"].asBool());
        EXPECT_EQ(std::make_pair(200, std::string()),
                  refusal(*g, "{\"patterns\": [{\"dna\": \"AACG\"}]}"));
    } else {
        EXPECT_EQ("alphabet_untested", expected);
        std::string error;
        EXPECT_EQ(std::make_pair(400, expected),
                  refusal(*g, "{\"patterns\": [{\"dna\": \"AACG\"}]}", &error));
        EXPECT_NE(std::string::npos, error.find("not served on the $ACGTN alphabet until a "
                                                "DNA5 build passes the pattern tests")) << error;
        EXPECT_FALSE(caps["available"].asBool());
        EXPECT_EQ(expected, caps["unavailable_reason"].asString());
        // a recognised graph: what it is stays stated
        EXPECT_EQ("basic", caps["graph_mode"].asString());
        EXPECT_EQ(engine.alphabet, caps["alphabet"].asString());
        EXPECT_EQ("file", caps["mask"].asString());
        EXPECT_EQ(2u, caps["scopes"].size());
    }
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
    // increment 3: the projection "all" is served; predicate_only is still to come
    ASSERT_EQ(2u, caps["projections"].size());
    EXPECT_EQ("none", caps["projections"][0].asString());
    EXPECT_EQ("all", caps["projections"][1].asString());
    ASSERT_EQ(1u, caps["projections_later_increment"].size());
    EXPECT_EQ("predicate_only", caps["projections_later_increment"][0].asString());
    EXPECT_TRUE(caps["default_occurrences"].asBool());
    EXPECT_EQ(4.0, caps["caps"]["min_information_bits"].asDouble());
    EXPECT_EQ(64u, caps["caps"]["max_labels_per_anchor"].asUInt64());
    EXPECT_EQ(100000000u, caps["caps"]["max_annotation_work"].asUInt64());
    EXPECT_EQ(256u, caps["caps"]["max_memory_mb"].asUInt64());
    EXPECT_EQ(1000u, caps["caps"]["max_labels"].asUInt64());
    EXPECT_EQ(16u, caps["caps"]["max_occurrences_per_label"].asUInt64());
    // increment 4: paths opt-in (owner decision #13); without the option a long pattern keeps
    // its anchor-only answer, which long_patterns describes
    EXPECT_EQ(1000u, caps["caps"]["max_paths"].asUInt64());
    EXPECT_EQ("anchors_counted", caps["long_patterns"].asString());
    ASSERT_EQ(2u, caps["long_search"].size());
    EXPECT_EQ("anchors", caps["long_search"][0].asString());
    EXPECT_EQ("paths", caps["long_search"][1].asString());
    EXPECT_EQ("anchors", caps["default_long_search"].asString());
    EXPECT_NE(std::string::npos, caps["caps_rule"].asString().find("max_paths are the maxima"));
    EXPECT_NE(std::string::npos, caps["caps_rule"].asString().find("long_search \"paths\""));
    PatternLimits other = limits();
    other.max_paths = 7;
    EXPECT_EQ(7u, pattern_capabilities_json(g.get(), other, false)["caps"]["max_paths"].asUInt64());
    EXPECT_EQ("none", caps["placement"].asString());
    // no coordinates: a path's labels can only be the intersection of its k-mers'
    EXPECT_EQ("label_intersection", caps["support"].asString());
    // a column annotation: no budget-aware decode
    EXPECT_EQ("unbudgeted", caps["annotation"].asString());
    // increment 5 (owner decision #15): the kind protein, its residues, the genetic codes
    // (every NCBI translation table) and the default, the standard code
    ASSERT_EQ(3u, caps["kinds"].size());
    EXPECT_EQ("protein", caps["kinds"][2].asString());
    EXPECT_TRUE(caps["kinds_later_increment"].isArray());
    EXPECT_EQ(0u, caps["kinds_later_increment"].size());
    std::string residues;
    for (const Json::Value &r : caps["protein_residues"]) {
        ASSERT_EQ(1u, r.asString().size());
        residues += r.asString();
    }
    EXPECT_EQ("ACDEFGHIKLMNPQRSTVWYXBZJ", residues);
    // the list is what the engine parses: each residue alone is a peptide, every other
    // character of the alphabet, the stop '*' included, is not
    for (char c = 'A'; c <= 'Z'; ++c) {
        const bool listed = residues.find(c) != std::string::npos;
        bool parsed = true;
        try {
            pattern::Pattern::parse(pattern::PatternKind::PROTEIN, std::string(1, c));
        } catch (const pattern::PatternError &e) {
            parsed = false;
            EXPECT_EQ("bad_alphabet", e.code()) << c;
        }
        EXPECT_EQ(listed, parsed) << c;
    }
    std::vector<int> codes;
    for (const Json::Value &id : caps["genetic_codes"]) {
        codes.push_back(id.asInt());
    }
    EXPECT_EQ(std::vector<int>({ 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23,
                                 24, 25, 26, 27, 28, 29, 30, 31, 32, 33 }), codes);
    EXPECT_EQ(1, caps["default_genetic_code"].asInt());
    EXPECT_NE(std::string::npos, caps["protein_rule"].asString().find("stop_unsupported"));
    // owner decision #16: a graph with its mask counts exactly, and has no dummy fraction
    EXPECT_EQ("exact", caps["counting"].asString());
    EXPECT_TRUE(caps.isMember("dummy_fraction"));
    EXPECT_TRUE(caps["dummy_fraction"].isNull());
    // no prose of its own: the document a service's MCP tool returns has a 32 KiB ceiling
    EXPECT_FALSE(caps.isMember("counting_rule"));

    // loading: nothing about the graph is known yet
    caps = pattern_capabilities_json(nullptr, limits(), false);
    EXPECT_TRUE(caps["available"].isNull());
    EXPECT_TRUE(caps["graph_mode"].isNull());
    EXPECT_TRUE(caps["k"].isNull());
    EXPECT_TRUE(caps.isMember("counting"));
    EXPECT_TRUE(caps["counting"].isNull());
    EXPECT_TRUE(caps["dummy_fraction"].isNull());

    caps = pattern_capabilities_json(g.get(), limits(), true);
    EXPECT_EQ(3u, caps.size());
    EXPECT_FALSE(caps["available"].asBool());
    EXPECT_EQ("multi_graph_later_increment", caps["unavailable_reason"].asString());
}


// ---------------------------------------------------------------- without the mask (#16)

// |masked|'s graph without its dummy-edge mask (owner decision #16): the same BOSS written
// and loaded again without the mask, so the same node ids and rows; a PRIMARY graph wrapped in
// CanonicalDBG again (the wrapper reads at construction whether its graph has a mask); an
// annotation of as many rows, without labels (these tests read none)
std::unique_ptr<AnnotatedDBG> unmasked(const AnnotatedDBG &masked, const std::string &name) {
    const DeBruijnGraph &graph = masked.get_graph();
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    const auto *dbg = dynamic_cast<const DBGSuccinct*>(canonical ? &canonical->get_graph()
                                                                 : &graph);
    EXPECT_TRUE(dbg);
    const std::string base = test_dump_dir() + "/pattern_route_unmasked_" + name;
    dbg->serialize(base);
    auto loaded = std::make_shared<DBGSuccinct>(2);
    EXPECT_TRUE(loaded->load_without_mask(base + ".dbg"));
    EXPECT_EQ(nullptr, loaded->get_mask());
    std::shared_ptr<DeBruijnGraph> served = loaded;
    if (canonical)
        served = std::make_shared<CanonicalDBG>(served);
    return std::make_unique<AnnotatedDBG>(
            served, std::make_unique<annot::ColumnCompressed<>>(loaded->max_index()));
}

bool iupac_matches(char code, char base) {
    static const std::map<char, std::string> sets = {
        { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" }, { 'R', "AG" }, { 'Y', "CT" },
        { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" }, { 'B', "CGT" },
        { 'D', "AGT" }, { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
    };
    return sets.at(code).find(base) != std::string::npos;
}

std::string iupac_rc(const std::string &p) {
    static const std::map<char, char> c = {
        { 'A', 'T' }, { 'C', 'G' }, { 'G', 'C' }, { 'T', 'A' }, { 'R', 'Y' }, { 'Y', 'R' },
        { 'S', 'S' }, { 'W', 'W' }, { 'K', 'M' }, { 'M', 'K' }, { 'B', 'V' }, { 'V', 'B' },
        { 'D', 'H' }, { 'H', 'D' }, { 'N', 'N' },
    };
    std::string r(p.rbegin(), p.rend());
    for (char &x : r) {
        x = c.at(x);
    }
    return r;
}

/**
 * The graph-walk oracle of a succinct graph (BASIC or native CANONICAL), never the engine:
 * every edge of its BOSS spelled (its node's symbols and W), the edges with W = $ (sink
 * dummies, the main dummy edge) left out as no pattern base matches there; for every
 * orientation of |p| searched (P and rc(P), or P once when palindromic) and every offset of the
 * scope, the contexts (orientation key, edge, offset) whose k-mer holds the oriented pattern
 * there, each marked real or a source dummy ($ in its node).
 */
struct EdgeContext {
    std::string orientation;
    uint64_t edge;
    uint64_t offset;
    bool dummy;
    std::string kmer;
};

std::vector<EdgeContext> edge_oracle(const DBGSuccinct &graph, const std::string &p,
                                     bool strand_stated, bool suffix_only) {
    const auto &boss = graph.get_boss();
    const size_t k = graph.get_k();
    std::vector<std::pair<std::string, std::string>> oriented;
    if (iupac_rc(p) == p) {
        oriented.emplace_back(strand_stated ? "both" : "palindromic", p);
    } else {
        oriented.emplace_back(strand_stated ? "+" : "forward", p);
        oriented.emplace_back(strand_stated ? "-" : "reverse", iupac_rc(p));
    }
    std::vector<EdgeContext> out;
    for (uint64_t e = 1; e <= boss.num_edges(); ++e) {
        const auto w = boss.get_W(e) % boss.alph_size;
        if (!w)
            continue;
        const std::string kmer = boss.get_node_str(e) + boss.decode(w);
        const bool dummy = kmer.find('$') != std::string::npos;
        for (const auto &[key, q] : oriented) {
            for (size_t o = suffix_only ? k - q.size() : 0; o + q.size() <= k; ++o) {
                bool ok = true;
                for (size_t i = 0; i < q.size() && ok; ++i) {
                    ok = kmer[o + i] != '$' && iupac_matches(q[i], kmer[o + i]);
                }
                if (ok)
                    out.push_back({ key, e, o, dummy, kmer });
            }
        }
    }
    return out;
}

/**
 * A count of a graph without its mask against the oracle's |real| contexts and |all| (real and
 * source dummies): exact only as the real count (and exact 0 whenever nothing is a candidate:
 * an empty block), else bounds with the real count inside, upper the candidates, and the
 * estimate upper x f rounded into the bounds; nothing else
 */
void expect_unmasked_count(const Json::Value &c, uint64_t real, uint64_t all,
                           const DummyFraction &f, const std::string &where) {
    SCOPED_TRACE(where);
    ASSERT_TRUE(c.isObject()) << c;
    if (!all) {
        EXPECT_EQ("exact", c["relation"].asString()) << c;
        EXPECT_EQ(0u, c["value"].asUInt64());
    }
    if (c["relation"].asString() == "exact") {
        EXPECT_EQ(real, c["value"].asUInt64()) << c;
        EXPECT_FALSE(c.isMember("estimate")) << c;
        EXPECT_FALSE(c.isMember("upper")) << c;
        return;
    }
    ASSERT_EQ("bounds", c["relation"].asString()) << c;
    const uint64_t lower = c["lower"].asUInt64(), upper = c["upper"].asUInt64();
    EXPECT_EQ(lower, c["value"].asUInt64());
    EXPECT_LE(lower, real) << c;
    EXPECT_EQ(all, upper) << c;
    ASSERT_TRUE(c.isMember("estimate")) << c;
    const uint64_t e = std::max<uint64_t>(lower, std::min<uint64_t>(
            upper, std::llround(static_cast<double>(upper) * f.value)));
    EXPECT_EQ(e, c["estimate"].asUInt64()) << c;
    // the count's own fields (a total also carries suffix, by_offset and by_strand)
    std::vector<std::string> own;
    for (const std::string &name : c.getMemberNames()) {
        if (name != "suffix" && name != "by_offset" && name != "by_strand"
                && name != "by_orientation") {
            own.push_back(name);
        }
    }
    EXPECT_EQ((std::vector<std::string>{ "estimate", "lower", "relation", "unit", "upper",
                                         "value" }), own);
}

const DBGSuccinct& succinct(const AnnotatedDBG &anno_graph) {
    const DeBruijnGraph &graph = anno_graph.get_graph();
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// patterns of the tiny records: ACGT and TTGG begin a record (their source dummies hold them:
// $ACGTTG, $$ACGTT, $$$ACGT), the others inside one, IUPAC ones, and absent ones
const std::vector<std::string> kUnmaskedPatterns = {
    "ACGT", "TTGG", "AACG", "CGTA", "GCTT", "RCGT", "TTNG", "GGGGG", "CCCCC", "ACGTTGC",
};

// Owner decision #16: a graph without its dummy-edge mask is served. Its counts against the
// graph-walk oracle over the BOSS's edges: upper the candidates, source dummies included (never
// a sink: W = $ matches no base), lower <= the real count, exact only as the real count, exact 0
// where nothing is a candidate; the estimate upper x f; the answer and the capabilities say so
TEST(PatternRoute, UnmaskedCountsAgainstTheEdgeOracle) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        const bool stated = mode == DeBruijnGraph::BASIC;
        SCOPED_TRACE(stated ? "basic" : "canonical");
        auto masked = tiny(mode);
        auto g = unmasked(*masked, stated ? "basic" : "canonical");
        const DBGSuccinct &dbg = succinct(*g);
        const std::optional<DummyFraction> f = dummy_fraction(*g);
        ASSERT_TRUE(f);
        EXPECT_FALSE(dummy_fraction(*masked));
        EXPECT_EQ(pattern::kRealFractionSamples, f->samples);
        EXPECT_LT(f->value, 1.0);

        // the capabilities: served, the mask absent, counting upper_bound, its fraction
        const Json::Value caps = pattern_capabilities_json(g.get(), limits(), false);
        EXPECT_TRUE(caps["available"].asBool());
        EXPECT_TRUE(caps["unavailable_reason"].isNull());
        EXPECT_EQ("absent", caps["mask"].asString());
        EXPECT_EQ("upper_bound", caps["counting"].asString());
        EXPECT_EQ(dummy_fraction_json(*f), caps["dummy_fraction"]);
        EXPECT_EQ("sampled", caps["dummy_fraction"]["source"].asString());
        // the rest of the block is the masked graph's
        Json::Value a = caps, b = pattern_capabilities_json(masked.get(), limits(), false);
        for (Json::Value *x : { &a, &b }) {
            x->removeMember("mask");
            x->removeMember("counting");
            x->removeMember("dummy_fraction");
        }
        EXPECT_EQ(b, a);

        for (const char *scope : { "any_offset", "suffix" }) {
            const bool suffix = std::string(scope) == "suffix";
            std::string body = "{\"patterns\": [";
            for (size_t i = 0; i < kUnmaskedPatterns.size(); ++i) {
                body += std::string(i ? ", " : "") + "{\"iupac\": \"" + kUnmaskedPatterns[i]
                        + "\"}";
            }
            body += "], \"mode\": \"count\", \"scope\": \"" + std::string(scope) + "\"}";
            const Json::Value out = run(*g, body);
            const Json::Value exact = run(*masked, body);
            // the answer's index states the counting and the fraction; the masked one does not
            EXPECT_EQ("upper_bound", out["index"]["counting"].asString());
            EXPECT_EQ(dummy_fraction_json(*f), out["index"]["dummy_fraction"]);
            EXPECT_FALSE(exact["index"].isMember("counting"));
            EXPECT_FALSE(exact["index"].isMember("dummy_fraction"));
            bool some_dummy = false;
            for (size_t i = 0; i < kUnmaskedPatterns.size(); ++i) {
                const std::string &p = kUnmaskedPatterns[i];
                SCOPED_TRACE(p + " " + scope);
                const Json::Value &e = out["patterns"][Json::ArrayIndex(i)];
                const Json::Value &m = exact["patterns"][Json::ArrayIndex(i)];
                ASSERT_FALSE(e.isMember("error")) << e;
                const auto oracle = edge_oracle(dbg, p, stated, suffix);
                auto tally = [&](const std::string &orientation, int64_t offset) {
                    uint64_t real = 0, all = 0;
                    for (const EdgeContext &x : oracle) {
                        if ((orientation.empty() || x.orientation == orientation)
                                && (offset < 0 || x.offset == static_cast<uint64_t>(offset))) {
                            ++all;
                            real += !x.dummy;
                        }
                    }
                    return std::make_pair(real, all);
                };
                const auto [real, all] = tally("", -1);
                some_dummy |= real < all;
                // the masked graph's exact count is the oracle's real one
                EXPECT_EQ("exact", m["counts"]["contexts"]["relation"].asString());
                EXPECT_EQ(real, m["counts"]["contexts"]["value"].asUInt64());

                const Json::Value &c = e["counts"]["contexts"];
                expect_unmasked_count(c, real, all, *f, "total");
                const auto [sreal, sall] = tally("", kK - p.size());
                expect_unmasked_count(c["suffix"], sreal, sall, *f, "suffix");
                for (const std::string &o : c["by_offset"].getMemberNames()) {
                    const auto [r, a] = tally("", std::stoll(o));
                    expect_unmasked_count(c["by_offset"][o], r, a, *f, "offset " + o);
                }
                EXPECT_EQ(m["counts"]["contexts"]["by_offset"].getMemberNames(),
                          c["by_offset"].getMemberNames());
                const Json::Value &by = c[stated ? "by_strand" : "by_orientation"];
                EXPECT_EQ(m["counts"]["contexts"][stated ? "by_strand" : "by_orientation"]
                                  .getMemberNames(), by.getMemberNames());
                for (const std::string &o : by.getMemberNames()) {
                    const auto [r, a] = tally(o, -1);
                    expect_unmasked_count(by[o], r, a, *f, "orientation " + o);
                }
                // the note says what an estimate rests on, exactly when one is stated
                bool estimated = false;
                std::function<void(const Json::Value&)> find = [&](const Json::Value &v) {
                    if (v.isObject()) {
                        estimated |= v.isMember("estimate");
                        for (const auto &name : v.getMemberNames()) {
                            find(v[name]);
                        }
                    }
                };
                find(e["counts"]);
                const auto &notes = e["notes"];
                const bool noted = std::any_of(notes.begin(), notes.end(), [](const auto &n) {
                    return n.asString() == "estimate_sampled_dummy_fraction";
                });
                EXPECT_EQ(estimated, noted) << e;
                // the masked answer is written as before: no estimate, no note
                Json::Value mc = m["counts"];
                estimated = false;
                find(mc);
                EXPECT_FALSE(estimated) << m;
                for (const auto &n : m["notes"]) {
                    EXPECT_NE("estimate_sampled_dummy_fraction", n.asString());
                }
                // everything but the counts, the work and the notes is the masked answer's
                Json::Value x = e, y = m;
                for (Json::Value *v : { &x, &y }) {
                    for (const char *field : { "counts", "work", "notes", "timing" }) {
                        v->removeMember(field);
                    }
                }
                EXPECT_EQ(y, x);
                // labels and occurrences are not read in mode count: unknown as before
                EXPECT_EQ(m["counts"]["labels"], e["counts"]["labels"]);
                EXPECT_EQ(m["counts"]["occurrences"], e["counts"]["occurrences"]);
            }
            // the records' starts put patterns in source dummies: the test sees some
            EXPECT_TRUE(some_dummy);
        }
    }
}

// Owner decision #16, the retrieval modes without the mask: the lists are exact (every
// released context a real k-mer: a dummy is dropped), equal to the masked graph's lists; a
// complete release makes the counts exact; all_or_count admits on the upper bound, stated
TEST(PatternRoute, UnmaskedRetrievalListsAreExact) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        const std::string name = mode == DeBruijnGraph::BASIC ? "basic"
                               : mode == DeBruijnGraph::CANONICAL ? "canonical" : "primary";
        SCOPED_TRACE(name);
        auto masked = tiny(mode);
        auto g = unmasked(*masked, "retrieval_" + name);
        std::string patterns;
        for (size_t i = 0; i < kUnmaskedPatterns.size(); ++i) {
            patterns += std::string(i ? ", " : "") + "{\"iupac\": \"" + kUnmaskedPatterns[i]
                        + "\"}";
        }
        for (const char *m : { "all_or_count", "partial" }) {
            SCOPED_TRACE(m);
            const std::string body = "{\"patterns\": [" + patterns + "], \"mode\": \""
                                     + std::string(m) + "\"}";
            const Json::Value out = run(*g, body);
            const Json::Value exact = run(*masked, body);
            for (size_t i = 0; i < kUnmaskedPatterns.size(); ++i) {
                SCOPED_TRACE(kUnmaskedPatterns[i]);
                const Json::Value &e = out["patterns"][Json::ArrayIndex(i)];
                const Json::Value &x = exact["patterns"][Json::ArrayIndex(i)];
                // every context released: the masked graph's list, node and row ids alike, and
                // the counts exact, the masked graph's
                EXPECT_EQ(x["results"], e["results"]);
                EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
                EXPECT_EQ(x["counts"], e["counts"]);
                EXPECT_EQ(x["withheld"], e["withheld"]);
                EXPECT_EQ(x["cut"], e["cut"]);
                for (const Json::Value &r : e["results"]) {
                    EXPECT_EQ(std::string::npos, r["kmer"].asString().find('$')) << r;
                }
            }
        }
        // partial, cut: the first contexts of the masked graph's list
        const Json::Value cut = run(*g, "{\"patterns\": [" + patterns + "], \"mode\": "
                                        "\"partial\", \"max_contexts\": 2}");
        const Json::Value cut_exact = run(*masked, "{\"patterns\": [" + patterns + "], "
                                                   "\"mode\": \"partial\", \"max_contexts\": 2}");
        for (size_t i = 0; i < kUnmaskedPatterns.size(); ++i) {
            SCOPED_TRACE(kUnmaskedPatterns[i]);
            const Json::Value &e = cut["patterns"][Json::ArrayIndex(i)];
            const Json::Value &x = cut_exact["patterns"][Json::ArrayIndex(i)];
            EXPECT_EQ(x["results"], e["results"]);
            EXPECT_EQ(x["returned"], e["returned"]);
            if (x["cut"].isNull()) {
                EXPECT_TRUE(e["cut"].isNull()) << e;
            } else {
                EXPECT_EQ(x["cut"], e["cut"]);
                EXPECT_FALSE(e["retrieval_complete"].asBool());
                // a count of a list cut: a true bound, at least what was returned
                const Json::Value &c = e["counts"]["contexts"];
                EXPECT_GE(c["value"].asUInt64(), e["returned"].asUInt64()) << c;
                EXPECT_LE(c["value"].asUInt64(), x["counts"]["contexts"]["value"].asUInt64());
                if (c["relation"].asString() == "bounds") {
                    EXPECT_GE(c["upper"].asUInt64(),
                              x["counts"]["contexts"]["value"].asUInt64());
                }
            }
        }
    }

    // all_or_count admits on the upper bound (conservative): ACGT's real contexts fit
    // max_contexts, its candidates (3 source dummies among them) do not. The masked graph
    // releases them; without the mask they are withheld, the count still a true bound, and the
    // engine says why (threshold_upper_bound)
    auto masked = tiny();
    auto g = unmasked(*masked, "admission");
    const uint64_t real = oracle("ACGT").size();
    const std::string body = "{\"patterns\": [{\"dna\": \"ACGT\"}], \"max_contexts\": "
                             + std::to_string(real) + "}";
    const Json::Value x = run(*masked, body)["patterns"][0];
    EXPECT_TRUE(x["retrieval_complete"].asBool());
    EXPECT_EQ(real, x["returned"].asUInt64());
    const Json::Value e = run(*g, body)["patterns"][0];
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    EXPECT_EQ("count_above_threshold", e["withheld"]["reason"].asString()) << e;
    EXPECT_EQ(0u, e["returned"].asUInt64());
    const Json::Value &c = e["counts"]["contexts"];
    EXPECT_EQ("bounds", c["relation"].asString()) << c;
    EXPECT_EQ(real + 3, c["upper"].asUInt64()) << c;
    EXPECT_LE(c["lower"].asUInt64(), real);
    bool noted = false;
    for (const auto &n : e["notes"]) {
        noted |= n.asString() == "threshold_upper_bound";
    }
    EXPECT_TRUE(noted) << e;
    // with room for the candidates: released, every one real, the count exact
    const Json::Value ok = run(*g, "{\"patterns\": [{\"dna\": \"ACGT\"}], \"max_contexts\": "
                                   + std::to_string(real + 3) + "}")["patterns"][0];
    EXPECT_TRUE(ok["retrieval_complete"].asBool());
    EXPECT_EQ(real, ok["returned"].asUInt64());
    EXPECT_EQ("exact", ok["counts"]["contexts"]["relation"].asString());
    EXPECT_EQ(real, ok["counts"]["contexts"]["value"].asUInt64());
}

// Owner decision #16: the dummy fraction is sampled once per graph and kept (the loader samples
// it in its thread; a graph not loaded that way at its first use), the same for every request
// and every reader; a graph with its mask has none
TEST(PatternRoute, UnmaskedDummyFractionIsKeptPerGraph) {
    auto masked = tiny();
    auto g = unmasked(*masked, "kept");
    const auto a = dummy_fraction(*g);
    const auto b = dummy_fraction(*g);
    ASSERT_TRUE(a && b);
    EXPECT_EQ(a->real, b->real);
    EXPECT_EQ(a->value, b->value);
    const DummyFraction direct = pattern::sample_real_fraction(succinct(*g));
    EXPECT_EQ(direct.real, a->real);
    EXPECT_EQ(direct.samples, a->samples);
    EXPECT_EQ(direct.lower, a->lower);
    EXPECT_EQ(direct.upper, a->upper);
    // the same graph loaded again: the same draws
    auto again = unmasked(*masked, "kept_again");
    EXPECT_EQ(a->real, dummy_fraction(*again)->real);
    EXPECT_FALSE(dummy_fraction(*masked));
}

} // namespace
