#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// Bounded mismatches on POST /pattern (SPEC §21): the request field max_mismatches and its
// refusals, the per-kind caps and the slot error mismatches_above_cap, the answer's fields per
// level (by_mismatches, steps_by_mismatches, stop.mismatch_level), each result's level and
// substitutions in the submitted pattern's coordinates, the floor slot's
// max_mismatches_within_floor, the capabilities' mismatches object, and the identity of a
// request with max_mismatches 0 with one without the field. Expected contexts and levels come
// from a brute force over the records' k-mers with the test's own IUPAC table, never from the
// engine.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;

// a floor low enough for the short patterns of tiny graphs, and the caps as configured
PatternLimits limits(uint64_t max_mismatches = 3, uint64_t max_mismatches_protein = 3) {
    PatternLimits l;
    l.min_information_bits = 4;
    l.max_mismatches = max_mismatches;
    l.max_mismatches_protein = max_mismatches_protein;
    return l;
}

std::unique_ptr<AnnotatedDBG> build(size_t k, const std::vector<std::string> &records) {
    std::vector<std::string> labels;
    for (size_t i = 0; i < records.size(); ++i) {
        labels.push_back("r" + std::to_string(i));
    }
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, records, labels, DeBruijnGraph::BASIC);
}

Json::Value run(const AnnotatedDBG &anno_graph, const std::string &body,
                const PatternLimits &caps = limits()) {
    return process_pattern_request(parse_pattern_body(body), anno_graph, caps, "rel",
                                   nullptr, nullptr, nullptr, nullptr);
}

// the refusal of |body|: {status, code}, and its message in |error|; {200, ""} when answered
std::pair<int, std::string> refusal(const AnnotatedDBG &anno_graph, const std::string &body,
                                    const PatternLimits &caps = limits(),
                                    std::string *error = nullptr) {
    try {
        run(anno_graph, body, caps);
    } catch (const PatternRefusal &e) {
        if (error)
            *error = e.what();
        return { e.status(), e.code() };
    }
    return { 200, "" };
}

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

// the answer without its timings (which vary between runs)
Json::Value untimed(Json::Value v) {
    v.removeMember("timing");
    for (Json::Value &e : v["patterns"]) {
        e.removeMember("timing");
    }
    return v;
}

// ---------------------------------------------------------------- the oracle

// the test's own IUPAC table and complement (never the engine's)
const std::map<char, std::string> kIupac = {
    { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" }, { 'R', "AG" }, { 'Y', "CT" },
    { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" }, { 'B', "CGT" }, { 'D', "AGT" },
    { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
};

char complement(char c) {
    const std::string from = "ACGTRYSWKMBDHVN";
    const std::string to = "TGCAYRSWMKVHDBN";
    return to[from.find(c)];
}

std::string rc(const std::string &s) {
    std::string out(s.rbegin(), s.rend());
    for (char &c : out) {
        c = complement(c);
    }
    return out;
}

std::string random_dna(std::mt19937 &rng, size_t n) {
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s.push_back("ACGT"[rng() % 4]);
    }
    return s;
}

// |s| with the base at |position| replaced by another one
std::string substituted(std::string s, size_t position, std::mt19937 &rng) {
    char c = s[position];
    while (c == s[position]) {
        c = "ACGT"[rng() % 4];
    }
    s[position] = c;
    return s;
}

// one substitution of an instance, in the submitted pattern's coordinates (SPEC §21.8)
using Sub = std::tuple<uint32_t, char, char>;  // position, expected, found

// a graph context of the pattern within m substitutions: (strand, k-mer, offset, level, subs)
using Ctx = std::tuple<std::string, std::string, uint64_t, uint32_t, std::vector<Sub>>;

// every context of |P| within |m| substitutions on a BASIC graph of |records| (every k-mer as
// deposited), scope any_offset, both strands
std::set<Ctx> oracle(const std::vector<std::string> &records, size_t k, const std::string &P,
                     uint32_t m) {
    std::set<std::string> kmers;
    for (const std::string &r : records) {
        for (size_t i = 0; i + k <= r.size(); ++i) {
            kmers.insert(r.substr(i, k));
        }
    }
    const size_t L = P.size();
    std::vector<std::pair<std::string, std::string>> oriented;
    if (rc(P) == P) {
        oriented.emplace_back("=", P);
    } else {
        oriented.emplace_back("+", P);
        oriented.emplace_back("-", rc(P));
    }
    std::set<Ctx> out;
    for (const auto &[strand, q] : oriented) {
        for (const std::string &kmer : kmers) {
            for (size_t o = 0; o + L <= k; ++o) {
                std::vector<Sub> subs;
                for (size_t i = 0; i < L; ++i) {
                    const char base = kmer[o + i];
                    if (kIupac.at(q[i]).find(base) != std::string::npos)
                        continue;
                    if (strand == "-") {
                        subs.emplace_back(L - 1 - i, P[L - 1 - i], complement(base));
                    } else {
                        subs.emplace_back(i, P[i], base);
                    }
                }
                if (subs.size() > m)
                    continue;
                std::sort(subs.begin(), subs.end());
                out.emplace(strand, kmer, o, subs.size(), subs);
            }
        }
    }
    return out;
}

// the answer's results as the oracle states them
std::set<Ctx> contexts_of(const Json::Value &entry) {
    std::set<Ctx> out;
    for (const Json::Value &r : entry["results"]) {
        std::vector<Sub> subs;
        for (const Json::Value &s : r["substitutions"]) {
            EXPECT_EQ(3u, s.size());
            subs.emplace_back(s["position"].asUInt(), s["expected"].asString().at(0),
                              s["found"].asString().at(0));
        }
        EXPECT_EQ(subs.size(), r["mismatches"].asUInt());
        out.emplace(r["strand"].asString(), r["kmer"].asString(), r["offset"].asUInt64(),
                    r["mismatches"].asUInt(), subs);
    }
    return out;
}

size_t at_level(const std::set<Ctx> &contexts, uint32_t level) {
    return std::count_if(contexts.begin(), contexts.end(),
                         [&](const Ctx &c) { return std::get<3>(c) == level; });
}

// an exact count of |value|
void expect_exact(const Json::Value &c, uint64_t value, const std::string &what) {
    EXPECT_EQ("exact", c["relation"].asString()) << what;
    EXPECT_EQ(value, c["value"].asUInt64()) << what;
}

// the records of the DNA cases: a random 70-mer R, its copies with one, two and three
// substitutions planted inside P = R[20, 51) (the 31-mer, one offset at k = 31), and the
// reverse complement of the one-substitution copy (a context on the other strand)
struct Planted {
    std::string P;
    std::vector<std::string> records;
};

Planted planted(std::mt19937 &rng) {
    const std::string R = random_dna(rng, 70);
    const std::string one = substituted(R, 25, rng);
    const std::string two = substituted(substituted(R, 23, rng), 48, rng);
    const std::string three = substituted(substituted(substituted(R, 21, rng), 35, rng), 50, rng);
    return { R.substr(20, 31), { R, one, two, three, rc(one) } };
}

const size_t kK = 31;

std::string request(const std::string &patterns, const std::string &rest = "") {
    return "{\"patterns\": [" + patterns + "], \"output\": {\"labels\": \"none\"}" + rest + "}";
}


// ---------------------------------------------------------------- the field

TEST(PatternMismatchRoute, FieldParsedAndRefused) {
    std::mt19937 rng(21);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    const std::string p = "{\"dna\": \"" + c.P + "\"}";
    const std::vector<std::tuple<std::string, int, std::string>> cases = {
        // the values (SPEC §21.2): 0 to 3, an integer; a fraction, a string, null, a boolean
        // and a negative or larger integer are invalid_request
        { request(p, ", \"max_mismatches\": 0"), 200, "" },
        { request(p, ", \"max_mismatches\": 1"), 200, "" },
        { request(p, ", \"max_mismatches\": 3"), 200, "" },
        { request(p, ", \"max_mismatches\": 3.0"), 200, "" },
        { request(p, ", \"max_mismatches\": 4"), 400, "invalid_request" },
        { request(p, ", \"max_mismatches\": -1"), 400, "invalid_request" },
        { request(p, ", \"max_mismatches\": 1.5"), 400, "invalid_request" },
        { request(p, ", \"max_mismatches\": \"1\""), 400, "invalid_request" },
        { request(p, ", \"max_mismatches\": null"), 400, "invalid_request" },
        { request(p, ", \"max_mismatches\": true"), 400, "invalid_request" },
        // with every mode and scope, and with long_search "paths" (a pattern of at most k
        // bases is answered alike)
        { request(p, ", \"max_mismatches\": 2, \"mode\": \"count\""), 200, "" },
        { request(p, ", \"max_mismatches\": 2, \"mode\": \"partial\""), 200, "" },
        { request(p, ", \"max_mismatches\": 2, \"scope\": \"suffix\""), 200, "" },
        { request(p, ", \"max_mismatches\": 2, \"long_search\": \"paths\""), 200, "" },
        { request(p, ", \"max_mismatches\": 2, \"strands\": \"reverse\""), 200, "" },
        // what is not served with substitutions yet (SPEC §21.16): the labels, a predicate,
        // the supported paths -- told what to wait for, before the annotation is looked at
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 1, \"output\": {\"labels\": "
          "\"all\"}}", 400, "later_increment" },
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 1, \"mode\": \"count\", \"output\": "
          "{\"labels\": \"all\"}}", 400, "later_increment" },
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 1, \"predicate\": {\"any\": "
          "[\"r0\"]}}", 400, "later_increment" },
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 1, \"predicate\": {\"any\": "
          "[\"r0\"]}, \"output\": {\"labels\": \"predicate_only\"}}", 400, "later_increment" },
        { request(p, ", \"max_mismatches\": 1, \"long_search\": \"supported_paths\""), 400,
          "later_increment" },
        // with 0 the same combinations are answered (or refused) as without the field
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 0, \"predicate\": {\"any\": "
          "[\"r0\"]}}", 400, "annotation_unbudgeted" },
        { "{\"patterns\": [" + p + "], \"max_mismatches\": 0, \"mode\": \"count\", \"output\": "
          "{\"labels\": \"all\"}}", 200, "" },
    };
    for (const auto &[body, status, code] : cases) {
        EXPECT_EQ(std::make_pair(status, code), refusal(*g, body)) << body;
    }
    std::string error;
    refusal(*g, request(p, ", \"max_mismatches\": 4"), limits(), &error);
    EXPECT_EQ("request.max_mismatches: expected an integer from 0 to 3 (substitutions per "
              "instance)", error);
    refusal(*g, request(p, ", \"max_mismatches\": 1, \"long_search\": \"supported_paths\""),
            limits(), &error);
    EXPECT_NE(std::string::npos, error.find("mismatches.later_increment")) << error;
    // the order of the checks (SPEC §5 step 8): max_mismatches after predicate_strands, before
    // predicate_scope; the later-increment combinations after every single field
    refusal(*g, request(p, ", \"max_mismatches\": 4, \"predicate_scope\": \"x\""), limits(),
            &error);
    EXPECT_EQ(0u, error.find("request.max_mismatches")) << error;
    refusal(*g, request(p, ", \"max_mismatches\": 4, \"predicate_strands\": \"x\""), limits(),
            &error);
    EXPECT_EQ(0u, error.find("request.predicate_strands")) << error;
    refusal(*g, "{\"patterns\": [" + p + "], \"max_mismatches\": 1, \"output\": {\"labels\": "
                "\"all\"}, \"predicate_scope\": \"motif\"}", limits(), &error);
    EXPECT_EQ(0u, error.find("request.predicate_scope")) << error;
}

TEST(PatternMismatchRoute, ZeroIsTheFieldsAbsence) {
    std::mt19937 rng(22);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    const std::string p = "{\"dna\": \"" + c.P + "\"}, {\"iupac\": \"" + c.P.substr(3, 20)
            + "N\"}, {\"protein\": \"MK\"}, {\"dna\": \"" + c.P + "ACGT\"}";
    for (const std::string &rest : { std::string(), std::string(", \"mode\": \"count\""),
                                     std::string(", \"mode\": \"partial\", \"max_steps\": 20") }) {
        const Json::Value without = untimed(run(*g, request(p, rest)));
        const Json::Value with = untimed(run(*g, request(p, rest + ", \"max_mismatches\": 0")));
        EXPECT_EQ(compact(without), compact(with)) << rest;
        // nothing of §21 in the answer
        EXPECT_EQ(std::string::npos, compact(with).find("mismatch")) << rest;
    }
}


// ---------------------------------------------------------------- the caps

TEST(PatternMismatchRoute, CapsPerKindAndTheSlotError) {
    std::mt19937 rng(23);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    const std::string dna = "{\"dna\": \"" + c.P + "\", \"id\": \"dna\"}";
    const std::string pep = "{\"protein\": \"MKLV\", \"id\": \"pep\"}";
    const std::string lng = "{\"dna\": \"" + c.P + "ACGT\", \"id\": \"long\"}";
    // a server with --pattern-max-mismatches 2 --pattern-max-mismatches-protein 1
    const PatternLimits caps = limits(2, 1);
    Json::Value a = run(*g, request(dna + ", " + pep + ", " + lng, ", \"max_mismatches\": 2"),
                        caps);
    ASSERT_EQ(3u, a["patterns"].size());
    EXPECT_EQ(2u, a["limits"]["max_mismatches"].asUInt64());
    // the dna pattern within its cap: answered, with the fields of §21.8
    const Json::Value &answered = a["patterns"][0];
    EXPECT_FALSE(answered.isMember("error"));
    EXPECT_EQ(2u, answered["max_mismatches"].asUInt64());
    EXPECT_NEAR(62 - std::log2(1 + 93 + 4185), answered["mismatch_information_bits"].asDouble(),
                1e-9);
    EXPECT_TRUE(answered["min_anchor_mismatch_information_bits"].isNull());
    EXPECT_FALSE(answered.isMember("max_mismatches_cap"));
    EXPECT_EQ(3u, answered["counts"]["contexts"]["by_mismatches"].size());
    // the peptide above the protein cap: its slot, the description as the engine would state
    // it, the cap applied, no counts; nothing charged
    const Json::Value &slot = a["patterns"][1];
    ASSERT_TRUE(slot.isMember("error"));
    EXPECT_EQ("mismatches_above_cap", slot["error"]["code"].asString());
    EXPECT_NE(std::string::npos, slot["error"]["message"].asString().find("mismatches.max.protein"));
    EXPECT_EQ(1u, slot["max_mismatches_cap"].asUInt64());
    EXPECT_EQ(2u, slot["max_mismatches"].asUInt64());
    EXPECT_EQ("protein", slot["kind"].asString());
    EXPECT_EQ("MKLV", slot["pattern"].asString());
    EXPECT_EQ(12u, slot["length"].asUInt64());
    EXPECT_EQ(4u, slot["residues"].asUInt64());
    EXPECT_EQ(1, slot["genetic_code"].asInt());
    EXPECT_TRUE(slot["mismatch_information_bits"].isDouble());
    EXPECT_LT(slot["mismatch_information_bits"].asDouble(), slot["information_bits"].asDouble());
    EXPECT_TRUE(slot["anchor_information_bits"].isNull());
    EXPECT_TRUE(slot["min_anchor_information_bits"].isNull());
    EXPECT_TRUE(slot["min_anchor_mismatch_information_bits"].isNull());
    EXPECT_FALSE(slot.isMember("counts"));
    EXPECT_FALSE(slot.isMember("work"));
    EXPECT_FALSE(slot.isMember("max_mismatches_within_floor"));
    std::set<std::string> fields;
    for (const std::string &f : slot.getMemberNames()) {
        fields.insert(f);
    }
    EXPECT_EQ(std::set<std::string>({ "id", "kind", "pattern", "length", "residues",
                                      "genetic_code", "information_bits",
                                      "anchor_information_bits", "min_anchor_information_bits",
                                      "max_mismatches", "mismatch_information_bits",
                                      "min_anchor_mismatch_information_bits",
                                      "max_mismatches_cap", "error" }), fields);
    // a pattern longer than k: the cap applied is 0 (a later increment), its anchor windows'
    // bits stated as for every long pattern (both orientations searched: the least of P[0, k)
    // and rc(P)[0, k), here the same exact 31-mers' bits)
    const Json::Value &lslot = a["patterns"][2];
    ASSERT_TRUE(lslot.isMember("error"));
    EXPECT_EQ("mismatches_above_cap", lslot["error"]["code"].asString());
    EXPECT_NE(std::string::npos, lslot["error"]["message"].asString().find("longer than k"));
    EXPECT_EQ(0u, lslot["max_mismatches_cap"].asUInt64());
    EXPECT_EQ(62.0, lslot["anchor_information_bits"].asDouble());
    EXPECT_EQ(62.0, lslot["min_anchor_information_bits"].asDouble());
    EXPECT_EQ(70.0, lslot["information_bits"].asDouble());
    EXPECT_NEAR(62 - std::log2(1 + 93 + 4185),
                lslot["min_anchor_mismatch_information_bits"].asDouble(), 1e-9);
    EXPECT_NEAR(70 - std::log2(1 + 3 * 35 + 9 * 35 * 34 / 2),
                lslot["mismatch_information_bits"].asDouble(), 1e-9);
    EXPECT_FALSE(lslot.isMember("residues"));
    // the dna pattern above the general cap
    a = run(*g, request(dna, ", \"max_mismatches\": 3"), caps);
    EXPECT_EQ("mismatches_above_cap", a["patterns"][0]["error"]["code"].asString());
    EXPECT_EQ(2u, a["patterns"][0]["max_mismatches_cap"].asUInt64());
    EXPECT_NE(std::string::npos,
              a["patterns"][0]["error"]["message"].asString().find("mismatches.max.dna"));
    // a kind whose cap is 0: not served with substitutions on this server
    a = run(*g, request(pep + ", " + dna, ", \"max_mismatches\": 1"), limits(3, 0));
    EXPECT_EQ("mismatches_above_cap", a["patterns"][0]["error"]["code"].asString());
    EXPECT_EQ(0u, a["patterns"][0]["max_mismatches_cap"].asUInt64());
    EXPECT_NE(std::string::npos,
              a["patterns"][0]["error"]["message"].asString().find("not served for protein"));
    EXPECT_FALSE(a["patterns"][1].isMember("error"));
    // a request with 0 is never above a cap, whatever the caps
    a = run(*g, request(dna + ", " + pep + ", " + lng, ", \"max_mismatches\": 0"), limits(0, 0));
    for (const Json::Value &e : a["patterns"]) {
        EXPECT_FALSE(e.isMember("error"));
    }
    // the slot of a peptide refused for its alphabet stays what it is (id, kind, error)
    a = run(*g, request("{\"protein\": \"MKLU\"}", ", \"max_mismatches\": 1"));
    EXPECT_EQ("bad_alphabet", a["patterns"][0]["error"]["code"].asString());
    EXPECT_EQ(3u, a["patterns"][0].size());
}


// ---------------------------------------------------------------- the answer

TEST(PatternMismatchRoute, LevelsAndSubstitutionsAgainstTheOracle) {
    std::mt19937 rng(24);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    for (uint32_t m : { 1u, 2u, 3u }) {
        const std::set<Ctx> expected = oracle(c.records, kK, c.P, m);
        ASSERT_GE(at_level(expected, 1), 2u);
        ASSERT_GE(at_level(expected, std::min(m, 2u)), 1u);
        const Json::Value a = run(*g, request("{\"dna\": \"" + c.P + "\"}",
                                              ", \"max_mismatches\": " + std::to_string(m)));
        const Json::Value &e = a["patterns"][0];
        ASSERT_FALSE(e.isMember("error")) << e;
        EXPECT_EQ(m, a["limits"]["max_mismatches"].asUInt());
        EXPECT_TRUE(e["retrieval_complete"].asBool());
        EXPECT_TRUE(e["stop"].isNull());
        // the counts per level, exact, keys "0" .. "m", zeros included; their sum the total
        const Json::Value &counts = e["counts"]["contexts"];
        const Json::Value &by_level = counts["by_mismatches"];
        ASSERT_EQ(m + 1, by_level.size());
        uint64_t sum = 0;
        for (uint32_t t = 0; t <= m; ++t) {
            const Json::Value &level = by_level[std::to_string(t)];
            expect_exact(level, at_level(expected, t), "level " + std::to_string(t));
            EXPECT_EQ("graph_contexts", level["unit"].asString());
            sum += level["value"].asUInt64();
        }
        expect_exact(counts, expected.size(), "total");
        EXPECT_EQ(expected.size(), sum);
        expect_exact(counts["by_offset"]["0"], expected.size(), "offset 0");
        // the results: every context within m, each with its level and its substitutions in
        // P's coordinates (the '-' strand's complemented, positions mirrored)
        EXPECT_EQ(expected, contexts_of(e));
        EXPECT_EQ(expected.size(), e["returned"].asUInt64());
        // in (mismatches, node, offset, orientation) order
        std::vector<std::tuple<uint32_t, uint64_t, uint64_t, int>> order;
        for (const Json::Value &r : e["results"]) {
            order.emplace_back(r["mismatches"].asUInt(), r["node"].asUInt64(),
                               r["offset"].asUInt64(), r["strand"].asString() == "+" ? 0 : 1);
        }
        EXPECT_TRUE(std::is_sorted(order.begin(), order.end()));
        // the work per level sums to the steps; a bits field per the ball
        const Json::Value &steps = e["work"]["steps_by_mismatches"];
        ASSERT_EQ(m + 1, steps.size());
        uint64_t steps_sum = 0;
        for (uint32_t t = 0; t <= m; ++t) {
            steps_sum += steps[std::to_string(t)].asUInt64();
        }
        EXPECT_EQ(e["work"]["steps"].asUInt64(), steps_sum);
        EXPECT_EQ(m, e["max_mismatches"].asUInt());
        EXPECT_LT(e["mismatch_information_bits"].asDouble(), e["information_bits"].asDouble());
        EXPECT_TRUE(e["min_anchor_mismatch_information_bits"].isNull());
        // mode count: the same counts, no results
        const Json::Value cnt = run(*g, request("{\"dna\": \"" + c.P + "\"}",
                                                ", \"max_mismatches\": " + std::to_string(m)
                                                + ", \"mode\": \"count\""))["patterns"][0];
        EXPECT_EQ(compact(by_level), compact(cnt["counts"]["contexts"]["by_mismatches"]));
        EXPECT_FALSE(cnt.isMember("results"));
    }
    // an IUPAC pattern: a substitution is a base outside the position's set, never at N
    std::string iupac = c.P.substr(0, 24);
    iupac[2] = kIupac.at(iupac[2]) == "A" || kIupac.at(iupac[2]) == "G" ? 'R' : 'Y';
    iupac[10] = 'N';
    iupac[20] = kIupac.at(iupac[20]) == "C" || kIupac.at(iupac[20]) == "G" ? 'S' : 'W';
    const std::set<Ctx> expected = oracle(c.records, kK, iupac, 2);
    ASSERT_GE(expected.size(), 1u);
    const Json::Value e = run(*g, request("{\"iupac\": \"" + iupac + "\"}",
                                          ", \"max_mismatches\": 2"))["patterns"][0];
    ASSERT_FALSE(e.isMember("error")) << e;
    EXPECT_EQ(expected, contexts_of(e));
    for (const Json::Value &r : e["results"]) {
        for (const Json::Value &s : r["substitutions"]) {
            EXPECT_NE('N', s["expected"].asString().at(0));
            EXPECT_EQ(iupac[s["position"].asUInt()], s["expected"].asString().at(0));
        }
    }
    for (uint32_t t = 0; t <= 2; ++t) {
        expect_exact(e["counts"]["contexts"]["by_mismatches"][std::to_string(t)],
                     at_level(expected, t), "iupac level " + std::to_string(t));
    }
}

TEST(PatternMismatchRoute, PeptideResidueSubstitutions) {
    // MKLV in the standard code: ATG AAR CTN|TTR GTN. Records: an instance, a synonymous one
    // (level 0), K replaced by R (AGA, a sense codon of another residue: level 1), K replaced
    // by a stop (TAA: no instance at any level), the reverse complement of the substituted
    // record (a '-' context at level 1), and V replaced by A beside the K substitution (level 2)
    const std::string flank_l = "CCGTA", flank_r = "GGATC";
    const std::string exact = flank_l + "ATGAAACTGGTT" + flank_r;
    const std::string synonymous = flank_l + "ATGAAGCTGGTT" + flank_r;
    const std::string one = flank_l + "ATGAGACTGGTT" + flank_r;
    const std::string stop = flank_l + "ATGTAACTGGTT" + flank_r;
    const std::string two = flank_l + "ATGAGACTGGCT" + flank_r;
    const std::vector<std::string> records = { exact, synonymous, one, stop, rc(one), two };
    const size_t k = 17;
    auto g = build(k, records);
    const std::string body = "{\"patterns\": [{\"protein\": \"MKLV\"}], \"output\": "
                             "{\"labels\": \"none\"}, \"max_mismatches\": ";
    Json::Value e = run(*g, body + "1}")["patterns"][0];
    ASSERT_FALSE(e.isMember("error")) << e;
    const Json::Value &by_level = e["counts"]["contexts"]["by_mismatches"];
    ASSERT_EQ(2u, by_level.size());
    // every k-mer holding an instance (k = 17 over 22-base records: the instance of 12 bases
    // sits at offsets 0 .. 5 of 6 k-mers per record) is a context
    expect_exact(by_level["0"], 2 * 6, "level 0");
    expect_exact(by_level["1"], 2 * 6, "level 1");
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    std::set<std::tuple<std::string, uint32_t, std::string, std::string, std::string>> subs;
    for (const Json::Value &r : e["results"]) {
        const Json::Value &list = r["substitutions"];
        ASSERT_EQ(r["mismatches"].asUInt(), list.size());
        for (const Json::Value &s : list) {
            EXPECT_EQ(4u, s.size());
            subs.emplace(r["strand"].asString(), s["residue"].asUInt(), s["expected"].asString(),
                         s["found"].asString(), s["codon"].asString());
        }
    }
    // one substituted residue, on both strands: residue 1 (K), found R, codon AGA in P's
    // reading direction
    EXPECT_EQ((std::set<std::tuple<std::string, uint32_t, std::string, std::string,
                                   std::string>>({ { "+", 1, "K", "R", "AGA" },
                                                    { "-", 1, "K", "R", "AGA" } })), subs);
    // the stop codon's record is no instance at any level; the two-substitution record is
    // found at level 2 only
    e = run(*g, body + "2}")["patterns"][0];
    expect_exact(e["counts"]["contexts"]["by_mismatches"]["0"], 12, "level 0 at m 2");
    expect_exact(e["counts"]["contexts"]["by_mismatches"]["1"], 12, "level 1 at m 2");
    expect_exact(e["counts"]["contexts"]["by_mismatches"]["2"], 6, "level 2 at m 2");
    for (const Json::Value &r : e["results"]) {
        // a stop codon substitutes no residue: the record holding TAA is no instance
        EXPECT_EQ(std::string::npos, r["instance"].asString().find("TAA")) << r["instance"];
        if (r["mismatches"].asUInt() == 2) {
            EXPECT_EQ("V", r["substitutions"][1]["expected"].asString());
            EXPECT_EQ("A", r["substitutions"][1]["found"].asString());
            EXPECT_EQ("GCT", r["substitutions"][1]["codon"].asString());
            EXPECT_EQ(3u, r["substitutions"][1]["residue"].asUInt());
        }
    }
    // the bits left: per residue the own codons and the sense substitutes
    EXPECT_LT(e["mismatch_information_bits"].asDouble(), e["information_bits"].asDouble());
    EXPECT_EQ(2u, e["max_mismatches"].asUInt());
}

TEST(PatternMismatchRoute, AStopLeavesTheLowerLevelsComplete) {
    std::mt19937 rng(25);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    const std::string p = "{\"dna\": \"" + c.P + "\"}";
    const Json::Value full = run(*g, request(p, ", \"max_mismatches\": 2"))["patterns"][0];
    ASSERT_TRUE(full["stop"].isNull());
    const Json::Value &steps = full["work"]["steps_by_mismatches"];
    const uint64_t s0 = steps["0"].asUInt64(), s1 = steps["1"].asUInt64(),
                   s2 = steps["2"].asUInt64();
    ASSERT_GT(s1, 1u);
    ASSERT_GT(s2, 1u);
    const std::set<Ctx> expected = oracle(c.records, kK, c.P, 2);
    // a budget spent inside the work of level 2: levels 0 and 1 exact, level 2 at_least, the
    // total at_least, stop {discovery, max_steps, mismatch_level 2}
    for (const std::string &mode : { std::string("all_or_count"), std::string("partial"),
                                     std::string("count") }) {
        const Json::Value e = run(*g, request(p, ", \"max_mismatches\": 2, \"max_steps\": "
                                              + std::to_string(s0 + s1 + 1) + ", \"mode\": \""
                                              + mode + "\""))["patterns"][0];
        ASSERT_FALSE(e.isMember("error")) << e;
        const Json::Value &by_level = e["counts"]["contexts"]["by_mismatches"];
        expect_exact(by_level["0"], at_level(expected, 0), mode + " level 0");
        expect_exact(by_level["1"], at_level(expected, 1), mode + " level 1");
        EXPECT_EQ("at_least", by_level["2"]["relation"].asString()) << mode;
        EXPECT_LE(by_level["2"]["value"].asUInt64(), at_level(expected, 2)) << mode;
        EXPECT_EQ("at_least", e["counts"]["contexts"]["relation"].asString()) << mode;
        EXPECT_EQ("max_steps", e["stop"]["reason"].asString()) << mode;
        EXPECT_EQ("discovery", e["stop"]["phase"].asString()) << mode;
        EXPECT_EQ(2u, e["stop"]["mismatch_level"].asUInt()) << mode;
        EXPECT_EQ(3u, e["stop"].size()) << mode;
        EXPECT_EQ(s0 + s1 + 1, e["work"]["steps"].asUInt64()) << mode;
        EXPECT_EQ(s0, e["work"]["steps_by_mismatches"]["0"].asUInt64()) << mode;
        EXPECT_EQ(s1, e["work"]["steps_by_mismatches"]["1"].asUInt64()) << mode;
        EXPECT_EQ(1u, e["work"]["steps_by_mismatches"]["2"].asUInt64()) << mode;
        if (mode == "partial") {
            // the levels below the stop released whole, the closest first
            const std::set<Ctx> got = contexts_of(e);
            for (const Ctx &x : expected) {
                if (std::get<3>(x) < 2) {
                    EXPECT_TRUE(got.count(x)) << std::get<1>(x);
                }
            }
            EXPECT_EQ("max_steps", e["cut"]["reason"].asString());
        } else if (mode == "all_or_count") {
            EXPECT_EQ("discovery_budget", e["withheld"]["reason"].asString());
        }
    }
    // a budget spent inside the work of level 1: level 2 never started (unknown)
    const Json::Value e = run(*g, request(p, ", \"max_mismatches\": 2, \"max_steps\": "
                                          + std::to_string(s0 + 1)))["patterns"][0];
    const Json::Value &by_level = e["counts"]["contexts"]["by_mismatches"];
    expect_exact(by_level["0"], at_level(expected, 0), "level 0");
    EXPECT_EQ("at_least", by_level["1"]["relation"].asString());
    EXPECT_EQ("unknown", by_level["2"]["relation"].asString());
    EXPECT_TRUE(by_level["2"]["value"].isNull());
    EXPECT_EQ(1u, e["stop"]["mismatch_level"].asUInt());
    EXPECT_EQ(0u, e["work"]["steps_by_mismatches"]["2"].asUInt64());
    // without mismatches a stop names no level
    const Json::Value plain = run(*g, request(p, ", \"max_steps\": 1"))["patterns"][0];
    EXPECT_EQ(2u, plain["stop"].size());
    EXPECT_FALSE(plain["stop"].isMember("mismatch_level"));
}

TEST(PatternMismatchRoute, TheFloorGatesTheBitsLeft) {
    std::mt19937 rng(26);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    PatternLimits floor = limits();
    floor.min_information_bits = 24;
    // an exact 16-mer: 32 bits, 26.4 within one substitution, 21.9 within two: refused at 2,
    // the slot stating 1 as the largest max_mismatches the floor admits
    const std::string p16 = "{\"dna\": \"" + c.P.substr(0, 16) + "\", \"id\": \"p16\"}";
    const std::string p10 = "{\"dna\": \"" + c.P.substr(0, 10) + "\", \"id\": \"p10\"}";
    Json::Value a = run(*g, request(p16 + ", " + p10, ", \"max_mismatches\": 2"), floor);
    const Json::Value &slot = a["patterns"][0];
    ASSERT_TRUE(slot.isMember("error")) << slot;
    EXPECT_EQ("information_below_floor", slot["error"]["code"].asString());
    EXPECT_EQ(32.0, slot["information_bits"].asDouble());
    EXPECT_NEAR(32 - std::log2(1 + 48 + 1080), slot["mismatch_information_bits"].asDouble(), 1e-9);
    EXPECT_EQ(1u, slot["max_mismatches_within_floor"].asUInt());
    EXPECT_EQ(2u, slot["max_mismatches"].asUInt());
    EXPECT_FALSE(slot.isMember("max_mismatches_cap"));
    EXPECT_NE(std::string::npos, slot["error"]["message"].asString().find("2 substitutions"));
    // 20 bits: below the floor at every max_mismatches, 0 included (null)
    const Json::Value &slot10 = a["patterns"][1];
    EXPECT_EQ("information_below_floor", slot10["error"]["code"].asString());
    EXPECT_TRUE(slot10.isMember("max_mismatches_within_floor"));
    EXPECT_TRUE(slot10["max_mismatches_within_floor"].isNull());
    // at 1 the 16-mer is served; a floor slot without mismatches has none of the fields
    a = run(*g, request(p16 + ", " + p10, ", \"max_mismatches\": 1"), floor);
    EXPECT_FALSE(a["patterns"][0].isMember("error"));
    a = run(*g, request(p10), floor);
    EXPECT_EQ("information_below_floor", a["patterns"][0]["error"]["code"].asString());
    EXPECT_FALSE(a["patterns"][0].isMember("max_mismatches_within_floor"));
    EXPECT_FALSE(a["patterns"][0].isMember("max_mismatches"));
    // no exemption: an exact pattern in suffix scope is gated at 1 (a 12-mer: 18.8 bits)
    a = run(*g, request("{\"dna\": \"" + c.P.substr(0, 12) + "\"}",
                        ", \"scope\": \"suffix\", \"max_mismatches\": 1"), floor);
    EXPECT_EQ("information_below_floor", a["patterns"][0]["error"]["code"].asString());
    EXPECT_EQ(0u, a["patterns"][0]["max_mismatches_within_floor"].asUInt());
    a = run(*g, request("{\"dna\": \"" + c.P.substr(0, 12) + "\"}", ", \"scope\": \"suffix\""),
            floor);
    EXPECT_FALSE(a["patterns"][0].isMember("error"));
}


// ---------------------------------------------------------------- the capabilities

TEST(PatternMismatchRoute, Capabilities) {
    std::mt19937 rng(27);
    const Planted c = planted(rng);
    auto g = build(kK, c.records);
    const std::vector<const AnnotatedDBG*> graphs = { g.get(), nullptr };
    for (const AnnotatedDBG *graph : graphs) {
        const Json::Value caps = pattern_capabilities_json(graph, limits(2, 1), false);
        ASSERT_TRUE(caps.isMember("mismatches"));
        const Json::Value &m = caps["mismatches"];
        std::set<std::string> fields;
        for (const std::string &f : m.getMemberNames()) {
            fields.insert(f);
        }
        EXPECT_EQ(std::set<std::string>({ "max", "default", "unit", "protein_substitutes",
                                          "edits", "order", "levels_rule", "later_increment" }),
                  fields);
        EXPECT_EQ(2u, m["max"]["dna"].asUInt64());
        EXPECT_EQ(2u, m["max"]["iupac"].asUInt64());
        EXPECT_EQ(1u, m["max"]["protein"].asUInt64());
        EXPECT_EQ(3u, m["max"].size());
        EXPECT_EQ(0, m["default"].asInt());
        EXPECT_EQ("base", m["unit"]["dna"].asString());
        EXPECT_EQ("base", m["unit"]["iupac"].asString());
        EXPECT_EQ("residue", m["unit"]["protein"].asString());
        EXPECT_EQ("sense_codons", m["protein_substitutes"].asString());
        ASSERT_EQ(1u, m["edits"].size());
        EXPECT_EQ("substitution", m["edits"][0].asString());
        EXPECT_EQ("mismatches_first", m["order"].asString());
        EXPECT_EQ("SPEC-pattern-search.md section 21.4", m["levels_rule"].asString());
        std::vector<std::string> later;
        for (const Json::Value &x : m["later_increment"]) {
            later.push_back(x.asString());
        }
        EXPECT_EQ(std::vector<std::string>({ "long_patterns", "labels", "predicate",
                                             "supported_paths" }), later);
        // not among the gate keys (SPEC §23: the full block only), not among the caps
        for (const std::string &gate : pattern_gate_keys()) {
            EXPECT_NE("mismatches", gate);
        }
        EXPECT_FALSE(caps["caps"].isMember("max_mismatches"));
    }
    // the defaults: the contract's maximum for every kind
    const Json::Value d = pattern_capabilities_json(g.get(), PatternLimits(), false)["mismatches"];
    EXPECT_EQ(3u, d["max"]["dna"].asUInt64());
    EXPECT_EQ(3u, d["max"]["protein"].asUInt64());
    // the multi-graph three-field block has none
    EXPECT_FALSE(pattern_capabilities_json(g.get(), limits(), true).isMember("mismatches"));
}

} // namespace
