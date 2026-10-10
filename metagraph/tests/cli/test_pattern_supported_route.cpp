#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// POST /pattern with long_search "supported_paths" (SPEC-pattern-search.md §20), with a
// predicate on the supported paths (§20.9) and with predicate_scope "motif" (§25), driven
// through THE ROUTE (process_pattern_request, src/cli/pattern.cpp) on tiny indexes built from
// explicit records in labelled columns, with their coordinates and record mapping. The
// expectations come from ORACLES over the records that never ask the engine, the trackers or
// the selection:
//  - a WALK oracle: every walk of the graph (the records' k-mers as deposited, consecutive
//    k-mers overlapping by k - 1) spelling the oriented pattern (this file's own IUPAC table),
//    and its support: the columns carrying every k-mer of it (label level), the columns one of
//    whose records holds it whole (record level); a walk is supported when its support is not
//    empty, and the search follows a branch only while its prefix is supported;
//  - a recursive evaluator of the request's predicate JSON (an unknown name is absent), on a
//    walk's support and, for predicate_strands "either", its reverse-complement walk's;
//  - the labels of the paths of long_search "paths" with output.labels "all" (the labelled
//    retrieval of every walk, its labels verified by a join of the k-mers' coordinate lists):
//    an independent implementation that the supported paths must agree with;
//  - for the motif: the union of the predicate's labels over the pattern's contexts, scanned
//    from the records' k-mers (and their reverse complements' rows for "either").

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;
namespace gp = mtg::graph::pattern;
using Clock = gp::Deadline::Clock;


// ------------------------------------------------------------------ the index

struct Record {
    std::string column;
    std::string header;
    std::string seq;
};

struct Index {
    size_t k = 0;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
};

// each record's k-mers numbered from its column's running count (as `annotate --coordinates`
// numbers them), the record mapping built from the same records; a row-diff annotation with
// coordinates (budgeted, placement record), or a column one without (unbudgeted, no
// coordinates)
template <class Annotation>
Index build(size_t k, const std::vector<Record> &records, bool coordinates) {
    Index idx;
    idx.k = k;
    idx.records = records;
    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    std::map<std::string, uint64_t> next;
    for (const Record &r : records) {
        EXPECT_GE(r.seq.size(), k);
        seqs.push_back(r.seq);
        labels.push_back(r.column);
        starts.push_back(next[r.column]);
        next[r.column] += r.seq.size() - k + 1;
    }
    idx.anno = test::build_anno_graph<DBGSuccinct, Annotation>(
            k, seqs, labels, DeBruijnGraph::BASIC, coordinates,
            coordinates ? starts : std::vector<uint64_t>{});
    if (coordinates) {
        const auto &encoder = idx.anno->get_annotator().get_label_encoder();
        std::vector<std::vector<std::string>> headers(encoder.size());
        std::vector<std::vector<uint64_t>> num_kmers(encoder.size());
        for (const Record &r : records) {
            const size_t c = encoder.encode(r.column);
            headers[c].push_back(r.header);
            num_kmers[c].push_back(r.seq.size() - k + 1);
        }
        idx.cth = std::make_unique<annot::CoordToHeader>(std::move(headers),
                                                         std::move(num_kmers));
    }
    return idx;
}

Index build_rd(size_t k, const std::vector<Record> &records) {
    return build<annot::RowDiffColumnAnnotator>(k, records, true);
}


// ------------------------------------------------------------------ the test's own tables

const std::map<char, std::string> kCodes = {
    { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" }, { 'R', "AG" }, { 'Y', "CT" },
    { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" }, { 'B', "CGT" },
    { 'D', "AGT" }, { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
};

char complement_code(char c) {
    const std::string set = kCodes.at(c);
    std::string out;
    for (char b : set) {
        out += b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
    }
    std::sort(out.begin(), out.end());
    for (const auto &[code, s] : kCodes) {
        std::string t = s;
        std::sort(t.begin(), t.end());
        if (t == out)
            return code;
    }
    return 'N';
}

std::string rc(const std::string &s) {
    std::string out(s.rbegin(), s.rend());
    for (char &c : out) {
        c = complement_code(c);
    }
    return out;
}

bool admits(char code, char base) {
    return kCodes.at(code).find(base) != std::string::npos;
}

bool matches(const std::string &pattern, const std::string &s) {
    for (size_t i = 0; i < s.size(); ++i) {
        if (!admits(pattern[i], s[i]))
            return false;
    }
    return true;
}

// the request predicate on a set of present column names (unknown names simply absent)
bool o_eval(const Json::Value &p, const std::set<std::string> &s) {
    const std::string op = p.getMemberNames().front();
    const Json::Value &v = p[op];
    if (op == "any" || op == "all" || op == "none") {
        size_t in = 0;
        for (const Json::Value &n : v) {
            in += s.count(n.asString());
        }
        return op == "any" ? in > 0 : op == "all" ? in == v.size() : in == 0;
    }
    if (op == "at_least") {
        size_t in = 0;
        for (const Json::Value &n : v["labels"]) {
            in += s.count(n.asString());
        }
        return in >= v["n"].asUInt64();
    }
    if (op == "not")
        return !o_eval(v, s);
    bool all = true, any = false;
    for (const Json::Value &q : v) {
        const bool x = o_eval(q, s);
        all &= x;
        any |= x;
    }
    return op == "and" ? all : any;
}

void o_names(const Json::Value &p, std::set<std::string> *out) {
    const std::string op = p.getMemberNames().front();
    const Json::Value &v = p[op];
    if (op == "any" || op == "all" || op == "none") {
        for (const Json::Value &n : v) {
            out->insert(n.asString());
        }
    } else if (op == "at_least") {
        for (const Json::Value &n : v["labels"]) {
            out->insert(n.asString());
        }
    } else if (op == "not") {
        o_names(v, out);
    } else {
        for (const Json::Value &q : v) {
            o_names(q, out);
        }
    }
}

// a predicate without none and not: false on a set stays false on its subsets
bool o_monotone(const Json::Value &p) {
    const std::string op = p.getMemberNames().front();
    if (op == "none" || op == "not")
        return false;
    if (op == "and" || op == "or") {
        for (const Json::Value &q : p[op]) {
            if (!o_monotone(q))
                return false;
        }
    }
    return true;
}


// ------------------------------------------------------------------ the walk oracle

struct Walk {
    std::string strand;      // "+", "-", "="
    std::string sequence;
};

struct Oracle {
    const Index &idx;
    const size_t k;
    std::set<std::string> kmers;
    std::map<std::string, std::set<std::string>> carried;
    std::map<std::string, std::vector<const Record*>> columns;

    explicit Oracle(const Index &idx) : idx(idx), k(idx.k) {
        for (const Record &r : idx.records) {
            columns[r.column].push_back(&r);
            for (size_t i = 0; i + k <= r.seq.size(); ++i) {
                kmers.insert(r.seq.substr(i, k));
                carried[r.column].insert(r.seq.substr(i, k));
            }
        }
    }

    // the columns carrying every k-mer of |s|
    std::set<std::string> carriers(const std::string &s) const {
        std::set<std::string> out;
        for (const auto &[column, set] : carried) {
            bool all = true;
            for (size_t i = 0; i + k <= s.size() && all; ++i) {
                all = set.count(s.substr(i, k));
            }
            if (all)
                out.insert(column);
        }
        return out;
    }

    // the columns one of whose records holds |s| whole
    std::set<std::string> holders(const std::string &s) const {
        std::set<std::string> out;
        for (const auto &[column, recs] : columns) {
            for (const Record *r : recs) {
                if (r->seq.find(s) != std::string::npos)
                    out.insert(column);
            }
        }
        return out;
    }

    std::set<std::string> support(const std::string &s, bool records) const {
        return records ? holders(s) : carriers(s);
    }

    // the oriented patterns searched: P (+) and rc(P) (-), or one (=) for a palindrome
    static std::vector<std::pair<std::string, std::string>> oriented(const std::string &p,
                                                                     const std::string &strands) {
        if (p == rc(p))
            return { { "=", p } };
        std::vector<std::pair<std::string, std::string>> out;
        if (strands != "reverse")
            out.emplace_back("+", p);
        if (strands != "forward")
            out.emplace_back("-", rc(p));
        return out;
    }

    /**
     * The walks of the graph spelling |pattern| in |strands|: every one (|records| ignored)
     * when |all|, else the supported ones at the level, followed only while their prefix is
     * supported. |walks|: the complete walks the supported search completes (a branch pruned at
     * its last k-mer included); |pruned_early|: a branch was pruned before L.
     */
    std::vector<Walk> walks(const std::string &pattern, const std::string &strands,
                            bool records, bool all, uint64_t *completed = nullptr,
                            bool *pruned_early = nullptr) const {
        std::vector<Walk> out;
        uint64_t done = 0;
        bool early = false;
        for (const auto &oriented_pattern : oriented(pattern, strands)) {
            const std::string &strand = oriented_pattern.first;
            const std::string &q = oriented_pattern.second;
            const size_t L = q.size();
            std::function<void(const std::string&)> dfs = [&](const std::string &s) {
                for (char b : std::string("ACGT")) {
                    const std::string t = s + b;
                    if (!admits(q[t.size() - 1], b) || !kmers.count(t.substr(t.size() - k)))
                        continue;
                    const bool alive = all || !support(t, records).empty();
                    if (t.size() == L) {
                        ++done;
                        if (alive)
                            out.push_back(Walk { strand, t });
                    } else if (alive) {
                        dfs(t);
                    } else {
                        early = true;
                    }
                }
            };
            for (const std::string &x : kmers) {
                if (!matches(q.substr(0, k), x))
                    continue;
                if (!all && support(x, records).empty()) {
                    early = true;
                    continue;
                }
                dfs(x);
            }
        }
        if (completed)
            *completed = done;
        if (pruned_early)
            *pruned_early = early;
        return out;
    }
};

std::set<std::pair<std::string, std::string>> walk_set(const std::vector<Walk> &w) {
    std::set<std::pair<std::string, std::string>> out;
    for (const Walk &x : w) {
        out.emplace(x.strand, x.sequence);
    }
    return out;
}


// ------------------------------------------------------------------ the route

struct Ask {
    std::vector<std::string> patterns;
    std::string kind = "dna";
    std::string mode = "all_or_count";
    std::string strands = "both";
    std::string long_search = "supported_paths";
    // "" (not named), "best", "label_intersection"
    std::string level;
    std::string labels = "none";
    bool occurrences = true;
    uint64_t max_paths = 1'000;
    bool stop_at_threshold = false;
    // "" (no predicate), else its JSON text
    std::string predicate;
    // "" (not named), "either", "context"
    std::string predicate_strands;
    // "" (not named), "context", "motif"
    std::string predicate_scope;
    uint64_t max_predicate_contexts = 100'000;
    uint64_t max_annotation_work = 0;
    std::string require_support;
    // the account in bytes (0: the request's), a virtual clock, the record mapping
    uint64_t max_memory_bytes = 0;
    std::function<Clock::time_point()> clock;
    bool records = true;
};

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

std::string body_of(const Ask &r) {
    Json::Value b;
    for (const std::string &p : r.patterns) {
        Json::Value q;
        q[r.kind] = p;
        b["patterns"].append(q);
    }
    b["mode"] = r.mode;
    b["strands"] = r.strands;
    b["long_search"] = r.long_search;
    if (!r.level.empty())
        b["supported_paths_level"] = r.level;
    b["output"]["labels"] = r.labels;
    if (r.labels != "none")
        b["output"]["occurrences"] = r.occurrences;
    b["max_paths"] = static_cast<Json::UInt64>(r.max_paths);
    b["stop_at_threshold"] = r.stop_at_threshold;
    if (!r.predicate.empty()) {
        b["predicate"] = parse_pattern_body(r.predicate);
        b["max_predicate_contexts"] = static_cast<Json::UInt64>(r.max_predicate_contexts);
    }
    if (!r.predicate_strands.empty())
        b["predicate_strands"] = r.predicate_strands;
    if (!r.predicate_scope.empty())
        b["predicate_scope"] = r.predicate_scope;
    if (r.max_annotation_work)
        b["max_annotation_work"] = static_cast<Json::UInt64>(r.max_annotation_work);
    if (!r.require_support.empty())
        b["require_support"] = r.require_support;
    return compact(b);
}

PatternLimits route_limits() {
    PatternLimits l;
    l.min_information_bits = 0;
    return l;
}

Json::Value answer_of(const Index &idx, const Ask &r) {
    RetrievalHooks hooks;
    if (r.records)
        hooks.coord_to_header = idx.cth.get();
    hooks.max_memory_bytes = r.max_memory_bytes;
    return process_pattern_request(parse_pattern_body(body_of(r)), *idx.anno, route_limits(),
                                   "rel", nullptr, nullptr, r.clock, &hooks);
}

std::set<std::pair<std::string, std::string>> result_set(const Json::Value &e) {
    std::set<std::pair<std::string, std::string>> out;
    for (const Json::Value &r : e["results"]) {
        out.emplace(r["strand"].asString(), r["sequence"].asString());
    }
    return out;
}

// the results are in answer order: anchor node, orientation (+, -, =), sequence
void check_answer_order(const Json::Value &e) {
    auto rank = [](const std::string &s) { return s == "+" ? 0 : s == "-" ? 1 : 2; };
    for (Json::ArrayIndex i = 1; i < e["results"].size(); ++i) {
        const Json::Value &a = e["results"][i - 1], &b = e["results"][i];
        EXPECT_LE(std::make_tuple(a["nodes"][0].asUInt64(), rank(a["strand"].asString()),
                                  a["sequence"].asString()),
                  std::make_tuple(b["nodes"][0].asUInt64(), rank(b["strand"].asString()),
                                  b["sequence"].asString()));
    }
}

// the supported paths' fields every answered entry of a pattern longer than k carries
void check_supported_shape(const Json::Value &e, const std::string &level) {
    const Json::Value &sp = e["counts"]["supported_paths"];
    ASSERT_TRUE(sp.isObject()) << compact(e);
    EXPECT_EQ("paths", sp["unit"].asString());
    EXPECT_EQ(level, sp["level"].asString());
    EXPECT_TRUE(sp.isMember("by_strand"));
    EXPECT_TRUE(sp.isMember("search"));
    EXPECT_TRUE(sp.isMember("candidates_examined"));
    EXPECT_TRUE(sp.isMember("branches_pruned"));
    // counts.paths is a plain count, not long_search "paths"' paths_count
    const Json::Value &paths = e["counts"]["paths"];
    EXPECT_EQ(3u, paths.size()) << compact(paths);
    for (const char *f : { "annotation_rows", "annotation_units", "memory_bytes", "anchor_rows",
                           "row_cache_hits", "row_cache_evictions", "extension_edges" }) {
        EXPECT_TRUE(e["work"].isMember(f)) << f;
    }
    EXPECT_TRUE(e["timing"].isMember("support_ms"));
    EXPECT_TRUE(e["rows_refused"].isArray());
}

// shared prefixes and a planted repeat: walks through the repeat are mosaics of the records
std::vector<Record> random_records(std::mt19937 &rng, size_t k, size_t columns) {
    auto bases = [&](size_t n) {
        std::string s;
        for (size_t i = 0; i < n; ++i) {
            s += "ACGT"[rng() % 4];
        }
        return s;
    };
    const std::string repeat = bases(k - 1);
    const std::string common = bases(k + 2);
    std::vector<Record> records;
    for (size_t c = 0; c < columns; ++c) {
        const size_t n = 1 + rng() % 2;
        for (size_t j = 0; j < n; ++j) {
            std::string s = bases(2 + rng() % 4);
            if (rng() % 2)
                s += common;
            s += bases(1 + rng() % 3) + repeat + bases(2 + rng() % 4);
            if (rng() % 3 == 0)
                s += repeat + bases(2 + rng() % 3);
            if (rng() % 4 == 0)
                s = rc(s);
            records.push_back(Record { "c" + std::to_string(c),
                                       "c" + std::to_string(c) + "_r" + std::to_string(j), s });
        }
    }
    // a column holding another's record on its other strand (a reverse-complement walk then
    // has a support of its own: predicate_strands "either" reads it)
    const size_t n = records.size();
    for (size_t i = 0; i < n; ++i) {
        if (rng() % 3)
            continue;
        const std::string column = "c" + std::to_string(rng() % columns);
        records.push_back(Record { column, column + "_rc" + std::to_string(i),
                                   rc(records[i].seq) });
    }
    return records;
}

// a pattern longer than k from a record (or its reverse complement), a position or two made
// degenerate so that it branches
std::string random_pattern(std::mt19937 &rng, const std::vector<Record> &records, size_t k) {
    for (;;) {
        const Record &r = records[rng() % records.size()];
        std::string s = rng() % 3 ? r.seq : rc(r.seq);
        const size_t max = std::min<size_t>(s.size(), 2 * k + 3);
        if (max <= k + 1)
            continue;
        const size_t L = k + 1 + rng() % (max - k);
        const size_t start = rng() % (s.size() - L + 1);
        std::string p = s.substr(start, L);
        const size_t degenerate = rng() % 3;
        for (size_t i = 0; i < degenerate; ++i) {
            const size_t at = rng() % L;
            const char code = "RYSWKMN"[rng() % 7];
            if (admits(code, p[at]))
                p[at] = code;
        }
        return p;
    }
}

// predicates over the columns c0.. and a name no index has
std::string random_predicate(std::mt19937 &rng, size_t columns) {
    auto name = [&]() {
        return "\"" + (rng() % 8 ? "c" + std::to_string(rng() % columns) : std::string("zz"))
                + "\"";
    };
    auto leaf = [&]() {
        const char *ops[] = { "any", "all", "none" };
        std::string a = name(), b = name();
        std::string list = a == b ? a : a + ", " + b;
        return std::string("{\"") + ops[rng() % 3] + "\": [" + list + "]}";
    };
    switch (rng() % 5) {
        case 0: return leaf();
        case 1: return "{\"and\": [" + leaf() + ", " + leaf() + "]}";
        case 2: return "{\"or\": [" + leaf() + ", " + leaf() + "]}";
        case 3: return "{\"not\": " + leaf() + "}";
        default: return "{\"any\": [" + name() + "]}";
    }
}


// ------------------------------------------------------------------ supported paths

// The design's k = 3 example: column A holds the records ACG and CGT; ACGT is a walk of the
// graph that no record holds whole. Label level: supported by A; record level: not supported
// (pruned at its last k-mer, a complete walk all the same)
TEST(PatternSupportedRoute, TwoRecordsOneLabel) {
    const Index idx = build_rd(3, { { "A", "a0", "ACG" }, { "A", "a1", "CGT" } });
    Ask r;
    r.patterns = { "ACGT" };
    r.strands = "forward";
    r.labels = "all";
    for (const std::string level : { "best", "label_intersection" }) {
        r.level = level;
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        const bool records = level == "best";
        check_supported_shape(e, records ? "record_verified" : "label_intersection");
        EXPECT_EQ(1u, e["counts"]["paths"]["value"].asUInt64());
        EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
        const Json::Value &sp = e["counts"]["supported_paths"];
        EXPECT_EQ("exact", sp["relation"].asString());
        EXPECT_EQ(records ? 0u : 1u, sp["value"].asUInt64());
        EXPECT_EQ(records ? 1u : 0u, sp["branches_pruned"].asUInt64());
        EXPECT_TRUE(e["retrieval_complete"].asBool());
        ASSERT_EQ(records ? 0u : 1u, e["results"].size());
        // what the answer places: the label level reads no coordinate, so nothing is placed
        // on this index that could place (and the labels are listed without occurrences)
        EXPECT_EQ(records ? "record" : "none", e["placement"].asString());
        EXPECT_EQ(records ? "exact" : "unknown", e["counts"]["occurrences"]["relation"].asString());
        bool only = false;
        for (const Json::Value &n : e["notes"]) {
            only |= n.asString() == "label_intersection_only";
        }
        EXPECT_EQ(!records, only);
        if (!records) {
            EXPECT_EQ("ACGT", e["results"][0]["sequence"].asString());
            EXPECT_EQ("A", e["results"][0]["labels"][0]["column"].asString());
            EXPECT_EQ("label_intersection", e["results"][0]["labels"][0]["support"].asString());
            EXPECT_FALSE(e["results"][0]["labels"][0].isMember("occurrences"));
            EXPECT_FALSE(e["results"][0]["labels"][0].isMember("occurrence_list"));
        }
        // a predicate: none(A) rejects the walk at the label level and has no candidate at
        // the record level
        Ask q = r;
        q.predicate = "{\"none\": [\"A\"]}";
        const Json::Value f = answer_of(idx, q)["patterns"][0];
        EXPECT_EQ(records ? 0u : 1u, f["counts"]["tested"]["value"].asUInt64());
        EXPECT_EQ(0u, f["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ("exact", f["counts"]["selected"]["relation"].asString());
        EXPECT_EQ("completed", f["selection"]["pass"].asString());
        EXPECT_EQ(0u, f["results"].size());
        EXPECT_TRUE(f["retrieval_complete"].asBool());
    }
}

// Random indexes with mosaics, random patterns: the supported paths are the oracle's at both
// levels, in answer order, each with the labels the oracle names; the walks' count is exact
// exactly when no branch was pruned before L; and the supported paths are the paths of
// long_search "paths" that list a label of the level (an independent implementation)
TEST(PatternSupportedRoute, RandomIndexesAgainstTheWalkOracle) {
    std::mt19937 rng(20261009);
    size_t checked = 0, supported_total = 0, pruned_cases = 0;
    for (size_t round = 0; round < 40; ++round) {
        const size_t k = 5 + rng() % 3;
        const std::vector<Record> records = random_records(rng, k, 2 + rng() % 3);
        const Index idx = build_rd(k, records);
        const Oracle oracle(idx);
        for (size_t t = 0; t < 4; ++t) {
            const std::string p = random_pattern(rng, records, k);
            for (const std::string level : { "best", "label_intersection" }) {
                const bool rec = level == "best";
                for (const std::string strands : { "both", "forward" }) {
                    Ask r;
                    r.patterns = { p };
                    r.kind = "iupac";
                    r.level = level;
                    r.strands = strands;
                    r.labels = "all";
                    const Json::Value e = answer_of(idx, r)["patterns"][0];
                    ASSERT_FALSE(e.isMember("error")) << compact(e);
                    check_supported_shape(e, rec ? "record_verified" : "label_intersection");
                    uint64_t walks = 0;
                    bool early = false;
                    const std::vector<Walk> expected
                            = oracle.walks(p, strands, rec, false, &walks, &early);
                    const Json::Value &sp = e["counts"]["supported_paths"];
                    EXPECT_EQ("exact", sp["relation"].asString()) << p;
                    EXPECT_EQ(expected.size(), sp["value"].asUInt64()) << p << " " << level;
                    EXPECT_EQ(walk_set(expected), result_set(e)) << p << " " << level;
                    EXPECT_TRUE(e["retrieval_complete"].asBool());
                    check_answer_order(e);
                    // counts.paths: the complete walks, exact iff nothing was pruned before L
                    const Json::Value &paths = e["counts"]["paths"];
                    EXPECT_EQ(walks, paths["value"].asUInt64()) << p;
                    EXPECT_EQ(early ? "at_least" : "exact", paths["relation"].asString()) << p;
                    if (!early) {
                        EXPECT_EQ(oracle.walks(p, strands, rec, true).size(),
                                  paths["value"].asUInt64());
                    }
                    pruned_cases += early;
                    // each path's labels: the columns carrying it, verified where a record
                    // holds it whole (record level)
                    for (const Json::Value &res : e["results"]) {
                        const std::string s = res["sequence"].asString();
                        std::set<std::string> listed, verified;
                        for (const Json::Value &l : res["labels"]) {
                            listed.insert(l["column"].asString());
                            if (l["support"].asString() == "record_verified")
                                verified.insert(l["column"].asString());
                        }
                        EXPECT_EQ(oracle.carriers(s), listed) << s;
                        EXPECT_EQ(rec ? oracle.holders(s) : std::set<std::string>(), verified)
                                << s;
                    }
                    // differential: long_search "paths" with labels "all", the paths listing a
                    // label of the level
                    Ask d = r;
                    d.long_search = "paths";
                    d.level.clear();
                    const Json::Value g = answer_of(idx, d)["patterns"][0];
                    std::set<std::pair<std::string, std::string>> from_paths;
                    for (const Json::Value &res : g["results"]) {
                        bool any = false;
                        for (const Json::Value &l : res["labels"]) {
                            any |= !rec || l["support"].asString() == "record_verified";
                        }
                        if (any)
                            from_paths.emplace(res["strand"].asString(),
                                               res["sequence"].asString());
                    }
                    EXPECT_EQ(from_paths, result_set(e)) << p << " " << level;
                    ++checked;
                    supported_total += expected.size();
                }
            }
        }
    }
    // the cases are not all trivial
    EXPECT_GT(supported_total, checked / 2);
    EXPECT_GT(pruned_cases, 10u);
}

// Short patterns are answered as without the option; the request reads the annotation in every
// mode (limits echo it), and a count states the supported paths without listing any
TEST(PatternSupportedRoute, ShortPatternsAlikeAndCountMode) {
    std::mt19937 rng(7);
    const size_t k = 6;
    const std::vector<Record> records = random_records(rng, k, 3);
    const Index idx = build_rd(k, records);
    const std::string shortp = records[0].seq.substr(0, 4);
    const std::string longp = records[0].seq.substr(0, k + 3);
    for (const std::string mode : { "count", "all_or_count", "partial" }) {
        for (const std::string labels : { "none", "all" }) {
            Ask r;
            r.patterns = { shortp, longp };
            r.mode = mode;
            r.labels = labels;
            const Json::Value a = answer_of(idx, r);
            Ask d = r;
            d.long_search = "paths";
            const Json::Value b = answer_of(idx, d);
            Json::Value x = a["patterns"][0], y = b["patterns"][0];
            x.removeMember("timing");
            y.removeMember("timing");
            EXPECT_EQ(compact(y), compact(x)) << mode << " " << labels;
            EXPECT_EQ("supported_paths", a["limits"]["long_search"].asString());
            EXPECT_EQ("best", a["limits"]["supported_paths_level"].asString());
            EXPECT_TRUE(a["limits"].isMember("max_annotation_work"));
            const Json::Value &e = a["patterns"][1];
            check_supported_shape(e, "record_verified");
            const Oracle oracle(idx);
            EXPECT_EQ(oracle.walks(longp, "both", true, false).size(),
                      e["counts"]["supported_paths"]["value"].asUInt64());
            if (mode == "count") {
                EXPECT_FALSE(e.isMember("results"));
                EXPECT_FALSE(e["retrieval_complete"].asBool());
                // named and not built: no projection (the search read rows all the same)
                const bool named = labels == "all";
                bool noted = false;
                for (const Json::Value &n : e["notes"]) {
                    noted |= n.asString() == "projection_not_read";
                    EXPECT_NE("annotation_not_read", n.asString());
                }
                EXPECT_EQ(named, noted);
            }
        }
    }
}

// projection_not_read under supported paths names a projection asked for and not built: the
// fields only the labels of the listed paths read (output.labels, max_labels,
// max_occurrences_per_label, require_support); max_annotation_work and max_memory_mb bound the
// search itself, which read its rows, so naming them says nothing; with a predicate alike
TEST(PatternSupportedRoute, ProjectionNotReadNamesTheProjectionsFields) {
    std::mt19937 rng(11);
    const size_t k = 6;
    const std::vector<Record> records = random_records(rng, k, 3);
    const Index idx = build_rd(k, records);
    const std::string longp = records[0].seq.substr(0, k + 3);
    auto noted = [&](const std::string &mode, const std::string &predicate,
                     const std::string &key, const Json::Value &value) {
        Ask r;
        r.patterns = { longp };
        r.mode = mode;
        r.predicate = predicate;
        Json::Value body = parse_pattern_body(body_of(r));
        if (!key.empty())
            body[key] = value;
        RetrievalHooks hooks;
        hooks.coord_to_header = idx.cth.get();
        const Json::Value a = process_pattern_request(body, *idx.anno, route_limits(), "rel",
                                                      nullptr, nullptr, nullptr, &hooks);
        const Json::Value &e = a["patterns"][0];
        EXPECT_TRUE(e.isMember("counts")) << compact(a);
        bool found = false;
        for (const Json::Value &n : e["notes"]) {
            found |= n.asString() == "projection_not_read";
            EXPECT_NE("annotation_not_read", n.asString());
        }
        return found;
    };
    for (const std::string mode : { "count", "all_or_count", "partial" }) {
        for (const std::string &predicate : { std::string(),
                                              std::string("{\"any\": [\"c0\"]}") }) {
            const std::string what = mode + " " + predicate;
            EXPECT_FALSE(noted(mode, predicate, "", Json::Value())) << what;
            EXPECT_FALSE(noted(mode, predicate, "max_annotation_work", 1'000'000)) << what;
            EXPECT_FALSE(noted(mode, predicate, "max_memory_mb", 64)) << what;
            EXPECT_TRUE(noted(mode, predicate, "max_labels", 5)) << what;
            EXPECT_TRUE(noted(mode, predicate, "max_occurrences_per_label", 5)) << what;
            EXPECT_TRUE(noted(mode, predicate, "require_support", "label_intersection")) << what;
        }
    }
}

// max_paths on the supported paths: all_or_count withholds above it, partial lists the first
// ones, stop_at_threshold stops the search once more are complete
TEST(PatternSupportedRoute, ThresholdsOnTheSupportedPaths) {
    // four records holding one 9-mer each with a different base at one position: an IUPAC
    // pattern with an N there has 4 supported walks on one strand
    const std::vector<Record> records = {
        { "a", "a0", "TTACGTAACGGAT" }, { "b", "b0", "TTACGTCACGGAT" },
        { "c", "c0", "TTACGTGACGGAT" }, { "d", "d0", "TTACGTTACGGAT" },
    };
    const Index idx = build_rd(5, records);
    Ask r;
    r.patterns = { "ACGTNACGG" };
    r.kind = "iupac";
    r.strands = "forward";
    r.max_paths = 2;
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ(4u, e["counts"]["supported_paths"]["value"].asUInt64());
        EXPECT_EQ("count_above_threshold", e["withheld"]["reason"].asString());
        EXPECT_EQ(0u, e["results"].size());
    }
    r.mode = "partial";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ(2u, e["results"].size());
        EXPECT_EQ("max_paths", e["cut"]["reason"].asString());
        EXPECT_FALSE(e["retrieval_complete"].asBool());
        check_answer_order(e);
    }
    r.stop_at_threshold = true;
    for (const std::string mode : { "all_or_count", "partial" }) {
        r.mode = mode;
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("extension", e["stop"]["phase"].asString());
        EXPECT_EQ("max_paths", e["stop"]["reason"].asString());
        EXPECT_EQ("at_least", e["counts"]["supported_paths"]["relation"].asString());
        EXPECT_EQ(3u, e["counts"]["supported_paths"]["value"].asUInt64());
        if (mode == "all_or_count") {
            EXPECT_EQ("threshold_crossed", e["withheld"]["reason"].asString());
        } else {
            EXPECT_EQ(2u, e["results"].size());
            EXPECT_EQ("max_paths", e["cut"]["reason"].asString());
        }
    }
    // with a predicate the threshold is on the selected paths: any(a, b) selects 2 of the 4
    r.stop_at_threshold = false;
    r.mode = "all_or_count";
    r.predicate = "{\"any\": [\"a\", \"b\"]}";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ(2u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_TRUE(e["withheld"].isNull()) << compact(e);
        EXPECT_EQ(2u, e["results"].size());
        EXPECT_TRUE(e["retrieval_complete"].asBool());
    }
    r.predicate = "{\"any\": [\"a\", \"b\", \"c\"]}";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ(3u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ("selected_above_threshold", e["withheld"]["reason"].asString());
    }
    r.stop_at_threshold = true;
    r.predicate_strands = "context";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("extension", e["stop"]["phase"].asString());
        EXPECT_EQ("max_paths", e["stop"]["reason"].asString());
        EXPECT_EQ("threshold_crossed", e["withheld"]["reason"].asString());
        EXPECT_EQ("stopped", e["selection"]["pass"].asString());
        EXPECT_EQ(3u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ("at_least", e["counts"]["selected"]["relation"].asString());
    }
}

// The note low_complexity_pattern beside a threshold the predicate's selection raises
// (stop_at_threshold: max_paths on the selected paths, max_predicate_contexts on the held
// walks): diagnosed, as beside the engine's own thresholds; beside the search's own budget
// (max_annotation_work): not stated
TEST(PatternSupportedRoute, LowComplexityNoteBesideTheSelectionsThresholds) {
    std::string repeat;
    for (int i = 0; i < 14; ++i) {
        repeat += "GCC";
    }
    const std::vector<Record> records = {
        { "a", "a0", "TTAT" + repeat + "TATT" }, { "b", "b0", "ATTA" + repeat + "ATAT" },
    };
    const Index idx = build_rd(15, records);
    Ask r;
    // flagged by sdust with the seeder's parameters (PatternSupport's oracle)
    r.patterns = { "GCCGCCGCCGCCGCCGCCGC" };
    r.strands = "forward";
    r.predicate = "{\"any\": [\"a\"]}";
    r.stop_at_threshold = true;
    auto noted = [](const Json::Value &e) {
        for (const Json::Value &n : e["notes"]) {
            if (n.asString() == "low_complexity_pattern")
                return true;
        }
        return false;
    };
    {
        Ask full = r;
        full.stop_at_threshold = false;
        const Json::Value e = answer_of(idx, full)["patterns"][0];
        ASSERT_TRUE(e["stop"].isNull()) << compact(e);
        ASSERT_LE(1u, e["counts"]["selected"]["value"].asUInt64()) << compact(e);
        EXPECT_TRUE(noted(e)) << compact(e);
    }
    {
        Ask t = r;
        t.predicate_strands = "context";
        t.max_paths = 0;
        const Json::Value e = answer_of(idx, t)["patterns"][0];
        EXPECT_EQ("extension", e["stop"]["phase"].asString()) << compact(e);
        EXPECT_EQ("max_paths", e["stop"]["reason"].asString());
        EXPECT_TRUE(noted(e)) << compact(e);
    }
    {
        Ask t = r;
        t.max_predicate_contexts = 0;
        const Json::Value e = answer_of(idx, t)["patterns"][0];
        EXPECT_EQ("extension", e["stop"]["phase"].asString()) << compact(e);
        EXPECT_EQ("max_predicate_contexts", e["stop"]["reason"].asString());
        EXPECT_TRUE(noted(e)) << compact(e);
    }
    {
        Ask b = r;
        b.stop_at_threshold = false;
        b.max_annotation_work = 1;
        const Json::Value e = answer_of(idx, b)["patterns"][0];
        EXPECT_EQ("extension", e["stop"]["phase"].asString()) << compact(e);
        EXPECT_EQ("max_annotation_work", e["stop"]["reason"].asString());
        EXPECT_FALSE(noted(e)) << compact(e);
    }
}


// ------------------------------------------------------------------ predicates (§20.9)

// The selection on random indexes against the oracle: a supported walk is selected when the
// predicate holds on its support ("context") or on its support and its reverse-complement
// walk's ("either"), whatever the strands searched (one strand reads the mirror walks); the
// selected paths are listed with their selection labels and strands; the counts are exact;
// a monotone predicate under "context" prunes, its supported paths then at_least
TEST(PatternSupportedRoute, PredicatesAgainstTheOracle) {
    std::mt19937 rng(91);
    size_t checked = 0, selected_total = 0, mirror_reads = 0, pruned = 0;
    for (size_t round = 0; round < 30; ++round) {
        const size_t k = 5 + rng() % 2;
        const size_t columns = 2 + rng() % 3;
        const std::vector<Record> records = random_records(rng, k, columns);
        const Index idx = build_rd(k, records);
        const Oracle oracle(idx);
        for (size_t t = 0; t < 3; ++t) {
            const std::string p = random_pattern(rng, records, k);
            const std::string predicate = random_predicate(rng, columns);
            const Json::Value pj = parse_pattern_body(predicate);
            std::set<std::string> names;
            o_names(pj, &names);
            for (const std::string level : { "best", "label_intersection" }) {
                const bool rec = level == "best";
                for (const std::string ps : { "either", "context" }) {
                    for (const std::string strands : { "both", "forward", "reverse" }) {
                        Ask r;
                        r.patterns = { p };
                        r.kind = "iupac";
                        r.level = level;
                        r.strands = strands;
                        r.predicate = predicate;
                        r.predicate_strands = ps;
                        r.labels = "predicate_only";
                        const Json::Value a = answer_of(idx, r);
                        const Json::Value &e = a["patterns"][0];
                        ASSERT_FALSE(e.isMember("error")) << compact(e);
                        if (e["selection"]["pass"].asString() == "constant")
                            continue;
                        // the oracle's selection
                        std::set<std::pair<std::string, std::string>> expected;
                        std::map<std::pair<std::string, std::string>,
                                 std::map<std::string, std::string>> strands_of;
                        const std::vector<Walk> supported = oracle.walks(p, strands, rec, false);
                        for (const Walk &w : supported) {
                            std::set<std::string> own = oracle.support(w.sequence, rec);
                            std::set<std::string> mirror;
                            if (ps == "either")
                                mirror = oracle.support(rc(w.sequence), rec);
                            std::set<std::string> present = own;
                            present.insert(mirror.begin(), mirror.end());
                            if (!o_eval(pj, present))
                                continue;
                            expected.emplace(w.strand, w.sequence);
                            auto &m = strands_of[{ w.strand, w.sequence }];
                            for (const std::string &n : present) {
                                if (!names.count(n))
                                    continue;
                                m[n] = own.count(n) && mirror.count(n) ? "both"
                                     : own.count(n) ? "context" : "reverse_complement";
                            }
                        }
                        const std::string what = p + " " + predicate + " " + level + " " + ps
                                                    + " " + strands;
                        EXPECT_EQ("completed", e["selection"]["pass"].asString()) << what;
                        EXPECT_EQ(rec ? "record_verified" : "label_intersection",
                                  e["selection"]["support"].asString());
                        EXPECT_EQ("predicate", e["absence_filter"].asString());
                        EXPECT_EQ(expected.size(), e["counts"]["selected"]["value"].asUInt64())
                                << what;
                        EXPECT_EQ("exact", e["counts"]["selected"]["relation"].asString());
                        EXPECT_EQ(expected, result_set(e)) << what;
                        EXPECT_TRUE(e["retrieval_complete"].asBool()) << what;
                        check_answer_order(e);
                        const uint64_t by_predicate
                                = e["counts"]["supported_paths"]["branches_pruned_by_predicate"]
                                        .asUInt64();
                        if (ps == "either" || !o_monotone(pj)) {
                            // no pruning but the support's: every supported walk tested
                            EXPECT_EQ(0u, by_predicate) << what;
                            EXPECT_EQ(supported.size(),
                                      e["counts"]["supported_paths"]["value"].asUInt64()) << what;
                            EXPECT_EQ(supported.size(), e["counts"]["tested"]["value"].asUInt64())
                                    << what;
                        } else {
                            EXPECT_LE(e["counts"]["tested"]["value"].asUInt64(), supported.size());
                            EXPECT_EQ(by_predicate ? "at_least" : "exact",
                                      e["counts"]["supported_paths"]["relation"].asString())
                                    << what;
                            pruned += by_predicate > 0;
                        }
                        // the mirror walks are read when the search did not reach them: rows
                        // of the annotation work, none of the selection's own
                        const bool reads = ps == "either" && strands != "both" && p != rc(p);
                        const uint64_t mirror_rows = e["work"]["mirror_rows"].asUInt64();
                        if (!reads) {
                            EXPECT_EQ(0u, mirror_rows) << what;
                        }
                        EXPECT_EQ(0u, e["work"]["predicate_rows"].asUInt64()) << what;
                        EXPECT_LE(mirror_rows, e["work"]["annotation_rows"].asUInt64()) << what;
                        mirror_reads += mirror_rows > 0;
                        // each result: why it was selected, and on which walk
                        for (const Json::Value &res : e["results"]) {
                            const auto key = std::make_pair(res["strand"].asString(),
                                                            res["sequence"].asString());
                            std::set<std::string> why;
                            std::map<std::string, std::string> on;
                            for (Json::ArrayIndex i = 0; i < res["selection_labels"].size(); ++i) {
                                why.insert(res["selection_labels"][i].asString());
                                on[res["selection_labels"][i].asString()]
                                        = res["selection_strands"][i].asString();
                            }
                            EXPECT_TRUE(o_eval(pj, why)) << what;
                            EXPECT_EQ(strands_of[key], on) << what;
                            // predicate_only: the predicate's labels carrying the walk itself
                            std::set<std::string> listed;
                            for (const Json::Value &l : res["labels"]) {
                                listed.insert(l["column"].asString());
                            }
                            std::set<std::string> own_names;
                            for (const std::string &n : oracle.carriers(key.second)) {
                                if (names.count(n))
                                    own_names.insert(n);
                            }
                            EXPECT_EQ(own_names, listed) << what;
                        }
                        ++checked;
                        selected_total += expected.size();
                    }
                }
            }
        }
    }
    EXPECT_GT(selected_total, checked / 4);
    EXPECT_GT(mirror_reads, 5u);
    EXPECT_GT(pruned, 2u);
}

// At the record level a predicate sees the labels a record of which holds the walk whole, not
// the labels that merely carry its k-mers: B carries every k-mer of W in two records, neither
// holding W; A holds it. any(B) selects W at the label level only, and under "context" (a
// monotone predicate) the record level prunes W where B's chains end, at its fifth k-mer
TEST(PatternSupportedRoute, RecordLevelSupportIsARecordHoldingTheWalk) {
    const std::string w = "ACGTTGCAAT";
    const Index idx = build_rd(5, { { "A", "a0", w },
                                    { "B", "b0", w.substr(0, 8) }, { "B", "b1", w.substr(4) } });
    Ask r;
    r.patterns = { w };
    r.strands = "forward";
    r.labels = "all";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        ASSERT_EQ(1u, e["results"].size());
        EXPECT_EQ("mixed", e["results"][0]["support"].asString());
        std::map<std::string, std::string> support;
        for (const Json::Value &l : e["results"][0]["labels"]) {
            support[l["column"].asString()] = l["support"].asString();
        }
        EXPECT_EQ((std::map<std::string, std::string> { { "A", "record_verified" },
                                                        { "B", "label_intersection" } }),
                  support);
    }
    r.labels = "none";
    r.predicate = "{\"any\": [\"B\"]}";
    for (const std::string ps : { "either", "context" }) {
        r.predicate_strands = ps;
        for (const std::string level : { "best", "label_intersection" }) {
            r.level = level;
            const Json::Value e = answer_of(idx, r)["patterns"][0];
            const bool records = level == "best";
            EXPECT_EQ(records ? 0u : 1u, e["counts"]["selected"]["value"].asUInt64())
                    << ps << " " << level;
            EXPECT_EQ("exact", e["counts"]["selected"]["relation"].asString());
            EXPECT_EQ(records && ps == "context" ? 1u : 0u,
                      e["counts"]["supported_paths"]["branches_pruned_by_predicate"].asUInt64())
                    << ps << " " << level;
        }
    }
}

// A monotone predicate under "context" prunes the branches it is already false on: fewer rows
// read than without it, the same selection, the supported paths of the pruned orientation
// at_least; "either" prunes nothing (a walk's decision needs its mirror)
TEST(PatternSupportedRoute, MonotonePruningReadsFewerRows) {
    const std::vector<Record> records = {
        { "a", "a0", "TTACGTAACGGATTCCAGT" }, { "b", "b0", "TTACGTCACGGATTCCAGT" },
        { "c", "c0", "TTACGTGACGGATTCCAGT" },
    };
    const Index idx = build_rd(5, records);
    Ask r;
    r.patterns = { "ACGTNACGGATTCC" };
    r.kind = "iupac";
    r.strands = "forward";
    const Json::Value plain = answer_of(idx, r)["patterns"][0];
    r.predicate = "{\"any\": [\"a\"]}";
    r.predicate_strands = "context";
    const Json::Value pruned = answer_of(idx, r)["patterns"][0];
    EXPECT_EQ(1u, pruned["counts"]["selected"]["value"].asUInt64());
    EXPECT_EQ("exact", pruned["counts"]["selected"]["relation"].asString());
    EXPECT_EQ(2u, pruned["counts"]["supported_paths"]["branches_pruned_by_predicate"].asUInt64());
    EXPECT_EQ("at_least", pruned["counts"]["supported_paths"]["relation"].asString());
    EXPECT_EQ(1u, pruned["counts"]["supported_paths"]["value"].asUInt64());
    EXPECT_EQ("at_least", pruned["counts"]["paths"]["relation"].asString());
    EXPECT_LT(pruned["work"]["annotation_rows"].asUInt64(),
              plain["work"]["annotation_rows"].asUInt64());
    EXPECT_LT(0u, pruned["work"]["predicate_units"].asUInt64());
    r.predicate_strands = "either";
    const Json::Value either = answer_of(idx, r)["patterns"][0];
    EXPECT_EQ(0u, either["counts"]["supported_paths"]["branches_pruned_by_predicate"].asUInt64());
    EXPECT_EQ("exact", either["counts"]["supported_paths"]["relation"].asString());
    EXPECT_EQ(3u, either["counts"]["supported_paths"]["value"].asUInt64());
    EXPECT_EQ(1u, either["counts"]["selected"]["value"].asUInt64());
    // (the mirror walks are not in this graph: nothing read for them)
    EXPECT_EQ(0u, either["work"]["mirror_rows"].asUInt64());
    EXPECT_EQ(plain["work"]["annotation_rows"].asUInt64(),
              either["work"]["annotation_rows"].asUInt64());
}

// The mirror walks of "either" with one strand searched are read by a tracker of their own,
// after the search: their rows are annotation work (work.mirror_rows, part of
// work.annotation_rows; the selection's own predicate_rows stay 0), and their row cache
// replaces the search's in the account instead of adding to it. Column d holds the reverse
// complement of a's record: the forward walk's record-level support is {a}, its mirror's {d}
TEST(PatternSupportedRoute, MirrorWalksAreAnnotationWork) {
    const std::string a = "TTACGTAACGGATTCCAGT";
    const Index idx = build_rd(5, { { "a", "a0", a }, { "d", "d0", rc(a) } });
    Ask r;
    r.patterns = { "ACGTAACGGATTCC" };
    r.strands = "forward";
    const Json::Value plain = answer_of(idx, r)["patterns"][0];
    ASSERT_EQ(1u, plain["counts"]["supported_paths"]["value"].asUInt64());
    r.predicate = "{\"any\": [\"d\"]}";
    r.predicate_strands = "context";
    const Json::Value context = answer_of(idx, r)["patterns"][0];
    EXPECT_EQ(0u, context["counts"]["selected"]["value"].asUInt64());
    EXPECT_EQ(0u, context["work"]["mirror_rows"].asUInt64());
    r.predicate_strands = "either";
    r.labels = "predicate_only";
    const Json::Value either = answer_of(idx, r)["patterns"][0];
    EXPECT_EQ(1u, either["counts"]["selected"]["value"].asUInt64());
    EXPECT_EQ("exact", either["counts"]["selected"]["relation"].asString());
    ASSERT_EQ(1u, either["results"].size());
    EXPECT_EQ("[\"d\"]", compact(either["results"][0]["selection_labels"]));
    EXPECT_EQ("[\"reverse_complement\"]", compact(either["results"][0]["selection_strands"]));
    EXPECT_EQ(1u, either["work"]["predicate_lookups"].asUInt64());
    const uint64_t mirror_rows = either["work"]["mirror_rows"].asUInt64();
    // the mirror's 10 k-mers, its anchor's row read whole and with coordinates
    EXPECT_LE(10u, mirror_rows);
    EXPECT_GE(12u, mirror_rows);
    EXPECT_EQ(0u, either["work"]["predicate_rows"].asUInt64());
    EXPECT_EQ(plain["work"]["annotation_rows"].asUInt64() + mirror_rows,
              either["work"]["annotation_rows"].asUInt64());
    EXPECT_LT(plain["work"]["annotation_units"].asUInt64(),
              either["work"]["annotation_units"].asUInt64());
    // the row cache's allotment is most of the peak here: read twice at once, it would double it
    EXPECT_LT(either["work"]["memory_bytes"].asUInt64(),
              plain["work"]["memory_bytes"].asUInt64() * 5 / 4);
}

// A listed path's selection_labels are charged with their names' bytes (SPEC §20.9), not as ids
// alone: on the graph of all 64 3-mers (one record of every 3-mer in a row) NNNN has 256
// supported paths at the label level (the record holds only 123 of the 4-mers whole), each
// selected by any(<the one label>) and listing its name; a name of 64 KiB under an account of
// 1 MiB (max_labels 0: the paths' own label lists copy nothing) cannot be listed 256 times
// (16 MiB), so partial ends the list at the account, stop {output, max_memory}, with an answer
// below the account, and all_or_count withholds; a short name lists every path
TEST(PatternSupportedRoute, SelectionLabelNamesAreCharged) {
    std::string record;
    for (char a : std::string("ACGT")) {
        for (char b : std::string("ACGT")) {
            for (char c : std::string("ACGT")) {
                record += a;
                record += b;
                record += c;
            }
        }
    }
    const std::string long_name(64 * 1024, 'x');
    auto body = [&](const std::string &name, const std::string &mode) {
        return "{\"patterns\": [{\"iupac\": \"NNNN\"}], \"mode\": \"" + mode + "\", \"strands\": "
               "\"forward\", \"long_search\": \"supported_paths\", \"supported_paths_level\": "
               "\"label_intersection\", \"predicate\": {\"any\": [\"" + name + "\"]}, "
               "\"predicate_strands\": \"context\", \"output\": {\"labels\": \"all\", "
               "\"occurrences\": false}, \"max_labels\": 0}";
    };
    auto answer = [&](const Index &idx, const std::string &name, const std::string &mode) {
        RetrievalHooks hooks;
        hooks.coord_to_header = idx.cth.get();
        hooks.max_memory_bytes = 1 << 20;
        return process_pattern_request(parse_pattern_body(body(name, mode)), *idx.anno,
                                       route_limits(), "rel", nullptr, nullptr, nullptr, &hooks);
    };
    {
        const Index idx = build_rd(3, { { long_name, "h0", record } });
        Json::Value out = answer(idx, long_name, "partial");
        const Json::Value &e = out["patterns"][0];
        ASSERT_FALSE(e.isMember("error")) << compact(e["error"]);
        EXPECT_EQ("exact", e["counts"]["selected"]["relation"].asString());
        EXPECT_EQ(256u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ("output", e["stop"]["phase"].asString());
        EXPECT_EQ("max_memory", e["stop"]["reason"].asString());
        EXPECT_EQ("max_memory", e["cut"]["reason"].asString());
        const size_t listed = e["results"].size();
        EXPECT_LT(0u, listed);
        EXPECT_GT(256u, listed);
        EXPECT_EQ(listed, e["returned"].asUInt64());
        for (const Json::Value &r : e["results"]) {
            ASSERT_EQ(1u, r["selection_labels"].size());
            EXPECT_EQ(long_name, r["selection_labels"][0].asString());
            EXPECT_EQ("context", r["selection_strands"][0].asString());
        }
        // the account holds the listed copies of the name, and the answer is within it
        EXPECT_LE(listed * long_name.size(), e["work"]["memory_bytes"].asUInt64());
        EXPECT_GE(1u << 20, e["work"]["memory_bytes"].asUInt64());
        EXPECT_GT(1u << 20, compact(out).size());

        out = answer(idx, long_name, "all_or_count");
        const Json::Value &w = out["patterns"][0];
        EXPECT_EQ("output_budget", w["withheld"]["reason"].asString());
        EXPECT_EQ("max_memory", w["stop"]["reason"].asString());
        EXPECT_EQ(0u, w["results"].size());
        EXPECT_EQ(256u, w["counts"]["selected"]["value"].asUInt64());
        EXPECT_GT(256u * 1024, compact(out).size());
    }
    {
        const Index idx = build_rd(3, { { "a", "h0", record } });
        for (const char *mode : { "partial", "all_or_count" }) {
            const Json::Value e = answer(idx, "a", mode)["patterns"][0];
            EXPECT_TRUE(e["stop"].isNull()) << mode;
            EXPECT_TRUE(e["cut"].isNull()) << mode;
            EXPECT_TRUE(e["withheld"].isNull()) << mode;
            EXPECT_EQ(256u, e["results"].size()) << mode;
        }
    }
}

// The mirror of a held walk ("either" on a BASIC graph) is read only when it decides (SPEC
// §20.9). On the graph of all 64 3-mers (one record of every 3-mer in a row) with one label
// shared on every k-mer, ANNNNN has 1,024 forward walks at the label level, each selected by
// any(shared) on its own support alone (a definite true, and count mode asks no label
// evidence), so "either" reads no mirror: the same 1,024 exact as "context", no lookup, no
// mirror row, the same annotation_rows -- 64, the 64 k-mers' rows each read once through the
// pattern's row cache (the anchors' whole rows are among them) -- and under max_annotation_work
// 4,000 both finish exact, where reading the 1,024 mirrors stopped "either" with bounds. A
// predicate the own support cannot decide still reads the mirror: a label other carried by the
// record TTTTTTTT only, any(other), selects AAAAAA alone (the one forward walk whose mirror
// TTTTTT that record holds), found on the mirror's rows.
TEST(PatternSupportedRoute, TheMirrorIsReadOnlyWhenItDecides) {
    std::string record;
    for (char a : std::string("ACGT")) {
        for (char b : std::string("ACGT")) {
            for (char c : std::string("ACGT")) {
                record += a;
                record += b;
                record += c;
            }
        }
    }
    const Index idx = build_rd(3, { { "shared", "h0", record } });
    Ask r;
    r.patterns = { "ANNNNN" };
    r.kind = "iupac";
    r.mode = "count";
    r.strands = "forward";
    r.level = "label_intersection";
    r.predicate = "{\"any\": [\"shared\"]}";
    for (uint64_t work : { uint64_t(0), uint64_t(4000) }) {
        r.max_annotation_work = work;
        r.predicate_strands = "context";
        const Json::Value context = answer_of(idx, r)["patterns"][0];
        r.predicate_strands = "either";
        const Json::Value either = answer_of(idx, r)["patterns"][0];
        for (const Json::Value *e : { &context, &either }) {
            ASSERT_FALSE(e->isMember("error")) << compact(*e);
            EXPECT_EQ("completed", (*e)["selection"]["pass"].asString()) << work;
            EXPECT_EQ("exact", (*e)["counts"]["selected"]["relation"].asString()) << work;
            EXPECT_EQ(1024u, (*e)["counts"]["selected"]["value"].asUInt64()) << work;
            EXPECT_EQ(1024u, (*e)["counts"]["tested"]["value"].asUInt64()) << work;
            EXPECT_TRUE((*e)["stop"].isNull()) << compact((*e)["stop"]);
        }
        EXPECT_EQ(64u, context["work"]["annotation_rows"].asUInt64());
        EXPECT_EQ(context["work"]["annotation_rows"].asUInt64(),
                  either["work"]["annotation_rows"].asUInt64());
        EXPECT_EQ(0u, either["work"]["mirror_rows"].asUInt64());
        EXPECT_EQ(0u, either["work"]["predicate_lookups"].asUInt64());
        // the decisions' units: each walk's eval3 on its own support, 1 and 1 per leaf its
        // present label is listed in ("context" charges more: its monotone pruning asks the
        // predicate of every branch too)
        EXPECT_EQ(2u * 1024, either["work"]["predicate_units"].asUInt64());
        EXPECT_EQ(context["work"]["annotation_units"].asUInt64(),
                  either["work"]["annotation_units"].asUInt64());
    }

    const Index two = build_rd(3, { { "shared", "h0", record }, { "other", "o0", "TTTTTTTT" } });
    r.max_annotation_work = 0;
    r.predicate_strands.clear();
    r.predicate.clear();
    const Json::Value plain = answer_of(two, r)["patterns"][0];
    r.predicate = "{\"any\": [\"other\"]}";
    r.predicate_strands = "context";
    const Json::Value context = answer_of(two, r)["patterns"][0];
    EXPECT_EQ(0u, context["counts"]["selected"]["value"].asUInt64());
    EXPECT_EQ("exact", context["counts"]["selected"]["relation"].asString());
    r.predicate_strands = "either";
    const Json::Value either = answer_of(two, r)["patterns"][0];
    EXPECT_EQ("completed", either["selection"]["pass"].asString());
    EXPECT_EQ("exact", either["counts"]["selected"]["relation"].asString());
    EXPECT_EQ(1u, either["counts"]["selected"]["value"].asUInt64());
    // every walk's own support {shared} leaves any(other) undecided: every mirror is read --
    // the 256 walks ending in T have theirs, a forward walk, among the held walks (no lookup),
    // the other 768 are looked up and read
    EXPECT_EQ(768u, either["work"]["predicate_lookups"].asUInt64());
    EXPECT_LT(0u, either["work"]["mirror_rows"].asUInt64());
    // the search's rows are a plain search's ("context" prunes on the monotone any(other),
    // "either" holds every supported walk and prunes nothing), the mirrors' on top
    EXPECT_EQ(plain["work"]["annotation_rows"].asUInt64()
                      + either["work"]["mirror_rows"].asUInt64(),
              either["work"]["annotation_rows"].asUInt64());
}

// "either" holds the supported walks for their decisions, at most max_predicate_contexts:
// all_or_count withholds above it (not admitted), count does not admit them, partial decides
// the first ones
TEST(PatternSupportedRoute, HeldAboveMaxPredicateContexts) {
    const std::vector<Record> records = {
        { "a", "a0", "TTACGTAACGGAT" }, { "b", "b0", "TTACGTCACGGAT" },
        { "c", "c0", "TTACGTGACGGAT" }, { "d", "d0", "TTACGTTACGGAT" },
    };
    const Index idx = build_rd(5, records);
    Ask r;
    r.patterns = { "ACGTNACGG" };
    r.kind = "iupac";
    r.strands = "forward";
    r.predicate = "{\"none\": [\"a\"]}";
    r.max_predicate_contexts = 2;
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("not_admitted", e["selection"]["pass"].asString());
        EXPECT_EQ("predicate_above_threshold", e["withheld"]["reason"].asString());
        EXPECT_EQ("unknown", e["counts"]["selected"]["relation"].asString());
        EXPECT_EQ(4u, e["counts"]["supported_paths"]["value"].asUInt64());
    }
    r.mode = "count";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("not_admitted", e["selection"]["pass"].asString());
    }
    r.mode = "partial";
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("stopped", e["selection"]["pass"].asString());
        EXPECT_EQ(2u, e["counts"]["tested"]["value"].asUInt64());
        EXPECT_EQ("max_predicate_contexts", e["cut"]["reason"].asString());
        // of the first two walks (A..., C...: a's and b's) b's is selected; the other two are
        // untested: bounds over the exact supported count
        EXPECT_EQ("bounds", e["counts"]["selected"]["relation"].asString());
        EXPECT_EQ(1u, e["counts"]["selected"]["lower"].asUInt64());
        EXPECT_EQ(3u, e["counts"]["selected"]["upper"].asUInt64());
        EXPECT_EQ(1u, e["results"].size());
    }
    r.mode = "all_or_count";
    r.max_predicate_contexts = 4;
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("completed", e["selection"]["pass"].asString());
        EXPECT_EQ(3u, e["results"].size());
    }
}

// A constant normal form selects without deciding: true is the request without the predicate,
// false lists nothing (complete); an unbound or spent selection does not start
TEST(PatternSupportedRoute, ConstantPredicates) {
    const std::vector<Record> records = {
        { "a", "a0", "TTACGTAACGGAT" }, { "b", "b0", "TTACGTCACGGAT" },
    };
    const Index idx = build_rd(5, records);
    Ask r;
    r.patterns = { "ACGTNACGG" };
    r.kind = "iupac";
    r.strands = "forward";
    r.labels = "predicate_only";
    r.predicate = "{\"none\": [\"zz\"]}";       // folds to true
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("constant", e["selection"]["pass"].asString());
        EXPECT_EQ(2u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ(2u, e["results"].size());
        for (const Json::Value &res : e["results"]) {
            EXPECT_EQ(0u, res["selection_labels"].size());
            EXPECT_EQ(0u, res["labels"].size());
        }
    }
    r.predicate = "{\"any\": [\"zz\"]}";        // folds to false
    {
        const Json::Value e = answer_of(idx, r)["patterns"][0];
        EXPECT_EQ("constant", e["selection"]["pass"].asString());
        EXPECT_EQ(0u, e["counts"]["selected"]["value"].asUInt64());
        EXPECT_EQ(2u, e["counts"]["tested"]["value"].asUInt64());
        EXPECT_EQ(0u, e["results"].size());
        EXPECT_TRUE(e["retrieval_complete"].asBool());
        EXPECT_EQ("exact", e["counts"]["labels"]["relation"].asString());
    }
}


// ------------------------------------------------------------------ budgets and stops

// Every max_annotation_work, every reading of a virtual clock, every account size: the answer
// is one of the complete answer's prefixes (partial) or withheld (all_or_count), its stop
// named, its supported paths at_least the ones it lists, never more than the oracle's
TEST(PatternSupportedRoute, BudgetAndClockSweeps) {
    std::mt19937 rng(5);
    const size_t k = 5;
    const std::vector<Record> records = random_records(rng, k, 3);
    const Index idx = build_rd(k, records);
    const Oracle oracle(idx);
    std::string p;
    std::vector<Walk> full_walks;
    for (size_t t = 0; t < 50 && full_walks.size() < 3; ++t) {
        p = random_pattern(rng, records, k);
        full_walks = oracle.walks(p, "both", true, false);
    }
    ASSERT_GE(full_walks.size(), 1u);
    for (const char *predicate : { "", "{\"none\": [\"zz\", \"c0\"]}" }) {
        Ask r;
        r.patterns = { p };
        r.kind = "iupac";
        r.labels = "all";
        r.mode = "partial";
        r.predicate = predicate;
        const Json::Value full = answer_of(idx, r)["patterns"][0];
        const uint64_t total = full["work"]["annotation_units"].asUInt64();
        const auto all = result_set(full);
        auto check = [&](const Json::Value &e, const std::string &what) {
            ASSERT_FALSE(e.isMember("error"));
            const auto got = result_set(e);
            for (const auto &x : got) {
                EXPECT_TRUE(all.count(x)) << what;
            }
            const Json::Value &sp = e["counts"]["supported_paths"];
            if (sp["relation"].asString() == "exact" || sp["relation"].asString() == "at_least") {
                EXPECT_LE(sp["value"].asUInt64(), full_walks.size()) << what;
            }
            if (got.size() < all.size()) {
                EXPECT_FALSE(e["retrieval_complete"].asBool()) << what;
                EXPECT_TRUE(e["cut"].isObject() || !e["stop"].isNull()) << what << compact(e);
            }
        };
        for (uint64_t w = 1; w <= total + 1; w += std::max<uint64_t>(1, total / 60)) {
            Ask q = r;
            q.max_annotation_work = w;
            const Json::Value e = answer_of(idx, q)["patterns"][0];
            check(e, "work " + std::to_string(w));
            if (w <= total && e["stop"].isObject() && e["stop"]["phase"] == "extension") {
                EXPECT_EQ("max_annotation_work", e["stop"]["reason"].asString());
            }
            Ask s = q;
            s.mode = "all_or_count";
            const Json::Value f = answer_of(idx, s)["patterns"][0];
            if (f["stop"].isObject() && f["stop"]["reason"] == "max_annotation_work") {
                EXPECT_EQ("annotation_budget", f["withheld"]["reason"].asString());
            }
        }
        // a virtual clock: the work time passes at the n-th reading
        for (size_t n = 1; n < 400; n += 7) {
            Ask q = r;
            auto readings = std::make_shared<size_t>(0);
            const Clock::time_point t0 = Clock::now();
            q.clock = [readings, n, t0]() {
                return ++*readings < n ? t0 : t0 + std::chrono::hours(1);
            };
            const Json::Value e = answer_of(idx, q)["patterns"][0];
            check(e, "clock " + std::to_string(n));
            if (!e["stop"].isNull() && e["stop"]["reason"] == "time") {
                EXPECT_EQ("time_limited", e["determinism"].asString());
            }
        }
        // the account, from the smallest
        for (uint64_t bytes = 1 << 12; bytes < (uint64_t(1) << 28); bytes *= 2) {
            Ask q = r;
            q.max_memory_bytes = bytes;
            check(answer_of(idx, q)["patterns"][0], "memory " + std::to_string(bytes));
        }
    }
}


// ------------------------------------------------------------------ the motif (§25)

// predicate_scope "motif": each pattern of at most k bases also answered as a whole, on the
// union of its contexts' labels (each context's row and, with "either", its reverse
// complement's); the oracle scans the records' k-mers
TEST(PatternSupportedRoute, MotifAgainstTheUnion) {
    // the motif AC in A's TACGG and C's CACCC: one context each; and(any A, none C) is true for
    // A's context and false for the motif (C carries AC elsewhere)
    const size_t k = 5;
    const std::vector<Record> records = {
        { "A", "a0", "TACGG" }, { "C", "c0", "CACCC" }, { "D", "d0", "GTTTC" },
    };
    const Index idx = build_rd(k, records);
    auto union_of = [&](const std::string &motif, bool either) {
        std::set<std::string> u;
        for (const Record &r : records) {
            for (size_t i = 0; i + k <= r.seq.size(); ++i) {
                const std::string x = r.seq.substr(i, k);
                // a context of the motif or of its reverse complement in this k-mer, or (either)
                // in its reverse complement
                if (x.find(motif) != std::string::npos || x.find(rc(motif)) != std::string::npos
                        || (either && (rc(x).find(motif) != std::string::npos
                                       || rc(x).find(rc(motif)) != std::string::npos))) {
                    u.insert(r.column);
                }
            }
        }
        return u;
    };
    for (const std::string predicate : { "{\"and\": [{\"any\": [\"A\"]}, {\"none\": [\"C\"]}]}",
                                         "{\"any\": [\"A\"]}", "{\"none\": [\"D\"]}",
                                         "{\"at_least\": {\"n\": 2, \"labels\": [\"A\", \"C\", \"D\"]}}" }) {
        for (const std::string ps : { "either", "context" }) {
            Ask r;
            r.patterns = { "AC", "GT" };
            r.long_search = "anchors";
            r.mode = "count";
            r.predicate = predicate;
            r.predicate_strands = ps;
            r.predicate_scope = "motif";
            const Json::Value a = answer_of(idx, r);
            EXPECT_EQ("shard_motif", a["predicate"]["motif_scope"].asString());
            EXPECT_EQ("motif", a["limits"]["predicate_scope"].asString());
            const Json::Value pj = parse_pattern_body(predicate);
            std::set<std::string> names;
            o_names(pj, &names);
            for (Json::ArrayIndex i = 0; i < 2; ++i) {
                const Json::Value &e = a["patterns"][i];
                const std::string motif = r.patterns[i];
                const std::set<std::string> u = union_of(motif, ps == "either");
                const Json::Value &m = e["motif"];
                ASSERT_TRUE(m.isObject()) << compact(e);
                EXPECT_EQ("every_context", m["decided_by"].asString()) << compact(m);
                EXPECT_EQ(o_eval(pj, u), m["selected"].asBool()) << predicate << " " << motif;
                std::set<std::string> present;
                for (const Json::Value &l : m["labels_present"]) {
                    present.insert(l["column"].asString());
                }
                std::set<std::string> expected;
                for (const std::string &n : u) {
                    if (names.count(n))
                        expected.insert(n);
                }
                EXPECT_EQ(expected, present) << predicate << " " << motif << " " << ps;
                EXPECT_EQ(m["labels"].asUInt64() - present.size(),
                          m["labels_absent"].asUInt64());
                EXPECT_TRUE(e["work"].isMember("motif_units"));
            }
        }
    }
    // without predicate_scope the answer has no motif block, and no new limits field
    Ask r;
    r.patterns = { "AC" };
    r.long_search = "anchors";
    r.mode = "count";
    r.predicate = "{\"any\": [\"A\"]}";
    const Json::Value a = answer_of(idx, r);
    EXPECT_FALSE(a["patterns"][0].isMember("motif"));
    EXPECT_FALSE(a["predicate"].isMember("motif_scope"));
    EXPECT_FALSE(a["limits"].isMember("predicate_scope"));
}

// A pattern without an instance on the graph (no context of a short one, no anchor of a long
// one) is decided on the empty, complete union and says so: no_instance, selected the normal
// form on the empty set, every label absent, nothing untested, under both long searches; a
// constant normal form is constant there too; a pattern with a context stays every_context.
TEST(PatternSupportedRoute, MotifNoInstance) {
    const size_t k = 5;
    const Index idx = build_rd(k, { { "A", "a0", "TACGG" }, { "C", "c0", "CACCC" },
                                    { "D", "d0", "GTTTC" } });
    for (const auto &[predicate, value] : std::vector<std::pair<std::string, bool>>{
            { "{\"any\": [\"A\"]}", false }, { "{\"none\": [\"A\"]}", true },
            { "{\"and\": [{\"any\": [\"A\"]}, {\"none\": [\"C\"]}]}", false },
            { "{\"at_least\": {\"n\": 1, \"labels\": [\"A\", \"C\", \"D\"]}}", false } }) {
        for (const std::string long_search : { "anchors", "supported_paths" }) {
            Ask r;
            // no k-mer of the records holds GAG or its reverse complement CTC
            r.patterns = { "GAG", "GAGAGAGAGA", "AC" };
            r.long_search = long_search;
            r.mode = "count";
            r.predicate = predicate;
            r.predicate_scope = "motif";
            const Json::Value a = answer_of(idx, r);
            SCOPED_TRACE(predicate + " " + long_search);
            std::set<std::string> names;
            o_names(parse_pattern_body(predicate), &names);
            const Json::Value &s = a["patterns"][0];
            const Json::Value &l = a["patterns"][1];
            EXPECT_EQ("exact", s["counts"]["contexts"]["relation"].asString()) << compact(s);
            EXPECT_EQ(0u, s["counts"]["contexts"]["value"].asUInt64());
            EXPECT_EQ("exact", l["counts"]["anchors"]["relation"].asString()) << compact(l);
            EXPECT_EQ(0u, l["counts"]["anchors"]["value"].asUInt64());
            for (const Json::Value *e : { &s, &l }) {
                const Json::Value &m = (*e)["motif"];
                ASSERT_TRUE(m.isObject()) << compact(*e);
                EXPECT_EQ(Json::Value(value), m["selected"]) << compact(m);
                EXPECT_EQ("no_instance", m["decided_by"].asString()) << compact(m);
                EXPECT_TRUE(m["untested"].isNull()) << compact(m);
                EXPECT_EQ(names.size(), m["labels"].asUInt64());
                EXPECT_TRUE(m["labels_present"].isArray() && m["labels_present"].empty());
                EXPECT_EQ(names.size(), m["labels_absent"].asUInt64());
                EXPECT_TRUE(m["stop"].isNull());
            }
            EXPECT_EQ("every_context", a["patterns"][2]["motif"]["decided_by"].asString())
                    << compact(a["patterns"][2]);
        }
    }
    // a constant normal form (Z is no column) stays constant on the empty union
    Ask r;
    r.patterns = { "GAG", "GAGAGAGAGA" };
    r.mode = "count";
    r.predicate = "{\"any\": [\"Z\"]}";
    r.predicate_scope = "motif";
    const Json::Value a = answer_of(idx, r);
    ASSERT_EQ(2u, a["patterns"].size());
    for (const Json::Value &e : a["patterns"]) {
        EXPECT_EQ("constant", e["motif"]["decided_by"].asString()) << compact(e);
        EXPECT_EQ(false, e["motif"]["selected"].asBool());
        EXPECT_EQ(0u, e["motif"]["labels_absent"].asUInt64());
    }
}

} // namespace
