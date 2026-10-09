#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../graph/alignment/pattern_test_support.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// Patterns longer than k as paths (long_search "paths", increment 4 of
// docs/DESIGN-pattern-search.md, §4.2, §4.3 "Label consistency for long"; owner decisions #13
// and #14): the route's path results (sequence, anchor_kmer, nodes, rows; never kmer), the two
// thresholds (max_anchors admits the extension, max_paths the release), the counts with their
// relations, and the labels of each path with its per-label support (label_intersection:
// every k-mer of the path annotated; record_verified: one record holds the whole path) and
// require_support. Tiny graphs built from explicit records in labelled columns with their
// coordinates and record mapping. The expectations come from two oracles that never ask the
// engine:
//  - a GRAPH-WALK oracle over the k-mer set of the records (as deposited on BASIC; with their
//    reverse complements on CANONICAL and PRIMARY graphs): every string of L bases that
//    instantiates the oriented pattern (the test's own IUPAC table) and whose every k-window
//    is such a k-mer, found by a depth-first walk over strings;
//  - a RECORD-SCAN oracle over the records: a column carries a path when every k-mer of it is
//    in one of its records; it verifies it when one of its records holds the whole path, at
//    the 1-based positions of the scan.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;

// ------------------------------------------------------------------ the test's own tables

const std::map<char, std::string> kCodes = {
    { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" }, { 'R', "AG" }, { 'Y', "CT" },
    { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" }, { 'B', "CGT" },
    { 'D', "AGT" }, { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
};

char complement_code(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'V': return 'B';
        case 'D': return 'H';
        case 'H': return 'D';
        default: return c;     // S, W, N
    }
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

// |s| instantiates |q| position by position
bool instantiates(const std::string &q, const std::string &s) {
    if (q.size() != s.size())
        return false;
    for (size_t i = 0; i < q.size(); ++i) {
        if (!admits(q[i], s[i]))
            return false;
    }
    return true;
}

// the oriented patterns of |p|: (+, P) and (-, rc(P)), or (=, P) for a palindrome; on a graph
// without strands, their orientation names
std::vector<std::pair<std::string, std::string>> oriented(const std::string &p,
                                                          bool strand_stated = true) {
    if (rc(p) == p)
        return { { strand_stated ? "=" : "palindromic", p } };
    return { { strand_stated ? "+" : "forward", p }, { strand_stated ? "-" : "reverse", rc(p) } };
}

// ------------------------------------------------------------------ the index

struct Record {
    std::string column;
    std::string header;
    std::string seq;
};

struct Index {
    size_t k = 0;
    DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
    std::vector<uint64_t> starts;
};

// as test_pattern_retrieval.cpp builds its indexes: each record's k-mers numbered from its
// column's running count, the record mapping from the same records
template <class Annotation>
Index build(size_t k, const std::vector<Record> &records, bool coordinates,
            DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    Index idx;
    idx.k = k;
    idx.mode = mode;
    idx.records = records;
    std::vector<std::string> seqs, labels;
    std::map<std::string, uint64_t> next;
    for (const Record &r : records) {
        EXPECT_GE(r.seq.size(), k);
        seqs.push_back(r.seq);
        labels.push_back(r.column);
        idx.starts.push_back(next[r.column]);
        next[r.column] += r.seq.size() - k + 1;
    }
    idx.anno = test::build_anno_graph<DBGSuccinct, Annotation>(
            k, seqs, labels, mode, coordinates,
            coordinates ? idx.starts : std::vector<uint64_t>{});
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

// ------------------------------------------------------------------ the oracles

/**
 * The graph-walk oracle: the k-mers of the records (with their reverse complements when the
 * graph holds both orientations) and the paths of an oriented pattern over them.
 */
struct Walk {
    size_t k = 0;
    std::set<std::string> kmers;

    explicit Walk(const Index &idx) : k(idx.k) {
        const bool both = idx.mode != DeBruijnGraph::BASIC;
        for (const Record &r : idx.records) {
            for (size_t i = 0; i + k <= r.seq.size(); ++i) {
                kmers.insert(r.seq.substr(i, k));
                if (both)
                    kmers.insert(rc(r.seq.substr(i, k)));
            }
        }
    }

    // the k-mers instantiating q[0, k)
    std::vector<std::string> anchors(const std::string &q) const {
        std::vector<std::string> out;
        for (const std::string &x : kmers) {
            if (instantiates(q.substr(0, k), x))
                out.push_back(x);
        }
        return out;
    }

    /**
     * Every string of |q|'s length that instantiates q and whose every k-window is a k-mer,
     * found by extending each anchor one base at a time; |candidates| counts the prefixes of
     * k + 1 .. L bases so formed (the branches the engine's DFS enters, §4.2).
     */
    std::vector<std::string> paths(const std::string &q, uint64_t *candidates = nullptr) const {
        std::vector<std::string> out;
        uint64_t entered = 0;
        std::function<void(const std::string&)> dfs = [&](const std::string &s) {
            if (s.size() == q.size()) {
                out.push_back(s);
                return;
            }
            for (char b : std::string("ACGT")) {
                if (!admits(q[s.size()], b))
                    continue;
                if (!kmers.count(s.substr(s.size() - k + 1) + b))
                    continue;
                ++entered;
                dfs(s + b);
            }
        };
        for (const std::string &x : anchors(q)) {
            dfs(x);
        }
        if (candidates)
            *candidates = entered;
        return out;
    }
};

/**
 * The record-scan oracle: per column its records (seq_id order); which columns carry a path
 * (every k-mer of it in one of their records, either orientation on a graph without strands),
 * and where a column's records hold the whole path (1-based starts).
 */
struct Scan {
    size_t k = 0;
    bool both = false;
    std::map<std::string, std::vector<const Record*>> columns;

    explicit Scan(const Index &idx) : k(idx.k), both(idx.mode != DeBruijnGraph::BASIC) {
        for (const Record &r : idx.records) {
            columns[r.column].push_back(&r);
        }
    }

    bool holds_kmer(const std::string &column, const std::string &kmer) const {
        for (const Record *r : columns.at(column)) {
            if (r->seq.find(kmer) != std::string::npos)
                return true;
            if (both && r->seq.find(rc(kmer)) != std::string::npos)
                return true;
        }
        return false;
    }

    std::set<std::string> carriers(const std::string &path) const {
        std::set<std::string> out;
        for (const auto &[column, recs] : columns) {
            bool all = true;
            for (size_t i = 0; i + k <= path.size() && all; ++i) {
                all = holds_kmer(column, path.substr(i, k));
            }
            if (all)
                out.insert(column);
        }
        return out;
    }

    // {(seq_id, 1-based start)} of |path| in the records of |column|
    std::set<std::pair<uint64_t, uint64_t>> occurrences(const std::string &column,
                                                        const std::string &path) const {
        std::set<std::pair<uint64_t, uint64_t>> out;
        const auto &recs = columns.at(column);
        for (uint64_t id = 0; id < recs.size(); ++id) {
            const std::string &seq = recs[id]->seq;
            for (size_t i = 0; i + path.size() <= seq.size(); ++i) {
                if (seq.compare(i, path.size(), path) == 0)
                    out.emplace(id, i + 1);
            }
        }
        return out;
    }
};

// ------------------------------------------------------------------ running requests

PatternLimits limits() {
    PatternLimits l;
    // the short anchor windows of the tiny graphs are below any useful floor
    l.min_information_bits = 0;
    return l;
}

Json::Value run(const Index &idx, const std::string &body, RetrievalHooks hooks = {},
                bool records = true, const PatternLimits &caps = limits()) {
    if (records)
        hooks.coord_to_header = idx.cth.get();
    return process_pattern_request(parse_pattern_body(body), *idx.anno, caps, "rel",
                                   nullptr, nullptr, nullptr, &hooks);
}

using Clock = pattern::Deadline::Clock;

// as run(), the request's clock |clock| (a virtual one the test moves)
Json::Value run_clocked(const Index &idx, const std::string &body, RetrievalHooks hooks,
                        const std::function<Clock::time_point()> &clock) {
    hooks.coord_to_header = idx.cth.get();
    return process_pattern_request(parse_pattern_body(body), *idx.anno, limits(), "rel",
                                   nullptr, nullptr, clock, &hooks);
}

std::pair<int, std::string> refusal(const Index &idx, const std::string &body,
                                    bool records = true) {
    try {
        run(idx, body, {}, records);
    } catch (const PatternRefusal &e) {
        return { e.status(), e.code() };
    }
    return { 200, "" };
}

std::string body(const std::string &patterns, const std::string &rest = "",
                 const std::string &labels = "all") {
    return "{\"patterns\": [" + patterns + "], \"long_search\": \"paths\", "
           "\"output\": {\"labels\": \"" + labels + "\"}" + (rest.empty() ? "" : ", " + rest)
           + "}";
}

// an answer without its timings
Json::Value untimed(Json::Value v) {
    if (v.isObject()) {
        v.removeMember("timing");
        for (const std::string &name : v.getMemberNames()) {
            v[name] = untimed(v[name]);
        }
    } else if (v.isArray()) {
        for (Json::ArrayIndex i = 0; i < v.size(); ++i) {
            v[i] = untimed(v[i]);
        }
    }
    return v;
}

// (strand or orientation, sequence) of paths
using PathSet = std::set<std::pair<std::string, std::string>>;

uint8_t orientation_rank(const std::string &s) {
    if (s == "+" || s == "forward")
        return 0;
    if (s == "-" || s == "reverse")
        return 1;
    return 2;
}

/**
 * The label-free shape of every path result against the graph-walk oracle: exactly the
 * oracle's paths, each with sequence = instance, anchor_kmer its first k bases as the graph
 * spells its first node, offset 0, its n nodes spelling its k-windows, its rows those of the
 * stored (or canonical) k-mers, no kmer; in the answer order (anchor node, orientation,
 * sequence). Returns the (strand, sequence) of every result, in order.
 */
std::vector<std::pair<std::string, std::string>>
check_path_results(const Index &idx, const Json::Value &e, const std::string &p) {
    const size_t k = idx.k, L = p.size(), n = L - k + 1;
    const bool strand_stated = idx.mode == DeBruijnGraph::BASIC;
    const std::string key = strand_stated ? "strand" : "orientation";
    const DeBruijnGraph &graph = idx.anno->get_graph();
    Walk walk(idx);
    PathSet expected;
    for (const auto &[strand, q] : oriented(p, strand_stated)) {
        for (const std::string &s : walk.paths(q)) {
            expected.emplace(strand, s);
        }
    }
    std::vector<std::pair<std::string, std::string>> got;
    std::vector<std::tuple<uint64_t, uint8_t, std::string>> order;
    for (const Json::Value &r : e["results"]) {
        EXPECT_FALSE(r.isMember("kmer")) << r;
        EXPECT_FALSE(r.isMember("node")) << r;
        EXPECT_FALSE(r.isMember("row")) << r;
        const std::string s = r["sequence"].asString();
        EXPECT_EQ(L, s.size());
        EXPECT_EQ(s, r["instance"].asString());
        EXPECT_EQ(s.substr(0, k), r["anchor_kmer"].asString());
        EXPECT_EQ(0u, r["offset"].asUInt64());
        EXPECT_TRUE(r.isMember(key)) << r;
        EXPECT_FALSE(r.isMember(strand_stated ? "orientation" : "strand")) << r;
        const Json::Value &nodes = r["nodes"];
        const Json::Value &rows = r["rows"];
        EXPECT_EQ(n, nodes.size());
        EXPECT_EQ(n, rows.size());
        for (Json::ArrayIndex j = 0; j < nodes.size() && j < n; ++j) {
            const uint64_t node = nodes[j].asUInt64();
            EXPECT_EQ(s.substr(j, k), graph.get_node_sequence(node)) << j << " " << r;
            EXPECT_TRUE(rows[j].isUInt64()) << r;
            if (!rows[j].isUInt64())
                continue;
            if (idx.mode == DeBruijnGraph::BASIC) {
                EXPECT_EQ(node - 1, rows[j].asUInt64());
            } else if (idx.mode == DeBruijnGraph::CANONICAL) {
                // the canonical k-mer's row, shared with the reverse complement
                DeBruijnGraph::node_index canonical = DeBruijnGraph::npos;
                graph.map_to_nodes(s.substr(j, k), [&](DeBruijnGraph::node_index x) {
                    canonical = x;
                });
                EXPECT_EQ(canonical - 1, rows[j].asUInt64());
            }
            EXPECT_LT(rows[j].asUInt64(), idx.anno->get_annotator().num_objects());
        }
        got.emplace_back(r[key].asString(), s);
        order.emplace_back(nodes[0].asUInt64(), orientation_rank(r[key].asString()), s);
    }
    EXPECT_EQ(expected, PathSet(got.begin(), got.end()));
    EXPECT_EQ(got.size(), PathSet(got.begin(), got.end()).size());
    EXPECT_TRUE(std::is_sorted(order.begin(), order.end()));
    EXPECT_EQ(got.size(), e["returned"].asUInt64());
    return got;
}

/**
 * counts.anchors and counts.paths of a completed extension against the graph-walk oracle:
 * per strand, exact, with the branches the extension entered and what it did.
 */
void check_path_counts(const Index &idx, const Json::Value &e, const std::string &p) {
    const bool strand_stated = idx.mode == DeBruijnGraph::BASIC;
    const std::string by = strand_stated ? "by_strand" : "by_orientation";
    Walk walk(idx);
    uint64_t anchors = 0, paths = 0, candidates = 0;
    const Json::Value &a = e["counts"]["anchors"];
    const Json::Value &c = e["counts"]["paths"];
    for (const auto &[strand, q] : oriented(p, strand_stated)) {
        const std::string name = strand == "=" ? "both" : strand;
        uint64_t entered = 0;
        const uint64_t found = walk.paths(q, &entered).size();
        anchors += walk.anchors(q).size();
        paths += found;
        candidates += entered;
        EXPECT_EQ(walk.anchors(q).size(), a[by][name]["value"].asUInt64()) << name;
        EXPECT_EQ(found, c[by][name]["value"].asUInt64()) << name << " " << c;
        EXPECT_EQ("exact", c[by][name]["relation"].asString());
        EXPECT_EQ("paths", c[by][name]["unit"].asString());
    }
    // a palindromic anchor window on a wrapped PRIMARY graph is one anchor, not two
    if (idx.mode != DeBruijnGraph::PRIMARY) {
        EXPECT_EQ(anchors, a["value"].asUInt64());
    }
    EXPECT_EQ("exact", a["relation"].asString());
    EXPECT_EQ(paths, c["value"].asUInt64()) << c;
    EXPECT_EQ("exact", c["relation"].asString());
    EXPECT_EQ("paths", c["unit"].asString());
    if (idx.mode == DeBruijnGraph::BASIC) {
        EXPECT_EQ(candidates, c["candidates_examined"].asUInt64()) << c;
    }
    EXPECT_EQ(anchors ? "completed" : "no_anchors", c["extension"].asString());
    EXPECT_FALSE(e.isMember("error"));
    EXPECT_EQ("long", e["scope"].asString());
    EXPECT_EQ("long", e["absence_scope"].asString());
    for (const Json::Value &note : e["notes"]) {
        EXPECT_NE("paths_later_increment", note.asString());
    }
    EXPECT_TRUE(e["work"].isMember("extension_edges"));
    EXPECT_LE(e["work"]["extension_edges"].asUInt64(), e["work"]["steps"].asUInt64());
    EXPECT_TRUE(e["timing"].isMember("extension_ms"));
}

// {(column, seq_id, start, strand)}: every occurrence of the oriented pattern in every record
std::set<std::tuple<std::string, uint64_t, uint64_t, std::string>>
placed_scan(const Index &idx, const std::string &p) {
    std::set<std::tuple<std::string, uint64_t, uint64_t, std::string>> out;
    std::map<std::string, uint64_t> seq_id;
    for (const Record &r : idx.records) {
        const uint64_t id = seq_id[r.column]++;
        for (const auto &[strand, q] : oriented(p)) {
            for (size_t i = 0; i + q.size() <= r.seq.size(); ++i) {
                if (instantiates(q, r.seq.substr(i, q.size())))
                    out.emplace(r.column, id, i + 1, strand);
            }
        }
    }
    return out;
}

/**
 * A complete labelled answer with record placement (BASIC) against both oracles: each path's
 * labels exactly the columns carrying it (or, with |require_verified|, those verifying it),
 * each with its support and its placed occurrences of the whole path; by_label and every
 * count exact, the support split, and the occurrences the record scan's.
 */
void check_labelled_paths(const Index &idx, const Json::Value &e, const std::string &p,
                          bool require_verified) {
    const size_t L = p.size();
    Scan scan(idx);
    check_path_counts(idx, e, p);
    const auto results = check_path_results(idx, e, p);
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    EXPECT_TRUE(e["withheld"].isNull());
    EXPECT_TRUE(e["stop"].isNull());
    EXPECT_TRUE(e["cut"].isNull());
    EXPECT_EQ("record", e["placement"].asString());
    EXPECT_EQ("budgeted", e["annotation"].asString());
    EXPECT_EQ(0u, e["rows_refused"].size());
    EXPECT_EQ(0u, e["anchors_truncated"].size());

    std::map<std::string, uint64_t> label_paths, label_verified;
    std::map<std::string, std::set<std::tuple<uint64_t, uint64_t, std::string>>> unions;
    std::set<std::string> carriers_all, verified_any;
    for (Json::ArrayIndex i = 0; i < e["results"].size(); ++i) {
        const Json::Value &r = e["results"][i];
        const std::string strand = results[i].first;
        const std::string s = results[i].second;
        const std::set<std::string> carriers = scan.carriers(s);
        std::set<std::string> verified;
        for (const std::string &c : carriers) {
            if (!scan.occurrences(c, s).empty())
                verified.insert(c);
        }
        carriers_all.insert(carriers.begin(), carriers.end());
        const std::set<std::string> &listed = require_verified ? verified : carriers;
        EXPECT_EQ("complete", r["labels_status"].asString()) << r;
        EXPECT_EQ(carriers.size(), r["labels_total"].asUInt64()) << r;
        if (require_verified) {
            EXPECT_EQ(carriers.size() - verified.size(),
                      r["labels_excluded_unverified"].asUInt64()) << r;
        } else {
            EXPECT_FALSE(r.isMember("labels_excluded_unverified"));
        }
        const std::string summary = listed.empty() ? ""
                : verified.size() == listed.size() ? "record_verified"
                : verified.empty() ? "label_intersection" : "mixed";
        if (summary.empty()) {
            EXPECT_TRUE(r["support"].isNull()) << r;
        } else {
            EXPECT_EQ(summary, r["support"].asString()) << r;
        }
        std::set<std::string> got;
        for (const Json::Value &l : r["labels"]) {
            const std::string column = l["column"].asString();
            got.insert(column);
            const auto occ = scan.occurrences(column, s);
            EXPECT_EQ(occ.empty() ? "label_intersection" : "record_verified",
                      l["support"].asString()) << s << " " << column;
            EXPECT_EQ(occ.size(), l["occurrences"]["value"].asUInt64()) << s << " " << column;
            EXPECT_EQ("exact", l["occurrences"]["relation"].asString());
            std::set<std::pair<uint64_t, uint64_t>> listed_occ;
            for (const Json::Value &o : l["occurrence_list"]) {
                const uint64_t seq_id = o["seq_id"].asUInt64();
                const Record &rec = *scan.columns.at(column).at(seq_id);
                EXPECT_EQ(rec.header, o["record"].asString());
                EXPECT_EQ(rec.seq.size(), o["nt_length"].asUInt64());
                EXPECT_EQ(strand, o["strand"].asString());
                const std::string coords = o["nt_coords"].asString();
                const uint64_t start = std::stoull(coords.substr(0, coords.find('-')));
                EXPECT_EQ(std::to_string(start) + "-" + std::to_string(start + L - 1), coords);
                listed_occ.emplace(seq_id, start);
                unions[column].emplace(seq_id, start, strand);
            }
            EXPECT_EQ(occ, listed_occ) << s << " " << column;
            label_paths[column]++;
            if (!occ.empty()) {
                label_verified[column]++;
                verified_any.insert(column);
            }
        }
        EXPECT_EQ(listed, got) << s;
    }

    // by_label: every listed column once, (paths desc, column asc), exact
    const Json::Value &by_label = e["by_label"];
    ASSERT_EQ(label_paths.size(), by_label.size()) << by_label;
    std::vector<std::pair<int64_t, std::string>> keys;
    uint64_t occurrences = 0;
    for (const Json::Value &b : by_label) {
        const std::string column = b["column"].asString();
        keys.emplace_back(-b["paths"]["value"].asInt64(), column);
        EXPECT_FALSE(b.isMember("contexts"));
        EXPECT_EQ(label_paths[column], b["paths"]["value"].asUInt64()) << column;
        EXPECT_EQ("exact", b["paths"]["relation"].asString());
        EXPECT_EQ("paths", b["paths"]["unit"].asString());
        EXPECT_EQ(label_verified[column], b["paths_record_verified"]["value"].asUInt64());
        EXPECT_EQ("exact", b["paths_record_verified"]["relation"].asString());
        EXPECT_EQ(unions[column].size(), b["occurrences"]["value"].asUInt64()) << column;
        EXPECT_EQ("exact", b["occurrences"]["relation"].asString());
        occurrences += unions[column].size();
    }
    EXPECT_TRUE(std::is_sorted(keys.begin(), keys.end()));

    const Json::Value &labels = e["counts"]["labels"];
    EXPECT_EQ(label_paths.size(), labels["value"].asUInt64());
    EXPECT_EQ("exact", labels["relation"].asString());
    EXPECT_EQ(verified_any.size(), labels["by_support"]["record_verified"]["value"].asUInt64());
    EXPECT_EQ("exact", labels["by_support"]["record_verified"]["relation"].asString());
    EXPECT_EQ(label_paths.size() - verified_any.size(),
              labels["by_support"]["label_intersection"]["value"].asUInt64());
    EXPECT_EQ("exact", labels["by_support"]["label_intersection"]["relation"].asString());
    EXPECT_EQ(occurrences, e["counts"]["occurrences"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["occurrences"]["relation"].asString());
    if (require_verified) {
        uint64_t excluded = 0;
        for (const std::string &c : carriers_all) {
            excluded += !verified_any.count(c);
        }
        EXPECT_EQ(excluded, e["labels_excluded_unverified"]["value"].asUInt64());
        EXPECT_EQ("exact", e["labels_excluded_unverified"]["relation"].asString());
    } else {
        EXPECT_FALSE(e.isMember("labels_excluded_unverified"));
    }

    // the placed occurrences: the record scan's, each once
    std::set<std::tuple<std::string, uint64_t, uint64_t, std::string>> placed;
    for (const auto &[column, set] : unions) {
        for (const auto &[seq_id, start, strand] : set) {
            placed.emplace(column, seq_id, start, strand);
        }
    }
    EXPECT_EQ(placed_scan(idx, p), placed);
}

// k = 5. The pattern ACGTACC (n = 3 k-mers: ACGTA, CGTAC, GTACC):
//  - column A holds it whole in its record a0 (verified, at 3-9);
//  - column B holds ACGTA as the last k-mer of b0 (local 3) and CGTAC, GTACC as the first
//    two of b1: consecutive column coordinates 3, 4, 5 crossing from b0 into b1, which the
//    record bounds must reject (label_intersection), as DESIGN §4.3 "Label consistency";
//  - CGTACC (n = 2) is held whole by a0 (at 4-9) and b1 (at 1-6): verified in both;
//  - ACCAAG: its anchor ACCAA (a0's last k-mer) has no outgoing k-mer: no path;
//  - GGACGTACCA (n = 6), a0's first ten bases: in A only.
const size_t kK = 5;
const std::vector<Record> kRecords = {
    { "A", "a0", "GGACGTACCAA" },
    { "B", "b0", "TTTACGTA" },
    { "B", "b1", "CGTACCGG" },
    { "C", "c0", "TTGCATGCAT" },
};


TEST(PatternPaths, TwoAndThreeKmersAgainstBothOracles) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    for (const std::string p : { "CGTACC", "ACGTACC", "GGACGTACCA" }) {
        Json::Value out = run(idx, body("{\"dna\": \"" + p + "\"}"));
        EXPECT_EQ("paths", out["limits"]["long_search"].asString());
        EXPECT_EQ(1000u, out["limits"]["max_paths"].asUInt64());
        EXPECT_EQ("label_intersection", out["limits"]["require_support"].asString());
        const Json::Value &e = out["patterns"][0];
        ASSERT_GT(e["results"].size(), 0u) << p;
        check_labelled_paths(idx, e, p, false);
        // the same paths without labels: nothing read, the same results' shape
        Json::Value none = run(idx, body("{\"dna\": \"" + p + "\"}", "", "none"));
        const Json::Value &f = none["patterns"][0];
        check_path_counts(idx, f, p);
        EXPECT_EQ(e["returned"], f["returned"]);
        EXPECT_TRUE(f["retrieval_complete"].asBool());
        check_path_results(idx, f, p);
        for (const Json::Value &r : f["results"]) {
            EXPECT_EQ(7u, r.size()) << r;     // sequence, anchor_kmer, instance, offset,
                                              // strand, nodes, rows
        }
        EXPECT_FALSE(none["limits"].isMember("require_support"));
        EXPECT_EQ("unknown", f["counts"]["labels"]["relation"].asString());
    }
}

TEST(PatternPaths, CrossRecordPathIsNotRecordVerified) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"ACGTACC\"}"));
    const Json::Value &e = out["patterns"][0];
    check_labelled_paths(idx, e, "ACGTACC", false);
    ASSERT_EQ(1u, e["results"].size());
    const Json::Value &r = e["results"][0];
    EXPECT_EQ("mixed", r["support"].asString());
    ASSERT_EQ(2u, r["labels"].size());
    std::map<std::string, std::string> support;
    for (const Json::Value &l : r["labels"]) {
        support[l["column"].asString()] = l["support"].asString();
    }
    const std::map<std::string, std::string> expected {
        { "A", "record_verified" }, { "B", "label_intersection" },
    };
    EXPECT_EQ(expected, support);
    // B's coordinates are consecutive (3, 4, 5) but cross from b0 into b1: no occurrence
    for (const Json::Value &l : r["labels"]) {
        if (l["column"].asString() == "B") {
            EXPECT_EQ(0u, l["occurrences"]["value"].asUInt64());
            EXPECT_EQ(0u, l["occurrence_list"].size());
        } else {
            ASSERT_EQ(1u, l["occurrence_list"].size());
            EXPECT_EQ("3-9", l["occurrence_list"][0]["nt_coords"].asString());
            EXPECT_EQ("a0", l["occurrence_list"][0]["record"].asString());
        }
    }
    EXPECT_EQ(1u, e["counts"]["labels"]["by_support"]["record_verified"]["value"].asUInt64());
    EXPECT_EQ(1u, e["counts"]["labels"]["by_support"]["label_intersection"]["value"].asUInt64());
    EXPECT_EQ(1u, e["counts"]["occurrences"]["value"].asUInt64());
}

TEST(PatternPaths, RequireSupportListsTheVerifiedLabels) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    for (const std::string p : { "ACGTACC", "CGTACC", "GGACGTACCA" }) {
        Json::Value out = run(idx, body("{\"dna\": \"" + p + "\"}",
                                        "\"require_support\": \"record_verified\""));
        EXPECT_EQ("record_verified", out["limits"]["require_support"].asString());
        check_labelled_paths(idx, out["patterns"][0], p, true);
    }
    Json::Value out = run(idx, body("{\"dna\": \"ACGTACC\"}",
                                    "\"require_support\": \"record_verified\""));
    const Json::Value &e = out["patterns"][0];
    ASSERT_EQ(1u, e["results"].size());
    // B carries the path but verifies none: left out, counted
    ASSERT_EQ(1u, e["results"][0]["labels"].size());
    EXPECT_EQ("A", e["results"][0]["labels"][0]["column"].asString());
    EXPECT_EQ("record_verified", e["results"][0]["support"].asString());
    EXPECT_EQ(1u, e["results"][0]["labels_excluded_unverified"].asUInt64());
    EXPECT_EQ(2u, e["results"][0]["labels_total"].asUInt64());
    ASSERT_EQ(1u, e["by_label"].size());
    EXPECT_EQ("A", e["by_label"][0]["column"].asString());
    EXPECT_EQ(1u, e["labels_excluded_unverified"]["value"].asUInt64());
    EXPECT_EQ(1u, e["counts"]["labels"]["value"].asUInt64());
    EXPECT_EQ(0u, e["counts"]["labels"]["by_support"]["label_intersection"]["value"].asUInt64());
    // an L <= k pattern in the same request is not touched (its labels are the k-mer's)
    out = run(idx, body("{\"dna\": \"ACG\"}, {\"dna\": \"ACGTACC\"}",
                        "\"require_support\": \"record_verified\""));
    Json::Value plain = run(idx, "{\"patterns\": [{\"dna\": \"ACG\"}], \"output\": "
                                 "{\"labels\": \"all\"}}");
    EXPECT_EQ(untimed(plain["patterns"][0]), untimed(out["patterns"][0]));
}

TEST(PatternPaths, AnAnchorWithoutAPath) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    // ACCAA is a0's last k-mer: an anchor of ACCAAG with no outgoing k-mer
    for (const std::string mode : { "all_or_count", "partial", "count" }) {
        Json::Value out = run(idx, body("{\"dna\": \"ACCAAG\"}", "\"mode\": \"" + mode + "\""));
        const Json::Value &e = out["patterns"][0];
        check_path_counts(idx, e, "ACCAAG");
        EXPECT_EQ(1u, e["counts"]["anchors"]["value"].asUInt64());
        EXPECT_EQ(0u, e["counts"]["paths"]["value"].asUInt64());
        EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
        EXPECT_EQ("completed", e["counts"]["paths"]["extension"].asString());
        EXPECT_EQ(0u, e["counts"]["paths"]["candidates_examined"].asUInt64());
        if (mode == "count") {
            EXPECT_FALSE(e["retrieval_complete"].asBool());
            EXPECT_FALSE(e.isMember("results"));
            continue;
        }
        // an exact 0: the empty answer is complete, its labels exact zeros
        EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
        EXPECT_EQ(0u, e["results"].size());
        EXPECT_TRUE(e["withheld"].isNull());
        EXPECT_EQ(0u, e["counts"]["labels"]["value"].asUInt64());
        EXPECT_EQ("exact", e["counts"]["labels"]["relation"].asString());
        EXPECT_EQ("exact", e["counts"]["labels"]["by_support"]["record_verified"]["relation"]
                                   .asString());
        EXPECT_EQ("exact", e["counts"]["occurrences"]["relation"].asString());
        EXPECT_EQ(0u, e["by_label"].size());
        EXPECT_EQ(0u, e["work"]["annotation_rows"].asUInt64());
    }
    // without the option: anchors counted, paths unknown, as before
    Json::Value out = run(idx, "{\"patterns\": [{\"dna\": \"ACCAAG\"}]}");
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("unknown", e["counts"]["paths"]["relation"].asString());
    EXPECT_EQ("paths_later_increment", e["withheld"]["reason"].asString());
}

// Owner decision #13: the paths are opt-in; a request without long_search "paths" — or with
// "anchors", its default — is answered as before, and "paths" changes nothing for a pattern of
// at most k bases
TEST(PatternPaths, WithoutTheOptionTheAnswerIsUnchanged) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    for (const std::string rest : { "", ", \"mode\": \"count\"", ", \"mode\": \"partial\"",
                                    ", \"output\": {\"labels\": \"all\"}",
                                    ", \"output\": {\"labels\": \"none\", \"paths\": false}" }) {
        const std::string patterns = "\"patterns\": [{\"dna\": \"ACGTACC\"}, {\"dna\": \"ACG\"}, "
                                     "{\"iupac\": \"GTACN\"}]";
        Json::Value plain = run(idx, "{" + patterns + rest + "}");
        Json::Value anchors = run(idx, "{" + patterns + rest + ", \"long_search\": "
                                       "\"anchors\"}");
        EXPECT_EQ(untimed(plain), untimed(anchors)) << rest;
        const Json::Value &e = plain["patterns"][0];
        EXPECT_EQ("unknown", e["counts"]["paths"]["relation"].asString());
        EXPECT_FALSE(e["counts"]["paths"].isMember("extension"));
        EXPECT_FALSE(e["work"].isMember("extension_edges"));
        EXPECT_FALSE(plain["limits"].isMember("long_search"));
        EXPECT_FALSE(plain["limits"].isMember("max_paths"));
        bool later = false;
        for (const Json::Value &note : e["notes"]) {
            later |= note.asString() == "paths_later_increment";
        }
        EXPECT_TRUE(later);
        // long_search "paths": the patterns of at most k bases answer alike (in a request of
        // their own: beside a long pattern's path labels, the request's one memory account
        // would state its own peak in their memory_bytes)
        const std::string shorter = "\"patterns\": [{\"dna\": \"ACG\"}, {\"iupac\": \"GTACN\"}]";
        Json::Value without = run(idx, "{" + shorter + rest + "}");
        Json::Value paths = run(idx, "{" + shorter + rest + ", \"long_search\": \"paths\"}");
        EXPECT_EQ(untimed(without["patterns"]), untimed(paths["patterns"])) << rest;
        EXPECT_EQ(untimed(without["index"]), untimed(paths["index"]));
        EXPECT_EQ(untimed(without["output"]), untimed(paths["output"]));
        // the limits gain the option and its threshold (and the support with labels read)
        Json::Value limits = paths["limits"];
        EXPECT_EQ("paths", limits["long_search"].asString());
        limits.removeMember("long_search");
        limits.removeMember("max_paths");
        limits.removeMember("require_support");
        EXPECT_EQ(without["limits"], limits) << rest;
    }
    // max_paths without the option: bounded, clamped and listed like the other caps, no effect
    PatternLimits caps = limits();
    caps.max_paths = 5;
    Json::Value out = run(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"max_paths\": 9}", {},
                          true, caps);
    ASSERT_EQ(1u, out["limits"]["clamped"].size());
    EXPECT_EQ("max_paths", out["limits"]["clamped"][0]["field"].asString());
    EXPECT_EQ(9u, out["limits"]["clamped"][0]["requested"].asUInt64());
    EXPECT_EQ(5u, out["limits"]["clamped"][0]["effective"].asUInt64());
    EXPECT_FALSE(out["limits"].isMember("max_paths"));
    out = run(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"max_paths\": 9, \"long_search\": "
                   "\"paths\"}", {}, true, caps);
    EXPECT_EQ(5u, out["limits"]["max_paths"].asUInt64());
    // require_support in a request that reads no labels: stated as not read
    out = run(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"require_support\": "
                   "\"record_verified\", \"long_search\": \"paths\"}");
    const Json::Value &notes = out["patterns"][0]["notes"];
    ASSERT_GT(notes.size(), 0u);
    EXPECT_EQ("annotation_not_read", notes[notes.size() - 1].asString());
}

// Patterns with several paths: random records of three columns, IUPAC patterns cut from them
// and widened, every answer against both oracles
TEST(PatternPaths, RandomRecordsAgainstTheOracles) {
    std::mt19937 rng(4711);
    auto random_seq = [&](size_t length) {
        std::string s;
        for (size_t i = 0; i < length; ++i) {
            s += "ACGT"[rng() % 4];
        }
        return s;
    };
    std::vector<Record> records;
    // a segment of 14 bases held whole by two records of column c1 (r1, r4), and split across
    // the adjacent records r3 and r6 of column c0 (seq_ids 1 and 2): r3 ends with its first 9
    // bases (its k-mers 0-2 at k = 7), r6 starts with the rest from base 3 (its k-mers 3-7),
    // so that c0's column coordinates of its k-mers are consecutive across the two records —
    // carried by c0 (label_intersection), verified by c1 only
    const std::string shared = random_seq(14);
    for (int i = 0; i < 9; ++i) {
        std::string seq = random_seq(20 + rng() % 30);
        if (i == 1 || i == 4)
            seq = seq.substr(0, 8) + shared + seq.substr(8);
        if (i == 3)
            seq += shared.substr(0, 9);
        if (i == 6)
            seq = shared.substr(3) + seq;
        records.push_back({ std::string("c") + char('0' + i % 3), "r" + std::to_string(i), seq });
    }
    Index idx = build<annot::RowDiffColumnAnnotator>(7, records, true);
    std::vector<std::string> patterns = { shared, shared.substr(0, 10), shared.substr(2, 9) };
    for (const Record &r : records) {
        const size_t L = 8 + rng() % 7;
        if (r.seq.size() < L + 1)
            continue;
        std::string p = r.seq.substr(rng() % (r.seq.size() - L), L);
        patterns.push_back(p);
        // widened at two positions
        p[1] = "NRYS"[rng() % 4];
        p[L - 2] = 'N';
        patterns.push_back(p);
    }
    {
        // the shared segment: verified in c1, carried across two records by c0
        Json::Value out = run(idx, body("{\"dna\": \"" + shared + "\"}"));
        const Json::Value &e = out["patterns"][0];
        bool c0 = false;
        for (const Json::Value &r : e["results"]) {
            if (r["sequence"].asString() != shared)
                continue;
            for (const Json::Value &l : r["labels"]) {
                if (l["column"].asString() == "c0") {
                    c0 = true;
                    EXPECT_EQ("label_intersection", l["support"].asString()) << r;
                } else if (l["column"].asString() == "c1") {
                    EXPECT_EQ("record_verified", l["support"].asString()) << r;
                    EXPECT_EQ(2u, l["occurrences"]["value"].asUInt64()) << r;
                }
            }
        }
        EXPECT_TRUE(c0) << e;
    }
    for (const std::string &p : patterns) {
        const bool iupac = p.find_first_not_of("ACGT") != std::string::npos;
        const std::string spec = "{\"" + std::string(iupac ? "iupac" : "dna") + "\": \"" + p
                                 + "\"}";
        for (bool require : { false, true }) {
            Json::Value out = run(idx, body(spec, require ? "\"require_support\": "
                                                            "\"record_verified\"" : ""));
            const Json::Value &e = out["patterns"][0];
            SCOPED_TRACE(p + (require ? " record_verified" : ""));
            check_labelled_paths(idx, e, p, require);
        }
    }
}

TEST(PatternPaths, ReverseHitOnAWrappedPrimaryGraph) {
    // k = 3: the record CGTT; the PRIMARY graph stores one orientation of each k-mer, the
    // wrapper exposes both. AACG's forward path AAC, ACG lies on the reverse complements of
    // the stored GTT and CGT (anchored through rc(AAC), DESIGN §4.1), its reverse one is the
    // record itself
    const std::vector<Record> records = { { "x", "x0", "CGTT" }, { "y", "y0", "GTTA" } };
    Index idx = build<annot::ColumnCompressed<>>(3, records, false, DeBruijnGraph::PRIMARY);
    Json::Value out = run(idx, body("{\"dna\": \"AACG\"}", "\"allow_unbudgeted_annotation\": true"),
                          {}, false);
    EXPECT_EQ("primary", out["index"]["graph_mode"].asString());
    const Json::Value &e = out["patterns"][0];
    check_path_counts(idx, e, "AACG");
    const auto results = check_path_results(idx, e, "AACG");
    const PathSet expected {
        { "forward", "AACG" }, { "reverse", "CGTT" },
    };
    EXPECT_EQ(expected, PathSet(results.begin(), results.end()));
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    EXPECT_EQ("none_canonical", e["placement"].asString());
    std::set<std::string> notes;
    for (const Json::Value &n : e["notes"]) {
        notes.insert(n.asString());
    }
    EXPECT_TRUE(notes.count("strand_unknown_canonical"));
    EXPECT_TRUE(notes.count("label_intersection_only"));
    Scan scan(idx);
    for (Json::ArrayIndex i = 0; i < e["results"].size(); ++i) {
        const Json::Value &r = e["results"][i];
        std::set<std::string> got;
        for (const Json::Value &l : r["labels"]) {
            got.insert(l["column"].asString());
            EXPECT_EQ("label_intersection", l["support"].asString());
            EXPECT_FALSE(l.isMember("occurrence_list"));
        }
        // x holds CGT and GTT; y holds GTT only: x carries both paths, y neither
        EXPECT_EQ(scan.carriers(results[i].second), got) << r;
        EXPECT_EQ(std::set<std::string>({ "x" }), got);
        EXPECT_EQ("label_intersection", r["support"].asString());
    }
    // nothing can be verified there
    EXPECT_EQ("unknown", e["counts"]["labels"]["by_support"]["record_verified"]["relation"]
                             .asString());
    EXPECT_EQ(1u, e["counts"]["labels"]["by_support"]["label_intersection"]["value"].asUInt64());
    EXPECT_EQ("unknown", e["by_label"][0]["paths_record_verified"]["relation"].asString());
    EXPECT_EQ(std::make_pair(400, std::string("support_unavailable")),
              refusal(idx, body("{\"dna\": \"AACG\"}", "\"allow_unbudgeted_annotation\": true, "
                                "\"require_support\": \"record_verified\""), false));
}

TEST(PatternPaths, NativeCanonicalGraph) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, false,
                                                     DeBruijnGraph::CANONICAL);
    for (const std::string p : { "ACGTACC", "GGACGTACCA", "CGTACC" }) {
        Json::Value out = run(idx, body("{\"dna\": \"" + p + "\"}"), {}, false);
        const Json::Value &e = out["patterns"][0];
        check_path_counts(idx, e, p);
        const auto results = check_path_results(idx, e, p);
        EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
        EXPECT_EQ("none_canonical", e["placement"].asString());
        Scan scan(idx);
        for (Json::ArrayIndex i = 0; i < e["results"].size(); ++i) {
            std::set<std::string> got;
            for (const Json::Value &l : e["results"][i]["labels"]) {
                got.insert(l["column"].asString());
            }
            EXPECT_EQ(scan.carriers(results[i].second), got) << results[i].second;
        }
    }
}

// The two thresholds (§4.2, §5.2): max_anchors admits the extension, max_paths the release;
// stop_at_threshold, max_steps in the extension, and partial's cut, each with its relation
TEST(PatternPaths, ThresholdsAndStops) {
    std::mt19937 rng(17);
    std::vector<Record> records;
    for (int i = 0; i < 6; ++i) {
        std::string seq;
        for (int j = 0; j < 120; ++j) {
            seq += "ACGT"[rng() % 4];
        }
        records.push_back({ i % 2 ? "odd" : "even", "r" + std::to_string(i), seq });
    }
    Index idx = build<annot::RowDiffColumnAnnotator>(5, records, true);
    const std::string p = "ACNNNNNT";
    const std::string spec = "{\"iupac\": \"" + p + "\"}";
    Walk walk(idx);
    uint64_t anchors = 0, paths = 0;
    for (const auto &[strand, q] : oriented(p)) {
        anchors += walk.anchors(q).size();
        paths += walk.paths(q).size();
    }
    ASSERT_GT(anchors, 3u);
    ASSERT_GT(paths, 3u);

    // everything within the thresholds: complete
    Json::Value full = run(idx, body(spec, "", "none"));
    const Json::Value &all = full["patterns"][0];
    check_path_counts(idx, all, p);
    check_path_results(idx, all, p);
    EXPECT_TRUE(all["retrieval_complete"].asBool());

    // max_anchors below the anchors: the extension is not admitted, in every mode
    for (const std::string mode : { "all_or_count", "partial", "count" }) {
        Json::Value out = run(idx, body(spec, "\"mode\": \"" + mode + "\", \"max_anchors\": "
                                              + std::to_string(anchors - 1), "none"));
        const Json::Value &e = out["patterns"][0];
        EXPECT_EQ(anchors, e["counts"]["anchors"]["value"].asUInt64());
        EXPECT_EQ("exact", e["counts"]["anchors"]["relation"].asString());
        EXPECT_EQ("unknown", e["counts"]["paths"]["relation"].asString());
        EXPECT_EQ("not_admitted", e["counts"]["paths"]["extension"].asString());
        EXPECT_EQ(0u, e["counts"]["paths"]["candidates_examined"].asUInt64());
        EXPECT_EQ(0u, e["work"]["extension_edges"].asUInt64());
        EXPECT_TRUE(e["stop"].isNull());
        if (mode != "count") {
            EXPECT_EQ("anchors_above_threshold", e["withheld"]["reason"].asString()) << mode;
            EXPECT_EQ(0u, e["results"].size());
            EXPECT_FALSE(e["retrieval_complete"].asBool());
        }
    }

    // max_paths below the paths: all_or_count withholds with the exact count, partial
    // returns the first max_paths in answer order
    Json::Value out = run(idx, body(spec, "\"max_paths\": " + std::to_string(paths - 1), "none"));
    const Json::Value &over = out["patterns"][0];
    EXPECT_EQ("count_above_threshold", over["withheld"]["reason"].asString());
    EXPECT_EQ(paths, over["counts"]["paths"]["value"].asUInt64());
    EXPECT_EQ("exact", over["counts"]["paths"]["relation"].asString());
    EXPECT_EQ(paths - 1, out["limits"]["max_paths"].asUInt64());
    out = run(idx, body(spec, "\"max_paths\": 2, \"mode\": \"partial\"", "none"));
    const Json::Value &cut = out["patterns"][0];
    EXPECT_EQ("max_paths", cut["cut"]["reason"].asString());
    EXPECT_EQ(2u, cut["returned"].asUInt64());
    EXPECT_FALSE(cut["retrieval_complete"].asBool());
    EXPECT_EQ(paths, cut["counts"]["paths"]["value"].asUInt64());
    EXPECT_EQ("exact", cut["counts"]["paths"]["relation"].asString());
    for (Json::ArrayIndex i = 0; i < 2; ++i) {
        EXPECT_EQ(all["results"][i], cut["results"][i]);
    }

    // stop_at_threshold on the paths: the extension stops past max_paths, its count at_least
    out = run(idx, body(spec, "\"max_paths\": 2, \"stop_at_threshold\": true", "none"));
    const Json::Value &crossed = out["patterns"][0];
    EXPECT_EQ("extension", crossed["stop"]["phase"].asString());
    EXPECT_EQ("max_paths", crossed["stop"]["reason"].asString());
    EXPECT_EQ("at_least", crossed["counts"]["paths"]["relation"].asString());
    EXPECT_EQ(3u, crossed["counts"]["paths"]["value"].asUInt64());
    EXPECT_EQ("stopped", crossed["counts"]["paths"]["extension"].asString());
    EXPECT_EQ("exact", crossed["counts"]["anchors"]["relation"].asString());
    EXPECT_EQ("threshold_crossed", crossed["withheld"]["reason"].asString());
    out = run(idx, body(spec, "\"max_paths\": 2, \"stop_at_threshold\": true, \"mode\": "
                              "\"partial\"", "none"));
    EXPECT_EQ("max_paths", out["patterns"][0]["cut"]["reason"].asString());
    EXPECT_EQ(2u, out["patterns"][0]["returned"].asUInt64());

    // stop_at_threshold on the anchors: the anchors stop, nothing is extended
    out = run(idx, body(spec, "\"max_anchors\": 1, \"stop_at_threshold\": true, \"mode\": "
                              "\"partial\"", "none"));
    const Json::Value &early = out["patterns"][0];
    EXPECT_EQ("discovery", early["stop"]["phase"].asString());
    EXPECT_EQ("max_anchors", early["stop"]["reason"].asString());
    EXPECT_EQ("at_least", early["counts"]["anchors"]["relation"].asString());
    EXPECT_EQ("unknown", early["counts"]["paths"]["relation"].asString());
    EXPECT_EQ("not_started", early["counts"]["paths"]["extension"].asString());
    EXPECT_EQ("max_anchors", early["cut"]["reason"].asString());
    EXPECT_EQ(0u, early["returned"].asUInt64());

    // max_steps in the extension: the paths completed before the stop, at_least
    Json::Value count = run(idx, "{\"patterns\": [" + spec + "], \"mode\": \"count\"}");
    const uint64_t discovery = count["patterns"][0]["work"]["steps"].asUInt64();
    const uint64_t steps = all["work"]["steps"].asUInt64();
    ASSERT_GT(steps, discovery + 2);
    EXPECT_EQ(steps - discovery, all["work"]["extension_edges"].asUInt64());
    Json::Value stopped = run(idx, body(spec, "\"max_steps\": "
                                              + std::to_string(discovery + (steps - discovery) / 2),
                                        "none"));
    const Json::Value &budget = stopped["patterns"][0];
    EXPECT_EQ("extension", budget["stop"]["phase"].asString());
    EXPECT_EQ("max_steps", budget["stop"]["reason"].asString());
    EXPECT_EQ("at_least", budget["counts"]["paths"]["relation"].asString());
    EXPECT_LT(budget["counts"]["paths"]["value"].asUInt64(), paths);
    EXPECT_EQ("discovery_budget", budget["withheld"]["reason"].asString());
    out = run(idx, body(spec, "\"mode\": \"partial\", \"max_steps\": "
                              + std::to_string(discovery + (steps - discovery) / 2), "none"));
    const Json::Value &partial = out["patterns"][0];
    EXPECT_EQ("max_steps", partial["cut"]["reason"].asString());
    EXPECT_EQ(budget["counts"]["paths"]["value"], partial["returned"]);
    for (Json::ArrayIndex i = 0; i < partial["results"].size(); ++i) {
        // the paths completed before the stop are a prefix of the answer order
        EXPECT_EQ(all["results"][i], partial["results"][i]);
    }
    // a step stop is sticky: the next pattern answers unknown
    out = run(idx, body(spec + ", " + spec, "\"max_steps\": "
                        + std::to_string(discovery + (steps - discovery) / 2), "none"));
    EXPECT_EQ("unknown", out["patterns"][1]["counts"]["anchors"]["relation"].asString());
    EXPECT_EQ("unknown", out["patterns"][1]["counts"]["paths"]["relation"].asString());
    EXPECT_EQ("not_started", out["patterns"][1]["counts"]["paths"]["extension"].asString());

    // labelled: the same thresholds, the labels read only for an admitted release
    out = run(idx, body(spec, "\"max_paths\": " + std::to_string(paths - 1)));
    EXPECT_EQ("count_above_threshold", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ(0u, out["patterns"][0]["work"]["annotation_rows"].asUInt64());
    EXPECT_TRUE(out["patterns"][0]["by_label"].isNull());
    Json::Value labelled = run(idx, body(spec));
    check_labelled_paths(idx, labelled["patterns"][0], p, false);
}

// The memory account charges each released path (its descriptor, sequence and node and row
// arrays: path_descriptor_bytes) before its result object is built (§5.3)
TEST(PatternPaths, TheMemoryAccountStopsPathRetention) {
    std::mt19937 rng(17);
    std::vector<Record> records;
    for (int i = 0; i < 6; ++i) {
        std::string seq;
        for (int j = 0; j < 120; ++j) {
            seq += "ACGT"[rng() % 4];
        }
        records.push_back({ i % 2 ? "odd" : "even", "r" + std::to_string(i), seq });
    }
    Index idx = build<annot::RowDiffColumnAnnotator>(5, records, true);
    const std::string spec = "{\"iupac\": \"ACNNNNNT\"}";
    const uint64_t descriptor = path_descriptor_bytes(5, 8);
    EXPECT_EQ(512u + 2 * 5 + 3 * 8 + 192 * 4, descriptor);
    Json::Value full = run(idx, body(spec, "\"mode\": \"partial\""));
    const uint64_t paths = full["patterns"][0]["returned"].asUInt64();
    ASSERT_GT(paths, 3u);

    // room for the descriptors of six paths: partial's take at most half, three
    RetrievalHooks hooks;
    hooks.max_memory_bytes = 6 * descriptor;
    Json::Value out = run(idx, body(spec, "\"mode\": \"partial\""), hooks);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ(3u, e["returned"].asUInt64());
    EXPECT_EQ("max_memory", e["cut"]["reason"].asString());
    EXPECT_EQ("output", e["stop"]["phase"].asString());
    EXPECT_EQ("max_memory", e["stop"]["reason"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    EXPECT_EQ(paths, e["counts"]["paths"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
    EXPECT_LE(e["work"]["memory_bytes"].asUInt64(), hooks.max_memory_bytes);
    for (Json::ArrayIndex i = 0; i < 3; ++i) {
        EXPECT_EQ(full["patterns"][0]["results"][i]["sequence"], e["results"][i]["sequence"]);
    }
    // all_or_count: all or nothing
    out = run(idx, body(spec), hooks);
    EXPECT_EQ("output_budget", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ(0u, out["patterns"][0]["results"].size());

    // any account smaller than a complete run's peak: something is cut, refused or stopped,
    // said, and the account never passes its maximum
    Json::Value complete = run(idx, body(spec));
    ASSERT_TRUE(complete["patterns"][0]["retrieval_complete"].asBool());
    const uint64_t peak = complete["patterns"][0]["work"]["memory_bytes"].asUInt64();
    for (uint64_t max = descriptor; max < peak; max += (peak - descriptor) / 23 + 1) {
        hooks.max_memory_bytes = max;
        for (const std::string mode : { "all_or_count", "partial" }) {
            Json::Value x = run(idx, body(spec, "\"mode\": \"" + mode + "\""), hooks);
            const Json::Value &w = x["patterns"][0];
            EXPECT_FALSE(w["retrieval_complete"].asBool()) << max << " " << mode;
            EXPECT_LE(w["work"]["memory_bytes"].asUInt64(), max) << mode;
            EXPECT_TRUE(!w["withheld"].isNull() || !w["cut"].isNull() || !w["stop"].isNull()
                        || w["rows_refused"].size()) << max << " " << mode << " " << w;
        }
    }
}

// The reads and the output of the paths' labels are work under the deadline (§5.3), as for
// contexts: a virtual clock moved past the work time at the first read, or when the second
// path's labels are to be built
TEST(PatternPaths, DeadlineInTheReadsAndInTheOutput) {
    std::mt19937 rng(17);
    std::vector<Record> records;
    for (int i = 0; i < 6; ++i) {
        std::string seq;
        for (int j = 0; j < 120; ++j) {
            seq += "ACGT"[rng() % 4];
        }
        records.push_back({ i % 2 ? "odd" : "even", "r" + std::to_string(i), seq });
    }
    Index idx = build<annot::RowDiffColumnAnnotator>(5, records, true);
    const std::string spec = "{\"iupac\": \"ACNNNNNT\"}";
    const Clock::time_point start = Clock::now();
    auto virtual_ms = std::make_shared<double>(0);
    auto clock = [start, virtual_ms]() {
        return start + std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double, std::milli>(*virtual_ms));
    };

    // the reads: every read moves the clock far past the budget
    RetrievalHooks reads;
    reads.read_hook = [virtual_ms](size_t) { *virtual_ms += 1e6; };
    Json::Value out = run_clocked(idx, body(spec + ", " + spec, "\"mode\": \"partial\""),
                                  reads, clock);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("label_discovery", e["stop"]["phase"].asString());
    EXPECT_EQ("time", e["stop"]["reason"].asString());
    EXPECT_EQ("time_limited", e["determinism"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    // the extension completed before the clock moved
    EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
    EXPECT_GT(e["results"].size(), 1u);
    for (const Json::Value &r : e["results"]) {
        EXPECT_NE("complete", r["labels_status"].asString());
    }
    // the next pattern: after a time stop, as after any
    const Json::Value &next = out["patterns"][1];
    EXPECT_EQ("discovery", next["stop"]["phase"].asString());
    EXPECT_EQ("time", next["stop"]["reason"].asString());
    EXPECT_EQ("unknown", next["counts"]["anchors"]["relation"].asString());
    *virtual_ms = 0;
    out = run_clocked(idx, body(spec), reads, clock);
    EXPECT_EQ("deadline", out["patterns"][0]["withheld"]["reason"].asString());

    // the output: past the work time when the second path's labels are to be built
    RetrievalHooks output;
    output.output_hook = [virtual_ms](size_t path) {
        if (path >= 1)
            *virtual_ms = 1e9;
    };
    *virtual_ms = 0;
    out = run_clocked(idx, body(spec, "\"mode\": \"partial\""), output, clock);
    const Json::Value &o = out["patterns"][0];
    EXPECT_EQ("output", o["stop"]["phase"].asString());
    EXPECT_EQ("time", o["stop"]["reason"].asString());
    EXPECT_EQ("time_limited", o["determinism"].asString());
    EXPECT_FALSE(o["retrieval_complete"].asBool());
    const Json::Value &results = o["results"];
    ASSERT_GT(results.size(), 1u);
    EXPECT_EQ("complete", results[0]["labels_status"].asString());
    EXPECT_TRUE(results[0]["labels"].isArray());
    for (Json::ArrayIndex i = 1; i < results.size(); ++i) {
        EXPECT_EQ("output_budget", results[i]["labels_status"].asString()) << i;
        EXPECT_TRUE(results[i]["labels"].isNull()) << i;
        EXPECT_TRUE(results[i]["support"].isNull()) << i;
    }
    *virtual_ms = 0;
    out = run_clocked(idx, body(spec), output, clock);
    EXPECT_EQ("deadline", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ(0u, out["patterns"][0]["results"].size());
    // the clock left alone: complete
    *virtual_ms = 0;
    out = run_clocked(idx, body(spec), {}, clock);
    EXPECT_TRUE(out["patterns"][0]["retrieval_complete"].asBool());
}

TEST(PatternPaths, TruncatedAndRefusedRows) {
    // GGACGTACCAA twice, in columns A and D: every row of its paths carries two labels
    std::vector<Record> records = kRecords;
    records.push_back({ "D", "d0", "GGACGTACCAA" });
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, records, true);
    Json::Value out = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"max_labels_per_anchor\": 1"));
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("anchor_labels_truncated", e["withheld"]["reason"].asString());
    EXPECT_EQ(0u, e["results"].size());
    // the six rows of the path, each stated once with its k-mer and its total (A and D; B
    // too for ACGTA, CGTAC and GTACC)
    Scan scan(idx);
    std::set<std::string> truncated;
    for (const Json::Value &t : e["anchors_truncated"]) {
        const std::string kmer = t["kmer"].asString();
        uint64_t total = 0;
        for (const auto &[column, recs] : scan.columns) {
            total += scan.holds_kmer(column, kmer);
        }
        EXPECT_EQ(total, t["total"].asUInt64()) << kmer;
        EXPECT_EQ(1u, t["cap"].asUInt64());
        truncated.insert(kmer);
    }
    const std::set<std::string> expected { "GGACG", "GACGT", "ACGTA", "CGTAC", "GTACC",
                                           "TACCA" };
    EXPECT_EQ(expected, truncated);
    // partial: the path's labels a true but partial list, their total unknown
    out = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"max_labels_per_anchor\": 1, \"mode\": "
                        "\"partial\""));
    const Json::Value &p = out["patterns"][0];
    ASSERT_EQ(1u, p["results"].size());
    EXPECT_EQ("truncated", p["results"][0]["labels_status"].asString());
    EXPECT_TRUE(p["results"][0]["labels_total"].isNull());
    EXPECT_LE(p["results"][0]["labels"].size(), 1u);
    EXPECT_EQ("at_least", p["counts"]["labels"]["relation"].asString());
    EXPECT_FALSE(p["retrieval_complete"].asBool());

    // every read refused: the rows stated, no label claimed
    RetrievalHooks hooks;
    hooks.deny_decode = [](uint64_t) { return true; };
    out = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &r = out["patterns"][0];
    ASSERT_EQ(1u, r["results"].size());
    EXPECT_EQ("refused", r["results"][0]["labels_status"].asString());
    EXPECT_TRUE(r["results"][0]["labels"].isNull());
    EXPECT_TRUE(r["results"][0]["support"].isNull());
    EXPECT_EQ(6u, r["rows_refused"].size());
    EXPECT_EQ("at_least", r["counts"]["labels"]["relation"].asString());
    out = run(idx, body("{\"dna\": \"GGACGTACCA\"}"), hooks);
    EXPECT_EQ("annotation_budget", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ(1u, out["patterns"][0]["rows_refused"].size());

    // the placement reads refused (after the discovery's): every label unverified, stated
    Json::Value plain = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"mode\": \"partial\""));
    // six rows, each read in both steps
    const uint64_t rows = plain["patterns"][0]["work"]["annotation_rows"].asUInt64() / 2;
    ASSERT_EQ(6u, rows);
    // the first charge of every read has the ordinal 0 (a DecodeBudget per read)
    auto reads = std::make_shared<uint64_t>(0);
    RetrievalHooks place;
    place.deny_decode = [reads, rows](uint64_t ordinal) {
        if (!ordinal)
            ++*reads;
        return *reads > rows;
    };
    out = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"mode\": \"partial\""), place);
    const Json::Value &v = out["patterns"][0];
    EXPECT_FALSE(v["retrieval_complete"].asBool());
    EXPECT_GT(v["rows_refused"].size(), 0u);
    for (const Json::Value &x : v["rows_refused"]) {
        EXPECT_EQ("placement", x["phase"].asString());
    }
    ASSERT_EQ(1u, v["results"].size());
    for (const Json::Value &l : v["results"][0]["labels"]) {
        EXPECT_EQ("label_intersection", l["support"].asString());
        EXPECT_EQ("unknown", l["occurrences"]["relation"].asString());
        EXPECT_TRUE(l["occurrence_list"].isNull());
    }
    EXPECT_EQ("at_least", v["counts"]["labels"]["by_support"]["record_verified"]["relation"]
                              .asString());
    EXPECT_EQ("unknown", v["counts"]["labels"]["by_support"]["label_intersection"]["relation"]
                             .asString());
    // with require_support the unverified labels are left out, their count unknown
    *reads = 0;
    out = run(idx, body("{\"dna\": \"GGACGTACCA\"}", "\"mode\": \"partial\", "
                        "\"require_support\": \"record_verified\""), place);
    const Json::Value &q = out["patterns"][0];
    EXPECT_EQ(0u, q["results"][0]["labels"].size());
    // neither label was verified nor refuted: excluded, their number not stated as known
    // (review of increments 4 and 5: it was 2, as if both had been refuted)
    EXPECT_EQ("complete", q["results"][0]["labels_status"].asString());
    EXPECT_EQ(2u, q["results"][0]["labels_total"].asUInt64());
    EXPECT_TRUE(q["results"][0]["labels_excluded_unverified"].isNull()) << q["results"][0];
    EXPECT_EQ("unknown", q["labels_excluded_unverified"]["relation"].asString());
    EXPECT_EQ("at_least", q["counts"]["labels"]["relation"].asString());
}

// A path's labels_excluded_unverified is an integer only when it is the true number: every row
// of the path read completely and every label carrying it verified or refuted; null otherwise,
// as labels_total (review of increments 4 and 5, finding 1: with max_labels_per_anchor below a
// row's label count the path's truncated rows gave a definite integer over the labels kept,
// 0 where B, carrying the path unverified, had been cut). The expectations: the record scan.
TEST(PatternPaths, ExcludedUnverifiedOnlyWhenEveryLabelIsDecided) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    const std::string p = "ACGTACC";
    Scan scan(idx);
    const std::set<std::string> carriers = scan.carriers(p);
    std::set<std::string> verified;
    for (const std::string &c : carriers) {
        if (!scan.occurrences(c, p).empty())
            verified.insert(c);
    }
    // A verifies it, B carries it across b0 -> b1: every row of the path has two labels
    ASSERT_EQ((std::set<std::string> { "A", "B" }), carriers);
    ASSERT_EQ(std::set<std::string> { "A" }, verified);
    const std::string require = "\"require_support\": \"record_verified\", \"mode\": \"partial\"";

    // the rows truncated (one label of two kept): the count unknown, for the path and the entry
    for (const std::string &rest : { require, std::string("\"mode\": \"partial\"") }) {
        const bool verified_only = rest == require;
        Json::Value out = run(idx, body("{\"dna\": \"" + p + "\"}",
                                        rest + ", \"max_labels_per_anchor\": 1"));
        const Json::Value &e = out["patterns"][0];
        ASSERT_EQ(1u, e["results"].size()) << e;
        const Json::Value &r = e["results"][0];
        EXPECT_EQ("truncated", r["labels_status"].asString()) << r;
        EXPECT_TRUE(r["labels_total"].isNull()) << r;
        EXPECT_EQ(3u, e["anchors_truncated"].size());
        EXPECT_FALSE(e["retrieval_complete"].asBool());
        if (verified_only) {
            EXPECT_TRUE(r["labels_excluded_unverified"].isNull()) << r;
            EXPECT_EQ("unknown", e["labels_excluded_unverified"]["relation"].asString());
            for (const Json::Value &l : r["labels"]) {
                EXPECT_EQ("record_verified", l["support"].asString());
            }
        } else {
            EXPECT_FALSE(r.isMember("labels_excluded_unverified")) << r;
            EXPECT_FALSE(e.isMember("labels_excluded_unverified"));
        }
    }
    // the cap at the rows' label count: read completely, every label decided, the count exact
    for (const char *cap : { "2", "3" }) {
        Json::Value out = run(idx, body("{\"dna\": \"" + p + "\"}",
                                        require + ", \"max_labels_per_anchor\": "
                                        + std::string(cap)));
        const Json::Value &e = out["patterns"][0];
        ASSERT_EQ(1u, e["results"].size()) << e;
        const Json::Value &r = e["results"][0];
        EXPECT_EQ("complete", r["labels_status"].asString()) << r;
        EXPECT_EQ(carriers.size(), r["labels_total"].asUInt64()) << r;
        EXPECT_EQ(carriers.size() - verified.size(),
                  r["labels_excluded_unverified"].asUInt64()) << r;
        EXPECT_EQ(carriers.size() - verified.size(),
                  e["labels_excluded_unverified"]["value"].asUInt64());
        EXPECT_EQ("exact", e["labels_excluded_unverified"]["relation"].asString());
        check_labelled_paths(idx, e, p, true);
    }

    // the placement stopped by the clock after its first read: the labels neither verified nor
    // refuted, the count unknown (stop {placement, time})
    Json::Value plain = run(idx, body("{\"dna\": \"" + p + "\"}", require));
    const uint64_t rows = plain["patterns"][0]["work"]["annotation_rows"].asUInt64() / 2;
    ASSERT_EQ(3u, rows);
    const Clock::time_point start = Clock::now();
    auto virtual_ms = std::make_shared<double>(0);
    auto clock = [start, virtual_ms]() {
        return start + std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double, std::milli>(*virtual_ms));
    };
    auto reads = std::make_shared<uint64_t>(0);
    RetrievalHooks late;
    late.read_hook = [reads, rows, virtual_ms](size_t) {
        if (++*reads > rows)
            *virtual_ms = 1e9;
    };
    Json::Value out = run_clocked(idx, body("{\"dna\": \"" + p + "\"}", require), late, clock);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("placement", e["stop"]["phase"].asString()) << e["stop"];
    EXPECT_EQ("time", e["stop"]["reason"].asString());
    ASSERT_EQ(1u, e["results"].size());
    const Json::Value &r = e["results"][0];
    EXPECT_EQ(carriers.size(), r["labels_total"].asUInt64()) << r;
    EXPECT_TRUE(r["labels_excluded_unverified"].isNull()) << r;
    EXPECT_EQ("unknown", e["labels_excluded_unverified"]["relation"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
}

TEST(PatternPaths, GlobalPlacementListsTheChainsUnverified) {
    // coordinates without the record mapping: B's chain across b0 and b1 cannot be told from
    // a record's, nothing is verified (record_bounds_unknown)
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"ACGTACC\"}"), {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("global", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    std::set<std::string> notes;
    for (const Json::Value &n : e["notes"]) {
        notes.insert(n.asString());
    }
    EXPECT_TRUE(notes.count("record_bounds_unknown"));
    EXPECT_FALSE(notes.count("label_intersection_only"));
    ASSERT_EQ(1u, e["results"].size());
    std::map<std::string, std::set<uint64_t>> chains;
    for (const Json::Value &l : e["results"][0]["labels"]) {
        EXPECT_EQ("label_intersection", l["support"].asString());
        EXPECT_FALSE(l.isMember("occurrences"));
        for (const Json::Value &o : l["occurrence_list"]) {
            EXPECT_EQ(0u, o["offset"].asUInt64());
            EXPECT_EQ("+", o["strand"].asString());
            chains[l["column"].asString()].insert(o["kmer_coord"].asUInt64());
        }
    }
    // A: a0's ACGTA at its column coordinate 2; B: b0's ACGTA at 3 (then b1's at 4, 5)
    const std::map<std::string, std::set<uint64_t>> expected { { "A", { 2 } }, { "B", { 3 } } };
    EXPECT_EQ(expected, chains);
    EXPECT_EQ("unknown", e["counts"]["occurrences"]["relation"].asString());
    EXPECT_EQ("unknown", e["counts"]["labels"]["by_support"]["record_verified"]["relation"]
                             .asString());
    EXPECT_EQ(2u, e["counts"]["labels"]["by_support"]["label_intersection"]["value"].asUInt64());
    EXPECT_EQ(std::make_pair(400, std::string("support_unavailable")),
              refusal(idx, body("{\"dna\": \"ACGTACC\"}", "\"require_support\": "
                                "\"record_verified\""), false));
}

TEST(PatternPaths, WithoutCoordinatesTheIntersectionOnly) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, false);
    Json::Value out = run(idx, body("{\"dna\": \"ACGTACC\"}"), {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("none", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    std::set<std::string> notes;
    for (const Json::Value &n : e["notes"]) {
        notes.insert(n.asString());
    }
    EXPECT_TRUE(notes.count("label_intersection_only"));
    ASSERT_EQ(1u, e["results"].size());
    std::set<std::string> got;
    for (const Json::Value &l : e["results"][0]["labels"]) {
        got.insert(l["column"].asString());
        EXPECT_EQ("label_intersection", l["support"].asString());
        EXPECT_FALSE(l.isMember("occurrences"));
        EXPECT_FALSE(l.isMember("occurrence_list"));
    }
    EXPECT_EQ(std::set<std::string>({ "A", "B" }), got);
    EXPECT_EQ(std::make_pair(400, std::string("support_unavailable")),
              refusal(idx, body("{\"dna\": \"ACGTACC\"}", "\"require_support\": "
                                "\"record_verified\""), false));
}

TEST(PatternPaths, OccurrencesNotRequested) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"long_search\": "
                               "\"paths\", \"output\": {\"labels\": \"all\", \"occurrences\": "
                               "false}}");
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("not_requested", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    bool note = false;
    for (const Json::Value &n : e["notes"]) {
        note |= n.asString() == "label_intersection_only";
    }
    EXPECT_TRUE(note);
    for (const Json::Value &l : e["results"][0]["labels"]) {
        EXPECT_EQ("label_intersection", l["support"].asString());
    }
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"long_search\": \"paths\", "
                           "\"require_support\": \"record_verified\", \"output\": {\"labels\": "
                           "\"all\", \"occurrences\": false}}"));
}

TEST(PatternPaths, CountMode) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, "{\"patterns\": [{\"dna\": \"ACGTACC\"}, {\"dna\": \"ACCAAG\"}], "
                               "\"long_search\": \"paths\", \"mode\": \"count\"}");
    EXPECT_TRUE(out["output"].isNull());
    for (int i : { 0, 1 }) {
        const Json::Value &e = out["patterns"][i];
        check_path_counts(idx, e, i ? "ACCAAG" : "ACGTACC");
        EXPECT_FALSE(e.isMember("results"));
        EXPECT_FALSE(e["retrieval_complete"].asBool());
    }
    // the CLI answers alike
    std::ostringstream text;
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    RetrievalHooks hooks;
    hooks.coord_to_header = idx.cth.get();
    ASSERT_TRUE(write_pattern_answer("{\"patterns\": [{\"dna\": \"ACGTACC\"}], \"long_search\": "
                                     "\"paths\", \"output\": {\"labels\": \"all\"}}",
                                     *idx.anno, limits(), "rel", nullptr, builder, text,
                                     "the request", &hooks));
    Json::Value cli;
    std::istringstream in(text.str());
    in >> cli;
    Json::Value server = run(idx, body("{\"dna\": \"ACGTACC\"}"));
    // (as text: a parsed number is a signed JSON value, the route's unsigned)
    EXPECT_EQ(Json::writeString(builder, untimed(server)), Json::writeString(builder, untimed(cli)));
}


// ------------------------------------------------------------------ repeats (review GPT-3)
//
// GPT review 3, finding 1: the verification of a path binary-searched every coordinate of its
// first k-mer in every other k-mer's list, made it all again for the output, and read no clock
// (a 30,000-base homopolymer and a 1,500-base path: 1.5 s past a budget of 500 ms, 503). The
// chains of a path are now the intersection of its k-mers' coordinate lists shifted (a leapfrog
// join, consecutive chains kept as runs), made once, under the clock. Records of repeats, where
// a path's chains form long runs: a homopolymer in two records of one column (its chains run
// on across the two records' coordinates, which the record bounds must cut), a dinucleotide
// repeat (no two chains consecutive), and runs of A broken by single C's.

std::string repeated(const std::string &unit, size_t times) {
    std::string s;
    for (size_t i = 0; i < times; ++i) {
        s += unit;
    }
    return s;
}

const std::vector<Record> kRepeats = {
    { "h", "h0", std::string(300, 'A') },
    { "h", "h1", std::string(120, 'A') },
    { "d", "d0", repeated("AC", 150) },
    { "m", "m0", std::string(40, 'A') + "C" + std::string(25, 'A') + "C" + std::string(60, 'A') },
};

TEST(PatternPaths, RepeatsAgainstTheOracles) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRepeats, true);
    const std::vector<std::string> patterns = {
        std::string(50, 'A'), std::string(25, 'A'), repeated("AC", 20), repeated("CA", 12),
        "C" + std::string(25, 'A') + "C", std::string(10, 'A') + "W" + std::string(9, 'A'),
        "N" + std::string(18, 'A') + "N",
    };
    for (const std::string &p : patterns) {
        const bool iupac = p.find_first_not_of("ACGT") != std::string::npos;
        const std::string spec = "{\"" + std::string(iupac ? "iupac" : "dna") + "\": \"" + p
                                 + "\"}";
        for (bool require : { false, true }) {
            SCOPED_TRACE(p + (require ? " record_verified" : ""));
            Json::Value out = run(idx, body(spec, require ? "\"require_support\": "
                                                            "\"record_verified\"" : ""));
            check_labelled_paths(idx, out["patterns"][0], p, require);
        }
    }
    // A^50: 251 starts in h0 and 71 in h1, none across the two (their chains are one run of
    // consecutive column coordinates); 11 in m0's run of 60
    Json::Value out = run(idx, body("{\"dna\": \"" + std::string(50, 'A') + "\"}"));
    const Json::Value &e = out["patterns"][0];
    ASSERT_EQ(1u, e["results"].size());
    std::map<std::string, uint64_t> occurrences;
    for (const Json::Value &l : e["results"][0]["labels"]) {
        occurrences[l["column"].asString()] = l["occurrences"]["value"].asUInt64();
    }
    EXPECT_EQ((std::map<std::string, uint64_t> { { "h", 251 + 71 }, { "m", 11 } }), occurrences);
}

// partial: each label lists the first max_occurrences_per_label occurrences of its union (over
// the paths returned), its counts exact over all of them (the record scan), and the cut stated
TEST(PatternPaths, TheCapListsTheFirstOccurrencesOfEachUnion) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRepeats, true);
    Scan scan(idx);
    const uint64_t cap = 3;
    for (const std::string &p : { std::string(25, 'A'), "N" + std::string(18, 'A') + "N",
                                 repeated("AC", 20) }) {
        SCOPED_TRACE(p);
        const bool iupac = p.find_first_not_of("ACGT") != std::string::npos;
        Json::Value out = run(idx, body("{\"" + std::string(iupac ? "iupac" : "dna") + "\": \""
                                        + p + "\"}", "\"mode\": \"partial\", "
                                        "\"max_occurrences_per_label\": " + std::to_string(cap)));
        const Json::Value &e = out["patterns"][0];
        ASSERT_TRUE(e["stop"].isNull()) << e["stop"];
        // per column: every occurrence of every path returned (the scan), in the union's order
        std::map<std::string, std::set<std::tuple<uint64_t, uint64_t, uint8_t>>> unions;
        std::map<std::string, std::set<std::tuple<uint64_t, uint64_t, uint8_t>>> listed;
        for (const Json::Value &r : e["results"]) {
            const std::string s = r["sequence"].asString();
            const uint8_t strand = orientation_rank(r["strand"].asString());
            for (const Json::Value &l : r["labels"]) {
                const std::string column = l["column"].asString();
                const auto occ = scan.occurrences(column, s);
                for (const auto &[seq_id, start] : occ) {
                    unions[column].emplace(seq_id, start, strand);
                }
                // the path's own count, exact, whatever is listed
                EXPECT_EQ(occ.size(), l["occurrences"]["value"].asUInt64()) << column;
                EXPECT_EQ("exact", l["occurrences"]["relation"].asString());
                for (const Json::Value &o : l["occurrence_list"]) {
                    const std::string coords = o["nt_coords"].asString();
                    const uint64_t start = std::stoull(coords.substr(0, coords.find('-')));
                    EXPECT_TRUE(occ.count({ o["seq_id"].asUInt64(), start })) << o;
                    listed[column].emplace(o["seq_id"].asUInt64(), start, strand);
                }
            }
        }
        ASSERT_FALSE(unions.empty());
        uint64_t cut = 0, total = 0;
        for (const auto &[column, u] : unions) {
            std::set<std::tuple<uint64_t, uint64_t, uint8_t>> first;
            for (auto it = u.begin(); it != u.end() && first.size() < cap; ++it) {
                first.insert(*it);
            }
            EXPECT_EQ(first, listed[column]) << column;
            cut += u.size() > cap;
            total += u.size();
        }
        for (const Json::Value &b : e["by_label"]) {
            EXPECT_EQ(unions[b["column"].asString()].size(), b["occurrences"]["value"].asUInt64());
            EXPECT_EQ("exact", b["occurrences"]["relation"].asString());
        }
        EXPECT_EQ(total, e["counts"]["occurrences"]["value"].asUInt64());
        EXPECT_EQ("exact", e["counts"]["occurrences"]["relation"].asString());
        if (cut) {
            EXPECT_EQ("max_occurrences_per_label", e["occurrences_cut"]["reason"].asString());
            EXPECT_EQ(cut, e["occurrences_cut"]["labels"].asUInt64());
        } else {
            EXPECT_TRUE(e["occurrences_cut"].isNull());
        }
    }
}

// The work of the verification, seen through the clock's readings (deterministic: every
// reading after kClockStride units of work): a 1,000-base path through a 3,000-base
// homopolymer has 2,001 chains, one run, found with a few units per k-mer — not a search per
// chain and k-mer (two million units, some 500 readings). The readings between the path's
// verification and the output of its labels: those of a few thousand units at most. (A
// column-compressed index, read unbudgeted: a row-diff conversion of the homopolymer's row
// takes seconds)
TEST(PatternPaths, AHomopolymersChainsAreOneRun) {
    Index idx = build<annot::ColumnCompressed<>>(kK, { { "h", "h0", std::string(3000, 'A') } },
                                                 true);
    const Clock::time_point start = Clock::now();
    auto armed = std::make_shared<bool>(false);
    auto readings = std::make_shared<uint64_t>(0);
    auto clock = [start, armed, readings]() {
        *readings += *armed;
        return start;
    };
    RetrievalHooks hooks;
    hooks.occurrences_hook = [armed](size_t) { *armed = true; };
    hooks.output_hook = [armed](size_t) { *armed = false; };
    Json::Value out = run_clocked(idx, body("{\"dna\": \"" + std::string(1000, 'A') + "\"}",
                                            "\"mode\": \"partial\", \"max_occurrences_per_label\": "
                                            "2, \"allow_unbudgeted_annotation\": true"), hooks,
                                  clock);
    const Json::Value &e = out["patterns"][0];
    ASSERT_EQ(1u, e["results"].size());
    const Json::Value &l = e["results"][0]["labels"][0];
    EXPECT_EQ("record_verified", l["support"].asString());
    EXPECT_EQ(2001u, l["occurrences"]["value"].asUInt64());
    ASSERT_EQ(2u, l["occurrence_list"].size());
    EXPECT_EQ("1-1000", l["occurrence_list"][0]["nt_coords"].asString());
    EXPECT_EQ("2-1001", l["occurrence_list"][1]["nt_coords"].asString());
    EXPECT_LE(*readings, 3u);
}

// The verification reads the clock inside a path (review GPT-3, finding 1): a dinucleotide
// repeat, whose 1,901 chains of a 200-base path are no run, each sought in 195 lists; the clock
// passes the work time once the verification began (after the path's own reading): stop
// {placement, time}, its labels neither verified nor refuted, stated — not finished first and
// stopped in the output after it (or 503 after the budget)
TEST(PatternPaths, TheVerificationReadsTheClock) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, { { "d", "d0", repeated("AC", 2000) } },
                                                     true);
    const std::string p = repeated("AC", 100);
    const Clock::time_point start = Clock::now();
    auto virtual_ms = std::make_shared<double>(0);
    auto clock = [start, virtual_ms]() {
        return start + std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double, std::milli>(*virtual_ms));
    };
    RetrievalHooks hooks;
    hooks.occurrences_hook = [virtual_ms](size_t) { *virtual_ms = 1e9; };
    for (const std::string mode : { "partial", "all_or_count" }) {
        *virtual_ms = 0;
        Json::Value out = run_clocked(idx, body("{\"dna\": \"" + p + "\"}",
                                                "\"mode\": \"" + mode + "\""), hooks, clock);
        const Json::Value &e = out["patterns"][0];
        EXPECT_EQ("placement", e["stop"]["phase"].asString()) << mode << " " << e["stop"];
        EXPECT_EQ("time", e["stop"]["reason"].asString());
        EXPECT_EQ("time_limited", e["determinism"].asString());
        EXPECT_FALSE(e["retrieval_complete"].asBool());
        EXPECT_EQ("exact", e["counts"]["paths"]["relation"].asString());
        if (mode == "partial") {
            ASSERT_EQ(1u, e["results"].size());
            EXPECT_TRUE(e["results"][0]["labels"].isNull());
            EXPECT_EQ("unknown", e["counts"]["occurrences"]["relation"].asString());
        } else {
            EXPECT_EQ("deadline", e["withheld"]["reason"].asString());
        }
    }
    // the clock left alone: every chain verified (the record scan's 1,901)
    *virtual_ms = 0;
    Json::Value out = run_clocked(idx, body("{\"dna\": \"" + p + "\"}"), {}, clock);
    check_labelled_paths(idx, out["patterns"][0], p, false);
    EXPECT_EQ(1901u, out["patterns"][0]["counts"]["occurrences"]["value"].asUInt64());
}

// Over a sweep of memory maxima on the repeats (the verification's runs, the occurrences'
// deduplication and the labels built all in the account), both modes: the account never passes
// its maximum, an incomplete entry says why, a path whose labels are listed lists them all
// (labels_total of them, unless max_labels cut them) with its true occurrences (the record
// scan's), and one whose labels were not built says so (output_budget, labels null) — not an
// empty list beside labels_total
TEST(PatternPaths, MemorySweepOverRepeats) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRepeats, true);
    Scan scan(idx);
    const std::string spec = "{\"iupac\": \"N" + std::string(18, 'A') + "N\"}, {\"dna\": \""
                             + repeated("AC", 20) + "\"}";
    for (const std::string mode : { "partial", "all_or_count" }) {
        Json::Value full = run(idx, body(spec, "\"mode\": \"" + mode + "\""));
        ASSERT_TRUE(full["patterns"][0]["retrieval_complete"].asBool() || mode == "partial");
        ASSERT_GT(full["patterns"][0]["results"].size(), 2u);
        const uint64_t peak = full["patterns"][1]["work"]["memory_bytes"].asUInt64();
        uint64_t stated = 0;
        for (uint64_t max = 2000; max <= peak; max += std::max<uint64_t>(1, (peak - 2000) / 100)) {
            SCOPED_TRACE(mode + " max " + std::to_string(max));
            RetrievalHooks hooks;
            hooks.max_memory_bytes = max;
            Json::Value out = run(idx, body(spec, "\"mode\": \"" + mode + "\""), hooks);
            for (const Json::Value &e : out["patterns"]) {
                ASSERT_LE(e["work"]["memory_bytes"].asUInt64(), max);
                if (!e["retrieval_complete"].asBool() && e["occurrences_cut"].isNull()) {
                    stated++;
                    EXPECT_TRUE(!e["withheld"].isNull() || !e["cut"].isNull()
                                || !e["stop"].isNull() || e["rows_refused"].size()) << e;
                }
                for (const Json::Value &r : e["results"]) {
                    const std::string status = r["labels_status"].asString();
                    if (status != "complete") {
                        if (status == "output_budget") {
                            EXPECT_TRUE(r["labels"].isNull()) << r;
                        }
                        continue;
                    }
                    ASSERT_TRUE(r["labels"].isArray()) << r;
                    if (e["labels_cut"].isNull()) {
                        EXPECT_EQ(r["labels_total"].asUInt64(), r["labels"].size()) << r;
                    }
                    const std::string s = r["sequence"].asString();
                    for (const Json::Value &l : r["labels"]) {
                        if (l["occurrences"]["relation"].asString() != "exact")
                            continue;
                        const auto occ = scan.occurrences(l["column"].asString(), s);
                        EXPECT_EQ(occ.size(), l["occurrences"]["value"].asUInt64()) << l;
                        EXPECT_LE(l["occurrence_list"].size(), occ.size());
                    }
                }
            }
        }
        EXPECT_GT(stated, 0u);
    }
}


// As MemorySweepOverRepeats over the clock: the clock passes the work time at its R-th reading,
// for every R of the request (the readings are a deterministic function of the request: every
// stride of work, every path, every row read), so that every time stop the request can meet —
// in the extension, the reads, the paths' label lists, their verification, the output — is
// met once. Each answer states it, a path returned with labels_status complete lists all its
// labels (not an empty list after a stop of the lists), and none is refused
TEST(PatternPaths, TimeSweepOverRepeats) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRepeats, true);
    Scan scan(idx);
    const std::string spec = "{\"iupac\": \"N" + std::string(18, 'A') + "N\"}";
    const Clock::time_point start = Clock::now();
    auto readings = std::make_shared<uint64_t>(0);
    auto jump = std::make_shared<uint64_t>(0);
    auto clock = [start, readings, jump]() {
        return ++*readings > *jump ? start + std::chrono::hours(1) : start;
    };
    for (const std::string mode : { "partial", "all_or_count" }) {
        const std::string b = body(spec, "\"mode\": \"" + mode + "\"");
        *readings = 0;
        *jump = std::numeric_limits<uint64_t>::max();
        Json::Value full = run_clocked(idx, b, {}, clock);
        ASSERT_GT(full["patterns"][0]["results"].size(), 2u);
        const uint64_t total = *readings;
        std::set<std::string> stops;
        for (uint64_t r = 0; r <= total; ++r) {
            SCOPED_TRACE(mode + " reading " + std::to_string(r));
            *readings = 0;
            *jump = r;
            Json::Value out = run_clocked(idx, b, {}, clock);
            const Json::Value &e = out["patterns"][0];
            if (!e["stop"].isNull())
                stops.insert(e["stop"]["phase"].asString() + "/" + e["stop"]["reason"].asString());
            if (!e["retrieval_complete"].asBool() && e["occurrences_cut"].isNull()) {
                EXPECT_TRUE(!e["withheld"].isNull() || !e["cut"].isNull() || !e["stop"].isNull())
                        << e;
            }
            for (const Json::Value &res : e["results"]) {
                const std::string status = res["labels_status"].asString();
                if (status == "output_budget") {
                    EXPECT_TRUE(res["labels"].isNull()) << res;
                }
                if (status != "complete")
                    continue;
                ASSERT_TRUE(res["labels"].isArray()) << res;
                EXPECT_EQ(res["labels_total"].asUInt64(), res["labels"].size()) << res;
                for (const Json::Value &l : res["labels"]) {
                    if (l["occurrences"]["relation"].asString() == "exact") {
                        EXPECT_EQ(scan.occurrences(l["column"].asString(),
                                                   res["sequence"].asString()).size(),
                                  l["occurrences"]["value"].asUInt64()) << l;
                    }
                }
            }
        }
        // the stops met on the way include the output's and the placement's
        EXPECT_TRUE(stops.count("output/time")) << mode;
    }
}


// The retrieval's counters (review GPT-3, for the route to state), read from PatternRetrieval
// itself on paths built by hand: the distinct rows read (one for a homopolymer's path, however
// long), the verification's work (a few units per k-mer for its one run of chains, not one per
// chain and k-mer), and the request's sums over its patterns
TEST(PatternPaths, TheRetrievalCounters) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRepeats, true);
    const DeBruijnGraph &graph = idx.anno->get_graph();
    // a path of |s| as the route collects it: the annotation key (row + 1) of every k-mer
    auto path_of = [&](const std::string &s) {
        RetrievalPath p;
        p.sequence = s;
        for (size_t j = 0; j + kK <= s.size(); ++j) {
            DeBruijnGraph::node_index node = DeBruijnGraph::npos;
            graph.map_to_nodes(s.substr(j, kK), [&](DeBruijnGraph::node_index x) { node = x; });
            EXPECT_NE(DeBruijnGraph::npos, node) << s.substr(j, kK);
            p.keys.push_back(AnnotatedDBG::anno_to_graph_index(
                    AnnotatedDBG::graph_to_anno_index(node)));
        }
        return p;
    };
    RetrievalLimits limits;
    pattern::Budget budget(1'000'000, mtg::test::unbounded_deadline());
    RetrievalHooks hooks;
    hooks.coord_to_header = idx.cth.get();
    PatternRetrieval retrieval(*idx.anno, pattern::GraphMode::BASIC, limits, budget, &hooks);
    pattern::Extraction x;
    x.complete = true;
    auto labels = [&](const std::string &s) {
        const std::vector<RetrievalPath> paths = { path_of(s) };
        retrieval.begin_release(pattern::Mode::PARTIAL);
        EXPECT_TRUE(retrieval.admit_path(s.size()));
        x.returned = 1;
        return retrieval.retrieve_paths(paths, 1, s.size(), pattern::Mode::PARTIAL, x,
                                        Json::Value(), false);
    };
    const std::string a50(50, 'A');
    LabelsAnswer h = labels(a50);
    ASSERT_EQ(1u, h.result_fields.size());
    EXPECT_TRUE(h.complete || !h.fields["occurrences_cut"].isNull()) << h.fields;
    EXPECT_EQ(1u, h.counters.rows_distinct);
    // 46 k-mers, two labels: h's 412 coordinates are one run (its rows looked up, the lists
    // ordered, sought and extended once each, a unit per record), m's three runs broken by C's
    // a few candidates each; a search per coordinate and k-mer would be 412 x 45 for h alone
    const uint64_t n = 50 - kK + 1;
    EXPECT_GT(h.counters.verification_steps, 4 * n);
    EXPECT_LT(h.counters.verification_steps, 412 * (n - 1) / 10) << h.counters.verification_steps;
    EXPECT_GE(h.counters.label_intersection_ms, 0);
    EXPECT_GE(h.counters.verification_ms, 0);

    // m0's C A^25 C: the rows CAAAA, AAAAA and AAAAC
    LabelsAnswer c = labels("C" + std::string(25, 'A') + "C");
    EXPECT_EQ(3u, c.counters.rows_distinct);
    const RetrievalCounters &total = retrieval.counters();
    EXPECT_EQ(4u, total.rows_distinct);
    EXPECT_EQ(h.counters.verification_steps + c.counters.verification_steps,
              total.verification_steps);

    // a work budget that one row passes: the others not read, not counted
    limits.max_annotation_work = 1;
    pattern::Budget small(1'000'000, mtg::test::unbounded_deadline());
    PatternRetrieval one(*idx.anno, pattern::GraphMode::BASIC, limits, small, &hooks);
    const std::string s = "C" + std::string(25, 'A') + "C";
    const std::vector<RetrievalPath> paths = { path_of(s) };
    one.begin_release(pattern::Mode::PARTIAL);
    ASSERT_TRUE(one.admit_path(s.size()));
    LabelsAnswer w = one.retrieve_paths(paths, 1, s.size(), pattern::Mode::PARTIAL, x,
                                        Json::Value(), false);
    ASSERT_TRUE(w.stop);
    EXPECT_EQ("max_annotation_work", w.stop->second);
    EXPECT_EQ(1u, w.work["annotation_rows"].asUInt64());
    EXPECT_EQ(1u, w.counters.rows_distinct);
}


// ------------------------------------------------------------------ peptides (increment 5)
//
// The route's protein kind (owner decision #15 of 2026-10-08, wired by the integration of
// increments 4 and 5): patterns[i].protein read in the request's genetic_code, the entry's
// kind, residues and genetic_code, its slot errors, and its contexts (3m <= k) and its paths
// (3m > k, long_search "paths") with their labels, against oracles that never ask the engine
// nor its tables:
//  - the genetic codes are the test's own copy of NCBI's ncbieaa strings (gc.prt 4.6) for the
//    tables used here (1, 2, 11; 27 and 31, whose stops code a residue unless in context:
//    owner decisions #19 and #21);
//  - a GRAPH-WALK oracle over the k-mers of the records (both orientations on CANONICAL and
//    PRIMARY graphs): every k-mer window, and every walk of k-mers, whose bases translate to
//    the peptide (forward), or whose reverse complement does (reverse);
//  - a SIX-FRAME oracle over the records: each record translated in its three frames and the
//    three of its reverse complement; every match of the peptide is a placed occurrence
//    (column, seq_id, 1-based start on the record, strand).

// NCBI's ncbieaa: the residue of each codon in TCAG order (TTT TTC TTA TTG TCT ... GGG)
const std::map<int, std::string> kTables = {
    { 1, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    { 2, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG" },
    { 11, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    // Karyorelict: TAA, TAG Q and TGA W (a stop in context only)
    { 27, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
    // Blastocrithidia: TAA, TAG E (a stop in context only), TGA W
    { 31, "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG" },
};

char translate(int table, const std::string &codon) {
    const std::string order = "TCAG";
    size_t i = 0;
    for (char b : codon) {
        i = 4 * i + order.find(b);
    }
    return kTables.at(table).at(i);
}

const std::vector<std::string>& all_codons() {
    static const std::vector<std::string> codons = [] {
        std::vector<std::string> out;
        for (char a : std::string("ACGT")) {
            for (char b : std::string("ACGT")) {
                for (char c : std::string("ACGT")) {
                    out.push_back({ a, b, c });
                }
            }
        }
        return out;
    }();
    return codons;
}

// a residue of a peptide (X, B, Z, J the ambiguity codes) admits the amino acid |aa|; a stop
// is admitted by the stop '*' only (owner decision #19), and '*' admits nothing else
bool residue_admits(char residue, char aa) {
    if (residue == '*' || aa == '*')
        return residue == aa;
    switch (residue) {
        case 'X': return true;
        case 'B': return aa == 'D' || aa == 'N';
        case 'Z': return aa == 'E' || aa == 'Q';
        case 'J': return aa == 'I' || aa == 'L';
        default: return aa == residue;
    }
}

struct Peptide {
    std::string residues;
    int table = 1;

    size_t length() const { return 3 * residues.size(); }

    /**
     * |s| (at most length() bases) is the prefix of an instance of the oriented peptide: each
     * complete codon translates to its residue and a partial last codon is the prefix of one
     * that does. The reverse orientation is rc(P): its i-th codon is the reverse complement of
     * a codon of residue m - 1 - i.
     */
    bool prefix(const std::string &s, bool reverse) const {
        const size_t m = residues.size();
        if (s.size() > 3 * m)
            return false;
        for (size_t i = 0; 3 * i < s.size(); ++i) {
            const std::string part = s.substr(3 * i, 3);
            const char residue = residues[reverse ? m - 1 - i : i];
            bool any = false;
            for (const std::string &codon : all_codons()) {
                const std::string oriented = reverse ? rc(codon) : codon;
                if (oriented.compare(0, part.size(), part) == 0
                        && residue_admits(residue, translate(table, codon))) {
                    any = true;
                    break;
                }
            }
            if (!any)
                return false;
        }
        return true;
    }

    bool instance(const std::string &s, bool reverse) const {
        return s.size() == length() && prefix(s, reverse);
    }

    // log2 of 64 / (the codons of the residue), summed: the bits of the whole peptide; a
    // residue without a codon (a stop '*' in a table without one) counts as one codon, 6 bits
    double bits() const {
        double out = 0;
        for (char r : residues) {
            size_t n = 0;
            for (const std::string &codon : all_codons()) {
                n += residue_admits(r, translate(table, codon));
            }
            out += std::log2(64.0 / std::max<size_t>(n, 1));
        }
        return out;
    }

    std::string json() const {
        return "{\"protein\": \"" + residues + "\"}";
    }
};

// the request of |p| with long_search "paths" (as body()) in its genetic code
std::string peptide_body(const Peptide &p, const std::string &rest = "",
                         const std::string &labels = "all") {
    return body(p.json(), "\"genetic_code\": " + std::to_string(p.table)
                          + (rest.empty() ? "" : ", " + rest), labels);
}

using ContextSet = std::set<std::tuple<std::string, std::string, uint64_t>>;

// the graph-walk oracle's contexts of a peptide of at most k bases: (strand or orientation,
// k-mer, offset)
ContextSet peptide_contexts(const Walk &walk, const Peptide &p, bool strand_stated) {
    ContextSet out;
    const size_t L = p.length();
    for (const std::string &kmer : walk.kmers) {
        for (size_t o = 0; o + L <= walk.k; ++o) {
            const std::string s = kmer.substr(o, L);
            if (p.instance(s, false))
                out.emplace(strand_stated ? "+" : "forward", kmer, o);
            if (p.instance(s, true))
                out.emplace(strand_stated ? "-" : "reverse", kmer, o);
        }
    }
    return out;
}

// the graph-walk oracle's paths of a peptide longer than k in one orientation: every walk of
// k-mers spelling an instance, extended base by base from each anchor; |entered| counts the
// prefixes of k + 1 .. L bases so formed (the branches the engine's DFS enters, §4.2)
std::vector<std::string> peptide_paths(const Walk &walk, const Peptide &p, bool reverse,
                                       uint64_t *entered = nullptr,
                                       uint64_t *anchors = nullptr) {
    std::vector<std::string> out;
    uint64_t branches = 0, starts = 0;
    std::function<void(const std::string&)> dfs = [&](const std::string &s) {
        if (s.size() == p.length()) {
            out.push_back(s);
            return;
        }
        for (char b : std::string("ACGT")) {
            const std::string t = s + b;
            if (!p.prefix(t, reverse) || !walk.kmers.count(t.substr(t.size() - walk.k)))
                continue;
            ++branches;
            dfs(t);
        }
    };
    for (const std::string &x : walk.kmers) {
        if (p.prefix(x, reverse)) {
            ++starts;
            dfs(x);
        }
    }
    if (entered)
        *entered = branches;
    if (anchors)
        *anchors = starts;
    return out;
}

using Placed = std::set<std::tuple<std::string, uint64_t, uint64_t, std::string>>;

// the six-frame oracle: {(column, seq_id, 1-based start, strand)} of every match of |p| in
// the translation of every record, three frames on each strand
Placed six_frames(const Index &idx, const Peptide &p) {
    Placed out;
    std::map<std::string, uint64_t> seq_id;
    const size_t m = p.residues.size(), L = p.length();
    for (const Record &r : idx.records) {
        const uint64_t id = seq_id[r.column]++;
        for (const std::string strand : { "+", "-" }) {
            const std::string t = strand == "+" ? r.seq : rc(r.seq);
            for (size_t frame = 0; frame < 3; ++frame) {
                std::string aa;
                for (size_t i = frame; i + 3 <= t.size(); i += 3) {
                    aa += translate(p.table, t.substr(i, 3));
                }
                for (size_t j = 0; j + m <= aa.size(); ++j) {
                    bool match = true;
                    for (size_t x = 0; x < m && match; ++x) {
                        match = residue_admits(p.residues[x], aa[j + x]);
                    }
                    if (!match)
                        continue;
                    const size_t start = frame + 3 * j;
                    out.emplace(r.column, id,
                                strand == "+" ? start + 1 : t.size() - (start + L) + 1, strand);
                }
            }
        }
    }
    return out;
}

// the entry of a peptide: its kind, residues, genetic code and bits (the test's own count of
// the codons of each residue)
void check_peptide_entry(const Json::Value &e, const Peptide &p) {
    EXPECT_EQ("protein", e["kind"].asString());
    EXPECT_EQ(p.residues, e["pattern"].asString());
    EXPECT_EQ(p.length(), e["length"].asUInt64());
    EXPECT_EQ(p.residues.size(), e["residues"].asUInt64());
    EXPECT_EQ(p.table, e["genetic_code"].asInt());
    EXPECT_NEAR(p.bits(), e["information_bits"].asDouble(), 1e-9) << e;
    EXPECT_FALSE(e.isMember("error")) << e;
    EXPECT_FALSE(e["palindromic"].asBool());
}

/**
 * A peptide of at most k bases on a BASIC index with record placement, labels "all": the
 * contexts are the graph-walk oracle's, each with the columns whose records hold its k-mer,
 * and the placed occurrences, over all contexts, the six-frame oracle's matches.
 */
void check_peptide_contexts(const Index &idx, const Json::Value &e, const Peptide &p) {
    check_peptide_entry(e, p);
    const Walk walk(idx);
    const Scan scan(idx);
    const ContextSet expected = peptide_contexts(walk, p, true);
    ContextSet got;
    Placed placed;
    for (const Json::Value &r : e["results"]) {
        const std::string kmer = r["kmer"].asString();
        const uint64_t offset = r["offset"].asUInt64();
        got.emplace(r["strand"].asString(), kmer, offset);
        EXPECT_EQ(kmer.substr(offset, p.length()), r["instance"].asString());
        std::set<std::string> columns;
        for (const Json::Value &l : r["labels"]) {
            columns.insert(l["column"].asString());
            for (const Json::Value &o : l["occurrence_list"]) {
                const std::string coords = o["nt_coords"].asString();
                placed.emplace(l["column"].asString(), o["seq_id"].asUInt64(),
                               std::stoull(coords.substr(0, coords.find('-'))),
                               o["strand"].asString());
            }
        }
        std::set<std::string> holding;
        for (const auto &[column, recs] : scan.columns) {
            if (scan.holds_kmer(column, kmer))
                holding.insert(column);
        }
        EXPECT_EQ(holding, columns) << kmer;
    }
    EXPECT_EQ(expected, got);
    EXPECT_EQ(expected.size(), e["counts"]["contexts"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["contexts"]["relation"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    const Placed frames = six_frames(idx, p);
    EXPECT_EQ(frames, placed);
    EXPECT_EQ(frames.size(), e["counts"]["occurrences"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["occurrences"]["relation"].asString());
}

/**
 * A peptide longer than k under long_search "paths" on a BASIC index with record placement,
 * labels "all": the paths and their counts the graph-walk oracle's, each path's labels the
 * columns carrying it (or, with |require_verified|, those holding it whole in a record), each
 * label's support and occurrences the record scan's, and the placed occurrences over all paths
 * the six-frame oracle's matches.
 */
void check_peptide_paths(const Index &idx, const Json::Value &e, const Peptide &p,
                         bool require_verified) {
    check_peptide_entry(e, p);
    const Walk walk(idx);
    const Scan scan(idx);
    const size_t k = idx.k, L = p.length();
    std::set<std::pair<std::string, std::string>> expected;
    uint64_t candidates = 0, anchors = 0;
    for (bool reverse : { false, true }) {
        uint64_t entered = 0, starts = 0;
        const std::string strand = reverse ? "-" : "+";
        const auto paths = peptide_paths(walk, p, reverse, &entered, &starts);
        for (const std::string &s : paths) {
            expected.emplace(strand, s);
        }
        candidates += entered;
        anchors += starts;
        EXPECT_EQ(paths.size(), e["counts"]["paths"]["by_strand"][strand]["value"].asUInt64());
        EXPECT_EQ(starts, e["counts"]["anchors"]["by_strand"][strand]["value"].asUInt64());
    }
    const Json::Value &c = e["counts"]["paths"];
    EXPECT_EQ(expected.size(), c["value"].asUInt64()) << c;
    EXPECT_EQ("exact", c["relation"].asString());
    EXPECT_EQ(candidates, c["candidates_examined"].asUInt64()) << c;
    EXPECT_EQ(anchors ? "completed" : "no_anchors", c["extension"].asString());
    EXPECT_EQ(anchors, e["counts"]["anchors"]["value"].asUInt64());
    EXPECT_EQ("long", e["scope"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;

    std::set<std::pair<std::string, std::string>> got;
    Placed placed;
    for (const Json::Value &r : e["results"]) {
        const std::string s = r["sequence"].asString();
        const std::string strand = r["strand"].asString();
        got.emplace(strand, s);
        EXPECT_FALSE(r.isMember("kmer"));
        EXPECT_EQ(s, r["instance"].asString());
        EXPECT_EQ(s.substr(0, k), r["anchor_kmer"].asString());
        EXPECT_EQ(L - k + 1, r["nodes"].size());
        const std::set<std::string> carriers = scan.carriers(s);
        std::set<std::string> verified;
        for (const std::string &column : carriers) {
            if (!scan.occurrences(column, s).empty())
                verified.insert(column);
        }
        EXPECT_EQ(carriers.size(), r["labels_total"].asUInt64()) << r;
        std::set<std::string> listed;
        for (const Json::Value &l : r["labels"]) {
            const std::string column = l["column"].asString();
            listed.insert(column);
            const auto occ = scan.occurrences(column, s);
            EXPECT_EQ(occ.empty() ? "label_intersection" : "record_verified",
                      l["support"].asString()) << s << " " << column;
            EXPECT_EQ(occ.size(), l["occurrences"]["value"].asUInt64());
            for (const Json::Value &o : l["occurrence_list"]) {
                const std::string coords = o["nt_coords"].asString();
                const uint64_t start = std::stoull(coords.substr(0, coords.find('-')));
                EXPECT_TRUE(occ.count({ o["seq_id"].asUInt64(), start })) << coords;
                placed.emplace(column, o["seq_id"].asUInt64(), start, o["strand"].asString());
            }
        }
        EXPECT_EQ(require_verified ? verified : carriers, listed) << s;
    }
    EXPECT_EQ(expected, got);
    // every placed occurrence is a whole path in one record: the six-frame matches
    EXPECT_EQ(six_frames(idx, p), placed);
    EXPECT_EQ(placed.size(), e["counts"]["occurrences"]["value"].asUInt64());
}

// random records of four columns with planted codings of one peptide, in synonymous codons,
// on both strands and split across two adjacent records of one column
struct PeptideRecords {
    std::vector<Record> records;
    std::string peptide;
};

PeptideRecords peptide_records(std::mt19937 &rng, size_t m, size_t k) {
    auto random_seq = [&](size_t length) {
        std::string s;
        for (size_t i = 0; i < length; ++i) {
            s += "ACGT"[rng() % 4];
        }
        return s;
    };
    // a peptide with single-, two-, four- and six-codon residues
    const std::string pool = "MWKHLSRAGVEDNQ";
    std::string peptide;
    for (size_t i = 0; i < m; ++i) {
        peptide += pool[rng() % pool.size()];
    }
    auto coding = [&]() {
        std::string out;
        for (char residue : peptide) {
            std::vector<std::string> codons;
            for (const std::string &codon : all_codons()) {
                if (translate(1, codon) == residue)
                    codons.push_back(codon);
            }
            out += codons[rng() % codons.size()];
        }
        return out;
    };
    PeptideRecords out;
    out.peptide = peptide;
    for (int i = 0; i < 10; ++i) {
        std::string seq = random_seq(30 + rng() % 40);
        if (i % 3 == 0) {
            const std::string planted = i % 2 ? rc(coding()) : coding();
            seq = seq.substr(0, 10) + planted + seq.substr(10);
        }
        if (i == 4) {
            // the record's end and the next record's start in column c1 (seq_ids 1 and 2)
            // spell one coding: the first record ends with its k-mers 0 .. h - 1, the next
            // starts with its k-mers h .., so that their column coordinates are consecutive
            // across the two records, and neither holds it whole
            const std::string planted = coding();
            const size_t h = (planted.size() - k + 1) / 2;
            seq += planted.substr(0, h + k - 1);
            out.records.push_back({ "c1", "r" + std::to_string(i), seq });
            out.records.push_back({ "c1", "r" + std::to_string(i) + "b",
                                    planted.substr(h) + random_seq(20) });
            continue;
        }
        out.records.push_back({ std::string("c") + char('0' + i % 4), "r" + std::to_string(i),
                                seq });
    }
    return out;
}

TEST(PatternRoutePeptide, ShortPeptidesAgainstTheOracles) {
    std::mt19937 rng(815);
    for (int trial = 0; trial < 6; ++trial) {
        const PeptideRecords data = peptide_records(rng, 4, 13);
        Index idx = build<annot::RowDiffColumnAnnotator>(13, data.records, true);
        std::vector<Peptide> peptides = { { data.peptide, 1 }, { data.peptide, 11 },
                                          { data.peptide.substr(1), 1 } };
        // widened by an ambiguity code
        Peptide x { data.peptide, 1 };
        x.residues[2] = "XBZJ"[rng() % 4];
        peptides.push_back(x);
        // windows of the records' translations (whatever they code), in table 2 too
        for (const Record &r : data.records) {
            const size_t start = rng() % (r.seq.size() - 12);
            std::string residues;
            for (size_t i = 0; i < 4; ++i) {
                residues += translate(trial % 2 ? 2 : 1, r.seq.substr(start + 3 * i, 3));
            }
            if (residues.find('*') == std::string::npos)
                peptides.push_back({ residues, trial % 2 ? 2 : 1 });
        }
        for (const Peptide &p : peptides) {
            SCOPED_TRACE(p.residues + " table " + std::to_string(p.table));
            Json::Value out = run(idx, "{\"patterns\": [" + p.json() + "], \"genetic_code\": "
                                       + std::to_string(p.table) + ", \"output\": {\"labels\": "
                                       "\"all\"}}");
            const Json::Value &e = out["patterns"][0];
            check_peptide_contexts(idx, e, p);
            // a peptide of at most k bases: long_search changes nothing
            Json::Value again = run(idx, peptide_body(p));
            EXPECT_EQ(untimed(out["patterns"]), untimed(again["patterns"]));
        }
        // the planted peptide is found on both strands
        Json::Value out = run(idx, "{\"patterns\": [{\"protein\": \"" + data.peptide + "\"}], "
                                   "\"mode\": \"count\"}");
        const Json::Value &by = out["patterns"][0]["counts"]["contexts"]["by_strand"];
        EXPECT_GT(by["+"]["value"].asUInt64(), 0u) << out;
        EXPECT_GT(by["-"]["value"].asUInt64(), 0u) << out;
    }
}

TEST(PatternRoutePeptide, LongPeptidesAsPathsAgainstTheOracles) {
    std::mt19937 rng(4242);
    size_t verified_seen = 0, intersection_seen = 0;
    for (int trial = 0; trial < 6; ++trial) {
        const PeptideRecords data = peptide_records(rng, 7, 11);
        Index idx = build<annot::RowDiffColumnAnnotator>(11, data.records, true);
        std::vector<Peptide> peptides = { { data.peptide, 1 }, { data.peptide.substr(1), 11 },
                                          { data.peptide.substr(0, 5), 1 } };
        Peptide x { data.peptide, 1 };
        x.residues[4] = 'X';
        peptides.push_back(x);
        for (const Peptide &p : peptides) {
            ASSERT_GT(p.length(), idx.k);
            for (bool require : { false, true }) {
                SCOPED_TRACE(p.residues + " table " + std::to_string(p.table)
                             + (require ? " record_verified" : ""));
                Json::Value out = run(idx, peptide_body(p, require ? "\"require_support\": "
                                                                     "\"record_verified\"" : ""));
                const Json::Value &e = out["patterns"][0];
                check_peptide_paths(idx, e, p, require);
                for (const Json::Value &r : e["results"]) {
                    for (const Json::Value &l : r["labels"]) {
                        verified_seen += l["support"].asString() == "record_verified";
                        intersection_seen += l["support"].asString() == "label_intersection";
                    }
                }
            }
            // without long_search "paths": the anchors only, as any pattern longer than k
            Json::Value plain = run(idx, "{\"patterns\": [" + p.json() + "], \"genetic_code\": "
                                         + std::to_string(p.table) + "}");
            const Json::Value &e = plain["patterns"][0];
            ASSERT_GT(e["counts"]["anchors"]["value"].asUInt64(), 0u) << e;
            EXPECT_EQ("unknown", e["counts"]["paths"]["relation"].asString());
            EXPECT_EQ("paths_later_increment", e["withheld"]["reason"].asString());
            EXPECT_EQ(0u, e["results"].size());
            EXPECT_FALSE(plain["limits"].isMember("long_search"));
        }
    }
    // both supports were seen: the planted coding held whole, and the one split across two
    // records of c1
    EXPECT_GT(verified_seen, 0u);
    EXPECT_GT(intersection_seen, 0u);
}

TEST(PatternRoutePeptide, CanonicalAndPrimaryGraphs) {
    std::mt19937 rng(97);
    const PeptideRecords data = peptide_records(rng, 6, 11);
    for (auto mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        Index idx = build<annot::ColumnCompressed<>>(11, data.records, false, mode);
        const Walk walk(idx);
        for (const Peptide &p : std::vector<Peptide> { { data.peptide.substr(0, 3), 1 },
                                                       { data.peptide, 1 },
                                                       { data.peptide.substr(1), 11 } }) {
            SCOPED_TRACE(p.residues + (mode == DeBruijnGraph::PRIMARY ? " primary"
                                                                       : " canonical"));
            Json::Value out = run(idx, peptide_body(p, "", "none"), {}, false);
            const Json::Value &e = out["patterns"][0];
            check_peptide_entry(e, p);
            EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
            if (p.length() <= idx.k) {
                ContextSet got;
                for (const Json::Value &r : e["results"]) {
                    got.emplace(r["orientation"].asString(), r["kmer"].asString(),
                                r["offset"].asUInt64());
                }
                EXPECT_EQ(peptide_contexts(walk, p, false), got);
                continue;
            }
            std::set<std::pair<std::string, std::string>> expected, got;
            for (bool reverse : { false, true }) {
                for (const std::string &s : peptide_paths(walk, p, reverse)) {
                    expected.emplace(reverse ? "reverse" : "forward", s);
                }
            }
            for (const Json::Value &r : e["results"]) {
                got.emplace(r["orientation"].asString(), r["sequence"].asString());
            }
            EXPECT_GT(expected.size(), 0u);
            EXPECT_EQ(expected, got);
            EXPECT_EQ(expected.size(), e["counts"]["paths"]["value"].asUInt64());
        }
    }
}

TEST(PatternRoutePeptide, TheEntryItsSlotsAndTheGeneticCode) {
    // ATG TGA AAA: M W K in table 2 (TGA a W), a stop in tables 1 and 11
    const std::vector<Record> records = {
        { "a", "a0", "CCCATGTGAAAAGGGTTTCCC" },
        { "b", "b0", "GGGATGTGGAAATTTGGG" },
    };
    Index idx = build<annot::RowDiffColumnAnnotator>(13, records, true);
    PatternLimits floor = limits();
    floor.min_information_bits = 16;
    Json::Value out = run(idx, "{\"patterns\": [{\"id\": \"mwk\", \"protein\": \"mwk\"}, "
                               "{\"id\": \"stop\", \"protein\": \"M*K\"}, "
                               "{\"id\": \"bad\", \"protein\": \"M*U\"}, {\"protein\": \"\"}, "
                               "{\"id\": \"floor\", \"protein\": \"MK\"}], "
                               "\"mode\": \"count\"}", {}, true, floor);
    const Json::Value &mwk = out["patterns"][0];
    check_peptide_entry(mwk, { "MWK", 1 });
    EXPECT_EQ("mwk", mwk["id"].asString());
    const Walk walk(idx);
    // ATG TGG AAA in b0, in the 4 k-mers holding it (offsets 0 to 3)
    EXPECT_EQ(4u, peptide_contexts(walk, { "MWK", 1 }, true).size());
    EXPECT_EQ(4u, mwk["counts"]["contexts"]["value"].asUInt64()) << mwk;
    for (int i : { 2, 3 }) {
        const Json::Value &slot = out["patterns"][i];
        EXPECT_EQ(std::vector<std::string>({ "error", "id", "kind" }), slot.getMemberNames());
        EXPECT_EQ("protein", slot["kind"].asString());
        EXPECT_EQ("bad_alphabet", slot["error"]["code"].asString()) << slot;
    }
    EXPECT_NE(std::string::npos, out["patterns"][2]["error"]["message"].asString()
                                         .find("'U' at position 2"));
    // owner decision #19: the stop '*' is a residue (stop_unsupported is gone): M*K is read,
    // its bits those of M, a stop codon (3 in table 1) and K, below this floor of 16
    const Json::Value &stop = out["patterns"][1];
    EXPECT_EQ("information_below_floor", stop["error"]["code"].asString()) << stop;
    EXPECT_EQ("M*K", stop["pattern"].asString());
    EXPECT_EQ(3u, stop["residues"].asUInt64());
    EXPECT_NEAR(6 + std::log2(64.0 / 3) + 5, stop["information_bits"].asDouble(), 1e-9);
    // below the floor (M 6 bits + K 5): the slot keeps the peptide's description
    const Json::Value &low = out["patterns"][4];
    EXPECT_EQ("information_below_floor", low["error"]["code"].asString());
    EXPECT_EQ(2u, low["residues"].asUInt64());
    EXPECT_EQ(6u, low["length"].asUInt64());
    EXPECT_EQ(1, low["genetic_code"].asInt());
    EXPECT_NEAR(11.0, low["information_bits"].asDouble(), 1e-9);

    // the genetic code changes the hit: MWK in table 2 also reads ATG TGA AAA
    for (int table : { 1, 2, 11 }) {
        const Peptide p { "MWK", table };
        Json::Value answer = run(idx, "{\"patterns\": [{\"protein\": \"MWK\"}], \"genetic_code\": "
                                      + std::to_string(table) + ", \"output\": {\"labels\": "
                                      "\"all\"}}");
        check_peptide_contexts(idx, answer["patterns"][0], p);
        EXPECT_EQ(table == 2 ? 2u : 1u, six_frames(idx, p).size()) << table;
    }
    // an unknown table refuses the request
    EXPECT_EQ(std::make_pair(400, std::string("genetic_code_unknown")),
              refusal(idx, "{\"patterns\": [{\"protein\": \"MWK\"}], \"genetic_code\": 19}"));
}


// Owner decisions #19 and #21: the stop '*' of a peptide is a stop codon of the request's
// genetic code at that position (the codons the table translates to '*'), against the
// graph-walk and six-frame oracles in tables 1, 2 and 11 (whose stops differ: table 2 stops at
// AGA and AGG and reads TGA as W); X never matches a stop; in tables 27 and 31, whose stop
// codons code a residue unless in context, those codons match as the residue and '*' matches
// nothing: the counts are exact 0 and the entry says why (no_stop_codon)
TEST(PatternRoutePeptide, TheStopIsAStopCodonOfTheTable) {
    const std::vector<Record> records = {
        // ATG TGA AAA: M * K (tables 1, 11), M W K (2, 27, 31)
        { "a", "a0", "CCCATGTGAAAAGGGTTTCCC" },
        // ATG AGG AAA: M R K (1, 11, 27, 31), M * K (2)
        { "b", "b0", "GGGATGAGGAAATTTGGG" },
        // ATG TAA AAA: M * K (1, 2, 11), M Q K (27), M E K (31)
        { "c", "c0", "TTTATGTAAAAACCCAAAT" },
        // ATG TAG AAA on the minus strand: M * K (1, 2, 11)
        { "d", "d0", "GGGTTTCTACATGGG" },
    };
    Index idx = build<annot::RowDiffColumnAnnotator>(13, records, true);
    const Walk walk(idx);
    const std::vector<Peptide> peptides = {
        { "M*K", 1 }, { "M*K", 2 }, { "M*K", 11 }, { "M*", 1 }, { "*K", 11 }, { "MX*", 11 },
        { "M*X", 2 }, { "MXK", 1 }, { "MXK", 27 }, { "MQK", 27 }, { "MWK", 27 }, { "MEK", 31 },
        { "M*K", 27 }, { "M*K", 31 },
    };
    for (const Peptide &p : peptides) {
        SCOPED_TRACE(p.residues + " in table " + std::to_string(p.table));
        const bool none = (p.table == 27 || p.table == 31)
                            && p.residues.find('*') != std::string::npos;
        const Json::Value answer = run(idx, "{\"patterns\": [" + p.json() + "], "
                                            "\"genetic_code\": " + std::to_string(p.table)
                                            + ", \"output\": {\"labels\": \"all\"}}");
        const Json::Value &e = answer["patterns"][0];
        check_peptide_contexts(idx, e, p);
        const ContextSet expected = peptide_contexts(walk, p, true);
        if (none) {
            EXPECT_TRUE(expected.empty());
        }
        bool noted = false;
        for (const Json::Value &n : e["notes"]) {
            noted |= n.asString() == "no_stop_codon";
        }
        EXPECT_EQ(none, noted) << e;
        if (none) {
            EXPECT_EQ(0u, e["counts"]["contexts"]["value"].asUInt64());
            EXPECT_EQ(0u, e["returned"].asUInt64());
        }
        // the count mode states the same count
        const Json::Value count = run(idx, "{\"patterns\": [" + p.json() + "], \"mode\": "
                                           "\"count\", \"genetic_code\": "
                                           + std::to_string(p.table) + "}")["patterns"][0];
        EXPECT_EQ(e["counts"]["contexts"], count["counts"]["contexts"]);
    }
    // the stop alone in a table without one: no instance at all (palindromic vacuously, as
    // its reverse complement has none either), exact 0, the note
    const Json::Value lone = run(idx, "{\"patterns\": [{\"protein\": \"*\"}], "
                                      "\"genetic_code\": 27, \"output\": {\"labels\": "
                                      "\"all\"}}")["patterns"][0];
    EXPECT_EQ("exact", lone["counts"]["contexts"]["relation"].asString()) << lone;
    EXPECT_EQ(0u, lone["counts"]["contexts"]["value"].asUInt64());
    EXPECT_TRUE(lone["retrieval_complete"].asBool());
    EXPECT_EQ(0u, lone["returned"].asUInt64());
    EXPECT_NE(lone["notes"].end(), std::find(lone["notes"].begin(), lone["notes"].end(),
                                             Json::Value("no_stop_codon")));
    // the hits the records were written for: M*K at TGA, TAA, TAG in tables 1 and 11 and at
    // AGG, TAA, TAG in table 2; X never a stop (MXK in table 1: MRK only), every codon a
    // residue in table 27 (MXK four times); MQK and MEK at the context stops TAA and TAG
    EXPECT_EQ(3u, six_frames(idx, { "M*K", 1 }).size());
    EXPECT_EQ(3u, six_frames(idx, { "M*K", 2 }).size());
    EXPECT_EQ(1u, six_frames(idx, { "MXK", 1 }).size());
    EXPECT_EQ(4u, six_frames(idx, { "MXK", 27 }).size());
    EXPECT_EQ(2u, six_frames(idx, { "MQK", 27 }).size());
    EXPECT_EQ(2u, six_frames(idx, { "MEK", 31 }).size());
}

} // namespace
