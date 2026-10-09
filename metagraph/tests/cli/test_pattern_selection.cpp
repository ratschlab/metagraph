#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"
#include "../graph/alignment/pattern_test_support.hpp"
#include "../test_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "cli/pattern_predicate.hpp"
#include "cli/pattern_retrieval.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// A predicate's selection for patterns of L <= k (increments 5b-3 and 5b-4; SPEC-pattern-search.md
// §19, TESTS §3), driven through THE ROUTE (process_pattern_request, src/cli/pattern.cpp): the
// request's predicate bound once, each pattern's raw contexts released into the selection pass
// (src/cli/pattern_selection.cpp), its rows read under max_predicate_work, the memory account
// and the deadline, the relations of §19.7, the selection admission, the results of the
// selected contexts and their projection (none, predicate_only, all), as the answer states
// them. The expectations come from ORACLES over the records (never the modules under test):
// the columns whose records hold a k-mer as deposited (and its reverse complement for
// "either"), the graph contexts of a pattern scanned from the records' k-mers, the placed
// occurrences from the records, and a recursive evaluator of the request's predicate JSON
// (unknown names are simply absent from a set, which is what the folding of §19.4 computes).
// The order in which the pass reads rows is the engine's answer order of the raw contexts,
// taken from the route's own release of the pattern without a predicate (labels none).

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;
namespace gp = mtg::graph::pattern;
using Clock = gp::Deadline::Clock;

struct Record {
    std::string column;
    std::string header;
    std::string seq;
};

std::string rc(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

// An index of |records| (as test_pattern_retrieval.cpp builds them): a BASIC (or |mode|)
// succinct graph at |k| with the dummy-edge mask, an annotation of type |Annotation| (with
// coordinates when |coordinates|) and the record mapping built from the same records
struct Index {
    size_t k = 0;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
};

template <class Annotation>
Index build(size_t k, const std::vector<Record> &records, bool coordinates,
            DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
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
            k, seqs, labels, mode, coordinates, coordinates ? starts : std::vector<uint64_t>{});
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

// k = 7, the records of test_pattern_retrieval.cpp's kRecords and c5, which holds c1/r0's
// motif on the - strand only (its reverse complement): AC's + contexts of GGACGACTTTG are not
// on c5's rows, their reverse complements are
const size_t kK = 7;
const std::vector<Record> kR1 = {
    { "c1", "c1_r0", "GGACGACTTTG" },
    { "c1", "c1_r1", "CCCTGTCCCCA" },
    { "c2", "c2_r0", "TTTTACTTTTG" },
    { "c2", "c2_r1", "GGGGACGGGGC" },
    { "c3", "c3_r0", "ACGTACGTACGTACG" },
    { "c4", "c4_r0", "GGACGACTTTG" },
    { "c5", "c5_r0", "CAAAGTCGTCC" },
};

// |columns| columns r00, r01, ... of one or two random records each (seeded)
std::vector<Record> random_records(size_t columns, unsigned seed) {
    std::mt19937 rng(seed);
    std::vector<Record> out;
    for (size_t c = 0; c < columns; ++c) {
        const std::string name = (c < 10 ? "r0" : "r") + std::to_string(c);
        const size_t n = 1 + rng() % 2;
        for (size_t r = 0; r < n; ++r) {
            std::string s(12 + rng() % 19, 'A');
            for (char &x : s) {
                x = "ACGT"[rng() % 4];
            }
            out.push_back({ name, name + "_" + std::to_string(r), s });
        }
    }
    return out;
}

// ---------------------------------------------------------------- the oracles

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
    std::string out(p.rbegin(), p.rend());
    for (char &x : out) {
        x = c.at(x);
    }
    return out;
}

// O-eval: the request's predicate on a set of column names, by recursion over its JSON
bool o_eval(const Json::Value &p, const std::set<std::string> &s) {
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &v = p[op];
    auto in = [&](const Json::Value &n) { return s.count(n.asString()) > 0; };
    if (op == "any" || op == "all" || op == "none") {
        uint64_t hits = 0;
        for (const Json::Value &n : v) {
            hits += in(n);
        }
        return op == "any" ? hits > 0 : op == "all" ? hits == v.size() : hits == 0;
    }
    if (op == "at_least") {
        uint64_t hits = 0;
        for (const Json::Value &n : v["labels"]) {
            hits += in(n);
        }
        return hits >= v["n"].asUInt64();
    }
    if (op == "and") {
        for (const Json::Value &q : v) {
            if (!o_eval(q, s))
                return false;
        }
        return true;
    }
    if (op == "or") {
        for (const Json::Value &q : v) {
            if (o_eval(q, s))
                return true;
        }
        return false;
    }
    if (op == "not")
        return !o_eval(v, s);
    throw std::invalid_argument("o_eval: " + op);
}

/**
 * O-eval's own folding (SPEC §19.4's table), to know which names the normal form keeps (the
 * labels the pass reads and a selection_labels list may hold): unknown names (not in
 * |columns|) folded away, constants through and, or and not. Returns the folded predicate, or
 * a constant (|value|) as a JSON boolean.
 */
Json::Value o_fold(const Json::Value &p, const std::set<std::string> &columns) {
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &v = p[op];
    Json::Value out;
    if (op == "any" || op == "all" || op == "none" || op == "at_least") {
        const Json::Value &list = op == "at_least" ? v["labels"] : v;
        Json::Value known(Json::arrayValue);
        for (const Json::Value &n : list) {
            if (columns.count(n.asString()))
                known.append(n);
        }
        if (op == "any" && known.empty())
            return false;
        if (op == "all" && known.size() < list.size())
            return false;
        if (op == "none" && known.empty())
            return true;
        if (op == "at_least" && known.size() < v["n"].asUInt64())
            return false;
        if (op == "at_least") {
            out[op]["n"] = v["n"];
            out[op]["labels"] = known;
        } else {
            out[op] = known;
        }
        return out;
    }
    if (op == "not") {
        const Json::Value q = o_fold(v, columns);
        if (q.isBool())
            return !q.asBool();
        out[op] = q;
        return out;
    }
    // and, or
    const bool absorbing = op == "or";
    Json::Value kept(Json::arrayValue);
    for (const Json::Value &q : v) {
        const Json::Value f = o_fold(q, columns);
        if (f.isBool()) {
            if (f.asBool() == absorbing)
                return absorbing;
            continue;
        }
        kept.append(f);
    }
    if (kept.empty())
        return !absorbing;
    if (kept.size() == 1)
        return kept[0];
    out[op] = kept;
    return out;
}

void names_of(const Json::Value &p, std::set<std::string> *out) {
    if (p.isBool())
        return;
    const std::string op = p.getMemberNames().at(0);
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
        names_of(v, out);
    } else {
        for (const Json::Value &q : v) {
            names_of(q, out);
        }
    }
}

/**
 * O-ctx over the records: the columns whose records hold a k-mer as deposited (a BASIC row),
 * the set a context is evaluated on (with "either", or on CANONICAL and PRIMARY graphs, the
 * reverse complement's too), and the graph contexts of a pattern (BASIC: the records' k-mers).
 */
struct Oracle {
    size_t k = 0;
    std::vector<Record> records;
    std::map<std::string, std::set<std::string>> cols;

    explicit Oracle(const Index &idx) : k(idx.k), records(idx.records) {
        for (const Record &r : idx.records) {
            for (size_t s = 0; s + k <= r.seq.size(); ++s) {
                cols[r.seq.substr(s, k)].insert(r.column);
            }
        }
    }
    bool present(const std::string &kmer) const { return cols.count(kmer); }
    std::set<std::string> columns() const {
        std::set<std::string> out;
        for (const Record &r : records) {
            out.insert(r.column);
        }
        return out;
    }
    std::set<std::string> own(const std::string &kmer) const {
        auto it = cols.find(kmer);
        return it == cols.end() ? std::set<std::string>() : it->second;
    }
    std::set<std::string> evaluated(const std::string &kmer, bool either) const {
        std::set<std::string> s = own(kmer);
        if (either) {
            for (const std::string &c : own(rc(kmer))) {
                s.insert(c);
            }
        }
        return s;
    }
    // (strand, k-mer, offset) of every context of |p| (BASIC)
    std::set<std::tuple<std::string, std::string, uint64_t>>
    contexts(const std::string &p, gp::Scope scope = gp::Scope::ANY_OFFSET,
             gp::Strands strands = gp::Strands::BOTH) const {
        std::vector<std::pair<std::string, std::string>> oriented;
        const std::string q = iupac_rc(p);
        if (q == p) {
            oriented = { { "=", p } };
        } else {
            if (strands != gp::Strands::REVERSE)
                oriented.emplace_back("+", p);
            if (strands != gp::Strands::FORWARD)
                oriented.emplace_back("-", q);
        }
        std::set<std::tuple<std::string, std::string, uint64_t>> out;
        for (const auto &[kmer, columns] : cols) {
            for (const auto &[strand, o] : oriented) {
                for (size_t off = 0; off + o.size() <= k; ++off) {
                    if (scope == gp::Scope::SUFFIX && off + o.size() != k)
                        continue;
                    bool match = true;
                    for (size_t i = 0; i < o.size() && match; ++i) {
                        match = iupac_matches(o[i], kmer[off + i]);
                    }
                    if (match)
                        out.emplace(strand, kmer, off);
                }
            }
        }
        return out;
    }
    // the placed occurrences (seq_id, 1-based start) of a context of |kmer| at |offset| in
    // |column|'s records (each record's k-mers as deposited)
    std::set<std::pair<uint64_t, uint64_t>> occurrences(const std::string &column,
                                                        const std::string &kmer,
                                                        uint64_t offset) const {
        std::set<std::pair<uint64_t, uint64_t>> out;
        uint64_t seq_id = 0;
        for (const Record &r : records) {
            if (r.column != column)
                continue;
            for (size_t s = 0; s + k <= r.seq.size(); ++s) {
                if (r.seq.compare(s, k, kmer) == 0)
                    out.emplace(seq_id, s + offset + 1);
            }
            ++seq_id;
        }
        return out;
    }
};

// ---------------------------------------------------------------- the route

// One request: its body's fields, and the test hooks of the route
struct Ask {
    // the predicate's JSON; empty: a request without one
    std::string predicate = "{\"any\": [\"c1\"]}";
    std::vector<std::string> patterns = { "AC" };
    // the patterns' key: dna, iupac or protein
    std::string kind = "dna";
    std::string mode = "all_or_count";
    std::string scope = "any_offset";
    std::string strands = "both";
    uint64_t max_contexts = 10'000;
    uint64_t max_predicate_contexts = 100'000;
    bool stop_at_threshold = false;
    // output.labels: none, predicate_only or all
    std::string labels = "none";
    // predicate_strands either (true) or context; not sent when |send_strands| is false (the
    // default applies)
    bool either = true;
    bool send_strands = true;
    uint64_t max_predicate_work = 100'000'000;
    // 0: not named (the cap)
    uint64_t max_steps = 0;
    bool allow_unbudgeted = false;
    // the index's record mapping given to the retrieval (placement record)
    bool records = true;
    // the account in bytes (0: the request's 256 MiB), the decode's deny, the read hook
    uint64_t max_memory_bytes = 0;
    std::function<bool(uint64_t)> deny_decode;
    std::function<void(size_t)> read_hook;
    // a virtual clock (null: the steady clock)
    std::function<Clock::time_point()> clock;
    // the caller left (PatternDelivery::set_abort; null: never asked)
    std::function<bool()> abort;
    // the server's policy on a graph without its mask (null: the default)
    std::optional<uint64_t> max_checked_entries;
};

// the server's caps, with no information floor (the tiny patterns have few bits)
PatternLimits route_limits(const Ask &r) {
    PatternLimits l;
    l.min_information_bits = 0;
    if (r.max_checked_entries)
        l.max_checked_entries = *r.max_checked_entries;
    return l;
}

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

std::string body_of(const Ask &r) {
    Json::Value b;
    Json::Value patterns(Json::arrayValue);
    for (const std::string &p : r.patterns) {
        Json::Value q;
        q[r.kind] = p;
        patterns.append(q);
    }
    b["patterns"] = patterns;
    b["mode"] = r.mode;
    b["scope"] = r.scope;
    b["strands"] = r.strands;
    b["max_contexts"] = static_cast<Json::UInt64>(r.max_contexts);
    b["stop_at_threshold"] = r.stop_at_threshold;
    if (r.max_steps)
        b["max_steps"] = static_cast<Json::UInt64>(r.max_steps);
    if (!r.predicate.empty()) {
        b["predicate"] = parse_pattern_body(r.predicate);
        b["max_predicate_contexts"] = static_cast<Json::UInt64>(r.max_predicate_contexts);
        b["max_predicate_work"] = static_cast<Json::UInt64>(r.max_predicate_work);
        if (r.send_strands)
            b["predicate_strands"] = r.either ? "either" : "context";
    }
    b["output"]["labels"] = r.labels;
    if (r.allow_unbudgeted)
        b["allow_unbudgeted_annotation"] = true;
    return compact(b);
}

// The route's answer to |r| on |idx|
Json::Value answer_of(const Index &idx, const Ask &r) {
    RetrievalHooks hooks;
    if (r.records)
        hooks.coord_to_header = idx.cth.get();
    hooks.max_memory_bytes = r.max_memory_bytes;
    hooks.deny_decode = r.deny_decode;
    hooks.read_hook = r.read_hook;
    PatternDelivery delivery;
    if (r.abort)
        delivery.set_abort(r.abort);
    return process_pattern_request(parse_pattern_body(body_of(r)), *idx.anno, route_limits(r),
                                   "rel", nullptr, r.abort ? &delivery : nullptr, r.clock,
                                   &hooks);
}

struct Cnt {
    std::string relation;
    uint64_t value = 0;
    uint64_t lower = 0;
    uint64_t upper = 0;
};

Cnt cnt(const Json::Value &c) {
    Cnt out;
    out.relation = c["relation"].asString();
    out.value = c["value"].isNull() ? 0 : c["value"].asUInt64();
    out.lower = c.isMember("lower") ? c["lower"].asUInt64() : out.value;
    out.upper = c.isMember("upper") ? c["upper"].asUInt64() : out.value;
    return out;
}

bool contains(const Cnt &c, uint64_t truth) {
    if (c.relation == "exact")
        return c.value == truth;
    if (c.relation == "at_least")
        return c.value <= truth;
    if (c.relation == "bounds")
        return c.lower <= truth && truth <= c.upper;
    return c.relation == "unknown";
}

// one context of a pattern as an answer names it: its strand (+, -, =; on CANONICAL and
// PRIMARY graphs its orientation written so), k-mer, offset, node, and the result's JSON
struct Ctx {
    std::string strand;
    std::string kmer;
    uint64_t offset = 0;
    uint64_t node = 0;
    Json::Value json;
};

Ctx ctx_of(const Json::Value &r) {
    Ctx c;
    if (r.isMember("strand")) {
        c.strand = r["strand"].asString();
    } else {
        const std::string o = r["orientation"].asString();
        c.strand = o == "forward" ? "+" : o == "reverse" ? "-" : "=";
    }
    c.kmer = r["kmer"].asString();
    c.offset = r["offset"].asUInt64();
    c.node = r["node"].asUInt64();
    c.json = r;
    return c;
}

uint8_t strand_rank(const std::string &s) {
    return s == "+" ? 0 : s == "-" ? 1 : 2;
}

// What the answer says of one pattern's selection
struct Outcome {
    Json::Value answer;
    Json::Value entry;
    // the predicate block: a constant normal form, the binding stopped (normal_form null)
    bool constant = false;
    bool bound = true;
    std::vector<std::string> unknown;
    std::string strands;
    // the entry's
    std::string pass, access, support, withheld, cut, determinism;
    Cnt raw, tested, selected;
    std::vector<Ctx> results;
    // per result, when the results carry them
    std::vector<std::vector<std::string>> selection_labels;
    // beside them, per label the orientation whose row carries it (the owner's answer to P11)
    std::vector<std::vector<std::string>> selection_strands;
    bool has_selection_labels = false;
    bool has_selection_strands = false;
    uint64_t rows = 0, units = 0, lookups = 0, memory = 0;
    std::optional<std::pair<std::string, std::string>> stop;
    std::vector<std::string> notes;
};

Outcome outcome(const Json::Value &answer, size_t i) {
    Outcome o;
    o.answer = answer;
    o.entry = answer["patterns"][static_cast<Json::ArrayIndex>(i)];
    const Json::Value &p = answer["predicate"];
    if (p.isObject()) {
        o.constant = p["normal_form"].isBool();
        o.bound = !p["normal_form"].isNull();
        for (const Json::Value &u : p["unknown_labels"]) {
            o.unknown.push_back(u.asString());
        }
        o.strands = p["strands"].asString();
    }
    const Json::Value &e = o.entry;
    if (e.isMember("error"))
        return o;
    o.pass = e["selection"]["pass"].asString();
    o.access = e["selection"]["access"].asString();
    o.support = e["selection"]["support"].asString();
    o.withheld = e["withheld"].isObject() ? e["withheld"]["reason"].asString() : "";
    o.cut = e["cut"].isObject() ? e["cut"]["reason"].asString() : "";
    o.determinism = e["determinism"].asString();
    o.raw = cnt(e["counts"]["contexts"]);
    if (e["counts"].isMember("tested")) {
        o.tested = cnt(e["counts"]["tested"]);
        o.selected = cnt(e["counts"]["selected"]);
    }
    for (const Json::Value &r : e["results"]) {
        o.results.push_back(ctx_of(r));
        if (r.isMember("selection_labels")) {
            o.has_selection_labels = true;
            std::vector<std::string> names;
            for (const Json::Value &n : r["selection_labels"]) {
                names.push_back(n.asString());
            }
            o.selection_labels.push_back(std::move(names));
        }
        if (r.isMember("selection_strands")) {
            o.has_selection_strands = true;
            std::vector<std::string> strands;
            for (const Json::Value &n : r["selection_strands"]) {
                strands.push_back(n.asString());
            }
            o.selection_strands.push_back(std::move(strands));
        }
    }
    o.rows = e["work"]["predicate_rows"].asUInt64();
    o.units = e["work"]["predicate_units"].asUInt64();
    o.lookups = e["work"]["predicate_lookups"].asUInt64();
    o.memory = e["work"]["memory_bytes"].asUInt64();
    if (e["stop"].isObject())
        o.stop = std::make_pair(e["stop"]["phase"].asString(), e["stop"]["reason"].asString());
    for (const Json::Value &n : e["notes"]) {
        o.notes.push_back(n.asString());
    }
    return o;
}

std::vector<Outcome> run_all(const Index &idx, const Ask &r) {
    const Json::Value answer = answer_of(idx, r);
    std::vector<Outcome> out;
    for (Json::ArrayIndex i = 0; i < answer["patterns"].size(); ++i) {
        out.push_back(outcome(answer, i));
    }
    return out;
}

Outcome run(const Index &idx, const Ask &r) {
    return run_all(idx, r).at(0);
}

// The raw contexts of |r|'s first pattern in answer order: the route's release of the pattern
// without a predicate (the engine's, labels none), the order the pass reads its rows in
std::vector<Ctx> raw_order(const Index &idx, Ask r) {
    r.predicate.clear();
    r.labels = "none";
    r.mode = "all_or_count";
    r.max_contexts = 10'000;
    r.stop_at_threshold = false;
    r.max_steps = 0;
    r.max_memory_bytes = 0;
    r.deny_decode = nullptr;
    r.read_hook = nullptr;
    r.clock = nullptr;
    r.abort = nullptr;
    const Outcome o = run(idx, r);
    EXPECT_EQ("exact", o.raw.relation);
    EXPECT_EQ(o.raw.value, o.results.size());
    return o.results;
}

std::set<std::tuple<std::string, std::string, uint64_t>> ctx_set(const std::vector<Ctx> &v) {
    std::set<std::tuple<std::string, std::string, uint64_t>> s;
    for (const Ctx &c : v) {
        s.emplace(c.strand, c.kmer, c.offset);
    }
    return s;
}

// the contexts of |raw| (in its order) the oracle selects
std::vector<Ctx> oracle_selected(const Oracle &oracle, const Json::Value &pred,
                                 const std::vector<Ctx> &raw, bool either) {
    std::vector<Ctx> out;
    for (const Ctx &c : raw) {
        if (o_eval(pred, oracle.evaluated(c.kmer, either)))
            out.push_back(c);
    }
    return out;
}

std::vector<std::tuple<std::string, std::string, uint64_t>> keys(const std::vector<Ctx> &v) {
    std::vector<std::tuple<std::string, std::string, uint64_t>> out;
    for (const Ctx &c : v) {
        out.emplace_back(c.strand, c.kmer, c.offset);
    }
    return out;
}

/**
 * The invariants of TESTS §1 on one answered entry: every result selected by the oracle, in
 * answer order (node, offset, orientation); tested exact once the pass ran; selected contains
 * |truth| (the oracle's selected among all raw contexts); a completed pass tested every raw
 * context and lists the oracle's selected ones (the first max_contexts in partial; all of them
 * or none in all_or_count); a withheld answer and mode count list nothing; the results'
 * selection_labels are the set they were evaluated on (restricted to the normal form's names),
 * satisfy the predicate and come in label order (contexts desc, column asc).
 */
void check_invariants(const Oracle &oracle, const Json::Value &pred, const Outcome &o,
                      bool either, std::optional<uint64_t> truth, const Ask &r,
                      const std::vector<Ctx> *raw = nullptr) {
    ASSERT_FALSE(o.entry.isMember("error")) << compact(o.entry["error"]);
    EXPECT_EQ("predicate", o.entry["absence_filter"].asString());
    for (const Ctx &c : o.results) {
        EXPECT_TRUE(o_eval(pred, oracle.evaluated(c.kmer, either)))
                << c.kmer << " " << c.strand << c.offset << ": not selected by the oracle";
    }
    for (size_t j = 1; j < o.results.size(); ++j) {
        const Ctx &a = o.results[j - 1], &b = o.results[j];
        EXPECT_TRUE(std::make_tuple(a.node, a.offset, strand_rank(a.strand))
                        < std::make_tuple(b.node, b.offset, strand_rank(b.strand)))
                << "answer order";
    }
    if (o.pass == "completed" || o.pass == "stopped") {
        EXPECT_EQ("exact", o.tested.relation);
    }
    if (truth) {
        EXPECT_TRUE(contains(o.selected, *truth)) << *truth << " " << compact(o.entry["counts"]);
    }
    if (!o.withheld.empty() || r.mode == "count") {
        EXPECT_TRUE(o.results.empty());
    }
    if (r.mode == "count") {
        EXPECT_FALSE(o.entry.isMember("results"));
        EXPECT_FALSE(o.entry["retrieval_complete"].asBool());
    }
    if (o.pass == "completed") {
        EXPECT_EQ("exact", o.selected.relation);
        // after a completed pass only the output and the projection's reads can stop
        EXPECT_TRUE(!o.stop || o.stop->first == "output" || o.stop->first == "placement"
                    || o.stop->first == "label_discovery") << o.stop->first;
        if (raw) {
            EXPECT_EQ(raw->size(), o.tested.value);
            const std::vector<Ctx> want = oracle_selected(oracle, pred, *raw, either);
            EXPECT_EQ(want.size(), o.selected.value);
            if (r.mode != "count" && o.withheld.empty() && !o.stop) {
                std::vector<Ctx> first(want.begin(),
                                       want.begin() + std::min<uint64_t>(want.size(),
                                                                         r.max_contexts));
                EXPECT_EQ(keys(first), keys(o.results)) << "the first selected, in order";
            }
        }
    }
    if (r.mode == "all_or_count" && o.pass == "completed" && truth) {
        if (*truth > r.max_contexts) {
            EXPECT_EQ("selected_above_threshold", o.withheld);
        } else if (!o.stop) {
            EXPECT_EQ("", o.withheld);
            EXPECT_EQ(*truth, o.results.size());
            EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
        }
    }
    // selection_labels: with a projection that reads labels only
    if (r.labels == "none" || o.results.empty()) {
        EXPECT_FALSE(o.has_selection_labels);
        EXPECT_FALSE(o.has_selection_strands);
        return;
    }
    std::set<std::string> names;
    names_of(o_fold(pred, oracle.columns()), &names);
    ASSERT_EQ(o.results.size(), o.selection_labels.size());
    // selection_strands (the owner's answer to P11): per label the orientation whose row
    // carries it, from the records: on a BASIC graph "context" when the context's k-mer x
    // holds it as deposited, "reverse_complement" when only rc(x) does (with "either"),
    // "both" when both do (x == rc(x) included); "either" on CANONICAL and PRIMARY graphs (one
    // row for x and rc(x)). Mutations tried: the flags swapped, every label "context", the
    // mirror's flag lost in the merge, a palindrome "context": each fails here
    ASSERT_EQ(o.results.size(), o.selection_strands.size());
    for (size_t j = 0; j < o.results.size(); ++j) {
        ASSERT_EQ(o.selection_labels[j].size(), o.selection_strands[j].size());
        const Ctx &c = o.results[j];
        const bool basic = c.json.isMember("strand");
        const std::set<std::string> on_x = oracle.own(c.kmer);
        const std::set<std::string> on_rc = oracle.own(rc(c.kmer));
        for (size_t t = 0; t < o.selection_labels[j].size(); ++t) {
            const std::string &label = o.selection_labels[j][t];
            std::string want = "either";
            if (basic) {
                const bool x = on_x.count(label);
                const bool y = either && on_rc.count(label);
                want = x && y ? "both" : y ? "reverse_complement" : "context";
            }
            EXPECT_EQ(want, o.selection_strands[j][t]) << c.kmer << " " << label;
        }
    }
    std::map<std::string, uint64_t> count;
    for (size_t j = 0; j < o.results.size(); ++j) {
        std::set<std::string> expected;
        for (const std::string &c : oracle.evaluated(o.results[j].kmer, either)) {
            if (names.count(c))
                expected.insert(c);
        }
        const std::set<std::string> got(o.selection_labels[j].begin(),
                                        o.selection_labels[j].end());
        EXPECT_EQ(expected, got) << o.results[j].kmer;
        EXPECT_EQ(got.size(), o.selection_labels[j].size());
        EXPECT_TRUE(o_eval(pred, got)) << "the selection_labels satisfy the predicate";
        for (const std::string &c : got) {
            count[c]++;
        }
    }
    for (const auto &list : o.selection_labels) {
        for (size_t t = 1; t < list.size(); ++t) {
            const auto &a = list[t - 1], &b = list[t];
            EXPECT_TRUE(count[a] > count[b] || (count[a] == count[b] && a < b))
                    << a << " before " << b << ": contexts desc, column asc";
        }
    }
}

uint64_t oracle_truth(const Oracle &oracle, const Json::Value &pred, const std::string &p,
                      bool either, gp::Scope scope = gp::Scope::ANY_OFFSET) {
    uint64_t n = 0;
    for (const auto &[strand, kmer, offset] : oracle.contexts(p, scope)) {
        n += o_eval(pred, oracle.evaluated(kmer, either));
    }
    return n;
}

// a random predicate over |universe| (depth <= 3, lists of 1-3 distinct names)
Json::Value random_predicate(std::mt19937 &rng, const std::vector<std::string> &universe,
                             int depth = 0) {
    auto list = [&]() {
        std::vector<std::string> u = universe;
        std::shuffle(u.begin(), u.end(), rng);
        Json::Value l(Json::arrayValue);
        const size_t n = 1 + rng() % 3;
        for (size_t i = 0; i < n && i < u.size(); ++i) {
            l.append(u[i]);
        }
        return l;
    };
    Json::Value p;
    const unsigned op = depth >= 3 ? rng() % 4 : rng() % 7;
    switch (op) {
        case 0: p["any"] = list(); break;
        case 1: p["all"] = list(); break;
        case 2: p["none"] = list(); break;
        case 3: {
            Json::Value a;
            a["labels"] = list();
            a["n"] = static_cast<Json::UInt64>(1 + rng() % a["labels"].size());
            p["at_least"] = a;
            break;
        }
        case 4:
        case 5: {
            Json::Value l(Json::arrayValue);
            const size_t n = 1 + rng() % 3;
            for (size_t i = 0; i < n; ++i) {
                l.append(random_predicate(rng, universe, depth + 1));
            }
            p[op == 4 ? "and" : "or"] = l;
            break;
        }
        default: p["not"] = random_predicate(rng, universe, depth + 1); break;
    }
    return p;
}

// the rows of the pass in their order of first appearance, from the raw contexts in answer
// order (BASIC): a context's own k-mer, then (either) its reverse complement's when the graph
// holds it
std::vector<std::string> row_order(const Oracle &oracle, const std::vector<Ctx> &raw,
                                   bool either) {
    std::vector<std::string> order;
    std::set<std::string> seen;
    for (const Ctx &c : raw) {
        if (seen.insert(c.kmer).second)
            order.push_back(c.kmer);
        const std::string y = rc(c.kmer);
        if (either && y != c.kmer && oracle.present(y) && seen.insert(y).second)
            order.push_back(y);
    }
    return order;
}

// the raw contexts all of whose rows are among |read| (decided by a pass that read them)
std::vector<Ctx> ready_among(const Oracle &oracle, const std::vector<Ctx> &raw,
                             const std::set<std::string> &read, bool either) {
    std::vector<Ctx> out;
    for (const Ctx &c : raw) {
        const std::string y = rc(c.kmer);
        if (read.count(c.kmer)
                && (!either || y == c.kmer || !oracle.present(y) || read.count(y))) {
            out.push_back(c);
        }
    }
    return out;
}

// |answer| without its timing (the only part of an answer that may differ between two runs)
Json::Value without_timing(Json::Value answer) {
    answer.removeMember("timing");
    for (Json::Value &e : answer["patterns"]) {
        e.removeMember("timing");
    }
    return answer;
}

bool has_note(const Outcome &o, const std::string &note) {
    return std::find(o.notes.begin(), o.notes.end(), note) != o.notes.end();
}


// ---------------------------------------------------------------- tests

// Mutations tried (5b-3 mut/, through the route since 5b-4): "either" ignored, the mirror's
// labels left out of the evaluated set, the label order by name only, rows not deduplicated by
// key, the kept rows dropped, selection_labels of the own row only: each fails here.
TEST(PatternSelection, PerRowAgainstTheOracle) {
    for (const bool records_rowdiff : { true, false }) {
        Index idx = records_rowdiff ? build<annot::RowDiffColumnAnnotator>(kK, kR1, true)
                                    : build<annot::RowDiffColumnAnnotator>(7, random_records(12, 7),
                                                                           true);
        Oracle oracle(idx);
        std::vector<std::string> universe;
        for (const Record &rec : idx.records) {
            universe.push_back(rec.column);
        }
        std::sort(universe.begin(), universe.end());
        universe.erase(std::unique(universe.begin(), universe.end()), universe.end());
        universe.push_back("zz");
        universe.push_back("c1x");
        std::mt19937 rng(records_rowdiff ? 11 : 13);
        const auto truth_contexts = oracle.contexts("AC");
        Ask base;
        const std::vector<Ctx> raw = raw_order(idx, base);
        ASSERT_EQ(truth_contexts, ctx_set(raw));
        for (int t = 0; t < 50; ++t) {
            Json::Value pred = random_predicate(rng, universe);
            for (const bool either : { true, false }) {
                for (const std::string labels : { "none", "predicate_only", "all" }) {
                    Ask r;
                    r.predicate = compact(pred);
                    r.either = either;
                    r.labels = labels;
                    Outcome o = run(idx, r);
                    SCOPED_TRACE(r.predicate + (either ? " either " : " context ") + labels);
                    const Json::Value folded = o_fold(pred, oracle.columns());
                    EXPECT_EQ(folded.isBool(), o.constant) << compact(folded);
                    // the echo: the normal form as the oracle folds it
                    EXPECT_EQ(compact(folded), compact(o.answer["predicate"]["normal_form"]));
                    if (o.constant) {
                        // folded to a constant (unknown names): no pass, no read
                        const bool value = o_eval(pred, {});
                        EXPECT_EQ("constant", o.pass);
                        EXPECT_EQ(0u, o.rows);
                        EXPECT_TRUE(has_note(o, "predicate_constant"));
                        bool all_same = true;
                        for (const auto &[s, kmer, off] : truth_contexts) {
                            all_same &= o_eval(pred, oracle.evaluated(kmer, either)) == value;
                        }
                        EXPECT_TRUE(all_same) << "a constant decides every context alike";
                        EXPECT_EQ(value ? truth_contexts.size() : 0u, o.results.size());
                        continue;
                    }
                    EXPECT_EQ("completed", o.pass);
                    check_invariants(oracle, pred, o, either,
                                     oracle_truth(oracle, pred, "AC", either), r, &raw);
                    EXPECT_EQ(either ? "either" : "context", o.strands);
                    EXPECT_EQ("rows", o.access);
                    EXPECT_EQ("kmer", o.support);
                    EXPECT_EQ(row_order(oracle, raw, either).size(), o.rows);
                    if (labels == "none") {
                        EXPECT_FALSE(o.entry.isMember("by_label"));
                        for (const Ctx &c : o.results) {
                            EXPECT_FALSE(c.json.isMember("labels"));
                        }
                        continue;
                    }
                    EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
                    std::set<std::string> names;
                    names_of(folded, &names);
                    for (const Ctx &c : o.results) {
                        std::set<std::string> expected;
                        for (const std::string &col : oracle.own(c.kmer)) {
                            if (labels == "all" || names.count(col))
                                expected.insert(col);
                        }
                        std::set<std::string> got;
                        for (const Json::Value &l : c.json["labels"]) {
                            const std::string column = l["column"].asString();
                            got.insert(column);
                            // placed as the record scan says
                            std::set<std::pair<uint64_t, uint64_t>> occ;
                            for (const Json::Value &e : l["occurrence_list"]) {
                                const std::string nt = e["nt_coords"].asString();
                                occ.emplace(e["seq_id"].asUInt64(),
                                            std::stoull(nt.substr(0, nt.find('-'))));
                                EXPECT_EQ(c.strand, e["strand"].asString());
                            }
                            EXPECT_EQ(oracle.occurrences(column, c.kmer, c.offset), occ)
                                    << column << " " << c.kmer;
                        }
                        EXPECT_EQ(expected, got) << c.kmer;
                        EXPECT_EQ(expected.size(), c.json["labels_total"].asUInt64());
                        EXPECT_EQ("complete", c.json["labels_status"].asString());
                    }
                }
            }
        }
    }
}

// The determinism of the route with a predicate: the same request twice, the same answer apart
// from timing (mutation tried: the kept rows dropped, retrieve_given refuses: fails)
TEST(PatternSelection, Deterministic) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Ask r;
    r.predicate = "{\"or\": [{\"any\": [\"c1\", \"c5\"]}, {\"none\": [\"c3\"]}]}";
    r.labels = "predicate_only";
    r.mode = "partial";
    r.max_contexts = 3;
    r.patterns = { "AC", "GA" };
    const Json::Value a = answer_of(idx, r), b = answer_of(idx, r);
    EXPECT_EQ(compact(without_timing(a)), compact(without_timing(b)));
    ASSERT_EQ(3u, a["patterns"][0]["results"].size());
    EXPECT_EQ("max_contexts", a["patterns"][0]["cut"]["reason"].asString());
}

// P11: c5 holds GGACGACTTTG's motif on the - strand only. Mutation tried: the mirror's row
// not read ("either" as "context"): fails.
TEST(PatternSelection, StrandsEitherAndContextOnBasic) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    // none(c5): "context" selects the + contexts of c1/r0's k-mers (their own rows lack c5),
    // "either" none of them (their reverse complements are c5's)
    for (const bool either : { false, true }) {
        Ask r;
        r.predicate = "{\"none\": [\"c5\"]}";
        r.either = either;
        Outcome o = run(idx, r);
        const Json::Value pred = parse_pattern_body(r.predicate);
        check_invariants(oracle, pred, o, either, oracle_truth(oracle, pred, "AC", either), r,
                         &raw);
        EXPECT_TRUE(o.answer["predicate"]["vacuous"].asBool());
        auto c1_plus = [&](const Ctx &c) {
            return c.strand == "+" && oracle.own(c.kmer).count("c1")
                    && oracle.own(rc(c.kmer)).count("c5");
        };
        const uint64_t all = std::count_if(raw.begin(), raw.end(), c1_plus);
        const uint64_t selected = std::count_if(o.results.begin(), o.results.end(), c1_plus);
        ASSERT_GT(all, 0u);
        EXPECT_EQ(either ? 0u : all, selected);
    }
    // the default of predicate_strands is "either" (P11): the same selection, echoed as such
    {
        Ask r;
        r.predicate = "{\"none\": [\"c5\"]}";
        Ask d = r;
        d.send_strands = false;
        const Outcome a = run(idx, r), b = run(idx, d);
        EXPECT_EQ(keys(a.results), keys(b.results));
        EXPECT_EQ("either", b.strands);
        EXPECT_EQ("either", b.answer["limits"]["predicate_strands"].asString());
        EXPECT_EQ(a.lookups, b.lookups);
    }
    // any(c5) with "either": the + contexts of c1/r0 are selected by their mirror's label,
    // which their selection_labels show and their own labels (predicate_only) do not
    Ask r;
    r.predicate = "{\"any\": [\"c5\"]}";
    r.labels = "predicate_only";
    Outcome o = run(idx, r);
    const Json::Value pred = parse_pattern_body(r.predicate);
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r, &raw);
    bool seen = false;
    for (size_t j = 0; j < o.results.size(); ++j) {
        if (oracle.own(o.results[j].kmer).count("c5"))
            continue;
        seen = true;
        EXPECT_EQ(std::vector<std::string>({ "c5" }), o.selection_labels[j]);
        // the owner's answer to P11: the answer says which orientation supported it
        EXPECT_EQ(std::vector<std::string>({ "reverse_complement" }), o.selection_strands[j]);
        EXPECT_EQ(0u, o.results[j].json["labels"].size()) << o.results[j].kmer;
    }
    EXPECT_TRUE(seen);
    // and every kind of row a label can be found on, over R1's contexts of AC
    std::set<std::string> kinds;
    for (const auto &strands : o.selection_strands) {
        kinds.insert(strands.begin(), strands.end());
    }
    EXPECT_TRUE(kinds.count("reverse_complement"));
    EXPECT_TRUE(kinds.count("context") || kinds.count("both"));
    // with "context" every label is the context's own
    r.either = false;
    r.predicate = "{\"any\": [\"c1\", \"c5\"]}";
    o = run(idx, r);
    check_invariants(oracle, parse_pattern_body(r.predicate), o, false, std::nullopt, r, &raw);
    ASSERT_FALSE(o.results.empty());
    for (const auto &strands : o.selection_strands) {
        for (const std::string &strand : strands) {
            EXPECT_EQ("context", strand);
        }
    }
}

// On CANONICAL and PRIMARY graphs one row serves both orientations: "either" whatever is
// asked, no lookup; with "predicate_only" each result's labels are the predicate's labels on
// that shared row (retrieve_given finds every returned context's row among the pass's kept
// ones: the route's key and the pass's are computed alike; the review of 5b, L4) and its
// selection_strands are "either". Mutations tried: the lookups run on every graph mode; the
// pass keying a CANONICAL context by its own node instead of the canonical k-mer's: fail.
TEST(PatternSelection, CanonicalAndPrimaryShareTheRow) {
    for (const bool primary : { false, true }) {
        Index idx = primary
                ? build<annot::ColumnCompressed<>>(kK, kR1, false, DeBruijnGraph::PRIMARY)
                : build<annot::RowDiffColumnAnnotator>(kK, kR1, false, DeBruijnGraph::CANONICAL);
        Oracle oracle(idx);
        Ask base;
        base.records = false;
        const std::vector<Ctx> raw = raw_order(idx, base);
        for (const std::string &predicate : { std::string("{\"none\": [\"c5\"]}"),
                                              std::string("{\"all\": [\"c1\", \"c5\"]}"),
                                              std::string("{\"any\": [\"c2\", \"c3\"]}") }) {
            for (const bool either : { false, true }) {
                Ask r;
                r.predicate = predicate;
                r.either = either;
                r.records = false;
                r.allow_unbudgeted = primary;
                Outcome o = run(idx, r);
                SCOPED_TRACE(predicate + (primary ? " primary" : " canonical"));
                EXPECT_EQ("either", o.strands);
                EXPECT_EQ(0u, o.lookups);
                EXPECT_EQ("completed", o.pass);
                ASSERT_GT(o.tested.value, 0u);
                check_invariants(oracle, parse_pattern_body(predicate), o, true, std::nullopt,
                                 r, &raw);
                EXPECT_EQ(primary ? "columns" : "rows", o.access);
                // the request echoes what it asked for, the block what was evaluated
                EXPECT_EQ(either ? "either" : "context",
                          o.answer["limits"]["predicate_strands"].asString());

                // predicate_only on the shared row
                r.labels = "predicate_only";
                const Outcome q = run(idx, r);
                check_invariants(oracle, parse_pattern_body(predicate), q, true, std::nullopt,
                                 r, &raw);
                EXPECT_EQ(keys(o.results), keys(q.results));
                std::set<std::string> names;
                names_of(o_fold(parse_pattern_body(predicate), oracle.columns()), &names);
                for (const Ctx &c : q.results) {
                    std::set<std::string> want;
                    for (const std::string &column : oracle.evaluated(c.kmer, true)) {
                        if (names.count(column))
                            want.insert(column);
                    }
                    std::set<std::string> got;
                    for (const Json::Value &l : c.json["labels"]) {
                        got.insert(l["column"].asString());
                    }
                    EXPECT_EQ(want, got) << c.kmer;
                    EXPECT_EQ(want.size(), c.json["labels_total"].asUInt64()) << c.kmer;
                }
            }
        }
    }
}

// NA: two contexts at one (k-mer, offset), one row read for both, decided alike. Mutation
// tried: rows not deduplicated by key: predicate_rows exceeds the distinct rows, fails.
TEST(PatternSelection, IupacTwoOrientationsOneRow) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    Ask base;
    base.patterns = { "NA" };
    base.kind = "iupac";
    const std::vector<Ctx> raw = raw_order(idx, base);
    EXPECT_EQ(oracle.contexts("NA"), ctx_set(raw));
    bool shared = false;
    for (size_t i = 0; i + 1 < raw.size(); ++i) {
        shared |= raw[i].kmer == raw[i + 1].kmer && raw[i].offset == raw[i + 1].offset;
    }
    EXPECT_TRUE(shared) << "NA has two contexts at one (k-mer, offset) here";
    for (const bool either : { false, true }) {
        Ask r = base;
        r.predicate = "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"c4\", \"c5\"]}}";
        r.either = either;
        Outcome o = run(idx, r);
        const Json::Value pred = parse_pattern_body(r.predicate);
        check_invariants(oracle, pred, o, either, std::nullopt, r, &raw);
        EXPECT_EQ(row_order(oracle, raw, either).size(), o.rows);
        std::map<std::pair<std::string, uint64_t>, std::set<std::string>> listed;
        for (const Ctx &c : o.results) {
            listed[{ c.kmer, c.offset }].insert(c.strand);
        }
        for (size_t i = 0; i + 1 < raw.size(); ++i) {
            if (raw[i].kmer != raw[i + 1].kmer || raw[i].offset != raw[i + 1].offset)
                continue;
            // both orientations listed, or neither
            const auto it = listed.find({ raw[i].kmer, raw[i].offset });
            EXPECT_TRUE(it == listed.end() || it->second.size() == 2) << raw[i].kmer;
        }
    }
}

// at_least(n, cohort) for every n, on R1 and on a 30-column random index; n = |cohort| is all
// (mutation tried: "either" ignored: fails)
TEST(PatternSelection, AtLeastOnCohorts) {
    for (const bool random : { false, true }) {
        Index idx = random ? build<annot::RowDiffColumnAnnotator>(7, random_records(30, 5), true)
                           : build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
        Oracle oracle(idx);
        const std::vector<Ctx> raw = raw_order(idx, Ask());
        std::vector<std::string> cohort;
        for (const Record &rec : idx.records) {
            if (cohort.empty() || cohort.back() != rec.column)
                cohort.push_back(rec.column);
        }
        if (cohort.size() > 12)
            cohort.resize(12);
        Json::Value labels(Json::arrayValue), all;
        for (const std::string &c : cohort) {
            labels.append(c);
        }
        all["all"] = labels;
        for (size_t n = 1; n <= cohort.size(); ++n) {
            Json::Value pred;
            pred["at_least"]["n"] = static_cast<Json::UInt64>(n);
            pred["at_least"]["labels"] = labels;
            Ask r;
            r.predicate = compact(pred);
            Outcome o = run(idx, r);
            SCOPED_TRACE(r.predicate);
            check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r,
                             &raw);
            if (n == cohort.size()) {
                Ask q = r;
                q.predicate = compact(all);
                Outcome a = run(idx, q);
                EXPECT_EQ(keys(a.results), keys(o.results));
                EXPECT_EQ(a.selected.value, o.selected.value);
            }
        }
    }
}

// A typo is an unknown name, folded away: and(any(c1), none(c1x)) selects as any(c1) and as
// the oracle says, the echo naming the typo (mutations tried: "either" ignored, the mirror's
// labels dropped: fail)
TEST(PatternSelection, TypoReportedUnknown) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Ask r;
    r.predicate = "{\"and\": [{\"any\": [\"c1\"]}, {\"none\": [\"c1x\"]}]}";
    Outcome o = run(idx, r);
    EXPECT_EQ(std::vector<std::string>({ "c1x" }), o.unknown);
    EXPECT_EQ("{\"any\":[\"c1\"]}", compact(o.answer["predicate"]["normal_form"]));
    EXPECT_EQ(2u, o.answer["predicate"]["names"].asUInt64());
    EXPECT_EQ(1u, o.answer["predicate"]["known"].asUInt64());
    EXPECT_FALSE(o.answer["predicate"]["vacuous"].asBool());
    Ask q = r;
    q.predicate = "{\"any\": [\"c1\"]}";
    Outcome p = run(idx, q);
    EXPECT_EQ(keys(p.results), keys(o.results));
    EXPECT_EQ(p.selected.value, o.selected.value);
    EXPECT_EQ("completed", o.pass);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body(r.predicate);
    const std::vector<Ctx> raw = raw_order(idx, r);
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r, &raw);
}

// P18: a predicate folded to a constant reads nothing (the read hook never called): pass
// constant, tested the raw count, selected 0 (false) or the raw count (true, the unfiltered
// answer); a constant false after a discovery stop still selects exactly 0. select() refuses a
// constant. Mutation tried: constant true answered selected exact 0: fails.
TEST(PatternSelection, AllUnknownNeedsNoRead) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const uint64_t raw = oracle.contexts("AC").size();
    for (const auto &[predicate, value] : std::vector<std::pair<std::string, bool>>{
            { "{\"any\": [\"zz\"]}", false }, { "{\"none\": [\"zz\"]}", true },
            { "{\"all\": [\"c1\", \"zz\"]}", false },
            { "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"zz\"]}}", false } }) {
        for (const std::string mode : { "all_or_count", "partial", "count" }) {
            for (const std::string labels : { "none", "predicate_only", "all" }) {
                auto reads = std::make_shared<uint64_t>(0);
                Ask r;
                r.predicate = predicate;
                r.mode = mode;
                r.labels = labels;
                r.read_hook = [reads](size_t) { ++*reads; };
                Outcome o = run(idx, r);
                SCOPED_TRACE(predicate + " " + mode + " " + labels);
                EXPECT_TRUE(o.constant);
                EXPECT_EQ(value, o.answer["predicate"]["normal_form"].asBool());
                EXPECT_EQ("constant", o.pass);
                EXPECT_TRUE(has_note(o, "predicate_constant"));
                // nothing read for the selection; "all" on a constant true reads the labels
                // of the unfiltered release as without a predicate
                if (!(value && labels == "all" && mode != "count")) {
                    EXPECT_EQ(0u, *reads);
                }
                EXPECT_EQ(0u, o.rows);
                EXPECT_EQ("exact", o.tested.relation);
                EXPECT_EQ(raw, o.tested.value);
                EXPECT_EQ("exact", o.selected.relation);
                EXPECT_EQ(value ? raw : 0u, o.selected.value);
                if (mode == "count")
                    continue;
                EXPECT_EQ(value ? raw : 0u, o.results.size());
                EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
                if (labels != "none") {
                    // nothing of the predicate on any context: empty lists, exact zeros
                    for (const Ctx &c : o.results) {
                        EXPECT_EQ(0u, c.json["selection_labels"].size());
                        if (labels == "predicate_only") {
                            EXPECT_TRUE(c.json["labels"].isArray());
                            EXPECT_EQ(0u, c.json["labels"].size());
                            EXPECT_TRUE(c.json["labels_total"].isUInt64());
                            EXPECT_EQ(0u, c.json["labels_total"].asUInt64());
                            EXPECT_EQ("complete", c.json["labels_status"].asString());
                        }
                    }
                    EXPECT_TRUE(o.entry["by_label"].isArray());
                    if (labels == "predicate_only" || !value) {
                        EXPECT_EQ("exact", o.entry["counts"]["labels"]["relation"].asString());
                        EXPECT_EQ(0u, o.entry["counts"]["labels"]["value"].asUInt64());
                    }
                }
            }
        }
    }
    // a constant false after a discovery stop: selected exact 0, the raw count at_least
    Ask r;
    r.predicate = "{\"any\": [\"zz\"]}";
    r.max_steps = 1;
    Outcome o = run(idx, r);
    EXPECT_EQ("at_least", o.raw.relation);
    EXPECT_EQ("at_least", o.tested.relation);
    EXPECT_EQ("exact", o.selected.relation);
    EXPECT_EQ(0u, o.selected.value);
    EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
    EXPECT_EQ("", o.withheld);
    ASSERT_TRUE(o.stop);
    EXPECT_EQ("max_steps", o.stop->second);
    // a constant true after it: the unfiltered answer's (discovery_budget)
    r.predicate = "{\"none\": [\"zz\"]}";
    o = run(idx, r);
    EXPECT_EQ("at_least", o.selected.relation);
    EXPECT_EQ("discovery_budget", o.withheld);

    // the API the route builds on: select() refuses a constant, the counts of the passes
    // that did not run
    {
        const gp::GraphSupport support = gp::PatternSearch::support(idx.anno->get_graph());
        gp::Budget budget(1'000'000, mtg::test::unbounded_deadline());
        PatternRetrieval retrieval(*idx.anno, support.mode, RetrievalLimits(), budget);
        retrieval.bind(predicate::Predicate::parse(parse_pattern_body("{\"any\": [\"zz\"]}"),
                                                   10'000),
                       SelectionLimits());
        std::vector<TestedContext> none;
        gp::Extraction x;
        EXPECT_THROW(retrieval.select(none, 0, gp::Count::exact(gp::Unit::GRAPH_CONTEXTS, 0), x,
                                      SelectionRequest()),
                     std::logic_error);
    }
    const gp::Count raw_stopped = gp::Count::at_least(gp::Unit::GRAPH_CONTEXTS, 3);
    SelectionAnswer t = constant_selection(true, raw_stopped);
    EXPECT_EQ(gp::Relation::AT_LEAST, t.selected.relation);
    EXPECT_EQ(3u, t.selected.value);
    SelectionAnswer n = selection_without_pass(SelectionPass::NOT_ADMITTED, raw_stopped);
    EXPECT_EQ(gp::Relation::UNKNOWN, n.selected.relation);
    n = selection_without_pass(SelectionPass::NOT_STARTED,
                               gp::Count::exact(gp::Unit::GRAPH_CONTEXTS, 0));
    EXPECT_EQ(gp::Relation::EXACT, n.selected.relation);
    EXPECT_THROW(selection_without_pass(SelectionPass::COMPLETED, raw_stopped),
                 std::logic_error);
}

/**
 * Every max_predicate_work from 1 to the pass's total + 1 (both strand settings, the three
 * modes): below the total the pass stops {selection, max_predicate_work} with selected bounds
 * [S, S + R - T], the decided contexts exactly those whose rows are among the first rows read
 * in the order of first appearance (a prefix: their number is tested, their selected ones are
 * S and, in partial, the results), more units never fewer of them; at the total + 1 it
 * completes. SPEC §19.13 (d): with "either" the first lookup comes first and no row is read at
 * a budget of 1. Mutations tried: the bounds' upper end without R - T (= S + R), the final
 * sweep removed (contexts after the cursor with read rows undecided), the gate checked after
 * the read (one row more), the refusal's units not charged: each fails here.
 */
TEST(PatternSelection, PassStoppedByItsBudget) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::string predicate = "{\"or\": [{\"all\": [\"c1\", \"c4\"]}, {\"any\": [\"c3\"]}]}";
    const Json::Value pred = parse_pattern_body(predicate);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    for (const bool either : { true, false }) {
        const std::vector<std::string> order = row_order(oracle, raw, either);
        for (const std::string mode : { "all_or_count", "partial", "count" }) {
            Ask r;
            r.predicate = predicate;
            r.either = either;
            r.mode = mode;
            const Outcome full = run(idx, r);
            ASSERT_EQ("completed", full.pass);
            const uint64_t total = full.units;
            ASSERT_GT(total, 0u);
            const uint64_t truth = oracle_truth(oracle, pred, "AC", either);
            ASSERT_GT(full.rows, 1u);
            // the budget at which exactly one row is read: every lookup (k units each) passes
            // its gate, the first read's gate passes (units < W), the second's does not
            const uint64_t one_row = either ? full.lookups * kK + 1 : 1;
            uint64_t last_tested = 0;
            bool completed_before = false;
            for (uint64_t w = 1; w <= total + 1; ++w) {
                r.max_predicate_work = w;
                const Outcome o = run(idx, r);
                SCOPED_TRACE("work " + std::to_string(w) + (either ? " either " : " context ")
                             + mode);
                check_invariants(oracle, pred, o, either, truth, r, &raw);
                EXPECT_GE(o.tested.value, last_tested) << "monotone in the budget";
                last_tested = o.tested.value;
                EXPECT_LE(o.units, total);
                if (w == total + 1) {
                    EXPECT_EQ("completed", o.pass);
                }
                if (w == one_row) {
                    EXPECT_EQ("stopped", o.pass);
                    EXPECT_EQ(1u, o.rows) << "the gate is checked before a read";
                }
                // completed at a budget, completed at every larger one
                EXPECT_TRUE(!completed_before || o.pass == "completed");
                completed_before |= o.pass == "completed";
                if (o.pass == "completed")
                    continue;
                ASSERT_EQ("stopped", o.pass);
                ASSERT_TRUE(o.stop);
                EXPECT_EQ("selection", o.stop->first);
                EXPECT_EQ("max_predicate_work", o.stop->second);
                EXPECT_EQ("full", o.determinism);
                // it stopped at a gate: the units had reached the budget
                EXPECT_GE(o.units, w);
                EXPECT_EQ("bounds", o.selected.relation);
                EXPECT_EQ(o.selected.lower + raw.size() - o.tested.value, o.selected.upper);
                // the decided contexts: those whose rows are among the first rows read
                ASSERT_LE(o.rows, order.size());
                const std::set<std::string> read(order.begin(), order.begin() + o.rows);
                const std::vector<Ctx> ready = ready_among(oracle, raw, read, either);
                EXPECT_EQ(ready.size(), o.tested.value);
                const std::vector<Ctx> sel = oracle_selected(oracle, pred, ready, either);
                EXPECT_EQ(sel.size(), o.selected.lower);
                if (mode == "all_or_count") {
                    EXPECT_EQ("predicate_budget", o.withheld);
                } else if (mode == "partial") {
                    EXPECT_EQ("max_predicate_work", o.cut);
                    EXPECT_EQ("", o.withheld);
                    EXPECT_EQ(keys(sel), keys(o.results));
                    EXPECT_FALSE(o.entry["retrieval_complete"].asBool());
                } else {
                    EXPECT_FALSE(o.entry.isMember("withheld"));
                    EXPECT_FALSE(o.entry.isMember("cut"));
                }
                if (w == 1 && either) {
                    // the first context's lookup (k units) comes first: no row read
                    EXPECT_EQ(0u, o.rows);
                    EXPECT_EQ(0u, o.tested.value);
                    EXPECT_EQ(kK, o.units);
                }
            }
        }
    }
}

// P17: a read is charged its decoded row's whole size, whatever the predicate restricts it
// to. Two predicates over the same rows: their units differ by the decisions' only (each
// context 1 + its present labels, every name in one leaf), computed by the oracle. Mutation
// tried: charging the hits (8 + hits + dependency units) instead of KeyCost::entries: fails.
TEST(PatternSelection, ReadsChargeTheWholeRow) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    std::vector<uint64_t> reads;
    for (const std::string &predicate : { std::string("{\"any\": [\"c1\"]}"),
                                          std::string("{\"any\": [\"c1\", \"c2\", \"c3\", "
                                                      "\"c4\", \"c5\"]}") }) {
        Ask r;
        r.predicate = predicate;
        r.either = false;
        Outcome o = run(idx, r);
        std::set<std::string> names;
        names_of(parse_pattern_body(predicate), &names);
        uint64_t decisions = 0;
        for (const Ctx &c : raw) {
            decisions += 1;
            for (const std::string &col : oracle.own(c.kmer)) {
                decisions += names.count(col);
            }
        }
        reads.push_back(o.units - decisions);
    }
    EXPECT_EQ(reads[0], reads[1]);
}

// One reverse-complement lookup per distinct row, its result kept, and none for a row found
// as another's mirror (rc is an involution). Mutation tried: a lookup per context: fails.
TEST(PatternSelection, MirrorLookupsOncePerRow) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    for (const std::string &p : { std::string("AC"), std::string("GTC"), std::string("GA") }) {
        Ask r;
        r.patterns = { p };
        r.predicate = "{\"any\": [\"c1\", \"c5\"]}";
        const std::vector<Ctx> raw = raw_order(idx, r);
        Outcome o = run(idx, r);
        std::set<std::string> known;
        uint64_t lookups = 0;
        for (const Ctx &c : raw) {
            if (known.count(c.kmer))
                continue;
            ++lookups;
            known.insert(c.kmer);
            if (oracle.present(rc(c.kmer)))
                known.insert(rc(c.kmer));
        }
        EXPECT_EQ(lookups, o.lookups) << p;
        // k units a lookup, inside the request's units
        EXPECT_GE(o.units, lookups * kK);
        check_invariants(oracle, parse_pattern_body(r.predicate), o, true, std::nullopt, r,
                         &raw);
    }
}

// A selective predicate on a pattern above max_contexts: the selected fit (all_or_count),
// where the same request without the predicate is withheld count_above_threshold (mutation
// tried: the list one short, selected < max_contexts: fails)
TEST(PatternSelection, SelectiveFilterAboveMaxContexts) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    ASSERT_GT(oracle.contexts("AC").size(), 10u);
    Ask r;
    r.predicate = "{\"all\": [\"c1\", \"c4\", \"c5\"]}";
    const Json::Value pred = parse_pattern_body(r.predicate);
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    ASSERT_GE(truth, 1u);
    // the raw contexts are above the threshold, the selected within it
    r.max_contexts = truth;
    ASSERT_GT(oracle.contexts("AC").size(), r.max_contexts);
    Outcome o = run(idx, r);
    EXPECT_EQ("completed", o.pass);
    EXPECT_EQ("", o.withheld);
    EXPECT_EQ(truth, o.results.size());
    EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
    const std::vector<Ctx> raw = raw_order(idx, r);
    check_invariants(oracle, pred, o, true, truth, r, &raw);
    Ask plain = r;
    plain.predicate.clear();
    const Outcome p = run(idx, plain);
    EXPECT_EQ("count_above_threshold", p.withheld);
    EXPECT_FALSE(p.entry.isMember("selection"));
    EXPECT_FALSE(p.answer.isMember("predicate"));
}

// The compute admission (§19.6 step 2): above max_predicate_contexts all_or_count and count
// read nothing (not_admitted; all_or_count withholds predicate_above_threshold), partial tests
// the first max_predicate_contexts (cut max_predicate_contexts). Mutations tried: not_admitted
// answered not_started; the bounds' upper end S + R: fail.
TEST(PatternSelection, ComputeAdmission) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    const Json::Value pred = parse_pattern_body("{\"any\": [\"c1\", \"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const std::string mode : { "all_or_count", "count", "partial" }) {
        auto reads = std::make_shared<uint64_t>(0);
        Ask r;
        r.predicate = compact(pred);
        r.mode = mode;
        r.max_predicate_contexts = raw.size() - 3;
        r.read_hook = [reads](size_t) { ++*reads; };
        Outcome o = run(idx, r);
        SCOPED_TRACE(mode);
        check_invariants(oracle, pred, o, true, truth, r, &raw);
        EXPECT_EQ(raw.size() - 3, o.answer["limits"]["max_predicate_contexts"].asUInt64());
        if (mode != "partial") {
            EXPECT_EQ("not_admitted", o.pass);
            EXPECT_EQ(0u, *reads);
            EXPECT_EQ("unknown", o.selected.relation);
            EXPECT_EQ("unknown", o.tested.relation);
            EXPECT_EQ(mode == "all_or_count" ? "predicate_above_threshold" : "", o.withheld);
            continue;
        }
        EXPECT_EQ(raw.size() - 3, o.tested.value);
        EXPECT_EQ("stopped", o.pass);
        EXPECT_FALSE(o.stop);
        EXPECT_EQ("bounds", o.selected.relation);
        EXPECT_EQ(o.selected.lower + 3, o.selected.upper);
        EXPECT_EQ("max_predicate_contexts", o.cut);
        // the first raw contexts tested, their selected ones listed
        const std::vector<Ctx> first(raw.begin(), raw.end() - 3);
        EXPECT_EQ(keys(oracle_selected(oracle, pred, first, true)), keys(o.results));
        // the pass's own stop comes before the raw release's cut (§19.8)
        r.max_predicate_work = 1;
        const Outcome w = run(idx, r);
        EXPECT_EQ("max_predicate_work", w.cut);
        EXPECT_EQ("stopped", w.pass);
    }
}

// More selected than max_contexts: all_or_count withholds selected_above_threshold, the counts
// kept, nothing of the list held after it (the next identical pattern peaks alike); partial
// lists the first max_contexts (cut max_contexts). Mutations tried: the list's condition
// selected < max_contexts (one fewer listed); the rows' entries never released: fail.
TEST(PatternSelection, SelectedAboveThreshold) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    const Json::Value pred = parse_pattern_body("{\"none\": [\"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    ASSERT_GT(truth, 3u);
    for (const std::string mode : { "all_or_count", "partial" }) {
        for (const std::string labels : { "none", "predicate_only" }) {
            Ask r;
            r.predicate = compact(pred);
            r.mode = mode;
            r.max_contexts = 3;
            r.labels = labels;
            r.patterns = { "AC", "AC" };
            const std::vector<Outcome> all = run_all(idx, r);
            const Outcome &o = all[0];
            SCOPED_TRACE(mode + " " + labels);
            check_invariants(oracle, pred, o, true, truth, r, &raw);
            EXPECT_EQ("exact", o.selected.relation);
            EXPECT_EQ(truth, o.selected.value);
            if (mode == "all_or_count") {
                EXPECT_EQ("selected_above_threshold", o.withheld);
                EXPECT_TRUE(o.results.empty());
                EXPECT_EQ(o.memory, all[1].memory) << "nothing of the list is held";
            } else {
                EXPECT_EQ(3u, o.results.size());
                EXPECT_EQ("max_contexts", o.cut);
                EXPECT_FALSE(o.entry["retrieval_complete"].asBool()) << "not every selected";
            }
        }
    }
}

// stop_at_threshold: on the selected count an early exit of the pass (stop {selection,
// max_contexts}); on the raw count, with max_predicate_contexts below it, a stop of the raw
// discovery (stop {discovery, max_predicate_contexts}). Mutation tried: the check removed:
// fails.
TEST(PatternSelection, StopAtThresholdRawAndSelected) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    const Json::Value pred = parse_pattern_body("{\"none\": [\"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    Ask r0;
    r0.predicate = compact(pred);
    const uint64_t all_rows = run(idx, r0).rows;
    for (const std::string mode : { "all_or_count", "partial" }) {
        Ask r = r0;
        r.mode = mode;
        r.max_contexts = 2;
        r.stop_at_threshold = true;
        Outcome o = run(idx, r);
        SCOPED_TRACE(mode);
        check_invariants(oracle, pred, o, true, truth, r, &raw);
        ASSERT_TRUE(o.stop);
        EXPECT_EQ("selection", o.stop->first);
        EXPECT_EQ("max_contexts", o.stop->second);
        EXPECT_EQ("stopped", o.pass);
        EXPECT_LT(o.rows, all_rows) << "an early exit";
        if (mode == "all_or_count") {
            EXPECT_EQ("threshold_crossed", o.withheld);
        } else {
            EXPECT_EQ("max_contexts", o.cut);
            EXPECT_EQ(2u, o.results.size());
        }
        // on the raw count
        r.max_contexts = 10'000;
        r.max_predicate_contexts = 3;
        o = run(idx, r);
        check_invariants(oracle, pred, o, true, truth, r, &raw);
        ASSERT_TRUE(o.stop);
        EXPECT_EQ("discovery", o.stop->first);
        EXPECT_EQ("max_predicate_contexts", o.stop->second);
        EXPECT_EQ("at_least", o.raw.relation);
        if (mode == "all_or_count") {
            EXPECT_EQ("threshold_crossed", o.withheld);
            EXPECT_EQ("not_started", o.pass);
        } else {
            EXPECT_EQ("max_predicate_contexts", o.cut);
        }
    }
}

// Mode count reads the annotation with a predicate: the counts, no list, the note
// projection_not_read for the projection it named; the pass's rows freed after it (the next
// identical pattern peaks alike; mutation tried: the rows' entries never released: fails)
TEST(PatternSelection, CountModeReads) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    Ask r;
    r.mode = "count";
    r.predicate = "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"c4\", \"c5\"]}}";
    r.labels = "predicate_only";
    r.patterns = { "AC", "AC" };
    const Json::Value pred = parse_pattern_body(r.predicate);
    const std::vector<Outcome> all = run_all(idx, r);
    const Outcome &o = all[0];
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r, &raw);
    EXPECT_EQ("completed", o.pass);
    EXPECT_GT(o.rows, 0u);
    EXPECT_TRUE(has_note(o, "projection_not_read"));
    EXPECT_FALSE(has_note(o, "annotation_not_read"));
    EXPECT_TRUE(o.answer["output"].isNull());
    EXPECT_EQ(o.memory, all[1].memory) << "the pass's rows freed";
    EXPECT_EQ(o.rows, all[1].rows);
}

// The scopes: suffix contexts only, the absence claim narrowed by the predicate (mutation
// tried: the mirror's labels dropped: fails)
TEST(PatternSelection, ScopesSuffixAndAnyOffset) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    for (const std::string scope : { "suffix", "any_offset" }) {
        Ask r;
        r.scope = scope;
        r.predicate = "{\"or\": [{\"any\": [\"c3\"]}, {\"none\": [\"c1\"]}]}";
        const Json::Value pred = parse_pattern_body(r.predicate);
        const std::vector<Ctx> raw = raw_order(idx, r);
        const gp::Scope s = scope == "suffix" ? gp::Scope::SUFFIX : gp::Scope::ANY_OFFSET;
        EXPECT_EQ(oracle.contexts("AC", s), ctx_set(raw));
        Outcome o = run(idx, r);
        check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true, s), r,
                         &raw);
        EXPECT_EQ(scope == "suffix" ? "suffix_only" : "any_offset",
                  o.entry["absence_scope"].asString());
        EXPECT_EQ("predicate", o.entry["absence_filter"].asString());
    }
}

// A peptide's contexts are contexts: MK (ATG, AAA/AAG) on the + and - strands (mutation
// tried: the mirror's labels dropped: fails)
TEST(PatternSelection, ProteinKind) {
    const std::vector<Record> records = {
        { "p1", "p1_r0", "CCATGAAAGG" },
        { "p2", "p2_r0", "TTCTTTCATT" },
        { "p3", "p3_r0", "GATGAAGTCC" },
    };
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, records, true);
    Oracle oracle(idx);
    for (const bool either : { true, false }) {
        Ask r;
        r.patterns = { "MK" };
        r.kind = "protein";
        r.either = either;
        r.predicate = "{\"none\": [\"p2\"]}";
        const std::vector<Ctx> raw = raw_order(idx, r);
        ASSERT_GT(raw.size(), 0u);
        for (const Ctx &c : raw) {
            std::string inst = c.kmer.substr(c.offset, 6);
            if (c.strand == "-")
                inst = rc(inst);
            EXPECT_EQ("ATG", inst.substr(0, 3));
            EXPECT_TRUE(inst.substr(3) == "AAA" || inst.substr(3) == "AAG") << inst;
        }
        Outcome o = run(idx, r);
        check_invariants(oracle, parse_pattern_body(r.predicate), o, either, std::nullopt, r,
                         &raw);
        EXPECT_EQ("completed", o.pass);
    }
}

/**
 * The memory account from the smallest up (every 1/150 of what the full run peaked at): each
 * stop stated where it falls (the binding, the descriptors, the rows, the list, the results),
 * never past the account, the decisions made right, complete once the account holds the full
 * run's peak. Mutation tried: a refused row's statement not charged (the reserve kept but not
 * held): the full run's peak drops and the answers at that peak change; and the descriptors'
 * allowance ignored in partial: an answer past the account.
 */
TEST(PatternSelection, MemoryAccount) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    const Json::Value pred = parse_pattern_body("{\"or\": [{\"any\": [\"c1\"]}, "
                                                "{\"none\": [\"c2\", \"c5\"]}]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const std::string mode : { "all_or_count", "partial" }) {
        for (const std::string labels : { "none", "predicate_only" }) {
            Ask r;
            r.predicate = compact(pred);
            r.mode = mode;
            r.labels = labels;
            r.max_memory_bytes = uint64_t(1) << 30;
            const Outcome full = run(idx, r);
            ASSERT_EQ("completed", full.pass);
            bool some_stopped = false, output_stopped = false;
            for (uint64_t step = 1; step <= 160; ++step) {
                r.max_memory_bytes = std::max<uint64_t>(1, full.memory * step / 150);
                const Outcome o = run(idx, r);
                SCOPED_TRACE(mode + " " + labels + " " + std::to_string(r.max_memory_bytes)
                             + " bytes");
                EXPECT_LE(o.memory, r.max_memory_bytes);
                check_invariants(oracle, pred, o, true, truth, r, &raw);
                if (!o.bound) {
                    // the bound predicate did not fit: no pass
                    EXPECT_EQ("not_started", o.pass);
                    ASSERT_TRUE(o.stop);
                    EXPECT_EQ("selection", o.stop->first);
                    EXPECT_EQ("max_memory", o.stop->second);
                }
                const bool complete = o.entry["retrieval_complete"].asBool();
                if (o.pass != "completed" || !complete) {
                    some_stopped = true;
                    const bool stated = (o.stop && o.stop->second == "max_memory")
                                            || o.entry["rows_refused"].size();
                    EXPECT_TRUE(stated);
                    if (mode == "all_or_count") {
                        // the pass's, the results', or the projection's reads
                        EXPECT_TRUE(o.withheld == "predicate_budget"
                                    || o.withheld == "output_budget"
                                    || o.withheld == "annotation_budget") << o.withheld;
                    } else if (o.pass != "completed") {
                        EXPECT_EQ("max_memory", o.cut);
                    } else {
                        // a projection's refused row alone leaves the list whole, no cut
                        EXPECT_TRUE(o.cut.empty() || o.cut == "max_memory") << o.cut;
                    }
                }
                for (const Json::Value &x : o.entry["rows_refused"]) {
                    EXPECT_EQ("max_memory", x["reason"].asString());
                    EXPECT_EQ(kK, x["kmer"].asString().size());
                }
                // the result objects of the selected contexts are charged too (512 + 2k each)
                output_stopped |= o.pass == "completed" && o.stop
                        && *o.stop == std::make_pair(std::string("output"),
                                                     std::string("max_memory"));
            }
            EXPECT_TRUE(some_stopped);
            EXPECT_TRUE(output_stopped) << "an account the pass fits and its results do not";
            // the account holds the run's peak and the reads' own decoding (a read's demand
            // is admitted against what is left, beside what the account holds)
            r.max_memory_bytes = 4 * full.memory;
            const Outcome o = run(idx, r);
            EXPECT_EQ("completed", o.pass);
            EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
        }
    }
}

// The descriptors the account could not hold end the pass's set (§19.9): partial tests what
// was admitted (cut max_memory), all_or_count reads nothing. 31 contexts of A on one row (a
// poly-A k-mer at every offset), the account sized to admit 30 (an unbudgeted annotation: its
// read holds only its hits, charged after it, so that the arithmetic is the descriptors' and
// the row's). Mutation tried: the admission cut ending the pass in every mode: fails.
TEST(PatternSelection, AdmissionCutTestsTheAdmitted) {
    const size_t k = 31;
    Index idx = build<annot::ColumnCompressed<>>(
            k, { { "c1", "c1_r0", std::string(40, 'A') }, { "c2", "c2_r0", std::string(40, 'C') } },
            false);
    Oracle oracle(idx);
    Ask r;
    r.patterns = { "A" };
    r.predicate = "{\"any\": [\"c1\"]}";
    r.either = false;
    r.records = false;
    r.allow_unbudgeted = true;
    r.mode = "partial";
    const Outcome full = run(idx, r);
    ASSERT_EQ(k, full.tested.value);
    ASSERT_EQ("completed", full.pass);
    // the bound predicate's model (SPEC §19.9): its one label (192 + 2 x 2), one node (64),
    // one name listed in a leaf (8)
    const uint64_t bound_bytes = 192 + 4 + 64 + 8;
    // partial's allowance is half of what the bound predicate left: 30 descriptors of 64
    r.max_memory_bytes = bound_bytes + 2 * (64 * 30 + 32);
    const Outcome o = run(idx, r);
    EXPECT_EQ(30u, o.tested.value);
    EXPECT_EQ("stopped", o.pass);
    EXPECT_EQ("bounds", o.selected.relation);
    EXPECT_EQ(30u, o.selected.lower);
    EXPECT_EQ(31u, o.selected.upper);
    ASSERT_TRUE(o.stop);
    EXPECT_EQ("selection", o.stop->first);
    EXPECT_EQ("max_memory", o.stop->second);
    EXPECT_EQ("max_memory", o.cut);
    EXPECT_EQ(1u, o.rows);
    // the results' own objects (512 + 2k each) are what is left of the account: fewer
    EXPECT_GE(30u, o.results.size());
    check_invariants(oracle, parse_pattern_body(r.predicate), o, false, 31, r);
    // one byte less: 29 admitted
    r.max_memory_bytes -= 2 * 32 + 2;
    EXPECT_EQ(29u, run(idx, r).tested.value);
    // all_or_count: nothing can be published, nothing read
    r.mode = "all_or_count";
    r.max_memory_bytes = bound_bytes + 64 * 30 + 32;
    const Outcome w = run(idx, r);
    EXPECT_EQ(0u, w.rows);
    EXPECT_EQ("predicate_budget", w.withheld);
    EXPECT_EQ("stopped", w.pass);
    EXPECT_EQ(0u, w.tested.value);
    // count: its descriptors may take the whole account, as all_or_count's, so a cut leaves no
    // room for a row: stopped at the first, stated, the bounds over every raw context
    r.mode = "count";
    const Outcome c = run(idx, r);
    EXPECT_EQ(0u, c.tested.value);
    EXPECT_EQ(0u, c.rows);
    EXPECT_EQ("stopped", c.pass);
    EXPECT_EQ("bounds", c.selected.relation);
    EXPECT_EQ(31u, c.selected.upper);
    ASSERT_TRUE(c.stop);
    EXPECT_EQ(std::make_pair(std::string("selection"), std::string("max_memory")), *c.stop);
}

// A row the account cannot hold: stated (rows_refused, phase selection); all_or_count ends
// the pass there, partial goes on, the contexts needing it untested. Mutation tried: the
// statement's reserve not kept (a refused row unstated): fails.
TEST(PatternSelection, RefusedRowsStated) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::vector<Ctx> raw = raw_order(idx, Ask());
    const std::vector<std::string> order = row_order(oracle, raw, true);
    const Json::Value pred = parse_pattern_body("{\"any\": [\"c1\", \"c3\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const std::string mode : { "all_or_count", "partial", "count" }) {
        // the second read refused
        auto reads = std::make_shared<uint64_t>(0);
        auto denied = std::make_shared<bool>(false);
        Ask r;
        r.predicate = compact(pred);
        r.mode = mode;
        r.read_hook = [reads](size_t) { ++*reads; };
        // one charge of the second read refused: that row only
        r.deny_decode = [reads, denied](uint64_t) {
            if (*reads != 2 || *denied)
                return false;
            *denied = true;
            return true;
        };
        Outcome o = run(idx, r);
        SCOPED_TRACE(mode);
        check_invariants(oracle, pred, o, true, truth, r, &raw);
        const Json::Value &refused = o.entry["rows_refused"];
        ASSERT_EQ(1u, refused.size());
        EXPECT_EQ("selection", refused[0]["phase"].asString());
        EXPECT_EQ("stopped", o.pass);
        EXPECT_EQ(order.at(1), refused[0]["kmer"].asString());
        if (mode == "all_or_count") {
            EXPECT_EQ("predicate_budget", o.withheld);
            EXPECT_EQ(1u, o.rows);
        } else {
            // partial and count go on (count publishes nothing: no withheld, no cut)
            EXPECT_EQ(mode == "partial" ? "max_memory" : "", o.cut);
            EXPECT_LE(o.rows, order.size() - 1);
            // the contexts needing the refused row are the untested ones; the others decided
            std::vector<Ctx> decided;
            for (const Ctx &c : raw) {
                if (c.kmer != order[1] && rc(c.kmer) != order[1])
                    decided.push_back(c);
            }
            EXPECT_EQ(decided.size(), o.tested.value);
            if (mode == "partial") {
                EXPECT_EQ(keys(oracle_selected(oracle, pred, decided, true)), keys(o.results));
            } else {
                EXPECT_EQ(oracle_selected(oracle, pred, decided, true).size(), o.selected.lower);
            }
        }
    }
}

/**
 * A virtual clock past the work time at each reading in turn (the binding's, the engine's,
 * the lookups', each read and its paced pieces, each 64 decisions, the results'): a time stop
 * stated, time_limited, never an exception; at the reading after the last, the pass
 * completes. Mutations tried: an interrupted read's stop unstated; the decisions' time stop
 * unstated: fail. (The clock before a read removed is equivalent here, the paced read seeing
 * the deadline itself: the abort test catches it.)
 */
TEST(PatternSelection, DeadlineAtEveryReading) {
    // R1 with AC; and k = 33 with one poly-A k-mer carrying 66 contexts of D (A, G or T: its
    // reverse complement H matches A too), decided after one read: the decisions' own clock
    // (every 64) is read there, and the 66 results' (every 64)
    struct Case {
        Index idx;
        std::string pattern;
        std::string kind;
        std::string predicate;
    };
    std::vector<Case> cases;
    cases.push_back({ build<annot::RowDiffColumnAnnotator>(kK, kR1, true), "AC", "dna",
                      "{\"none\": [\"c4\"]}" });
    cases.push_back({ build<annot::RowDiffColumnAnnotator>(
                              33, { { "c1", "c1_r0", std::string(50, 'A') },
                                    { "c2", "c2_r0", std::string(50, 'C') } }, false),
                      "D", "iupac", "{\"any\": [\"c1\"]}" });
    for (const Case &c : cases) {
        Oracle oracle(c.idx);
        const Json::Value pred = parse_pattern_body(c.predicate);
        for (const auto &[mode, labels] : std::vector<std::pair<std::string, std::string>>{
                { "all_or_count", "predicate_only" }, { "partial", "predicate_only" },
                { "all_or_count", "none" }, { "partial", "none" } }) {
            auto readings = std::make_shared<uint64_t>(0);
            auto at = std::make_shared<uint64_t>(UINT64_MAX);
            Ask r;
            r.patterns = { c.pattern };
            r.kind = c.kind;
            r.records = c.idx.cth != nullptr;
            r.predicate = c.predicate;
            r.mode = mode;
            r.labels = labels;
            r.clock = [readings, at]() {
                ++*readings;
                return Clock::time_point() + (*readings >= *at ? std::chrono::hours(1)
                                                               : std::chrono::hours(0));
            };
            const Outcome full = run(c.idx, r);
            ASSERT_EQ("completed", full.pass);
            const uint64_t total = *readings;
            ASSERT_GT(total, 3u);
            bool pass_stopped = false, output_stopped = false;
            for (uint64_t n = 1; n <= total + 1; ++n) {
                *readings = 0;
                *at = n;
                Outcome o;
                ASSERT_NO_THROW(o = run(c.idx, r));
                SCOPED_TRACE(c.pattern + " " + mode + " " + labels + " reading "
                             + std::to_string(n));
                check_invariants(oracle, pred, o, true, std::nullopt, r);
                EXPECT_TRUE(contains(o.selected, full.selected.value));
                if (n == total + 1) {
                    EXPECT_EQ("completed", o.pass);
                    EXPECT_TRUE(o.entry["retrieval_complete"].asBool());
                    continue;
                }
                // an incomplete answer says why: a time stop, the engine's or its own
                const bool complete = o.entry["retrieval_complete"].asBool();
                if (o.pass != "completed" || !complete) {
                    EXPECT_EQ("time_limited", o.determinism);
                    ASSERT_TRUE(o.stop) << compact(o.entry);
                    EXPECT_EQ("time", o.stop->second);
                    EXPECT_EQ(mode == "all_or_count" ? "deadline" : "", o.withheld);
                    // the list cut by the time stop of the search, the pass or the results;
                    // a time stop of the projection's reads or of its labels' output alone
                    // leaves the list whole (no cut, increment 3)
                    if (mode == "partial") {
                        if (o.stop->first == "output" || o.stop->first == "placement") {
                            EXPECT_TRUE(o.cut.empty() || o.cut == "time") << o.cut;
                        } else {
                            EXPECT_EQ("time", o.cut) << o.stop->first;
                        }
                    }
                }
                if (o.stop && o.stop->second == "time") {
                    pass_stopped |= o.stop->first == "selection";
                    output_stopped |= o.stop->first == "output";
                    EXPECT_EQ("time_limited", o.determinism);
                }
            }
            EXPECT_TRUE(pass_stopped);
            EXPECT_TRUE(output_stopped);
        }
    }
}

// A caller that left is not answered: the pass reads the clock (and so asks the abort
// predicate) before every row it reads, and the work ends there (Aborted). Mutation tried: the
// clock before a read removed (the paced read still sees the deadline, not the abort): fails.
TEST(PatternSelection, AbortedCallerEndsThePass) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    for (const bool either : { false, true }) {
        auto left = std::make_shared<bool>(false);
        auto reads = std::make_shared<uint64_t>(0);
        Ask r;
        r.predicate = "{\"any\": [\"c1\"]}";
        r.either = either;
        // the caller leaves during the first read
        r.read_hook = [left, reads](size_t) {
            ++*reads;
            *left = true;
        };
        r.abort = [left]() { return *left; };
        EXPECT_THROW(run(idx, r), gp::Aborted);
        EXPECT_EQ(1u, *reads) << "no row read after the caller left";
    }
}

// A max_predicate_work stop is sticky for the later patterns' selection (not_started,
// predicate_budget); their discovery still runs (mutation tried: no sticky check: fails)
TEST(PatternSelection, StickyWorkStopForLaterPatterns) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    // (A: more raw contexts than max_predicate_contexts, which AC's fit: after the sticky stop
    // its selection is not_started, never not_admitted)
    Ask ac;
    const uint64_t raw_ac = raw_order(idx, ac).size();
    for (const std::string mode : { "all_or_count", "partial", "count" }) {
        Ask r;
        r.patterns = { "AC", "GA", "TT", "A" };
        r.predicate = "{\"any\": [\"c1\"]}";
        r.max_predicate_work = 40;
        r.max_predicate_contexts = raw_ac;
        r.mode = mode;
        std::vector<Outcome> all = run_all(idx, r);
        ASSERT_EQ(4u, all.size());
        ASSERT_TRUE(all[0].stop);
        EXPECT_EQ("max_predicate_work", all[0].stop->second);
        ASSERT_GT(oracle.contexts("A").size(), raw_ac);
        for (size_t p = 1; p < 4; ++p) {
            SCOPED_TRACE(mode + " pattern " + std::to_string(p));
            EXPECT_EQ("not_started", all[p].pass);
            ASSERT_TRUE(all[p].stop);
            EXPECT_EQ("selection", all[p].stop->first);
            EXPECT_EQ("max_predicate_work", all[p].stop->second);
            EXPECT_EQ(mode == "all_or_count" ? "predicate_budget" : "", all[p].withheld);
            EXPECT_EQ(mode == "partial" ? "max_predicate_work" : "", all[p].cut);
            EXPECT_EQ(0u, all[p].rows);
            EXPECT_EQ("unknown", all[p].selected.relation);
            EXPECT_EQ("exact", all[p].raw.relation);
            EXPECT_EQ(oracle.contexts(r.patterns[p]).size(), all[p].raw.value);
        }
    }
}

// An unbudgeted annotation with direct access: single cells for at most 16 known labels
// ("columns"), rows above; stated by the note annotation_unbudgeted. Mutation tried: always
// rows: fails.
TEST(PatternSelection, UnbudgetedAnnotation) {
    Index small = build<annot::ColumnCompressed<>>(kK, kR1, true);
    Oracle so(small);
    Ask r;
    r.allow_unbudgeted = true;
    r.predicate = "{\"or\": [{\"any\": [\"c1\"]}, {\"none\": [\"c2\", \"c5\"]}]}";
    Outcome o = run(small, r);
    EXPECT_EQ("columns", o.access);
    EXPECT_TRUE(has_note(o, "annotation_unbudgeted"));
    const std::vector<Ctx> raw = raw_order(small, r);
    check_invariants(so, parse_pattern_body(r.predicate), o, true,
                     oracle_truth(so, parse_pattern_body(r.predicate), "AC", true), r, &raw);
    // without the opt-in: refused in every mode
    for (const std::string mode : { "count", "all_or_count", "partial" }) {
        Ask q = r;
        q.allow_unbudgeted = false;
        q.mode = mode;
        try {
            run(small, q);
            ADD_FAILURE() << "not refused";
        } catch (const PatternRefusal &e) {
            EXPECT_EQ("annotation_unbudgeted", e.code());
        }
    }
    Index wide = build<annot::ColumnCompressed<>>(7, random_records(30, 3), true);
    Oracle wo(wide);
    Json::Value list(Json::arrayValue);
    for (size_t c = 0; c < 20; ++c) {
        list.append((c < 10 ? "r0" : "r") + std::to_string(c));
    }
    Json::Value pred;
    pred["at_least"]["n"] = 2;
    pred["at_least"]["labels"] = list;
    r.predicate = compact(pred);
    o = run(wide, r);
    EXPECT_EQ("rows", o.access);
    const std::vector<Ctx> wraw = raw_order(wide, r);
    check_invariants(wo, pred, o, true, oracle_truth(wo, pred, "AC", true), r, &wraw);
}

// A graph without its dummy-edge mask (owner decision #16): the admitted release enumerates
// the candidates and drops the dummies; the selection is the masked twin's and the oracle's.
// The compute admission on such a graph compares the raw count's upper bound U (stated by the
// note threshold_upper_bound when its lower bound fits). Mutation tried: the mirror's labels
// dropped: fails.
TEST(PatternSelection, UnmaskedGraphAsTheMaskedTwin) {
    Index masked = build<annot::ColumnCompressed<>>(kK, kR1, false);
    const auto &dbg = dynamic_cast<const DBGSuccinct&>(masked.anno->get_graph());
    const std::string base = test_dump_dir() + "/pattern_selection_unmasked";
    dbg.serialize(base);
    masked.anno->get_annotator().serialize(base);
    auto loaded = std::make_shared<DBGSuccinct>(2);
    ASSERT_TRUE(loaded->load_without_mask(base + ".dbg"));
    ASSERT_EQ(nullptr, loaded->get_mask());
    auto annotation = std::make_unique<annot::ColumnCompressed<>>();
    ASSERT_TRUE(annotation->load(base));
    Index unmasked;
    unmasked.k = kK;
    unmasked.records = kR1;
    unmasked.anno = std::make_unique<AnnotatedDBG>(loaded, std::move(annotation));
    Oracle oracle(masked);
    for (const std::string &p : { std::string("AC"), std::string("GGA") }) {
        Ask r;
        r.patterns = { p };
        r.records = false;
        r.allow_unbudgeted = true;
        r.predicate = "{\"or\": [{\"any\": [\"c1\"]}, {\"none\": [\"c2\", \"c5\"]}]}";
        const Outcome a = run(masked, r), b = run(unmasked, r);
        EXPECT_EQ("exact", b.raw.relation) << p;
        EXPECT_EQ(keys(a.results), keys(b.results)) << p;
        EXPECT_EQ(a.selected.value, b.selected.value);
        EXPECT_EQ("completed", b.pass);
        const std::vector<Ctx> raw = raw_order(masked, r);
        check_invariants(oracle, parse_pattern_body(r.predicate), b, true, std::nullopt, r,
                         &raw);
    }
    // no unchecked candidate tested (--pattern-max-checked-entries 0): the raw count of AC is
    // bounds [lower, U] (mode count without a predicate); the compute admission compares U and
    // says so when the lower bound fits (all_or_count and count, which runs the admission too)
    Ask r;
    r.patterns = { "AC" };
    r.records = false;
    r.allow_unbudgeted = true;
    r.mode = "count";
    r.max_checked_entries = 0;
    r.predicate.clear();
    const Outcome c = run(unmasked, r);
    ASSERT_EQ("bounds", c.raw.relation);
    ASSERT_LT(c.raw.lower, c.raw.upper);
    r.predicate = "{\"any\": [\"c1\"]}";
    r.max_predicate_contexts = c.raw.upper - 1;
    for (const std::string mode : { "all_or_count", "count" }) {
        r.mode = mode;
        const Outcome u = run(unmasked, r);
        SCOPED_TRACE(mode);
        EXPECT_EQ(mode == "count" ? "" : "predicate_above_threshold", u.withheld);
        EXPECT_EQ("not_admitted", u.pass);
        EXPECT_EQ("bounds", u.raw.relation);
        EXPECT_TRUE(has_note(u, "threshold_upper_bound"));
    }
    // admitted on U: the release makes the raw count exact, the selection the oracle's
    r.max_predicate_contexts = c.raw.upper;
    r.mode = "all_or_count";
    const Outcome a = run(unmasked, r);
    EXPECT_EQ("exact", a.raw.relation);
    EXPECT_EQ(oracle.contexts("AC").size(), a.raw.value);
    EXPECT_EQ("completed", a.pass);
    EXPECT_EQ(oracle_truth(oracle, parse_pattern_body(r.predicate), "AC", true),
              a.selected.value);
}

// The three projections select the same contexts; predicate_only lists the predicate's labels
// on each context's own row, placed as the record scan says; all lists every label, as the
// request without a predicate lists them for the same contexts (increment 3); the answer's
// shape: output, the predicate block, limits (TESTS §3 ProjectionNonePredicateOnlyAll).
// Mutation tried: retrieve() for predicate_only (every label listed): fails.
TEST(PatternSelection, ProjectionNonePredicateOnlyAll) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body("{\"and\": [{\"any\": [\"c1\", \"c3\"]}, "
                                                "{\"none\": [\"c2\"]}]}");
    std::map<std::string, Outcome> by;
    for (const std::string labels : { "none", "predicate_only", "all" }) {
        Ask r;
        r.predicate = compact(pred);
        r.labels = labels;
        by[labels] = run(idx, r);
        const Outcome &o = by[labels];
        EXPECT_EQ(labels, o.answer["output"]["labels"].asString());
        EXPECT_EQ(labels != "none", o.answer["output"].isMember("occurrences"));
        const Json::Value &p = o.answer["predicate"];
        EXPECT_EQ("shard_context", p["scope"].asString());
        EXPECT_EQ("either", p["strands"].asString());
        EXPECT_EQ(3u, p["names"].asUInt64());
        EXPECT_EQ(3u, p["known"].asUInt64());
        const Json::Value &l = o.answer["limits"];
        EXPECT_EQ(100000u, l["max_predicate_contexts"].asUInt64());
        EXPECT_EQ(100000000u, l["max_predicate_work"].asUInt64());
        EXPECT_EQ(10000u, l["max_predicate_labels"].asUInt64());
        EXPECT_EQ(256u, l["max_memory_mb"].asUInt64());
        EXPECT_EQ(labels != "none", o.entry.isMember("by_label"));
        EXPECT_TRUE(o.entry["rows_refused"].isArray());
    }
    EXPECT_EQ(keys(by["none"].results), keys(by["predicate_only"].results));
    EXPECT_EQ(keys(by["none"].results), keys(by["all"].results));
    ASSERT_GT(by["none"].results.size(), 0u);
    // all: the labels of increment 3 for the same contexts
    Ask plain;
    plain.predicate.clear();
    plain.labels = "all";
    const Outcome p = run(idx, plain);
    // (as sets: a result lists its labels in the pattern's label order, which depends on the
    // contexts returned)
    auto label_set = [](const Json::Value &labels) {
        std::set<std::string> out;
        for (const Json::Value &l : labels) {
            out.insert(compact(l));
        }
        return out;
    };
    std::map<std::tuple<std::string, std::string, uint64_t>, std::set<std::string>> labels_of;
    for (const Ctx &c : p.results) {
        labels_of[{ c.strand, c.kmer, c.offset }] = label_set(c.json["labels"]);
    }
    for (const Ctx &c : by["all"].results) {
        EXPECT_EQ(labels_of.at({ c.strand, c.kmer, c.offset }), label_set(c.json["labels"]));
    }
    for (const Ctx &c : by["predicate_only"].results) {
        for (const Json::Value &l : c.json["labels"]) {
            const std::string column = l["column"].asString();
            EXPECT_TRUE(column == "c1" || column == "c3") << column;
            EXPECT_TRUE(oracle.own(c.kmer).count(column));
        }
    }
    // the same request without a predicate states none of it
    EXPECT_FALSE(p.answer.isMember("predicate"));
    EXPECT_FALSE(p.answer["limits"].isMember("max_predicate_contexts"));
    EXPECT_FALSE(p.entry.isMember("selection"));
    EXPECT_FALSE(p.entry.isMember("absence_filter"));
    EXPECT_FALSE(p.entry["counts"].isMember("tested"));
}

// A pattern longer than k under long_search "anchors": its anchors' answer as without the
// predicate, its selection not_started (SPEC §19.5), selected exact 0 when it has no anchor
TEST(PatternSelection, LongPatternNotStarted) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    for (const std::string mode : { "all_or_count", "count" }) {
        Ask r;
        r.patterns = { "GGACGACTTT", "GGGGGGGGGGG", "AC" };
        r.mode = mode;
        r.labels = mode == "count" ? "none" : "predicate_only";
        const Json::Value a = answer_of(idx, r);
        Ask plain = r;
        plain.predicate.clear();
        plain.labels = "none";
        const Json::Value b = answer_of(idx, plain);
        for (Json::ArrayIndex i = 0; i < 2; ++i) {
            const Json::Value &e = a["patterns"][i], &f = b["patterns"][i];
            SCOPED_TRACE(mode + " " + e["pattern"].asString());
            EXPECT_EQ("not_started", e["selection"]["pass"].asString());
            EXPECT_EQ("predicate", e["absence_filter"].asString());
            EXPECT_EQ(compact(f["counts"]["anchors"]), compact(e["counts"]["anchors"]));
            EXPECT_EQ(compact(f["withheld"]), compact(e["withheld"]));
            EXPECT_EQ("unknown", e["counts"]["tested"]["relation"].asString());
            EXPECT_EQ("paths", e["counts"]["tested"]["unit"].asString());
            EXPECT_EQ("paths", e["counts"]["selected"]["unit"].asString());
            const bool none = f["counts"]["anchors"]["value"].asUInt64() == 0;
            EXPECT_EQ(none ? "exact" : "unknown", e["counts"]["selected"]["relation"].asString());
            EXPECT_EQ(0u, e["work"]["predicate_rows"].asUInt64());
        }
        EXPECT_EQ("completed", a["patterns"][2]["selection"]["pass"].asString());
    }
}

// A palindromic k-mer (k = 6, GAATTC) is its own reverse complement: under "either" a label
// on its row supports it on both strands ("both"), under "context" on its own ("context"),
// against the records (check_invariants). Mutation tried: a palindrome's own labels marked
// "context" under "either": fails.
TEST(PatternSelection, PalindromicKmerOnBothStrands) {
    Index idx = build<annot::RowDiffColumnAnnotator>(
            6, { { "c1", "c1_r0", "TTGAATTCAA" }, { "c2", "c2_r0", "CCGAATTGCC" } }, true);
    Oracle oracle(idx);
    ASSERT_EQ("GAATTC", rc("GAATTC"));
    ASSERT_TRUE(oracle.own("GAATTC").count("c1"));
    for (const bool either : { true, false }) {
        Ask r;
        r.patterns = { "GAAT" };
        r.predicate = "{\"any\": [\"c1\", \"c2\"]}";
        r.either = either;
        r.labels = "predicate_only";
        const std::vector<Ctx> raw = raw_order(idx, r);
        const Outcome o = run(idx, r);
        const Json::Value pred = parse_pattern_body(r.predicate);
        check_invariants(oracle, pred, o, either, oracle_truth(oracle, pred, "GAAT", either), r,
                         &raw);
        bool palindrome = false;
        for (size_t j = 0; j < o.results.size(); ++j) {
            if (o.results[j].kmer != "GAATTC")
                continue;
            palindrome = true;
            EXPECT_EQ(std::vector<std::string>({ "c1" }), o.selection_labels[j]);
            EXPECT_EQ(std::vector<std::string>({ either ? "both" : "context" }),
                      o.selection_strands[j]);
        }
        EXPECT_TRUE(palindrome) << compact(o.entry);
    }
}

// The engine stopped before it released a context (max_steps 1): nothing reached the pass, so
// every mode says not_started, its counts unknown (§19.7; the review of 5b, L1: partial said
// stopped with tested exact 0 while all_or_count and count said not_started). Mutation tried:
// the partial branch removed: fails.
TEST(PatternSelection, EngineStopBeforeTheReleaseNotStarted) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    for (const std::string mode : { "partial", "all_or_count", "count" }) {
        for (const std::string labels : { "none", "predicate_only" }) {
            Ask r;
            r.predicate = "{\"any\": [\"c1\"]}";
            r.mode = mode;
            r.labels = labels;
            r.max_steps = 1;
            const Outcome o = run(idx, r);
            SCOPED_TRACE(mode + " " + labels);
            ASSERT_FALSE(o.entry.isMember("error")) << compact(o.entry);
            ASSERT_TRUE(o.stop) << compact(o.entry);
            EXPECT_EQ("discovery", o.stop->first);
            EXPECT_EQ("max_steps", o.stop->second);
            EXPECT_EQ("not_started", o.pass);
            EXPECT_EQ("unknown", o.tested.relation);
            EXPECT_EQ("unknown", o.selected.relation);
            EXPECT_EQ(0u, o.rows);
            EXPECT_TRUE(o.results.empty());
            if (mode == "partial") {
                EXPECT_EQ("max_steps", o.cut);
            } else if (mode == "all_or_count") {
                EXPECT_EQ("discovery_budget", o.withheld);
            }
        }
    }
}

} // namespace
