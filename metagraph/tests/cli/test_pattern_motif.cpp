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
#include "cli/pattern.hpp"
#include "cli/pattern_predicate.hpp"
#include "cli/pattern_retrieval.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// The motif-level predicate of a pattern of L <= k (SelectionRequest::motif,
// src/cli/pattern_selection.cpp): the request's normal form evaluated once on the union of the
// predicate's labels over the pattern's graph contexts, from the selection pass's own reads.
// Driven through PatternRetrieval (bind, the engine's release into the pass, select) as the
// route drives it. The expectations come from ORACLES over the records, never the modules under
// test: the graph contexts of a pattern scanned from the records' k-mers (and their reverse
// complements on a CANONICAL graph), the columns whose records hold a k-mer as deposited, the
// union over the contexts (with "either" the reverse complement's columns too), a recursive
// two-valued evaluator of the request's predicate JSON and a strong-Kleene three-valued one.

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

struct Index {
    size_t k = 0;
    DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
};

// a succinct graph of |records| at |k| (with its dummy-edge mask), a row-diff annotation with
// coordinates on BASIC (as refseq33m's), without on CANONICAL
Index build(size_t k, const std::vector<Record> &records,
            DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    Index idx;
    idx.k = k;
    idx.mode = mode;
    idx.records = records;
    const bool coordinates = mode == DeBruijnGraph::BASIC;
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
    idx.anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            k, seqs, labels, mode, coordinates, coordinates ? starts : std::vector<uint64_t>{});
    return idx;
}

// k = 7: c5 holds c1/r0's motif on the other strand only (its reverse complement), c4 the same
// record as c1/r0, c3 a palindromic repeat (ACGT is its own reverse complement)
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

// |columns| columns r00, r01, ... of one or two random records each, some of them the reverse
// complement of an earlier column's record (seeded)
std::vector<Record> random_records(size_t columns, unsigned seed) {
    std::mt19937 rng(seed);
    std::vector<Record> out;
    for (size_t c = 0; c < columns; ++c) {
        const std::string name = (c < 10 ? "r0" : "r") + std::to_string(c);
        const size_t n = 1 + rng() % 2;
        for (size_t r = 0; r < n; ++r) {
            std::string s;
            if (!out.empty() && rng() % 4 == 0) {
                s = rc(out[rng() % out.size()].seq);
            } else {
                s.assign(12 + rng() % 19, 'A');
                for (char &x : s) {
                    x = "ACGT"[rng() % 4];
                }
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

// the request's predicate on a set of column names, by recursion over its JSON (an unknown name
// is simply absent from every set)
bool o_eval(const Json::Value &p, const std::set<std::string> &s) {
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &v = p[op];
    auto hits = [&](const Json::Value &list) {
        uint64_t n = 0;
        for (const Json::Value &x : list) {
            n += s.count(x.asString());
        }
        return n;
    };
    if (op == "any")
        return hits(v) > 0;
    if (op == "all")
        return hits(v) == v.size();
    if (op == "none")
        return hits(v) == 0;
    if (op == "at_least")
        return hits(v["labels"]) >= v["n"].asUInt64();
    if (op == "not")
        return !o_eval(v, s);
    const bool is_and = op == "and";
    for (const Json::Value &q : v) {
        if (o_eval(q, s) != is_and)
            return !is_and;
    }
    return is_and;
}

// strong Kleene over the JSON: the names |sure| present, the names |maybe| possibly present,
// every other name absent; nullopt for "cannot tell"
std::optional<bool> o_kleene(const Json::Value &p, const std::set<std::string> &sure,
                             const std::set<std::string> &maybe) {
    const std::string op = p.getMemberNames().at(0);
    const Json::Value &v = p[op];
    auto counts = [&](const Json::Value &list) {
        uint64_t s = 0, m = 0;
        for (const Json::Value &x : list) {
            if (sure.count(x.asString())) {
                ++s;
            } else if (maybe.count(x.asString())) {
                ++m;
            }
        }
        return std::make_pair(s, m);
    };
    if (op == "any" || op == "all" || op == "none" || op == "at_least") {
        const Json::Value &list = op == "at_least" ? v["labels"] : v;
        const auto [s, m] = counts(list);
        const uint64_t n = op == "any" || op == "none" ? 1
                         : op == "all" ? list.size() : v["n"].asUInt64();
        std::optional<bool> at_least_n;
        if (s >= n) {
            at_least_n = true;
        } else if (s + m < n) {
            at_least_n = false;
        }
        if (op == "none" && at_least_n)
            return !*at_least_n;
        return at_least_n;
    }
    if (op == "not") {
        const std::optional<bool> q = o_kleene(v, sure, maybe);
        return q ? std::optional<bool>(!*q) : std::nullopt;
    }
    const bool is_and = op == "and";
    bool undecided = false;
    for (const Json::Value &q : v) {
        const std::optional<bool> r = o_kleene(q, sure, maybe);
        if (!r) {
            undecided = true;
        } else if (*r != is_and) {
            return !is_and;
        }
    }
    return undecided ? std::nullopt : std::optional<bool>(is_and);
}

/**
 * The folding of SPEC §19.4 by its table, to know which names the normal form keeps (the labels
 * the pass reads and the motif's union is over): unknown names (not in |columns|) folded away,
 * constants through and, or and not. Returns the folded predicate, or a JSON boolean.
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

using CtxKey = std::tuple<std::string, std::string, uint64_t>;   // strand, k-mer, offset

// The records' view of an index: the columns holding a k-mer as deposited, the k-mers of the
// graph (BASIC: the records'; CANONICAL: those and their reverse complements), the contexts of a
// pattern, and a pattern's motif union
struct Oracle {
    size_t k = 0;
    bool canonical = false;
    std::map<std::string, std::set<std::string>> cols;
    std::set<std::string> kmers;

    explicit Oracle(const Index &idx)
          : k(idx.k), canonical(idx.mode == DeBruijnGraph::CANONICAL) {
        for (const Record &r : idx.records) {
            for (size_t s = 0; s + k <= r.seq.size(); ++s) {
                cols[r.seq.substr(s, k)].insert(r.column);
                kmers.insert(r.seq.substr(s, k));
                if (canonical)
                    kmers.insert(rc(r.seq.substr(s, k)));
            }
        }
    }
    std::set<std::string> own(const std::string &kmer) const {
        auto it = cols.find(kmer);
        return it == cols.end() ? std::set<std::string>() : it->second;
    }
    std::set<std::string> columns() const {
        std::set<std::string> out;
        for (const auto &[kmer, c] : cols) {
            out.insert(c.begin(), c.end());
        }
        return out;
    }
    std::set<CtxKey> contexts(const std::string &p, gp::Scope scope, gp::Strands strands) const {
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
        std::set<CtxKey> out;
        for (const std::string &kmer : kmers) {
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
    // the rows a context of |kmer| is evaluated on, restricted to |names|: per name the rows it
    // was found on (1 the k-mer's, 2 its reverse complement's). CANONICAL: one row for both
    // (flagged 1); BASIC "context": the k-mer's row only
    std::map<std::string, uint8_t> evaluated(const std::string &kmer, bool either,
                                             const std::set<std::string> &names) const {
        std::map<std::string, uint8_t> out;
        for (const std::string &c : own(kmer)) {
            if (names.count(c))
                out[c] |= 1;
        }
        if (either || canonical) {
            for (const std::string &c : own(rc(kmer))) {
                if (names.count(c))
                    out[c] |= canonical ? 1 : 2;
            }
        }
        return out;
    }
};

// The motif union of a set of contexts: per label the contexts carrying it and the rows (1, 2)
struct Union {
    std::map<std::string, uint64_t> contexts;
    std::map<std::string, uint8_t> on;
    std::set<std::string> names() const {
        std::set<std::string> out;
        for (const auto &[n, c] : contexts) {
            out.insert(n);
        }
        return out;
    }
};

template <class Kmers>
Union union_of(const Oracle &o, const Kmers &kmers, bool either,
               const std::set<std::string> &names) {
    Union u;
    for (const std::string &kmer : kmers) {
        for (const auto &[c, on] : o.evaluated(kmer, either, names)) {
            ++u.contexts[c];
            u.on[c] |= on;
        }
    }
    return u;
}

std::string strands_name(uint8_t on, bool canonical) {
    return canonical ? "either" : on == 3 ? "both" : on == 2 ? "reverse_complement" : "context";
}

// ---------------------------------------------------------------- the pass

struct Ask {
    std::string pattern = "AC";
    gp::PatternKind kind = gp::PatternKind::DNA;
    std::string predicate = "{\"any\": [\"c1\"]}";
    gp::Mode mode = gp::Mode::COUNT;
    gp::Scope scope = gp::Scope::ANY_OFFSET;
    gp::Strands strands = gp::Strands::BOTH;
    bool either = true;
    bool motif = true;
    uint64_t max_predicate_contexts = 100'000;
    uint64_t max_contexts = 10'000;
    bool stop_at_threshold = false;
    Projection projection = Projection::NONE;
    uint64_t max_predicate_work = 100'000'000;
    uint64_t max_steps = 1'000'000;
    // the account in bytes (0: the default 256 MiB), the decode's deny, a hook before each read
    uint64_t max_memory_bytes = 0;
    std::function<bool(uint64_t)> deny_decode;
    std::function<void(size_t)> read_hook;
    // a virtual clock (null: never expires)
    std::function<Clock::time_point()> clock;
};

struct Pass {
    SelectionAnswer answer;
    Json::Value motif;
    // the raw count, the engine's extraction, and the released contexts the pass tested
    gp::Count raw;
    gp::Extraction x;
    uint64_t released = 0;
    std::vector<TestedContext> tested;
    // per tested context its (strand, k-mer, offset)
    std::vector<CtxKey> keys;
    // the bound predicate's labels (names, by LabelId) and its constant, if any
    std::vector<std::string> labels;
    std::optional<bool> constant;
    // the binding stopped (the account could not hold the predicate, or the time passed)
    bool unbound = false;
    uint64_t memory_peak = 0;
};

const char* strand_of(gp::Orientation o) {
    return o == gp::Orientation::FORWARD ? "+" : o == gp::Orientation::REVERSE ? "-" : "=";
}

// one pattern's selection as the route runs it (SPEC §19.6 steps 1-4): the predicate bound, the
// engine's release into the pass (max_contexts = max_predicate_contexts, ALL_OR_COUNT for a
// count), each context's descriptor admitted, then select() with |r|.motif
Pass run(const Index &idx, const Ask &r) {
    const DeBruijnGraph &graph = idx.anno->get_graph();
    const gp::GraphSupport support = gp::PatternSearch::support(graph);
    EXPECT_TRUE(support.supported) << support.reason;
    gp::Budget budget(r.max_steps,
                      r.clock ? gp::Deadline(Clock::time_point(), 1000, 0, r.clock)
                              : test::unbounded_deadline());
    RetrievalLimits limits;
    RetrievalHooks hooks;
    hooks.max_memory_bytes = r.max_memory_bytes;
    hooks.deny_decode = r.deny_decode;
    hooks.read_hook = r.read_hook;
    PatternRetrieval retrieval(*idx.anno, support.mode, limits, budget, &hooks);
    SelectionLimits sl;
    sl.either = r.either;
    sl.max_predicate_work = r.max_predicate_work;
    Pass out;
    retrieval.bind(predicate::Predicate::parse(parse_pattern_body(r.predicate), 10'000), sl);
    const predicate::Bound *bound = retrieval.bound();
    if (!bound) {
        out.unbound = true;
        return out;
    }
    for (const auto &l : bound->labels()) {
        out.labels.push_back(l.name);
    }
    out.constant = bound->constant();
    if (out.constant)
        return out;

    const gp::PatternSearch search(graph);
    const gp::Pattern pattern = gp::Pattern::parse(r.kind, r.pattern);
    gp::Request request;
    request.mode = r.mode == gp::Mode::COUNT ? gp::Mode::ALL_OR_COUNT : r.mode;
    request.scope = r.scope;
    request.strands = r.strands;
    request.max_contexts = r.max_predicate_contexts;
    request.stop_at_threshold = r.stop_at_threshold;
    request.min_information_bits = 0;
    retrieval.begin_selection(r.mode);
    const gp::Result result = search.enumerate(pattern, request, budget,
                                               [&](const gp::Context &c) {
        ++out.released;
        if (!retrieval.admit_tested())
            return;
        out.tested.push_back(retrieval.tested_context(c));
        out.keys.emplace_back(strand_of(c.orientation), graph.get_node_sequence(c.node),
                              c.offset);
    });
    EXPECT_FALSE(result.refusal);
    EXPECT_TRUE(result.extraction);
    if (result.refusal || !result.extraction)
        return out;
    out.raw = result.contexts->total;
    out.x = *result.extraction;
    SelectionRequest sr;
    sr.mode = r.mode;
    sr.max_contexts = r.max_contexts;
    sr.stop_at_threshold = r.stop_at_threshold;
    sr.projection = r.projection;
    sr.motif = r.motif;
    out.answer = retrieval.select(out.tested, out.released, out.raw, out.x, sr);
    if (out.answer.motif)
        out.motif = retrieval.motif_json(*out.answer.motif);
    retrieval.end_selection();
    out.memory_peak = retrieval.memory_peak();
    return out;
}

// the labels of the request's predicate's normal form (the index's columns the folding keeps)
std::set<std::string> known_names(const Oracle &o, const std::string &predicate) {
    std::set<std::string> names;
    names_of(o_fold(parse_pattern_body(predicate), o.columns()), &names);
    return names;
}

/**
 * Every motif answer, against the oracles. |memory_cut|: the account may have refused a union
 * entry (the union then a subset of the decided contexts'). Returns the oracle's value over
 * every context of the pattern.
 */
bool check_motif(const Oracle &o, const Ask &r, const Pass &run, bool memory_cut = false) {
    const Json::Value pred = parse_pattern_body(r.predicate);
    const std::set<std::string> known = known_names(o, r.predicate);
    const std::set<CtxKey> contexts = o.contexts(r.pattern, r.scope, r.strands);
    std::vector<std::string> all_kmers;
    for (const auto &[strand, kmer, offset] : contexts) {
        all_kmers.push_back(kmer);
    }
    const Union truth = union_of(o, all_kmers, r.either, known);
    const bool value = o_eval(pred, truth.names());

    EXPECT_TRUE(run.answer.motif);
    if (!run.answer.motif)
        return value;
    const MotifAnswer &m = *run.answer.motif;
    const Json::Value &j = run.motif;
    EXPECT_EQ(known.size(), m.labels);
    EXPECT_EQ(known.size(), j["labels"].asUInt64());
    std::set<std::string> label_set(run.labels.begin(), run.labels.end());
    EXPECT_EQ(known, label_set) << "the normal form's labels";

    // the released contexts are the pattern's (an independent scan), each tested one among them
    for (const CtxKey &key : run.keys) {
        EXPECT_TRUE(contexts.count(key)) << std::get<1>(key);
    }
    const bool completed = run.answer.pass == SelectionPass::COMPLETED;
    if (completed) {
        EXPECT_EQ(contexts.size(), run.keys.size());
        EXPECT_EQ(gp::Relation::EXACT, run.answer.tested.relation);
        EXPECT_EQ(contexts.size(), run.answer.tested.value);
    }

    // a value is stated only when it is the motif's: definite values are never contradicted
    if (m.value) {
        EXPECT_EQ(value, *m.value) << "a definite motif value is the oracle's";
        EXPECT_NE(MotifBasis::UNDECIDED, m.basis);
        EXPECT_EQ(Json::Value(*m.value), j["selected"]);
    } else {
        EXPECT_EQ(MotifBasis::UNDECIDED, m.basis);
        EXPECT_TRUE(j["selected"].isNull());
        EXPECT_TRUE(j["decided_by"].isNull());
    }
    EXPECT_EQ(completed && !m.stop, m.basis == MotifBasis::EVERY_CONTEXT);
    EXPECT_EQ(completed, m.untested == MotifUntested::NONE);
    EXPECT_EQ(completed, j["untested"].isNull());
    if (m.stop) {
        EXPECT_FALSE(m.value);
        EXPECT_EQ(m.stop, j["stop"].asString());
    } else {
        EXPECT_TRUE(j["stop"].isNull());
    }
    if (m.stop && std::string(m.stop) == "time") {
        EXPECT_FALSE(m.present);
        EXPECT_TRUE(j["labels_present"].isNull());
        EXPECT_TRUE(run.answer.time_limited);
        return value;
    }
    EXPECT_TRUE(m.present);
    if (!m.present)
        return value;

    // the labels found: a subset of the motif's union, every one on a decided context
    std::vector<std::string> decided_kmers;
    for (size_t i = 0; i < run.tested.size(); ++i) {
        if (run.tested[i].decided)
            decided_kmers.push_back(std::get<1>(run.keys[i]));
    }
    const Union seen = union_of(o, decided_kmers, r.either, known);
    std::set<std::string> found;
    EXPECT_EQ(m.present->size(), j["labels_present"].size());
    if (m.present->size() != j["labels_present"].size())
        return value;
    std::map<std::string, std::string> json_strands;
    for (size_t i = 0; i < m.present->size(); ++i) {
        const MotifLabel &l = (*m.present)[i];
        const std::string name = run.labels.at(l.label);
        const Json::Value &e = j["labels_present"][static_cast<Json::ArrayIndex>(i)];
        EXPECT_EQ(name, e["column"].asString());
        found.insert(name);
        json_strands[name] = e["strands"].asString();
        EXPECT_TRUE(truth.contexts.count(name)) << name << " is not on any context";
        EXPECT_TRUE(seen.contexts.count(name)) << name << " is not on a decided context";
        EXPECT_LE(l.contexts, seen.contexts.count(name) ? seen.contexts.at(name) : 0) << name;
        EXPECT_GE(l.contexts, 1u);
        EXPECT_EQ(l.contexts, e["contexts"]["value"].asUInt64());
        EXPECT_EQ(m.basis == MotifBasis::EVERY_CONTEXT ? "exact" : "at_least",
                  e["contexts"]["relation"].asString());
        // its first context carries it
        EXPECT_LT(l.first, run.tested.size());
        if (l.first >= run.tested.size())
            continue;
        EXPECT_TRUE(run.tested[l.first].decided);
        EXPECT_TRUE(o.evaluated(std::get<1>(run.keys[l.first]), r.either, known).count(name));
        if (!memory_cut) {
            // every decided context joined the union: its counts and rows are theirs exactly
            EXPECT_EQ(seen.contexts.at(name), l.contexts) << name;
            EXPECT_EQ(strands_name(seen.on.at(name), o.canonical), e["strands"].asString())
                    << name;
            size_t first = 0;
            while (first < run.tested.size()
                    && !(run.tested[first].decided
                         && o.evaluated(std::get<1>(run.keys[first]), r.either, known)
                                 .count(name))) {
                ++first;
            }
            EXPECT_EQ(first, l.first) << name;
        }
        // in label order: contexts desc, column asc
        if (i) {
            const MotifLabel &p = (*m.present)[i - 1];
            EXPECT_TRUE(p.contexts > l.contexts
                        || (p.contexts == l.contexts && run.labels.at(p.label) < name));
        }
    }
    if (!memory_cut) {
        EXPECT_EQ(seen.names(), found) << "the union is the decided contexts' labels";
    }

    if (m.basis == MotifBasis::EVERY_CONTEXT) {
        EXPECT_EQ(truth.names(), found);
        for (const auto &[name, n] : truth.contexts) {
            EXPECT_EQ(strands_name(truth.on.at(name), o.canonical), json_strands[name]) << name;
            for (const MotifLabel &l : *m.present) {
                if (run.labels.at(l.label) == name) {
                    EXPECT_EQ(n, l.contexts) << name;
                }
            }
        }
        EXPECT_EQ(known.size() - found.size(), j["labels_absent"].asUInt64());
        EXPECT_EQ("every_context", j["decided_by"].asString());
    } else {
        EXPECT_TRUE(j["labels_absent"].isNull());
        // three-valued on what was found: the oracle's strong Kleene over the same sets
        std::set<std::string> maybe;
        for (const std::string &n : known) {
            if (!found.count(n))
                maybe.insert(n);
        }
        if (!m.stop) {
            EXPECT_EQ(o_kleene(pred, found, maybe), m.value);
        }
        if (m.value) {
            EXPECT_EQ(MotifBasis::TESTED_CONTEXTS, m.basis);
            EXPECT_EQ("tested_contexts", j["decided_by"].asString());
        }
    }
    return value;
}

// predicates over |names| (seeded): the motif question and the other operators, unknown names
// among them
std::string random_predicate(std::mt19937 &rng, const std::vector<std::string> &names) {
    auto pick = [&](size_t n) {
        std::vector<std::string> out;
        std::set<std::string> used;
        while (out.size() < n) {
            const std::string s = rng() % 9 == 0 ? std::string("zz") : names[rng() % names.size()];
            if (used.insert(s).second)
                out.push_back(s);
        }
        Json::Value list(Json::arrayValue);
        for (const std::string &s : out) {
            list.append(s);
        }
        return list;
    };
    Json::Value p;
    switch (rng() % 7) {
        case 0: p["any"] = pick(1 + rng() % 2); break;
        case 1: p["none"] = pick(1 + rng() % 3); break;
        case 2: p["all"] = pick(2); break;
        case 3: {
            Json::Value a, c;
            a["any"] = pick(1 + rng() % 2);
            c["none"] = pick(1 + rng() % 3);
            p["and"].append(a);
            p["and"].append(c);
            break;
        }
        case 4: {
            const Json::Value list = pick(3);
            p["at_least"]["n"] = 2;
            p["at_least"]["labels"] = list;
            break;
        }
        case 5: {
            Json::Value a, n, q;
            a["all"] = pick(2);
            q["any"] = pick(1);
            n["not"] = q;
            p["or"].append(a);
            p["or"].append(n);
            break;
        }
        default: {
            Json::Value a, c;
            a["all"] = pick(2);
            c["none"] = pick(1);
            p["and"].append(a);
            p["and"].append(c);
        }
    }
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, p);
}

// ---------------------------------------------------------------- the tests

// The design's example (DESIGN §5.6): with k = 5 and the pattern AC, a sample A holding TACGG
// and a sample C holding CACCC. "any(A) and none(C)" selects A's context per context, but the
// motif is in C too: false at the motif level, decided on every context, both labels found.
// Mutation tried: the motif value as "some context selected": fails.
TEST(PatternMotif, MotifIsNotTheContextSelection) {
    const Index idx = build(5, { { "A", "A_r0", "TACGG" }, { "C", "C_r0", "CACCC" } });
    const Oracle oracle(idx);
    Ask r;
    r.predicate = "{\"and\": [{\"any\": [\"A\"]}, {\"none\": [\"C\"]}]}";
    for (const gp::Mode mode : { gp::Mode::COUNT, gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
        r.mode = mode;
        const Pass o = run(idx, r);
        ASSERT_EQ(SelectionPass::COMPLETED, o.answer.pass);
        EXPECT_EQ(1u, o.answer.selected.value) << "per context: A's flank only";
        ASSERT_TRUE(o.answer.motif);
        EXPECT_EQ(false, o.answer.motif->value);
        EXPECT_EQ(MotifBasis::EVERY_CONTEXT, o.answer.motif->basis);
        EXPECT_EQ(false, check_motif(oracle, r, o));
        ASSERT_EQ(2u, o.motif["labels_present"].size());
        EXPECT_EQ(0u, o.motif["labels_absent"].asUInt64());
        EXPECT_EQ("every_context", o.motif["decided_by"].asString());
    }
    // "present in A, absent throughout D" (D not a column of this index: folded to any(A))
    r.predicate = "{\"and\": [{\"any\": [\"A\"]}, {\"none\": [\"D\"]}]}";
    const Pass o = run(idx, r);
    EXPECT_EQ(true, o.answer.motif->value);
    EXPECT_EQ(true, check_motif(oracle, r, o));
}

// Every pattern, scope, strand set and predicate_strands over random indexes: the motif value
// is the predicate on the union of the contexts' labels (both strands with "either"), every
// label's context count and strands the oracle's. Mutations tried: the mirror row left out of
// the union under "either"; a context counted once per row instead of once; the union taken
// over the selected contexts only; the label order by id: each fails.
TEST(PatternMotif, MotifIsThePredicateOnTheUnion) {
    uint64_t differs = 0, checked = 0;
    std::mt19937 rng(7);
    for (unsigned seed = 1; seed <= 6; ++seed) {
        const std::vector<Record> records = seed == 1 ? kR1 : random_records(8, seed);
        const Index idx = build(kK, records);
        const Oracle oracle(idx);
        std::vector<std::string> names;
        for (const std::string &c : oracle.columns()) {
            names.push_back(c);
        }
        for (const std::string pattern : { "AC", "GT", "TA", "ACG", "CNG", "RY", "GGAC" }) {
            for (int t = 0; t < 6; ++t) {
                Ask r;
                r.pattern = pattern;
                r.kind = pattern.find_first_not_of("ACGT") == std::string::npos
                        ? gp::PatternKind::DNA : gp::PatternKind::IUPAC;
                r.predicate = random_predicate(rng, names);
                r.either = rng() % 2;
                r.scope = rng() % 4 ? gp::Scope::ANY_OFFSET : gp::Scope::SUFFIX;
                r.strands = std::vector<gp::Strands>{ gp::Strands::BOTH, gp::Strands::FORWARD,
                                                      gp::Strands::REVERSE }[rng() % 3];
                r.mode = std::vector<gp::Mode>{ gp::Mode::COUNT, gp::Mode::ALL_OR_COUNT,
                                                gp::Mode::PARTIAL }[rng() % 3];
                r.projection = r.mode == gp::Mode::COUNT ? Projection::NONE
                             : rng() % 2 ? Projection::PREDICATE_ONLY : Projection::NONE;
                const Pass o = run(idx, r);
                SCOPED_TRACE(pattern + " " + r.predicate + (r.either ? " either" : " context")
                             + " seed " + std::to_string(seed));
                if (o.constant) {
                    // folded to a constant: no pass (the route answers constant_motif)
                    continue;
                }
                ASSERT_EQ(SelectionPass::COMPLETED, o.answer.pass);
                ASSERT_EQ(MotifBasis::EVERY_CONTEXT, o.answer.motif->basis);
                const bool value = check_motif(oracle, r, o);
                ++checked;
                differs += value != (o.answer.selected.value > 0);
            }
        }
    }
    EXPECT_GT(checked, 150u);
    // the motif and the context-level selection do differ in this panel
    EXPECT_GT(differs, 0u);
}

// On a CANONICAL graph one row serves a k-mer and its reverse complement: the union is the
// canonical rows', every label's strands "either", whatever predicate_strands says.
TEST(PatternMotif, CanonicalRowsServeBothStrands) {
    const Index idx = build(kK, kR1, DeBruijnGraph::CANONICAL);
    const Oracle oracle(idx);
    for (const std::string predicate : { "{\"and\": [{\"any\": [\"c1\"]}, {\"none\": [\"c5\"]}]}",
                                         "{\"at_least\": {\"n\": 2, \"labels\": [\"c2\", "
                                         "\"c3\", \"c4\"]}}",
                                         "{\"none\": [\"c2\"]}" }) {
        for (const bool either : { false, true }) {
            for (const std::string pattern : { "AC", "TA", "GTC" }) {
                Ask r;
                r.pattern = pattern;
                r.predicate = predicate;
                r.either = either;
                const Pass o = run(idx, r);
                SCOPED_TRACE(pattern + " " + predicate);
                ASSERT_EQ(SelectionPass::COMPLETED, o.answer.pass);
                EXPECT_EQ(0u, o.answer.lookups);
                check_motif(oracle, r, o);
                for (const Json::Value &l : o.motif["labels_present"]) {
                    EXPECT_EQ("either", l["strands"].asString());
                }
            }
        }
    }
}

// A label holding the motif on the other strand only (c5 holds rc(c1/r0)): with "either" it
// carries the motif (strands "reverse_complement" where no context's own row has it), with
// "context" (the deposited strand only, BASIC) it does not, and none(c5) flips.
TEST(PatternMotif, OtherStrandOnlyWithEither) {
    const Index idx = build(kK, kR1);
    const Oracle oracle(idx);
    Ask r;
    r.pattern = "GGACG";
    r.strands = gp::Strands::FORWARD;
    r.predicate = "{\"none\": [\"c5\"]}";
    r.either = false;
    Pass o = run(idx, r);
    EXPECT_EQ(true, o.answer.motif->value);
    EXPECT_EQ(true, check_motif(oracle, r, o));
    r.either = true;
    o = run(idx, r);
    EXPECT_EQ(false, o.answer.motif->value);
    EXPECT_EQ(false, check_motif(oracle, r, o));
    ASSERT_EQ(1u, o.motif["labels_present"].size());
    EXPECT_EQ("c5", o.motif["labels_present"][0]["column"].asString());
    EXPECT_EQ("reverse_complement", o.motif["labels_present"][0]["strands"].asString());
}

// A palindromic k-mer (even k) is its own reverse complement: under "either" a label on its row
// carries the context on both strands ("both"), under "context" on the deposited one only. X
// holds only the palindrome GAATTC, Y only a flank of it on one strand.
TEST(PatternMotif, PalindromicKmerIsBothStrands) {
    const Index idx = build(6, { { "X", "X_r0", "GAATTC" }, { "Y", "Y_r0", "TGAATT" } });
    const Oracle oracle(idx);
    for (const bool either : { false, true }) {
        Ask r;
        r.pattern = "AATT";
        r.predicate = "{\"and\": [{\"any\": [\"X\"]}, {\"none\": [\"Y\"]}]}";
        r.either = either;
        const Pass o = run(idx, r);
        SCOPED_TRACE(either ? "either" : "context");
        ASSERT_EQ(SelectionPass::COMPLETED, o.answer.pass);
        EXPECT_EQ(false, check_motif(oracle, r, o));
        std::map<std::string, std::string> strands;
        for (const Json::Value &l : o.motif["labels_present"]) {
            strands[l["column"].asString()] = l["strands"].asString();
        }
        EXPECT_EQ(either ? "both" : "context", strands["X"]);
        EXPECT_EQ("context", strands["Y"]);
    }
}

/**
 * Every max_predicate_work from 1 to the complete pass's total + 1: the motif is exact once the
 * pass completes, and before that three-valued on the labels the decided contexts carry — a
 * definite value never contradicted, untested "selection". A label found decides any(A) true
 * and none(C) false while contexts are still untested (tested_contexts). Mutations tried: the
 * incomplete motif evaluated two-valued on what was found (an absence read as false): fails;
 * the incomplete union taken as every label undecided (no early decision): fails.
 */
TEST(PatternMotif, IncompletePassIsThreeValued) {
    const Index idx = build(kK, random_records(8, 11));
    const Oracle oracle(idx);
    const std::set<std::string> columns = oracle.columns();
    const std::vector<std::string> names(columns.begin(), columns.end());
    for (const bool either : { false, true }) {
        // a label on the first context in answer order (decided after the first rows read):
        // any(A) true and none(A) false before the other contexts are tested
        Ask probe;
        probe.pattern = "AC";
        probe.predicate = std::string("{\"any\": [\"") + names[0] + "\"]}";
        probe.either = either;
        const Pass first = run(idx, probe);
        ASSERT_GT(first.keys.size(), 2u);
        const auto on_first = oracle.evaluated(std::get<1>(first.keys[0]), either, columns);
        ASSERT_FALSE(on_first.empty());
        const std::string a = on_first.begin()->first;
        const std::vector<std::pair<std::string, std::optional<bool>>> predicates = {
            { std::string("{\"any\": [\"") + a + "\"]}", true },
            { std::string("{\"none\": [\"") + a + "\", \"" + (a == names[1] ? names[2] : names[1])
                      + "\"]}", false },
            { std::string("{\"and\": [{\"any\": [\"") + names[2] + "\"]}, {\"none\": [\""
                      + names[3] + "\", \"" + names[4] + "\"]}]}", std::nullopt },
            { std::string("{\"at_least\": {\"n\": 2, \"labels\": [\"") + names[0] + "\", \""
                      + names[5] + "\", \"" + names[6] + "\"]}}", std::nullopt },
        };
        for (const auto &[predicate, early_value] : predicates) {
            for (const gp::Mode mode : { gp::Mode::COUNT, gp::Mode::PARTIAL,
                                         gp::Mode::ALL_OR_COUNT }) {
                Ask r;
                r.pattern = "AC";
                r.predicate = predicate;
                r.either = either;
                r.mode = mode;
                const Pass full = run(idx, r);
                ASSERT_EQ(SelectionPass::COMPLETED, full.answer.pass);
                const bool truth = check_motif(oracle, r, full);
                const uint64_t total = full.answer.units;
                bool early = false;
                // every budget up to 64, then about 200 more up to the total + 1 (once a label
                // decides the value it does so at every larger budget until the pass completes)
                const uint64_t stride = std::max<uint64_t>(1, total / 200);
                for (uint64_t w = 1; w <= total + 1; w += w < 64 ? 1 : stride) {
                    r.max_predicate_work = w;
                    const Pass o = run(idx, r);
                    SCOPED_TRACE(predicate + (either ? " either" : " context") + " work "
                                 + std::to_string(w));
                    EXPECT_EQ(truth, check_motif(oracle, r, o));
                    if (o.answer.pass == SelectionPass::COMPLETED) {
                        EXPECT_EQ(MotifBasis::EVERY_CONTEXT, o.answer.motif->basis);
                        continue;
                    }
                    EXPECT_EQ(MotifUntested::SELECTION, o.answer.motif->untested);
                    EXPECT_EQ("selection", o.motif["untested"].asString());
                    early |= o.answer.motif->basis == MotifBasis::TESTED_CONTEXTS;
                }
                if (early_value) {
                    EXPECT_EQ(*early_value, truth);
                    EXPECT_TRUE(early) << predicate << " decided before every context";
                }
            }
        }
    }
}

// Why not every context was tested: the raw count above max_predicate_contexts (nothing read:
// not_admitted), partial's release cut there (release), the engine's step budget (discovery),
// rows the decode refused (selection); each undecided or decided soundly.
TEST(PatternMotif, UntestedContextsAreNamed) {
    const Index idx = build(kK, random_records(8, 5));
    const Oracle oracle(idx);
    const std::set<std::string> columns = oracle.columns();
    const std::vector<std::string> names(columns.begin(), columns.end());
    Ask base;
    base.pattern = "AC";
    base.predicate = std::string("{\"and\": [{\"any\": [\"") + names[0] + "\"]}, {\"none\": [\""
            + names[1] + "\"]}]}";
    const Pass full = run(idx, base);
    ASSERT_EQ(SelectionPass::COMPLETED, full.answer.pass);
    const uint64_t raw = full.raw.value;
    ASSERT_GT(raw, 4u);

    for (const gp::Mode mode : { gp::Mode::COUNT, gp::Mode::ALL_OR_COUNT }) {
        Ask r = base;
        r.mode = mode;
        r.max_predicate_contexts = raw - 1;
        const Pass o = run(idx, r);
        EXPECT_EQ(SelectionPass::NOT_ADMITTED, o.answer.pass);
        EXPECT_EQ(MotifUntested::NOT_ADMITTED, o.answer.motif->untested);
        EXPECT_EQ(MotifBasis::UNDECIDED, o.answer.motif->basis);
        EXPECT_EQ(0u, o.motif["labels_present"].size());
        check_motif(oracle, r, o);
    }
    {
        Ask r = base;
        r.mode = gp::Mode::PARTIAL;
        r.max_predicate_contexts = raw / 2;
        const Pass o = run(idx, r);
        EXPECT_EQ(SelectionPass::STOPPED, o.answer.pass);
        EXPECT_EQ(MotifUntested::RELEASE, o.answer.motif->untested);
        EXPECT_EQ("release", o.motif["untested"].asString());
        check_motif(oracle, r, o);
    }
    {
        // the engine's step budget: the raw count at_least
        bool seen = false;
        for (uint64_t steps = 1; steps < 400 && !seen; ++steps) {
            Ask r = base;
            r.mode = gp::Mode::PARTIAL;
            r.max_steps = steps;
            const Pass o = run(idx, r);
            if (o.raw.relation != gp::Relation::AT_LEAST || o.tested.empty())
                continue;
            seen = true;
            EXPECT_EQ(MotifUntested::DISCOVERY, o.answer.motif->untested);
            check_motif(oracle, r, o);
        }
        EXPECT_TRUE(seen);
    }
    for (const gp::Mode mode : { gp::Mode::COUNT, gp::Mode::PARTIAL }) {
        // every second read refused by the decode: rows_refused, the others decided
        Ask r = base;
        r.mode = mode;
        r.deny_decode = [](uint64_t ordinal) { return ordinal % 2 == 1; };
        const Pass o = run(idx, r);
        SCOPED_TRACE(std::string(gp::to_string(mode)));
        if (!o.answer.rows_refused.empty()) {
            EXPECT_EQ(SelectionPass::STOPPED, o.answer.pass);
            EXPECT_EQ(MotifUntested::SELECTION, o.answer.motif->untested);
        }
        check_motif(oracle, r, o);
    }
}

// The motif needs no read of its own: with it the pass reads the same rows, makes the same
// lookups, tests and selects the same contexts and lists the same ones; its units are the
// context-level pass's plus the final evaluation's. Without it the answer has no motif. A
// mutation reading a row again for the union, or charging it as a read: fails.
TEST(PatternMotif, SameReadsAsTheContextPass) {
    const Index idx = build(kK, random_records(10, 3));
    const Oracle oracle(idx);
    const std::set<std::string> columns = oracle.columns();
    const std::vector<std::string> names(columns.begin(), columns.end());
    std::mt19937 rng(3);
    for (const std::string pattern : { "AC", "CNG", "TA" }) {
        for (int t = 0; t < 8; ++t) {
            Ask r;
            r.pattern = pattern;
            r.kind = pattern == "CNG" ? gp::PatternKind::IUPAC : gp::PatternKind::DNA;
            r.predicate = random_predicate(rng, names);
            r.either = t % 2;
            r.mode = t % 3 == 0 ? gp::Mode::PARTIAL : gp::Mode::ALL_OR_COUNT;
            r.projection = t % 4 < 2 ? Projection::PREDICATE_ONLY : Projection::NONE;
            r.max_contexts = 5;
            auto reads_on = std::make_shared<uint64_t>(0), reads_off = std::make_shared<uint64_t>(0);
            r.read_hook = [reads_on](size_t) { ++*reads_on; };
            const Pass on = run(idx, r);
            r.motif = false;
            r.read_hook = [reads_off](size_t) { ++*reads_off; };
            const Pass off = run(idx, r);
            SCOPED_TRACE(pattern + " " + r.predicate);
            if (on.constant)
                continue;
            EXPECT_FALSE(off.answer.motif);
            ASSERT_TRUE(on.answer.motif);
            EXPECT_EQ(*reads_off, *reads_on);
            EXPECT_EQ(off.answer.rows, on.answer.rows);
            EXPECT_EQ(off.answer.lookups, on.answer.lookups);
            EXPECT_EQ(off.answer.pass, on.answer.pass);
            EXPECT_EQ(off.answer.tested.value, on.answer.tested.value);
            EXPECT_EQ(off.answer.selected.value, on.answer.selected.value);
            EXPECT_EQ(off.answer.chosen, on.answer.chosen);
            EXPECT_EQ(off.answer.selection_labels, on.answer.selection_labels);
            EXPECT_EQ(off.answer.withheld, on.answer.withheld);
            EXPECT_EQ(off.answer.stop, on.answer.stop);
            EXPECT_EQ(off.answer.units + on.answer.motif->units, on.answer.units);
            EXPECT_GE(on.answer.motif->units, 1u);
            EXPECT_GE(on.memory_peak, off.memory_peak);
            r.motif = true;
            check_motif(oracle, r, on);
        }
    }
}

// The union's entries are charged to the request's account before they are made: at every
// account size the motif is sound (the union a subset of the decided contexts'), and some size
// that completes the pass without the motif stops it with the motif ({selection, max_memory}).
// Mutation tried: the entries not charged: the last expectation fails.
TEST(PatternMotif, AccountHoldsTheUnion) {
    // many labels on one motif: 40 columns of one random record each, all with AC
    std::vector<Record> records;
    std::mt19937 rng(17);
    for (size_t c = 0; c < 40; ++c) {
        std::string s(16, 'A');
        for (char &x : s) {
            x = "ACGT"[rng() % 4];
        }
        s.replace(4, 2, "AC");
        records.push_back({ "label_with_a_long_name_" + std::to_string(c), "h" + std::to_string(c),
                            s });
    }
    const Index idx = build(kK, records);
    const Oracle oracle(idx);
    Json::Value list(Json::arrayValue);
    for (size_t c = 0; c < 40; ++c) {
        list.append("label_with_a_long_name_" + std::to_string(c));
    }
    Json::Value p;
    p["none"] = list;
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    Ask r;
    r.pattern = "AC";
    r.predicate = Json::writeString(b, p);
    const Pass full = run(idx, r);
    ASSERT_EQ(SelectionPass::COMPLETED, full.answer.pass);
    ASSERT_EQ(40u, full.answer.motif->present->size());
    bool stopped_by_motif = false;
    uint64_t motif_stops = 0, motif_stops_stated = 0;
    for (uint64_t bytes = 2'000; bytes <= 200'000; bytes += 1'000) {
        r.max_memory_bytes = bytes;
        r.motif = false;
        const Pass off = run(idx, r);
        r.motif = true;
        const Pass on = run(idx, r);
        SCOPED_TRACE("bytes " + std::to_string(bytes));
        EXPECT_EQ(off.unbound, on.unbound);
        if (on.unbound)
            continue;
        check_motif(oracle, r, on, true);
        if (on.answer.motif && on.answer.motif->stop
                && std::string(on.answer.motif->stop) == "max_memory") {
            // the motif's own stop is the entry's too, unless an earlier one is
            EXPECT_TRUE(on.answer.stop.has_value());
            ++motif_stops;
            if (on.answer.stop && *on.answer.stop
                    == std::make_pair(std::string("output"), std::string("max_memory"))) {
                ++motif_stops_stated;
            }
        }
        if (off.answer.pass == SelectionPass::COMPLETED
                && on.answer.pass != SelectionPass::COMPLETED) {
            // the union's entries took what a row or the next entry needed
            EXPECT_TRUE((on.answer.stop && *on.answer.stop
                            == std::make_pair(std::string("selection"),
                                              std::string("max_memory")))
                        || !on.answer.rows_refused.empty());
            stopped_by_motif = true;
        }
    }
    EXPECT_TRUE(stopped_by_motif);
    std::cerr << "motif stops " << motif_stops << " stated " << motif_stops_stated << std::endl;
}

// A virtual clock past the work time at each reading in turn: never an exception, a definite
// motif value always the oracle's; a stop of the motif's own evaluation (the clock read before
// it, after a pass that completed) leaves it undecided and unlisted, {output, time}. At the
// reading after the last the motif is exact. Mutation tried: no clock before the evaluation:
// fails.
TEST(PatternMotif, DeadlineAtEveryReading) {
    const Index idx = build(kK, kR1);
    const Oracle oracle(idx);
    for (const bool either : { false, true }) {
        auto readings = std::make_shared<uint64_t>(0);
        auto at = std::make_shared<uint64_t>(UINT64_MAX);
        Ask r;
        r.pattern = "AC";
        r.predicate = "{\"and\": [{\"any\": [\"c1\"]}, {\"none\": [\"c2\"]}]}";
        r.either = either;
        r.clock = [readings, at]() {
            ++*readings;
            return Clock::time_point() + (*readings >= *at ? std::chrono::hours(1)
                                                           : std::chrono::hours(0));
        };
        const Pass full = run(idx, r);
        ASSERT_EQ(SelectionPass::COMPLETED, full.answer.pass);
        const uint64_t total = *readings;
        bool motif_stopped = false;
        for (uint64_t n = 1; n <= total + 1; ++n) {
            *readings = 0;
            *at = n;
            Pass o;
            ASSERT_NO_THROW(o = run(idx, r));
            SCOPED_TRACE("reading " + std::to_string(n));
            if (o.unbound)
                continue;
            ASSERT_TRUE(o.answer.motif);
            if (o.tested.empty() && o.answer.pass != SelectionPass::COMPLETED) {
                // stopped before the pass: undecided
                EXPECT_FALSE(o.answer.motif->value);
                continue;
            }
            check_motif(oracle, r, o);
            if (n == total + 1) {
                EXPECT_EQ(MotifBasis::EVERY_CONTEXT, o.answer.motif->basis);
            }
            if (o.answer.pass == SelectionPass::COMPLETED && o.answer.motif->stop) {
                motif_stopped = true;
                EXPECT_EQ("time", std::string(o.answer.motif->stop));
                ASSERT_TRUE(o.answer.stop);
                EXPECT_EQ("output", o.answer.stop->first);
                EXPECT_EQ("time", o.answer.stop->second);
            }
        }
        EXPECT_TRUE(motif_stopped);
    }
}

// Without contexts the union is empty and complete: the value is the normal form's on the
// empty set (none(A) true, any(A) false); a constant normal form decides without a context; a
// pass that did not run is undecided. The motif block's shape.
TEST(PatternMotif, EmptyConstantAndWithoutPass) {
    const Index idx = build(kK, kR1);
    const Oracle oracle(idx);
    for (const auto &[predicate, value] : std::vector<std::pair<std::string, bool>>{
            { "{\"none\": [\"c1\"]}", true }, { "{\"any\": [\"c1\"]}", false } }) {
        for (const gp::Mode mode : { gp::Mode::COUNT, gp::Mode::ALL_OR_COUNT }) {
            Ask r;
            r.pattern = "CCCCCC";   // no k-mer of R1 holds it
            r.predicate = predicate;
            r.mode = mode;
            const Pass o = run(idx, r);
            ASSERT_EQ(gp::Relation::EXACT, o.raw.relation);
            ASSERT_EQ(0u, o.raw.value);
            EXPECT_EQ(value, o.answer.motif->value);
            EXPECT_EQ(MotifBasis::EVERY_CONTEXT, o.answer.motif->basis);
            EXPECT_EQ(value, check_motif(oracle, r, o));
            EXPECT_EQ(1u, o.motif["labels_absent"].asUInt64());
        }
    }
    const MotifAnswer t = constant_motif(true);
    EXPECT_EQ(true, t.value);
    EXPECT_EQ(MotifBasis::CONSTANT, t.basis);
    EXPECT_EQ(MotifUntested::NONE, t.untested);
    ASSERT_TRUE(t.present);
    EXPECT_TRUE(t.present->empty());
    EXPECT_STREQ("constant", to_string(t.basis));

    gp::Budget budget(1'000'000, test::unbounded_deadline());
    PatternRetrieval retrieval(*idx.anno, gp::GraphMode::BASIC, RetrievalLimits(), budget);
    retrieval.bind(predicate::Predicate::parse(parse_pattern_body("{\"none\": [\"c1\"]}"), 10),
                   SelectionLimits());
    const MotifAnswer n = retrieval.motif_without_pass(
            SelectionPass::NOT_STARTED, gp::Count::at_least(gp::Unit::GRAPH_CONTEXTS, 3));
    EXPECT_FALSE(n.value);
    EXPECT_EQ(MotifUntested::NOT_STARTED, n.untested);
    const Json::Value j = retrieval.motif_json(n);
    EXPECT_TRUE(j["selected"].isNull());
    EXPECT_TRUE(j["decided_by"].isNull());
    EXPECT_EQ("not_started", j["untested"].asString());
    EXPECT_EQ(1u, j["labels"].asUInt64());
    EXPECT_EQ(0u, j["labels_present"].size());
    EXPECT_TRUE(j["labels_absent"].isNull());
    EXPECT_TRUE(j["stop"].isNull());
    const MotifAnswer z = retrieval.motif_without_pass(
            SelectionPass::NOT_ADMITTED, gp::Count::exact(gp::Unit::GRAPH_CONTEXTS, 0));
    EXPECT_EQ(true, z.value);
    EXPECT_EQ(MotifBasis::EVERY_CONTEXT, z.basis);
    EXPECT_THROW(retrieval.motif_without_pass(SelectionPass::COMPLETED,
                                              gp::Count::exact(gp::Unit::GRAPH_CONTEXTS, 0)),
                 std::logic_error);
    const Json::Value c = retrieval.motif_json(constant_motif(false));
    EXPECT_EQ(false, c["selected"].asBool());
    EXPECT_EQ("constant", c["decided_by"].asString());
    EXPECT_EQ(0u, c["labels_absent"].asUInt64());
}

} // namespace
