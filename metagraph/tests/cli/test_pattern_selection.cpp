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


// The selection pass of a predicate for patterns of L <= k (src/cli/pattern_selection.cpp,
// increment 5b-3; SPEC-DRAFT §19.5-§19.9, TESTS §3): which raw contexts of a pattern satisfy
// the bound predicate, read from their annotation rows, under max_predicate_work, the memory
// account and the deadline, with the relations of §19.7, the selection admission and the
// projection predicate_only (retrieve_given). Driven by a TEST-ONLY ENTRY POINT (run_all
// below) that runs a pattern as PLAN §2.1's pipeline does (the route wires the real one in
// 5b-4): the engine's release into the pass, select(), the chosen contexts' projection.
// The expectations come from ORACLES over the records (never the module under test): the
// columns whose records hold a k-mer as deposited (and its reverse complement for "either"),
// the graph contexts of a pattern scanned from the records' k-mers, and a recursive
// evaluator of the request's predicate JSON (unknown names are simply absent from a set,
// which is what the folding of §19.4 computes).

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

// ---------------------------------------------------------------- the test-only entry point

struct Ask {
    std::string predicate = "{\"any\": [\"c1\"]}";
    std::vector<std::string> patterns = { "AC" };
    gp::PatternKind kind = gp::PatternKind::DNA;
    gp::Mode mode = gp::Mode::ALL_OR_COUNT;
    gp::Scope scope = gp::Scope::ANY_OFFSET;
    gp::Strands strands = gp::Strands::BOTH;
    uint64_t max_contexts = 10'000;
    uint64_t max_predicate_contexts = 100'000;
    bool stop_at_threshold = false;
    Projection projection = Projection::NONE;
    bool either = true;
    uint64_t max_predicate_work = 100'000'000;
    uint64_t max_steps = 1'000'000'000;
    bool allow_unbudgeted = false;
    bool records = true;
    // the account in bytes (0: the request's 256 MiB), the decode's deny, the read hook
    uint64_t max_memory_bytes = 0;
    std::function<bool(uint64_t)> deny_decode;
    std::function<void(size_t)> read_hook;
    // a virtual clock (null: no deadline): the work time is 1,000 ms from the epoch
    std::function<Clock::time_point()> clock;
    // the budget's abort predicate (Budget::set_abort; null: never asked)
    std::function<bool()> abort;
    // call select() also for a constant predicate (it must refuse)
    bool force_select = false;
};

struct Outcome {
    bool constant = false;
    bool select_threw = false;
    std::vector<std::string> unknown;
    const char *bind_stop = nullptr;
    uint64_t bound_bytes = 0;
    std::optional<gp::Result> result;
    uint64_t released = 0;
    std::vector<TestedContext> tested;
    // the tested contexts' k-mers (spelled by the test, for the oracles)
    std::vector<std::string> kmers;
    SelectionAnswer sel;
    std::optional<LabelsAnswer> labels;
    std::vector<std::vector<std::string>> selection_labels;
    std::string strands, access;
    // the request's units after this pattern; the account after end_selection
    uint64_t predicate_units = 0;
    uint64_t held = 0;
    uint64_t peak = 0;
    uint64_t length = 0;
};

// PLAN §2.1 per pattern: bind once; the engine's release into the pass (admit_tested before
// each descriptor, mode ALL_OR_COUNT for count, max_contexts = max_predicate_contexts);
// select(); the chosen contexts' projection (admit_context, then retrieve_given or
// retrieve); end_selection()
std::vector<Outcome> run_all(const Index &idx, const Ask &r) {
    const DeBruijnGraph &graph = idx.anno->get_graph();
    const gp::GraphSupport support = gp::PatternSearch::support(graph);
    gp::Budget budget(r.max_steps, r.clock ? gp::Deadline(Clock::time_point(), 1000, 0, r.clock)
                                           : gp::Deadline::unbounded());
    if (r.abort)
        budget.set_abort(r.abort);
    RetrievalLimits limits;
    limits.allow_unbudgeted = r.allow_unbudgeted;
    RetrievalHooks hooks;
    if (r.records)
        hooks.coord_to_header = idx.cth.get();
    hooks.max_memory_bytes = r.max_memory_bytes;
    hooks.deny_decode = r.deny_decode;
    hooks.read_hook = r.read_hook;
    PatternRetrieval retrieval(*idx.anno, support.mode, limits, budget, &hooks);
    const predicate::Predicate p = predicate::Predicate::parse(parse_pattern_body(r.predicate),
                                                               10'000);
    SelectionLimits sl;
    sl.max_predicate_work = r.max_predicate_work;
    sl.either = r.either;
    const predicate::Binding &binding = retrieval.bind(p, sl);
    const gp::PatternSearch search(graph);

    std::vector<Outcome> out;
    for (const std::string &text : r.patterns) {
        Outcome o;
        o.bind_stop = binding.stop;
        o.bound_bytes = binding.admitted;
        const predicate::Bound *bound = retrieval.bound();
        if (bound) {
            o.constant = bound->constant().has_value();
            o.unknown = bound->unknown_labels();
        }
        const gp::Pattern pattern = gp::Pattern::parse(r.kind, text);
        o.length = pattern.length();
        gp::Request rq;
        rq.mode = r.mode == gp::Mode::COUNT ? gp::Mode::ALL_OR_COUNT : r.mode;
        rq.scope = r.scope;
        rq.strands = r.strands;
        rq.max_contexts = r.max_predicate_contexts;
        rq.stop_at_threshold = r.stop_at_threshold;
        rq.min_information_bits = 0;
        SelectionRequest sr;
        sr.mode = r.mode;
        sr.max_contexts = r.max_contexts;
        sr.stop_at_threshold = r.stop_at_threshold;
        sr.projection = r.projection;
        if (o.constant && !r.force_select) {
            // a constant: no pass (P18), the engine counts
            o.result = search.count(pattern, rq, budget);
            o.sel = constant_selection(*bound->constant(), o.result->contexts->total);
            out.push_back(std::move(o));
            continue;
        }
        retrieval.begin_selection(r.mode);
        o.result = search.enumerate(pattern, rq, budget, [&](const gp::Context &c) {
            ++o.released;
            if (!retrieval.admit_tested())
                return;
            o.tested.push_back(retrieval.tested_context(c));
        });
        for (const TestedContext &t : o.tested) {
            o.kmers.push_back(graph.get_node_sequence(t.node));
        }
        if (o.result->refusal) {
            out.push_back(std::move(o));
            continue;
        }
        try {
            o.sel = retrieval.select(o.tested, o.released, o.result->contexts->total,
                                     *o.result->extraction, sr);
        } catch (const std::logic_error &) {
            o.select_threw = true;
            retrieval.end_selection();
            out.push_back(std::move(o));
            continue;
        }
        for (size_t j = 0; j < o.sel.selection_labels.size(); ++j) {
            std::vector<std::string> names;
            for (const Json::Value &n : retrieval.selection_labels_json(o.sel, j)) {
                names.push_back(n.asString());
            }
            o.selection_labels.push_back(std::move(names));
        }
        if (r.projection != Projection::NONE && r.mode != gp::Mode::COUNT
                && !o.sel.chosen.empty()) {
            retrieval.begin_release(r.mode);
            std::vector<RetrievalContext> collected;
            for (uint32_t j : o.sel.chosen) {
                if (!retrieval.admit_context())
                    break;
                RetrievalContext c;
                c.orientation = o.tested[j].orientation;
                c.offset = o.tested[j].offset;
                c.kmer = o.kmers[j];
                c.key = o.tested[j].key;
                collected.push_back(std::move(c));
            }
            gp::Extraction x;
            x.returned = collected.size();
            x.complete = o.sel.selected.relation == gp::Relation::EXACT
                            && o.sel.chosen.size() == o.sel.selected.value;
            o.labels = r.projection == Projection::PREDICATE_ONLY
                    ? retrieval.retrieve_given(collected, o.sel.chosen.size(), o.length, r.mode,
                                               x, Json::Value())
                    : retrieval.retrieve(collected, o.sel.chosen.size(), o.length, r.mode, x,
                                         Json::Value());
        }
        retrieval.end_selection();
        o.predicate_units = retrieval.predicate_units();
        o.held = retrieval.memory_held();
        o.peak = retrieval.memory_peak();
        o.strands = retrieval.selection_strands();
        o.access = retrieval.selection_access();
        out.push_back(std::move(o));
    }
    return out;
}

Outcome run(const Index &idx, const Ask &r) {
    std::vector<Outcome> all = run_all(idx, r);
    return std::move(all.at(0));
}

std::string strand_of(const TestedContext &c) {
    return gp::strand_symbol(c.orientation);
}

// the oracle's decision of tested context |i|
bool want(const Oracle &oracle, const Json::Value &pred, const Outcome &o, size_t i,
          bool either) {
    return o_eval(pred, oracle.evaluated(o.kmers[i], either));
}

bool contains(const gp::Count &c, uint64_t truth) {
    switch (c.relation) {
        case gp::Relation::EXACT: return c.value == truth;
        case gp::Relation::AT_LEAST: return c.value <= truth;
        case gp::Relation::BOUNDS: return c.lower <= truth && truth <= c.upper;
        case gp::Relation::UNKNOWN: return true;
    }
    return false;
}

/**
 * The invariants of TESTS §1 on one outcome: every decided context decided as the oracle
 * says; tested exact = the decisions; selected contains |truth| (the oracle's selected among
 * all raw contexts); the chosen are selected, in answer order, the first of them; their
 * selection_labels the set they were evaluated on (restricted to the predicate's names), in
 * label order; a withheld answer lists nothing.
 */
void check_invariants(const Oracle &oracle, const Json::Value &pred, const Outcome &o,
                      bool either, std::optional<uint64_t> truth, const Ask &r) {
    // the names the pass reads: the normal form's (a known name the folding dropped is not)
    std::set<std::string> names;
    names_of(o_fold(pred, oracle.columns()), &names);
    uint64_t decided = 0;
    std::vector<uint32_t> selected;
    for (size_t i = 0; i < o.tested.size(); ++i) {
        if (!o.tested[i].decided) {
            EXPECT_FALSE(o.tested[i].selected);
            continue;
        }
        ++decided;
        EXPECT_EQ(want(oracle, pred, o, i, either), o.tested[i].selected)
                << o.kmers[i] << " " << strand_of(o.tested[i]) << o.tested[i].offset;
        if (o.tested[i].selected)
            selected.push_back(i);
    }
    if (o.sel.pass == SelectionPass::COMPLETED || o.sel.pass == SelectionPass::STOPPED) {
        EXPECT_EQ(gp::Relation::EXACT, o.sel.tested.relation);
        EXPECT_EQ(decided, o.sel.tested.value);
        EXPECT_EQ(selected.size(), o.sel.selected.value) << "S counts the decided selected";
    }
    if (truth) {
        EXPECT_TRUE(contains(o.sel.selected, *truth)) << *truth;
    }
    if (o.sel.pass == SelectionPass::COMPLETED) {
        EXPECT_EQ(gp::Relation::EXACT, o.sel.selected.relation);
        EXPECT_EQ(o.tested.size(), decided);
        EXPECT_FALSE(o.sel.stop);
    }
    // the chosen: selected, ascending, the first of the decided selected
    EXPECT_TRUE(std::is_sorted(o.sel.chosen.begin(), o.sel.chosen.end()));
    if (o.sel.withheld || r.mode == gp::Mode::COUNT) {
        EXPECT_TRUE(o.sel.chosen.empty());
    }
    for (size_t j = 0; j < o.sel.chosen.size(); ++j) {
        ASSERT_LT(j, selected.size());
        EXPECT_EQ(selected[j], o.sel.chosen[j]);
    }
    if (r.mode == gp::Mode::PARTIAL && !o.sel.stop) {
        EXPECT_EQ(std::min<uint64_t>(selected.size(), r.max_contexts), o.sel.chosen.size());
    }
    if (r.mode == gp::Mode::ALL_OR_COUNT && o.sel.pass == SelectionPass::COMPLETED) {
        if (selected.size() > r.max_contexts) {
            EXPECT_EQ("selected_above_threshold", o.sel.withheld.value_or(""));
        } else {
            EXPECT_FALSE(o.sel.withheld);
            EXPECT_EQ(selected, o.sel.chosen);
        }
    }
    // selection_labels: the set it was evaluated on, the predicate's names only, label order
    if (r.projection == Projection::NONE || o.sel.chosen.empty()) {
        EXPECT_TRUE(o.selection_labels.empty());
        return;
    }
    ASSERT_EQ(o.sel.chosen.size(), o.selection_labels.size());
    std::map<std::string, uint64_t> count;
    for (size_t j = 0; j < o.sel.chosen.size(); ++j) {
        std::set<std::string> expected;
        for (const std::string &c : oracle.evaluated(o.kmers[o.sel.chosen[j]], either)) {
            if (names.count(c))
                expected.insert(c);
        }
        const std::set<std::string> got(o.selection_labels[j].begin(),
                                        o.selection_labels[j].end());
        EXPECT_EQ(expected, got) << o.kmers[o.sel.chosen[j]];
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

// the tested contexts as (strand, k-mer, offset)
std::set<std::tuple<std::string, std::string, uint64_t>> tested_set(const Outcome &o) {
    std::set<std::tuple<std::string, std::string, uint64_t>> s;
    for (size_t i = 0; i < o.tested.size(); ++i) {
        s.emplace(strand_of(o.tested[i]), o.kmers[i], o.tested[i].offset);
    }
    return s;
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

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

// the rows of the pass in their order of first appearance, from the tested k-mers (BASIC):
// a context's own k-mer, then (either) its reverse complement's when the graph holds it
std::vector<std::string> row_order(const Oracle &oracle, const Outcome &o, bool either) {
    std::vector<std::string> order;
    std::set<std::string> seen;
    for (const std::string &x : o.kmers) {
        if (seen.insert(x).second)
            order.push_back(x);
        const std::string y = rc(x);
        if (either && y != x && oracle.present(y) && seen.insert(y).second)
            order.push_back(y);
    }
    return order;
}


// ---------------------------------------------------------------- tests

// Mutations tried (5b-3 mut/): "either" ignored, the mirror's labels left out of the evaluated
// set, the label order by name only, rows not deduplicated by key, the rows' entries never
// released, the kept rows dropped, selection_labels of the own row only: each fails here.
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
        for (int t = 0; t < 50; ++t) {
            Json::Value pred = random_predicate(rng, universe);
            for (const bool either : { true, false }) {
                for (const Projection proj : { Projection::NONE, Projection::PREDICATE_ONLY,
                                               Projection::ALL }) {
                    Ask r;
                    r.predicate = compact(pred);
                    r.either = either;
                    r.projection = proj;
                    Outcome o = run(idx, r);
                    SCOPED_TRACE(r.predicate + (either ? " either" : " context"));
                    const Json::Value folded = o_fold(pred, oracle.columns());
                    EXPECT_EQ(folded.isBool(), o.constant) << compact(folded);
                    if (o.constant) {
                        // folded to a constant (unknown names): no pass, no read
                        const bool value = o_eval(pred, {});
                        EXPECT_EQ(SelectionPass::CONSTANT, o.sel.pass);
                        bool all_same = true;
                        for (const auto &[s, kmer, off] : truth_contexts) {
                            all_same &= o_eval(pred, oracle.evaluated(kmer, either)) == value;
                        }
                        EXPECT_TRUE(all_same) << "a constant decides every context alike";
                        continue;
                    }
                    EXPECT_EQ(truth_contexts, tested_set(o));
                    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
                    check_invariants(oracle, pred, o, either,
                                     oracle_truth(oracle, pred, "AC", either), r);
                    EXPECT_EQ(either ? "either" : "context", o.strands);
                    EXPECT_EQ("rows", o.access);
                    EXPECT_EQ(row_order(oracle, o, either).size(), o.sel.rows);
                    if (proj == Projection::NONE) {
                        // nothing of the pass is held after it but the bound predicate
                        EXPECT_EQ(o.bound_bytes, o.held);
                        continue;
                    }
                    if (o.sel.chosen.empty())
                        continue;
                    ASSERT_TRUE(o.labels);
                    ASSERT_EQ(o.sel.chosen.size(), o.labels->result_fields.size());
                    EXPECT_TRUE(o.labels->complete);
                    std::set<std::string> names;
                    names_of(o_fold(pred, oracle.columns()), &names);
                    for (size_t j = 0; j < o.sel.chosen.size(); ++j) {
                        const uint32_t i = o.sel.chosen[j];
                        const Json::Value &f = o.labels->result_fields[j];
                        std::set<std::string> expected;
                        for (const std::string &c : oracle.own(o.kmers[i])) {
                            if (proj == Projection::ALL || names.count(c))
                                expected.insert(c);
                        }
                        std::set<std::string> got;
                        for (const Json::Value &l : f["labels"]) {
                            const std::string column = l["column"].asString();
                            got.insert(column);
                            // placed as the record scan says
                            std::set<std::pair<uint64_t, uint64_t>> occ;
                            for (const Json::Value &e : l["occurrence_list"]) {
                                const std::string nt = e["nt_coords"].asString();
                                occ.emplace(e["seq_id"].asUInt64(),
                                            std::stoull(nt.substr(0, nt.find('-'))));
                                EXPECT_EQ(strand_of(o.tested[i]), e["strand"].asString());
                            }
                            EXPECT_EQ(oracle.occurrences(column, o.kmers[i],
                                                         o.tested[i].offset), occ)
                                    << column << " " << o.kmers[i];
                        }
                        EXPECT_EQ(expected, got) << o.kmers[i];
                        EXPECT_EQ(expected.size(), f["labels_total"].asUInt64());
                        EXPECT_EQ("complete", f["labels_status"].asString());
                    }
                }
            }
        }
    }
}

// The determinism of the pass: the same request twice, the same answer (mutation tried: the
// kept rows dropped, retrieve_given refuses: fails)
TEST(PatternSelection, Deterministic) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Ask r;
    r.predicate = "{\"or\": [{\"any\": [\"c1\", \"c5\"]}, {\"none\": [\"c3\"]}]}";
    r.projection = Projection::PREDICATE_ONLY;
    r.mode = gp::Mode::PARTIAL;
    r.max_contexts = 3;
    Outcome a = run(idx, r), b = run(idx, r);
    EXPECT_EQ(a.sel.chosen, b.sel.chosen);
    EXPECT_EQ(a.selection_labels, b.selection_labels);
    EXPECT_EQ(a.sel.units, b.sel.units);
    EXPECT_EQ(a.sel.rows, b.sel.rows);
    EXPECT_EQ(a.peak, b.peak);
    ASSERT_TRUE(a.labels && b.labels);
    ASSERT_EQ(a.labels->result_fields.size(), b.labels->result_fields.size());
    for (size_t j = 0; j < a.labels->result_fields.size(); ++j) {
        EXPECT_EQ(a.labels->result_fields[j], b.labels->result_fields[j]);
    }
}

// P11: c5 holds GGACGACTTTG's motif on the - strand only. Mutation tried: the mirror's row
// not read (either as context): fails.
TEST(PatternSelection, StrandsEitherAndContextOnBasic) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    // none(c5): "context" selects the + contexts of c1/r0's k-mers (their own rows lack c5),
    // "either" none of them (their reverse complements are c5's)
    for (const bool either : { false, true }) {
        Ask r;
        r.predicate = "{\"none\": [\"c5\"]}";
        r.either = either;
        Outcome o = run(idx, r);
        const Json::Value pred = parse_pattern_body(r.predicate);
        check_invariants(oracle, pred, o, either, oracle_truth(oracle, pred, "AC", either), r);
        uint64_t c1_plus = 0, c1_plus_selected = 0;
        for (size_t i = 0; i < o.tested.size(); ++i) {
            if (strand_of(o.tested[i]) == "+" && oracle.own(o.kmers[i]).count("c1")
                    && oracle.own(rc(o.kmers[i])).count("c5")) {
                ++c1_plus;
                c1_plus_selected += o.tested[i].selected;
            }
        }
        ASSERT_GT(c1_plus, 0u);
        EXPECT_EQ(either ? 0u : c1_plus, c1_plus_selected);
    }
    // any(c5) with "either": the + contexts of c1/r0 are selected by their mirror's label,
    // which their selection_labels show and their own labels (predicate_only) do not
    Ask r;
    r.predicate = "{\"any\": [\"c5\"]}";
    r.projection = Projection::PREDICATE_ONLY;
    Outcome o = run(idx, r);
    const Json::Value pred = parse_pattern_body(r.predicate);
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r);
    bool seen = false;
    for (size_t j = 0; j < o.sel.chosen.size(); ++j) {
        const uint32_t i = o.sel.chosen[j];
        if (oracle.own(o.kmers[i]).count("c5"))
            continue;
        seen = true;
        EXPECT_EQ(std::vector<std::string>({ "c5" }), o.selection_labels[j]);
        ASSERT_TRUE(o.labels);
        EXPECT_EQ(0u, o.labels->result_fields[j]["labels"].size()) << o.kmers[i];
    }
    EXPECT_TRUE(seen);
}

// On CANONICAL and PRIMARY graphs one row serves both orientations: "either" whatever is
// asked, no lookup. Mutation tried: the lookups run on every graph mode: fails.
TEST(PatternSelection, CanonicalAndPrimaryShareTheRow) {
    for (const bool primary : { false, true }) {
        Index idx = primary
                ? build<annot::ColumnCompressed<>>(kK, kR1, false, DeBruijnGraph::PRIMARY)
                : build<annot::RowDiffColumnAnnotator>(kK, kR1, false, DeBruijnGraph::CANONICAL);
        Oracle oracle(idx);
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
                EXPECT_EQ(0u, o.sel.lookups);
                EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
                ASSERT_GT(o.tested.size(), 0u);
                check_invariants(oracle, parse_pattern_body(predicate), o, true, std::nullopt,
                                 r);
                EXPECT_EQ(primary ? "columns" : "rows", o.access);
            }
        }
    }
}

// NA: two contexts at one (k-mer, offset), one row read for both, decided alike. Mutation
// tried: rows not deduplicated by key: predicate_rows exceeds the distinct rows, fails.
TEST(PatternSelection, IupacTwoOrientationsOneRow) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    for (const bool either : { false, true }) {
        Ask r;
        r.patterns = { "NA" };
        r.kind = gp::PatternKind::IUPAC;
        r.predicate = "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"c4\", \"c5\"]}}";
        r.either = either;
        Outcome o = run(idx, r);
        const Json::Value pred = parse_pattern_body(r.predicate);
        EXPECT_EQ(oracle.contexts("NA"), tested_set(o));
        check_invariants(oracle, pred, o, either, std::nullopt, r);
        EXPECT_EQ(row_order(oracle, o, either).size(), o.sel.rows);
        std::map<std::pair<std::string, uint64_t>, std::set<bool>> at;
        for (size_t i = 0; i < o.tested.size(); ++i) {
            at[{ o.kmers[i], o.tested[i].offset }].insert(o.tested[i].selected);
        }
        bool shared = false;
        for (size_t i = 0; i + 1 < o.tested.size(); ++i) {
            shared |= o.kmers[i] == o.kmers[i + 1]
                        && o.tested[i].offset == o.tested[i + 1].offset;
        }
        EXPECT_TRUE(shared) << "NA has two contexts at one (k-mer, offset) here";
        for (const auto &[where, decisions] : at) {
            EXPECT_EQ(1u, decisions.size()) << where.first;
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
            check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r);
            if (n == cohort.size()) {
                Ask q = r;
                q.predicate = compact(all);
                Outcome a = run(idx, q);
                EXPECT_EQ(a.sel.chosen, o.sel.chosen);
                EXPECT_EQ(a.sel.selected.value, o.sel.selected.value);
            }
        }
    }
}

// A typo is an unknown name, folded away: and(any(c1), none(c1x)) selects as any(c1) and as
// the oracle says (mutations tried: "either" ignored, the mirror's labels dropped: fail)
TEST(PatternSelection, TypoReportedUnknown) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Ask r;
    r.predicate = "{\"and\": [{\"any\": [\"c1\"]}, {\"none\": [\"c1x\"]}]}";
    Outcome o = run(idx, r);
    EXPECT_EQ(std::vector<std::string>({ "c1x" }), o.unknown);
    Ask q = r;
    q.predicate = "{\"any\": [\"c1\"]}";
    Outcome p = run(idx, q);
    EXPECT_EQ(p.sel.chosen, o.sel.chosen);
    EXPECT_EQ(p.sel.selected.value, o.sel.selected.value);
    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body(r.predicate);
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r);
}

// P18: a predicate folded to a constant reads nothing; select() refuses it; the counts of
// constant_selection. Mutation tried: constant true answered selected exact 0: fails.
TEST(PatternSelection, AllUnknownNeedsNoRead) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const uint64_t raw = oracle.contexts("AC").size();
    for (const auto &[predicate, value] : std::vector<std::pair<std::string, bool>>{
            { "{\"any\": [\"zz\"]}", false }, { "{\"none\": [\"zz\"]}", true },
            { "{\"all\": [\"c1\", \"zz\"]}", false },
            { "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"zz\"]}}", false } }) {
        auto reads = std::make_shared<uint64_t>(0);
        Ask r;
        r.predicate = predicate;
        r.read_hook = [reads](size_t) { ++*reads; };
        Outcome o = run(idx, r);
        SCOPED_TRACE(predicate);
        EXPECT_TRUE(o.constant);
        EXPECT_EQ(SelectionPass::CONSTANT, o.sel.pass);
        EXPECT_EQ(0u, *reads);
        EXPECT_EQ(gp::Relation::EXACT, o.sel.tested.relation);
        EXPECT_EQ(raw, o.sel.tested.value);
        EXPECT_EQ(gp::Relation::EXACT, o.sel.selected.relation);
        EXPECT_EQ(value ? raw : 0u, o.sel.selected.value);
        r.force_select = true;
        EXPECT_TRUE(run(idx, r).select_threw);
    }
    // a constant false after a discovery stop: selected exact 0, tested the raw count's
    const gp::Count raw_stopped = gp::Count::at_least(gp::Unit::GRAPH_CONTEXTS, 3);
    SelectionAnswer f = constant_selection(false, raw_stopped);
    EXPECT_EQ(gp::Relation::EXACT, f.selected.relation);
    EXPECT_EQ(0u, f.selected.value);
    EXPECT_EQ(gp::Relation::AT_LEAST, f.tested.relation);
    SelectionAnswer t = constant_selection(true, raw_stopped);
    EXPECT_EQ(gp::Relation::AT_LEAST, t.selected.relation);
    EXPECT_EQ(3u, t.selected.value);
    // passes that did not run
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
 * in the order of first appearance (a prefix), more units never fewer of them; at the total +
 * 1 it completes. SPEC §19.13 (d): with "either" the first lookup comes first and no row is
 * read at a budget of 1. Mutations tried: the bounds' upper end without R - T (= S + R), the
 * final sweep removed (contexts after the cursor with read rows undecided), the gate checked
 * after the read (one row more), the refusal's units not charged: each fails here.
 */
TEST(PatternSelection, PassStoppedByItsBudget) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const std::string predicate = "{\"or\": [{\"all\": [\"c1\", \"c4\"]}, {\"any\": [\"c3\"]}]}";
    const Json::Value pred = parse_pattern_body(predicate);
    const uint64_t raw = oracle.contexts("AC").size();
    for (const bool either : { true, false }) {
        for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL,
                                     gp::Mode::COUNT }) {
            Ask r;
            r.predicate = predicate;
            r.either = either;
            r.mode = mode;
            const Outcome full = run(idx, r);
            ASSERT_EQ(SelectionPass::COMPLETED, full.sel.pass);
            const uint64_t total = full.sel.units;
            ASSERT_GT(total, 0u);
            const uint64_t truth = oracle_truth(oracle, pred, "AC", either);
            ASSERT_GT(full.sel.rows, 1u);
            // the budget at which exactly one row is read: every lookup (k units each) passes
            // its gate, the first read's gate passes (units < W), the second's does not
            const uint64_t one_row = either ? full.sel.lookups * kK + 1 : 1;
            uint64_t last_tested = 0;
            bool completed_before = false;
            for (uint64_t w = 1; w <= total + 1; ++w) {
                r.max_predicate_work = w;
                const Outcome o = run(idx, r);
                SCOPED_TRACE("work " + std::to_string(w) + (either ? " either " : " context ")
                             + gp::to_string(mode));
                check_invariants(oracle, pred, o, either, truth, r);
                EXPECT_GE(o.sel.tested.value, last_tested) << "monotone in the budget";
                last_tested = o.sel.tested.value;
                EXPECT_LE(o.sel.units, total);
                if (w == total + 1) {
                    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
                }
                if (w == one_row) {
                    EXPECT_EQ(SelectionPass::STOPPED, o.sel.pass);
                    EXPECT_EQ(1u, o.sel.rows) << "the gate is checked before a read";
                }
                // completed at a budget, completed at every larger one
                EXPECT_TRUE(!completed_before || o.sel.pass == SelectionPass::COMPLETED);
                completed_before |= o.sel.pass == SelectionPass::COMPLETED;
                if (o.sel.pass == SelectionPass::COMPLETED)
                    continue;
                ASSERT_EQ(SelectionPass::STOPPED, o.sel.pass);
                ASSERT_TRUE(o.sel.stop);
                EXPECT_EQ("selection", o.sel.stop->first);
                EXPECT_EQ("max_predicate_work", o.sel.stop->second);
                // it stopped at a gate: the units had reached the budget
                EXPECT_GE(o.sel.units, w);
                EXPECT_EQ(gp::Relation::BOUNDS, o.sel.selected.relation);
                EXPECT_EQ(o.sel.selected.lower + raw - o.sel.tested.value,
                          o.sel.selected.upper);
                // the decided contexts: those whose rows are among the first rows read
                const std::vector<std::string> order = row_order(oracle, o, either);
                ASSERT_LE(o.sel.rows, order.size());
                const std::set<std::string> read(order.begin(), order.begin() + o.sel.rows);
                for (size_t i = 0; i < o.tested.size(); ++i) {
                    const std::string y = rc(o.kmers[i]);
                    const bool ready = read.count(o.kmers[i])
                            && (!either || y == o.kmers[i] || !oracle.present(y)
                                || read.count(y));
                    EXPECT_EQ(ready, o.tested[i].decided) << o.kmers[i];
                }
                if (mode == gp::Mode::ALL_OR_COUNT) {
                    EXPECT_EQ("predicate_budget", o.sel.withheld.value_or(""));
                } else if (mode == gp::Mode::PARTIAL) {
                    EXPECT_EQ("max_predicate_work", o.sel.cut.value_or(""));
                    EXPECT_FALSE(o.sel.withheld);
                } else {
                    EXPECT_FALSE(o.sel.withheld);
                    EXPECT_FALSE(o.sel.cut);
                }
                if (w == 1 && either) {
                    // the first context's lookup (k units) comes first: no row read
                    EXPECT_EQ(0u, o.sel.rows);
                    EXPECT_EQ(0u, o.sel.tested.value);
                    EXPECT_EQ(kK, o.sel.units);
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
        for (size_t i = 0; i < o.tested.size(); ++i) {
            decisions += 1;
            for (const std::string &c : oracle.own(o.kmers[i])) {
                decisions += names.count(c);
            }
        }
        reads.push_back(o.sel.units - decisions);
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
        Outcome o = run(idx, r);
        std::set<std::string> known;
        uint64_t lookups = 0;
        for (const std::string &x : o.kmers) {
            if (known.count(x))
                continue;
            ++lookups;
            known.insert(x);
            if (oracle.present(rc(x)))
                known.insert(rc(x));
        }
        EXPECT_EQ(lookups, o.sel.lookups) << p;
        check_invariants(oracle, parse_pattern_body(r.predicate), o, true, std::nullopt, r);
    }
}

// A selective predicate on a pattern above max_contexts: the selected fit (all_or_count;
// mutation tried: the list one short, selected < max_contexts: fails)
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
    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
    EXPECT_FALSE(o.sel.withheld);
    EXPECT_EQ(truth, o.sel.chosen.size());
    check_invariants(oracle, pred, o, true, truth, r);
}

// The compute admission (§19.6 step 2): above max_predicate_contexts all_or_count and count
// read nothing (not_admitted), partial tests the first max_predicate_contexts (mutations
// tried: not_admitted answered not_started; the bounds' upper end S + R: fail)
TEST(PatternSelection, ComputeAdmission) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const uint64_t raw = oracle.contexts("AC").size();
    const Json::Value pred = parse_pattern_body("{\"any\": [\"c1\", \"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::COUNT, gp::Mode::PARTIAL }) {
        auto reads = std::make_shared<uint64_t>(0);
        Ask r;
        r.predicate = compact(pred);
        r.mode = mode;
        r.max_predicate_contexts = raw - 3;
        r.read_hook = [reads](size_t) { ++*reads; };
        Outcome o = run(idx, r);
        SCOPED_TRACE(gp::to_string(mode));
        check_invariants(oracle, pred, o, true, truth, r);
        if (mode != gp::Mode::PARTIAL) {
            EXPECT_EQ(SelectionPass::NOT_ADMITTED, o.sel.pass);
            EXPECT_EQ(0u, *reads);
            EXPECT_EQ(gp::Relation::UNKNOWN, o.sel.selected.relation);
            continue;
        }
        EXPECT_EQ(raw - 3, o.tested.size());
        EXPECT_EQ(raw - 3, o.sel.tested.value);
        EXPECT_EQ(SelectionPass::STOPPED, o.sel.pass);
        EXPECT_FALSE(o.sel.stop);
        EXPECT_EQ(gp::Relation::BOUNDS, o.sel.selected.relation);
        EXPECT_EQ(o.sel.selected.lower + 3, o.sel.selected.upper);
    }
}

// More selected than max_contexts: all_or_count withholds selected_above_threshold, the
// counts kept; partial lists the first max_contexts. Mutations tried: the list's condition
// selected < max_contexts (one fewer listed); the rows' entries never released: fail.
TEST(PatternSelection, SelectedAboveThreshold) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body("{\"none\": [\"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    ASSERT_GT(truth, 3u);
    for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
        for (const Projection proj : { Projection::NONE, Projection::PREDICATE_ONLY }) {
            Ask r;
            r.predicate = compact(pred);
            r.mode = mode;
            r.max_contexts = 3;
            r.projection = proj;
            Outcome o = run(idx, r);
            check_invariants(oracle, pred, o, true, truth, r);
            EXPECT_EQ(gp::Relation::EXACT, o.sel.selected.relation);
            EXPECT_EQ(truth, o.sel.selected.value);
            if (mode == gp::Mode::ALL_OR_COUNT) {
                EXPECT_EQ("selected_above_threshold", o.sel.withheld.value_or(""));
                EXPECT_TRUE(o.sel.chosen.empty());
                EXPECT_EQ(o.bound_bytes, o.held) << "nothing of the list is held";
            } else {
                EXPECT_EQ(3u, o.sel.chosen.size());
                EXPECT_TRUE(o.sel.list_cut);
                EXPECT_FALSE(o.sel.cut);
            }
        }
    }
}

// stop_at_threshold on the selected count: an early exit, stated. Mutation tried: the check
// removed: fails.
TEST(PatternSelection, StopAtThresholdSelected) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body("{\"none\": [\"c2\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    Ask r0;
    r0.predicate = compact(pred);
    const uint64_t all_rows = run(idx, r0).sel.rows;
    for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
        Ask r = r0;
        r.mode = mode;
        r.max_contexts = 2;
        r.stop_at_threshold = true;
        Outcome o = run(idx, r);
        check_invariants(oracle, pred, o, true, truth, r);
        ASSERT_TRUE(o.sel.stop);
        EXPECT_EQ("selection", o.sel.stop->first);
        EXPECT_EQ("max_contexts", o.sel.stop->second);
        EXPECT_EQ(SelectionPass::STOPPED, o.sel.pass);
        EXPECT_LT(o.sel.rows, all_rows) << "an early exit";
        if (mode == gp::Mode::ALL_OR_COUNT) {
            EXPECT_EQ("threshold_crossed", o.sel.withheld.value_or(""));
        } else {
            EXPECT_EQ("max_contexts", o.sel.cut.value_or(""));
            EXPECT_EQ(2u, o.sel.chosen.size());
        }
    }
}

// Mode count reads the annotation with a predicate: the counts, no list (mutation tried: the
// rows' entries never released: fails)
TEST(PatternSelection, CountModeReads) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    Ask r;
    r.mode = gp::Mode::COUNT;
    r.predicate = "{\"at_least\": {\"n\": 2, \"labels\": [\"c1\", \"c4\", \"c5\"]}}";
    r.projection = Projection::PREDICATE_ONLY;
    const Json::Value pred = parse_pattern_body(r.predicate);
    Outcome o = run(idx, r);
    check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true), r);
    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
    EXPECT_GT(o.sel.rows, 0u);
    EXPECT_TRUE(o.sel.chosen.empty());
    EXPECT_FALSE(o.labels);
    EXPECT_EQ(o.bound_bytes, o.held);
}

// The scopes: suffix contexts only (mutation tried: the mirror's labels dropped: fails)
TEST(PatternSelection, ScopesSuffixAndAnyOffset) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    for (const gp::Scope scope : { gp::Scope::SUFFIX, gp::Scope::ANY_OFFSET }) {
        Ask r;
        r.scope = scope;
        r.predicate = "{\"or\": [{\"any\": [\"c3\"]}, {\"none\": [\"c1\"]}]}";
        const Json::Value pred = parse_pattern_body(r.predicate);
        Outcome o = run(idx, r);
        EXPECT_EQ(oracle.contexts("AC", scope), tested_set(o));
        check_invariants(oracle, pred, o, true, oracle_truth(oracle, pred, "AC", true, scope),
                         r);
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
        r.kind = gp::PatternKind::PROTEIN;
        r.either = either;
        r.predicate = "{\"none\": [\"p2\"]}";
        Outcome o = run(idx, r);
        ASSERT_GT(o.tested.size(), 0u);
        for (size_t i = 0; i < o.tested.size(); ++i) {
            std::string inst = o.kmers[i].substr(o.tested[i].offset, 6);
            if (o.tested[i].orientation == gp::Orientation::REVERSE)
                inst = rc(inst);
            EXPECT_EQ("ATG", inst.substr(0, 3));
            EXPECT_TRUE(inst.substr(3) == "AAA" || inst.substr(3) == "AAG") << inst;
        }
        check_invariants(oracle, parse_pattern_body(r.predicate), o, either, std::nullopt, r);
        EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
    }
}

/**
 * The memory account from the smallest up (every 1/150 of what the full run peaked at): each
 * stop stated where it falls (the binding, the descriptors, the rows, the list), never past
 * the account, the decisions made right, complete once the account holds the full run's peak.
 * Mutation tried: a refused row's statement not charged (the reserve kept but not held): the
 * full run's peak drops and the answers at that peak change; and the descriptors' allowance
 * ignored in partial: an answer past the account.
 */
TEST(PatternSelection, MemoryAccount) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body("{\"or\": [{\"any\": [\"c1\"]}, "
                                                "{\"none\": [\"c2\", \"c5\"]}]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
        for (const Projection proj : { Projection::NONE, Projection::PREDICATE_ONLY }) {
            Ask r;
            r.predicate = compact(pred);
            r.mode = mode;
            r.projection = proj;
            r.max_memory_bytes = uint64_t(1) << 30;
            const Outcome full = run(idx, r);
            ASSERT_EQ(SelectionPass::COMPLETED, full.sel.pass);
            bool some_stopped = false;
            for (uint64_t step = 1; step <= 160; ++step) {
                r.max_memory_bytes = std::max<uint64_t>(1, full.peak * step / 150);
                const Outcome o = run(idx, r);
                SCOPED_TRACE(std::to_string(r.max_memory_bytes) + " bytes");
                EXPECT_LE(o.peak, r.max_memory_bytes);
                check_invariants(oracle, pred, o, true, truth, r);
                if (o.bind_stop) {
                    EXPECT_EQ(std::string("max_memory"), o.bind_stop);
                    EXPECT_EQ(SelectionPass::NOT_STARTED, o.sel.pass);
                }
                if (o.sel.pass != SelectionPass::COMPLETED) {
                    some_stopped = true;
                    const bool stated = (o.sel.stop && o.sel.stop->second == "max_memory")
                                            || o.sel.rows_refused.size();
                    EXPECT_TRUE(stated);
                    if (mode == gp::Mode::ALL_OR_COUNT) {
                        EXPECT_EQ("predicate_budget", o.sel.withheld.value_or(""));
                    } else {
                        EXPECT_EQ("max_memory", o.sel.cut.value_or(""));
                    }
                }
                for (const Json::Value &x : o.sel.rows_refused) {
                    EXPECT_EQ("selection", x["phase"].asString());
                    EXPECT_EQ("max_memory", x["reason"].asString());
                    EXPECT_EQ(kK, x["kmer"].asString().size());
                }
            }
            EXPECT_TRUE(some_stopped);
            // the account holds the run's peak and the reads' own decoding (a read's demand
            // is admitted against what is left, beside what the account holds)
            r.max_memory_bytes = 4 * full.peak;
            EXPECT_EQ(SelectionPass::COMPLETED, run(idx, r).sel.pass);
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
    r.mode = gp::Mode::PARTIAL;
    const Outcome full = run(idx, r);
    ASSERT_EQ(k, full.tested.size());
    ASSERT_EQ(SelectionPass::COMPLETED, full.sel.pass);
    // partial's allowance is half of what the bound predicate left: 30 descriptors of 64
    r.max_memory_bytes = full.bound_bytes + 2 * (64 * 30 + 32);
    const Outcome o = run(idx, r);
    EXPECT_EQ(k, o.released);
    EXPECT_EQ(30u, o.tested.size());
    EXPECT_EQ(SelectionPass::STOPPED, o.sel.pass);
    EXPECT_EQ(30u, o.sel.tested.value);
    EXPECT_EQ(gp::Relation::BOUNDS, o.sel.selected.relation);
    EXPECT_EQ(30u, o.sel.selected.lower);
    EXPECT_EQ(31u, o.sel.selected.upper);
    ASSERT_TRUE(o.sel.stop);
    EXPECT_EQ("selection", o.sel.stop->first);
    EXPECT_EQ("max_memory", o.sel.stop->second);
    EXPECT_EQ("max_memory", o.sel.cut.value_or(""));
    EXPECT_EQ(30u, o.sel.chosen.size());
    EXPECT_EQ(1u, o.sel.rows);
    check_invariants(oracle, parse_pattern_body(r.predicate), o, false, 31, r);
    // all_or_count: nothing can be published, nothing read
    r.mode = gp::Mode::ALL_OR_COUNT;
    r.max_memory_bytes = full.bound_bytes + 64 * 30 + 32;
    const Outcome w = run(idx, r);
    EXPECT_EQ(30u, w.tested.size());
    EXPECT_EQ(0u, w.sel.rows);
    EXPECT_EQ("predicate_budget", w.sel.withheld.value_or(""));
    EXPECT_EQ(SelectionPass::STOPPED, w.sel.pass);
}

// A row the account cannot hold: stated (phase selection); all_or_count ends the pass there,
// partial goes on, the contexts needing it untested. Mutation tried: the statement's reserve
// not kept (a refused row unstated): fails.
TEST(PatternSelection, RefusedRowsStated) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kR1, true);
    Oracle oracle(idx);
    const Json::Value pred = parse_pattern_body("{\"any\": [\"c1\", \"c3\"]}");
    const uint64_t truth = oracle_truth(oracle, pred, "AC", true);
    for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
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
        check_invariants(oracle, pred, o, true, truth, r);
        ASSERT_EQ(1u, o.sel.rows_refused.size());
        EXPECT_EQ("selection", o.sel.rows_refused[0]["phase"].asString());
        EXPECT_EQ(SelectionPass::STOPPED, o.sel.pass);
        const std::vector<std::string> order = row_order(oracle, o, true);
        EXPECT_EQ(order.at(1), o.sel.rows_refused[0]["kmer"].asString());
        if (mode == gp::Mode::ALL_OR_COUNT) {
            EXPECT_EQ("predicate_budget", o.sel.withheld.value_or(""));
            EXPECT_EQ(1u, o.sel.rows);
        } else {
            EXPECT_EQ("max_memory", o.sel.cut.value_or(""));
            EXPECT_LE(o.sel.rows, order.size() - 1);
            // the contexts needing the refused row are the untested ones, and the rows of the
            // others were all read
            std::set<std::string> needed;
            for (size_t i = 0; i < o.tested.size(); ++i) {
                const bool needs = o.kmers[i] == order[1] || rc(o.kmers[i]) == order[1];
                EXPECT_EQ(!needs, o.tested[i].decided) << o.kmers[i];
                if (!needs) {
                    needed.insert(o.kmers[i]);
                    if (oracle.present(rc(o.kmers[i])))
                        needed.insert(rc(o.kmers[i]));
                }
            }
            EXPECT_GE(o.sel.rows, needed.size());
        }
    }
}

/**
 * A virtual clock past the work time at each reading in turn (the engine's, the lookups', each
 * read and its paced pieces, each 64 decisions): a time stop stated, time_limited, never an
 * exception; at the reading after the last, the pass completes. Mutations tried: an
 * interrupted read's stop unstated; the decisions' time stop unstated: fail. (The clock before
 * a read removed is equivalent here, the paced read seeing the deadline itself: the abort test
 * catches it.)
 */
TEST(PatternSelection, DeadlineAtEveryReading) {
    // R1 with AC; and k = 33 with one poly-A k-mer carrying 66 contexts of D (A, G or T: its
    // reverse complement H matches A too), decided after one read: the decisions' own clock
    // (every 64) is read there
    struct Case {
        Index idx;
        std::string pattern;
        gp::PatternKind kind;
        std::string predicate;
    };
    std::vector<Case> cases;
    cases.push_back({ build<annot::RowDiffColumnAnnotator>(kK, kR1, true), "AC",
                      gp::PatternKind::DNA, "{\"none\": [\"c4\"]}" });
    cases.push_back({ build<annot::RowDiffColumnAnnotator>(
                              33, { { "c1", "c1_r0", std::string(50, 'A') },
                                    { "c2", "c2_r0", std::string(50, 'C') } }, false),
                      "D", gp::PatternKind::IUPAC, "{\"any\": [\"c1\"]}" });
    for (const Case &c : cases) {
        Oracle oracle(c.idx);
        const Json::Value pred = parse_pattern_body(c.predicate);
        for (const gp::Mode mode : { gp::Mode::ALL_OR_COUNT, gp::Mode::PARTIAL }) {
            auto readings = std::make_shared<uint64_t>(0);
            auto at = std::make_shared<uint64_t>(UINT64_MAX);
            Ask r;
            r.patterns = { c.pattern };
            r.kind = c.kind;
            r.records = c.idx.cth != nullptr;
            r.predicate = c.predicate;
            r.mode = mode;
            r.clock = [readings, at]() {
                ++*readings;
                return Clock::time_point() + (*readings >= *at ? std::chrono::hours(1)
                                                               : std::chrono::hours(0));
            };
            const Outcome full = run(c.idx, r);
            ASSERT_EQ(SelectionPass::COMPLETED, full.sel.pass);
            const uint64_t total = *readings;
            ASSERT_GT(total, 3u);
            bool pass_stopped = false;
            for (uint64_t n = 1; n <= total + 1; ++n) {
                *readings = 0;
                *at = n;
                Outcome o;
                ASSERT_NO_THROW(o = run(c.idx, r));
                SCOPED_TRACE(c.pattern + " reading " + std::to_string(n));
                check_invariants(oracle, pred, o, true, std::nullopt, r);
                EXPECT_TRUE(contains(o.sel.selected, full.sel.selected.value));
                const bool engine_time = o.result->stop
                                            && o.result->stop->reason == gp::StopReason::TIME;
                if (n == total + 1) {
                    EXPECT_EQ(SelectionPass::COMPLETED, o.sel.pass);
                    continue;
                }
                // an incomplete pass says why: a time stop, the engine's or its own
                if (o.sel.pass != SelectionPass::COMPLETED) {
                    EXPECT_TRUE(engine_time || o.sel.time_limited);
                }
                if (o.sel.stop && o.sel.stop->second == "time") {
                    pass_stopped = true;
                    EXPECT_TRUE(o.sel.time_limited);
                    EXPECT_EQ(mode == gp::Mode::ALL_OR_COUNT ? "deadline" : "",
                              o.sel.withheld.value_or(""));
                    if (mode == gp::Mode::PARTIAL) {
                        EXPECT_EQ("time", o.sel.cut.value_or(""));
                    }
                }
            }
            EXPECT_TRUE(pass_stopped);
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
    Ask r;
    r.patterns = { "AC", "GA", "TT" };
    r.predicate = "{\"any\": [\"c1\"]}";
    r.max_predicate_work = 40;
    std::vector<Outcome> all = run_all(idx, r);
    ASSERT_EQ(3u, all.size());
    ASSERT_TRUE(all[0].sel.stop);
    EXPECT_EQ("max_predicate_work", all[0].sel.stop->second);
    for (size_t p = 1; p < 3; ++p) {
        EXPECT_EQ(SelectionPass::NOT_STARTED, all[p].sel.pass);
        ASSERT_TRUE(all[p].sel.stop);
        EXPECT_EQ("max_predicate_work", all[p].sel.stop->second);
        EXPECT_EQ("predicate_budget", all[p].sel.withheld.value_or(""));
        EXPECT_EQ(0u, all[p].sel.rows);
        EXPECT_EQ(gp::Relation::EXACT, all[p].result->contexts->total.relation);
        EXPECT_EQ(oracle.contexts(r.patterns[p]).size(), all[p].result->contexts->total.value);
    }
}

// An unbudgeted annotation with direct access: single cells for at most 16 known labels
// ("columns"), rows above. Mutation tried: always rows: fails.
TEST(PatternSelection, UnbudgetedAnnotation) {
    Index small = build<annot::ColumnCompressed<>>(kK, kR1, true);
    Oracle so(small);
    Ask r;
    r.allow_unbudgeted = true;
    r.predicate = "{\"or\": [{\"any\": [\"c1\"]}, {\"none\": [\"c2\", \"c5\"]}]}";
    Outcome o = run(small, r);
    EXPECT_EQ("columns", o.access);
    check_invariants(so, parse_pattern_body(r.predicate), o, true,
                     oracle_truth(so, parse_pattern_body(r.predicate), "AC", true), r);
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
    check_invariants(wo, pred, o, true, oracle_truth(wo, pred, "AC", true), r);
}

// A graph without its dummy-edge mask (owner decision #16): the admitted release enumerates
// the candidates and drops the dummies; the selection is the masked twin's and the oracle's
// (mutation tried: the mirror's labels dropped: fails)
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
        EXPECT_EQ(gp::Relation::EXACT, b.result->contexts->total.relation) << p;
        EXPECT_EQ(tested_set(a), tested_set(b)) << p;
        ASSERT_EQ(a.tested.size(), b.tested.size());
        for (size_t i = 0; i < a.tested.size(); ++i) {
            EXPECT_EQ(a.tested[i].selected, b.tested[i].selected) << a.kmers[i];
        }
        EXPECT_EQ(a.sel.selected.value, b.sel.selected.value);
        check_invariants(oracle, parse_pattern_body(r.predicate), b, true, std::nullopt, r);
    }
}

} // namespace
