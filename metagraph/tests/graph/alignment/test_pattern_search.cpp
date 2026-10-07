/**
 * The oracle suite of the pattern search (docs/DESIGN-pattern-search.md §13, increments 0-2):
 * tiny graphs built from explicit records, checked against two brute-force oracles that
 * share no code with the engine:
 *  - the graph-walk oracle enumerates the built graph's k-mers (every valid edge of the
 *    DBGSuccinct; on a wrapped PRIMARY graph also the wrapper's reverse complements, numbered
 *    by CanonicalDBG::reverse_complement) and lists every (node, offset, orientation) whose
 *    k-mer instantiates the oriented pattern at that offset: the expected graph contexts,
 *    with their node ids and in answer order;
 *  - the record-scan oracle scans the records' retained islands (the stretches over A, C, G,
 *    T) as strings on both strands and lists the distinct (k-mer, offset, orientation) they
 *    contain: what the records say the contexts must be, independently of the graph.
 * count() must equal the first in every unit and relation, enumerate() must return exactly
 * its list in its order, and the first must equal the second wherever no k-mer was pruned.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/aligner_seeder_methods.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

#if ! _PROTEIN_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;
using mtg::test::build_graph_batch;

typedef DeBruijnGraph::node_index node_index;

constexpr uint64_t kManySteps = 1'000'000'000;


// ---------------------------------------------------------------- helpers

std::string rev_comp(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

BaseSet base_bit(char c) {
    switch (c) {
        case 'A': return kBaseA;
        case 'C': return kBaseC;
        case 'G': return kBaseG;
        case 'T': return kBaseT;
        default: return 0;  // N, $: never matched by a pattern position (§3)
    }
}

bool matches(const std::vector<BaseSet> &q, std::string_view s) {
    if (s.size() != q.size())
        return false;
    for (size_t i = 0; i < q.size(); ++i) {
        if (!(q[i] & base_bit(s[i])))
            return false;
    }
    return true;
}

std::shared_ptr<DeBruijnGraph> build(size_t k, const std::vector<std::string> &records,
                                     DeBruijnGraph::Mode mode, bool batch = false) {
    return batch ? build_graph_batch<DBGSuccinct>(k, records, mode)
                 : build_graph<DBGSuccinct>(k, records, mode);
}

const DBGSuccinct& base_dbg(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// a request with no information floor (tiny k) and the given scope and strands
Request make_request(Scope scope = Scope::ANY_OFFSET, Strands strands = Strands::BOTH) {
    Request request;
    request.scope = scope;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = 1'000'000'000;
    request.max_anchors = 1'000'000'000;
    return request;
}

Budget unbounded_budget(uint64_t max_steps = kManySteps) {
    return Budget(max_steps, Deadline::unbounded());
}

struct Ctx {
    node_index node;
    uint32_t offset;
    Orientation orientation;

    bool operator<(const Ctx &o) const {
        return std::tie(node, offset, orientation) < std::tie(o.node, o.offset, o.orientation);
    }
    bool operator==(const Ctx &o) const {
        return std::tie(node, offset, orientation) == std::tie(o.node, o.offset, o.orientation);
    }
};

std::ostream& operator<<(std::ostream &out, const Ctx &c) {
    return out << "(" << c.node << ", " << c.offset << ", "
               << orientation_key(c.orientation) << ")";
}

// the orientations the engine searches, with the oriented pattern of each (§3)
std::vector<std::pair<Orientation, Pattern>> orientations(const Pattern &pattern,
                                                          Strands strands) {
    if (pattern.is_palindromic())
        return { { Orientation::PALINDROMIC, pattern } };
    std::vector<std::pair<Orientation, Pattern>> result;
    if (strands != Strands::REVERSE)
        result.emplace_back(Orientation::FORWARD, pattern);
    if (strands != Strands::FORWARD)
        result.emplace_back(Orientation::REVERSE, pattern.reverse_complement());
    return result;
}

std::vector<uint32_t> scope_offsets(size_t L, size_t k, Scope scope) {
    if (L > k)
        return { 0 };
    if (scope == Scope::SUFFIX)
        return { static_cast<uint32_t>(k - L) };
    std::vector<uint32_t> offsets;
    for (uint32_t p = 0; p + L <= k; ++p) {
        offsets.push_back(p);
    }
    return offsets;
}

std::vector<BaseSet> window(const Pattern &q, size_t k) {
    return std::vector<BaseSet>(q.positions().begin(),
                                q.positions().begin() + std::min(q.length(), k));
}


// ---------------------------------------------------------------- the oracles

// every k-mer of the served graph with its node id, directly from the graph
std::vector<std::pair<node_index, std::string>> graph_kmers(const DeBruijnGraph &graph) {
    std::vector<std::pair<node_index, std::string>> kmers;
    const DBGSuccinct &dbg_succ = base_dbg(graph);
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (!dbg_succ.in_graph(y))
            continue;
        kmers.emplace_back(y, graph.get_node_sequence(y));
        if (canonical) {
            node_index z = canonical->reverse_complement(y);
            if (z != y)
                kmers.emplace_back(z, graph.get_node_sequence(z));
        }
    }
    return kmers;
}

// the graph-walk oracle: every context of the served graph, in answer order
std::vector<Ctx> walk_oracle(const DeBruijnGraph &graph, const Pattern &pattern,
                             const Request &request) {
    const size_t k = graph.get_k();
    std::vector<Ctx> result;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
            std::vector<BaseSet> w = window(q, k);
            for (uint32_t p : scope_offsets(pattern.length(), k, request.scope)) {
                if (matches(w, std::string_view(kmer).substr(p, w.size())))
                    result.push_back(Ctx { node, p, orientation });
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

typedef std::set<std::tuple<std::string, uint32_t, Orientation>> SpelledContexts;

// the record-scan oracle: the k-mers of the records' islands (both orientations unless
// BASIC) and the contexts they spell
SpelledContexts record_oracle(const std::vector<std::string> &records, size_t k,
                              DeBruijnGraph::Mode mode, const Pattern &pattern,
                              const Request &request) {
    std::set<std::string> kmers;
    for (const std::string &record : records) {
        for (size_t i = 0; i + k <= record.size(); ++i) {
            std::string kmer = record.substr(i, k);
            if (!std::all_of(kmer.begin(), kmer.end(), [](char c) { return base_bit(c); }))
                continue;
            kmers.insert(kmer);
            if (mode != DeBruijnGraph::BASIC)
                kmers.insert(rev_comp(kmer));
        }
    }
    SpelledContexts result;
    for (const std::string &kmer : kmers) {
        for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
            std::vector<BaseSet> w = window(q, k);
            for (uint32_t p : scope_offsets(pattern.length(), k, request.scope)) {
                if (matches(w, std::string_view(kmer).substr(p, w.size())))
                    result.emplace(kmer, p, orientation);
            }
        }
    }
    return result;
}

SpelledContexts spelled(const DeBruijnGraph &graph, const std::vector<Ctx> &contexts) {
    SpelledContexts result;
    for (const Ctx &c : contexts) {
        result.emplace(graph.get_node_sequence(c.node), c.offset, c.orientation);
    }
    return result;
}


// ---------------------------------------------------------------- the checks

std::vector<Ctx> run_enumerate(const PatternSearch &engine, const Pattern &pattern,
                               const Request &request, Result *result,
                               uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    std::vector<Ctx> contexts;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        contexts.push_back(Ctx { c.node, c.offset, c.orientation });
    });
    return contexts;
}

/**
 * Every check of one (graph, pattern, request) against the oracles: counts in every unit
 * and relation, the released list and its order, the spelled k-mers, the PARTIAL prefix and
 * the ALL_OR_COUNT threshold. |records| enables the record-scan comparison.
 */
void check_against_oracles(const DeBruijnGraph &graph, const Pattern &pattern,
                           const Request &request,
                           const std::vector<std::string> *records = nullptr,
                           DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC,
                           size_t *num_expected = nullptr) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    PatternSearch engine(graph);

    std::vector<Ctx> expected = walk_oracle(graph, pattern, request);
    if (num_expected)
        *num_expected = expected.size();
    if (records && L <= k) {
        // the graph holds exactly the records' k-mers (no cleaning in these builds)
        ASSERT_EQ(record_oracle(*records, k, mode, pattern, request), spelled(graph, expected))
            << pattern.text();
    }

    Budget budget = unbounded_budget();
    Result counted = engine.count(pattern, request, budget);
    ASSERT_FALSE(counted.refusal) << counted.refusal->message;
    EXPECT_FALSE(counted.stop);
    EXPECT_FALSE(counted.time_limited);
    EXPECT_FALSE(counted.extraction);

    std::map<uint32_t, uint64_t> by_offset;
    std::map<Orientation, uint64_t> by_orientation;
    for (const Ctx &c : expected) {
        ++by_offset[c.offset];
        ++by_orientation[c.orientation];
    }

    const Count &total = L > k ? counted.anchors->total : counted.contexts->total;
    EXPECT_EQ(Relation::EXACT, total.relation) << pattern.text();
    EXPECT_EQ(expected.size(), total.value) << pattern.text();
    EXPECT_EQ(L > k ? Unit::ANCHORS : Unit::GRAPH_CONTEXTS, total.unit);

    const auto &orientation_counts = L > k ? counted.anchors->by_orientation
                                           : counted.contexts->by_orientation;
    ASSERT_EQ(orientations(pattern, request.strands).size(), orientation_counts.size());
    for (const auto &[orientation, count] : orientation_counts) {
        EXPECT_EQ(Relation::EXACT, count.relation);
        EXPECT_EQ(by_orientation[orientation], count.value)
            << pattern.text() << " " << orientation_key(orientation);
    }
    if (L <= k) {
        std::vector<uint32_t> offsets = scope_offsets(L, k, request.scope);
        ASSERT_EQ(offsets.size(), counted.contexts->by_offset.size());
        for (uint32_t p : offsets) {
            const Count &count = counted.contexts->by_offset.at(p);
            EXPECT_EQ(Relation::EXACT, count.relation);
            EXPECT_EQ(by_offset[p], count.value) << pattern.text() << " offset " << p;
        }
        EXPECT_EQ(by_offset[k - L], counted.contexts->suffix.value);
    } else {
        // paths are a later increment: unknown, unless no anchor exists
        EXPECT_EQ(expected.empty() ? Relation::EXACT : Relation::UNKNOWN,
                  counted.anchors->paths.relation);
    }

    // enumerate(): discovery identical to count(), then every context in answer order
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result enumerated;
    std::vector<Ctx> released = run_enumerate(engine, pattern, all, &enumerated);
    EXPECT_EQ(counted.work.steps, enumerated.work.steps);
    EXPECT_EQ(counted.work.ranges_visited, enumerated.work.ranges_visited);
    ASSERT_TRUE(enumerated.extraction);
    if (L > k && !request.release_anchors) {
        // the results of a long pattern are its paths (a later increment): withheld
        EXPECT_TRUE(released.empty());
        if (expected.empty()) {
            EXPECT_TRUE(enumerated.extraction->complete);
        } else {
            EXPECT_EQ(Withheld::PATHS_LATER_INCREMENT, enumerated.extraction->withheld);
        }
        // on request, the anchors are released like contexts (offset 0), and checked so
        Request anchors = request;
        anchors.release_anchors = true;
        check_against_oracles(graph, pattern, anchors);
        return;
    }
    EXPECT_TRUE(enumerated.extraction->complete);
    EXPECT_FALSE(enumerated.extraction->withheld);
    EXPECT_EQ(expected.size(), enumerated.extraction->returned);
    ASSERT_EQ(expected, released) << pattern.text();

    for (const Ctx &c : released) {
        // every released context spells the oriented pattern at its offset
        std::string kmer = graph.get_node_sequence(c.node);
        Pattern q = c.orientation == Orientation::REVERSE ? pattern.reverse_complement()
                                                          : pattern;
        EXPECT_TRUE(matches(window(q, k), std::string_view(kmer).substr(c.offset, L)));
    }

    // ALL_OR_COUNT: all or nothing at the threshold (max_anchors for released anchors)
    auto set_cap = [&](Request *r, uint64_t cap) {
        (L > k ? r->max_anchors : r->max_contexts) = cap;
    };
    if (!expected.empty()) {
        Request tight = all;
        set_cap(&tight, expected.size() - 1);
        Result withheld;
        EXPECT_TRUE(run_enumerate(engine, pattern, tight, &withheld).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
        EXPECT_FALSE(withheld.extraction->complete);
        EXPECT_EQ(Relation::EXACT, L > k ? withheld.anchors->total.relation
                                         : withheld.contexts->total.relation);
    }

    // PARTIAL: the first max_contexts in answer order, the cut stated
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()) }) {
        set_cap(&partial, cap);
        Result cut;
        std::vector<Ctx> prefix = run_enumerate(engine, pattern, partial, &cut);
        ASSERT_EQ(std::min<uint64_t>(cap, expected.size()), prefix.size());
        EXPECT_TRUE(std::equal(prefix.begin(), prefix.end(), expected.begin()));
        if (cap >= expected.size()) {
            EXPECT_TRUE(cut.extraction->complete);
            EXPECT_FALSE(cut.extraction->cut);
        } else {
            EXPECT_FALSE(cut.extraction->complete);
            EXPECT_EQ(L > k ? StopReason::MAX_ANCHORS : StopReason::MAX_CONTEXTS,
                      cut.extraction->cut);
        }
    }
}

// the contexts the engine releases, for the named cases
std::vector<Ctx> contexts_of(const DeBruijnGraph &graph, const Pattern &pattern,
                             const Request &request) {
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result result;
    return run_enumerate(PatternSearch(graph), pattern, all, &result);
}

Result count_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    return PatternSearch(graph).count(pattern, request, budget);
}


// ---------------------------------------------------------------- patterns and counts

TEST(PatternSearch, ParsePatterns) {
    Pattern dna = Pattern::parse(PatternKind::DNA, "acgT");
    EXPECT_EQ("ACGT", dna.text());
    EXPECT_EQ(4u, dna.length());
    EXPECT_TRUE(dna.is_exact());
    EXPECT_TRUE(dna.is_palindromic());
    EXPECT_DOUBLE_EQ(8.0, dna.information_bits());

    for (std::string bad : { "", "ACGU", "AC-G", "AC G", "ACGN", "ACGR" }) {
        try {
            Pattern::parse(PatternKind::DNA, bad);
            FAIL() << "accepted " << bad;
        } catch (const PatternError &e) {
            EXPECT_EQ("bad_alphabet", e.code());
        }
    }
    EXPECT_THROW(Pattern::parse(PatternKind::IUPAC, "ACGU"), PatternError);
    EXPECT_THROW(Pattern::parse(PatternKind::IUPAC, "AC.G"), PatternError);

    Pattern iupac = Pattern::parse(PatternKind::IUPAC, "ACGTRYSWKMBDHVN");
    EXPECT_FALSE(iupac.is_exact());
    EXPECT_EQ("NBDHVKMWSRYACGT", iupac.reverse_complement().text());
    EXPECT_EQ(iupac.text(), iupac.reverse_complement().reverse_complement().text());
    EXPECT_DOUBLE_EQ(0.0, Pattern::parse(PatternKind::IUPAC, "NNN").information_bits());
    EXPECT_DOUBLE_EQ(1.0, Pattern::parse(PatternKind::IUPAC, "R").information_bits());
    EXPECT_NEAR(std::log2(4.0 / 3), Pattern::parse(PatternKind::IUPAC, "B").information_bits(),
                1e-12);
    EXPECT_DOUBLE_EQ(3.0, iupac.information_bits(0, 2) - 1.0);

    EXPECT_TRUE(Pattern::parse(PatternKind::IUPAC, "RY").is_palindromic());
    EXPECT_TRUE(Pattern::parse(PatternKind::IUPAC, "NN").is_palindromic());
    EXPECT_FALSE(Pattern::parse(PatternKind::IUPAC, "NA").is_palindromic());
    EXPECT_EQ("TN", Pattern::parse(PatternKind::IUPAC, "NA").reverse_complement().text());
    EXPECT_TRUE(Pattern::parse(PatternKind::DNA, "A").reverse_complement().is_exact());
    EXPECT_EQ(PatternKind::DNA, dna.reverse_complement().kind());
    EXPECT_EQ(kBaseC | kBaseT, Pattern::parse(PatternKind::IUPAC, "AY").allowed(1, "A"));
}

TEST(PatternSearch, CountAlgebra) {
    const Unit u = Unit::GRAPH_CONTEXTS;
    Count c = Count::exact(u, 3);
    c += Count::bounds(u, 1, 4);
    EXPECT_EQ(Relation::BOUNDS, c.relation);
    EXPECT_EQ(4u, c.lower);
    EXPECT_EQ(7u, c.upper);

    Count unknown = Count::unknown(u);
    unknown += Count::unknown(u);
    EXPECT_EQ(Relation::UNKNOWN, unknown.relation);
    unknown += Count::exact(u, 5);
    EXPECT_EQ(Relation::AT_LEAST, unknown.relation);
    EXPECT_EQ(5u, unknown.value);

    Count same = Count::bounds(u, 0, 0);
    same += Count::exact(u, 0);
    EXPECT_EQ(Relation::BOUNDS, same.relation);  // never promoted

    Count at_least = Count::at_least(u, 2);
    at_least += Count::bounds(u, 1, 9);
    EXPECT_EQ(Relation::AT_LEAST, at_least.relation);
    EXPECT_EQ(3u, at_least.value);

    Count exact = Count::exact(u, 2);
    exact += Count::exact(u, 3);
    EXPECT_EQ(Relation::EXACT, exact.relation);
    EXPECT_EQ(5u, exact.value);
}

TEST(PatternSearch, BudgetAndDeadline) {
    Budget budget(10, Deadline::unbounded());
    EXPECT_TRUE(budget.charge(4));
    EXPECT_TRUE(budget.charge(6));
    EXPECT_FALSE(budget.charge(1));
    EXPECT_EQ(StopReason::MAX_STEPS, budget.stopped());
    EXPECT_EQ(10u, budget.steps_used());
    EXPECT_FALSE(budget.charge(0));  // sticky

    // an injected clock: 1 ms per reading; the clock is read at stride crossings only
    auto start = Deadline::Clock::now();
    auto now = std::make_shared<Deadline::Clock::time_point>(start);
    auto clock = [now]() { return *now += std::chrono::milliseconds(1); };
    Budget timed(1'000'000, Deadline(start, 10, 5, clock));
    EXPECT_TRUE(timed.check_time());  // 1 ms
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(timed.charge(Budget::kClockStride));  // 2, 3, 4 ms
    }
    EXPECT_TRUE(timed.charge(Budget::kClockStride - 1));  // no reading: no crossing
    EXPECT_FALSE(timed.charge(1));  // 5 ms: the work time (10 - 5) has passed
    EXPECT_EQ(StopReason::TIME, timed.stopped());
    EXPECT_FALSE(timed.deadline().respond_expired());  // 6 ms < 10
}


// ---------------------------------------------------------------- graphs

TEST(PatternSearch, GraphSupport) {
    std::vector<std::string> records { "ACGTACGGTTAC" };
    auto basic = build(5, records, DeBruijnGraph::BASIC);
    auto canonical = build(5, records, DeBruijnGraph::CANONICAL);
    auto primary = build(5, records, DeBruijnGraph::PRIMARY);

    GraphSupport s = PatternSearch::support(*basic);
    EXPECT_TRUE(s.supported);
    EXPECT_EQ(GraphMode::BASIC, s.mode);
    EXPECT_TRUE(s.strand_stated);
    EXPECT_TRUE(s.mask_present);
    EXPECT_EQ("$ACGT", s.alphabet);
    EXPECT_EQ(5u, s.k);
    EXPECT_EQ(2u, s.scopes.size());

    s = PatternSearch::support(*canonical);
    EXPECT_TRUE(s.supported);
    EXPECT_EQ(GraphMode::CANONICAL, s.mode);
    EXPECT_FALSE(s.strand_stated);

    s = PatternSearch::support(*primary);
    EXPECT_TRUE(s.supported);
    EXPECT_EQ(GraphMode::PRIMARY, s.mode);
    EXPECT_EQ(std::vector<Scope>{ Scope::ANY_OFFSET }, s.scopes);

    // the PRIMARY graph without its wrapper: refused
    const auto &unwrapped = dynamic_cast<const CanonicalDBG&>(*primary).get_graph();
    s = PatternSearch::support(unwrapped);
    EXPECT_FALSE(s.supported);
    EXPECT_EQ("primary_unwrapped", s.reason);
    EXPECT_THROW(PatternSearch engine(unwrapped), std::invalid_argument);

    // another representation: refused
    auto hash = build_graph<DBGHashFast>(5, records);
    s = PatternSearch::support(*hash);
    EXPECT_FALSE(s.supported);
    EXPECT_EQ("representation_unsupported", s.reason);
}

TEST(PatternSearch, GraphWithoutMaskRefused) {
    // without the mask every dummy edge would count as a k-mer (§4): refused, never guessed
    auto graph = build(4, { "ACGTT", "TTGCA" }, DeBruijnGraph::BASIC);
    auto &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(*graph));
    dbg_succ.reset_mask();
    GraphSupport s = PatternSearch::support(*graph);
    EXPECT_FALSE(s.supported);
    EXPECT_FALSE(s.mask_present);
    EXPECT_EQ("mask_required", s.reason);
    EXPECT_THROW(PatternSearch engine(*graph), std::invalid_argument);
}

TEST(PatternSearch, MainDummySourceIsInvalid) {
    // count_edges_with_last_symbol discounts edges with W = $ from the invalid ones; that is
    // sound only if no such edge is valid, the main dummy source (edge 1) included
    for (bool batch : { false, true }) {
        auto graph = build(4, { "ACGTT", "TTGCA", "GGGGA" }, DeBruijnGraph::BASIC, batch);
        const DBGSuccinct &dbg_succ = base_dbg(*graph);
        const auto &boss = dbg_succ.get_boss();
        for (node_index e = 1; e <= dbg_succ.max_index(); ++e) {
            // braces: the gtest macro expands to an if/else (GCC's -Wdangling-else)
            if (boss.get_W(e) % boss.alph_size == 0) {
                EXPECT_FALSE(dbg_succ.in_graph(e)) << e;
            }
        }
    }
}


// ---------------------------------------------------------------- the range DFS

// the loop of suffix_to_prefix as it was before its symbol set became a parameter (copied
// from aligner_seeder_methods.cpp at 804731aa): the default must reproduce it exactly
template <class BOSSEdgeRange>
void suffix_to_prefix_804731aa(const DBGSuccinct &dbg_succ,
                               const BOSSEdgeRange &index_range,
                               const std::function<void(DBGSuccinct::node_index)> &callback) {
    const auto &boss = dbg_succ.get_boss();
    auto call_nodes_in_range = [&](const BOSSEdgeRange &final_range) {
        const auto &[first, last, seed_length] = final_range;
        for (boss::BOSS::edge_index i = first; i <= last; ++i) {
            DBGSuccinct::node_index node = dbg_succ.validate_edge(i);
            if (node)
                callback(node);
        }
    };
    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }
    std::vector<BOSSEdgeRange> range_stack { index_range };
    while (range_stack.size()) {
        BOSSEdgeRange cur_range = std::move(range_stack.back());
        range_stack.pop_back();
        ++std::get<2>(cur_range);
        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
            auto next_range = cur_range;
            auto &[first, last, seed_length] = next_range;
            if (boss.tighten_range(&first, &last, s)) {
                if (seed_length == boss.get_k()) {
                    call_nodes_in_range(next_range);
                } else {
                    range_stack.emplace_back(std::move(next_range));
                }
            }
        }
    }
}

TEST(PatternSearch, SuffixToPrefixDefaultUnchanged) {
    // the seeder's calls (no symbol set) must call the same nodes in the same order as
    // before, from every suffix range of length 1 .. k - 1 the seeder can start from
    std::mt19937 rng(42);
    for (size_t k : { 3, 4, 5, 7 }) {
        for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
            for (bool batch : { false, true }) {
                std::vector<std::string> records;
                for (int r = 0; r < 4; ++r) {
                    std::string record;
                    for (int i = 0; i < 25; ++i) {
                        record.push_back("ACGT"[rng() % 4]);
                    }
                    records.push_back(record);
                }
                auto graph = build(k, records, mode, batch);
                const DBGSuccinct &dbg_succ = base_dbg(*graph);
                const auto &boss = dbg_succ.get_boss();

                size_t calls = 0;
                for (size_t len = 1; len < k; ++len) {
                    for (size_t code = 0; code < (1u << (2 * len)); ++code) {
                        std::string s;
                        for (size_t i = 0; i < len; ++i) {
                            s.push_back("ACGT"[(code >> (2 * i)) & 3]);
                        }
                        auto encoded = boss.encode(s);
                        auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
                        size_t seed_length = end - encoded.begin();
                        if (!seed_length)
                            continue;
                        auto range = std::make_tuple(boss.pred_last(first - 1) + 1, last,
                                                     seed_length);
                        std::vector<node_index> before, after;
                        suffix_to_prefix_804731aa(dbg_succ, range,
                                                  [&](node_index n) { before.push_back(n); });
                        align::suffix_to_prefix(dbg_succ, range,
                                                [&](node_index n) { after.push_back(n); });
                        ASSERT_EQ(before, after) << s << " k=" << k;
                        calls += after.size();
                    }
                }
                EXPECT_LT(0u, calls);
            }
        }
    }
}

TEST(PatternSearch, SuffixToPrefixSymbolSet) {
    // a symbol set restricting the extension: only the k-mers whose prefix continues with
    // the chosen symbol are called
    auto graph = build(5, { "ACGTACCTGA", "ACATTTGACA" }, DeBruijnGraph::BASIC);
    const DBGSuccinct &dbg_succ = base_dbg(*graph);
    const auto &boss = dbg_succ.get_boss();
    auto encoded = boss.encode("AC");
    auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
    ASSERT_EQ(2, end - encoded.begin());
    auto range = std::make_tuple(boss.pred_last(first - 1) + 1, last, size_t(2));

    std::vector<std::string> called;
    align::suffix_to_prefix(dbg_succ, range,
        [&](node_index n) { called.push_back(dbg_succ.get_node_sequence(n)); },
        [&](const boss::BOSS &b, const auto &, const auto &try_symbol) {
            try_symbol(b.encode('G'));
        });
    // nodes ACG? → k-mers ACGG? only: "ACGT" continues with A in ACGTA
    std::sort(called.begin(), called.end());
    EXPECT_EQ(std::vector<std::string>{}, called);

    called.clear();
    align::suffix_to_prefix(dbg_succ, range,
        [&](node_index n) { called.push_back(dbg_succ.get_node_sequence(n)); },
        [&](const boss::BOSS &b, const auto &r, const auto &try_symbol) {
            // depth 3: G after AC; depth 4: every symbol
            if (std::get<2>(r) == 3) {
                try_symbol(b.encode('G'));
            } else {
                align::NonSentinelSymbols()(b, r, try_symbol);
            }
        });
    std::sort(called.begin(), called.end());
    EXPECT_EQ(std::vector<std::string>{ "ACGTA" }, called);
}


// ---------------------------------------------------------------- §13: the named cases

TEST(PatternSearch, PrunedDummyRecordStart) {
    // k = 4: ACGA is a k-mer of TACGAT, so the record ACGA starts without a dummy chain of
    // its own; AC at its start is found by any_offset, never through a dummy
    std::vector<std::string> records { "TACGAT", "ACGA" };
    for (bool batch : { false, true }) {
        auto graph = build(4, records, DeBruijnGraph::BASIC, batch);
        Pattern ac = Pattern::parse(PatternKind::DNA, "AC");

        // suffix: no k-mer ends with AC (a $$AC dummy, if any, is invalid)
        Result suffix = count_of(*graph, ac, make_request(Scope::SUFFIX));
        EXPECT_EQ(Relation::EXACT, suffix.contexts->total.relation);
        EXPECT_EQ(0u, suffix.contexts->total.value);

        // any_offset: TACG at 1 and ACGA at 0 (+); GT (-) nowhere
        Result any = count_of(*graph, ac, make_request());
        EXPECT_EQ(2u, any.contexts->total.value);
        EXPECT_EQ(2u, any.contexts->by_orientation.at(Orientation::FORWARD).value);
        EXPECT_EQ(0u, any.contexts->by_orientation.at(Orientation::REVERSE).value);
        EXPECT_EQ(1u, any.contexts->by_offset.at(0).value);
        EXPECT_EQ(1u, any.contexts->by_offset.at(1).value);
        EXPECT_EQ(0u, any.contexts->by_offset.at(2).value);

        check_against_oracles(*graph, ac, make_request(), &records);
        check_against_oracles(*graph, ac, make_request(Scope::SUFFIX), &records);
    }
}

TEST(PatternSearch, LengthKLookup) {
    // L = k: positions 0 .. k-2 on node ranges, the last on W (pick_edge's rule)
    std::vector<std::string> records { "AACT", "GAAC" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    Pattern aac = Pattern::parse(PatternKind::DNA, "AAC");
    for (Scope scope : { Scope::SUFFIX, Scope::ANY_OFFSET }) {
        Result r = count_of(*graph, aac, make_request(scope, Strands::FORWARD));
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
        EXPECT_EQ(1u, r.contexts->total.value);
        EXPECT_EQ(1u, r.contexts->by_offset.size());
        check_against_oracles(*graph, aac, make_request(scope), &records);
    }
}

TEST(PatternSearch, OnePlainTwoMarkedEdges) {
    // ACAT, CCAT, GCAT: three edges with W = T into the target CAT from nodes ending with CA:
    // one plain, two marked. A count of the plain rank alone would say 1
    std::vector<std::string> records { "ACAT", "CCAT", "GCAT" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);
    const DBGSuccinct &dbg_succ = base_dbg(*graph);
    const auto &boss = dbg_succ.get_boss();

    auto encoded = boss.encode("CA");
    auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
    ASSERT_EQ(2, end - encoded.begin());
    first = boss.pred_last(first - 1) + 1;
    auto t = boss.encode('T');
    EXPECT_EQ(1u, boss.rank_W(last, t) - boss.rank_W(first - 1, t));
    EXPECT_EQ(2u, boss.rank_W(last, t + boss.alph_size)
                    - boss.rank_W(first - 1, t + boss.alph_size));

    Pattern cat = Pattern::parse(PatternKind::DNA, "CAT");
    Result r = count_of(*graph, cat, make_request(Scope::SUFFIX, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
    EXPECT_EQ(3u, r.contexts->total.value);
    check_against_oracles(*graph, cat, make_request(Scope::SUFFIX), &records);
    check_against_oracles(*graph, cat, make_request(), &records);
}

TEST(PatternSearch, SinkDummyInsideRange) {
    // k = 3: the record GAC ends at the node AC, whose edge AC$ is a dummy sink inside the
    // range of nodes ending with AC; it is not a k-mer starting with AC
    std::vector<std::string> records { "GAC" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    const DBGSuccinct &dbg_succ = base_dbg(*graph);
    const auto &boss = dbg_succ.get_boss();
    auto encoded = boss.encode("AC");
    auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
    ASSERT_EQ(2, end - encoded.begin());
    first = boss.pred_last(first - 1) + 1;
    EXPECT_LT(dbg_succ.count_valid_edges_in_range(first, last), last - first + 1);

    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    Result r = count_of(*graph, ac, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(1u, r.contexts->total.value);  // GAC at offset 1
    EXPECT_EQ(0u, r.contexts->by_offset.at(0).value);
    check_against_oracles(*graph, ac, make_request(), &records);
}

TEST(PatternSearch, SourceDummyLeavingRange) {
    // k = 4: the record ACGT has the dummy $$AC; the leaf of AC (nodes ending with A) has a
    // candidate with W = C that is no k-mer. Its scan removes it: 0, exactly
    std::vector<std::string> records { "ACGT" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    Result r = count_of(*graph, ac, make_request(Scope::SUFFIX, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
    EXPECT_EQ(0u, r.contexts->total.value);
    EXPECT_EQ(1u, r.work.mask_scans);
    EXPECT_LT(r.work.ranges_visited, r.work.steps);
    check_against_oracles(*graph, ac, make_request(Scope::SUFFIX), &records);
}

TEST(PatternSearch, InterruptedMaskScanBounds) {
    // the leaf of AC (nodes ending with A) holds TTA -> C, GGA -> C and the source dummy
    // $$A with edges C and G (records starting with AC and AG): R = 3 candidates, 2 invalid
    // non-sentinel edges, 2 true contexts (TTAC, GGAC)
    std::vector<std::string> records { "ACGT", "AGTT", "TTAC", "GGAC" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    Request request = make_request(Scope::SUFFIX, Strands::FORWARD);

    Result full = count_of(*graph, ac, request);
    ASSERT_EQ(Relation::EXACT, full.contexts->total.relation);
    ASSERT_EQ(2u, full.contexts->total.value);
    ASSERT_EQ(1u, full.work.mask_scans);
    const uint64_t discovery = full.work.ranges_visited;
    const uint64_t scan = full.work.steps - discovery;
    ASSERT_EQ(2u, scan);

    // every range discovered, no scan: bounds [R - J, R] = [1, 3]
    Result none = count_of(*graph, ac, request, discovery);
    ASSERT_TRUE(none.stop);
    EXPECT_EQ(StopPhase::MASK_SCAN, none.stop->phase);
    EXPECT_EQ(StopReason::MAX_STEPS, none.stop->reason);
    const Count &bounds = none.contexts->total;
    EXPECT_EQ(Relation::BOUNDS, bounds.relation);
    EXPECT_EQ(1u, bounds.lower);
    EXPECT_EQ(3u, bounds.upper);
    EXPECT_EQ(bounds.lower, bounds.value);
    EXPECT_EQ(Relation::BOUNDS, none.contexts->by_offset.at(2).relation);

    // one edge scanned ($$AC, matching): bounds [3 - 1 - 1, 3 - 1] = [1, 2]
    Result half = count_of(*graph, ac, request, discovery + 1);
    EXPECT_EQ(StopPhase::MASK_SCAN, half.stop->phase);
    EXPECT_EQ(Relation::BOUNDS, half.contexts->total.relation);
    EXPECT_LE(half.contexts->total.lower, 2u);
    EXPECT_GE(half.contexts->total.upper, 2u);
    EXPECT_LE(half.contexts->total.upper, 3u);
    EXPECT_GE(half.contexts->total.lower, bounds.lower);

    // a discovery stop: at_least, never bounds (an undiscovered branch has no upper bound)
    Result early = count_of(*graph, ac, request, discovery - 1);
    EXPECT_EQ(StopPhase::DISCOVERY, early.stop->phase);
    EXPECT_EQ(Relation::AT_LEAST, early.contexts->total.relation);
    EXPECT_LE(early.contexts->total.value, 2u);
}

TEST(PatternSearch, DNA5Flank) {
    // §4.1: the flank admits every symbol of the build's alphabet. On a DNA5 build ACNTA is a
    // k-mer at k = 5 and contains AC at offset 0; on a DNA4 build the record splits into the
    // islands AC and TA, both shorter than k, so nothing is indexed and nothing is claimed
    std::vector<std::string> records { "ACNTA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    Result r = count_of(*graph, ac, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
#if _DNA5_GRAPH
    EXPECT_EQ(1u, r.contexts->total.value);
    EXPECT_EQ(1u, r.contexts->by_offset.at(0).value);
#else
    EXPECT_EQ(0u, r.contexts->total.value);
    check_against_oracles(*graph, ac, make_request(), &records);
#endif
}

TEST(PatternSearch, IslandShorterThanK) {
    // k = 5, DNA4: AAAAANACNCCCCC keeps AAAAA and CCCCC; AC sits in a two-base island
    std::vector<std::string> records { "AAAAANACNCCCCC" };
    for (bool batch : { false, true }) {
        auto graph = build(5, records, DeBruijnGraph::BASIC, batch);
        Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
        Result r = count_of(*graph, ac, make_request());
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
        EXPECT_EQ(0u, r.contexts->total.value);
        check_against_oracles(*graph, ac, make_request(), &records);
    }
}

TEST(PatternSearch, ZeroSuffixAnchorsPrefixContext) {
    // k = 5, ACGTT: AC ends no k-mer, but sits at offset 0; GT (the - strand) at offset 2
    std::vector<std::string> records { "ACGTT" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");

    Result suffix = count_of(*graph, ac, make_request(Scope::SUFFIX));
    EXPECT_EQ(Relation::EXACT, suffix.contexts->total.relation);
    EXPECT_EQ(0u, suffix.contexts->total.value);

    Result any = count_of(*graph, ac, make_request());
    EXPECT_EQ(2u, any.contexts->total.value);
    EXPECT_EQ(1u, any.contexts->by_offset.at(0).value);
    EXPECT_EQ(1u, any.contexts->by_offset.at(2).value);
    EXPECT_EQ(0u, any.contexts->suffix.value);
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::FORWARD).value);
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::REVERSE).value);
    EXPECT_LT(1u, any.work.ranges_visited);  // the flank is branching work
    check_against_oracles(*graph, ac, make_request(), &records);
}

TEST(PatternSearch, DoubleOccurrenceInOneKmer) {
    // ACGAC contains AC at offsets 0 and 3: two contexts of one k-mer, in offset order
    std::vector<std::string> records { "ACGAC" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    auto contexts = contexts_of(*graph, ac, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    ASSERT_EQ(2u, contexts.size());
    EXPECT_EQ(contexts[0].node, contexts[1].node);
    EXPECT_EQ(0u, contexts[0].offset);
    EXPECT_EQ(3u, contexts[1].offset);
    check_against_oracles(*graph, ac, make_request(), &records);
}

TEST(PatternSearch, RepeatedKmerOneContext) {
    // AAA occurs five times in the record and is one k-mer: AA has two contexts in it
    // (offsets 0 and 1), AAA one; graph contexts are not occurrences (§3)
    std::vector<std::string> records { "AAAAAAA" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    Result aa = count_of(*graph, Pattern::parse(PatternKind::DNA, "AA"),
                         make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(2u, aa.contexts->total.value);
    Result aaa = count_of(*graph, Pattern::parse(PatternKind::DNA, "AAA"),
                          make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(1u, aaa.contexts->total.value);

    // a homopolymer is what sdust flags with the seeder's parameters: stated, not refused
    Result poly = count_of(*graph, Pattern::parse(PatternKind::DNA, std::string(40, 'A')),
                           make_request());
    EXPECT_FALSE(poly.refusal);
    EXPECT_EQ(kNoteLowComplexity, poly.notes.front());
    check_against_oracles(*graph, Pattern::parse(PatternKind::DNA, "AA"), make_request(),
                          &records);
}

TEST(PatternSearch, PalindromeCountedOnce) {
    // ACGT = rc(ACGT): searched once, orientation palindromic, whatever strands say
    std::vector<std::string> records { "TTACGTAA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    Pattern acgt = Pattern::parse(PatternKind::DNA, "ACGT");
    for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
        Result r = count_of(*graph, acgt, make_request(Scope::ANY_OFFSET, strands));
        EXPECT_TRUE(r.palindromic);
        EXPECT_EQ(std::vector<Orientation>{ Orientation::PALINDROMIC }, r.searched);
        EXPECT_EQ(2u, r.contexts->total.value);  // TACGT at 1, ACGTA at 0
        EXPECT_EQ(2u, r.contexts->by_orientation.at(Orientation::PALINDROMIC).value);
        EXPECT_EQ(1u, r.contexts->by_orientation.size());
    }
    check_against_oracles(*graph, acgt, make_request(), &records);
}

TEST(PatternSearch, IUPACBothOrientationsOneOffset) {
    // NA and its reverse complement TN both match TA: two contexts at one (node, offset)
    // that differ in orientation (§3: orientation is part of a context's identity)
    std::vector<std::string> records { "GGTAGG" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);
    Pattern na = Pattern::parse(PatternKind::IUPAC, "NA");
    auto contexts = contexts_of(*graph, na, make_request());
    size_t pairs = 0;
    for (size_t i = 1; i < contexts.size(); ++i) {
        if (contexts[i].node == contexts[i - 1].node
                && contexts[i].offset == contexts[i - 1].offset) {
            EXPECT_EQ(Orientation::FORWARD, contexts[i - 1].orientation);
            EXPECT_EQ(Orientation::REVERSE, contexts[i].orientation);
            ++pairs;
        }
    }
    EXPECT_EQ(3u, pairs);  // GGTA at 2, GTAG at 1, TAGG at 0
    check_against_oracles(*graph, na, make_request(), &records);
}

TEST(PatternSearch, SuffixOnWrappedPrimary) {
    // the wrapped PRIMARY graph of ACGA exposes ACGA and TCGT; CGT is in TCGT, a virtual
    // node whose stored k-mer has ACG as its prefix, not CGT as a suffix: suffix is refused,
    // any_offset finds it (§4.1)
    std::vector<std::string> records { "ACGA" };
    auto graph = build(4, records, DeBruijnGraph::PRIMARY);
    Pattern cgt = Pattern::parse(PatternKind::DNA, "CGT");

    Result suffix = count_of(*graph, cgt, make_request(Scope::SUFFIX));
    ASSERT_TRUE(suffix.refusal);
    EXPECT_EQ("scope_unsupported", suffix.refusal->code);
    EXPECT_FALSE(suffix.contexts);

    Result any = count_of(*graph, cgt, make_request());
    EXPECT_EQ(Relation::EXACT, any.contexts->total.relation);
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::FORWARD).value);   // TCGT
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::REVERSE).value);   // ACGA
    EXPECT_EQ(std::vector<std::string>{ kNoteStrandUnknown }, any.notes);

    auto contexts = contexts_of(*graph, cgt, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    ASSERT_EQ(1u, contexts.size());
    EXPECT_EQ("TCGT", graph->get_node_sequence(contexts[0].node));
    EXPECT_EQ(1u, contexts[0].offset);
    check_against_oracles(*graph, cgt, make_request(), &records, DeBruijnGraph::PRIMARY);
}

TEST(PatternSearch, PalindromicAnchorOnWrappedPrimary) {
    // k = 4, ACGTA: the anchor window ACGT is a palindromic k-mer, found by the search of
    // ACGT and by that of rc(ACGT) = ACGT, mapped to itself: one anchor, not two (§4.1)
    std::vector<std::string> records { "ACGTA" };
    auto graph = build(4, records, DeBruijnGraph::PRIMARY);
    Pattern acgta = Pattern::parse(PatternKind::DNA, "ACGTA");
    Result r = count_of(*graph, acgta, make_request());
    ASSERT_TRUE(r.anchors);
    EXPECT_EQ(Scope::LONG, r.scope);
    EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
    EXPECT_EQ(1u, r.anchors->by_orientation.at(Orientation::FORWARD).value);   // ACGT
    EXPECT_EQ(1u, r.anchors->by_orientation.at(Orientation::REVERSE).value);   // TACG
    EXPECT_EQ(2u, r.anchors->total.value);
    EXPECT_EQ(Relation::UNKNOWN, r.anchors->paths.relation);
    check_against_oracles(*graph, acgta, make_request(), &records, DeBruijnGraph::PRIMARY);

    // an IUPAC window with one palindromic instance (ACGT) and one not (ACGA)
    std::vector<std::string> mixed { "ACGTA", "ACGAC" };
    auto graph2 = build(4, mixed, DeBruijnGraph::PRIMARY);
    Pattern acgwa = Pattern::parse(PatternKind::IUPAC, "ACGWA");
    Result m = count_of(*graph2, acgwa, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, m.anchors->total.relation);
    EXPECT_EQ(2u, m.anchors->total.value);  // ACGT once, ACGA
    check_against_oracles(*graph2, acgwa, make_request(), &mixed, DeBruijnGraph::PRIMARY);
    // and every short IUPAC pattern on it, to exercise the palindrome checks per offset
    for (std::string p : { "CG", "ACG", "GT", "NN", "WS", "CGT", "ACGT", "RCGY", "N" }) {
        check_against_oracles(*graph2, Pattern::parse(PatternKind::IUPAC, p), make_request(),
                              &mixed, DeBruijnGraph::PRIMARY);
    }
}

TEST(PatternSearch, ReverseAnchorOnWrappedPrimary) {
    // k = 3, AACG: its anchor window AAC exists only as the reverse complement of stored
    // GTT (records CGTT): found through rc(AAC), not through rc(AACG)[0, 3) = CGT
    std::vector<std::string> records { "CGTT" };
    auto graph = build(3, records, DeBruijnGraph::PRIMARY);
    Pattern aacg = Pattern::parse(PatternKind::DNA, "AACG");
    Result r = count_of(*graph, aacg, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(1u, r.anchors->total.value);
    check_against_oracles(*graph, aacg, make_request(), &records, DeBruijnGraph::PRIMARY);
}

TEST(PatternSearch, NativeCanonical) {
    // both orientations stored: P and rc(P) count alike, no mapping
    std::vector<std::string> records { "ACGTTGCAAGGCTTAC", "TTTACGGATC" };
    auto graph = build(5, records, DeBruijnGraph::CANONICAL);
    for (std::string p : { "AC", "GTT", "ACGG", "TAC", "G" }) {
        Pattern pattern = Pattern::parse(PatternKind::DNA, p);
        Result r = count_of(*graph, pattern, make_request());
        if (!pattern.is_palindromic()) {
            EXPECT_EQ(r.contexts->by_orientation.at(Orientation::FORWARD).value,
                      r.contexts->by_orientation.at(Orientation::REVERSE).value) << p;
        }
        for (Scope scope : { Scope::SUFFIX, Scope::ANY_OFFSET }) {
            check_against_oracles(*graph, pattern, make_request(scope), &records,
                                  DeBruijnGraph::CANONICAL);
        }
    }
}

TEST(PatternSearch, IUPACPatterns) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(6, records, mode);
        for (std::string p : { "ANR", "YNNA", "RY", "N", "NNNNNN", "GATYR", "ACGTTGC", "WSWS" }) {
            for (Scope scope : { Scope::SUFFIX, Scope::ANY_OFFSET }) {
                if (mode == DeBruijnGraph::PRIMARY && scope == Scope::SUFFIX)
                    continue;
                check_against_oracles(*graph, Pattern::parse(PatternKind::IUPAC, p),
                                      make_request(scope), &records, mode);
            }
        }
    }
}

TEST(PatternSearch, AnyOffsetVersusSuffix) {
    // the suffix count is any_offset's count at offset k - L
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(7, records, DeBruijnGraph::BASIC);
    for (std::string p : { "GAT", "AC", "TCGA", "C" }) {
        Pattern pattern = Pattern::parse(PatternKind::DNA, p);
        Result suffix = count_of(*graph, pattern, make_request(Scope::SUFFIX));
        Result any = count_of(*graph, pattern, make_request());
        EXPECT_EQ(suffix.contexts->total.value, any.contexts->suffix.value) << p;
        EXPECT_EQ(suffix.contexts->total.value,
                  any.contexts->by_offset.at(7 - pattern.length()).value) << p;
        EXPECT_GE(any.contexts->total.value, suffix.contexts->total.value);
        EXPECT_LE(suffix.work.steps, any.work.steps);
    }
}


// ---------------------------------------------------------------- stops and relations

TEST(PatternSearch, MaxStepsAtLeast) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    auto graph = build(7, records, DeBruijnGraph::BASIC);
    Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NA");
    Result full = count_of(*graph, pattern, make_request());
    ASSERT_EQ(Relation::EXACT, full.contexts->total.relation);

    for (uint64_t steps = 0; steps < full.work.steps; ++steps) {
        Result r = count_of(*graph, pattern, make_request(), steps);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);
        EXPECT_LE(r.work.steps, steps);
        const Count &total = r.contexts->total;
        // a discovery stop: at_least or (only scans interrupted) bounds, never exact
        EXPECT_NE(Relation::EXACT, total.relation);
        if (r.stop->phase == StopPhase::DISCOVERY) {
            EXPECT_TRUE(total.relation == Relation::AT_LEAST
                        || total.relation == Relation::UNKNOWN) << steps;
        }
        EXPECT_LE(total.value, full.contexts->total.value) << steps;
        // offsets whose discovery completed are exact and right; the others are lower bounds
        for (const auto &[p, count] : r.contexts->by_offset) {
            uint64_t truth = full.contexts->by_offset.at(p).value;
            if (count.relation == Relation::EXACT) {
                EXPECT_EQ(truth, count.value) << steps << " offset " << p;
            } else if (count.relation == Relation::BOUNDS) {
                EXPECT_LE(count.lower, truth);
                EXPECT_GE(count.upper, truth);
            } else if (count.relation == Relation::AT_LEAST) {
                EXPECT_LE(count.value, truth);
            }
        }
        for (const auto &[o, count] : r.contexts->by_orientation) {
            uint64_t truth = full.contexts->by_orientation.at(o).value;
            if (count.relation == Relation::EXACT) {
                EXPECT_EQ(truth, count.value);
            }
            EXPECT_LE(count.value, truth);
        }
    }
}

TEST(PatternSearch, RelationRuleAcrossPatterns) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern first = Pattern::parse(PatternKind::DNA, "GAT");
    Pattern second = Pattern::parse(PatternKind::DNA, "TTAC");

    Budget probe = unbounded_budget();
    Result alone = engine.count(first, make_request(), probe);
    Budget forward_budget = unbounded_budget();
    Result forward = engine.count(first, make_request(Scope::ANY_OFFSET, Strands::FORWARD),
                                  forward_budget);
    ASSERT_EQ(0u, forward.work.mask_scans);  // all its steps are discovery
    const uint64_t forward_steps = forward.work.steps;
    ASSERT_LT(forward_steps, alone.work.steps);

    // the budget ends right after the forward search: + exact, - not started (unknown),
    // the total at_least (§3: an unexplored strand has no upper bound)
    Budget budget = unbounded_budget(forward_steps);
    Result r = engine.count(first, make_request(), budget);
    EXPECT_EQ(Relation::EXACT, r.contexts->by_orientation.at(Orientation::FORWARD).relation);
    EXPECT_NE(Relation::EXACT, r.contexts->by_orientation.at(Orientation::REVERSE).relation);
    EXPECT_EQ(Relation::AT_LEAST, r.contexts->total.relation);
    EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);

    // the stop is sticky: the next pattern of the request answers unknown with that stop
    Result next = engine.count(second, make_request(), budget);
    ASSERT_TRUE(next.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, next.stop->phase);
    EXPECT_EQ(StopReason::MAX_STEPS, next.stop->reason);
    EXPECT_EQ(Relation::UNKNOWN, next.contexts->total.relation);
    for (const auto &[p, count] : next.contexts->by_offset) {
        EXPECT_EQ(Relation::UNKNOWN, count.relation);
    }
    EXPECT_EQ(0u, next.work.steps);

    Request partial = make_request();
    partial.mode = Mode::PARTIAL;
    std::vector<Ctx> got;
    Result withheld = engine.enumerate(second, partial, budget,
                                       [&](const Context &c) { got.push_back({ c.node, c.offset, c.orientation }); });
    EXPECT_TRUE(got.empty());
    EXPECT_EQ(StopReason::MAX_STEPS, withheld.extraction->cut);

    Request all = make_request();
    all.mode = Mode::ALL_OR_COUNT;
    Result none = engine.enumerate(second, all, budget, [&](const Context &) { FAIL(); });
    EXPECT_EQ(Withheld::DISCOVERY_BUDGET, none.extraction->withheld);
}

TEST(PatternSearch, StopAtThreshold) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    auto graph = build(7, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NA");
    std::vector<Ctx> expected = walk_oracle(*graph, pattern, make_request());
    ASSERT_LT(10u, expected.size());

    Request request = make_request();
    request.stop_at_threshold = true;
    request.max_contexts = 5;
    request.mode = Mode::ALL_OR_COUNT;
    Result all;
    EXPECT_TRUE(run_enumerate(engine, pattern, request, &all).empty());
    ASSERT_TRUE(all.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, all.stop->phase);
    EXPECT_EQ(StopReason::MAX_CONTEXTS, all.stop->reason);
    EXPECT_EQ(Relation::AT_LEAST, all.contexts->total.relation);
    EXPECT_GT(all.contexts->total.value, 5u);
    EXPECT_EQ(Withheld::THRESHOLD_CROSSED, all.extraction->withheld);

    request.mode = Mode::PARTIAL;
    Result partial;
    std::vector<Ctx> some = run_enumerate(engine, pattern, request, &partial);
    EXPECT_EQ(5u, some.size());
    EXPECT_EQ(StopReason::MAX_CONTEXTS, partial.extraction->cut);
    EXPECT_TRUE(std::is_sorted(some.begin(), some.end()));
    for (const Ctx &c : some) {
        EXPECT_TRUE(std::binary_search(expected.begin(), expected.end(), c));
    }

    // a threshold stop ends only its own pattern
    Budget budget = unbounded_budget();
    request.mode = Mode::COUNT;
    engine.count(pattern, request, budget);
    EXPECT_FALSE(budget.stopped());
    Result next = engine.count(Pattern::parse(PatternKind::DNA, "GGGCCC"), request, budget);
    EXPECT_FALSE(next.stop);
    EXPECT_EQ(Relation::EXACT, next.contexts->total.relation);
}

TEST(PatternSearch, PartialAfterStepStop) {
    // PARTIAL delivers what was discovered (a prefix-free subset of the truth, sorted), the
    // cut stated; ALL_OR_COUNT delivers nothing
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    auto graph = build(7, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NR");
    std::vector<Ctx> expected = walk_oracle(*graph, pattern, make_request());

    Request request = make_request();
    request.mode = Mode::PARTIAL;
    for (uint64_t steps : { 5, 20, 60 }) {
        Result r;
        std::vector<Ctx> some = run_enumerate(engine, pattern, request, &r, steps);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopReason::MAX_STEPS, r.extraction->cut);
        EXPECT_FALSE(r.extraction->complete);
        EXPECT_TRUE(std::is_sorted(some.begin(), some.end()));
        EXPECT_LE(some.size(), expected.size());
        for (const Ctx &c : some) {
            EXPECT_TRUE(std::binary_search(expected.begin(), expected.end(), c));
        }
    }
}

TEST(PatternSearch, DeadlineStops) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::DNA, "GAT");

    // a clock already past the work time: the pattern does not start
    auto start = Deadline::Clock::now();
    auto late = [start]() { return start + std::chrono::seconds(10); };
    {
        Budget budget(kManySteps, Deadline(start, 1000, 250, late));
        Result r = engine.count(pattern, make_request(), budget);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_TRUE(r.time_limited);
        EXPECT_EQ(Relation::UNKNOWN, r.contexts->total.relation);
    }
    {
        Budget budget(kManySteps, Deadline(start, 1000, 250, late));
        Request all = make_request();
        all.mode = Mode::ALL_OR_COUNT;
        Result r = engine.enumerate(pattern, all, budget, [](const Context &) { FAIL(); });
        EXPECT_EQ(Withheld::DEADLINE, r.extraction->withheld);
    }
    {
        Budget budget(kManySteps, Deadline(start, 1000, 250, late));
        Request partial = make_request();
        partial.mode = Mode::PARTIAL;
        Result r = engine.enumerate(pattern, partial, budget, [](const Context &) { FAIL(); });
        EXPECT_EQ(StopReason::TIME, r.extraction->cut);
        EXPECT_TRUE(r.time_limited);
    }
    {
        // the work time passes between discovery and the release: counts exact, the
        // release withheld with stop {extraction, time}
        int readings = 0;
        auto clock = [start, &readings]() {
            return start + std::chrono::milliseconds(++readings > 1 ? 10'000 : 0);
        };
        Budget budget(kManySteps, Deadline(start, 1000, 250, clock));
        Request all = make_request();
        all.mode = Mode::ALL_OR_COUNT;
        Result r = engine.enumerate(pattern, all, budget, [](const Context &) { FAIL(); });
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopPhase::EXTRACTION, r.stop->phase);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_EQ(Withheld::DEADLINE, r.extraction->withheld);
        EXPECT_TRUE(r.time_limited);
    }
}

TEST(PatternSearch, Refusals) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Request request = make_request();
    request.min_information_bits = 24;

    Budget budget = unbounded_budget();
    // an exact pattern in suffix scope is always admitted, however short
    Result suffix = engine.count(Pattern::parse(PatternKind::DNA, "GA"),
                                 make_request(Scope::SUFFIX), budget);
    EXPECT_FALSE(suffix.refusal);
    request.scope = Scope::SUFFIX;
    suffix = engine.count(Pattern::parse(PatternKind::DNA, "GA"), request, budget);
    EXPECT_FALSE(suffix.refusal);

    // any_offset, IUPAC and long patterns are gated (§5.3), and refused before any step
    request.scope = Scope::ANY_OFFSET;
    Result refused = engine.count(Pattern::parse(PatternKind::DNA, "GATC"), request, budget);
    ASSERT_TRUE(refused.refusal);
    EXPECT_EQ("information_below_floor", refused.refusal->code);
    EXPECT_EQ(0u, refused.work.steps);
    request.scope = Scope::SUFFIX;
    refused = engine.count(Pattern::parse(PatternKind::IUPAC, "GATN"), request, budget);
    ASSERT_TRUE(refused.refusal);
    refused = engine.count(Pattern::parse(PatternKind::DNA, "ACGTTGCAAGGCTTAC"), request, budget);
    ASSERT_TRUE(refused.refusal);  // 12 bits in its anchor window at k = 6
    EXPECT_NE(std::string::npos, refused.refusal->message.find("anchor window"));
    EXPECT_DOUBLE_EQ(12.0, *refused.anchor_information_bits);

    Request count = make_request();
    count.mode = Mode::COUNT;
    EXPECT_THROW(engine.enumerate(Pattern::parse(PatternKind::DNA, "GA"), count, budget,
                                  [](const Context &) {}),
                 std::invalid_argument);
}

TEST(PatternSearch, LongPatternAnchors) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(5, records, mode);
        for (std::string p : { "ACGTTGC", "GATCGATC", "NNNNNNA", "GGGCCCAA", "TTTTTTTT" }) {
            Pattern pattern = Pattern::parse(PatternKind::IUPAC, p);
            check_against_oracles(*graph, pattern, make_request(), &records, mode);
            Result r = count_of(*graph, pattern, make_request());
            EXPECT_EQ(Scope::LONG, r.scope);
            EXPECT_TRUE(r.anchor_information_bits);
            EXPECT_EQ(kNotePathsLater, r.notes.back());
        }
    }
}


// ---------------------------------------------------------------- even k, wrapped PRIMARY

TEST(PatternSearch, EvenKPrimaryPalindromeScanBounds) {
    // k = 6: palindromic k-mers (ACGCGT, AATATT, ...) are stored once and found by both
    // probes; the union subtracts them after a palindrome check charged as a scan.
    // Interrupting that check leaves valid bounds
    std::vector<std::string> records { "TTACGCGTAA", "GAATATTCCG", "ACGCGTACGT", "CCGGAATTCC" };
    auto graph = build(6, records, DeBruijnGraph::PRIMARY);
    for (std::string p : { "CG", "AT", "GCG", "N", "AATT" }) {
        Pattern pattern = Pattern::parse(PatternKind::IUPAC, p);
        check_against_oracles(*graph, pattern, make_request(), &records, DeBruijnGraph::PRIMARY);

        Result full = count_of(*graph, pattern, make_request());
        ASSERT_EQ(Relation::EXACT, full.contexts->total.relation);
        uint64_t truth = full.contexts->total.value;
        for (uint64_t steps = full.work.ranges_visited; steps < full.work.steps; ++steps) {
            Result r = count_of(*graph, pattern, make_request(), steps);
            ASSERT_TRUE(r.stop);
            EXPECT_EQ(StopPhase::MASK_SCAN, r.stop->phase);
            const Count &total = r.contexts->total;
            EXPECT_EQ(Relation::BOUNDS, total.relation) << p << " " << steps;
            EXPECT_LE(total.lower, truth) << p << " " << steps;
            EXPECT_GE(total.upper, truth) << p << " " << steps;
        }
    }
}


// ---------------------------------------------------------------- randomised

TEST(PatternSearch, RandomGraphsAgainstOracles) {
    // fixed seeds: a few hundred (graph, pattern, request) cases over every graph mode, both
    // builders, both scopes, every strand choice, DNA and IUPAC patterns of length 1 .. k + 2
    const char *iupac = "ACGTRYSWKMBDHVN";
    size_t cases = 0;
    size_t nonempty = 0;
    size_t multiple = 0;
    size_t max_expected = 0;
    for (uint32_t seed = 1; seed <= 60; ++seed) {
        std::mt19937 rng(seed);
        size_t k = 3 + rng() % 6;
        auto mode = static_cast<DeBruijnGraph::Mode>(rng() % 3);
        bool batch = rng() % 2;

        std::vector<std::string> records(1 + rng() % 5);
        for (std::string &record : records) {
            size_t length = 3 + rng() % 28;
            for (size_t i = 0; i < length; ++i) {
                // an occasional N splits a record into islands (DNA4)
                record.push_back(rng() % 25 ? "ACGT"[rng() % 4] : 'N');
            }
        }
        auto graph = build(k, records, mode, batch);

        for (int t = 0; t < 6; ++t) {
            size_t length = 1 + rng() % (k + 2);
            bool use_iupac = rng() % 2;
            std::string text;
            for (size_t i = 0; i < length; ++i) {
                text.push_back(use_iupac && rng() % 3 == 0 ? iupac[rng() % 15] : "ACGT"[rng() % 4]);
            }
            // two thirds of the patterns from the records, so that most have contexts
            if (rng() % 3) {
                const std::string &record = records[rng() % records.size()];
                size_t from = rng() % record.size();
                std::string piece = record.substr(from, length);
                if (piece.size() == length && piece.find('N') == std::string::npos)
                    text = piece;
            }
            Scope scope = mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0
                ? Scope::SUFFIX
                : Scope::ANY_OFFSET;
            auto strands = static_cast<Strands>(rng() % 3);
            Pattern pattern = Pattern::parse(PatternKind::IUPAC, text);
            SCOPED_TRACE("seed " + std::to_string(seed) + " k " + std::to_string(k)
                         + " mode " + std::to_string(mode) + " batch " + std::to_string(batch)
                         + " pattern " + text);
            size_t num_expected = 0;
            check_against_oracles(*graph, pattern, make_request(scope, strands), &records, mode,
                                  &num_expected);
            ++cases;
            nonempty += num_expected > 0;
            multiple += num_expected > 3;
            max_expected = std::max(max_expected, num_expected);
        }
    }
    EXPECT_EQ(360u, cases);
    // the cases are not vacuous: most have contexts, many several
    EXPECT_LT(180u, nonempty);
    EXPECT_LT(90u, multiple);
    std::cerr << "random cases: " << cases << ", with contexts: " << nonempty
              << ", with more than 3: " << multiple << ", most: " << max_expected << std::endl;
}

#endif // ! _PROTEIN_GRAPH

} // namespace
