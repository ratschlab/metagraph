/**
 * The oracle suite of the pattern search (docs/DESIGN-pattern-search.md, the counts and lists
 * and the engine half of the extension): tiny graphs built from explicit records, checked
 * against two brute-force oracles:
 *  - the graph-walk oracle enumerates the built graph's k-mers (every valid edge of the
 *    DBGSuccinct; on a wrapped PRIMARY graph also the wrapper's reverse complements, numbered
 *    by CanonicalDBG::reverse_complement) and lists every (node, offset, orientation) whose
 *    k-mer instantiates the oriented pattern at that offset: the expected graph contexts,
 *    with their node ids and in answer order;
 *  - the record-scan oracle scans the records' retained islands (the stretches over the
 *    symbols the build indexes: A, C, G, T, and N on a DNA5 build) as strings on both strands
 *    and lists the distinct (k-mer, offset, orientation) they contain: what the records say
 *    the contexts must be, independently of the graph.
 * count() must equal the first in every unit and relation, enumerate() must return exactly
 * its list in its order, and the first must equal the second wherever no k-mer was pruned.
 *
 * What the oracles share with the engine, and what they do not. They read a pattern only
 * through this file's own IUPAC table (iupac_bases: the bases of each code, its reverse
 * complement and palindromy derived from it), never through Pattern::positions(),
 * reverse_complement() or is_palindromic(); ParsePatterns checks those against the table,
 * code by code. They do share the graph: the graph-walk oracle reads the DBGSuccinct's
 * valid-edge mask through in_graph() and the node ids of the graph and of the CanonicalDBG
 * wrapper, as the engine does; the record-scan oracle, which reads neither, is the check on
 * that. Context::base_node is checked against the stored node found by spelling every stored
 * k-mer (StoredNodes), not against the engine's or CanonicalDBG's id arithmetic.
 *
 * Alphabet: the file runs on the DNA (DNA4) and DNA5 builds. Every DNA5-specific
 * expectation is under #if _DNA5_GRAPH. CI builds DNA and Protein only; the DNA5 branches
 * were run by hand on a DNA5 build of the engine and these tests (every PatternSearch and
 * PatternSearchFixes test passed), not in CI.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"
#include "pattern_test_support.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "graph/alignment/aligner_seeder_methods.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

// the nucleotide builds the engine serves ($ACGT and $ACGTN); the case-sensitive DNA build
// ($ACGTNacgt) is refused by the engine (alphabet_unsupported) and Protein has no pattern
// search, so neither runs this suite
#if _DNA_GRAPH || _DNA5_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;
using mtg::test::build_graph_batch;
using mtg::test::unbounded_deadline;

typedef DeBruijnGraph::node_index node_index;

constexpr uint64_t kManySteps = 1'000'000'000;


// ---------------------------------------------------------------- helpers

std::string rev_comp(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

// the record symbols the build indexes (§3): A, C, G, T, and N on a DNA5 build; any other
// symbol splits a record into islands
bool indexed_symbol(char c) {
#if _DNA5_GRAPH
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
#else
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
#endif
}

bool indexed(std::string_view s) {
    return std::all_of(s.begin(), s.end(), indexed_symbol);
}


// ---------------------------------------------------------------- the oracles' IUPAC table

/**
 * The IUPAC-IUB nucleotide codes, written out here as the bases each admits (sorted): the
 * oracles' only reading of a pattern. Nothing below derives a base set, a reverse complement
 * or palindromy from Pattern; Pattern is the engine's input, and ParsePatterns checks the
 * engine's tables against this one.
 */
constexpr char kIUPACCodes[] = "ACGTRYSWKMBDHVN";

std::string iupac_bases(char code) {
    switch (std::toupper(static_cast<unsigned char>(code))) {
        case 'A': return "A";
        case 'C': return "C";
        case 'G': return "G";
        case 'T': return "T";
        case 'R': return "AG";    // purine
        case 'Y': return "CT";    // pyrimidine
        case 'S': return "CG";    // strong
        case 'W': return "AT";    // weak
        case 'K': return "GT";    // keto
        case 'M': return "AC";    // amino
        case 'B': return "CGT";   // not A
        case 'D': return "AGT";   // not C
        case 'H': return "ACT";   // not G
        case 'V': return "ACG";   // not T
        case 'N': return "ACGT";  // any base, never the record symbol N (§3)
        default: return "";
    }
}

// the code admitting exactly |bases| (sorted)
char iupac_code(const std::string &bases) {
    for (const char *c = kIUPACCodes; *c; ++c) {
        if (iupac_bases(*c) == bases)
            return *c;
    }
    return '?';
}

char complement_base(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return '?';
    }
}

// an oriented pattern as the oracles read it: per position, the bases it admits (sorted)
typedef std::vector<std::string> Bases;

Bases oracle_pattern(std::string_view text) {
    Bases q;
    for (char c : text) {
        q.push_back(iupac_bases(c));
        EXPECT_FALSE(q.back().empty()) << "not an IUPAC code: " << c;
    }
    return q;
}

Bases oracle_pattern(const Pattern &pattern) {
    return oracle_pattern(pattern.text());
}

// rc(q): the positions reversed, every base of each complemented
Bases oracle_reverse_complement(const Bases &q) {
    Bases rc(q.rbegin(), q.rend());
    for (std::string &bases : rc) {
        for (char &base : bases) {
            base = complement_base(base);
        }
        std::sort(bases.begin(), bases.end());
    }
    return rc;
}

bool oracle_palindromic(const Bases &q) {
    return q == oracle_reverse_complement(q);
}

std::string oracle_text(const Bases &q) {
    std::string text;
    for (const std::string &bases : q) {
        text.push_back(iupac_code(bases));
    }
    return text;
}

// |s| instantiates |q|: every symbol is one its position admits. A position admits A, C, G,
// T only, never the record symbol N or the sentinel $ (§3); the flank is not a position
bool matches(const Bases &q, std::string_view s) {
    if (s.size() != q.size())
        return false;
    for (size_t i = 0; i < q.size(); ++i) {
        if (q[i].find(s[i]) == std::string::npos)
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
    return Budget(max_steps, unbounded_deadline());
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

// the orientations searched, with the oriented pattern of each (§3), from the oracles' table
std::vector<std::pair<Orientation, Bases>> orientations(const Pattern &pattern,
                                                        Strands strands) {
    const Bases q = oracle_pattern(pattern);
    if (oracle_palindromic(q))
        return { { Orientation::PALINDROMIC, q } };
    std::vector<std::pair<Orientation, Bases>> result;
    if (strands != Strands::REVERSE)
        result.emplace_back(Orientation::FORWARD, q);
    if (strands != Strands::FORWARD)
        result.emplace_back(Orientation::REVERSE, oracle_reverse_complement(q));
    return result;
}

// the oriented pattern a context or path of |orientation| instantiates
Bases oriented(const Pattern &pattern, Orientation orientation) {
    const Bases q = oracle_pattern(pattern);
    return orientation == Orientation::REVERSE ? oracle_reverse_complement(q) : q;
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

// the first min(L, k) positions: the whole pattern for L <= k, the anchor window for L > k
Bases window(const Bases &q, size_t k) {
    return Bases(q.begin(), q.begin() + std::min(q.size(), k));
}

/**
 * The stored node of every node of the served graph, found by spelling: on a wrapped
 * PRIMARY graph the stored k-mer whose sequence is the node's k-mer or its reverse complement
 * (every stored k-mer of the unwrapped DBGSuccinct is spelled once); on BASIC and CANONICAL
 * the node itself. Never the engine's base_node() or CanonicalDBG's id arithmetic.
 */
class StoredNodes {
  public:
    explicit StoredNodes(const DeBruijnGraph &graph)
          : graph_(graph), wrapped_(dynamic_cast<const CanonicalDBG*>(&graph)) {
        if (!wrapped_)
            return;
        const DBGSuccinct &dbg_succ = dynamic_cast<const DBGSuccinct&>(wrapped_->get_graph());
        for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
            if (dbg_succ.in_graph(y))
                stored_.emplace(dbg_succ.get_node_sequence(y), y);
        }
    }

    node_index of(node_index node) const {
        if (!wrapped_)
            return node;
        const std::string kmer = graph_.get_node_sequence(node);
        auto it = stored_.find(kmer);
        if (it == stored_.end())
            it = stored_.find(rev_comp(kmer));
        return it == stored_.end() ? DeBruijnGraph::npos : it->second;
    }

  private:
    const DeBruijnGraph &graph_;
    const CanonicalDBG *wrapped_;
    std::map<std::string, node_index> stored_;
};


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
            Bases w = window(q, k);
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
// BASIC) and the contexts they spell. An island is a stretch of indexed symbols: on DNA4 an
// N splits a record; on DNA5 the k-mers containing N are kept, and their N can lie in a flank
// but never at a pattern position (matches())
SpelledContexts record_oracle(const std::vector<std::string> &records, size_t k,
                              DeBruijnGraph::Mode mode, const Pattern &pattern,
                              const Request &request) {
    std::set<std::string> kmers;
    for (const std::string &record : records) {
        for (size_t i = 0; i + k <= record.size(); ++i) {
            std::string kmer = record.substr(i, k);
            if (!indexed(kmer))
                continue;
            kmers.insert(kmer);
            if (mode != DeBruijnGraph::BASIC)
                kmers.insert(rev_comp(kmer));
        }
    }
    SpelledContexts result;
    for (const std::string &kmer : kmers) {
        for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
            Bases w = window(q, k);
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

std::vector<Ctx> run_enumerate(const PatternSearch &engine, const StoredNodes &stored,
                               const Pattern &pattern, const Request &request, Result *result,
                               uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    std::vector<Ctx> contexts;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        // the node of the annotation row (the route's row is graph_to_anno_index(base_node)):
        // the stored k-mer of the context, found by spelling
        EXPECT_EQ(stored.of(c.node), c.base_node) << "node " << c.node;
        EXPECT_TRUE(c.path.empty());
        EXPECT_TRUE(c.sequence.empty());
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
    const StoredNodes stored(graph);

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
        // paths need extend_paths: unknown, unless no anchor exists
        EXPECT_EQ(expected.empty() ? Relation::EXACT : Relation::UNKNOWN,
                  counted.anchors->paths.relation);
    }

    // enumerate(): discovery identical to count(), then every context in answer order
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result enumerated;
    std::vector<Ctx> released = run_enumerate(engine, stored, pattern, all, &enumerated);
    EXPECT_EQ(counted.work.steps, enumerated.work.steps);
    EXPECT_EQ(counted.work.ranges_visited, enumerated.work.ranges_visited);
    ASSERT_TRUE(enumerated.extraction);
    if (L > k && !request.release_anchors) {
        // the results of a long pattern are its paths (extend_paths): withheld
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
        EXPECT_TRUE(matches(window(oriented(pattern, c.orientation), k),
                            std::string_view(kmer).substr(c.offset, L)));
    }

    // ALL_OR_COUNT: all or nothing at the threshold (max_anchors for released anchors)
    auto set_cap = [&](Request *r, uint64_t cap) {
        (L > k ? r->max_anchors : r->max_contexts) = cap;
    };
    if (!expected.empty()) {
        Request tight = all;
        set_cap(&tight, expected.size() - 1);
        Result withheld;
        EXPECT_TRUE(run_enumerate(engine, stored, pattern, tight, &withheld).empty());
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
        std::vector<Ctx> prefix = run_enumerate(engine, stored, pattern, partial, &cut);
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
    return run_enumerate(PatternSearch(graph), StoredNodes(graph), pattern, all, &result);
}

Result count_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    return PatternSearch(graph).count(pattern, request, budget);
}


// ---------------------------------------------------------------- relations against the truth

/**
 * The true counts of one (graph, pattern, request), from the graph-walk oracle: contexts
 * (L <= k) or anchors (L > k), in total, per offset of the scope and per orientation
 * searched, zeros included.
 */
struct Truth {
    uint64_t total = 0;
    std::map<uint32_t, uint64_t> by_offset;
    std::map<Orientation, uint64_t> by_orientation;
};

Truth truth_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
               const std::vector<Ctx> &expected) {
    Truth truth;
    for (uint32_t p : scope_offsets(pattern.length(), graph.get_k(), request.scope)) {
        truth.by_offset[p] = 0;
    }
    for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
        truth.by_orientation[orientation] = 0;
    }
    for (const Ctx &c : expected) {
        ++truth.total;
        ++truth.by_offset[c.offset];
        ++truth.by_orientation[c.orientation];
    }
    return truth;
}

// |count| states a true relation to |truth| (§3): EXACT equal to it, AT_LEAST at most it,
// BOUNDS around it with the lower bound as the value, UNKNOWN without a value
::testing::AssertionResult true_relation(const Count &count, uint64_t truth) {
    bool holds = false;
    switch (count.relation) {
        case Relation::EXACT:
            holds = count.value == truth;
            break;
        case Relation::AT_LEAST:
            holds = count.value <= truth;
            break;
        case Relation::BOUNDS:
            holds = count.lower == count.value && count.lower <= truth && truth <= count.upper;
            break;
        case Relation::UNKNOWN:
            holds = count.value == 0;
            break;
    }
    if (holds)
        return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure()
        << to_string(count.relation) << " " << count.value << " [" << count.lower << ", "
        << count.upper << "] does not hold for the true count " << truth;
}

// every count of |r| (total, suffix, per offset, per orientation) holds for |truth|
void expect_true_relations(const Result &r, const Truth &truth, size_t k, size_t L) {
    if (L > k) {
        ASSERT_TRUE(r.anchors);
        EXPECT_TRUE(true_relation(r.anchors->total, truth.total)) << "anchors";
        ASSERT_EQ(truth.by_orientation.size(), r.anchors->by_orientation.size());
        for (const auto &[o, count] : r.anchors->by_orientation) {
            EXPECT_TRUE(true_relation(count, truth.by_orientation.at(o)))
                << "anchors " << orientation_key(o);
        }
        return;
    }
    ASSERT_TRUE(r.contexts);
    EXPECT_TRUE(true_relation(r.contexts->total, truth.total)) << "total";
    EXPECT_TRUE(true_relation(r.contexts->suffix, truth.by_offset.at(k - L))) << "suffix";
    ASSERT_EQ(truth.by_offset.size(), r.contexts->by_offset.size());
    for (const auto &[p, count] : r.contexts->by_offset) {
        EXPECT_TRUE(true_relation(count, truth.by_offset.at(p))) << "offset " << p;
    }
    ASSERT_EQ(truth.by_orientation.size(), r.contexts->by_orientation.size());
    for (const auto &[o, count] : r.contexts->by_orientation) {
        EXPECT_TRUE(true_relation(count, truth.by_orientation.at(o)))
            << "orientation " << orientation_key(o);
    }
}

/**
 * Every step budget below the complete run's: count() stopped at max_steps = s for
 * every s in [0, steps of the complete run) states a true relation in every count against
 * the graph-walk oracle; the total is AT_LEAST after a stop in discovery (s below the
 * complete run's ranges_visited) and BOUNDS after a stop in a mask scan, never EXACT; the
 * stop is {that phase, max_steps} after exactly s steps. For L <= k, enumerate() at the same
 * budget counts alike and delivers, in PARTIAL, a sorted subset of the oracle's contexts with
 * the cut max_steps, and in ALL_OR_COUNT nothing (discovery_budget); enumerate() is run at
 * every third budget and at the budgets around the end of discovery. A case whose complete
 * run takes more than |max_sweep| steps is sampled: the first and last 64 budgets, the 64
 * around the end of discovery, and about |max_sweep| more evenly spread. Returns the number
 * of stopped runs checked.
 */
uint64_t check_halts_against_oracles(const DeBruijnGraph &graph, const Pattern &pattern,
                                     const Request &request, uint64_t max_sweep = 1500) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    PatternSearch engine(graph);
    const StoredNodes stored(graph);
    const std::vector<Ctx> expected = walk_oracle(graph, pattern, request);
    const Truth truth = truth_of(graph, pattern, request, expected);

    Budget unbounded = unbounded_budget();
    const Result full = engine.count(pattern, request, unbounded);
    EXPECT_FALSE(full.refusal);
    if (full.refusal)
        return 0;
    EXPECT_FALSE(full.stop);
    expect_true_relations(full, truth, k, L);
    EXPECT_EQ(Relation::EXACT, (L > k ? full.anchors->total : full.contexts->total).relation);
    const uint64_t discovery = full.work.ranges_visited;
    const uint64_t steps = full.work.steps;

    std::set<uint64_t> budgets;
    if (steps <= max_sweep) {
        for (uint64_t s = 0; s < steps; ++s) {
            budgets.insert(s);
        }
    } else {
        auto add = [&](uint64_t from, uint64_t to) {
            for (uint64_t s = from; s < std::min(to, steps); ++s) {
                budgets.insert(s);
            }
        };
        add(0, 64);
        add(steps - 64, steps);
        add(discovery > 32 ? discovery - 32 : 0, discovery + 32);
        for (uint64_t s = 0; s < steps; s += steps / max_sweep + 1) {
            budgets.insert(s);
        }
    }

    uint64_t runs = 0;
    for (uint64_t s : budgets) {
        SCOPED_TRACE("max_steps " + std::to_string(s) + " of " + std::to_string(steps)
                     + " (discovery " + std::to_string(discovery) + ")");
        Budget budget = unbounded_budget(s);
        const Result r = engine.count(pattern, request, budget);
        ++runs;
        EXPECT_FALSE(r.refusal);
        if (!r.stop) {
            ADD_FAILURE() << "no stop";
            continue;
        }
        const StopPhase phase = s < discovery ? StopPhase::DISCOVERY : StopPhase::MASK_SCAN;
        EXPECT_EQ(phase, r.stop->phase);
        EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);
        EXPECT_FALSE(r.time_limited);
        // every charge is one step, and the refused one is not counted
        EXPECT_EQ(s, r.work.steps);
        EXPECT_EQ(std::min(s, discovery), r.work.ranges_visited);
        const Count &total = L > k ? r.anchors->total : r.contexts->total;
        EXPECT_EQ(phase == StopPhase::DISCOVERY ? Relation::AT_LEAST : Relation::BOUNDS,
                  total.relation);
        expect_true_relations(r, truth, k, L);

        // enumerate() at every third budget and around the phase boundaries (its discovery
        // is count()'s, checked above at every budget)
        const bool boundary = (s + 1 >= discovery && s <= discovery + 1) || s + 1 == steps;
        if (L > k || !(s % 3 == 0 || boundary))
            continue;

        // enumerate(): the same discovery, then what the mode delivers after a step stop
        Request partial = request;
        partial.mode = Mode::PARTIAL;
        partial.max_contexts = kManySteps;
        Result cut;
        std::vector<Ctx> some = run_enumerate(engine, stored, pattern, partial, &cut, s);
        EXPECT_EQ(r.work.steps, cut.work.steps);
        EXPECT_EQ(total.relation, cut.contexts->total.relation);
        EXPECT_EQ(total.value, cut.contexts->total.value);
        if (!cut.extraction) {
            ADD_FAILURE() << "no extraction";
            continue;
        }
        EXPECT_FALSE(cut.extraction->complete);
        EXPECT_FALSE(cut.extraction->withheld);
        EXPECT_EQ(StopReason::MAX_STEPS, cut.extraction->cut);
        EXPECT_EQ(some.size(), cut.extraction->returned);
        EXPECT_TRUE(std::adjacent_find(some.begin(), some.end(),
                                       [](const Ctx &a, const Ctx &b) { return !(a < b); })
                        == some.end()) << "not strictly in answer order";
        for (const Ctx &c : some) {
            EXPECT_TRUE(std::binary_search(expected.begin(), expected.end(), c)) << c;
        }

        Request all = request;
        all.mode = Mode::ALL_OR_COUNT;
        Result none;
        EXPECT_TRUE(run_enumerate(engine, stored, pattern, all, &none, s).empty());
        EXPECT_TRUE(none.extraction && none.extraction->withheld == Withheld::DISCOVERY_BUDGET
                        && !none.extraction->complete);
    }
    return runs;
}

/**
 * stop_at_threshold and the thresholds of the release at every cap from 0 to one past the
 * truth (sampled above 64): a threshold stop implies cap < value <= truth, AT_LEAST
 * in discovery with the reason max_contexts (max_anchors for L > k), and in ALL_OR_COUNT
 * nothing (threshold_crossed); no stop implies EXACT and right (the stop is not promised
 * whenever the truth exceeds the cap: the running lower bound may stay below it). Without
 * stop_at_threshold,
 * ALL_OR_COUNT releases exactly the oracle's list when the truth is within the cap and
 * withholds count_above_threshold otherwise, and PARTIAL releases the first cap contexts.
 */
void check_thresholds_against_oracles(const DeBruijnGraph &graph, const Pattern &pattern,
                                      const Request &request) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    PatternSearch engine(graph);
    const StoredNodes stored(graph);
    const std::vector<Ctx> expected = walk_oracle(graph, pattern, request);
    const Truth truth = truth_of(graph, pattern, request, expected);
    const uint64_t n = truth.total;

    std::set<uint64_t> caps;
    for (uint64_t cap = 0; cap <= n + 1; cap += (cap < 32 || cap + 33 > n) ? 1 : n / 32 + 1) {
        caps.insert(cap);
    }
    for (uint64_t cap : caps) {
        SCOPED_TRACE("cap " + std::to_string(cap) + " truth " + std::to_string(n));
        auto with_cap = [&](Mode mode, bool stop) {
            Request r = request;
            r.mode = mode;
            r.stop_at_threshold = stop;
            (L > k ? r.max_anchors : r.max_contexts) = cap;
            return r;
        };

        Budget budget = unbounded_budget();
        const Result c = engine.count(pattern, with_cap(Mode::COUNT, true), budget);
        const Count &total = L > k ? c.anchors->total : c.contexts->total;
        expect_true_relations(c, truth, k, L);
        if (c.stop) {
            EXPECT_EQ(StopPhase::DISCOVERY, c.stop->phase);
            EXPECT_EQ(L > k ? StopReason::MAX_ANCHORS : StopReason::MAX_CONTEXTS,
                      c.stop->reason);
            EXPECT_EQ(Relation::AT_LEAST, total.relation);
            EXPECT_LT(cap, total.value);
            EXPECT_FALSE(budget.stopped());  // a threshold ends only its own pattern
            // ... the next pattern of the request answers on the same budget
            const Result next = engine.count(pattern, with_cap(Mode::COUNT, false), budget);
            EXPECT_FALSE(next.stop);
            EXPECT_EQ(Relation::EXACT,
                      (L > k ? next.anchors->total : next.contexts->total).relation);
        } else {
            EXPECT_EQ(Relation::EXACT, total.relation);
        }
        if (cap >= n) {
            EXPECT_FALSE(c.stop);
        }

        if (L > k)
            continue;

        Result all;
        std::vector<Ctx> released
            = run_enumerate(engine, stored, pattern, with_cap(Mode::ALL_OR_COUNT, true), &all);
        if (c.stop) {
            EXPECT_TRUE(released.empty());
            EXPECT_EQ(Withheld::THRESHOLD_CROSSED, all.extraction->withheld);
        } else if (n <= cap) {
            EXPECT_EQ(expected, released);
            EXPECT_TRUE(all.extraction->complete);
        } else {
            // the running lower bound need not cross the cap (a wrapped PRIMARY union counts
            // the larger part only): no stop, and the exact count above the cap withholds
            EXPECT_TRUE(released.empty());
            EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, all.extraction->withheld);
        }

        Result plain;
        released = run_enumerate(engine, stored, pattern, with_cap(Mode::ALL_OR_COUNT, false),
                                 &plain);
        EXPECT_EQ(Relation::EXACT, plain.contexts->total.relation);
        if (n <= cap) {
            EXPECT_TRUE(plain.extraction->complete);
            EXPECT_EQ(expected, released);
        } else {
            EXPECT_TRUE(released.empty());
            EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, plain.extraction->withheld);
        }

        Result partial;
        released = run_enumerate(engine, stored, pattern, with_cap(Mode::PARTIAL, true),
                                 &partial);
        ASSERT_LE(released.size(), std::min(cap, n));
        EXPECT_TRUE(std::is_sorted(released.begin(), released.end()));
        for (const Ctx &x : released) {
            EXPECT_TRUE(std::binary_search(expected.begin(), expected.end(), x)) << x;
        }
        if (!partial.stop) {
            EXPECT_EQ(std::min(cap, n), released.size());
            EXPECT_TRUE(std::equal(released.begin(), released.end(), expected.begin()));
        }
    }
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

// the engine's base set of a position, spelled as bases through the header's bit convention
// (bit 0 A, bit 1 C, bit 2 G, bit 3 T)
std::string spelled_set(BaseSet set) {
    std::string bases;
    if (set & kBaseA) bases.push_back('A');
    if (set & kBaseC) bases.push_back('C');
    if (set & kBaseG) bases.push_back('G');
    if (set & kBaseT) bases.push_back('T');
    return bases;
}

TEST(PatternSearch, IUPACTableAgainstTheOracles) {
    // the engine's IUPAC table, code by code, against the oracles' own; a consistent swap of
    // two codes in both of the engine's tables (K/M, B/V, D/H) passes every pin of
    // ParsePatterns but not this
    for (const char *c = kIUPACCodes; *c; ++c) {
        SCOPED_TRACE(std::string("code ") + *c);
        const std::string bases = iupac_bases(*c);
        for (char text : { *c, static_cast<char>(std::tolower(*c)) }) {
            Pattern p = Pattern::parse(PatternKind::IUPAC, std::string(1, text));
            ASSERT_EQ(1u, p.length());
            EXPECT_EQ(std::string(1, *c), p.text());
            EXPECT_EQ(bases, spelled_set(p.positions()[0]));
            EXPECT_EQ(bases, spelled_set(p.allowed(0, "")));
            EXPECT_DOUBLE_EQ(std::log2(4.0 / bases.size()), p.information_bits());
            EXPECT_EQ(bases.size() == 1, p.is_exact());
            const Bases rc = oracle_reverse_complement({ bases });
            EXPECT_EQ(oracle_text(rc), p.reverse_complement().text());
            EXPECT_EQ(rc[0], spelled_set(p.reverse_complement().positions()[0]));
            EXPECT_EQ(oracle_palindromic({ bases }), p.is_palindromic());
        }
        // a DNA pattern is A, C, G, T only
        if (bases.size() == 1) {
            EXPECT_EQ(std::string(1, *c), Pattern::parse(PatternKind::DNA, std::string(1, *c)).text());
        } else {
            EXPECT_THROW(Pattern::parse(PatternKind::DNA, std::string(1, *c)), PatternError);
        }
    }
    // single codes, by the IUPAC rules themselves: S, W, N are their own complements
    for (char c : std::string("SWN")) {
        EXPECT_TRUE(Pattern::parse(PatternKind::IUPAC, std::string(1, c)).is_palindromic()) << c;
    }
    for (char c : std::string("ACGTRYKMBDHV")) {
        EXPECT_FALSE(Pattern::parse(PatternKind::IUPAC, std::string(1, c)).is_palindromic()) << c;
    }
    EXPECT_FALSE(Pattern::parse(PatternKind::DNA, "CAG").is_palindromic());
    EXPECT_FALSE(Pattern::parse(PatternKind::DNA, "TTAGGA").is_palindromic());
    EXPECT_TRUE(Pattern::parse(PatternKind::DNA, "GAATTC").is_palindromic());

    // every IUPAC pattern of length 1 .. 3 (3 615 of them): the reverse complement, its text
    // and palindromy, and every position, against the oracles' table
    size_t checked = 0;
    std::function<void(std::string&)> all = [&](std::string &text) {
        if (text.size()) {
            Pattern p = Pattern::parse(PatternKind::IUPAC, text);
            const Bases q = oracle_pattern(text);
            const Bases rc = oracle_reverse_complement(q);
            ASSERT_EQ(q.size(), p.positions().size()) << text;
            for (size_t i = 0; i < q.size(); ++i) {
                ASSERT_EQ(q[i], spelled_set(p.positions()[i])) << text << " " << i;
            }
            ASSERT_EQ(oracle_text(rc), p.reverse_complement().text()) << text;
            for (size_t i = 0; i < rc.size(); ++i) {
                ASSERT_EQ(rc[i], spelled_set(p.reverse_complement().positions()[i])) << text;
            }
            ASSERT_EQ(oracle_palindromic(q), p.is_palindromic()) << text;
            ASSERT_EQ(p.is_palindromic(), p.reverse_complement().is_palindromic()) << text;
            ++checked;
        }
        if (text.size() == 3)
            return;
        for (const char *c = kIUPACCodes; *c; ++c) {
            text.push_back(*c);
            all(text);
            text.pop_back();
        }
    };
    std::string text;
    all(text);
    EXPECT_EQ(15u + 15 * 15 + 15 * 15 * 15, checked);

    // the information of a long text, bit for bit the sum over its codes
    std::mt19937 rng(3);
    for (int i = 0; i < 200; ++i) {
        std::string t;
        for (size_t j = 0, n = 1 + rng() % 300; j < n; ++j) {
            t.push_back(kIUPACCodes[rng() % 15]);
        }
        double expected = 0;
        for (char c : t) {
            expected += std::log2(4.0 / iupac_bases(c).size());
        }
        EXPECT_EQ(expected, Pattern::parse(PatternKind::IUPAC, t).information_bits()) << t;
    }
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
    Budget budget(10, unbounded_deadline());
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
#if _DNA5_GRAPH
    EXPECT_EQ("$ACGTN", s.alphabet);
#else
    EXPECT_EQ("$ACGT", s.alphabet);
#endif
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

// the loop of suffix_to_prefix with its fixed symbol set (copied from
// aligner_seeder_methods.cpp): the default symbol set must reproduce it exactly
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

TEST(PatternSearch, FlankSymbolsAreTheBuildAlphabet) {
    // §4.1: the flank of any_offset tries every symbol of the graph's alphabet but the
    // sentinel, the seeder's NonSentinelSymbols (pattern_search.cpp, expand()): A, C, G, T on
    // a DNA4 build, and N as well on DNA5. On DNA4 that set equals the pattern alphabet, so
    // no DNA4 run can tell a flank restricted to A, C, G, T from the right one; only the
    // _DNA5_GRAPH cases (DNA5Flank, NCentredSelfComplementOddK) can, on a DNA5 build (run by
    // hand, not in CI). What DNA4 can check is the symbol set itself
    auto graph = build(4, { "ACGTACGT" }, DeBruijnGraph::BASIC);
    const boss::BOSS &boss = base_dbg(*graph).get_boss();
    std::string symbols;
    align::NonSentinelSymbols()(boss, std::make_tuple(uint64_t(1), uint64_t(1), size_t(1)),
                                [&](boss::BOSS::TAlphabet s) {
                                    symbols.push_back(boss.decode(s));
                                });
#if _DNA5_GRAPH
    EXPECT_EQ("ACGTN", symbols);
#else
    EXPECT_EQ("ACGT", symbols);
#endif
    EXPECT_EQ(boss.alphabet.substr(1), symbols);
}


// ---------------------------------------------------------------- §13: the named cases

// the k-mers of the graph's invalid (dummy or pruned) edges, spelled through the BOSS, $
// for the sentinel
std::set<std::string> invalid_edges(const DBGSuccinct &dbg_succ) {
    const boss::BOSS &boss = dbg_succ.get_boss();
    std::set<std::string> result;
    for (node_index e = 1; e <= dbg_succ.max_index(); ++e) {
        if (!dbg_succ.in_graph(e))
            result.insert(boss.get_node_str(e) + boss.decode(boss.get_W(e) % boss.alph_size));
    }
    return result;
}

// the edges, valid or not, of the nodes ending with |suffix| that can end with |c|: the
// W rule's leaf (§4.1)
DBGSuccinct::LastSymbolEdges leaf_edges(const DBGSuccinct &dbg_succ,
                                        const std::string &suffix, char c) {
    const boss::BOSS &boss = dbg_succ.get_boss();
    auto encoded = boss.encode(suffix);
    auto [first, last, end] = boss.index_range(encoded.begin(), encoded.end());
    EXPECT_EQ(suffix.size(), static_cast<size_t>(end - encoded.begin())) << suffix;
    first = boss.pred_last(first - 1) + 1;
    return dbg_succ.count_edges_with_last_symbol(first, last, boss.encode(c));
}

TEST(PatternSearch, PrunedDummyRecordStart) {
    // k = 4: ACGA is a k-mer of TACGAT, so the record ACGA starts without a dummy chain of
    // its own; AC at its start is found by any_offset, never through a dummy
    std::vector<std::string> records { "TACGAT", "ACGA" };
    for (bool batch : { false, true }) {
        auto graph = build(4, records, DeBruijnGraph::BASIC, batch);
        Pattern ac = Pattern::parse(PatternKind::DNA, "AC");

        // the premise, checked: no dummy chain $$$A, $$AC, $ACG of ACGA; TACGAT's
        // chain $$$T, $$TA, $TAC is there, and $TAC is the one dummy in the leaf of AC
        // (nodes ending with A, W = C), which only a scan can tell from a k-mer
        const DBGSuccinct &dbg_succ = base_dbg(*graph);
        const std::set<std::string> dummies = invalid_edges(dbg_succ);
        for (std::string chain : { "$$$A", "$$AC", "$ACG" }) {
            EXPECT_FALSE(dummies.count(chain)) << chain;
        }
        for (std::string chain : { "$$$T", "$$TA", "$TAC" }) {
            EXPECT_TRUE(dummies.count(chain)) << chain;
        }
        const auto leaf = leaf_edges(dbg_succ, "A", 'C');
        EXPECT_EQ(1u, leaf.candidates);
        EXPECT_EQ(1u, leaf.invalid_non_sentinel);

        // suffix: no k-mer ends with AC; the leaf's one candidate is the dummy $TAC
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

    // one edge scanned ($$AC, matching): bounds [3 - 1 - 1, 3 - 1] = [1, 2], exactly
    // ([1, 3] and [2, 2] would be wrong here)
    Result half = count_of(*graph, ac, request, discovery + 1);
    EXPECT_EQ(StopPhase::MASK_SCAN, half.stop->phase);
    EXPECT_EQ(Relation::BOUNDS, half.contexts->total.relation);
    EXPECT_EQ(1u, half.contexts->total.lower);
    EXPECT_EQ(2u, half.contexts->total.upper);
    EXPECT_EQ(1u, half.contexts->total.value);

    // a discovery stop: at_least, never bounds (an undiscovered branch has no upper bound)
    Result early = count_of(*graph, ac, request, discovery - 1);
    EXPECT_EQ(StopPhase::DISCOVERY, early.stop->phase);
    EXPECT_EQ(Relation::AT_LEAST, early.contexts->total.relation);
    EXPECT_LE(early.contexts->total.value, 2u);
}

TEST(PatternSearch, SentinelTighteningPinned) {
    // The bound of a W-rule leaf counts only the invalid edges whose W is not $ (J), since
    // a dummy sink (W = $) never carries the symbol c: the design's bound I tightened
    // (Item::lower). k = 3, GCA and TAC: the leaf of AC (nodes ending with A) holds the sink
    // CA$ and the k-mer TAC: one candidate (TAC), one invalid edge, none of it non-sentinel.
    // So AC has exactly one suffix context known by ranks alone: no scan, no scan step, and
    // a budget of exactly the discovery's steps answers exact. Counting the sink as a
    // possible loss (J = I) would scan it, charge a step and, at that budget, answer bounds
    std::vector<std::string> records { "GCA", "TAC" };
    for (bool batch : { false, true }) {
        auto graph = build(3, records, DeBruijnGraph::BASIC, batch);
        const auto leaf = leaf_edges(base_dbg(*graph), "A", 'C');
        ASSERT_EQ(1u, leaf.candidates);
        ASSERT_EQ(1u, leaf.invalid);
        ASSERT_EQ(0u, leaf.invalid_non_sentinel);
        ASSERT_TRUE(invalid_edges(base_dbg(*graph)).count("CA$"));

        Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
        Request request = make_request(Scope::SUFFIX, Strands::FORWARD);
        Result r = count_of(*graph, ac, request);
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
        EXPECT_EQ(1u, r.contexts->total.value);
        EXPECT_EQ(0u, r.work.mask_scans);
        EXPECT_EQ(r.work.ranges_visited, r.work.steps);

        Result tight = count_of(*graph, ac, request, r.work.ranges_visited);
        EXPECT_FALSE(tight.stop);
        EXPECT_EQ(Relation::EXACT, tight.contexts->total.relation);
        EXPECT_EQ(1u, tight.contexts->total.value);
        check_against_oracles(*graph, ac, request, &records);
        check_against_oracles(*graph, ac, make_request(Scope::SUFFIX), &records);
    }
}

TEST(PatternSearch, DNA5Flank) {
    // §4.1: the flank admits every symbol of the build's alphabet. On a DNA5 build ACNTA is a
    // k-mer at k = 5 and contains AC at offset 0; on a DNA4 build the record splits into the
    // islands AC and TA, both shorter than k, so nothing is indexed and nothing is claimed.
    // The DNA5 branch is not in CI (no DNA5 build there; run by hand)
    std::vector<std::string> records { "ACNTA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    Result r = count_of(*graph, ac, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
#if _DNA5_GRAPH
    EXPECT_EQ(1u, r.contexts->total.value);
    EXPECT_EQ(1u, r.contexts->by_offset.at(0).value);
    check_against_oracles(*graph, ac, make_request(), &records);

    // N in flanks, in every graph mode, against both oracles (the record-scan oracle keeps
    // the k-mers containing N on DNA5); a pattern position never matches a record's N. No
    // k-mer here equals its reverse complement (an odd k-mer with N at its centre between
    // complementary flanks would: NCentredSelfComplementOddK)
    std::vector<std::string> flanked { "ACNTA", "GGACNNTTACCA", "TTNACGTNAA", "CANNNGTAC" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        for (bool batch : { false, true }) {
            auto g = build(5, flanked, mode, batch);
            for (std::string p : { "AC", "GT", "TA", "N", "NN", "ACG", "TTAC", "NNNNN", "ACNTA" }) {
                SCOPED_TRACE(p + " mode " + std::to_string(mode));
                // ACNTA is not an IUPAC pattern with a record N: its N admits A, C, G, T
                Pattern pattern = Pattern::parse(PatternKind::IUPAC, p);
                for (Scope scope : { Scope::SUFFIX, Scope::ANY_OFFSET }) {
                    if (mode == DeBruijnGraph::PRIMARY && scope == Scope::SUFFIX)
                        continue;
                    check_against_oracles(*g, pattern, make_request(scope), &flanked, mode);
                }
            }
        }
    }
#else
    EXPECT_EQ(0u, r.contexts->total.value);
    check_against_oracles(*graph, ac, make_request(), &records);
#endif
}

#if _DNA5_GRAPH
TEST(PatternSearch, NCentredSelfComplementOddK) {
    // On DNA5 N complements to N, so at odd k = 5 the k-mer ACNGT is its own reverse
    // complement, which no DNA4 k-mer of odd length can be. On BASIC and CANONICAL one k-mer
    // is one context per (offset, orientation) (§3, SPEC §7.1): AC at offset 0 (forward) and
    // GT, rc(AC), at offset 3 (reverse), once each. On a wrapped PRIMARY graph this is the
    // KNOWN DNA5 LIMITATION the engine header, SPEC §8.2 and DESIGN §15 state: CanonicalDBG
    // serves the stored ACNGT at two wrapper ids (y and y + offset, the same spelling), the
    // engine unites palindromes for even k only, so each orientation counts and releases it
    // twice: 2 + 2, the same two spellings at four node ids. The expectation was taken from a
    // run on a DNA5 build.
    std::vector<std::string> records { "ACNGT" };
    Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        SCOPED_TRACE(mode);
        const uint64_t per_orientation = mode == DeBruijnGraph::PRIMARY ? 2 : 1;
        auto graph = build(5, records, mode);
        Result r = count_of(*graph, ac, make_request());
        ASSERT_TRUE(r.contexts);
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
        EXPECT_EQ(per_orientation, r.contexts->by_orientation.at(Orientation::FORWARD).value);
        EXPECT_EQ(per_orientation, r.contexts->by_orientation.at(Orientation::REVERSE).value);
        EXPECT_EQ(2 * per_orientation, r.contexts->total.value);
        // the distinct (k-mer, offset, orientation) among the released contexts: two in
        // every mode; their node ids: two, four on PRIMARY (the double count)
        std::set<std::string> kmers;
        std::set<std::pair<uint64_t, std::string>> nodes;
        for (const Ctx &c : contexts_of(*graph, ac, make_request())) {
            const std::string key = graph->get_node_sequence(c.node) + "@"
                    + std::to_string(c.offset) + orientation_key(c.orientation);
            kmers.insert(key);
            nodes.emplace(c.node, key);
        }
        EXPECT_EQ(2u, kmers.size());
        EXPECT_EQ(2 * per_orientation, nodes.size());
        check_against_oracles(*graph, ac, make_request(), &records, mode);
    }
}
#endif

TEST(PatternSearch, IslandShorterThanK) {
    // k = 5, DNA4: AAAAANACNCCCCC keeps AAAAA and CCCCC; AC sits in a two-base island. On
    // DNA5 the N k-mers are kept, and AC lies in four of them: AANAC at 3, ANACN at 2, NACNC
    // at 1, ACNCC at 0 (a branch not in CI: no DNA5 build there)
    std::vector<std::string> records { "AAAAANACNCCCCC" };
    for (bool batch : { false, true }) {
        auto graph = build(5, records, DeBruijnGraph::BASIC, batch);
        Pattern ac = Pattern::parse(PatternKind::DNA, "AC");
        Result r = count_of(*graph, ac, make_request());
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
#if _DNA5_GRAPH
        EXPECT_EQ(4u, r.contexts->total.value);
        for (uint32_t p : { 0, 1, 2, 3 }) {
            EXPECT_EQ(1u, r.contexts->by_offset.at(p).value) << p;
        }
#else
        EXPECT_EQ(0u, r.contexts->total.value);
#endif
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
    // the wrapped PRIMARY graph of ACGA exposes ACGA and TCGT, one stored and one virtual
    // (which one is the builder's choice: TCGT is stored in this build). The last
    // three bases of the virtual k-mer form a pattern that is a suffix of no stored k-mer,
    // nor is its reverse complement: it lies at the prefix of the stored k-mer, read on the
    // other strand. A stored-suffix lookup cannot find it, so suffix is refused on a wrapped
    // PRIMARY graph and any_offset finds it, through the search of its reverse complement
    // mapped to the virtual node (§4.1)
    std::vector<std::string> records { "ACGA" };
    auto graph = build(4, records, DeBruijnGraph::PRIMARY);
    const DBGSuccinct &stored_graph = base_dbg(*graph);
    const bool acga_stored = stored_graph.kmer_to_node("ACGA") != DeBruijnGraph::npos;
    const std::string stored = acga_stored ? "ACGA" : "TCGT";
    const std::string virtual_kmer = acga_stored ? "TCGT" : "ACGA";
    ASSERT_NE(DeBruijnGraph::npos, stored_graph.kmer_to_node(stored));
    ASSERT_EQ(DeBruijnGraph::npos, stored_graph.kmer_to_node(virtual_kmer));
    const node_index virtual_node = graph->kmer_to_node(virtual_kmer);
    ASSERT_NE(DeBruijnGraph::npos, virtual_node);
    ASSERT_LT(stored_graph.max_index(), virtual_node);

    // the pattern sitting only at the virtual k-mer's suffix (CGA in this build)
    const std::string tail = virtual_kmer.substr(1);
    Pattern q = Pattern::parse(PatternKind::DNA, tail);
    for (std::string s : { tail, rev_comp(tail) }) {
        EXPECT_NE(s, stored.substr(1)) << "a stored suffix";
    }
    Result suffix = count_of(*graph, q, make_request(Scope::SUFFIX));
    ASSERT_TRUE(suffix.refusal);
    EXPECT_EQ("scope_unsupported", suffix.refusal->code);
    EXPECT_FALSE(suffix.contexts);
    EXPECT_EQ(0u, suffix.work.steps);

    auto contexts = contexts_of(*graph, q, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    ASSERT_EQ(1u, contexts.size());
    EXPECT_EQ(virtual_node, contexts[0].node);
    EXPECT_EQ(virtual_kmer, graph->get_node_sequence(contexts[0].node));
    EXPECT_EQ(1u, contexts[0].offset);  // k - L: the suffix
    Result any = count_of(*graph, q, make_request());
    EXPECT_EQ(Relation::EXACT, any.contexts->total.relation);
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::FORWARD).value);
    EXPECT_EQ(1u, any.contexts->suffix.value);
    check_against_oracles(*graph, q, make_request(), &records, DeBruijnGraph::PRIMARY);

    // and CGT, the suffix of TCGT, for contrast: refused in suffix scope all the same
    Pattern cgt = Pattern::parse(PatternKind::DNA, "CGT");
    suffix = count_of(*graph, cgt, make_request(Scope::SUFFIX));
    ASSERT_TRUE(suffix.refusal);
    EXPECT_EQ("scope_unsupported", suffix.refusal->code);
    EXPECT_FALSE(suffix.contexts);

    any = count_of(*graph, cgt, make_request());
    EXPECT_EQ(Relation::EXACT, any.contexts->total.relation);
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::FORWARD).value);   // TCGT
    EXPECT_EQ(1u, any.contexts->by_orientation.at(Orientation::REVERSE).value);   // ACGA
    EXPECT_EQ(std::vector<std::string>{ kNoteStrandUnknown }, any.notes);

    contexts = contexts_of(*graph, cgt, make_request(Scope::ANY_OFFSET, Strands::FORWARD));
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


// ---------------------------------------------------------------- stops and relations

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
    EXPECT_DOUBLE_EQ(12.0, *refused.min_anchor_information_bits);

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
            EXPECT_TRUE(r.min_anchor_information_bits);
            EXPECT_EQ(kNotePathsLater, r.notes.back());
        }
    }
}


// ---------------------------------------------------------------- stops in every graph mode

struct HaltCase {
    size_t k;
    DeBruijnGraph::Mode mode;
    std::vector<std::string> records;
    std::vector<std::string> patterns;
};

const std::vector<HaltCase>& halt_cases() {
    static const std::vector<std::string> kLong {
        "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA"
    };
    static const std::vector<std::string> kShort { "ACGTTGCAAGGCTTAC", "TTTACGGATC" };
    static const std::vector<HaltCase> cases {
        // even k, wrapped PRIMARY, palindromic 4-mers (AATT, CGCG, GCGC, AGCT, ...): the
        // union of a completed and a stopped probe may share palindromes, so its lower bound
        // is the larger part, not the sum (T, strands forward, max_steps 28, offset 2 is
        // at_least 7 of 9; the sum of the parts would claim 10)
        { 4, DeBruijnGraph::PRIMARY, { "ACGAAATTATG", "AGCTGTCTCGCGCGC" },
          { "T", "A", "G", "C", "N", "AT", "CG", "TA", "W", "ACGT", "AATT", "NNNN", "CGCGC" } },
        { 6, DeBruijnGraph::PRIMARY, { "TTACGCGTAA", "GAATATTCCG", "ACGCGTACGT", "CCGGAATTCC" },
          { "G", "N", "CG", "AT", "GCG", "AATT", "ACGCGT", "TTACGC", "ACGCGTA" } },
        // L = 1 (offset k - 1 of a stopped search is open, not exact)
        { 6, DeBruijnGraph::PRIMARY, { "TTA", "NGAGGTCGTGATCTCTAGCGCNTCCGGG", "CCTAGCCGCTCAAG" },
          { "G", "C", "T", "CG", "GTCGTG" } },
        // odd k, wrapped PRIMARY: no palindromic k-mer, the parts are disjoint
        { 5, DeBruijnGraph::PRIMARY, kShort, { "A", "G", "N", "AC", "GTT", "ACGTT", "ACGTTG" } },
        { 3, DeBruijnGraph::CANONICAL, kShort, { "A", "T", "N", "AC", "ACG", "GCAA" } },
        { 5, DeBruijnGraph::CANONICAL, kShort, { "C", "N", "AC", "GTT", "TTTAC", "TTTACG" } },
        { 4, DeBruijnGraph::BASIC, kLong, { "A", "T", "N", "NA", "GAT", "GATC", "R", "GATCG" } },
        { 6, DeBruijnGraph::BASIC, kLong, { "C", "G", "AC", "GATCGA", "GATCGAT" } },
    };
    return cases;
}

TEST(PatternSearch, HaltsAgainstOraclesInEveryMode) {
    // A step stop at every budget below the complete run's, in every graph mode (BASIC,
    // native CANONICAL, wrapped PRIMARY at odd k and at even k with palindromic k-mers), both
    // builders, both scopes where served, every strand choice, patterns of length 1, 2, k and
    // k + 1: every count states a true relation to the graph-walk oracle
    // (check_halts_against_oracles), not only on BASIC graphs with L >= 2
    uint64_t runs = 0;
    for (const HaltCase &c : halt_cases()) {
        for (bool batch : { false, true }) {
            auto graph = build(c.k, c.records, c.mode, batch);
            for (const std::string &p : c.patterns) {
                Pattern pattern = Pattern::parse(PatternKind::IUPAC, p);
                for (Scope scope : { Scope::ANY_OFFSET, Scope::SUFFIX }) {
                    if (scope == Scope::SUFFIX
                            && (c.mode == DeBruijnGraph::PRIMARY || p.size() > c.k))
                        continue;
                    // every strand choice for L = 1, both strands otherwise
                    for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
                        if (strands != Strands::BOTH && p.size() > 1)
                            continue;
                        SCOPED_TRACE("k " + std::to_string(c.k) + " mode "
                                     + std::to_string(c.mode) + " batch "
                                     + std::to_string(batch) + " " + p + " "
                                     + to_string(scope) + " " + to_string(strands));
                        runs += check_halts_against_oracles(*graph, pattern,
                                                            make_request(scope, strands));
                    }
                }
            }
        }
    }
    std::cerr << "stopped runs checked: " << runs << std::endl;
    EXPECT_LT(5'000u, runs);
}

TEST(PatternSearch, ThresholdsAgainstOraclesInEveryMode) {
    // stop_at_threshold and the release thresholds at every cap up to one past the truth, in
    // every graph mode, the even-k wrapped PRIMARY union with palindromic k-mers included (its
    // running lower bound is the larger part, never the sum)
    for (const HaltCase &c : halt_cases()) {
        auto graph = build(c.k, c.records, c.mode);
        for (const std::string &p : c.patterns) {
            Pattern pattern = Pattern::parse(PatternKind::IUPAC, p);
            for (Strands strands : { Strands::BOTH, Strands::FORWARD }) {
                SCOPED_TRACE("k " + std::to_string(c.k) + " mode " + std::to_string(c.mode)
                             + " " + p + " " + to_string(strands));
                check_thresholds_against_oracles(*graph, pattern,
                                                 make_request(Scope::ANY_OFFSET, strands));
            }
        }
    }
}

TEST(PatternSearch, RandomGraphsHaltsAgainstOracles) {
    // Step stops on random graphs of every mode (the mode rotates with the seed), three
    // patterns each: one of length 1, one of length k, one of a random length up to k + 1;
    // sampled sweeps (check_halts_against_oracles with max_sweep 300)
    const char *iupac = "ACGTRYSWKMBDHVN";
    uint64_t runs = 0;
    for (uint32_t seed = 1; seed <= 24; ++seed) {
        std::mt19937 rng(7919 * seed + 13);
        size_t k = 3 + rng() % 6;
        auto mode = static_cast<DeBruijnGraph::Mode>(seed % 3);
        bool batch = rng() % 2;
        std::vector<std::string> records(1 + rng() % 4);
        for (std::string &record : records) {
            size_t length = k + rng() % 25;
            for (size_t i = 0; i < length; ++i) {
                record.push_back(rng() % 25 ? "ACGT"[rng() % 4] : 'N');
            }
        }
        auto graph = build(k, records, mode, batch);
        for (size_t length : { size_t(1), k, 1 + rng() % (k + 1) }) {
            std::string text;
            const std::string &record = records[rng() % records.size()];
            if (record.size() >= length) {
                text = record.substr(rng() % (record.size() - length + 1), length);
            }
            for (char &c : text) {
                if (c == 'N' || rng() % 4 == 0)
                    c = iupac[rng() % 15];
            }
            while (text.size() < length) {
                text.push_back(iupac[rng() % 15]);
            }
            Scope scope = mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0
                ? Scope::SUFFIX
                : Scope::ANY_OFFSET;
            auto strands = static_cast<Strands>(rng() % 3);
            SCOPED_TRACE("seed " + std::to_string(seed) + " k " + std::to_string(k) + " mode "
                         + std::to_string(mode) + " batch " + std::to_string(batch) + " "
                         + text + " " + to_string(scope) + " " + to_string(strands));
            runs += check_halts_against_oracles(*graph, Pattern::parse(PatternKind::IUPAC, text),
                                                make_request(scope, strands), 300);
        }
    }
    std::cerr << "random stopped runs checked: " << runs << std::endl;
    EXPECT_LT(1'000u, runs);
}

TEST(PatternSearch, Kmer31BothStrandsMultiLocus) {
    // The engine's half (the integration suite's patterns are test_pattern.py's): at k = 31
    // on a BASIC graph, an exact pattern with several loci on each strand, so that the -
    // strand, several k-mers per (strand, offset) and the answer order across strands are
    // compared with the oracles, and an IUPAC pattern whose reverse complement matches the
    // same instances, so that one (node, offset) carries both orientations
    std::mt19937 rng(31);
    auto random_bases = [&](size_t n) {
        std::string s;
        for (size_t i = 0; i < n; ++i) {
            s.push_back("ACGT"[rng() % 4]);
        }
        return s;
    };
    const std::string p14 = "AAGGCCATTTCCGG";
    const std::string site = "GAATTC";
    std::vector<std::string> records;
    for (int r = 0; r < 4; ++r) {
        std::string record = random_bases(60);
        // two loci of P and one of rc(P) per record, each with its own flanks
        record += p14 + random_bases(45) + rev_comp(p14) + random_bases(40) + p14;
        // EcoRI sites inside RNNNNGAATTCNNNNY instances
        record += random_bases(30) + "A" + random_bases(4) + site + random_bases(4) + "C"
                + random_bases(50);
        records.push_back(record);
    }
    auto graph = build(31, records, DeBruijnGraph::BASIC);

    Pattern exact = Pattern::parse(PatternKind::DNA, p14);
    check_against_oracles(*graph, exact, make_request(), &records);
    check_against_oracles(*graph, exact, make_request(Scope::SUFFIX), &records);
    check_against_oracles(*graph, exact, make_request(Scope::ANY_OFFSET, Strands::REVERSE),
                          &records);
    // the inputs are not trivial: at some offset, two distinct k-mers on each strand
    const std::vector<Ctx> contexts = walk_oracle(*graph, exact, make_request());
    for (Orientation o : { Orientation::FORWARD, Orientation::REVERSE }) {
        std::map<uint32_t, std::set<node_index>> nodes;
        for (const Ctx &c : contexts) {
            if (c.orientation == o)
                nodes[c.offset].insert(c.node);
        }
        size_t most = 0;
        for (const auto &[p, at] : nodes) {
            most = std::max(most, at.size());
        }
        EXPECT_LE(2u, most) << orientation_key(o);
    }

    // rc(RNNNNGAATTCNNNNN) = NNNNNGAATTCNNNNY: an instance with R first and Y last matches
    // both, at the same (node, offset)
    Pattern both = Pattern::parse(PatternKind::IUPAC, "RNNNNGAATTCNNNNN");
    ASSERT_FALSE(both.is_palindromic());
    check_against_oracles(*graph, both, make_request(), &records);
    auto released = contexts_of(*graph, both, make_request());
    size_t pairs = 0;
    for (size_t i = 1; i < released.size(); ++i) {
        if (released[i].node == released[i - 1].node
                && released[i].offset == released[i - 1].offset) {
            // one (node, offset), two orientations: forward first (§5.5)
            EXPECT_EQ(Orientation::FORWARD, released[i - 1].orientation);
            EXPECT_EQ(Orientation::REVERSE, released[i].orientation);
            ++pairs;
        }
    }
    EXPECT_LE(4u, pairs);
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


// ---------------------------------------------------------------- extension

/*
 * Phase 2 (§4.2): a pattern longer than k with Request::extend_paths. Two oracles, sharing
 * no code with the engine's DFS:
 *  - the graph-walk oracle extends strings over the served graph's k-mers (graph_kmers, a
 *    map from each k-mer to its node id): every walk of L - k + 1 k-mers spelling an
 *    instance of an oriented pattern, never through call_outgoing_kmers; it also counts the
 *    partial walks of k + 1 .. L bases the DFS must enter (candidates_examined);
 *  - the record-scan oracle lists the windows of L bases of the records' islands (and of
 *    their reverse complements unless BASIC) that instantiate an oriented pattern and whose
 *    k-mers were all retained (the records' k-mers less those a test pruned).
 * Every record occurrence is a graph path; the converse fails exactly where two records
 * join in the graph (a path is a graph context, not an occurrence: §3, §4.3).
 */

Request path_request(Scope scope = Scope::ANY_OFFSET, Strands strands = Strands::BOTH) {
    Request request = make_request(scope, strands);
    request.extend_paths = true;
    request.max_paths = 1'000'000'000;
    return request;
}

struct PathCtx {
    node_index anchor;
    Orientation orientation;
    std::string sequence;
    std::vector<node_index> path;

    // the answer order (§5.5): anchor node, orientation, then the DFS in symbol order
    bool operator<(const PathCtx &o) const {
        return std::tie(anchor, orientation, sequence)
                < std::tie(o.anchor, o.orientation, o.sequence);
    }
    bool operator==(const PathCtx &o) const {
        return std::tie(anchor, orientation, sequence, path)
                == std::tie(o.anchor, o.orientation, o.sequence, o.path);
    }
};

std::ostream& operator<<(std::ostream &out, const PathCtx &p) {
    out << "(" << p.anchor << ", " << orientation_key(p.orientation) << ", " << p.sequence
        << ", [";
    for (node_index n : p.path) {
        out << " " << n;
    }
    return out << " ])";
}

struct PathOracle {
    // in answer order
    std::vector<PathCtx> paths;
    // the anchors per orientation
    std::map<Orientation, uint64_t> anchors;
    // the partial walks of k + 1 .. L bases instantiating the oriented pattern's prefix
    uint64_t candidates = 0;
    // the walks of k .. L - 1 bases (anchors included) that two or more such walks extend by
    // one base: the DFS's branchings (Work::extension_branches)
    uint64_t branchings = 0;
};

PathOracle path_oracle(const DeBruijnGraph &graph, const Pattern &pattern,
                       const Request &request) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    std::map<std::string, node_index> node_of;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        node_of.emplace(kmer, node);
    }
    PathOracle result;
    for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
        const Orientation o = orientation;
        const Bases &oriented = q;
        std::function<void(std::string&, std::vector<node_index>&)> grow
                = [&](std::string &s, std::vector<node_index> &p) {
            if (s.size() == L) {
                result.paths.push_back(PathCtx { p.front(), o, s, p });
                return;
            }
            // a pattern position admits A, C, G, T only, never N
            uint64_t children = 0;
            for (char b : std::string("ACGT")) {
                if (oriented[s.size()].find(b) != std::string::npos
                        && node_of.count(s.substr(s.size() - k + 1) + b))
                    ++children;
            }
            result.branchings += children > 1;
            for (char b : std::string("ACGT")) {
                if (oriented[s.size()].find(b) == std::string::npos)
                    continue;
                auto it = node_of.find(s.substr(s.size() - k + 1) + b);
                if (it == node_of.end())
                    continue;
                ++result.candidates;
                s.push_back(b);
                p.push_back(it->second);
                grow(s, p);
                s.pop_back();
                p.pop_back();
            }
        };
        for (const auto &[kmer, node] : node_of) {
            if (!matches(window(q, k), kmer))
                continue;
            ++result.anchors[o];
            std::string s = kmer;
            std::vector<node_index> p { node };
            grow(s, p);
        }
    }
    std::sort(result.paths.begin(), result.paths.end());
    return result;
}

typedef std::set<std::pair<std::string, Orientation>> SpelledPaths;

SpelledPaths record_path_oracle(const std::vector<std::string> &records, size_t k,
                                DeBruijnGraph::Mode mode, const Pattern &pattern,
                                const Request &request,
                                const std::set<std::string> &pruned = {}) {
    std::vector<std::string> strands = records;
    if (mode != DeBruijnGraph::BASIC) {
        for (const std::string &record : records) {
            strands.push_back(rev_comp(record));
        }
    }
    // the k-mers of the islands (indexed symbols only: N splits a record on DNA4); an
    // occurrence's N, on DNA5, is never at a pattern position (matches())
    std::set<std::string> retained;
    for (const std::string &s : strands) {
        for (size_t i = 0; i + k <= s.size(); ++i) {
            std::string kmer = s.substr(i, k);
            if (indexed(kmer) && !pruned.count(kmer))
                retained.insert(kmer);
        }
    }
    const size_t L = pattern.length();
    SpelledPaths result;
    for (const std::string &s : strands) {
        for (size_t i = 0; i + L <= s.size(); ++i) {
            std::string occurrence = s.substr(i, L);
            if (!indexed(occurrence))
                continue;
            bool all_retained = true;
            for (size_t j = 0; j + k <= L; ++j) {
                all_retained &= retained.count(occurrence.substr(j, k)) > 0;
            }
            if (!all_retained)
                continue;
            for (const auto &[orientation, q] : orientations(pattern, request.strands)) {
                if (matches(q, occurrence))
                    result.emplace(occurrence, orientation);
            }
        }
    }
    return result;
}

std::vector<PathCtx> run_paths(const PatternSearch &engine, const StoredNodes &stored,
                               const Pattern &pattern, const Request &request, Result *result,
                               uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    std::vector<PathCtx> paths;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        EXPECT_EQ(0u, c.offset);
        EXPECT_FALSE(c.path.empty());
        // braces: the gtest macro expands to an if/else (GCC's -Wdangling-else)
        if (c.path.size()) {
            EXPECT_EQ(c.node, c.path.front());
        }
        // the anchor's stored node, and the one the route maps every node of the path to
        // (PatternSearch::base_node), against the stored k-mers spelled
        EXPECT_EQ(stored.of(c.node), c.base_node);
        for (node_index n : c.path) {
            EXPECT_EQ(stored.of(n), engine.base_node(n)) << n;
        }
        paths.push_back(PathCtx { c.node, c.orientation, c.sequence, c.path });
    });
    return paths;
}

Result count_paths(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                   uint64_t max_steps = kManySteps) {
    Budget budget = unbounded_budget(max_steps);
    return PatternSearch(graph).count(pattern, request, budget);
}

std::vector<PathCtx> paths_of(const DeBruijnGraph &graph, const Pattern &pattern,
                              const Request &request) {
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result result;
    return run_paths(PatternSearch(graph), StoredNodes(graph), pattern, all, &result);
}

// marks |kmer| pruned in the graph's valid-edge mask, as graph cleaning would (§3)
void prune(const DeBruijnGraph &graph, const std::string &kmer) {
    auto &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(graph));
    node_index node = dbg_succ.kmer_to_node(kmer);
    ASSERT_NE(DeBruijnGraph::npos, node) << kmer;
    auto *mask = dynamic_cast<bit_vector_dyn*>(const_cast<bit_vector*>(dbg_succ.get_mask()));
    ASSERT_TRUE(mask);
    mask->set(node, false);
    ASSERT_FALSE(dbg_succ.in_graph(node));
}

/**
 * Every check of one (graph, long pattern, request with extend_paths) against the oracles:
 * the anchor and path counts in every unit and relation, candidates_examined, the steps
 * (discovery as without the extension, plus one per outgoing edge examined), the released
 * paths and their order, each path's nodes and sequence, ALL_OR_COUNT at the threshold and
 * PARTIAL's prefixes. |records| enables the record-scan comparison (equality when
 * |records_exact|, else inclusion).
 */
void check_paths_against_oracles(const DeBruijnGraph &graph, const Pattern &pattern,
                                 const Request &request,
                                 const std::vector<std::string> *records = nullptr,
                                 DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC,
                                 bool records_exact = false,
                                 const std::set<std::string> &pruned = {},
                                 size_t *num_expected = nullptr) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    ASSERT_LT(k, L);
    ASSERT_TRUE(request.extend_paths);
    PatternSearch engine(graph);
    const StoredNodes stored(graph);

    const PathOracle oracle = path_oracle(graph, pattern, request);
    const std::vector<PathCtx> &expected = oracle.paths;
    if (num_expected)
        *num_expected = expected.size();
    uint64_t num_anchors = 0;
    for (const auto &[o, n] : oracle.anchors) {
        num_anchors += n;
    }

    if (records) {
        SpelledPaths in_records = record_path_oracle(*records, k, mode, pattern, request, pruned);
        SpelledPaths in_graph;
        for (const PathCtx &p : expected) {
            in_graph.emplace(p.sequence, p.orientation);
        }
        // every occurrence whose k-mers were all retained is a graph path (§5.1)
        for (const auto &occurrence : in_records) {
            EXPECT_TRUE(in_graph.count(occurrence))
                << occurrence.first << " " << orientation_key(occurrence.second);
        }
        if (records_exact) {
            EXPECT_EQ(in_records, in_graph) << pattern.text();
        }
    }

    // count(): the anchors as without the extension, the paths exact
    Budget budget = unbounded_budget();
    Result counted = engine.count(pattern, request, budget);
    ASSERT_FALSE(counted.refusal) << counted.refusal->message;
    ASSERT_TRUE(counted.anchors);
    EXPECT_FALSE(counted.contexts);
    EXPECT_EQ(Scope::LONG, counted.scope);
    EXPECT_FALSE(counted.stop);
    EXPECT_FALSE(counted.time_limited);
    EXPECT_FALSE(counted.extraction);
    const AnchorCounts &a = *counted.anchors;
    EXPECT_EQ(Relation::EXACT, a.total.relation);
    EXPECT_EQ(num_anchors, a.total.value) << pattern.text();
    EXPECT_EQ(num_anchors ? Extension::COMPLETED : Extension::NO_ANCHORS, a.extension);
    EXPECT_EQ(Unit::PATHS, a.paths.unit);
    EXPECT_EQ(Relation::EXACT, a.paths.relation);
    EXPECT_EQ(expected.size(), a.paths.value) << pattern.text();
    EXPECT_EQ(oracle.candidates, a.candidates_examined) << pattern.text();
    EXPECT_LE(expected.size(), a.candidates_examined);
    // the extension's counters beside its steps: every anchor spelled and extended, and the
    // walks the DFS branched at
    EXPECT_EQ(num_anchors, counted.work.extension_anchors) << pattern.text();
    EXPECT_EQ(oracle.branchings, counted.work.extension_branches) << pattern.text();
    std::map<Orientation, uint64_t> by_orientation;
    for (const PathCtx &p : expected) {
        ++by_orientation[p.orientation];
    }
    ASSERT_EQ(orientations(pattern, request.strands).size(), a.paths_by_orientation.size());
    for (const auto &[o, count] : a.paths_by_orientation) {
        EXPECT_EQ(Relation::EXACT, count.relation);
        EXPECT_EQ(Unit::PATHS, count.unit);
        EXPECT_EQ(by_orientation[o], count.value) << pattern.text() << " " << orientation_key(o);
    }
    EXPECT_TRUE(std::find(counted.notes.begin(), counted.notes.end(), kNotePathsLater)
                    == counted.notes.end());
    if (!num_anchors) {
        EXPECT_EQ(0u, counted.work.extension_edges);
    }

    // the discovery is that of the anchors alone; the extension adds one step per edge
    Request anchors_only = request;
    anchors_only.extend_paths = false;
    Result plain = count_of(graph, pattern, anchors_only);
    EXPECT_EQ(plain.anchors->total.value, a.total.value);
    EXPECT_EQ(Extension::NOT_REQUESTED, plain.anchors->extension);
    EXPECT_EQ(plain.work.ranges_visited, counted.work.ranges_visited);
    EXPECT_EQ(plain.work.steps + counted.work.extension_edges, counted.work.steps);

    // enumerate(), ALL_OR_COUNT: every path in answer order
    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result enumerated;
    std::vector<PathCtx> released = run_paths(engine, stored, pattern, all, &enumerated);
    EXPECT_EQ(counted.work.steps, enumerated.work.steps);
    EXPECT_EQ(counted.work.extension_edges, enumerated.work.extension_edges);
    EXPECT_EQ(counted.work.extension_anchors, enumerated.work.extension_anchors);
    EXPECT_EQ(counted.work.extension_branches, enumerated.work.extension_branches);
    EXPECT_EQ(a.candidates_examined, enumerated.anchors->candidates_examined);
    ASSERT_TRUE(enumerated.extraction);
    EXPECT_TRUE(enumerated.extraction->complete);
    EXPECT_FALSE(enumerated.extraction->withheld);
    EXPECT_FALSE(enumerated.extraction->cut);
    EXPECT_EQ(expected.size(), enumerated.extraction->returned);
    ASSERT_EQ(expected, released) << pattern.text();

    for (const PathCtx &p : released) {
        // n = L - k + 1 retained k-mers, each the next window of the spelled sequence, which
        // instantiates the oriented pattern
        ASSERT_EQ(L, p.sequence.size());
        ASSERT_EQ(L - k + 1, p.path.size());
        for (size_t i = 0; i < p.path.size(); ++i) {
            EXPECT_EQ(p.sequence.substr(i, k), graph.get_node_sequence(p.path[i]));
        }
        EXPECT_TRUE(matches(oriented(pattern, p.orientation), p.sequence)) << p;
    }

    // ALL_OR_COUNT: all or nothing at max_paths
    if (!expected.empty()) {
        Request tight = all;
        tight.max_paths = expected.size() - 1;
        Result withheld;
        EXPECT_TRUE(run_paths(engine, stored, pattern, tight, &withheld).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
        EXPECT_FALSE(withheld.extraction->complete);
        EXPECT_EQ(Relation::EXACT, withheld.anchors->paths.relation);
        EXPECT_EQ(expected.size(), withheld.anchors->paths.value);
    }

    // PARTIAL: the first max_paths in answer order, the cut stated
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()) }) {
        partial.max_paths = cap;
        Result cut;
        std::vector<PathCtx> prefix = run_paths(engine, stored, pattern, partial, &cut);
        ASSERT_EQ(std::min<uint64_t>(cap, expected.size()), prefix.size());
        EXPECT_TRUE(std::equal(prefix.begin(), prefix.end(), expected.begin()));
        EXPECT_EQ(Relation::EXACT, cut.anchors->paths.relation);
        if (cap >= expected.size()) {
            EXPECT_TRUE(cut.extraction->complete);
            EXPECT_FALSE(cut.extraction->cut);
        } else {
            EXPECT_FALSE(cut.extraction->complete);
            EXPECT_EQ(StopReason::MAX_PATHS, cut.extraction->cut);
        }
    }
}


TEST(PatternSearch, ExtensionNotRequestedUnchanged) {
    // without extend_paths a long pattern is answered by its anchors, step for step
    std::vector<std::string> records { "TTACGTACCA", "GGACGTTT" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);
    Pattern acgta = Pattern::parse(PatternKind::DNA, "ACGTA");

    Result plain = count_of(*graph, acgta, make_request());
    EXPECT_EQ(Extension::NOT_REQUESTED, plain.anchors->extension);
    EXPECT_EQ(Relation::UNKNOWN, plain.anchors->paths.relation);
    EXPECT_TRUE(plain.anchors->paths_by_orientation.empty());
    EXPECT_EQ(0u, plain.anchors->candidates_examined);
    EXPECT_EQ(0u, plain.work.extension_edges);
    EXPECT_EQ(0.0, plain.extension_ms);
    EXPECT_EQ(kNotePathsLater, plain.notes.back());

    Result extended = count_paths(*graph, acgta, path_request());
    EXPECT_EQ(Extension::COMPLETED, extended.anchors->extension);
    EXPECT_EQ(Relation::EXACT, extended.anchors->paths.relation);
    EXPECT_TRUE(extended.notes.empty());
    EXPECT_LT(0u, extended.work.extension_edges);
    EXPECT_EQ(plain.work.steps + extended.work.extension_edges, extended.work.steps);

    // the anchors are never results once paths are: the combination is refused
    Request both = path_request();
    both.release_anchors = true;
    Budget budget = unbounded_budget();
    EXPECT_THROW(PatternSearch(*graph).enumerate(acgta, both, budget, [](const Context &) {}),
                 std::invalid_argument);
    // count() reads no release field
    EXPECT_NO_THROW(PatternSearch(*graph).count(acgta, both, budget));

    // L <= k: extend_paths changes nothing
    Pattern cg = Pattern::parse(PatternKind::DNA, "CG");
    Result short_plain = count_of(*graph, cg, make_request());
    Result short_ext = count_paths(*graph, cg, path_request());
    EXPECT_EQ(short_plain.work.steps, short_ext.work.steps);
    EXPECT_EQ(short_plain.contexts->total.value, short_ext.contexts->total.value);
    EXPECT_EQ(contexts_of(*graph, cg, make_request()).size(),
              contexts_of(*graph, cg, path_request()).size());
}

TEST(PatternSearch, ExtensionTwoAndThreeKmers) {
    // k = 4: ACGTA spans two k-mers (ACGT, CGTA), ACGTAC three
    std::vector<std::string> records { "TTACGTACCA", "GGACGTTT" };
    auto graph = build(4, records, DeBruijnGraph::BASIC);

    Pattern two = Pattern::parse(PatternKind::DNA, "ACGTA");
    Result r = count_paths(*graph, two, path_request());
    ASSERT_TRUE(r.anchors);
    // + : ACGT is the anchor; of its outgoing CGTA and CGTT only CGTA continues the pattern
    // - : TACGT, anchored at TACG
    EXPECT_EQ(2u, r.anchors->total.value);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(2u, r.anchors->paths.value);
    EXPECT_EQ(1u, r.anchors->paths_by_orientation.at(Orientation::FORWARD).value);
    EXPECT_EQ(1u, r.anchors->paths_by_orientation.at(Orientation::REVERSE).value);
    EXPECT_EQ(2u, r.anchors->candidates_examined);
    // the outgoing k-mers of ACGT (CGTA, CGTT) and of TACG (ACGT), one step each
    EXPECT_EQ(3u, r.work.extension_edges);
    EXPECT_EQ(r.work.extension_edges + count_of(*graph, two, make_request()).work.steps,
              r.work.steps);

    auto paths = paths_of(*graph, two, path_request(Scope::ANY_OFFSET, Strands::FORWARD));
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ("ACGTA", paths[0].sequence);
    ASSERT_EQ(2u, paths[0].path.size());
    EXPECT_EQ("ACGT", graph->get_node_sequence(paths[0].path[0]));
    EXPECT_EQ("CGTA", graph->get_node_sequence(paths[0].path[1]));
    check_paths_against_oracles(*graph, two, path_request(), &records, DeBruijnGraph::BASIC,
                                true);

    // three k-mers: + ACGTAC lies in the first record; - GTACGT is a graph path through
    // GTAC -> TACG -> ACGT that no record contains (a context, not an occurrence)
    Pattern three = Pattern::parse(PatternKind::DNA, "ACGTAC");
    Result t = count_paths(*graph, three, path_request());
    EXPECT_EQ(Relation::EXACT, t.anchors->paths.relation);
    EXPECT_EQ(1u, t.anchors->paths_by_orientation.at(Orientation::FORWARD).value);
    EXPECT_EQ(1u, t.anchors->paths_by_orientation.at(Orientation::REVERSE).value);
    auto reverse = paths_of(*graph, three, path_request(Scope::ANY_OFFSET, Strands::REVERSE));
    ASSERT_EQ(1u, reverse.size());
    EXPECT_EQ("GTACGT", reverse[0].sequence);
    EXPECT_EQ(Orientation::REVERSE, reverse[0].orientation);
    EXPECT_EQ(0u, record_path_oracle(records, 4, DeBruijnGraph::BASIC, three,
                                     path_request(Scope::ANY_OFFSET, Strands::REVERSE)).size());
    check_paths_against_oracles(*graph, three, path_request(), &records);

    for (auto mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto other = build(4, records, mode);
        for (const Pattern &p : { two, three }) {
            check_paths_against_oracles(*other, p, path_request(), &records, mode);
        }
    }
}

TEST(PatternSearch, ExtensionAnchorWithoutPath) {
    // k = 3, AAAC: the anchor AAA exists, AAC does not: one anchor, zero paths, exactly. An
    // anchor count says nothing about paths (§3): without the extension they are unknown
    std::vector<std::string> records { "AAAA", "GAAT" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    Pattern aaac = Pattern::parse(PatternKind::DNA, "AAAC");
    Request request = path_request(Scope::ANY_OFFSET, Strands::FORWARD);

    Result r = count_paths(*graph, aaac, request);
    EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
    EXPECT_EQ(1u, r.anchors->total.value);
    EXPECT_EQ(Extension::COMPLETED, r.anchors->extension);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(0u, r.anchors->paths.value);
    EXPECT_EQ(0u, r.anchors->candidates_examined);
    // AAA's outgoing k-mers AAA and AAT were examined, neither allowed
    EXPECT_EQ(2u, r.work.extension_edges);

    Request plain = request;
    plain.extend_paths = false;
    EXPECT_EQ(Relation::UNKNOWN, count_of(*graph, aaac, plain).anchors->paths.relation);

    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result e;
    EXPECT_TRUE(run_paths(PatternSearch(*graph), StoredNodes(*graph), aaac, all, &e).empty());
    EXPECT_TRUE(e.extraction->complete);
    EXPECT_EQ(0u, e.extraction->returned);
    check_paths_against_oracles(*graph, aaac, path_request(), &records, DeBruijnGraph::BASIC,
                                true);

    // and no anchor at all: nothing to extend, paths exact 0 by derivation
    Pattern cccca = Pattern::parse(PatternKind::DNA, "CCCCA");
    Result none = count_paths(*graph, cccca, path_request());
    EXPECT_EQ(Extension::NO_ANCHORS, none.anchors->extension);
    EXPECT_EQ(Relation::EXACT, none.anchors->paths.relation);
    EXPECT_EQ(0u, none.anchors->paths.value);
    EXPECT_EQ(0u, none.work.extension_edges);
    check_paths_against_oracles(*graph, cccca, path_request(), &records);
}

TEST(PatternSearch, ExtensionPrunedMiddleKmer) {
    // k = 3, ACGTA with CGT pruned: ACG and GTA kept, so every base is covered, but no path
    // spells ACGTA (§3 "Covered sequence": a long pattern needs every k-mer retained)
    std::vector<std::string> records { "ACGTA" };
    Pattern acgta = Pattern::parse(PatternKind::DNA, "ACGTA");
    Request request = path_request(Scope::ANY_OFFSET, Strands::FORWARD);

    auto intact = build(3, records, DeBruijnGraph::BASIC);
    Result before = count_paths(*intact, acgta, request);
    EXPECT_EQ(1u, before.anchors->paths.value);
    check_paths_against_oracles(*intact, acgta, path_request(), &records, DeBruijnGraph::BASIC,
                                true);

    auto graph = build(3, records, DeBruijnGraph::BASIC);
    prune(*graph, "CGT");
    Result r = count_paths(*graph, acgta, request);
    EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
    EXPECT_EQ(1u, r.anchors->total.value);  // ACG
    EXPECT_EQ(Extension::COMPLETED, r.anchors->extension);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(0u, r.anchors->paths.value);
    // the pruned k-mer is not even an outgoing edge of ACG
    EXPECT_EQ(0u, r.work.extension_edges);
    EXPECT_EQ(0u, record_path_oracle(records, 3, DeBruijnGraph::BASIC, acgta, request,
                                     { "CGT" }).size());
    check_paths_against_oracles(*graph, acgta, path_request(), &records, DeBruijnGraph::BASIC,
                                true, { "CGT" });
    // both flanking k-mers are still found as contexts of their own
    EXPECT_EQ(1u, count_of(*graph, Pattern::parse(PatternKind::DNA, "GTA"),
                           make_request(Scope::ANY_OFFSET, Strands::FORWARD))
                      .contexts->total.value);
}

TEST(PatternSearch, ExtensionCrossRecordPath) {
    // k = 3, records ACG and CGT: the graph path ACGT, which no record contains (§13: the
    // graph-walk oracle expects one context, the record scan none). The engine counts graph
    // contexts; telling them from record occurrences is retrieval's per-label support (§4.3)
    std::vector<std::string> records { "ACG", "CGT" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    Pattern acgt = Pattern::parse(PatternKind::DNA, "ACGT");
    Result r = count_paths(*graph, acgt, path_request());
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(1u, r.anchors->paths.value);  // palindromic: searched once
    EXPECT_TRUE(r.palindromic);
    EXPECT_EQ(0u, record_path_oracle(records, 3, DeBruijnGraph::BASIC, acgt,
                                     path_request()).size());
    check_paths_against_oracles(*graph, acgt, path_request(), &records);
}

TEST(PatternSearch, ExtensionReverseHitOnWrappedPrimary) {
    // k = 3, AACG on the wrapped PRIMARY graph of CGTT: the anchor AAC is the virtual
    // reverse complement of stored GTT (found through rc(AAC), §4.1), and the path goes on
    // along the wrapper to ACG, the virtual reverse complement of CGT
    std::vector<std::string> records { "CGTT" };
    auto graph = build(3, records, DeBruijnGraph::PRIMARY);
    Pattern aacg = Pattern::parse(PatternKind::DNA, "AACG");

    auto paths = paths_of(*graph, aacg, path_request(Scope::ANY_OFFSET, Strands::FORWARD));
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ("AACG", paths[0].sequence);
    ASSERT_EQ(2u, paths[0].path.size());
    EXPECT_EQ("AAC", graph->get_node_sequence(paths[0].path[0]));
    EXPECT_EQ("ACG", graph->get_node_sequence(paths[0].path[1]));

    Result r = count_paths(*graph, aacg, path_request());
    EXPECT_EQ(2u, r.anchors->total.value);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(1u, r.anchors->paths_by_orientation.at(Orientation::FORWARD).value);   // AACG
    EXPECT_EQ(1u, r.anchors->paths_by_orientation.at(Orientation::REVERSE).value);   // CGTT
    check_paths_against_oracles(*graph, aacg, path_request(), &records,
                                DeBruijnGraph::PRIMARY, true);

    // on the BASIC graph of the same record only the reverse orientation exists
    auto basic = build(3, records, DeBruijnGraph::BASIC);
    Result b = count_paths(*basic, aacg, path_request());
    EXPECT_EQ(0u, b.anchors->paths_by_orientation.at(Orientation::FORWARD).value);
    EXPECT_EQ(1u, b.anchors->paths_by_orientation.at(Orientation::REVERSE).value);
    check_paths_against_oracles(*basic, aacg, path_request(), &records, DeBruijnGraph::BASIC,
                                true);
}

TEST(PatternSearch, ExtensionPalindromicAnchorWindow) {
    // k = 4, ACGTA on the wrapped PRIMARY graph of ACGTA: the anchor window ACGT is a
    // palindromic k-mer found by both probes and united: one anchor, one path (§4.1)
    std::vector<std::string> records { "ACGTA" };
    auto graph = build(4, records, DeBruijnGraph::PRIMARY);
    Pattern acgta = Pattern::parse(PatternKind::DNA, "ACGTA");
    Request forward = path_request(Scope::ANY_OFFSET, Strands::FORWARD);
    Result r = count_paths(*graph, acgta, forward);
    EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
    EXPECT_EQ(1u, r.anchors->total.value);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(1u, r.anchors->paths.value);
    auto paths = paths_of(*graph, acgta, forward);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ("ACGTA", paths[0].sequence);
    check_paths_against_oracles(*graph, acgta, path_request(), &records,
                                DeBruijnGraph::PRIMARY, true);

    // an IUPAC window with one palindromic instance (ACGT, which continues with A) and one
    // not (ACGA, which does not)
    std::vector<std::string> mixed { "ACGTA", "ACGAC" };
    auto graph2 = build(4, mixed, DeBruijnGraph::PRIMARY);
    Pattern acgwa = Pattern::parse(PatternKind::IUPAC, "ACGWA");
    Result m = count_paths(*graph2, acgwa, forward);
    EXPECT_EQ(2u, m.anchors->total.value);
    EXPECT_EQ(1u, m.anchors->paths.value);
    check_paths_against_oracles(*graph2, acgwa, path_request(), &mixed,
                                DeBruijnGraph::PRIMARY, true);

    // a palindromic k-mer inside a path: TACG -> ACGT -> CGTA
    std::vector<std::string> inside { "TTACGTAA" };
    auto graph3 = build(4, inside, DeBruijnGraph::PRIMARY);
    for (std::string p : { "TACGTA", "TTACGTAA", "NACGTN", "ACGTAA" }) {
        check_paths_against_oracles(*graph3, Pattern::parse(PatternKind::IUPAC, p),
                                    path_request(), &inside, DeBruijnGraph::PRIMARY);
    }
    EXPECT_EQ(1u, count_paths(*graph3, Pattern::parse(PatternKind::DNA, "TACGTA"),
                              forward).anchors->paths.value);
}

TEST(PatternSearch, ExtensionIUPAC) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(5, records, mode);
        for (std::string p : { "ACGTTGC", "GATCGATC", "NNNNNNA", "GGGCCCAA", "TTTTTTTT",
                               "RNNNYNG", "GATCNATC", "ACGNTGCAAGG", "NNNNNNNN", "WSWSWS" }) {
            size_t n = 0;
            check_paths_against_oracles(*graph, Pattern::parse(PatternKind::IUPAC, p),
                                        path_request(), &records, mode, false, {}, &n);
            Result r = count_paths(*graph, Pattern::parse(PatternKind::IUPAC, p),
                                   path_request());
            EXPECT_EQ(n, r.anchors->paths.value) << p;
        }
    }
}

TEST(PatternSearch, ExtensionCandidatesExamined) {
    // k = 3, records ACGTT and ACGAA: from the anchor ACG the DFS examines CGA and CGT
    // (two edges) and enters CGT (one branch), then examines GTT and enters it (complete)
    std::vector<std::string> records { "ACGTT", "ACGAA" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    Request forward = path_request(Scope::ANY_OFFSET, Strands::FORWARD);

    Result exact = count_paths(*graph, Pattern::parse(PatternKind::DNA, "ACGTT"), forward);
    EXPECT_EQ(1u, exact.anchors->paths.value);
    EXPECT_EQ(2u, exact.anchors->candidates_examined);
    EXPECT_EQ(3u, exact.work.extension_edges);

    // ACGNN: both branches entered and completed, in symbol order
    Pattern acgnn = Pattern::parse(PatternKind::IUPAC, "ACGNN");
    Result iupac = count_paths(*graph, acgnn, forward);
    EXPECT_EQ(2u, iupac.anchors->paths.value);
    EXPECT_EQ(4u, iupac.anchors->candidates_examined);
    EXPECT_EQ(4u, iupac.work.extension_edges);
    auto paths = paths_of(*graph, acgnn, forward);
    ASSERT_EQ(2u, paths.size());
    EXPECT_EQ("ACGAA", paths[0].sequence);
    EXPECT_EQ("ACGTT", paths[1].sequence);
    check_paths_against_oracles(*graph, acgnn, path_request(), &records, DeBruijnGraph::BASIC,
                                true);
}

TEST(PatternSearch, ExtensionMaxAnchors) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NNNNNNA");
    Result full = count_paths(*graph, pattern, path_request());
    ASSERT_EQ(Relation::EXACT, full.anchors->total.relation);
    const uint64_t anchors = full.anchors->total.value;
    ASSERT_LT(4u, anchors);
    ASSERT_LT(0u, full.anchors->paths.value);

    // above max_anchors: the extension is not admitted, in every mode (§4.2)
    Request request = path_request();
    request.max_anchors = anchors - 1;
    Result r = count_paths(*graph, pattern, request);
    EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
    EXPECT_EQ(anchors, r.anchors->total.value);
    EXPECT_EQ(Extension::NOT_ADMITTED, r.anchors->extension);
    EXPECT_EQ(Relation::UNKNOWN, r.anchors->paths.relation);
    EXPECT_EQ(0u, r.anchors->candidates_examined);
    EXPECT_EQ(0u, r.work.extension_edges);
    EXPECT_FALSE(r.stop);
    for (const auto &[o, count] : r.anchors->paths_by_orientation) {
        EXPECT_EQ(Relation::UNKNOWN, count.relation);
    }
    for (Mode mode : { Mode::ALL_OR_COUNT, Mode::PARTIAL }) {
        request.mode = mode;
        Result e;
        EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, request, &e).empty());
        EXPECT_EQ(Withheld::ANCHORS_ABOVE_THRESHOLD, e.extraction->withheld);
        EXPECT_FALSE(e.extraction->complete);
        EXPECT_EQ(Extension::NOT_ADMITTED, e.anchors->extension);
    }

    // at max_anchors: admitted
    request.max_anchors = anchors;
    EXPECT_EQ(Extension::COMPLETED, count_paths(*graph, pattern, request).anchors->extension);

    // stop_at_threshold: the anchors stop above max_anchors; nothing is extended
    request.stop_at_threshold = true;
    request.max_anchors = 2;
    Result s = count_paths(*graph, pattern, request);
    ASSERT_TRUE(s.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, s.stop->phase);
    EXPECT_EQ(StopReason::MAX_ANCHORS, s.stop->reason);
    EXPECT_EQ(Relation::AT_LEAST, s.anchors->total.relation);
    EXPECT_EQ(Extension::NOT_STARTED, s.anchors->extension);
    EXPECT_EQ(Relation::UNKNOWN, s.anchors->paths.relation);
    request.mode = Mode::ALL_OR_COUNT;
    Result t;
    EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, request, &t).empty());
    EXPECT_EQ(Withheld::THRESHOLD_CROSSED, t.extraction->withheld);
    request.mode = Mode::PARTIAL;
    Result u;
    EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, request, &u).empty());
    EXPECT_FALSE(u.extraction->withheld);
    EXPECT_EQ(StopReason::MAX_ANCHORS, u.extraction->cut);
}

TEST(PatternSearch, ExtensionMaxPaths) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NNNNNNA");
    const std::vector<PathCtx> expected = path_oracle(*graph, pattern, path_request()).paths;
    ASSERT_LT(4u, expected.size());

    // without stop_at_threshold the paths are counted exactly; the release is thresholded
    Request request = path_request();
    request.max_paths = expected.size() - 1;
    request.mode = Mode::ALL_OR_COUNT;
    Result all;
    EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, request, &all).empty());
    EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, all.extraction->withheld);
    EXPECT_EQ(Relation::EXACT, all.anchors->paths.relation);
    EXPECT_EQ(expected.size(), all.anchors->paths.value);
    EXPECT_FALSE(all.stop);

    request.mode = Mode::PARTIAL;
    Result partial;
    auto some = run_paths(engine, StoredNodes(*graph), pattern, request, &partial);
    ASSERT_EQ(expected.size() - 1, some.size());
    EXPECT_TRUE(std::equal(some.begin(), some.end(), expected.begin()));
    EXPECT_EQ(StopReason::MAX_PATHS, partial.extraction->cut);
    EXPECT_EQ(Relation::EXACT, partial.anchors->paths.relation);

    // stop_at_threshold: the extension stops as soon as max_paths is exceeded
    request.stop_at_threshold = true;
    request.max_paths = 2;
    Budget budget = unbounded_budget();
    Result s = engine.count(pattern, request, budget);
    ASSERT_TRUE(s.stop);
    EXPECT_EQ(StopPhase::EXTENSION, s.stop->phase);
    EXPECT_EQ(StopReason::MAX_PATHS, s.stop->reason);
    EXPECT_FALSE(s.time_limited);
    EXPECT_EQ(Relation::EXACT, s.anchors->total.relation);
    EXPECT_EQ(Extension::STOPPED, s.anchors->extension);
    EXPECT_EQ(Relation::AT_LEAST, s.anchors->paths.relation);
    EXPECT_EQ(3u, s.anchors->paths.value);
    for (const auto &[o, count] : s.anchors->paths_by_orientation) {
        EXPECT_NE(Relation::UNKNOWN, count.relation);
    }
    // a threshold stop ends only its own pattern
    EXPECT_FALSE(budget.stopped());
    Result next = engine.count(Pattern::parse(PatternKind::DNA, "GGGCCCAA"), path_request(),
                               budget);
    EXPECT_FALSE(next.stop);
    EXPECT_EQ(Relation::EXACT, next.anchors->paths.relation);

    request.mode = Mode::ALL_OR_COUNT;
    Result crossed;
    EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, request, &crossed).empty());
    EXPECT_EQ(Withheld::THRESHOLD_CROSSED, crossed.extraction->withheld);
    request.mode = Mode::PARTIAL;
    Result first;
    auto two = run_paths(engine, StoredNodes(*graph), pattern, request, &first);
    ASSERT_EQ(2u, two.size());
    EXPECT_TRUE(std::equal(two.begin(), two.end(), expected.begin()));
    EXPECT_EQ(StopReason::MAX_PATHS, first.extraction->cut);
    EXPECT_FALSE(first.extraction->complete);
}

TEST(PatternSearch, ExtensionMaxStepsAtLeast) {
    // every step budget: a stop in discovery leaves the anchors a true lower bound and the
    // paths unknown (nothing extended); a stop in the extension leaves the anchors exact and
    // the paths at_least, PARTIAL a prefix of the answer, ALL_OR_COUNT nothing. In every
    // graph mode, the even-k wrapped PRIMARY graph (palindromic anchors) included, with the
    // truth from the oracles
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA", "GGGCCCAATTGCA" };
    const std::vector<std::pair<size_t, DeBruijnGraph::Mode>> graphs {
        { 5, DeBruijnGraph::BASIC }, { 5, DeBruijnGraph::PRIMARY },
        { 5, DeBruijnGraph::CANONICAL }, { 4, DeBruijnGraph::PRIMARY },
    };
    for (const auto &[k, mode] : graphs) {
        SCOPED_TRACE("k " + std::to_string(k) + " mode " + std::to_string(mode));
        auto graph = build(k, records, mode);
        PatternSearch engine(*graph);
        Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NNNNNRN");
        check_paths_against_oracles(*graph, pattern, path_request(), &records, mode);
        const std::vector<PathCtx> expected = path_oracle(*graph, pattern, path_request()).paths;
        const Truth anchors = truth_of(*graph, pattern, make_request(),
                                       walk_oracle(*graph, pattern, make_request()));
        Result full = count_paths(*graph, pattern, path_request());
        ASSERT_EQ(anchors.total, full.anchors->total.value);
        ASSERT_EQ(Relation::EXACT, full.anchors->paths.relation);
        ASSERT_EQ(expected.size(), full.anchors->paths.value);
        Request plain = make_request();
        const uint64_t discovery = count_of(*graph, pattern, plain).work.steps;
        ASSERT_EQ(discovery + full.work.extension_edges, full.work.steps);

        for (uint64_t steps = 0; steps < full.work.steps; ++steps) {
            Result r = count_paths(*graph, pattern, path_request(), steps);
            ASSERT_TRUE(r.stop) << steps;
            EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);
            EXPECT_LE(r.work.steps, steps);
            const AnchorCounts &a = *r.anchors;
            expect_true_relations(r, anchors, k, pattern.length());
            if (steps < discovery) {
                EXPECT_NE(StopPhase::EXTENSION, r.stop->phase) << steps;
                EXPECT_NE(Relation::EXACT, a.total.relation);
                EXPECT_EQ(Extension::NOT_STARTED, a.extension);
                EXPECT_EQ(Relation::UNKNOWN, a.paths.relation);
                EXPECT_EQ(0u, r.work.extension_edges);
            } else {
                EXPECT_EQ(StopPhase::EXTENSION, r.stop->phase) << steps;
                EXPECT_EQ(Relation::EXACT, a.total.relation);
                EXPECT_EQ(full.anchors->total.value, a.total.value);
                EXPECT_EQ(Extension::STOPPED, a.extension);
                EXPECT_EQ(Relation::AT_LEAST, a.paths.relation);
                EXPECT_LE(a.paths.value, expected.size());
                EXPECT_EQ(steps, r.work.steps);
                EXPECT_EQ(steps - discovery, r.work.extension_edges);
                EXPECT_LE(a.candidates_examined, full.anchors->candidates_examined);
                for (const auto &[o, count] : a.paths_by_orientation) {
                    const Count &truth = full.anchors->paths_by_orientation.at(o);
                    if (count.relation == Relation::EXACT) {
                        EXPECT_EQ(truth.value, count.value) << steps;
                    } else {
                        EXPECT_EQ(Relation::AT_LEAST, count.relation);
                        EXPECT_LE(count.value, truth.value);
                    }
                }
            }

            Request partial = path_request();
            partial.mode = Mode::PARTIAL;
            Result p;
            auto prefix = run_paths(engine, StoredNodes(*graph), pattern, partial, &p, steps);
            EXPECT_EQ(StopReason::MAX_STEPS, p.extraction->cut);
            EXPECT_FALSE(p.extraction->complete);
            ASSERT_LE(prefix.size(), expected.size());
            EXPECT_TRUE(std::equal(prefix.begin(), prefix.end(), expected.begin())) << steps;
            EXPECT_EQ(steps < discovery ? 0u : a.paths.value, prefix.size());

            Request all = path_request();
            all.mode = Mode::ALL_OR_COUNT;
            Result w;
            EXPECT_TRUE(run_paths(engine, StoredNodes(*graph), pattern, all, &w, steps).empty());
            EXPECT_EQ(Withheld::DISCOVERY_BUDGET, w.extraction->withheld);
        }
    }
}

TEST(PatternSearch, ExtensionDeadline) {
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Pattern pattern = Pattern::parse(PatternKind::DNA, "GATCGATC");
    const auto start = Deadline::Clock::now();

    // the clock reads late from its |on_time|+1-th reading on: the pattern's start reads it
    // first, the listing of the anchors second, then the extension before each anchor (one
    // anchor here, GATCG), the release of the paths after them
    auto run = [&](int on_time, Mode mode, std::vector<PathCtx> *released) {
        auto readings = std::make_shared<int>(0);
        auto clock = [start, readings, on_time]() {
            return start + std::chrono::milliseconds(++*readings > on_time ? 10'000 : 0);
        };
        Budget budget(kManySteps, Deadline(start, 1000, 250, clock));
        Request request = path_request();
        request.mode = mode;
        if (!released)
            return engine.count(pattern, request, budget);
        return engine.enumerate(pattern, request, budget, [&](const Context &c) {
            released->push_back(PathCtx { c.node, c.orientation, c.sequence, c.path });
        });
    };

    // late before the pattern: nothing runs
    Result before = run(0, Mode::COUNT, nullptr);
    EXPECT_EQ(StopPhase::DISCOVERY, before.stop->phase);
    EXPECT_EQ(StopReason::TIME, before.stop->reason);
    EXPECT_EQ(Extension::NOT_STARTED, before.anchors->extension);
    EXPECT_EQ(Relation::UNKNOWN, before.anchors->paths.relation);
    EXPECT_TRUE(before.time_limited);

    // late at the extension: anchors exact, paths at_least 0
    Result extension = run(1, Mode::COUNT, nullptr);
    ASSERT_TRUE(extension.stop);
    EXPECT_EQ(StopPhase::EXTENSION, extension.stop->phase);
    EXPECT_EQ(StopReason::TIME, extension.stop->reason);
    EXPECT_EQ(Relation::EXACT, extension.anchors->total.relation);
    EXPECT_EQ(Extension::STOPPED, extension.anchors->extension);
    EXPECT_EQ(Relation::AT_LEAST, extension.anchors->paths.relation);
    EXPECT_EQ(0u, extension.anchors->paths.value);
    EXPECT_TRUE(extension.time_limited);
    std::vector<PathCtx> none;
    Result all = run(1, Mode::ALL_OR_COUNT, &none);
    EXPECT_EQ(Withheld::DEADLINE, all.extraction->withheld);
    Result partial = run(1, Mode::PARTIAL, &none);
    EXPECT_EQ(StopReason::TIME, partial.extraction->cut);
    EXPECT_TRUE(none.empty());

    // late before the anchor's extension: the anchors listed, none extended
    Result anchor = run(2, Mode::COUNT, nullptr);
    ASSERT_TRUE(anchor.stop);
    EXPECT_EQ(StopPhase::EXTENSION, anchor.stop->phase);
    EXPECT_EQ(StopReason::TIME, anchor.stop->reason);
    EXPECT_EQ(1u, anchor.anchors->total.value);
    EXPECT_EQ(Extension::STOPPED, anchor.anchors->extension);
    EXPECT_EQ(Relation::AT_LEAST, anchor.anchors->paths.relation);
    EXPECT_EQ(0u, anchor.anchors->paths.value);
    EXPECT_EQ(0u, anchor.work.extension_anchors);
    EXPECT_EQ(0u, anchor.work.extension_edges);

    // late at the release: paths exact, nothing released, stop {extraction, time}
    Result release = run(3, Mode::ALL_OR_COUNT, &none);
    EXPECT_EQ(Relation::EXACT, release.anchors->paths.relation);
    EXPECT_LT(0u, release.anchors->paths.value);
    ASSERT_TRUE(release.stop);
    EXPECT_EQ(StopPhase::EXTRACTION, release.stop->phase);
    EXPECT_EQ(Withheld::DEADLINE, release.extraction->withheld);
    Result streamed = run(3, Mode::PARTIAL, &none);
    EXPECT_EQ(StopReason::TIME, streamed.extraction->cut);
    EXPECT_TRUE(none.empty());

    // on time throughout
    std::vector<PathCtx> got;
    Result fine = run(100, Mode::ALL_OR_COUNT, &got);
    EXPECT_FALSE(fine.stop);
    EXPECT_TRUE(fine.extraction->complete);
    EXPECT_EQ(fine.anchors->paths.value, got.size());
}

TEST(PatternSearch, ExtensionDeterministicOrder) {
    // the same request on the same index gives the same paths in the same order (§5.5):
    // anchors by node, then orientation, then the DFS in symbol order
    std::mt19937 rng(7);
    std::vector<std::string> records(4);
    for (std::string &record : records) {
        for (int i = 0; i < 40; ++i) {
            record.push_back("ACGT"[rng() % 4]);
        }
    }
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(4, records, mode);
        Pattern pattern = Pattern::parse(PatternKind::IUPAC, "NNRNYN");
        auto first = paths_of(*graph, pattern, path_request());
        auto second = paths_of(*graph, pattern, path_request());
        ASSERT_LT(10u, first.size());
        EXPECT_EQ(first, second);
        EXPECT_TRUE(std::is_sorted(first.begin(), first.end()));
        // within one anchor and orientation, the paths come by their spelled sequence
        for (size_t i = 1; i < first.size(); ++i) {
            if (first[i].anchor == first[i - 1].anchor
                    && first[i].orientation == first[i - 1].orientation) {
                EXPECT_LT(first[i - 1].sequence, first[i].sequence);
            }
        }
        check_paths_against_oracles(*graph, pattern, path_request(), &records, mode);
    }
}

TEST(PatternSearch, ExtensionRandomGraphsAgainstOracles) {
    // fixed seeds: a few hundred (graph, long pattern, request) cases over every graph mode,
    // both builders, every strand choice, DNA and IUPAC patterns of length k + 1 .. k + 6
    const char *iupac = "ACGTRYSWKMBDHVN";
    size_t cases = 0;
    size_t nonempty = 0;
    size_t multiple = 0;
    size_t max_expected = 0;
    for (uint32_t seed = 1; seed <= 60; ++seed) {
        std::mt19937 rng(1000 + seed);
        size_t k = 3 + rng() % 6;
        auto mode = static_cast<DeBruijnGraph::Mode>(rng() % 3);
        bool batch = rng() % 2;

        std::vector<std::string> records(1 + rng() % 5);
        for (std::string &record : records) {
            size_t length = k + 2 + rng() % 30;
            for (size_t i = 0; i < length; ++i) {
                // an occasional N splits a record into islands (DNA4)
                record.push_back(rng() % 25 ? "ACGT"[rng() % 4] : 'N');
            }
        }
        auto graph = build(k, records, mode, batch);

        for (int t = 0; t < 6; ++t) {
            size_t length = k + 1 + rng() % 6;
            bool use_iupac = rng() % 2;
            std::string text;
            for (size_t i = 0; i < length; ++i) {
                text.push_back(use_iupac && rng() % 3 == 0 ? iupac[rng() % 15] : "ACGT"[rng() % 4]);
            }
            // three quarters of the patterns from the records' islands, so that most have
            // paths, and those of them that are IUPAC made degenerate at two positions
            if (rng() % 4) {
                for (int attempt = 0; attempt < 10; ++attempt) {
                    const std::string &record = records[rng() % records.size()];
                    if (record.size() < length)
                        continue;
                    std::string piece = record.substr(rng() % (record.size() - length + 1),
                                                      length);
                    if (piece.find('N') != std::string::npos)
                        continue;
                    text = piece;
                    for (int d = 0; use_iupac && d < 2; ++d) {
                        text[rng() % length] = iupac[rng() % 15];
                    }
                    break;
                }
            }
            auto strands = static_cast<Strands>(rng() % 3);
            Pattern pattern = Pattern::parse(PatternKind::IUPAC, text);
            SCOPED_TRACE("seed " + std::to_string(seed) + " k " + std::to_string(k)
                         + " mode " + std::to_string(mode) + " batch " + std::to_string(batch)
                         + " pattern " + text);
            size_t num_expected = 0;
            check_paths_against_oracles(*graph, pattern,
                                        path_request(Scope::ANY_OFFSET, strands),
                                        &records, mode, false, {}, &num_expected);
            ++cases;
            nonempty += num_expected > 0;
            multiple += num_expected > 1;
            max_expected = std::max(max_expected, num_expected);
        }
    }
    EXPECT_EQ(360u, cases);
    // the cases are not vacuous: many have paths, some several
    EXPECT_LT(150u, nonempty);
    EXPECT_LT(50u, multiple);
    std::cerr << "random long cases: " << cases << ", with paths: " << nonempty
              << ", with more than 1: " << multiple << ", most: " << max_expected << std::endl;
}


// ---------------------------------------------------------------- time stops inside a phase

TEST(PatternSearch, TimeStopMidDiscovery) {
    // A time stop inside discovery, not only before the pattern starts. An injected
    // clock is on time for its first n readings and late from then on, for every n below the
    // number of readings a complete run takes (whatever the engine's reading points are: the
    // pattern's start, the stride crossings of Budget::kClockStride steps, the release).
    // Every such run stops with reason time, time_limited, never an exact total, and every
    // count a true relation of the oracle's truth; a stop inside discovery (after steps were
    // charged) leaves the total at_least, ALL_OR_COUNT withholds deadline, PARTIAL returns
    // nothing (cut time), and the next pattern of the request answers unknown
    std::mt19937 rng(11);
    std::vector<std::string> records(3);
    for (std::string &record : records) {
        for (int i = 0; i < 1500; ++i) {
            record.push_back("ACGT"[rng() % 4]);
        }
    }
    auto graph = build(9, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    // A in any_offset at k = 9: the flank ranges of eight offsets, far more than a stride
    Pattern pattern = Pattern::parse(PatternKind::DNA, "A");
    const std::vector<Ctx> expected = walk_oracle(*graph, pattern, make_request());
    const Truth truth = truth_of(*graph, pattern, make_request(), expected);
    const Result full = count_of(*graph, pattern, make_request());
    ASSERT_EQ(truth.total, full.contexts->total.value);
    ASSERT_LT(2 * Budget::kClockStride, full.work.ranges_visited);  // discovery spans strides

    // a budget whose clock reads on time |on_time| times, then late; |readings| counts them
    const auto start = Deadline::Clock::now();
    auto late_from = [start](int on_time, std::shared_ptr<int> readings = nullptr) {
        if (!readings)
            readings = std::make_shared<int>(0);
        return Budget(kManySteps, Deadline(start, 1000, 250, [start, readings, on_time]() {
            return start + std::chrono::milliseconds(++*readings > on_time ? 10'000 : 0);
        }));
    };
    auto readings_of_complete_run = [&](const PatternSearch &search, const Pattern &p,
                                        const Request &request) {
        auto readings = std::make_shared<int>(0);
        Budget budget = late_from(std::numeric_limits<int>::max(), readings);
        Result r = search.count(p, request, budget);
        EXPECT_FALSE(r.stop);
        return *readings;
    };

    const int readings = readings_of_complete_run(engine, pattern, make_request());
    ASSERT_LT(2, readings);
    int inside = 0;
    int first_inside = -1;
    for (int on_time = 0; on_time < readings; ++on_time) {
        SCOPED_TRACE("on time for " + std::to_string(on_time) + " of " + std::to_string(readings)
                     + " readings");
        Budget budget = late_from(on_time);
        Result r = engine.count(pattern, make_request(), budget);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_TRUE(r.time_limited);
        EXPECT_NE(Relation::EXACT, r.contexts->total.relation);
        expect_true_relations(r, truth, 9, 1);
        if (r.stop->phase != StopPhase::DISCOVERY || !r.work.steps)
            continue;
        // inside discovery: a lower bound, and no offset exact (every offset sums both
        // orientations, and the stopped one is open)
        ++inside;
        if (first_inside < 0)
            first_inside = on_time;
        EXPECT_EQ(Relation::AT_LEAST, r.contexts->total.relation);
        EXPECT_LT(r.work.steps, full.work.steps);
        for (const auto &[p, count] : r.contexts->by_offset) {
            EXPECT_NE(Relation::EXACT, count.relation) << p;
        }

        // the stop is sticky: the next pattern of the request does not start
        Result next = engine.count(Pattern::parse(PatternKind::DNA, "ACGTAC"), make_request(),
                                   budget);
        ASSERT_TRUE(next.stop);
        EXPECT_EQ(StopPhase::DISCOVERY, next.stop->phase);
        EXPECT_EQ(StopReason::TIME, next.stop->reason);
        EXPECT_EQ(Relation::UNKNOWN, next.contexts->total.relation);
        EXPECT_EQ(0u, next.work.steps);
    }
    EXPECT_LE(2, inside);  // discovery spans at least two stride crossings
    ASSERT_LE(0, first_inside);

    // the release after a stop inside discovery: nothing, whatever the mode
    Request all = make_request();
    all.mode = Mode::ALL_OR_COUNT;
    Budget budget = late_from(first_inside);
    Result withheld = engine.enumerate(pattern, all, budget, [](const Context &) { FAIL(); });
    EXPECT_EQ(StopPhase::DISCOVERY, withheld.stop->phase);
    EXPECT_EQ(Withheld::DEADLINE, withheld.extraction->withheld);
    EXPECT_EQ(Relation::AT_LEAST, withheld.contexts->total.relation);
    EXPECT_FALSE(withheld.extraction->complete);

    Request partial = make_request();
    partial.mode = Mode::PARTIAL;
    budget = late_from(first_inside);
    Result cut = engine.enumerate(pattern, partial, budget, [](const Context &) { FAIL(); });
    EXPECT_EQ(StopPhase::DISCOVERY, cut.stop->phase);
    EXPECT_EQ(StopReason::TIME, cut.extraction->cut);
    EXPECT_EQ(0u, cut.extraction->returned);
    EXPECT_FALSE(cut.extraction->complete);
    EXPECT_TRUE(cut.time_limited);

    // a long pattern stopped by time inside its extension's DFS: the anchors exact, the
    // paths at_least and a lower bound of the oracle's, nothing released in ALL_OR_COUNT
    Pattern long_pattern = Pattern::parse(PatternKind::IUPAC, "AC" + std::string(12, 'N'));
    auto graph7 = build(7, records, DeBruijnGraph::BASIC);
    PatternSearch extender(*graph7);
    Request paths = make_request(Scope::ANY_OFFSET, Strands::FORWARD);
    paths.extend_paths = true;
    paths.max_paths = kManySteps;
    const PathOracle path_truth = path_oracle(*graph7, long_pattern, paths);
    Budget unbounded = unbounded_budget();
    const Result complete = extender.count(long_pattern, paths, unbounded);
    ASSERT_EQ(Extension::COMPLETED, complete.anchors->extension);
    ASSERT_EQ(path_truth.paths.size(), complete.anchors->paths.value);
    // the extension spans more than two stride crossings
    ASSERT_LT(complete.work.ranges_visited + 2 * Budget::kClockStride, complete.work.steps);

    const int long_readings = readings_of_complete_run(extender, long_pattern, paths);
    int in_extension = 0;
    for (int on_time = 0; on_time < long_readings; ++on_time) {
        SCOPED_TRACE("long pattern, on time for " + std::to_string(on_time) + " of "
                     + std::to_string(long_readings) + " readings");
        budget = late_from(on_time);
        Result r = extender.count(long_pattern, paths, budget);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_TRUE(r.time_limited);
        EXPECT_NE(Relation::EXACT, r.anchors->paths.relation);
        EXPECT_TRUE(true_relation(r.anchors->paths, path_truth.paths.size()));
        if (r.stop->phase != StopPhase::EXTENSION || !r.work.extension_edges)
            continue;
        ++in_extension;
        EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
        EXPECT_EQ(complete.anchors->total.value, r.anchors->total.value);
        EXPECT_EQ(Extension::STOPPED, r.anchors->extension);
        EXPECT_EQ(Relation::AT_LEAST, r.anchors->paths.relation);
        for (const auto &[o, count] : r.anchors->paths_by_orientation) {
            EXPECT_EQ(Relation::AT_LEAST, count.relation) << orientation_key(o);
            EXPECT_LE(count.value, complete.anchors->paths_by_orientation.at(o).value);
        }
        Request paths_all = paths;
        paths_all.mode = Mode::ALL_OR_COUNT;
        budget = late_from(on_time);
        Result none = extender.enumerate(long_pattern, paths_all, budget,
                                         [](const Context &) { FAIL(); });
        EXPECT_EQ(Withheld::DEADLINE, none.extraction->withheld);
        EXPECT_EQ(Relation::AT_LEAST, none.anchors->paths.relation);
    }
    EXPECT_LE(2, in_extension);
    std::cerr << "clock readings: discovery run " << readings << " (" << inside
              << " stops inside discovery), extension run " << long_readings << " ("
              << in_extension << " stops inside the DFS)" << std::endl;
}

#endif // _DNA_GRAPH || _DNA5_GRAPH

} // namespace
