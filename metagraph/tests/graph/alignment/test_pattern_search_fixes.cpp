/**
 * Regression tests of the pattern engine's fixes after the review of 2026-10-07 (milestones
 * 1/1b): one test (or more) per finding, each failing without its fix. The oracle here is
 * its own: the IUPAC table, the reverse complement and the matching are written from the
 * pattern TEXT with this file's tables (never Pattern::positions, is_palindromic or
 * reverse_complement), over the served graph's k-mers (the valid edges of the DBGSuccinct,
 * and on a wrapped PRIMARY graph their reverse complements at CanonicalDBG's ids). What stays
 * shared with the engine: the DBGSuccinct mask (in_graph) and CanonicalDBG's numbering.
 * The findings of GPT review 3 (2026-10-08) close the file: the low-complexity note against
 * sdust run over the whole text (the library the engine runs piece by piece), and the clock
 * readings of the diagnostic and of the extension.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
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

#include "common/vectors/bit_vector_dyn.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include <sdust.h>


namespace {

#if ! _PROTEIN_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;

typedef DeBruijnGraph::node_index node_index;
typedef Deadline::Clock Clock;

constexpr uint64_t kSteps = 1'000'000'000;


// ---------------------------------------------------------------- an oracle of its own

// the IUPAC set of a letter (A 1, C 2, G 4, T 8), written out here
uint8_t own_set(char c) {
    switch (c) {
        case 'A': return 1;  case 'C': return 2;  case 'G': return 4;  case 'T': return 8;
        case 'R': return 1 | 4;  case 'Y': return 2 | 8;  case 'S': return 2 | 4;
        case 'W': return 1 | 8;  case 'K': return 4 | 8;  case 'M': return 1 | 2;
        case 'B': return 2 | 4 | 8;  case 'D': return 1 | 4 | 8;  case 'H': return 1 | 2 | 8;
        case 'V': return 1 | 2 | 4;  case 'N': return 15;
        default: return 0;
    }
}

char own_complement(char c) {
    switch (c) {
        case 'A': return 'T';  case 'T': return 'A';  case 'C': return 'G';  case 'G': return 'C';
        case 'R': return 'Y';  case 'Y': return 'R';  case 'K': return 'M';  case 'M': return 'K';
        case 'B': return 'V';  case 'V': return 'B';  case 'D': return 'H';  case 'H': return 'D';
        case 'S': return 'S';  case 'W': return 'W';  case 'N': return 'N';
        default: return '?';
    }
}

std::string own_rc(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (char &c : r) {
        c = own_complement(c);
    }
    return r;
}

// a graph base (A, C, G, T only: never $ or the graph's N) admitted by a pattern letter
bool own_match(std::string_view pattern, std::string_view s) {
    if (pattern.size() != s.size())
        return false;
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t base = s[i] == 'A' ? 1 : s[i] == 'C' ? 2 : s[i] == 'G' ? 4 : s[i] == 'T' ? 8 : 0;
        if (!(own_set(pattern[i]) & base))
            return false;
    }
    return true;
}

const DBGSuccinct& base_of(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// the served graph's k-mers: the valid edges, and on a wrapped PRIMARY graph the reverse
// complement of each non-palindromic one at id + max_index (CanonicalDBG's numbering)
std::vector<std::pair<node_index, std::string>> served_kmers(const DeBruijnGraph &graph) {
    const DBGSuccinct &dbg = base_of(graph);
    const bool wrapped = dynamic_cast<const CanonicalDBG*>(&graph);
    std::vector<std::pair<node_index, std::string>> kmers;
    for (node_index y = 1; y <= dbg.max_index(); ++y) {
        if (!dbg.in_graph(y))
            continue;
        std::string s = dbg.get_node_sequence(y);
        kmers.emplace_back(y, s);
        if (wrapped && own_rc(s) != s)
            kmers.emplace_back(y + dbg.max_index(), own_rc(s));
    }
    return kmers;
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
    return out << "(" << c.node << ", " << c.offset << ", " << orientation_key(c.orientation)
               << ")";
}

// the oriented patterns of a request, from the text: a palindromic text once
std::vector<std::pair<Orientation, std::string>> own_orientations(const std::string &text,
                                                                  Strands strands) {
    if (own_rc(text) == text)
        return { { Orientation::PALINDROMIC, text } };
    std::vector<std::pair<Orientation, std::string>> result;
    if (strands != Strands::REVERSE)
        result.emplace_back(Orientation::FORWARD, text);
    if (strands != Strands::FORWARD)
        result.emplace_back(Orientation::REVERSE, own_rc(text));
    return result;
}

/**
 * Every context (L <= k) or anchor (L > k: the k-mers instantiating the oriented pattern's
 * first k positions, offset 0) of the served graph, in answer order.
 */
std::vector<Ctx> oracle(const std::vector<std::pair<node_index, std::string>> &kmers,
                        size_t k, const std::string &text, Scope scope, Strands strands) {
    const size_t L = text.size();
    std::vector<Ctx> result;
    for (const auto &[node, kmer] : kmers) {
        for (const auto &[orientation, q] : own_orientations(text, strands)) {
            if (L > k) {
                if (own_match(std::string_view(q).substr(0, k), kmer))
                    result.push_back(Ctx { node, 0, orientation });
                continue;
            }
            for (uint32_t p = 0; p + L <= k; ++p) {
                if (scope == Scope::SUFFIX && p + L != k)
                    continue;
                if (own_match(q, std::string_view(kmer).substr(p, L)))
                    result.push_back(Ctx { node, p, orientation });
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<Ctx> oracle(const DeBruijnGraph &graph, const std::string &text, Scope scope,
                        Strands strands) {
    return oracle(served_kmers(graph), graph.get_k(), text, scope, strands);
}

std::shared_ptr<DeBruijnGraph> build(size_t k, const std::vector<std::string> &records,
                                     DeBruijnGraph::Mode mode) {
    return build_graph<DBGSuccinct>(k, records, mode);
}

std::vector<std::string> random_records(size_t num, size_t length, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<std::string> records;
    for (size_t i = 0; i < num; ++i) {
        std::string r;
        for (size_t j = 0; j < length; ++j) {
            r.push_back("ACGT"[rng() % 4]);
        }
        records.push_back(r);
    }
    return records;
}

Request request_of(Mode mode, uint64_t max_contexts, Scope scope = Scope::ANY_OFFSET,
                   Strands strands = Strands::BOTH) {
    Request request;
    request.mode = mode;
    request.scope = scope;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = max_contexts;
    request.max_anchors = max_contexts;
    return request;
}

Budget budget_of(uint64_t max_steps = kSteps) {
    return Budget(max_steps, Deadline::unbounded());
}

Pattern iupac(const std::string &text) {
    return Pattern::parse(PatternKind::IUPAC, text);
}

std::vector<Ctx> released_by(const PatternSearch &engine, const Pattern &pattern,
                             const Request &request, Budget &budget, Result *result) {
    std::vector<Ctx> released;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        released.push_back(Ctx { c.node, c.offset, c.orientation });
    });
    return released;
}

// a clock that reads |start| plus |*now| milliseconds; tests move |*now|
struct VirtualClock {
    Clock::time_point start = Clock::now();
    std::shared_ptr<double> now = std::make_shared<double>(0);

    std::function<Clock::time_point()> clock() const {
        auto s = start;
        auto n = now;
        return [s, n]() {
            return s + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double, std::milli>(*n));
        };
    }
    Deadline deadline(double budget_ms, double reserve_ms) const {
        return Deadline(start, budget_ms, reserve_ms, clock());
    }
};


// ---------------------------------------------------------------- E4-01, X-EFFICIENCY-03,
// X-GUARANTEES-02: PARTIAL keeps only what its release can use

// one graph per mode with many ranges for a one-base pattern
struct Broad {
    std::vector<std::string> records = random_records(10, 2000, 101);
    std::shared_ptr<DeBruijnGraph> graph;
    std::vector<Ctx> expected;
};

TEST(PatternSearchFixes, PartialRetentionBoundedByWhatItReleases) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        for (size_t k : { 13, 12 }) {
            if (mode != DeBruijnGraph::PRIMARY && k == 12)
                continue;
            Broad b;
            b.graph = build(k, b.records, mode);
            PatternSearch engine(*b.graph);
            const std::string text = "A";
            b.expected = oracle(*b.graph, text, Scope::ANY_OFFSET, Strands::BOTH);
            ASSERT_GT(b.expected.size(), 50'000u);

            // everything retained: ALL_OR_COUNT under a threshold above the count
            Budget all_budget = budget_of();
            Result all;
            std::vector<Ctx> full = released_by(engine, iupac(text),
                                                request_of(Mode::ALL_OR_COUNT, 1'000'000'000),
                                                all_budget, &all);
            ASSERT_EQ(b.expected, full) << "mode " << mode << " k " << k;
            const uint64_t retained_all = all.work.spans_retained_peak;
            ASSERT_GT(retained_all, 20'000u);

            for (uint64_t cap : { 0, 1, 7, 100, 4097, 20'000 }) {
                Budget budget = budget_of();
                Result r;
                std::vector<Ctx> got = released_by(engine, iupac(text),
                                                   request_of(Mode::PARTIAL, cap), budget, &r);
                // the first |cap| contexts in answer order, exactly
                ASSERT_EQ(std::vector<Ctx>(full.begin(), full.begin() + cap), got)
                        << "mode " << mode << " k " << k << " cap " << cap;
                EXPECT_EQ(StopReason::MAX_CONTEXTS, r.extraction->cut);
                EXPECT_EQ(all.work.steps, r.work.steps);
                // bounded by the cap, not by the discovery (24 bytes each): the descriptors
                // kept are about cap plus the straddling ranges, never above twice the
                // compaction threshold
                EXPECT_LE(r.work.spans_retained_peak,
                          std::max<uint64_t>(4096, 2 * (cap + 500)) + 1)
                        << "mode " << mode << " k " << k << " cap " << cap;
                if (cap <= 100) {
                    EXPECT_LE(r.work.spans_retained_peak, 4097u);
                }
                if (!cap) {
                    EXPECT_EQ(0u, r.work.spans_retained_peak);
                }
            }

            // a cap at and above the count: everything, complete
            for (uint64_t cap : { b.expected.size(), b.expected.size() + 1 }) {
                Budget budget = budget_of();
                Result r;
                std::vector<Ctx> got = released_by(engine, iupac(text),
                                                   request_of(Mode::PARTIAL, cap), budget, &r);
                EXPECT_EQ(full, got);
                EXPECT_TRUE(r.extraction->complete);
            }
        }
    }
}

TEST(PatternSearchFixes, PartialRetentionUnderStepStops) {
    // after a max_steps stop, PARTIAL releases the first cap among the contexts discovered:
    // the pruned release equals the unpruned one cut at cap
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
        for (size_t k : { 11, 10 }) {
            if (mode == DeBruijnGraph::BASIC && k == 10)
                continue;
            auto graph = build(k, random_records(8, 1500, 7 + k), mode);
            PatternSearch engine(*graph);
            for (const std::string &text : std::vector<std::string> { "C", "RA", "NNG" }) {
                for (uint64_t steps : { 500, 3000, 9000, 20000 }) {
                    Budget unpruned_budget = budget_of(steps);
                    Result unpruned;
                    std::vector<Ctx> reference = released_by(
                            engine, iupac(text), request_of(Mode::PARTIAL, 1'000'000'000),
                            unpruned_budget, &unpruned);
                    for (uint64_t cap : { 1, 50, 5000 }) {
                        Budget budget = budget_of(steps);
                        Result r;
                        std::vector<Ctx> got = released_by(
                                engine, iupac(text), request_of(Mode::PARTIAL, cap), budget, &r);
                        std::vector<Ctx> prefix(reference.begin(),
                                                reference.begin()
                                                    + std::min<size_t>(cap, reference.size()));
                        ASSERT_EQ(prefix, got) << text << " steps " << steps << " cap " << cap
                                               << " mode " << mode << " k " << k;
                    }
                }
            }
        }
    }
}

TEST(PatternSearchFixes, ReleaseSetupReadsTheClock) {
    // many cursors (no pruning under a huge cap): a deadline that passes while the release
    // prepares them stops it before any context, {extraction, time}
    auto graph = build(13, random_records(12, 2500, 101), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const Pattern pattern = iupac("A");

    // the clock's readings before the release: the pattern's start, the discovery's
    // stride crossings, the release's start
    Budget probe = budget_of();
    Result counted = engine.count(pattern, request_of(Mode::COUNT, 0), probe);
    const uint64_t before = 1 + counted.work.steps / Budget::kClockStride;

    auto readings = std::make_shared<uint64_t>(0);
    const Clock::time_point start = Clock::now();
    auto clock = [start, readings, before]() {
        // on time through the release's start, late from the next reading on
        return start + std::chrono::seconds(++*readings > before + 1 ? 100 : 0);
    };
    Budget budget(kSteps, Deadline(start, 1000, 250, clock));
    Result r;
    std::vector<Ctx> got = released_by(engine, pattern,
                                       request_of(Mode::PARTIAL, 1'000'000'000), budget, &r);
    EXPECT_TRUE(got.empty());
    ASSERT_TRUE(r.stop);
    EXPECT_EQ(StopPhase::EXTRACTION, r.stop->phase);
    EXPECT_EQ(StopReason::TIME, r.stop->reason);
    EXPECT_EQ(StopReason::TIME, r.extraction->cut);
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
}


// ---------------------------------------------------------------- R1-01, X-EFFICIENCY-02,
// C1-02, D1-06, E4-02, X-GUARANTEES-03, E4-03: the caller's work per context is clocked

std::shared_ptr<DeBruijnGraph> medium_graph() {
    return build(9, random_records(6, 600, 11), DeBruijnGraph::BASIC);
}

TEST(PatternSearchFixes, AllOrCountDeliveryUnderTheClock) {
    auto graph = medium_graph();
    PatternSearch engine(*graph);
    const Pattern pattern = iupac("AC");
    std::vector<Ctx> expected = oracle(*graph, "AC", Scope::ANY_OFFSET, Strands::BOTH);
    ASSERT_GT(expected.size(), 300u);

    // every context costs the caller 1 ms of a virtual clock; the work time is 100 ms
    VirtualClock vc;
    Budget budget(kSteps, vc.deadline(350, 250));
    std::vector<Ctx> received;
    Result r = engine.enumerate(pattern, request_of(Mode::ALL_OR_COUNT, 1'000'000), budget,
                                [&](const Context &c) {
        received.push_back(Ctx { c.node, c.offset, c.orientation });
        *vc.now += 1;
    });
    // the counts stay exact; the release is withheld, nothing is claimed complete
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
    EXPECT_EQ(expected.size(), r.contexts->total.value);
    ASSERT_TRUE(r.extraction);
    EXPECT_EQ(Withheld::DEADLINE, r.extraction->withheld);
    EXPECT_FALSE(r.extraction->complete);
    EXPECT_EQ(0u, r.extraction->returned);
    ASSERT_TRUE(r.stop);
    EXPECT_EQ(StopPhase::EXTRACTION, r.stop->phase);
    EXPECT_EQ(StopReason::TIME, r.stop->reason);
    EXPECT_TRUE(r.time_limited);
    // the caller got a prefix (to discard), stopped within one stride of the work time
    EXPECT_LT(received.size(), 100u + Budget::kReleaseClockStride);
    EXPECT_GE(received.size(), 100u);
    EXPECT_EQ(std::vector<Ctx>(expected.begin(), expected.begin() + received.size()), received);

    // a delivery that fits the work time is complete, and reads nothing it must not
    VirtualClock fast;
    Budget enough(kSteps, fast.deadline(10'000, 250));
    received.clear();
    r = engine.enumerate(pattern, request_of(Mode::ALL_OR_COUNT, 1'000'000), enough,
                         [&](const Context &c) {
        received.push_back(Ctx { c.node, c.offset, c.orientation });
        *fast.now += 1;
    });
    EXPECT_TRUE(r.extraction->complete);
    EXPECT_EQ(expected, received);
    EXPECT_FALSE(r.stop);
    EXPECT_FALSE(r.time_limited);
}

TEST(PatternSearchFixes, PartialReleaseReadsTheClockPerReleasedContexts) {
    auto graph = medium_graph();
    PatternSearch engine(*graph);
    const Pattern pattern = iupac("AC");
    std::vector<Ctx> expected = oracle(*graph, "AC", Scope::ANY_OFFSET, Strands::BOTH);
    ASSERT_GT(expected.size(), 300u);

    VirtualClock vc;
    Budget budget(kSteps, vc.deadline(260, 250));
    std::vector<Ctx> received;
    Result r = engine.enumerate(pattern, request_of(Mode::PARTIAL, 1'000'000), budget,
                                [&](const Context &c) {
        received.push_back(Ctx { c.node, c.offset, c.orientation });
        *vc.now += 1;
    });
    // the work time (10 ms) passes at the 10th context; the next reading is at the 64th
    EXPECT_EQ(StopReason::TIME, r.extraction->cut);
    EXPECT_EQ(Budget::kReleaseClockStride, received.size());
    EXPECT_EQ(received.size(), r.extraction->returned);
    EXPECT_EQ(std::vector<Ctx>(expected.begin(), expected.begin() + received.size()), received);
    ASSERT_TRUE(r.stop);
    EXPECT_EQ(StopPhase::EXTRACTION, r.stop->phase);
    EXPECT_TRUE(r.time_limited);
}


// ---------------------------------------------------------------- X-CONCURRENCY-01, R2-02:
// a departed caller stops the work

TEST(PatternSearchFixes, AbortAtAClockReading) {
    auto graph = build(13, random_records(12, 2500, 101), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const Pattern pattern = iupac("A");

    auto asked = std::make_shared<uint64_t>(0);
    Budget budget = budget_of();
    budget.set_abort([asked]() { return ++*asked > 3; });
    EXPECT_THROW(engine.count(pattern, request_of(Mode::COUNT, 0), budget), Aborted);
    EXPECT_EQ(4u, *asked);
    // fewer steps than a run to the end: the stop came at the fourth reading
    Budget full = budget_of();
    Result r = engine.count(pattern, request_of(Mode::COUNT, 0), full);
    EXPECT_LT(budget.steps_used(), r.work.steps);

    // asked and never true: nothing changes
    Budget calm = budget_of();
    calm.set_abort([]() { return false; });
    Result same = engine.count(pattern, request_of(Mode::COUNT, 0), calm);
    EXPECT_EQ(r.contexts->total.value, same.contexts->total.value);
    EXPECT_EQ(r.work.steps, same.work.steps);

    // in the release: the reading at its start
    auto in_release = std::make_shared<bool>(false);
    Budget late = budget_of();
    late.set_abort([in_release]() { return *in_release; });
    bool thrown = false;
    try {
        Request partial = request_of(Mode::PARTIAL, 1'000'000'000);
        engine.enumerate(pattern, partial, late, [&](const Context &) { *in_release = true; });
    } catch (const Aborted &) {
        thrown = true;
    }
    EXPECT_TRUE(thrown);
}


// ---------------------------------------------------------------- X-EFFICIENCY-01: a leading
// pattern-N run is not branched

TEST(PatternSearchFixes, LeadingNRunIsNotSearched) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        const size_t k = 15;
        auto graph = build(k, random_records(10, 3000, 23), mode);
        PatternSearch engine(*graph);
        const auto kmers = served_kmers(*graph);
        // a core taken from the graph, so that it has contexts
        const std::string core = kmers[100].second.substr(3, 6);
        // (owner decision #8: runs at both ends at once too, L = k and L < k)
        for (const std::string &text : { "NNNNNNNN" + core, core + "NNNNNNNN",
                                        "NNNNN" + core, "NNNNNNNNN" + core,
                                        "NNNN" + core + "NNNNN", "NNN" + core + "NNN" }) {
            for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
                for (Scope scope : { Scope::ANY_OFFSET, Scope::SUFFIX }) {
                    if (scope == Scope::SUFFIX && mode == DeBruijnGraph::PRIMARY)
                        continue;
                    std::vector<Ctx> expected = oracle(kmers, k, text, scope, strands);
                    Budget budget = budget_of();
                    Result r;
                    std::vector<Ctx> got = released_by(
                            engine, iupac(text),
                            request_of(Mode::ALL_OR_COUNT, 1'000'000'000, scope, strands),
                            budget, &r);
                    ASSERT_FALSE(r.refusal);
                    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
                    EXPECT_EQ(expected.size(), r.contexts->total.value) << text;
                    EXPECT_EQ(expected, got) << text << " mode " << mode;
                    std::map<uint32_t, uint64_t> by_offset;
                    for (const Ctx &c : expected) {
                        ++by_offset[c.offset];
                    }
                    for (const auto &[p, count] : r.contexts->by_offset) {
                        EXPECT_EQ(by_offset[p], count.value) << text << " offset " << p;
                    }
#if !_DNA5_GRAPH
                    // an 8-N run branched four ways per position costs 4^8 = 65536 range
                    // evaluations at least in the orientation where it leads; its core costs
                    // a few thousand at most. Not on $ACGTN: a pattern N admits A, C, G, T
                    // but a flank also N, so the run is searched there (stated, SPEC §7.8;
                    // the first DNA5 run, fix1-integrator, showed ~1.1e5 steps)
                    EXPECT_LT(r.work.steps, 8'000u) << text << " mode " << mode
                                                   << " strands " << to_string(strands);
#endif
                }
            }
        }
        // a long pattern whose anchor window starts with an N run (anchors counted)
        const std::string longp = "NNNNNN" + kmers[50].second.substr(0, 9) + "ACGT";
        ASSERT_GT(longp.size(), k);
        std::vector<Ctx> anchors = oracle(kmers, k, longp, Scope::ANY_OFFSET, Strands::FORWARD);
        Budget budget = budget_of();
        Result r = engine.count(iupac(longp), request_of(Mode::COUNT, 0, Scope::ANY_OFFSET,
                                                         Strands::FORWARD), budget);
        ASSERT_FALSE(r.refusal);
        EXPECT_EQ(anchors.size(), r.anchors->total.value);
        EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
#if !_DNA5_GRAPH
        EXPECT_LT(r.work.steps, 2'000u);
#endif
    }
}


TEST(PatternSearchFixes, CheapOrientationFirst) {
    // GG + N8 + TTGGCGATCT: wide in its forward orientation (the run after two bases), narrow
    // in its reverse one (the run after ten): the reverse search runs first, so a budget that
    // fits it leaves it exact (it was not started: the forward search spent everything)
    const size_t k = 21;
    std::vector<std::string> records = random_records(8, 2500, 31);
    const std::string text = "GGNNNNNNNNTTGGCGATCT";
    records.push_back("ACGG" + std::string("ACGTTGCA") + "TTGGCGATCTAC");
    records.push_back(own_rc("TTGG" + std::string("CATGCATG") + "TTGGCGATCTGA"));
    auto graph = build(k, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);

    Budget reverse_budget = budget_of();
    Result reverse = engine.count(iupac(text), request_of(Mode::COUNT, 0, Scope::ANY_OFFSET,
                                                          Strands::REVERSE), reverse_budget);
    Budget forward_budget = budget_of();
    Result forward = engine.count(iupac(text), request_of(Mode::COUNT, 0, Scope::ANY_OFFSET,
                                                          Strands::FORWARD), forward_budget);
    ASSERT_LT(10 * reverse.work.steps, forward.work.steps);
    ASSERT_GT(reverse.contexts->total.value, 0u);

    Budget budget = budget_of(reverse.work.steps);
    Result both = engine.count(iupac(text), request_of(Mode::COUNT, 0), budget);
    ASSERT_TRUE(both.stop);
    const Count &rev = both.contexts->by_orientation.at(Orientation::REVERSE);
    EXPECT_EQ(Relation::EXACT, rev.relation);
    EXPECT_EQ(reverse.contexts->total.value, rev.value);
    EXPECT_NE(Relation::EXACT, both.contexts->by_orientation.at(Orientation::FORWARD).relation);
    // the answer still lists the strands in plan order
    EXPECT_EQ((std::vector<Orientation> { Orientation::FORWARD, Orientation::REVERSE }),
              both.searched);

    // an exact pattern keeps the plan's order: forward first
    const std::string exact = "TTGGCGATCTACGTTG";
    Budget exact_forward = budget_of();
    Result ef = engine.count(iupac(exact), request_of(Mode::COUNT, 0, Scope::ANY_OFFSET,
                                                      Strands::FORWARD), exact_forward);
    Budget exact_budget = budget_of(ef.work.steps);
    Result eb = engine.count(iupac(exact), request_of(Mode::COUNT, 0), exact_budget);
    EXPECT_EQ(Relation::EXACT, eb.contexts->by_orientation.at(Orientation::FORWARD).relation);
}


// ---------------------------------------------------------------- X-GUARANTEES-01: the
// floor reads every searched orientation's anchor window

TEST(PatternSearchFixes, LongFloorPerSearchedWindow) {
    const size_t k = 11;
    auto graph = build(k, random_records(4, 400, 5), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const std::string p = "ACGTTGCAAGG";  // 11 bases: 22 bits
    const std::string n(k, 'N');
    Request request = request_of(Mode::COUNT, 0);
    request.min_information_bits = 20;

    auto run = [&](const std::string &text, Strands strands) {
        Request r = request;
        r.strands = strands;
        Budget budget = budget_of();
        return engine.count(iupac(text), r, budget);
    };

    // P + N^k: the reverse orientation's window is rc(N^k), 0 bits: refused when searched
    for (Strands strands : { Strands::BOTH, Strands::REVERSE }) {
        Result r = run(p + n, strands);
        ASSERT_TRUE(r.refusal) << to_string(strands);
        EXPECT_EQ("information_below_floor", r.refusal->code);
        EXPECT_NE(std::string::npos, r.refusal->message.find("reverse orientation"));
        EXPECT_EQ(0u, r.work.steps);
        // anchor_information_bits keeps its version-1 meaning, P[0, k); the floor's operand,
        // the least searched window, is min_anchor_information_bits (the owner's decision of
        // 2026-10-07)
        ASSERT_TRUE(r.anchor_information_bits);
        EXPECT_EQ(22.0, *r.anchor_information_bits);
        ASSERT_TRUE(r.min_anchor_information_bits);
        EXPECT_EQ(0.0, *r.min_anchor_information_bits);
    }
    Result forward = run(p + n, Strands::FORWARD);
    EXPECT_FALSE(forward.refusal);
    EXPECT_EQ(22.0, *forward.anchor_information_bits);
    EXPECT_EQ(22.0, *forward.min_anchor_information_bits);
    EXPECT_EQ(1u, forward.anchor_window_bits.size());

    // N^k + P: the reverse orientation's window is rc(P), 22 bits: admitted alone
    Result reverse = run(n + p, Strands::REVERSE);
    EXPECT_FALSE(reverse.refusal);
    EXPECT_EQ(0.0, *reverse.anchor_information_bits);
    EXPECT_EQ(22.0, *reverse.min_anchor_information_bits);
    for (Strands strands : { Strands::BOTH, Strands::FORWARD }) {
        Result r = run(n + p, strands);
        ASSERT_TRUE(r.refusal);
        // the forward window is the least informative: the message as before
        EXPECT_EQ("pattern: 0.0 information bits in the anchor window, below the floor of "
                  "20.0 for scope long", r.refusal->message);
        EXPECT_EQ(0.0, *r.anchor_information_bits);
        EXPECT_EQ(0.0, *r.min_anchor_information_bits);
    }

    // both windows informative: the least of the two, each stated
    // the last 11 positions: 10 bases and an N, 20 bits
    Result both = run(p + "A" + "CCGTAGGTAC" + "N", Strands::BOTH);
    ASSERT_FALSE(both.refusal) << both.refusal->message;
    ASSERT_EQ(2u, both.anchor_window_bits.size());
    EXPECT_EQ(22.0, both.anchor_window_bits.at(Orientation::FORWARD));
    EXPECT_EQ(20.0, both.anchor_window_bits.at(Orientation::REVERSE));
    EXPECT_EQ(22.0, *both.anchor_information_bits);
    EXPECT_EQ(20.0, *both.min_anchor_information_bits);

    // L <= k: neither is set
    Result shorter = run(p, Strands::BOTH);
    EXPECT_FALSE(shorter.anchor_information_bits);
    EXPECT_FALSE(shorter.min_anchor_information_bits);
}


// ---------------------------------------------------------------- E2-01, E2-02: even-k
// wrapped PRIMARY

const std::vector<std::string> kPalindromeRich {
    "TTACGCGTAAGGATCCTTA", "GAATATTCCGACGTACGTTG", "CCGGAATTCCATGCATGGC",
    "ACGCGTACGTAGCTAGCTTAAGCTT", "GGGCCCAATTGCAACGTTGCA", "TATATATAGCGCGCGCA"
};

TEST(PatternSearchFixes, EvenPrimaryPartialReleasesWhatTheCountCredits) {
    for (size_t k : { 4, 6 }) {
        auto graph = build(k, kPalindromeRich, DeBruijnGraph::PRIMARY);
        PatternSearch engine(*graph);
        for (const std::string &text
                : std::vector<std::string> { "NNN", "NA", "ACG", "N", "WS", "CG" }) {
            if (text.size() > k)
                continue;
            std::vector<Ctx> truth = oracle(*graph, text, Scope::ANY_OFFSET, Strands::BOTH);
            std::set<Ctx> genuine(truth.begin(), truth.end());
            Budget full_budget = budget_of();
            Result full = engine.count(iupac(text), request_of(Mode::COUNT, 0), full_budget);
            ASSERT_EQ(truth.size(), full.contexts->total.value) << text << " k " << k;

            for (uint64_t steps = 1; steps <= full.work.steps; ++steps) {
                Budget budget = budget_of(steps);
                Result r;
                std::vector<Ctx> got = released_by(engine, iupac(text),
                                                   request_of(Mode::PARTIAL, 1'000'000'000),
                                                   budget, &r);
                std::map<Orientation, uint64_t> per;
                for (const Ctx &c : got) {
                    ++per[c.orientation];
                    EXPECT_TRUE(genuine.count(c)) << text << " " << c;
                }
                EXPECT_TRUE(std::adjacent_find(got.begin(), got.end()) == got.end());
                for (const auto &[o, count] : r.contexts->by_orientation) {
                    // the list holds at least what the count credits each orientation
                    EXPECT_GE(per[o], count.value) << text << " k " << k << " steps " << steps
                                                   << " " << orientation_key(o);
                }
            }

            // the threshold stop: cut max_contexts means max_contexts returned
            for (uint64_t cap = 1; cap < truth.size(); cap += std::max<uint64_t>(1, truth.size() / 7)) {
                Request thr = request_of(Mode::PARTIAL, cap);
                thr.stop_at_threshold = true;
                Budget budget = budget_of();
                Result r;
                std::vector<Ctx> got = released_by(engine, iupac(text), thr, budget, &r);
                if (r.extraction->cut == StopReason::MAX_CONTEXTS) {
                    EXPECT_EQ(cap, got.size()) << text << " k " << k << " cap " << cap;
                }
            }
        }
    }
}

TEST(PatternSearchFixes, EvenPrimaryThresholdStopsOnTime) {
    // an exact pattern no offset of which can lie in a palindromic k-mer: the two base
    // searches find disjoint contexts, so the running bound is their sum (was the larger
    // of the two, about half: the stop never fired between half and the whole count)
    const size_t k = 10;
    auto graph = build(k, random_records(6, 15'000, 77), DeBruijnGraph::PRIMARY);
    PatternSearch engine(*graph);
    const std::string text = "AACCTG";
    for (Strands strands : { Strands::FORWARD, Strands::BOTH }) {
        const uint64_t truth = oracle(*graph, text, Scope::ANY_OFFSET, strands).size();
        ASSERT_GT(truth, 50u);
        // (not up to truth - 1: the running bound lags the count by the masked edges not
        // yet scanned, stated)
        for (uint64_t cap : { truth * 6 / 10, truth * 8 / 10, truth * 9 / 10 }) {
            Request request = request_of(Mode::COUNT, cap, Scope::ANY_OFFSET, strands);
            request.stop_at_threshold = true;
            Budget budget = budget_of();
            Result r = engine.count(iupac(text), request, budget);
            ASSERT_TRUE(r.stop) << "cap " << cap << " of " << truth;
            EXPECT_EQ(StopReason::MAX_CONTEXTS, r.stop->reason);
            EXPECT_EQ(Relation::AT_LEAST, r.contexts->total.relation);
            EXPECT_GT(r.contexts->total.value, cap);
            EXPECT_LE(r.contexts->total.value, truth);
        }
    }
}


// ---------------------------------------------------------------- T1-02, E2-04, E1-02, E3-02:
// a release that disagrees with its exact count is never published as complete

TEST(PatternSearchFixes, ReleaseDisagreeingWithTheCountThrows) {
    // a masked graph extended in place: DBGSuccinct::add_sequence marks every inserted edge
    // valid, the new sink dummy (W = $) included, which breaks the discount of
    // count_edges_with_last_symbol (the mask is trusted): GA in suffix scope counts 2, the
    // release finds 1
    auto graph = build(4, { "GAAAT", "CCTGA" }, DeBruijnGraph::BASIC);
    auto &dbg = const_cast<DBGSuccinct&>(base_of(*graph));
    dbg.add_sequence("CCTTG", [](node_index) {});
    PatternSearch engine(*graph);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, "GA");
    Request request = request_of(Mode::ALL_OR_COUNT, 1000, Scope::SUFFIX, Strands::FORWARD);

    Budget count_budget = budget_of();
    Result counted = engine.count(pattern, request, count_budget);
    ASSERT_EQ(Relation::EXACT, counted.contexts->total.relation);
    // partial, its cap above the count, releases fewer than the exact count too: it fails the
    // same way rather than state the short list as cut at max_contexts (owner decision #9 of
    // 2026-10-07, also in partial)
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    Budget partial_budget = budget_of();
    Result listed;
    EXPECT_THROW(released_by(engine, pattern, partial, partial_budget, &listed), std::logic_error);

    Budget budget = budget_of();
    std::vector<Ctx> received;
    EXPECT_THROW(engine.enumerate(pattern, request, budget, [&](const Context &c) {
        received.push_back(Ctx { c.node, c.offset, c.orientation });
    }), std::logic_error);
    // nothing reached the caller
    EXPECT_TRUE(received.empty());
}


// ---------------------------------------------------------------- pins of stated behaviour
// (E1-01, C1-03, X-DETERMINISM-01, T1-07)

TEST(PatternSearchFixes, ThresholdIsNotCheckedInTheMaskScans) {
    // E1-01: the deferred scans do not consult stop_at_threshold; the pattern can answer
    // EXACT above its threshold with no stop. Records whose last bases are mostly T: the
    // W-rule leaf of T has as many dummy-source edges as candidates, so discovery's lower
    // bound stays below the threshold and the scan settles the count
    // (the review's t2_marked graph, k = 5, built as the CLI builds it)
    std::vector<std::string> records { "AACGT", "CACGT", "GACGT", "TACGT", "AACGTT", "ACCGTT",
                                       "AGCGTT", "ATCGTT", "GACGTT", "GCCGTT", "GGCGTT",
                                       "GTCGTT", "CGTAAC", "GTACGTAC", "TTTTTTT",
                                       "ACGTACGTA" };
    auto graph = mtg::test::build_graph_batch<DBGSuccinct>(5, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const std::string text = "T";
    const uint64_t truth = oracle(*graph, text, Scope::SUFFIX, Strands::FORWARD).size();
    Request request = request_of(Mode::ALL_OR_COUNT, 1, Scope::SUFFIX, Strands::FORWARD);
    request.stop_at_threshold = true;
    Budget budget = budget_of();
    request.max_contexts = 3;
    Result r = engine.enumerate(iupac(text), request, budget, [](const Context &) {});
    // the scan ran to its end: exact, above the threshold, no stop
    ASSERT_GT(r.work.mask_scans, 0u);
    EXPECT_GT(truth, 3u);
    EXPECT_FALSE(r.stop);
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
    EXPECT_EQ(truth, r.contexts->total.value);
    EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, r.extraction->withheld);
}

TEST(PatternSearchFixes, ReleaseTimeStopAfterAStepStopKeepsTheFirstStop) {
    // C1-03, X-DETERMINISM-01: PARTIAL, discovery stopped by max_steps, then the work time
    // passes before the release: stop names the first stop, the time stop shows in cut and
    // determinism; the next pattern carries the sticky max_steps stop with cut time
    auto graph = medium_graph();
    PatternSearch engine(*graph);
    const Pattern pattern = iupac("AC");
    const uint64_t steps = 200;
    Budget untimed = budget_of(steps);
    Result reference = engine.count(pattern, request_of(Mode::COUNT, 0), untimed);
    ASSERT_TRUE(reference.stop);
    ASSERT_EQ(StopReason::MAX_STEPS, reference.stop->reason);

    Request partial = request_of(Mode::PARTIAL, 1'000'000);
    // discovery stops at max_steps having read the clock once (the pattern's start; no
    // stride is crossed below 4096 steps); the clock is late from the release's reading on
    auto readings = std::make_shared<int>(0);
    const Clock::time_point start = Clock::now();
    Budget jumping(steps, Deadline(start, 1000, 250, [start, readings]() {
        return start + std::chrono::seconds(++*readings > 1 ? 100 : 0);
    }));
    std::vector<Ctx> none;
    Result first = engine.enumerate(pattern, partial, jumping, [&](const Context &c) {
        none.push_back(Ctx { c.node, c.offset, c.orientation });
    });
    ASSERT_TRUE(first.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, first.stop->phase);
    EXPECT_EQ(StopReason::MAX_STEPS, first.stop->reason);
    EXPECT_EQ(StopReason::TIME, first.extraction->cut);
    EXPECT_TRUE(first.time_limited);
    EXPECT_TRUE(none.empty());
    EXPECT_EQ(reference.contexts->total.value, first.contexts->total.value);
    EXPECT_EQ(reference.contexts->total.relation, first.contexts->total.relation);

    Result second = engine.enumerate(pattern, partial, jumping, [&](const Context &) {});
    ASSERT_TRUE(second.stop);
    EXPECT_EQ(StopReason::MAX_STEPS, second.stop->reason);
    EXPECT_EQ(Relation::UNKNOWN, second.contexts->total.relation);
    EXPECT_EQ(StopReason::TIME, second.extraction->cut);
    EXPECT_TRUE(second.time_limited);
}

TEST(PatternSearchFixes, StrandHaltedAtItsFirstStepIsAtLeastZero) {
    // T1-07: the stop that refuses the - search's very first step leaves it AT_LEAST 0 (its
    // discovery was entered); one step earlier it never starts: UNKNOWN
    std::vector<std::string> records { "ACGTTGCAAGGCTTACGATCGATCGGGATTACA" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const Pattern pattern = Pattern::parse(PatternKind::DNA, "GAT");
    Budget fwd = budget_of();
    Result forward = engine.count(pattern, request_of(Mode::COUNT, 0, Scope::ANY_OFFSET,
                                                      Strands::FORWARD), fwd);
    const uint64_t forward_steps = forward.work.steps;

    Budget at = budget_of(forward_steps);
    Result r = engine.count(pattern, request_of(Mode::COUNT, 0), at);
    EXPECT_EQ(Relation::EXACT, r.contexts->by_orientation.at(Orientation::FORWARD).relation);
    EXPECT_EQ(Relation::AT_LEAST, r.contexts->by_orientation.at(Orientation::REVERSE).relation);
    EXPECT_EQ(0u, r.contexts->by_orientation.at(Orientation::REVERSE).value);

    Budget before = budget_of(forward_steps - 1);
    r = engine.count(pattern, request_of(Mode::COUNT, 0), before);
    EXPECT_EQ(Relation::AT_LEAST, r.contexts->by_orientation.at(Orientation::FORWARD).relation);
    EXPECT_EQ(Relation::UNKNOWN, r.contexts->by_orientation.at(Orientation::REVERSE).relation);
}


// ---------------------------------------------------------------- T1-10, C1-01: the pattern
// algebra against this file's own tables

TEST(PatternSearchFixes, ReverseComplementAndPalindromesExhaustive) {
    const std::string letters = "ACGTRYSWKMBDHVN";
    std::vector<std::string> texts { "" };
    for (size_t length = 1; length <= 4; ++length) {
        std::vector<std::string> next;
        for (const std::string &t : texts) {
            for (char c : letters) {
                next.push_back(t + c);
            }
        }
        for (const std::string &t : next) {
            Pattern p = iupac(t);
            const std::string rc = own_rc(t);
            ASSERT_EQ(rc, p.reverse_complement().text()) << t;
            ASSERT_EQ(rc == t, p.is_palindromic()) << t;
        }
        texts = std::move(next);
    }
    EXPECT_FALSE(iupac("A").is_palindromic());
    for (const char *t : { "S", "W", "N" }) {
        EXPECT_TRUE(iupac(t).is_palindromic()) << t;
    }
    for (const char *t : { "AN", "CAG", "TTAGGA" }) {
        EXPECT_FALSE(iupac(t).is_palindromic()) << t;
    }
}

TEST(PatternSearchFixes, InformationBitsUnchangedByTheTable) {
    std::mt19937 rng(3);
    const std::string letters = "ACGTRYSWKMBDHVN";
    for (int i = 0; i < 200; ++i) {
        std::string t;
        for (size_t j = 0, n = 1 + rng() % 300; j < n; ++j) {
            t.push_back(letters[rng() % letters.size()]);
        }
        double expected = 0;
        for (char c : t) {
            expected += std::log2(4.0 / static_cast<uint32_t>(__builtin_popcount(own_set(c))));
        }
        // the same double, bit for bit
        EXPECT_EQ(expected, iupac(t).information_bits()) << t;
    }
}


// ---------------------------------------------------------------- X-TESTS-04: larger even k

TEST(PatternSearchFixes, LargeEvenKPrimaryAgainstTheOracle) {
    for (size_t k : { 32, 64 }) {
        // records holding palindromic k-mers: a half and its reverse complement
        std::vector<std::string> records = random_records(3, 3 * k, 1000 + k);
        std::string half = random_records(1, k / 2, 2000 + k)[0];
        records.push_back("ACG" + half + own_rc(half) + "TTG");
        records.push_back(half + own_rc(half));
        auto graph = build(k, records, DeBruijnGraph::PRIMARY);
        PatternSearch engine(*graph);
        const std::string mid = half.substr(half.size() - 3) + own_rc(half).substr(0, 3);
        for (const std::string &text : { mid, std::string("ACG"), half.substr(0, 5),
                                        half + own_rc(half).substr(0, 2), std::string("WS") }) {
            if (text.size() > k)
                continue;
            std::vector<Ctx> expected = oracle(*graph, text, Scope::ANY_OFFSET, Strands::BOTH);
            Budget budget = budget_of();
            Result r;
            std::vector<Ctx> got = released_by(engine, iupac(text),
                                               request_of(Mode::ALL_OR_COUNT, 1'000'000'000),
                                               budget, &r);
            ASSERT_EQ(expected.size(), r.contexts->total.value) << text << " k " << k;
            ASSERT_EQ(expected, got) << text << " k " << k;

            // interrupted releases: genuine, no duplicates
            std::set<Ctx> genuine(expected.begin(), expected.end());
            for (uint64_t steps = 1; steps < r.work.steps; steps += 1 + r.work.steps / 40) {
                Budget b = budget_of(steps);
                Result s;
                std::vector<Ctx> partial = released_by(
                        engine, iupac(text), request_of(Mode::PARTIAL, 1'000'000'000), b, &s);
                EXPECT_TRUE(std::adjacent_find(partial.begin(), partial.end()) == partial.end());
                for (const Ctx &c : partial) {
                    EXPECT_TRUE(genuine.count(c)) << text << " " << c;
                }
                EXPECT_GE(partial.size(), s.contexts->total.value);
            }
        }
    }
}


// ---------------------------------------------------------------- M1-02: the sink pass
// sets the same bits

TEST(PatternSearchFixes, SinkPassMarksTheSameEdges) {
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        auto graph = build(9, random_records(40, 120, 9), mode);
        const boss::BOSS &boss = base_of(*graph).get_boss();
        sdsl::bit_vector expected(boss.num_edges() + 1, 0);
        // the former loop: every edge from 2 on with W = $
        for (boss::BOSS::edge_index i = 2; i <= boss.num_edges(); ++i) {
            if (!boss.get_W(i))
                expected[i] = true;
        }
        sdsl::bit_vector marked(boss.num_edges() + 1, 0);
        EXPECT_EQ(boss.rank_W(boss.num_edges(), 0) - 1, boss.mark_sink_dummy_edges(&marked));
        EXPECT_EQ(expected, marked);
        EXPECT_GT(sdsl::util::cnt_one_bits(marked), 30u);
    }
}


// ---------------------------------------------------------------- GPT review 3, item 2: the
// low-complexity note, an optional diagnostic, under the clock and never after a stop

// the oracle: sdust over the whole text, as the seeder calls it (T = 20, W = 64)
bool sdust_flags(const std::string &text) {
    int n = 0;
    uint64_t *intervals = sdust(0, reinterpret_cast<const uint8_t*>(text.data()),
                                static_cast<int>(text.size()), 20, 64, &n);
    std::free(intervals);
    return n > 0;
}

bool noted_low_complexity(const Result &r) {
    return std::find(r.notes.begin(), r.notes.end(), std::string(kNoteLowComplexity))
            != r.notes.end();
}

// a budget whose clock counts its readings in |*readings| and reads late from the
// |on_time| + 1-th on
Budget counted_budget(std::shared_ptr<uint64_t> readings,
                      uint64_t on_time = std::numeric_limits<uint64_t>::max(),
                      uint64_t max_steps = kSteps) {
    const Clock::time_point start = Clock::now();
    return Budget(max_steps, Deadline(start, 1000, 250, [start, readings, on_time]() {
        return start + std::chrono::seconds(++*readings > on_time ? 100 : 0);
    }));
}

TEST(PatternSearchFixes, LowComplexityNoteAsSdustOverTheWholePattern) {
    // the engine reads a pattern in pieces of 128 bases (plus the 63 of the window before the
    // next): the note must be stated exactly when sdust flags the whole text. Random texts of
    // up to 1,500 bases, two in three with a short repeat (8 to 40 bases), half of these
    // across a multiple of 128, where a piece without the overlap would see only part of it
    auto graph = build(13, random_records(2, 300, 5), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    std::mt19937 rng(23);
    const std::vector<std::string> units { "A", "AT", "CA", "ACG", "GGT", "ATG" };
    uint64_t flagged = 0;
    uint64_t clean = 0;
    for (int t = 0; t < 300; ++t) {
        const size_t length = rng() % 2
            ? 1 + rng() % 1500
            : 128 * (1 + rng() % 8) + 64 - rng() % 128;
        std::string text;
        for (size_t i = 0; i < length; ++i) {
            text.push_back("ACGT"[rng() % 4]);
        }
        if (rng() % 3) {
            const std::string &unit = units[rng() % units.size()];
            const size_t run = 8 + rng() % 33;
            const size_t at = rng() % 2 && length > 128 + run
                ? 128 * (1 + rng() % (length / 128)) - run / 2
                : rng() % length;
            for (size_t i = 0; i < run && at + i < length; ++i) {
                text[at + i] = unit[i % unit.size()];
            }
        }
        Budget budget = budget_of();
        Result r = engine.count(iupac(text), request_of(Mode::COUNT, 1000), budget);
        ASSERT_FALSE(r.refusal) << text;
        ASSERT_FALSE(r.stop) << text;
        EXPECT_FALSE(r.time_limited) << text;
        const bool expected = sdust_flags(text);
        EXPECT_EQ(expected, noted_low_complexity(r)) << text;
        ++(expected ? flagged : clean);
    }
    EXPECT_LT(60u, flagged);
    EXPECT_LT(60u, clean);
}

TEST(PatternSearchFixes, LowComplexityNoteNotStatedAfterAStop) {
    // the review's repro: 10,000 Met (ATG x 10,000, 0.7 s of sdust read whole) with
    // max_steps 1 answered 503 after its 1 ms of work time. After any stop the diagnostic is
    // not run: no note, and no clock reading for it
    auto graph = build(31, random_records(2, 300, 5), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const Pattern met = Pattern::parse(PatternKind::PROTEIN, std::string(10'000, 'M'));
    ASSERT_TRUE(met.is_exact());
    Request request = request_of(Mode::COUNT, 1000);
    request.extend_paths = true;

    auto readings = std::make_shared<uint64_t>(0);
    Budget one = counted_budget(readings, std::numeric_limits<uint64_t>::max(), 1);
    Result stopped = engine.count(met, request, one);
    ASSERT_TRUE(stopped.stop);
    EXPECT_EQ(StopPhase::DISCOVERY, stopped.stop->phase);
    EXPECT_EQ(StopReason::MAX_STEPS, stopped.stop->reason);
    EXPECT_FALSE(noted_low_complexity(stopped));
    EXPECT_EQ(1u, *readings);  // the pattern's start only
    // the next pattern of the request, stopped before it starts
    Result next = engine.count(Pattern::parse(PatternKind::DNA, std::string(60, 'A')), request,
                               one);
    ASSERT_TRUE(next.stop);
    EXPECT_FALSE(noted_low_complexity(next));

    // completed: the note, from the first piece, which reads no clock (a twin whose bases
    // differ only after the anchor windows, and is not exact, reads it as often: every
    // further piece would read it)
    Budget completed = counted_budget(readings);
    *readings = 0;
    Result full = engine.count(met, request, completed);
    EXPECT_FALSE(full.stop);
    EXPECT_FALSE(full.time_limited);
    EXPECT_TRUE(noted_low_complexity(full));
    const uint64_t with_diagnostic = *readings;
    std::string inexact(10'000, 'M');
    inexact[5'000] = 'X';
    *readings = 0;
    Budget twin_budget = counted_budget(readings);
    Result twin = engine.count(Pattern::parse(PatternKind::PROTEIN, inexact), request,
                               twin_budget);
    EXPECT_FALSE(twin.stop);
    EXPECT_FALSE(noted_low_complexity(twin));
    EXPECT_EQ(twin.work.steps, full.work.steps);
    EXPECT_EQ(*readings, with_diagnostic);
}

TEST(PatternSearchFixes, LowComplexityDiagnosticReadsTheClock) {
    // a text sdust flags only at its end: 3,000 random bases it does not flag (the first seed
    // that gives such), then 40 A. Complete, the diagnostic reads the clock before each of its
    // pieces but the first; late from its third reading on, it ends there: the counts
    // complete, no stop, the note left out and the answer time_limited
    std::string text;
    for (uint32_t seed = 1; text.empty(); ++seed) {
        ASSERT_GT(100u, seed);
        std::string candidate = random_records(1, 3'000, seed).front();
        if (!sdust_flags(candidate))
            text = candidate;
    }
    text += std::string(40, 'A');
    ASSERT_TRUE(sdust_flags(text));
    std::string inexact = text;
    inexact[1'500] = 'N';

    auto graph = build(31, random_records(2, 300, 5), DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    const Request request = request_of(Mode::COUNT, 1000);
    auto readings = std::make_shared<uint64_t>(0);
    Budget twin_budget = counted_budget(readings);
    Result twin = engine.count(iupac(inexact), request, twin_budget);
    const uint64_t before = *readings;
    *readings = 0;
    Budget full_budget = counted_budget(readings);
    Result full = engine.count(iupac(text), request, full_budget);
    EXPECT_TRUE(noted_low_complexity(full));
    EXPECT_FALSE(full.time_limited);
    EXPECT_EQ(twin.work.steps, full.work.steps);
    // the 22 pieces inside the 3,000 random bases are not flagged (a piece is flagged only
    // where its text is): the 23rd or the 24th, which ends the text, is the last read
    EXPECT_LE(before + 22, *readings);
    EXPECT_GE(before + 23, *readings);

    *readings = 0;
    Budget late = counted_budget(readings, before + 2);
    Result cut = engine.count(iupac(text), request, late);
    EXPECT_FALSE(cut.stop);
    EXPECT_TRUE(cut.time_limited);
    EXPECT_FALSE(noted_low_complexity(cut));
    EXPECT_EQ(before + 3, *readings);
    ASSERT_TRUE(cut.anchors);
    EXPECT_EQ(Relation::EXACT, cut.anchors->total.relation);
    EXPECT_EQ(full.anchors->total.value, cut.anchors->total.value);
    EXPECT_EQ(full.work.steps, cut.work.steps);
}


// ---------------------------------------------------------------- GPT review 3, item 5: the
// extension reads the clock before every anchor, every 64 anchors it lists and every 64
// nodes its DFS expands

TEST(PatternSearchFixes, ExtensionReadsTheClockBeforeEveryAnchor) {
    // staging (refseq33m, count, long_search paths, a 40-base pattern half N): 674 anchors
    // listed, spelled and extended on a cold index ran 3 s past the work time without a clock
    // reading (their DFS charged 1,129 steps, within one stride), and the answer was 503
    // instead of a stated stop. Here every k-mer is an anchor (a window of N) and has an
    // outgoing k-mer (circular records), and the clock runs with the steps: on time through
    // discovery and the listing, late as soon as the extension charged an edge. Discovery is
    // made to end one step past a stride, so that no stride crossing reads the clock after it
    const size_t k = 11;
    std::vector<std::string> records = random_records(3, 400, 31);
    for (std::string &record : records) {
        record += record.substr(0, k - 1);
    }
    auto graph = build(k, records, DeBruijnGraph::BASIC);
    PatternSearch engine(*graph);
    Request request = request_of(Mode::COUNT, 1'000'000'000, Scope::ANY_OFFSET,
                                 Strands::FORWARD);
    request.extend_paths = true;
    request.max_paths = 1'000'000'000;

    // |tail| bases after the window: 10 (an exact tail, most anchors end at its first base)
    // or 200 (N: every anchor's DFS follows its record for 200 nodes)
    for (const std::string &tail : { std::string("ACGTACGTAC"), std::string(200, 'N') }) {
        SCOPED_TRACE(tail.size());
        const Pattern pattern = iupac(std::string(k, 'N') + tail);
        Budget unbounded = budget_of();
        const Result full = engine.count(pattern, request, unbounded);
        ASSERT_EQ(Extension::COMPLETED, full.anchors->extension);
        const uint64_t anchors = full.anchors->total.value;
        ASSERT_LT(2 * Budget::kReleaseClockStride, anchors);
        EXPECT_EQ(anchors, full.work.extension_anchors);
        ASSERT_LE(anchors, full.work.extension_edges);
        const uint64_t discovery = full.work.steps - full.work.extension_edges;

        const uint64_t stride = Budget::kClockStride;
        const uint64_t ahead = (stride + 1 - discovery % stride) % stride;
        const uint64_t at_end = ahead + discovery;
        auto budget_of_clock = std::make_shared<const Budget*>(nullptr);
        auto readings_at_end = std::make_shared<uint64_t>(0);
        const Clock::time_point start = Clock::now();
        Budget budget(kSteps, Deadline(start, 1000, 250,
                                       [start, budget_of_clock, readings_at_end, at_end]() {
            const uint64_t steps = (*budget_of_clock)->steps_used();
            *readings_at_end += steps == at_end;
            return start + std::chrono::seconds(steps > at_end ? 100 : 0);
        }));
        *budget_of_clock = &budget;
        ASSERT_TRUE(budget.charge(ahead));
        Result r = engine.count(pattern, request, budget);

        // a stated stop: the anchors exact, the paths a lower bound
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopPhase::EXTENSION, r.stop->phase);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_TRUE(r.time_limited);
        EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
        EXPECT_EQ(anchors, r.anchors->total.value);
        EXPECT_EQ(Extension::STOPPED, r.anchors->extension);
        EXPECT_EQ(Relation::AT_LEAST, r.anchors->paths.relation);
        EXPECT_LE(r.anchors->paths.value, full.anchors->paths.value);
        EXPECT_EQ(discovery, r.work.steps - r.work.extension_edges);
        // the first anchor only, and of its DFS (200 nodes deep for the N tail) at most the
        // nodes before the 64th: each entered as a candidate before it is expanded
        EXPECT_EQ(1u, r.work.extension_anchors);
        EXPECT_LE(1u, r.work.extension_edges);
        EXPECT_GT(Budget::kReleaseClockStride, r.anchors->candidates_examined);
        // on time: the listing's start, its reading every 64 anchors, the first anchor's
        EXPECT_EQ(2 + (anchors - 1) / Budget::kReleaseClockStride, *readings_at_end);
    }
}

#endif // ! _PROTEIN_GRAPH

} // namespace
