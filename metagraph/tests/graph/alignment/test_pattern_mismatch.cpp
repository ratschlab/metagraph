/**
 * Bounded mismatches in the pattern search (SPEC §21; the mismatch plan's increments 5m-0 and
 * 5m-1): DNA and IUPAC patterns of L <= k searched with Request::max_mismatches 1 to 3, checked
 * against the oracles of pattern_mismatch_oracles.hpp, which read the pattern through their own
 * IUPAC table and share only the graph with the engine:
 *  - the Hamming graph-walk (O1) gives every context within m with its level and substituted
 *    positions, in answer order: every count per level, offset and orientation, every release
 *    and its order are checked against it;
 *  - the Hamming record-scan (O2) checks O1 against the records, without the graph's ids;
 *  - the candidate enumeration (O3) checks O1 by a different algorithm: every string of the
 *    ball searched exactly by the engine without mismatches;
 *  - the brute-force ball (O5) checks the information bits left after m substitutions.
 * Requests without mismatches are the other suites' (their answers must not change); here the
 * pass of level 0 is checked to be that search, step for step.
 *
 * Alphabet: DNA (DNA4) builds only; the engine's DNA5 path is not served (alphabet_untested).
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../all/test_dbg_helpers.hpp"
#include "pattern_mismatch_oracles.hpp"
#include "pattern_test_support.hpp"

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

#if _DNA_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using namespace mtg::test::mismatch;
using mtg::test::build_graph;
using mtg::test::unbounded_deadline;

typedef DeBruijnGraph::node_index node_index;
typedef Deadline::Clock Clock;

constexpr uint64_t kManySteps = 1'000'000'000;


// ---------------------------------------------------------------- helpers

std::shared_ptr<DeBruijnGraph> build(size_t k, const std::vector<std::string> &records,
                                     DeBruijnGraph::Mode mode) {
    return build_graph<DBGSuccinct>(k, records, mode);
}

std::vector<std::string> random_records(std::mt19937 &rng, size_t num, size_t min_length,
                                        size_t max_length) {
    std::vector<std::string> records;
    for (size_t i = 0; i < num; ++i) {
        const size_t length = min_length + rng() % (max_length - min_length + 1);
        std::string r;
        for (size_t j = 0; j < length; ++j) {
            r.push_back("ACGT"[rng() % 4]);
        }
        records.push_back(r);
    }
    return records;
}

// a request with no information floor and no threshold in the way
Request make_request(Scope scope, Strands strands, uint8_t m, Mode mode = Mode::ALL_OR_COUNT) {
    Request request;
    request.mode = mode;
    request.scope = scope;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = 1'000'000'000;
    request.max_anchors = 1'000'000'000;
    request.max_mismatches = m;
    return request;
}

Budget budget_of(uint64_t max_steps = kManySteps) {
    return Budget(max_steps, unbounded_deadline());
}

Pattern parse(const std::string &text) {
    return Pattern::parse(text.find_first_not_of("ACGT") == std::string::npos
                                  ? PatternKind::DNA
                                  : PatternKind::IUPAC,
                          text);
}

Result count_of(const DeBruijnGraph &graph, const std::string &text, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget = budget_of(max_steps);
    return PatternSearch(graph).count(parse(text), request, budget);
}

/**
 * enumerate(): every released context as an MCtx, its level the engine's mismatches and its
 * positions those substitutions() reads off its instance, each substitution checked against
 * the oracles' table (expected is P's code, found the instance's base in P's orientation and
 * outside expected's bases).
 */
std::vector<MCtx> released_by(const DeBruijnGraph &graph, const std::string &text,
                              const Request &request, Budget &budget, Result *result) {
    const Pattern pattern = parse(text);
    const size_t L = text.size();
    std::vector<MCtx> released;
    *result = PatternSearch(graph).enumerate(pattern, request, budget, [&](const Context &c) {
        const std::string kmer = graph.get_node_sequence(c.node);
        const std::string instance = kmer.substr(c.offset, L);
        const std::vector<Substitution> subs = substitutions(pattern, c.orientation, instance);
        EXPECT_EQ(c.mismatches, subs.size()) << text << " " << kmer;
        MCtx m { c.node, c.offset, c.orientation, c.mismatches, {} };
        for (const Substitution &s : subs) {
            m.positions.push_back(s.position);
            EXPECT_EQ(text[s.position], s.expected);
            const char base = c.orientation == Orientation::REVERSE
                ? complement_base(instance[L - 1 - s.position])
                : instance[s.position];
            EXPECT_EQ(base, s.found);
            EXPECT_EQ(std::string::npos, iupac_bases(s.expected).find(s.found));
        }
        released.push_back(m);
    });
    return released;
}

std::vector<MCtx> released_by(const DeBruijnGraph &graph, const std::string &text,
                              const Request &request, Result *result,
                              uint64_t max_steps = kManySteps) {
    Budget budget = budget_of(max_steps);
    return released_by(graph, text, request, budget, result);
}

// |count| states a true relation to |truth| (§3)
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
        << count.upper << "] for the true " << truth;
}

// the oracle's contexts per level (zeros included), per offset and per orientation
struct Truth {
    std::vector<uint64_t> by_level;
    std::map<uint32_t, uint64_t> by_offset;
    std::map<Orientation, uint64_t> by_orientation;
    uint64_t total = 0;
};

Truth truth_of(const std::vector<MCtx> &expected, uint8_t m, size_t L, size_t k,
               const Request &request, const std::string &text) {
    Truth truth;
    truth.by_level.assign(m + 1, 0);
    for (uint32_t p : scope_offsets(L, k, request.scope)) {
        truth.by_offset[p] = 0;
    }
    for (const auto &[o, q] : orientations(text, request.strands)) {
        truth.by_orientation[o] = 0;
    }
    for (const MCtx &c : expected) {
        ++truth.by_level.at(c.level);
        ++truth.by_offset.at(c.offset);
        ++truth.by_orientation.at(c.orientation);
        ++truth.total;
    }
    return truth;
}

/**
 * Every check of one completed (graph, pattern, request) against O1 (and O2 with |records|):
 * every count EXACT per level, offset and orientation; the steps per pass summing to the
 * steps; enumerate() releasing O1's list in (level, node, offset, orientation) order with each
 * context's level and positions; ALL_OR_COUNT's threshold on the total; PARTIAL's prefixes.
 * Returns the number of contexts.
 */
size_t check_against_oracles(const DeBruijnGraph &graph, const std::string &text,
                             const Request &request,
                             const std::vector<std::string> *records = nullptr,
                             DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    const size_t k = graph.get_k();
    const size_t L = text.size();
    const uint8_t m = request.max_mismatches;
    const std::vector<MCtx> expected = hamming_walk(graph, text, request.scope,
                                                    request.strands, m);
    if (records) {
        // the graph holds exactly the records' k-mers (no cleaning in these builds)
        EXPECT_EQ(hamming_records(*records, k, mode, text, request.scope, request.strands, m),
                  spelled(graph, expected)) << text;
    }
    const Truth truth = truth_of(expected, m, L, k, request, text);

    Result counted = count_of(graph, text, request);
    EXPECT_FALSE(counted.refusal);
    if (counted.refusal)
        return 0;
    EXPECT_FALSE(counted.stop) << text;
    EXPECT_EQ(size_t(m) + 1, counted.contexts->by_mismatches.size());
    EXPECT_EQ(Relation::EXACT, counted.contexts->total.relation);
    EXPECT_EQ(truth.total, counted.contexts->total.value) << text << " m " << int(m);
    for (uint8_t t = 0; t <= m && t < counted.contexts->by_mismatches.size(); ++t) {
        const Count &c = counted.contexts->by_mismatches[t];
        EXPECT_EQ(Relation::EXACT, c.relation);
        EXPECT_EQ(truth.by_level[t], c.value) << text << " level " << int(t);
    }
    EXPECT_EQ(truth.by_offset.size(), counted.contexts->by_offset.size());
    for (const auto &[p, n] : truth.by_offset) {
        EXPECT_EQ(Relation::EXACT, counted.contexts->by_offset.at(p).relation);
        EXPECT_EQ(n, counted.contexts->by_offset.at(p).value) << text << " offset " << p;
    }
    EXPECT_EQ(truth.by_offset.at(k - L), counted.contexts->suffix.value);
    EXPECT_EQ(truth.by_orientation.size(), counted.contexts->by_orientation.size());
    for (const auto &[o, n] : truth.by_orientation) {
        EXPECT_EQ(n, counted.contexts->by_orientation.at(o).value) << orientation_key(o);
    }
    EXPECT_EQ(size_t(m) + 1, counted.work.steps_by_mismatches.size());
    EXPECT_EQ(counted.work.steps,
              std::accumulate(counted.work.steps_by_mismatches.begin(),
                              counted.work.steps_by_mismatches.end(), uint64_t(0)));
    EXPECT_TRUE(counted.mismatch_information_bits);

    // the release: O1's list, levels and positions included, in answer order
    Result all;
    const std::vector<MCtx> released = released_by(graph, text, request, &all);
    EXPECT_EQ(expected, released) << text << " m " << int(m);
    EXPECT_TRUE(all.extraction && all.extraction->complete);
    EXPECT_EQ(counted.work.steps, all.work.steps);
    EXPECT_EQ(counted.work.steps_by_mismatches, all.work.steps_by_mismatches);

    // ALL_OR_COUNT: all levels or none, on the total
    if (!expected.empty()) {
        Request tight = request;
        tight.max_contexts = expected.size() - 1;
        Result withheld;
        EXPECT_TRUE(released_by(graph, text, tight, &withheld).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
    }
    // PARTIAL: the first max_contexts of the closest levels
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()) }) {
        partial.max_contexts = cap;
        Result cut;
        const std::vector<MCtx> prefix = released_by(graph, text, partial, &cut);
        const size_t n = std::min<size_t>(cap, expected.size());
        EXPECT_EQ(std::vector<MCtx>(expected.begin(), expected.begin() + n), prefix)
            << text << " cap " << cap;
        EXPECT_EQ(cap >= expected.size(), cut.extraction->complete);
    }
    return expected.size();
}

// a random DNA or IUPAC pattern of |L| positions: IUPAC with every code, and N runs leading,
// internal or trailing
std::string random_pattern(std::mt19937 &rng, size_t L, bool iupac) {
    static const std::string kCodes = "RYSWKMBDHVN";
    std::string text;
    for (size_t i = 0; i < L; ++i) {
        text.push_back(iupac && rng() % 3 == 0 ? kCodes[rng() % kCodes.size()]
                                               : "ACGT"[rng() % 4]);
    }
    if (iupac && L >= 3 && rng() % 3 == 0) {
        // an N run at the start, inside or at the end
        const size_t run = 1 + rng() % std::max<size_t>(1, L / 3);
        const size_t where = rng() % 3;
        const size_t at = where == 0 ? 0 : where == 1 ? rng() % (L - run + 1) : L - run;
        for (size_t i = at; i < at + run; ++i) {
            text[i] = 'N';
        }
    }
    return text;
}

const char* mode_name(DeBruijnGraph::Mode mode) {
    switch (mode) {
        case DeBruijnGraph::BASIC: return "basic";
        case DeBruijnGraph::CANONICAL: return "canonical";
        case DeBruijnGraph::PRIMARY: return "primary";
    }
    return "?";
}

/**
 * One random case per seed: a graph of a few random records (k in {5, 6, 7, 9}, BASIC,
 * CANONICAL or wrapped PRIMARY, odd and even k), a random DNA or IUPAC pattern of 1 .. k
 * positions, m in {1, 2, 3}, a scope (any_offset only on PRIMARY) and strands.
 */
struct Case {
    size_t k;
    DeBruijnGraph::Mode mode;
    std::vector<std::string> records;
    std::shared_ptr<DeBruijnGraph> graph;
    std::string text;
    Request request;
};

Case random_case(uint32_t seed) {
    std::mt19937 rng(seed);
    Case c;
    c.k = std::vector<size_t> { 5, 6, 7, 9 }[rng() % 4];
    c.mode = std::vector<DeBruijnGraph::Mode> {
        DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[rng() % 3];
    c.records = random_records(rng, 2 + rng() % 4, 10, 60);
    if (rng() % 4 == 0) {
        // a record split by a symbol the build does not index, and a record shorter than k
        c.records[0][rng() % c.records[0].size()] = 'N';
        c.records.push_back(std::string("ACGT").substr(0, 1 + rng() % 3));
    }
    c.graph = build(c.k, c.records, c.mode);
    const size_t L = 1 + rng() % c.k;
    c.text = random_pattern(rng, L, rng() % 2);
    const Scope scope = c.mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0 ? Scope::SUFFIX
                                                                           : Scope::ANY_OFFSET;
    const Strands strands = std::vector<Strands> {
        Strands::BOTH, Strands::BOTH, Strands::FORWARD, Strands::REVERSE }[rng() % 4];
    c.request = make_request(scope, strands, 1 + rng() % 3);
    return c;
}

std::string describe(const Case &c) {
    return c.text + " k " + std::to_string(c.k) + " " + mode_name(c.mode) + " m "
         + std::to_string(c.request.max_mismatches) + " " + to_string(c.request.scope) + " "
         + to_string(c.request.strands);
}


// ---------------------------------------------------------------- 5m-0: the oracles

TEST(PatternMismatch, OraclesAgree) {
    // O1 = O2 on every case (the graph holds the records' k-mers); O1 = O3 where the ball is
    // small enough to enumerate; O1's level 0 = the search without mismatches
    size_t compared_o3 = 0;
    size_t with_substitutions = 0;
    for (uint32_t seed = 0; seed < 120; ++seed) {
        const Case c = random_case(seed);
        SCOPED_TRACE(describe(c));
        const uint8_t m = c.request.max_mismatches;
        const std::vector<MCtx> o1 = hamming_walk(*c.graph, c.text, c.request.scope,
                                                  c.request.strands, m);
        EXPECT_EQ(hamming_records(c.records, c.k, c.mode, c.text, c.request.scope,
                                  c.request.strands, m),
                  spelled(*c.graph, o1));
        if (c.text.size() <= 8 && m <= 2) {
            const PatternSearch engine(*c.graph);
            EXPECT_EQ(o1, candidate_enumeration(engine, c.text, c.request.scope,
                                                c.request.strands, m));
            ++compared_o3;
        }
        // level 0 is the search without mismatches
        Request exact = c.request;
        exact.max_mismatches = 0;
        Result r;
        const std::vector<MCtx> v1 = released_by(*c.graph, c.text, exact, &r);
        std::vector<MCtx> level0;
        for (const MCtx &x : o1) {
            if (!x.level)
                level0.push_back(x);
            with_substitutions += x.level > 0;
        }
        EXPECT_EQ(level0, v1);
    }
    EXPECT_LT(40u, compared_o3);
    EXPECT_LT(200u, with_substitutions);
}

TEST(PatternMismatch, OracleCountsTheNoteExample) {
    // the oracles' own arithmetic on a hand-made case: AACGTA and its reverse complement TACGTT
    // are each one substitution from the instance AACGTT (both orientations at one place)
    const std::vector<MCtx> both = hamming_walk(*build(6, { "AACGTT" }, DeBruijnGraph::BASIC),
                                                "AACGTA", Scope::ANY_OFFSET, Strands::BOTH, 1);
    ASSERT_EQ(2u, both.size());
    EXPECT_EQ(Orientation::FORWARD, both[0].orientation);
    EXPECT_EQ(Orientation::REVERSE, both[1].orientation);
    // P's last base A: T found (forward); rc(P)'s first base T, P's position 0 ... read back
    // in P's coordinates: the reverse instance AACGTT is rc(P) = TACGTT with A at 0, which is
    // P's position 5
    EXPECT_EQ(std::vector<uint32_t> { 5 }, both[0].positions);
    EXPECT_EQ(std::vector<uint32_t> { 5 }, both[1].positions);
    EXPECT_EQ(1u, both[0].level);
    EXPECT_EQ(1u, both[1].level);
    // the brute-force ball of an exact 3-mer: 1, 9, 27, 27 strings at 0 .. 3 substitutions
    EXPECT_EQ((std::vector<uint64_t> { 1, 9, 27, 27 }), ball_histogram("ACG"));
    // R admits 2 bases and 2 substitutes, N 4 and none
    EXPECT_EQ((std::vector<uint64_t> { 2, 2 }), ball_histogram("R"));
    EXPECT_EQ((std::vector<uint64_t> { 4, 0 }), ball_histogram("N"));
}


// ---------------------------------------------------------------- 5m-1: against the oracles

TEST(PatternMismatch, RandomCasesAgainstOracles) {
    // 300 random (graph, pattern, request): every count per level, offset and orientation,
    // the release with each context's level and positions, the thresholds
    size_t contexts = 0;
    std::vector<size_t> at_level(kMaxMismatches + 1, 0);
    for (uint32_t seed = 1000; seed < 1300; ++seed) {
        const Case c = random_case(seed);
        SCOPED_TRACE(describe(c));
        contexts += check_against_oracles(*c.graph, c.text, c.request, &c.records, c.mode);
        for (const MCtx &x : hamming_walk(*c.graph, c.text, c.request.scope, c.request.strands,
                                          c.request.max_mismatches)) {
            ++at_level[x.level];
        }
    }
    EXPECT_LT(2000u, contexts);
    // every level is released (level 3 survives the two bits a retained range keeps)
    for (size_t t = 0; t <= kMaxMismatches; ++t) {
        EXPECT_LT(50u, at_level[t]) << t;
    }
}


TEST(PatternMismatch, RandomCasesAtK31) {
    // the served k: records of a few hundred bases, patterns of 10 to 31 positions (a 31-mer
    // is the NOTE's k-mer neighbourhood: one offset), every mode, m 1 to 3
    size_t contexts = 0;
    size_t substituted = 0;
    for (uint32_t seed = 6000; seed < 6030; ++seed) {
        std::mt19937 rng(seed);
        const auto mode = std::vector<DeBruijnGraph::Mode> {
            DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[seed % 3];
        std::vector<std::string> records = random_records(rng, 3, 150, 300);
        // variants of a stretch of the first record, so that the neighbourhood is not empty
        const std::string stretch = records[0].substr(20, 60);
        for (int v = 0; v < 4; ++v) {
            std::string variant = stretch;
            for (int s = 0; s <= v % 4; ++s) {
                variant[rng() % variant.size()] = "ACGT"[rng() % 4];
            }
            records.push_back(variant);
        }
        auto graph = build(31, records, mode);
        const size_t L = 10 + rng() % 22;
        std::string text = records[0].substr(25 + rng() % 20, L);
        if (rng() % 2)
            text[rng() % L] = "RYSWKMBDHVN"[rng() % 11];
        const Request request = make_request(
                mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0 ? Scope::SUFFIX
                                                                 : Scope::ANY_OFFSET,
                Strands::BOTH, 1 + rng() % 3);
        SCOPED_TRACE(text + " " + mode_name(mode) + " m "
                     + std::to_string(request.max_mismatches));
        contexts += check_against_oracles(*graph, text, request, &records, mode);
        for (const MCtx &c : hamming_walk(*graph, text, request.scope, request.strands,
                                          request.max_mismatches)) {
            substituted += c.level > 0;
        }
    }
    EXPECT_LT(100u, contexts);
    EXPECT_LT(20u, substituted);
}

// ---------------------------------------------------------------- pass 0 and the search
// without mismatches

// the same graph built again, its dummy-edge mask removed (ids and BOSS identical)
std::shared_ptr<DeBruijnGraph> unmasked_copy(const Case &c) {
    auto graph = build(c.k, c.records, c.mode);
    const_cast<DBGSuccinct&>(stored_graph(*graph)).reset_mask();
    return graph;
}

TEST(PatternMismatch, PassZeroIsTheSearchWithoutMismatches) {
    // pass 0 charges the steps of the search without mismatches in the same order: under every
    // max_steps that stops that search, the pass stops at the same step, in the same phase,
    // with the same work, its level 0 counted as that search counted everything and its
    // PARTIAL release the same; under its full steps, pass 0 completes with that search's
    // count and steps. With and without the mask (the check of few unchecked candidates in
    // both)
    size_t stops = 0;
    for (uint32_t seed = 2000; seed < 2060; ++seed) {
        const Case c = random_case(seed);
        const bool unmasked = seed % 3 == 0;
        auto graph = unmasked ? unmasked_copy(c) : c.graph;
        SCOPED_TRACE(describe(c) + (unmasked ? " unmasked" : ""));
        Request with = c.request;
        with.max_checked_entries = unmasked ? 3 : 0;
        Request without = with;
        without.max_mismatches = 0;

        const Result v1 = count_of(*graph, c.text, without);
        const Result full = count_of(*graph, c.text, with);
        ASSERT_FALSE(v1.stop);
        ASSERT_FALSE(full.stop);
        const uint64_t steps = v1.work.steps;
        EXPECT_EQ(steps, full.work.steps_by_mismatches.at(0));
        const Count &level0 = full.contexts->by_mismatches.at(0);
        EXPECT_EQ(v1.contexts->total.relation, level0.relation);
        EXPECT_EQ(v1.contexts->total.value, level0.value);
        EXPECT_EQ(v1.contexts->total.upper, level0.upper);

        for (uint64_t s = 0; s <= steps; ++s) {
            const Result a = count_of(*graph, c.text, without, s);
            const Result b = count_of(*graph, c.text, with, s);
            ASSERT_EQ(bool(a.stop), s < steps);
            ASSERT_TRUE(b.stop || s == steps);
            if (s == steps) {
                // pass 0 completed; a stop, if any, came in a later pass
                if (b.stop) {
                    EXPECT_LT(0, *b.stop->mismatch_level);
                }
                EXPECT_EQ(steps, b.work.steps_by_mismatches.at(0));
                EXPECT_EQ(v1.contexts->total.value, b.contexts->by_mismatches.at(0).value);
                EXPECT_EQ(v1.contexts->total.relation, b.contexts->by_mismatches.at(0).relation);
                continue;
            }
            ++stops;
            ASSERT_TRUE(b.stop->mismatch_level);
            EXPECT_EQ(0, *b.stop->mismatch_level) << s;
            EXPECT_FALSE(a.stop->mismatch_level);
            EXPECT_EQ(a.stop->phase, b.stop->phase) << s;
            EXPECT_EQ(a.stop->reason, b.stop->reason) << s;
            EXPECT_EQ(a.work.ranges_visited, b.work.ranges_visited) << s;
            EXPECT_EQ(a.work.mask_scans, b.work.mask_scans) << s;
            EXPECT_EQ(a.work.steps, b.work.steps) << s;
            const Count &x = a.contexts->total;
            const Count &y = b.contexts->by_mismatches.at(0);
            EXPECT_EQ(x.relation, y.relation) << s;
            EXPECT_EQ(x.value, y.value) << s;
            EXPECT_EQ(x.upper, y.upper) << s;
            for (uint8_t t = 1; t <= with.max_mismatches; ++t) {
                EXPECT_EQ(Relation::UNKNOWN, b.contexts->by_mismatches.at(t).relation);
            }
            for (const auto &[p, count] : a.contexts->by_offset) {
                EXPECT_EQ(count.value, b.contexts->by_offset.at(p).value) << s << " " << p;
            }

            // the PARTIAL release after the stop: the same contexts, all of level 0
            Request partial_a = without;
            Request partial_b = with;
            partial_a.mode = partial_b.mode = Mode::PARTIAL;
            Result ra;
            Result rb;
            const std::vector<MCtx> la = released_by(*graph, c.text, partial_a, &ra, s);
            const std::vector<MCtx> lb = released_by(*graph, c.text, partial_b, &rb, s);
            EXPECT_EQ(la, lb) << s;
            EXPECT_EQ(ra.extraction->cut, rb.extraction->cut) << s;
        }
    }
    EXPECT_LT(1000u, stops);
}


// ---------------------------------------------------------------- stops

// the true counts behind a case, from O1
Truth truth_of(const Case &c, const DeBruijnGraph &masked) {
    return truth_of(hamming_walk(masked, c.text, c.request.scope, c.request.strands,
                                 c.request.max_mismatches),
                    c.request.max_mismatches, c.text.size(), c.k, c.request, c.text);
}

/**
 * The relations a stop leaves (SPEC §21.6): the levels below the stop's level as they
 * completed (|exact_below|: EXACT, on a masked graph; else a true relation), the stop's level
 * AT_LEAST (discovery), BOUNDS (its scans or check) or UNKNOWN (the clock reading that starts
 * its pass), the levels above UNKNOWN; every count a true relation; the steps per pass summing
 * to the steps. Returns the stop's level.
 */
uint8_t check_stopped(const Result &r, const Truth &truth, bool exact_below = true) {
    EXPECT_TRUE(r.stop);
    if (!r.stop)
        return 0;
    EXPECT_TRUE(r.stop->mismatch_level);
    const uint8_t t = r.stop->mismatch_level.value_or(0);
    const auto &levels = r.contexts->by_mismatches;
    EXPECT_EQ(truth.by_level.size(), levels.size());
    for (uint8_t u = 0; u < levels.size(); ++u) {
        EXPECT_TRUE(true_relation(levels[u], truth.by_level[u])) << "level " << int(u);
        if (u < t) {
            if (exact_below) {
                EXPECT_EQ(Relation::EXACT, levels[u].relation) << "level " << int(u);
            }
        } else if (u == t) {
            switch (levels[u].relation) {
                case Relation::AT_LEAST:
                    EXPECT_EQ(StopPhase::DISCOVERY, r.stop->phase);
                    break;
                case Relation::BOUNDS:
                    EXPECT_EQ(StopPhase::MASK_SCAN, r.stop->phase);
                    break;
                case Relation::UNKNOWN:
                    // the clock at the pass's start (or a stop before the pattern)
                    EXPECT_EQ(StopPhase::DISCOVERY, r.stop->phase);
                    EXPECT_EQ(0u, r.work.steps_by_mismatches.empty()
                                    ? 0 : r.work.steps_by_mismatches.at(u));
                    break;
                case Relation::EXACT:
                    ADD_FAILURE() << "the interrupted level " << int(u) << " is exact";
                    break;
            }
        } else {
            EXPECT_EQ(Relation::UNKNOWN, levels[u].relation) << "level " << int(u);
        }
    }
    EXPECT_NE(Relation::EXACT, r.contexts->total.relation);
    EXPECT_TRUE(true_relation(r.contexts->total, truth.total));
    for (const auto &[p, n] : truth.by_offset) {
        EXPECT_TRUE(true_relation(r.contexts->by_offset.at(p), n)) << "offset " << p;
    }
    for (const auto &[o, n] : truth.by_orientation) {
        EXPECT_TRUE(true_relation(r.contexts->by_orientation.at(o), n));
    }
    EXPECT_EQ(r.work.steps, std::accumulate(r.work.steps_by_mismatches.begin(),
                                            r.work.steps_by_mismatches.end(), uint64_t(0)));
    return t;
}

// the first pass whose cumulative steps exceed |s|: where a run with max_steps s stops
uint8_t pass_of_step(const std::vector<uint64_t> &steps_by_level, uint64_t s) {
    uint64_t cumulative = 0;
    for (uint8_t t = 0; t < steps_by_level.size(); ++t) {
        cumulative += steps_by_level[t];
        if (cumulative > s)
            return t;
    }
    return static_cast<uint8_t>(steps_by_level.size());
}

/**
 * A PARTIAL release after a stop in the pass of level |t|: in answer order, every context of
 * the levels below t, then some of level t (each one O1's), none above.
 */
void check_partial_after_stop(const std::vector<MCtx> &released,
                              const std::vector<MCtx> &expected, uint8_t t) {
    EXPECT_TRUE(std::is_sorted(released.begin(), released.end()));
    EXPECT_EQ(released.end(), std::adjacent_find(released.begin(), released.end()));
    const std::set<std::tuple<uint8_t, node_index, uint32_t, Orientation>> all = [&]() {
        std::set<std::tuple<uint8_t, node_index, uint32_t, Orientation>> s;
        for (const MCtx &c : expected) {
            s.emplace(c.level, c.node, c.offset, c.orientation);
        }
        return s;
    }();
    size_t below = 0;
    for (const MCtx &c : expected) {
        below += c.level < t;
    }
    ASSERT_LE(below, released.size());
    EXPECT_TRUE(std::equal(expected.begin(), expected.begin() + below, released.begin()));
    for (const MCtx &c : released) {
        EXPECT_LE(c.level, t);
        EXPECT_TRUE(all.count({ c.level, c.node, c.offset, c.orientation })) << c;
    }
}

TEST(PatternMismatch, StoppedAtEveryStep) {
    // every max_steps below a case's total: the stop in the pass the cumulative steps name,
    // the levels below it complete, its own AT_LEAST or BOUNDS, the ones above UNKNOWN;
    // PARTIAL a prefix of the closest levels and some of the stop's level; ALL_OR_COUNT
    // withheld; the same answer twice. Even-k PRIMARY cases put stops in the palindrome scans
    std::map<StopPhase, size_t> phases;
    std::set<uint8_t> levels;
    size_t cases = 0;
    for (uint32_t seed = 3000; cases < 24; ++seed) {
        Case c = random_case(seed);
        if (seed % 2) {
            c.k = 6;
            c.mode = DeBruijnGraph::PRIMARY;
            c.graph = build(c.k, c.records, c.mode);
            c.request.scope = Scope::ANY_OFFSET;
            if (c.text.size() > c.k)
                c.text.resize(c.k);
        }
        const Result full = count_of(*c.graph, c.text, c.request);
        if (full.work.steps < 50 || full.work.steps > 20'000
                || !full.contexts->by_mismatches.back().value)
            continue;
        ++cases;
        SCOPED_TRACE(describe(c));
        const Truth truth = truth_of(c, *c.graph);
        const std::vector<MCtx> expected = hamming_walk(*c.graph, c.text, c.request.scope,
                                                        c.request.strands,
                                                        c.request.max_mismatches);
        const uint64_t total = full.work.steps;
        const uint64_t stride = std::max<uint64_t>(1, total / 400);
        for (uint64_t s = 0; s < total; s += (s < 200 ? 1 : stride)) {
            const Result r = count_of(*c.graph, c.text, c.request, s);
            ASSERT_TRUE(r.stop) << s;
            EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);
            const uint8_t t = check_stopped(r, truth);
            EXPECT_EQ(pass_of_step(full.work.steps_by_mismatches, s), t) << s;
            ++phases[r.stop->phase];
            levels.insert(t);
            // deterministic
            const Result again = count_of(*c.graph, c.text, c.request, s);
            EXPECT_EQ(r.work.steps, again.work.steps);
            EXPECT_EQ(r.contexts->total.value, again.contexts->total.value);

            Request partial = c.request;
            partial.mode = Mode::PARTIAL;
            Result p;
            check_partial_after_stop(released_by(*c.graph, c.text, partial, &p, s), expected, t);
            EXPECT_EQ(r.work.steps, p.work.steps);
            EXPECT_EQ(StopReason::MAX_STEPS, p.extraction->cut);
            Result w;
            EXPECT_TRUE(released_by(*c.graph, c.text, c.request, &w, s).empty());
            EXPECT_EQ(Withheld::DISCOVERY_BUDGET, w.extraction->withheld);
        }
    }
    EXPECT_LT(100u, phases[StopPhase::DISCOVERY]);
    EXPECT_LT(5u, phases[StopPhase::MASK_SCAN]);
    EXPECT_EQ((std::set<uint8_t> { 0, 1, 2, 3 }), levels);
}

// a budget whose clock counts its readings in |*readings| and reads late from the
// |on_time| + 1-th on
Budget counted_budget(std::shared_ptr<uint64_t> readings, uint64_t on_time) {
    const Clock::time_point start = Clock::now();
    return Budget(kManySteps, Deadline(start, 1000, 250, [start, readings, on_time]() {
        return start + std::chrono::seconds(++*readings > on_time ? 100 : 0);
    }));
}

TEST(PatternMismatch, StoppedAtEveryClockReading) {
    // a virtual clock late from its n-th reading on, for every n of a case's readings: the
    // reading before the pattern, at the start of every pass after the first, and every 4,096
    // steps inside a pass (a graph large enough for several per pass). The stop is
    // {discovery or mask_scan, time}, time_limited, the levels as for a step stop; a stop at
    // a pass's own reading leaves that level UNKNOWN with no step charged in it
    std::mt19937 rng(41);
    const std::vector<std::string> records = random_records(rng, 6, 1500, 2500);
    size_t pass_start_stops = 0;
    size_t inside_stops = 0;
    for (const auto &[mode, k, text, m] : {
            std::make_tuple(DeBruijnGraph::BASIC, size_t(11), std::string("ACRTG"), uint8_t(2)),
            std::make_tuple(DeBruijnGraph::PRIMARY, size_t(12), std::string("GATCA"),
                            uint8_t(3)),
            std::make_tuple(DeBruijnGraph::CANONICAL, size_t(9), std::string("ANNGT"),
                            uint8_t(1)) }) {
        auto graph = build(k, records, mode);
        Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, m, Mode::COUNT);
        SCOPED_TRACE(text + " " + mode_name(mode));
        const PatternSearch engine(*graph);
        auto readings = std::make_shared<uint64_t>(0);
        Budget unlimited = counted_budget(readings, std::numeric_limits<uint64_t>::max());
        const Result full = engine.count(parse(text), request, unlimited);
        ASSERT_FALSE(full.stop);
        const uint64_t total_readings = *readings;
        ASSERT_LT(uint64_t(m) + 3, total_readings);
        const Truth truth = truth_of(hamming_walk(*graph, text, request.scope, request.strands,
                                                  m),
                                     m, text.size(), k, request, text);
        for (uint64_t n = 0; n < total_readings; ++n) {
            auto count = std::make_shared<uint64_t>(0);
            Budget budget = counted_budget(count, n);
            const Result r = engine.count(parse(text), request, budget);
            ASSERT_TRUE(r.stop) << n;
            EXPECT_EQ(StopReason::TIME, r.stop->reason);
            EXPECT_TRUE(r.time_limited);
            const uint8_t t = check_stopped(r, truth);
            if (r.contexts->by_mismatches.at(t).relation == Relation::UNKNOWN) {
                pass_start_stops += t > 0;
            } else {
                ++inside_stops;
            }
        }
    }
    EXPECT_LE(5u, pass_start_stops);
    EXPECT_LT(20u, inside_stops);
}

TEST(PatternMismatch, ThresholdStopsLeaveTheCloserLevelsComplete) {
    // stop_at_threshold compares the running total over the levels: a stop in the pass of
    // level t means the contexts of levels <= t exceed max_contexts, and the levels below t
    // are complete; ALL_OR_COUNT then withholds THRESHOLD_CROSSED
    size_t stopped = 0;
    std::set<uint8_t> levels;
    for (uint32_t seed = 4000; seed < 4040; ++seed) {
        Case c = random_case(seed);
        if (c.mode == DeBruijnGraph::PRIMARY)
            continue;
        SCOPED_TRACE(describe(c));
        const Truth truth = truth_of(c, *c.graph);
        for (uint64_t cap = 0; cap <= truth.total; ++cap) {
            Request request = c.request;
            request.stop_at_threshold = true;
            request.max_contexts = cap;
            const Result r = count_of(*c.graph, c.text, request);
            if (!r.stop) {
                EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
                EXPECT_EQ(truth.total, r.contexts->total.value);
                continue;
            }
            ++stopped;
            EXPECT_EQ(StopReason::MAX_CONTEXTS, r.stop->reason);
            const uint8_t t = check_stopped(r, truth);
            levels.insert(t);
            uint64_t up_to = 0;
            for (uint8_t u = 0; u <= t; ++u) {
                up_to += truth.by_level[u];
            }
            EXPECT_LT(cap, up_to);
            Result w;
            EXPECT_TRUE(released_by(*c.graph, c.text, request, &w).empty());
            EXPECT_EQ(Withheld::THRESHOLD_CROSSED, w.extraction->withheld);
        }
    }
    EXPECT_LT(100u, stopped);
    EXPECT_LE(3u, levels.size());
}


// ---------------------------------------------------------------- graphs without the mask

TEST(PatternMismatch, UnmaskedTwinsPerLevel) {
    // the same records with and without the mask: every EXACT level equals the masked one,
    // every BOUNDS level holds it; the check of few unchecked candidates runs after each pass
    // with a pattern-wide allowance X: level t is EXACT iff it has no unchecked candidate or
    // the levels up to t have at most X in all; ALL_OR_COUNT's release makes every level EXACT
    size_t bounded = 0;
    size_t checked = 0;
    for (uint32_t seed = 5000; seed < 5080; ++seed) {
        Case c = random_case(seed);
        if (c.mode == DeBruijnGraph::PRIMARY)
            continue;
        auto unmasked = unmasked_copy(c);
        SCOPED_TRACE(describe(c));
        const std::vector<MCtx> expected = hamming_walk(*c.graph, c.text, c.request.scope,
                                                        c.request.strands,
                                                        c.request.max_mismatches);
        const Truth truth = truth_of(expected, c.request.max_mismatches, c.text.size(), c.k,
                                     c.request, c.text);
        const Result plain = count_of(*unmasked, c.text, c.request);
        ASSERT_FALSE(plain.stop);
        std::vector<uint64_t> unchecked;
        for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
            const Count &level = plain.contexts->by_mismatches.at(t);
            EXPECT_TRUE(true_relation(level, truth.by_level[t])) << int(t);
            unchecked.push_back(level.relation == Relation::BOUNDS ? level.upper - level.lower
                                                                   : 0);
            bounded += level.relation == Relation::BOUNDS;
        }
        // the allowance at, below and above every level's cumulative unchecked candidates
        std::set<uint64_t> limits { 0 };
        uint64_t cumulative = 0;
        for (uint64_t u : unchecked) {
            cumulative += u;
            limits.insert({ cumulative, cumulative + 1 });
            if (cumulative)
                limits.insert(cumulative - 1);
        }
        for (uint64_t limit : limits) {
            Request request = c.request;
            request.max_checked_entries = limit;
            const Result r = count_of(*unmasked, c.text, request);
            ASSERT_FALSE(r.stop);
            uint64_t so_far = 0;
            for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
                so_far += unchecked[t];
                const Count &level = r.contexts->by_mismatches.at(t);
                const bool exact = !unchecked[t] || so_far <= limit;
                EXPECT_EQ(exact ? Relation::EXACT : Relation::BOUNDS, level.relation)
                    << "level " << int(t) << " limit " << limit;
                EXPECT_TRUE(true_relation(level, truth.by_level[t]));
                checked += exact && unchecked[t];
            }
        }
        // the release enumerates every candidate and drops the source dummies
        Result released;
        EXPECT_EQ(expected, released_by(*unmasked, c.text, c.request, &released));
        for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
            EXPECT_EQ(Relation::EXACT, released.contexts->by_mismatches.at(t).relation);
            EXPECT_EQ(truth.by_level[t], released.contexts->by_mismatches.at(t).value);
        }
    }
    EXPECT_LT(30u, bounded);
    EXPECT_LT(30u, checked);
}


TEST(PatternMismatch, UnmaskedReleasesPerLevel) {
    // without the mask, every mode (wrapped PRIMARY at odd and even k included): each level a
    // true relation of the masked count; a release that enumerated every candidate makes each
    // level EXACT (ALL_OR_COUNT, and PARTIAL with a cap above the total); a PARTIAL cut short
    // releases the masked list's prefix and raises each level's lower bound to what it
    // released there
    size_t cases = 0;
    for (uint32_t seed = 7000; seed < 7060; ++seed) {
        Case c = random_case(seed);
        auto unmasked = unmasked_copy(c);
        SCOPED_TRACE(describe(c));
        const std::vector<MCtx> expected = hamming_walk(*c.graph, c.text, c.request.scope,
                                                        c.request.strands,
                                                        c.request.max_mismatches);
        const Truth truth = truth_of(expected, c.request.max_mismatches, c.text.size(), c.k,
                                     c.request, c.text);
        const Result counted = count_of(*unmasked, c.text, c.request);
        ASSERT_FALSE(counted.stop);
        for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
            EXPECT_TRUE(true_relation(counted.contexts->by_mismatches.at(t), truth.by_level[t]));
        }
        Result all;
        EXPECT_EQ(expected, released_by(*unmasked, c.text, c.request, &all));
        Request partial = c.request;
        partial.mode = Mode::PARTIAL;
        for (uint64_t cap : { uint64_t(expected.size() / 2), uint64_t(expected.size() + 1) }) {
            partial.max_contexts = cap;
            Result r;
            const std::vector<MCtx> got = released_by(*unmasked, c.text, partial, &r);
            EXPECT_EQ(std::vector<MCtx>(expected.begin(),
                                        expected.begin() + std::min<size_t>(cap, expected.size())),
                      got);
            std::vector<uint64_t> released(c.request.max_mismatches + 1, 0);
            for (const MCtx &x : got) {
                ++released[x.level];
            }
            for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
                const Count &level = r.contexts->by_mismatches.at(t);
                EXPECT_TRUE(true_relation(level, truth.by_level[t])) << int(t) << " " << cap;
                if (cap > expected.size()) {
                    EXPECT_EQ(Relation::EXACT, level.relation) << int(t);
                } else {
                    EXPECT_LE(released[t], level.value) << int(t);
                }
            }
        }
        ++cases;
    }
    EXPECT_EQ(60u, cases);
}

// ---------------------------------------------------------------- order and retention

TEST(PatternMismatch, PartialKeepsTheClosestLevels) {
    // a short pattern with ~10^5 contexts, most of them at level 1: PARTIAL releases the first
    // max_contexts in (level, node, offset, orientation) order and keeps about that many
    // ranges (plus the straddling ones), not the discovery's
    std::mt19937 rng(7);
    const std::vector<std::string> records = random_records(rng, 10, 2000, 2000);
    for (const auto &[mode, k] : { std::make_pair(DeBruijnGraph::BASIC, size_t(13)),
                                   std::make_pair(DeBruijnGraph::PRIMARY, size_t(12)) }) {
        auto graph = build(k, records, mode);
        const std::string text = "ACG";
        const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 1);
        const std::vector<MCtx> expected = hamming_walk(*graph, text, request.scope,
                                                        request.strands, 1);
        ASSERT_GT(expected.size(), 50'000u);
        Result all;
        const std::vector<MCtx> full = released_by(*graph, text, request, &all);
        ASSERT_EQ(expected, full) << mode_name(mode);
        ASSERT_GT(all.work.spans_retained_peak, 20'000u);
        Request partial = request;
        partial.mode = Mode::PARTIAL;
        for (uint64_t cap : { 0, 1, 7, 100, 4097, 20'000 }) {
            partial.max_contexts = cap;
            Result r;
            const std::vector<MCtx> got = released_by(*graph, text, partial, &r);
            ASSERT_EQ(std::vector<MCtx>(full.begin(), full.begin() + cap), got)
                << mode_name(mode) << " cap " << cap;
            EXPECT_EQ(StopReason::MAX_CONTEXTS, r.extraction->cut);
            EXPECT_EQ(all.work.steps, r.work.steps);
            EXPECT_LE(r.work.spans_retained_peak, std::max<uint64_t>(4096, 2 * (cap + 500)) + 1)
                << mode_name(mode) << " cap " << cap;
            if (!cap) {
                EXPECT_EQ(0u, r.work.spans_retained_peak);
            }
        }
    }
}

TEST(PatternMismatch, CheaperOrientationFirstInEveryPass) {
    // an internal N run early in P is a wide fan for P and a narrow one for rc(P) (whose run
    // comes late): in pass 1 rc(P) runs first, so that a step stop there leaves the reverse
    // orientation complete before the forward one
    std::mt19937 rng(11);
    auto graph = build(13, random_records(rng, 6, 800, 800), DeBruijnGraph::BASIC);
    const std::string text = "ANNNNACGTACG";
    const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 1, Mode::COUNT);
    const Result full = count_of(*graph, text, request);
    ASSERT_FALSE(full.stop);
    const uint64_t pass0 = full.work.steps_by_mismatches.at(0);
    size_t reverse_only = 0;
    for (uint64_t s = pass0; s < full.work.steps; ++s) {
        const Result r = count_of(*graph, text, request, s);
        ASSERT_TRUE(r.stop);
        ASSERT_EQ(1, *r.stop->mismatch_level);
        const bool forward = r.contexts->by_orientation.at(Orientation::FORWARD).relation
                                == Relation::EXACT;
        const bool reverse = r.contexts->by_orientation.at(Orientation::REVERSE).relation
                                == Relation::EXACT;
        EXPECT_FALSE(forward && !reverse) << s;
        reverse_only += reverse && !forward;
    }
    EXPECT_LT(0u, reverse_only);
}


TEST(PatternMismatch, MatchingBaseFirst) {
    // within a pass the DFS follows the exact prefix and places the substitution from the
    // deepest position back: a pass 1 stopped step by step finds the graph's single mutants of
    // P in the order of their substituted position from the last to the first (7, 4, 1, 0),
    // the first position's last of all (the matching base at depth 0 is tried first too)
    std::mt19937 rng(29);
    std::string p;
    for (int i = 0; i < 8; ++i) {
        p.push_back("ACGT"[rng() % 4]);
    }
    std::vector<std::string> records;
    for (size_t position : { 0, 1, 4, 7 }) {
        std::string variant = p;
        // another base
        variant[position] = "ACGT"[(std::string("ACGT").find(p[position]) + 1 + rng() % 3) % 4];
        records.push_back(random_records(rng, 1, 6, 6)[0] + variant
                          + random_records(rng, 1, 6, 6)[0]);
    }
    auto graph = build(12, records, DeBruijnGraph::BASIC);
    Request request = make_request(Scope::SUFFIX, Strands::FORWARD, 1);
    check_against_oracles(*graph, p, request, &records, DeBruijnGraph::BASIC);
    const Result full = count_of(*graph, p, request);
    ASSERT_EQ(4u, full.contexts->by_mismatches.at(1).value);
    request.mode = Mode::PARTIAL;
    std::vector<uint32_t> order;
    for (uint64_t s = full.work.steps_by_mismatches.at(0); s <= full.work.steps; ++s) {
        Result r;
        for (const MCtx &c : released_by(*graph, p, request, &r, s)) {
            if (c.level != 1)
                continue;
            ASSERT_EQ(1u, c.positions.size());
            if (std::find(order.begin(), order.end(), c.positions[0]) == order.end())
                order.push_back(c.positions[0]);
        }
    }
    EXPECT_EQ((std::vector<uint32_t> { 7, 4, 1, 0 }), order);
}

// ---------------------------------------------------------------- named cases

std::vector<MCtx> released_all(const DeBruijnGraph &graph, const std::string &text,
                               uint8_t m, Strands strands = Strands::BOTH,
                               Scope scope = Scope::ANY_OFFSET) {
    Result r;
    return released_by(graph, text, make_request(scope, strands, m), &r);
}

TEST(PatternMismatch, DoubleMutantWithoutItsIntermediates) {
    // the NOTE's case: AACCA is two substitutions from AAAAA, and neither single mutant
    // (AACAA, AAACA) is in the graph; the shared DFS still reaches it at level 2
    const std::vector<std::string> records { "AAAAA", "AACCA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    const std::vector<MCtx> contexts = released_all(*graph, "AAAAA", 2, Strands::FORWARD);
    ASSERT_EQ(2u, contexts.size());
    EXPECT_EQ(0, contexts[0].level);
    EXPECT_EQ(2, contexts[1].level);
    EXPECT_EQ((std::vector<uint32_t> { 2, 3 }), contexts[1].positions);
    EXPECT_EQ("AACCA", graph->get_node_sequence(contexts[1].node));
    EXPECT_EQ(1u, released_all(*graph, "AAAAA", 1, Strands::FORWARD).size());
    check_against_oracles(*graph, "AAAAA", make_request(Scope::ANY_OFFSET, Strands::BOTH, 2),
                          &records, DeBruijnGraph::BASIC);
}

TEST(PatternMismatch, IupacSetsAndSubstitutions) {
    // at R (A, G): A and G keep the level, C and T substitute; at B (C, G, T) only A does; at
    // N nothing does
    const std::vector<std::string> records { "CAAAC", "CGAAC", "CCAAC", "CTAAC" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    auto level_of = [&](const std::string &text, const std::string &kmer) -> int {
        for (const MCtx &c : released_all(*graph, text, 1, Strands::FORWARD)) {
            if (graph->get_node_sequence(c.node) == kmer)
                return c.level;
        }
        return -1;
    };
    EXPECT_EQ(0, level_of("CRAAC", "CAAAC"));
    EXPECT_EQ(0, level_of("CRAAC", "CGAAC"));
    EXPECT_EQ(1, level_of("CRAAC", "CCAAC"));
    EXPECT_EQ(1, level_of("CRAAC", "CTAAC"));
    EXPECT_EQ(1, level_of("CBAAC", "CAAAC"));
    EXPECT_EQ(0, level_of("CBAAC", "CCAAC"));
    EXPECT_EQ(0, level_of("CBAAC", "CGAAC"));
    EXPECT_EQ(0, level_of("CNAAC", "CTAAC"));
    for (const std::string text : { "CRAAC", "CBAAC", "CNAAC", "NRNBN" }) {
        for (uint8_t m = 1; m <= 3; ++m) {
            check_against_oracles(*graph, text, make_request(Scope::ANY_OFFSET, Strands::BOTH, m),
                                  &records, DeBruijnGraph::BASIC);
        }
    }
    // a pattern of N only: no position can hold a substitution, every level above 0 EXACT 0
    // and the answer that of max_mismatches 0
    std::mt19937 rng(3);
    auto big = build(7, random_records(rng, 4, 50, 80), DeBruijnGraph::CANONICAL);
    for (uint8_t m = 1; m <= 3; ++m) {
        const Result r = count_of(*big, "NNNN", make_request(Scope::ANY_OFFSET, Strands::BOTH, m));
        const Result v1 = count_of(*big, "NNNN", make_request(Scope::ANY_OFFSET, Strands::BOTH, 0));
        EXPECT_EQ(v1.contexts->total.value, r.contexts->by_mismatches.at(0).value);
        EXPECT_EQ(v1.contexts->total.value, r.contexts->total.value);
        EXPECT_EQ(v1.work.steps, r.work.steps);
        for (uint8_t t = 1; t <= m; ++t) {
            EXPECT_EQ(Relation::EXACT, r.contexts->by_mismatches.at(t).relation);
            EXPECT_EQ(0u, r.contexts->by_mismatches.at(t).value);
            EXPECT_EQ(0u, r.work.steps_by_mismatches.at(t));
        }
    }
}

TEST(PatternMismatch, BothOrientationsAtOnePlace) {
    // AACGTA and its reverse complement TACGTT are each one substitution from AACGTT: two
    // contexts at one (k-mer, offset), + and -, both of level 1; none without mismatches
    const std::vector<std::string> records { "AACGTT" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    const std::vector<MCtx> contexts = released_all(*graph, "AACGTA", 1);
    ASSERT_EQ(2u, contexts.size());
    EXPECT_EQ(contexts[0].node, contexts[1].node);
    EXPECT_EQ(contexts[0].offset, contexts[1].offset);
    EXPECT_EQ(Orientation::FORWARD, contexts[0].orientation);
    EXPECT_EQ(Orientation::REVERSE, contexts[1].orientation);
    EXPECT_EQ(1, contexts[0].level);
    EXPECT_EQ(1, contexts[1].level);
    // P's position 5 in both: found T (forward), and A complemented to T (reverse)
    EXPECT_EQ(std::vector<uint32_t> { 5 }, contexts[0].positions);
    EXPECT_EQ(std::vector<uint32_t> { 5 }, contexts[1].positions);
    EXPECT_TRUE(released_all(*graph, "AACGTA", 0).empty());
}

TEST(PatternMismatch, SeveralOffsetsInOneKmer) {
    // AAAA in a k-mer of A's with one C: a context at every offset, at level 0 where the C
    // lies outside the window and level 1 where it lies inside
    const std::vector<std::string> records { "AAAACAAAA" };
    auto graph = build(9, records, DeBruijnGraph::BASIC);
    const std::vector<MCtx> contexts = released_all(*graph, "AAAA", 1, Strands::FORWARD);
    ASSERT_EQ(6u, contexts.size());
    std::vector<uint32_t> level0;
    std::vector<uint32_t> level1;
    for (const MCtx &c : contexts) {
        (c.level ? level1 : level0).push_back(c.offset);
    }
    EXPECT_EQ((std::vector<uint32_t> { 0, 5 }), level0);
    EXPECT_EQ((std::vector<uint32_t> { 1, 2, 3, 4 }), level1);
    check_against_oracles(*graph, "AAAA", make_request(Scope::ANY_OFFSET, Strands::BOTH, 1),
                          &records, DeBruijnGraph::BASIC);
}

TEST(PatternMismatch, PalindromesPerLevel) {
    // a palindromic P is searched once at every level; on an even-k wrapped PRIMARY graph a
    // palindromic stored k-mer is one context however both base searches meet it, also where
    // it is possible only from level 2 on (AAAAAT against the palindrome AAATTT: two of the
    // three mirrored pairs admit no complementary bases)
    std::mt19937 rng(5);
    std::vector<std::string> records = random_records(rng, 3, 40, 60);
    records.push_back("GGAAATTTCC");
    records.push_back("ACGCGT");
    for (size_t k : { 6, 8 }) {
        auto graph = build(k, records, DeBruijnGraph::PRIMARY);
        for (const std::string text : { "AAAAAT", "ACGCGT", "AATT", "ACGT", "AAAAAA" }) {
            for (uint8_t m = 1; m <= 3; ++m) {
                SCOPED_TRACE(text + " k " + std::to_string(k) + " m " + std::to_string(m));
                check_against_oracles(*graph, text,
                                      make_request(Scope::ANY_OFFSET, Strands::BOTH, m), &records,
                                      DeBruijnGraph::PRIMARY);
                const Result r = count_of(*graph, text,
                                          make_request(Scope::ANY_OFFSET, Strands::BOTH, m));
                const bool palindromic = text == reverse_complement_string(text);
                EXPECT_EQ(palindromic ? 1u : 2u, r.contexts->by_orientation.size());
            }
        }
        // AAATTT holds AAAAAT at level 2 and is its own reverse complement: one context
        bool found = false;
        for (const MCtx &c : released_all(*graph, "AAAAAT", 2, Strands::FORWARD)) {
            const std::string kmer = graph->get_node_sequence(c.node);
            if (k == 6 && kmer == "AAATTT") {
                EXPECT_FALSE(found);
                found = true;
                EXPECT_EQ(2, c.level);
            }
        }
        EXPECT_EQ(k == 6, found);
    }
}

TEST(PatternMismatch, GraphBoundariesAndLeadingN) {
    // instances at the records' first and last k-mers, records shorter than k, the leading N
    // run skipped on $ACGT (its offsets shifted back), an N run inside and at the end
    const std::vector<std::string> records { "TTGCAACGTAGG", "ACG", "GGCATT", "CCTANNAGGTCA" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        auto graph = build(6, records, mode);
        for (const std::string text : { "TTGCAA", "GTAGG", "NNNGCA", "NNTAGG", "GNNNAG",
                                        "AGGNN", "NNNNNC", "CATT" }) {
            for (uint8_t m = 1; m <= 3; ++m) {
                for (Scope scope : { Scope::ANY_OFFSET, Scope::SUFFIX }) {
                    SCOPED_TRACE(text + " m " + std::to_string(m) + " " + mode_name(mode));
                    check_against_oracles(*graph, text, make_request(scope, Strands::BOTH, m),
                                          &records, mode);
                }
            }
        }
    }
}

TEST(PatternMismatch, OnlyTheLastPositionsSubstitutable) {
    // reachability: a pass t places its t substitutions where positions are left that can hold
    // them; with N everywhere but the last t positions, every branch must survive to them, and
    // a pass above the substitutable positions is empty without a step
    std::mt19937 rng(9);
    std::vector<std::string> records = random_records(rng, 4, 30, 50);
    auto graph = build(7, records, DeBruijnGraph::BASIC);
    for (const std::string text : { "NNNNNAC", "NNNNNNA", "NNNNNGC", "ANNNNNC", "NNTNN" }) {
        for (uint8_t m = 1; m <= 3; ++m) {
            SCOPED_TRACE(text + " m " + std::to_string(m));
            check_against_oracles(*graph, text,
                                  make_request(Scope::ANY_OFFSET, Strands::BOTH, m), &records,
                                  DeBruijnGraph::BASIC);
            const Result r = count_of(*graph, text,
                                      make_request(Scope::ANY_OFFSET, Strands::BOTH, m));
            const size_t substitutable = std::count_if(text.begin(), text.end(),
                                                       [](char c) { return c != 'N'; });
            for (uint8_t t = substitutable + 1; t <= m; ++t) {
                EXPECT_EQ(0u, r.contexts->by_mismatches.at(t).value);
                EXPECT_EQ(Relation::EXACT, r.contexts->by_mismatches.at(t).relation);
                EXPECT_EQ(0u, r.work.steps_by_mismatches.at(t));
            }
        }
    }
}


// ---------------------------------------------------------------- the information floor

TEST(PatternMismatch, BallBitsAgainstTheBruteForceBall) {
    // mismatch_information_bits against every string enumerated (O5), for DNA and IUPAC
    // patterns of up to 10 positions and every m; 0 is information_bits itself
    std::mt19937 rng(17);
    for (int i = 0; i < 60; ++i) {
        const std::string text = random_pattern(rng, 1 + rng() % 10, i % 2);
        const Pattern pattern = parse(text);
        const std::vector<uint64_t> histogram = ball_histogram(text);
        EXPECT_EQ(pattern.information_bits(), pattern.mismatch_information_bits(0));
        for (uint8_t m = 0; m <= kMaxMismatches; ++m) {
            const double expected = ball_bits(histogram, m);
            const double bits = pattern.mismatch_information_bits(m);
            EXPECT_NEAR(expected, bits, 1e-9 * std::max(1.0, std::abs(expected)))
                << text << " m " << int(m);
        }
    }
    // an exact pattern: |ball| = sum of C(L, u) 3^u; a 31-mer loses 6.55, 12.06 and 16.94 bits
    // (94, 4,279 and 125,644 strings)
    const Pattern p31 = parse(std::string(31, 'A'));
    EXPECT_DOUBLE_EQ(94, exact_ball(31, 1));
    EXPECT_DOUBLE_EQ(4279, exact_ball(31, 2));
    EXPECT_DOUBLE_EQ(125644, exact_ball(31, 3));
    for (uint8_t m = 1; m <= 3; ++m) {
        EXPECT_NEAR(62 - std::log2(exact_ball(31, m)), p31.mismatch_information_bits(m), 1e-9);
    }
    EXPECT_NEAR(6.55, 62 - p31.mismatch_information_bits(1), 0.005);
    EXPECT_NEAR(12.06, 62 - p31.mismatch_information_bits(2), 0.005);
    EXPECT_NEAR(16.94, 62 - p31.mismatch_information_bits(3), 0.005);
    // a long pattern's ball overflows every integer type; the log-space DP stays finite
    const Pattern long_one = parse(std::string(5000, 'C'));
    EXPECT_NEAR(10000 - std::log2(exact_ball(5000, 3)), long_one.mismatch_information_bits(3),
                1e-6);
    // a window of positions, as for an anchor window
    EXPECT_NEAR(ball_bits(ball_histogram("ACGTR"), 2),
                parse("TTACGTRAA").mismatch_information_bits(2, 7, 2), 1e-9);
}

TEST(PatternMismatch, FloorOnTheBitsLeft) {
    // the floor (24 bits) gates the bits left after m substitutions, with no exemption:
    // exact DNA of 14/15 bases at m = 1, 17/18 at m = 2, 19/20 at m = 3 refused / served; a
    // refused pattern states the largest m' < m the floor admits
    std::mt19937 rng(23);
    auto graph = build(31, random_records(rng, 3, 200, 200), DeBruijnGraph::BASIC);
    auto exact = [&](size_t L) {
        std::string text;
        for (size_t i = 0; i < L; ++i) {
            text.push_back("ACGT"[rng() % 4]);
        }
        return text;
    };
    for (const auto &[L, m, served, within] : {
            std::make_tuple(14, 1, false, 0), std::make_tuple(15, 1, true, -1),
            std::make_tuple(17, 2, false, 1), std::make_tuple(18, 2, true, -1),
            std::make_tuple(19, 3, false, 2), std::make_tuple(20, 3, true, -1),
            std::make_tuple(11, 1, false, -1) }) {
        for (Scope scope : { Scope::ANY_OFFSET, Scope::SUFFIX }) {
            Request request = make_request(scope, Strands::BOTH, m);
            request.min_information_bits = 24;
            const std::string text = exact(L);
            const Result r = count_of(*graph, text, request);
            const double bits = 2.0 * L - std::log2(exact_ball(L, m));
            ASSERT_TRUE(r.mismatch_information_bits);
            EXPECT_NEAR(bits, *r.mismatch_information_bits, 1e-9);
            EXPECT_EQ(bits >= 24, served);
            if (served) {
                EXPECT_FALSE(r.refusal) << L << " " << m;
                continue;
            }
            ASSERT_TRUE(r.refusal) << L << " " << m;
            EXPECT_EQ("information_below_floor", r.refusal->code);
            EXPECT_FALSE(r.contexts);
            EXPECT_EQ(0u, r.work.steps);
            // without mismatches an exact pattern in suffix scope is exempt: m' = 0 admitted
            const int expected_within = within < 0 && scope == Scope::SUFFIX ? 0 : within;
            if (expected_within < 0) {
                EXPECT_FALSE(r.max_mismatches_within_floor) << L;
            } else {
                ASSERT_TRUE(r.max_mismatches_within_floor) << L;
                EXPECT_EQ(expected_within, *r.max_mismatches_within_floor) << L;
            }
        }
    }
    // no exemption: an exact 12-mer in suffix scope is served without mismatches, refused
    // with one (18.8 bits)
    Request suffix = make_request(Scope::SUFFIX, Strands::BOTH, 0);
    suffix.min_information_bits = 24;
    const std::string twelve = exact(12);
    EXPECT_FALSE(count_of(*graph, twelve, suffix).refusal);
    suffix.max_mismatches = 1;
    const Result refused = count_of(*graph, twelve, suffix);
    ASSERT_TRUE(refused.refusal);
    EXPECT_EQ(0, refused.max_mismatches_within_floor.value_or(99));
    // an exact 31-mer in suffix scope at m = 1 is gated and passes (55.4 bits)
    suffix.max_mismatches = 1;
    const Result passed = count_of(*graph, exact(31), suffix);
    EXPECT_FALSE(passed.refusal);
    EXPECT_NEAR(55.4, *passed.mismatch_information_bits, 0.05);
}


// ---------------------------------------------------------------- the engine API

TEST(PatternMismatch, RequestsNotServedYet) {
    // above the contract's 3, a peptide, a pattern longer than k: refused as arguments (the
    // route checks them before; 5m-2 and 5m-4 serve the latter two)
    auto graph = build(5, { "ACGTACGTAC" }, DeBruijnGraph::BASIC);
    const PatternSearch engine(*graph);
    Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 4);
    Budget budget = budget_of();
    EXPECT_THROW(engine.count(parse("ACG"), request, budget), std::invalid_argument);
    request.max_mismatches = 1;
    EXPECT_THROW(engine.count(Pattern::parse(PatternKind::PROTEIN, "MK"), request, budget),
                 std::invalid_argument);
    EXPECT_THROW(engine.count(parse("ACGTAC"), request, budget), std::invalid_argument);
    EXPECT_THROW(engine.enumerate(parse("ACGTAC"), request, budget, [](const Context&) {}),
                 std::invalid_argument);
    EXPECT_EQ(0u, budget.steps_used());
    // without mismatches nothing of the new answer is set
    request.max_mismatches = 0;
    const Result r = engine.count(parse("ACG"), request, budget);
    EXPECT_TRUE(r.contexts->by_mismatches.empty());
    EXPECT_TRUE(r.work.steps_by_mismatches.empty());
    EXPECT_FALSE(r.mismatch_information_bits);
    EXPECT_FALSE(r.max_mismatches_within_floor);
}

TEST(PatternMismatch, SubstitutionsInTheSubmittedPatternsCoordinates) {
    const Pattern p = parse("ACGRN");
    // forward: positions where the instance leaves the set; N never
    auto subs = substitutions(p, Orientation::FORWARD, "TCGCA");
    ASSERT_EQ(2u, subs.size());
    EXPECT_EQ(0u, subs[0].position);
    EXPECT_EQ('A', subs[0].expected);
    EXPECT_EQ('T', subs[0].found);
    EXPECT_EQ(3u, subs[1].position);
    EXPECT_EQ('R', subs[1].expected);
    EXPECT_EQ('C', subs[1].found);
    EXPECT_TRUE(substitutions(p, Orientation::FORWARD, "ACGGT").empty());
    // reverse: the instance instantiates rc(P) = NYCGT; P's position i is read at L - 1 - i,
    // complemented: the instance GGCGT holds G at rc's position 1 (Y = C, T): P's position 3
    // found complement(G) = C; GCCGT holds C there, inside Y
    EXPECT_TRUE(substitutions(p, Orientation::REVERSE, "GCCGT").empty());
    subs = substitutions(p, Orientation::REVERSE, "GGCGT");
    ASSERT_EQ(1u, subs.size());
    EXPECT_EQ(3u, subs[0].position);
    EXPECT_EQ('R', subs[0].expected);
    EXPECT_EQ('C', subs[0].found);
    // every substitution of a palindromic P reads forward
    EXPECT_EQ(1u, substitutions(parse("ACGT"), Orientation::PALINDROMIC, "ACGA").size());
    EXPECT_THROW(substitutions(p, Orientation::FORWARD, "ACG"), std::invalid_argument);
    EXPECT_THROW(substitutions(p, Orientation::FORWARD, "ACG$A"), std::invalid_argument);
    EXPECT_THROW(substitutions(p, Orientation::FORWARD, "ACGNA"), std::invalid_argument);
    EXPECT_THROW(substitutions(Pattern::parse(PatternKind::PROTEIN, "M"), Orientation::FORWARD,
                               "ATG"),
                 std::invalid_argument);
}

#endif // _DNA_GRAPH

} // namespace
