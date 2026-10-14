/**
 * Bounded mismatches in the pattern search (SPEC §21; the mismatch plan's increments 5m-0 to
 * 5m-2): DNA, IUPAC and peptide patterns searched with Request::max_mismatches 1 to 3, checked
 * against the oracles of pattern_mismatch_oracles.hpp, which read the pattern through their own
 * IUPAC table and genetic codes and share only the graph with the engine:
 *  - the Hamming graph-walk (O1) gives every context within m with its level and substituted
 *    positions, in answer order: every count per level, offset and orientation, every release
 *    and its order are checked against it;
 *  - the Hamming record-scan (O2) checks O1 against the records, without the graph's ids;
 *  - the candidate enumeration (O3) checks O1 by a different algorithm: every string of the
 *    ball searched exactly by the engine without mismatches;
 *  - the peptide residue oracle (O4) does the same for peptides: residue substitutions read
 *    codon by codon on the graph's k-mers (a long peptide's anchor windows with their cut
 *    codons included), and checked against the records' six-frame translations;
 *  - the brute-force ball (O5) checks the information bits left after m substitutions.
 * Requests without mismatches are the other suites' (their answers must not change); here the
 * pass of level 0 is checked to be that search, step for step. Patterns longer than k are
 * checked through their anchors (the extension with mismatches is not served yet).
 *
 * Alphabet: DNA (DNA4) builds only; the engine's DNA5 path is not served (alphabet_untested).
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../all/test_dbg_helpers.hpp"
#include "pattern_mismatch_oracles.hpp"
#include "pattern_test_support.hpp"

#include "graph/alignment/genetic_code.hpp"
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

Result count_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget = budget_of(max_steps);
    return PatternSearch(graph).count(pattern, request, budget);
}

Result count_of(const DeBruijnGraph &graph, const std::string &text, const Request &request,
                uint64_t max_steps = kManySteps) {
    return count_of(graph, parse(text), request, max_steps);
}

/**
 * The substituted positions of a released context, read off its instance by substitutions()
 * and checked against the oracles' table (expected is P's code, found the instance's base in
 * P's orientation and outside expected's bases); for a peptide its substituted residues by
 * residue_substitutions(), checked against the oracles' genetic code (expected is P's letter,
 * the codon the instance's in P's reading direction, found its residue there, one
 * substitution by O4's rule). |instance| is the context's L bases (a long pattern's anchor:
 * its k-mer, not read here).
 */
std::vector<uint32_t> substituted_in(const Pattern &pattern, Orientation orientation,
                                     const std::string &instance) {
    const std::string &text = pattern.text();
    const size_t L = pattern.length();
    std::vector<uint32_t> positions;
    if (pattern.kind() != PatternKind::PROTEIN) {
        for (const Substitution &s : substitutions(pattern, orientation, instance)) {
            positions.push_back(s.position);
            EXPECT_EQ(text[s.position], s.expected);
            const char base = orientation == Orientation::REVERSE
                ? complement_base(instance[L - 1 - s.position])
                : instance[s.position];
            EXPECT_EQ(base, s.found);
            EXPECT_EQ(std::string::npos, iupac_bases(s.expected).find(s.found));
        }
        return positions;
    }
    const OracleCode &code = oracle_code(pattern.genetic_code());
    for (const ResidueSubstitution &s : residue_substitutions(pattern, orientation, instance)) {
        positions.push_back(s.residue);
        EXPECT_EQ(text[s.residue], s.expected);
        const std::string codon = orientation == Orientation::REVERSE
            ? reverse_complement_string(instance.substr(L - 3 * s.residue - 3, 3))
            : instance.substr(3 * s.residue, 3);
        EXPECT_EQ(codon, s.codon);
        EXPECT_EQ(code.translate(codon), s.found) << codon;
        EXPECT_EQ(1, codon_level(s.expected, s.found).value_or(9)) << s.expected << s.found;
    }
    return positions;
}

/**
 * enumerate(): every released context as an MCtx, its level the engine's mismatches and its
 * positions (a peptide's residues) those substituted_in() reads off its instance, as many as
 * its level. A long pattern's anchors (release_anchors) have no whole instance: their
 * positions are left empty.
 */
std::vector<MCtx> released_by(const DeBruijnGraph &graph, const Pattern &pattern,
                              const Request &request, Budget &budget, Result *result) {
    const size_t L = pattern.length();
    const bool anchors = L > graph.get_k();
    std::vector<MCtx> released;
    *result = PatternSearch(graph).enumerate(pattern, request, budget, [&](const Context &c) {
        const std::string kmer = graph.get_node_sequence(c.node);
        MCtx m { c.node, c.offset, c.orientation, c.mismatches, {} };
        if (!anchors) {
            m.positions = substituted_in(pattern, c.orientation, kmer.substr(c.offset, L));
            EXPECT_EQ(c.mismatches, m.positions.size()) << pattern.text() << " " << kmer;
        }
        released.push_back(m);
    });
    return released;
}

std::vector<MCtx> released_by(const DeBruijnGraph &graph, const Pattern &pattern,
                              const Request &request, Result *result,
                              uint64_t max_steps = kManySteps) {
    Budget budget = budget_of(max_steps);
    return released_by(graph, pattern, request, budget, result);
}

std::vector<MCtx> released_by(const DeBruijnGraph &graph, const std::string &text,
                              const Request &request, Result *result,
                              uint64_t max_steps = kManySteps) {
    return released_by(graph, parse(text), request, result, max_steps);
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

// the oracle's counts of |expected| over the offsets and orientations searched
Truth truth_of(const std::vector<MCtx> &expected, uint8_t m,
               const std::vector<uint32_t> &offsets, const std::vector<Orientation> &searched) {
    Truth truth;
    truth.by_level.assign(m + 1, 0);
    for (uint32_t p : offsets) {
        truth.by_offset[p] = 0;
    }
    for (Orientation o : searched) {
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

std::vector<Orientation> searched_of(const std::string &text, Strands strands) {
    std::vector<Orientation> searched;
    for (const auto &[o, q] : orientations(text, strands)) {
        searched.push_back(o);
    }
    return searched;
}

Truth truth_of(const std::vector<MCtx> &expected, uint8_t m, size_t L, size_t k,
               const Request &request, const std::string &text) {
    return truth_of(expected, m, scope_offsets(L, k, request.scope),
                    searched_of(text, request.strands));
}

/**
 * Every check of one completed (graph, pattern, request) of L <= k against the oracle's list
 * |expected| (its orientations |searched|): every count EXACT per level, offset and
 * orientation; the steps per pass summing to the steps; enumerate() releasing the list in
 * (level, node, offset, orientation) order with each context's level and positions;
 * ALL_OR_COUNT's threshold on the total; PARTIAL's prefixes. Returns the number of contexts.
 */
size_t check_against(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                     const std::vector<MCtx> &expected, const std::vector<Orientation> &searched) {
    const size_t k = graph.get_k();
    const size_t L = pattern.length();
    const uint8_t m = request.max_mismatches;
    const std::string &text = pattern.text();
    const Truth truth = truth_of(expected, m, scope_offsets(L, k, request.scope), searched);

    Result counted = count_of(graph, pattern, request);
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

    EXPECT_EQ(searched, counted.searched);

    // the release: the oracle's list, levels and positions included, in answer order
    Result all;
    const std::vector<MCtx> released = released_by(graph, pattern, request, &all);
    EXPECT_EQ(expected, released) << text << " m " << int(m);
    EXPECT_TRUE(all.extraction && all.extraction->complete);
    EXPECT_EQ(counted.work.steps, all.work.steps);
    EXPECT_EQ(counted.work.steps_by_mismatches, all.work.steps_by_mismatches);

    // ALL_OR_COUNT: all levels or none, on the total
    if (!expected.empty()) {
        Request tight = request;
        tight.max_contexts = expected.size() - 1;
        Result withheld;
        EXPECT_TRUE(released_by(graph, pattern, tight, &withheld).empty());
        EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
    }
    // PARTIAL: the first max_contexts of the closest levels
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    for (uint64_t cap : { uint64_t(0), uint64_t(1), uint64_t(expected.size() / 2),
                          uint64_t(expected.size()) }) {
        partial.max_contexts = cap;
        Result cut;
        const std::vector<MCtx> prefix = released_by(graph, pattern, partial, &cut);
        const size_t n = std::min<size_t>(cap, expected.size());
        EXPECT_EQ(std::vector<MCtx>(expected.begin(), expected.begin() + n), prefix)
            << text << " cap " << cap;
        EXPECT_EQ(cap >= expected.size(), cut.extraction->complete);
    }
    return expected.size();
}

/**
 * check_against() for a DNA or IUPAC pattern against O1, and O1 against O2 with |records|.
 */
size_t check_against_oracles(const DeBruijnGraph &graph, const std::string &text,
                             const Request &request,
                             const std::vector<std::string> *records = nullptr,
                             DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    const uint8_t m = request.max_mismatches;
    const std::vector<MCtx> expected = hamming_walk(graph, text, request.scope,
                                                    request.strands, m);
    if (records) {
        // the graph holds exactly the records' k-mers (no cleaning in these builds)
        EXPECT_EQ(hamming_records(*records, graph.get_k(), mode, text, request.scope,
                                  request.strands, m),
                  spelled(graph, expected)) << text;
    }
    return check_against(graph, parse(text), request, expected,
                         searched_of(text, request.strands));
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
 * AT_LEAST (discovery), BOUNDS (its scans or check), UNKNOWN (the clock reading that starts
 * its pass) or, for a threshold crossed by the levels' counts once the pass's scans and check
 * are done, as it completed (the stop before the next pass); the levels above UNKNOWN; every
 * count a true relation; one entry of steps per level, summing to the steps. Returns the
 * stop's level.
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
                    if (r.stop->phase == StopPhase::DISCOVERY) {
                        // without the mask, its check not run: the threshold compared at the
                        // end of the pass, as for EXACT below
                        EXPECT_TRUE(r.stop->reason == StopReason::MAX_CONTEXTS
                                        || r.stop->reason == StopReason::MAX_ANCHORS)
                            << "the interrupted level " << int(u) << " is bounds";
                        EXPECT_GT(levels.size(), u + 1u) << "a stop after the last pass";
                        if (u + 1u < r.work.steps_by_mismatches.size()) {
                            EXPECT_EQ(0u, r.work.steps_by_mismatches[u + 1]);
                        }
                    } else {
                        EXPECT_EQ(StopPhase::MASK_SCAN, r.stop->phase);
                    }
                    break;
                case Relation::UNKNOWN:
                    // the clock at the pass's start (or a stop before the pattern)
                    EXPECT_EQ(StopPhase::DISCOVERY, r.stop->phase);
                    EXPECT_EQ(0u, r.work.steps_by_mismatches.at(u));
                    break;
                case Relation::EXACT:
                    // the threshold compared at the end of the pass: nothing of the next one
                    EXPECT_TRUE(r.stop->reason == StopReason::MAX_CONTEXTS
                                    || r.stop->reason == StopReason::MAX_ANCHORS)
                        << "the interrupted level " << int(u) << " is exact";
                    EXPECT_EQ(StopPhase::DISCOVERY, r.stop->phase);
                    EXPECT_GT(levels.size(), u + 1u) << "a stop after the last pass";
                    if (u + 1u < r.work.steps_by_mismatches.size()) {
                        EXPECT_EQ(0u, r.work.steps_by_mismatches[u + 1]);
                    }
                    break;
            }
        } else {
            EXPECT_EQ(Relation::UNKNOWN, levels[u].relation) << "level " << int(u);
        }
    }
    EXPECT_NE(Relation::EXACT, r.contexts->total.relation);
    EXPECT_TRUE(true_relation(r.contexts->total, truth.total));
    EXPECT_EQ(levels.size(), r.work.steps_by_mismatches.size());
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

// a clock that reads |start| plus |*now| milliseconds; the test moves |*now|
struct VirtualClock {
    Clock::time_point start = Clock::now();
    std::shared_ptr<double> now = std::make_shared<double>(0);

    Deadline deadline(double budget_ms, double reserve_ms) const {
        auto s = start;
        auto n = now;
        return Deadline(start, budget_ms, reserve_ms, [s, n]() {
            return s + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double, std::milli>(*n));
        });
    }
};

TEST(PatternMismatch, ExtractionStopCarriesNoLevel) {
    // a stop of a later phase carries no level (SPEC §21.7): the deadline during the PARTIAL
    // release of a request with mismatches (the caller's work per context clocked, as
    // PatternSearchFixes.PartialReleaseReadsTheClockPerReleasedContexts without them) stops
    // {EXTRACTION, TIME} with mismatch_level unset, the first kReleaseClockStride contexts in
    // answer order released, cut TIME, the steps per level and every level's count as the
    // count alone states them; without the mask a cut release raises each level's lower bound
    // to what it released there (§21.6), nothing else. ALL_OR_COUNT's delivery under the
    // clock withholds DEADLINE, again with no level, every level EXACT (the buffered release
    // enumerated every candidate)
    std::mt19937 rng(43);
    const std::vector<std::string> records = random_records(rng, 6, 600, 600);
    const std::string text = "ACGT";
    const uint8_t m = 2;
    const size_t k = 9;
    const Pattern pattern = parse(text);
    auto masked = build(k, records, DeBruijnGraph::BASIC);
    const std::vector<MCtx> expected = hamming_walk(*masked, text, Scope::ANY_OFFSET,
                                                    Strands::BOTH, m);
    ASSERT_GT(expected.size(), 300u);
    const Truth truth = truth_of(expected, m, text.size(), k,
                                 make_request(Scope::ANY_OFFSET, Strands::BOTH, m), text);
    for (bool unmasked : { false, true }) {
        SCOPED_TRACE(unmasked ? "unmasked" : "masked");
        std::shared_ptr<DeBruijnGraph> graph = masked;
        if (unmasked) {
            graph = build(k, records, DeBruijnGraph::BASIC);
            const_cast<DBGSuccinct&>(stored_graph(*graph)).reset_mask();
        }
        auto as_mctx = [&](const Context &c) {
            const std::string kmer = graph->get_node_sequence(c.node);
            return MCtx { c.node, c.offset, c.orientation, c.mismatches,
                          substituted_in(pattern, c.orientation,
                                         kmer.substr(c.offset, text.size())) };
        };
        Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, m, Mode::PARTIAL);
        const Result counted = count_of(*graph, text, request);
        ASSERT_FALSE(counted.stop);

        // the work time (10 ms) passes at the 10th context; the next reading is at the 64th
        VirtualClock vc;
        Budget budget(kManySteps, vc.deadline(260, 250));
        std::vector<MCtx> received;
        const Result r = PatternSearch(*graph).enumerate(pattern, request, budget,
                                                         [&](const Context &c) {
            received.push_back(as_mctx(c));
            *vc.now += 1;
        });
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopPhase::EXTRACTION, r.stop->phase);
        EXPECT_EQ(StopReason::TIME, r.stop->reason);
        EXPECT_FALSE(r.stop->mismatch_level);
        EXPECT_TRUE(r.time_limited);
        ASSERT_TRUE(r.extraction);
        EXPECT_EQ(StopReason::TIME, r.extraction->cut);
        EXPECT_FALSE(r.extraction->complete);
        EXPECT_FALSE(r.extraction->withheld);
        EXPECT_EQ(Budget::kReleaseClockStride, received.size());
        EXPECT_EQ(received.size(), r.extraction->returned);
        EXPECT_EQ(std::vector<MCtx>(expected.begin(), expected.begin() + received.size()),
                  received);
        EXPECT_EQ(counted.work.steps_by_mismatches, r.work.steps_by_mismatches);
        std::vector<uint64_t> released(m + 1, 0);
        for (const MCtx &c : received) {
            ++released[c.level];
        }
        for (uint8_t t = 0; t <= m; ++t) {
            const Count &before = counted.contexts->by_mismatches.at(t);
            const Count &after = r.contexts->by_mismatches.at(t);
            EXPECT_EQ(before.relation, after.relation) << "level " << int(t);
            EXPECT_TRUE(true_relation(after, truth.by_level[t])) << "level " << int(t);
            if (before.relation == Relation::EXACT) {
                EXPECT_EQ(before.value, after.value) << "level " << int(t);
            } else {
                EXPECT_EQ(Relation::BOUNDS, before.relation) << "level " << int(t);
                EXPECT_EQ(before.upper, after.upper) << "level " << int(t);
                EXPECT_LE(std::max(before.lower, released[t]), after.lower)
                    << "level " << int(t);
            }
        }
        EXPECT_EQ(counted.contexts->total.relation, r.contexts->total.relation);

        // ALL_OR_COUNT: every context costs the caller 1 ms; the work time is 100 ms
        request.mode = Mode::ALL_OR_COUNT;
        VirtualClock slow;
        Budget delivery(kManySteps, slow.deadline(350, 250));
        received.clear();
        const Result d = PatternSearch(*graph).enumerate(pattern, request, delivery,
                                                         [&](const Context &c) {
            received.push_back(as_mctx(c));
            *slow.now += 1;
        });
        ASSERT_TRUE(d.stop);
        EXPECT_EQ(StopPhase::EXTRACTION, d.stop->phase);
        EXPECT_EQ(StopReason::TIME, d.stop->reason);
        EXPECT_FALSE(d.stop->mismatch_level);
        EXPECT_TRUE(d.time_limited);
        ASSERT_TRUE(d.extraction);
        ASSERT_TRUE(d.extraction->withheld);
        EXPECT_EQ(Withheld::DEADLINE, *d.extraction->withheld);
        EXPECT_FALSE(d.extraction->complete);
        EXPECT_EQ(0u, d.extraction->returned);
        EXPECT_GE(received.size(), 100u);
        EXPECT_LT(received.size(), 100u + Budget::kReleaseClockStride);
        EXPECT_EQ(std::vector<MCtx>(expected.begin(), expected.begin() + received.size()),
                  received);
        for (uint8_t t = 0; t <= m; ++t) {
            const Count &level = d.contexts->by_mismatches.at(t);
            EXPECT_EQ(Relation::EXACT, level.relation) << "level " << int(t);
            EXPECT_EQ(truth.by_level[t], level.value) << "level " << int(t);
        }
    }
}

TEST(PatternMismatch, ThresholdStopsLeaveTheCloserLevelsComplete) {
    // stop_at_threshold compares the running total over the levels: a stop in the pass of
    // level t means the contexts of levels <= t exceed max_contexts, and the levels below t
    // are complete; ALL_OR_COUNT then withholds THRESHOLD_CROSSED. The comparison is repeated
    // at the end of each pass but the last, on the counts its scans and check completed: once
    // the levels 0 .. t alone exceed max_contexts, no later pass starts
    size_t stopped = 0;
    size_t after_its_pass = 0;
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
            // the first level up to which the closer levels alone exceed the cap, before the
            // last: the stop comes in its pass at the latest
            std::optional<uint8_t> first_above;
            uint64_t closer = 0;
            for (uint8_t u = 0; u + 1u < truth.by_level.size(); ++u) {
                closer += truth.by_level[u];
                if (closer > cap) {
                    first_above = u;
                    break;
                }
            }
            if (first_above) {
                EXPECT_TRUE(r.stop) << "cap " << cap;
                if (r.stop) {
                    EXPECT_GE(*first_above, r.stop->mismatch_level.value_or(kMaxMismatches + 1))
                        << "cap " << cap;
                }
            }
            if (!r.stop) {
                EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
                EXPECT_EQ(truth.total, r.contexts->total.value);
                continue;
            }
            ++stopped;
            EXPECT_EQ(StopReason::MAX_CONTEXTS, r.stop->reason);
            const uint8_t t = check_stopped(r, truth);
            levels.insert(t);
            after_its_pass += r.contexts->by_mismatches.at(t).relation == Relation::EXACT;
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
    EXPECT_LT(0u, after_its_pass);
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

// ---------------------------------------------------------------- the thresholds after the
// check of few unchecked candidates

// what a sweep of the thresholds after the check exercised
struct ThresholdSweep {
    // the allowances under which a checked level dropped source dummies before a later pass
    size_t stale_windows = 0;
    size_t releases = 0;
    size_t stops = 0;
    // stops at the end of a pass, its level complete: on an even-k wrapped PRIMARY graph the
    // ranges its palindrome scans counted crossed the cap
    size_t pass_end_stops = 0;
};

/**
 * A graph without the mask: the check after each pass turns the real candidates of its level
 * into exact contexts and drops the source dummies, and the passes after it compare the
 * counts as the check left them (SPEC §7.4, §21.6). For |text| on |unmasked| under |base|
 * (no threshold in the way), |expected| the masked list and |truth| its counts: every
 * allowance (Request::max_checked_entries) at, below and above each level's cumulative
 * unchecked candidates, and for each the caps at every decision boundary (the final upper
 * bound U, what a bound still carrying the dropped dummies would state, and per pass the
 * completed levels plus its candidates). ALL_OR_COUNT releases every context iff U is within
 * max_contexts, never failing an admitted release (its ranges are retained);
 * stop_at_threshold stops in the pass of level t iff the levels below t as completed plus
 * level t's candidates exceed the cap, and only then, the levels below t as the check left
 * them. |share|: the counted candidates one checked entry stands for (candidates_per_entry).
 * |palindrome_scans| (an even-k wrapped PRIMARY graph): the discovery of level t compares its
 * candidates without the ranges whose palindrome scan is deferred, which enter the running
 * bound scanned at the end of the pass, each palindrome once (SPEC §21.6): the stop comes in
 * that pass or a later one, but at the latest at the end of the first pass but the last whose
 * levels as completed exceed the cap, and a stop at the end of a pass leaves its level as
 * completed.
 */
void sweep_thresholds_after_the_check(const DeBruijnGraph &unmasked, const std::string &text,
                                      const Request &base, const std::vector<MCtx> &expected,
                                      const Truth &truth, ThresholdSweep *sweep,
                                      uint64_t share = 1, bool palindrome_scans = false) {
    const uint8_t m = base.max_mismatches;
    const size_t k = unmasked.get_k();
    // per level: its candidates U_t with no check, and the entries the check tests for those
    // of them unchecked: k - 1 steps each, what its pass charges with every entry tested above
    // what it charges with none
    Request none = base;
    none.max_checked_entries = 0;
    Request every = base;
    every.max_checked_entries = 1'000'000'000;
    const Result plain = count_of(unmasked, text, none);
    const Result checked = count_of(unmasked, text, every);
    ASSERT_FALSE(plain.stop);
    ASSERT_FALSE(checked.stop);
    std::vector<uint64_t> upper(m + 1, 0);
    std::vector<uint64_t> unchecked(m + 1, 0);
    for (uint8_t t = 0; t <= m; ++t) {
        const Count &level = plain.contexts->by_mismatches.at(t);
        EXPECT_TRUE(true_relation(level, truth.by_level[t])) << "level " << int(t);
        upper[t] = level.relation == Relation::BOUNDS ? level.upper : level.value;
        const uint64_t check_steps = checked.work.steps_by_mismatches.at(t)
                                        - plain.work.steps_by_mismatches.at(t);
        EXPECT_EQ(0u, check_steps % (k - 1)) << "level " << int(t);
        unchecked[t] = check_steps / (k - 1);
        // the level's bounds differ by its entries, each counted once per orientation counting
        // it; on an even-k wrapped PRIMARY graph a palindrome only the counting search spells
        // exactly can raise the lower bound of the union instead (§7.4: at least each part)
        if (!palindrome_scans) {
            EXPECT_EQ(level.relation == Relation::BOUNDS ? level.upper - level.lower : 0,
                      unchecked[t] * share) << "level " << int(t);
        }
    }
    std::set<uint64_t> limits { 0 };
    uint64_t cumulative = 0;
    for (uint64_t u : unchecked) {
        cumulative += u;
        limits.insert({ cumulative, cumulative + 1 });
        if (cumulative)
            limits.insert(cumulative - 1);
    }
    auto bound_of = [](const Count &count) {
        return count.relation == Relation::EXACT ? count.value : count.upper;
    };
    for (uint64_t limit : limits) {
        SCOPED_TRACE("allowance " + std::to_string(limit));
        // the counts the check leaves: level t is exact iff it has no unchecked candidate or
        // the levels up to t have at most |limit| in all; else U_t
        std::vector<bool> exact(m + 1, false);
        std::vector<uint64_t> completed(m + 1, 0);
        uint64_t so_far = 0;
        for (uint8_t t = 0; t <= m; ++t) {
            so_far += unchecked[t];
            exact[t] = !unchecked[t] || so_far <= limit;
            completed[t] = exact[t] ? truth.by_level[t] : upper[t];
        }
        const uint64_t final_upper = std::accumulate(completed.begin(), completed.end(),
                                                     uint64_t(0));
        // the dummies the check dropped from a level before the last pass
        uint64_t dropped = 0;
        for (uint8_t t = 0; t < m; ++t) {
            dropped += upper[t] - completed[t];
        }
        sweep->stale_windows += dropped > 0;

        Request request = base;
        request.max_checked_entries = limit;
        const Result counted = count_of(unmasked, text, request);
        ASSERT_FALSE(counted.stop);
        EXPECT_EQ(final_upper, bound_of(counted.contexts->total));
        for (uint8_t t = 0; t <= m; ++t) {
            const Count &level = counted.contexts->by_mismatches.at(t);
            EXPECT_EQ(exact[t] ? Relation::EXACT : Relation::BOUNDS, level.relation)
                << "level " << int(t);
            EXPECT_EQ(completed[t], bound_of(level)) << "level " << int(t);
        }

        std::set<uint64_t> caps { 0 };
        auto around = [&](uint64_t cap) {
            caps.insert({ cap, cap + 1 });
            if (cap)
                caps.insert(cap - 1);
        };
        around(final_upper);
        around(final_upper + dropped);
        uint64_t below = 0;
        uint64_t stale = 0;
        for (uint8_t t = 0; t <= m; ++t) {
            around(below + upper[t]);
            around(stale + upper[t]);
            below += completed[t];
            stale += upper[t];
        }
        for (uint64_t cap : caps) {
            SCOPED_TRACE("cap " + std::to_string(cap));
            request.max_contexts = cap;
            request.stop_at_threshold = false;
            Result r;
            std::vector<MCtx> got;
            ASSERT_NO_THROW(got = released_by(unmasked, text, request, &r));
            ASSERT_TRUE(r.extraction);
            EXPECT_FALSE(r.stop);
            if (final_upper <= cap) {
                ++sweep->releases;
                EXPECT_EQ(expected, got);
                EXPECT_TRUE(r.extraction->complete);
                EXPECT_FALSE(r.extraction->withheld);
                for (uint8_t t = 0; t <= m; ++t) {
                    const Count &level = r.contexts->by_mismatches.at(t);
                    EXPECT_EQ(Relation::EXACT, level.relation) << "level " << int(t);
                    EXPECT_EQ(truth.by_level[t], level.value) << "level " << int(t);
                }
            } else {
                EXPECT_TRUE(got.empty());
                ASSERT_TRUE(r.extraction->withheld);
                EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, *r.extraction->withheld);
            }

            // the stop's pass: the first level whose candidates on top of the completed
            // levels exceed the cap
            request.stop_at_threshold = true;
            const Result s = count_of(unmasked, text, request);
            std::optional<uint8_t> level;
            uint64_t done = 0;
            for (uint8_t t = 0; t <= m; ++t) {
                if (done + upper[t] > cap) {
                    level = t;
                    break;
                }
                done += completed[t];
            }
            // the comparison at the end of each pass but the last: the first whose levels as
            // completed exceed the cap
            std::optional<uint8_t> latest;
            uint64_t through = 0;
            for (uint8_t t = 0; t < m; ++t) {
                through += completed[t];
                if (through > cap) {
                    latest = t;
                    break;
                }
            }
            if (palindrome_scans) {
                ASSERT_TRUE(!latest || s.stop);
                ASSERT_TRUE(level || !s.stop);
            } else {
                ASSERT_EQ(bool(level), bool(s.stop));
            }
            if (!s.stop) {
                EXPECT_EQ(final_upper, bound_of(s.contexts->total));
                continue;
            }
            ++sweep->stops;
            EXPECT_EQ(StopReason::MAX_CONTEXTS, s.stop->reason);
            EXPECT_EQ(StopPhase::DISCOVERY, s.stop->phase);
            const uint8_t stopped = check_stopped(s, truth, false);
            if (palindrome_scans) {
                EXPECT_LE(*level, stopped);
                if (latest) {
                    EXPECT_GE(*latest, stopped);
                }
            } else {
                EXPECT_EQ(*level, stopped);
            }
            for (uint8_t u = 0; u < stopped; ++u) {
                const Count &c = s.contexts->by_mismatches.at(u);
                EXPECT_EQ(exact[u] ? Relation::EXACT : Relation::BOUNDS, c.relation)
                    << "level " << int(u);
                EXPECT_EQ(completed[u], bound_of(c)) << "level " << int(u);
            }
            const Count &at = s.contexts->by_mismatches.at(stopped);
            if (at.relation != Relation::AT_LEAST) {
                // stopped at the end of its pass: the level as completed, above the cap with
                // the levels below it
                ++sweep->pass_end_stops;
                EXPECT_EQ(exact[stopped] ? Relation::EXACT : Relation::BOUNDS, at.relation);
                EXPECT_EQ(completed[stopped], bound_of(at));
                uint64_t up_to = 0;
                for (uint8_t u = 0; u <= stopped; ++u) {
                    up_to += completed[u];
                }
                EXPECT_LT(cap, up_to);
            }
            Result w;
            EXPECT_TRUE(released_by(unmasked, text, request, &w).empty());
            ASSERT_TRUE(w.extraction && w.extraction->withheld);
            EXPECT_EQ(Withheld::THRESHOLD_CROSSED, *w.extraction->withheld);
        }
    }
}

// The counted candidates one entry the check tests stands for (SPEC §7.4: each range of a base
// search is tested once, and on a wrapped PRIMARY graph an entry can enter both orientations'
// counts). A wrapped PRIMARY orientation searches its window and the window's reverse
// complement, which is the other orientation's window: with both strands, or a palindromic
// pattern (one orientation whose two windows coincide), each entry is counted twice
uint64_t candidates_per_entry(DeBruijnGraph::Mode mode, const std::string &text,
                              Strands strands) {
    return mode == DeBruijnGraph::PRIMARY
            && (palindromic(read_pattern(text)) || strands == Strands::BOTH) ? 2 : 1;
}

TEST(PatternMismatch, UnmaskedThresholdsAfterTheCheck) {
    // random cases without the mask (BASIC, CANONICAL and wrapped PRIMARY at odd and even k;
    // at even k a wrapped PRIMARY graph's palindrome scans enter the running bound at the end
    // of their pass, SPEC §21.6): the thresholds after the check, every allowance and every
    // boundary cap
    ThresholdSweep sweep;
    size_t cases = 0;
    size_t primary = 0;
    size_t even_primary = 0;
    for (uint32_t seed = 8000; seed < 8080; ++seed) {
        Case c = random_case(seed);
        const bool palindrome_scans = c.mode == DeBruijnGraph::PRIMARY && c.k % 2 == 0;
        primary += c.mode == DeBruijnGraph::PRIMARY;
        even_primary += palindrome_scans;
        auto unmasked = unmasked_copy(c);
        SCOPED_TRACE(describe(c));
        const std::vector<MCtx> expected = hamming_walk(*c.graph, c.text, c.request.scope,
                                                        c.request.strands,
                                                        c.request.max_mismatches);
        const Truth truth = truth_of(expected, c.request.max_mismatches, c.text.size(), c.k,
                                     c.request, c.text);
        sweep_thresholds_after_the_check(*unmasked, c.text, c.request, expected, truth, &sweep,
                                         candidates_per_entry(c.mode, c.text,
                                                              c.request.strands),
                                         palindrome_scans);
        ++cases;
    }
    EXPECT_LT(40u, cases);
    EXPECT_LT(5u, primary);
    EXPECT_LT(3u, even_primary);
    EXPECT_LT(10u, sweep.pass_end_stops);
    EXPECT_LT(10u, sweep.stale_windows);
    EXPECT_LT(50u, sweep.releases);
    EXPECT_LT(50u, sweep.stops);
}

// The mini index (scripts/traversal/build_mini_refseq.sh) as the MiniRefSeq tests find it:
// mini_refseq/ in the working directory (build/), built without its .edgemask
TEST(PatternMismatch, UnmaskedMiniRefSeqThresholdsAfterTheCheck) {
    // the first 16 bases of an E. coli record (NZ_CP021206.1) whose first k-mer no k-mer
    // enters: the graph's source dummies hold them at offsets 1 to 15 (17 contexts with the
    // mask, 32 candidates without, pattern_fixtures.py's START16), so that pass 0's check
    // drops 15 dummies that no later pass's threshold may carry
    const std::string path = "mini_refseq/graph_k31.dbg";
    if (!std::ifstream(path).good()) {
        GTEST_SKIP() << "mini_refseq index not found in the working directory (build it with "
                        "scripts/traversal/build_mini_refseq.sh)";
    }
    auto graph = std::make_shared<DBGSuccinct>(2);
    ASSERT_TRUE(graph->load(path));
    graph->reset_mask();
    const std::string text = "GATGCCGGTGAACAAC";
    for (uint8_t m = 1; m <= 3; ++m) {
        SCOPED_TRACE("m " + std::to_string(m));
        const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, m);
        // the masked list: every candidate enumerated, the source dummies dropped
        Result all;
        const std::vector<MCtx> expected = released_by(*graph, text, request, &all);
        ASSERT_FALSE(all.stop);
        ASSERT_TRUE(all.extraction && all.extraction->complete);
        const Truth truth = truth_of(expected, m, text.size(), graph->get_k(), request, text);
        EXPECT_EQ(17u, truth.by_level[0]);
        ThresholdSweep sweep;
        sweep_thresholds_after_the_check(*graph, text, request, expected, truth, &sweep);
        EXPECT_LT(0u, sweep.stale_windows);
        EXPECT_LT(0u, sweep.releases);
        EXPECT_LT(0u, sweep.stops);
    }
}

// a k-mer holding the DNA pattern |q| at |p| can be its own reverse complement with |t|
// substitutions only where at most t mirrored pairs of pattern positions are not complementary
// (even k; no odd-length k-mer is its own reverse complement)
bool palindrome_possible(const std::string &q, size_t p, size_t k, uint8_t t) {
    if (k % 2) {
        return false;
    }
    size_t conflicts = 0;
    for (size_t i = 0; i < q.size(); ++i) {
        const size_t mirror = k - 1 - (p + i);
        if (mirror > p + i && mirror < p + q.size()
                && q[mirror - p] != complement_base(q[i])) {
            ++conflicts;
        }
    }
    return conflicts <= t;
}

// a wrapped PRIMARY graph whose stored k-mers are those of |contigs| as given (none the reverse
// complement of another), with the dummy-edge mask or without it
std::shared_ptr<DeBruijnGraph> primary_of(size_t k, const std::vector<std::string> &contigs,
                                          bool masked) {
    auto graph = std::make_shared<DBGSuccinct>(k, DeBruijnGraph::PRIMARY);
    for (const std::string &contig : contigs) {
        graph->add_sequence(contig);
    }
    graph->mask_dummy_kmers(1, false);
    if (!masked) {
        graph->reset_mask();
    }
    return std::make_shared<CanonicalDBG>(std::static_pointer_cast<DeBruijnGraph>(graph), 2);
}

TEST(PatternMismatch, EvenKPrimaryUnmaskedStopsAtTheEndOfThePass) {
    // k = 10, a wrapped PRIMARY graph without its mask, L = 8, strands forward: ACGAATTC holds
    // GAATTC, its own reverse complement, at 2, so a k-mer holding it at offset 0 can be
    // palindromic at every level, and the direct search's ranges there carry a palindrome scan,
    // deferred to the end of their pass. The comparison of stop_at_threshold at the end of
    // pass t (SPEC §21.6) counts those ranges as the scans counted them, each palindrome once:
    // a cap C between the total of the levels 0 .. t without them and with them stops there
    // ({discovery, max_contexts}, mismatch_level t, the levels 0 .. t exact, the later ones
    // unknown), not in a later pass
    const size_t k = 10;
    const std::string text = "ACGAATTC";
    // stored as given, between random flanks: the pattern at offset 0 at levels 0 and 1 (each
    // with a palindromic k-mer, ACGAATTCGT and TCGAATTCGA), and TGGAATTCGT, whose reverse
    // complement ACGAATTCCA only the reverse-complement search finds
    std::mt19937 rng(23);
    std::vector<std::string> contigs;
    for (const char *planted : { "ACGAATTCAA", "ACGAATTCTC", "ACGAATTCCC", "ACGAATTCAG",
                                 "ACGAATTCGT", "TGGAATTCGT", "ACGTATTCAA", "ACCAATTCCA",
                                 "AGGAATTCTT", "ACGAATACAG", "TCGAATTCGA" }) {
        const std::vector<std::string> flanks = random_records(rng, 2, 12, 12);
        contigs.push_back(flanks[0] + planted + flanks[1]);
    }
    std::set<std::string> kmers;
    for (const std::string &contig : contigs) {
        for (size_t i = 0; i + k <= contig.size(); ++i) {
            kmers.insert(contig.substr(i, k));
        }
    }
    for (const std::string &kmer : kmers) {
        const std::string rc = reverse_complement_string(kmer);
        ASSERT_TRUE(rc == kmer || !kmers.count(rc)) << kmer;
    }
    auto masked = primary_of(k, contigs, true);
    auto unmasked = primary_of(k, contigs, false);
    const node_index stored = stored_graph(*masked).max_index();

    for (uint8_t m : { 1, 2 }) {
        SCOPED_TRACE("m " + std::to_string(m));
        Request request = make_request(Scope::ANY_OFFSET, Strands::FORWARD, m);
        request.max_checked_entries = 1'000'000;
        const std::vector<MCtx> expected = hamming_walk(*masked, text, request.scope,
                                                        request.strands, m);
        const Truth truth = truth_of(expected, m, text.size(), k, request, text);
        const Result full = count_of(*unmasked, text, request);
        ASSERT_FALSE(full.stop);
        Request unchecked = request;
        unchecked.max_checked_entries = 0;
        const Result plain = count_of(*unmasked, text, unchecked);
        ASSERT_FALSE(plain.stop);
        // per level: the contexts of the ranges the palindrome scans count less their
        // palindromes (the stored k-mers holding the pattern at an offset where it can be
        // palindromic, each not its own reverse complement), and the source dummies among the
        // candidates
        std::vector<uint64_t> scanned(m + 1, 0);
        for (const MCtx &c : expected) {
            const std::string kmer = masked->get_node_sequence(c.node);
            scanned[c.level] += c.node <= stored
                && palindrome_possible(text, c.offset, k, c.level)
                && kmer != reverse_complement_string(kmer);
        }
        std::vector<uint64_t> dummies(m + 1, 0);
        for (uint8_t t = 0; t <= m; ++t) {
            ASSERT_EQ(Relation::EXACT, full.contexts->by_mismatches.at(t).relation);
            ASSERT_EQ(truth.by_level[t], full.contexts->by_mismatches.at(t).value);
            const Count &level = plain.contexts->by_mismatches.at(t);
            dummies[t] = (level.relation == Relation::BOUNDS ? level.upper : level.value)
                            - truth.by_level[t];
        }
        // the last pass compared at its end; C just below the levels 0 .. t, above them
        // without the scanned ranges, and above what the discovery of level t compares (the
        // levels below t as completed, level t's candidates without the scanned ranges)
        const uint8_t t = m - 1;
        uint64_t with = 0;
        uint64_t without = 0;
        for (uint8_t u = 0; u <= t; ++u) {
            with += truth.by_level[u];
            without += truth.by_level[u] - scanned[u];
        }
        const uint64_t cap = with - 1;
        ASSERT_LE(without, cap);
        ASSERT_LE(with - scanned[t] + dummies[t], cap);

        request.stop_at_threshold = true;
        request.max_contexts = cap;
        const Result s = count_of(*unmasked, text, request);
        ASSERT_TRUE(s.stop);
        EXPECT_EQ(StopPhase::DISCOVERY, s.stop->phase);
        EXPECT_EQ(StopReason::MAX_CONTEXTS, s.stop->reason);
        EXPECT_EQ(t, s.stop->mismatch_level.value_or(kMaxMismatches + 1));
        EXPECT_EQ(t, check_stopped(s, truth));
        for (uint8_t u = 0; u <= m; ++u) {
            const Count &level = s.contexts->by_mismatches.at(u);
            if (u <= t) {
                EXPECT_EQ(Relation::EXACT, level.relation) << "level " << int(u);
                EXPECT_EQ(truth.by_level[u], level.value) << "level " << int(u);
            } else {
                EXPECT_EQ(Relation::UNKNOWN, level.relation) << "level " << int(u);
                EXPECT_EQ(0u, s.work.steps_by_mismatches.at(u)) << "level " << int(u);
            }
        }
        Result w;
        EXPECT_TRUE(released_by(*unmasked, text, request, &w).empty());
        ASSERT_TRUE(w.extraction && w.extraction->withheld);
        EXPECT_EQ(Withheld::THRESHOLD_CROSSED, *w.extraction->withheld);

        // at the levels' total nothing is crossed up to the end of pass t
        request.max_contexts = with;
        const Result above = count_of(*unmasked, text, request);
        EXPECT_TRUE(!above.stop || *above.stop->mismatch_level > t);
    }
}

TEST(PatternMismatch, ThresholdStopsAgainstTheUnboundedRun) {
    // stop_at_threshold against the run without a threshold on random cases of every mode, odd
    // and even k, with and without the mask (no check of unchecked candidates: a level's final
    // bound is then what its pass leaves the running bound): a stop in the pass of level t needs
    // the levels 0 .. t, as the unbounded run counts them, above the cap (the running bound
    // never exceeds them), comes at the latest at the end of the first pass but the last whose
    // levels exceed it (SPEC §21.6), and leaves the levels below t as the unbounded run counts
    // them, level t too when it came at the end of the pass; without a stop the counts are the
    // unbounded run's. ALL_OR_COUNT under the caps around the final bound never drops a range
    // that its admitted release needs (the engine throws) and releases O1's list when it fits
    std::map<std::string, size_t> pass_end_stops;
    std::map<std::string, size_t> cases;
    size_t stops = 0;
    for (uint32_t seed = 9000; seed < 9200; ++seed) {
        const Case c = random_case(seed);
        const uint8_t m = c.request.max_mismatches;
        const std::vector<MCtx> expected = hamming_walk(*c.graph, c.text, c.request.scope,
                                                        c.request.strands, m);
        const Truth truth = truth_of(expected, m, c.text.size(), c.k, c.request, c.text);
        for (bool unmasked : { false, true }) {
            const std::string kind = std::string(mode_name(c.mode))
                + (c.k % 2 ? " odd" : " even") + (unmasked ? " unmasked" : " masked");
            SCOPED_TRACE(describe(c) + " " + kind);
            ++cases[kind];
            auto graph = unmasked ? unmasked_copy(c) : c.graph;
            const Result full = count_of(*graph, c.text, c.request);
            ASSERT_FALSE(full.stop);
            auto bound_of = [](const Count &count) {
                return count.relation == Relation::BOUNDS ? count.upper : count.value;
            };
            std::vector<uint64_t> bound(m + 1, 0);
            std::set<uint64_t> caps { 0 };
            uint64_t closer = 0;
            for (uint8_t t = 0; t <= m; ++t) {
                const Count &level = full.contexts->by_mismatches.at(t);
                if (!unmasked) {
                    EXPECT_EQ(Relation::EXACT, level.relation) << "level " << int(t);
                }
                EXPECT_TRUE(true_relation(level, truth.by_level[t])) << "level " << int(t);
                bound[t] = bound_of(level);
                closer += bound[t];
                caps.insert({ closer, closer + 1 });
                if (closer) {
                    caps.insert(closer - 1);
                }
            }
            const uint64_t total = bound_of(full.contexts->total);
            EXPECT_EQ(closer, total);

            for (uint64_t cap : caps) {
                SCOPED_TRACE("cap " + std::to_string(cap));
                Request request = c.request;
                request.stop_at_threshold = true;
                request.max_contexts = cap;
                const Result s = count_of(*graph, c.text, request);
                std::optional<uint8_t> latest;
                uint64_t through = 0;
                for (uint8_t t = 0; t < m; ++t) {
                    through += bound[t];
                    if (through > cap) {
                        latest = t;
                        break;
                    }
                }
                ASSERT_TRUE(!latest || s.stop);
                if (!s.stop) {
                    EXPECT_EQ(full.contexts->total.relation, s.contexts->total.relation);
                    EXPECT_EQ(total, bound_of(s.contexts->total));
                    continue;
                }
                ++stops;
                EXPECT_EQ(StopPhase::DISCOVERY, s.stop->phase);
                EXPECT_EQ(StopReason::MAX_CONTEXTS, s.stop->reason);
                const uint8_t t = check_stopped(s, truth, !unmasked);
                if (latest) {
                    EXPECT_GE(*latest, t);
                }
                uint64_t up_to = 0;
                for (uint8_t u = 0; u <= t; ++u) {
                    up_to += bound[u];
                }
                EXPECT_LT(cap, up_to);
                for (uint8_t u = 0; u < t; ++u) {
                    const Count &level = s.contexts->by_mismatches.at(u);
                    EXPECT_EQ(full.contexts->by_mismatches.at(u).relation, level.relation)
                        << "level " << int(u);
                    EXPECT_EQ(bound[u], bound_of(level)) << "level " << int(u);
                }
                const Count &at = s.contexts->by_mismatches.at(t);
                if (at.relation != Relation::AT_LEAST) {
                    ++pass_end_stops[kind];
                    EXPECT_EQ(full.contexts->by_mismatches.at(t).relation, at.relation);
                    EXPECT_EQ(bound[t], bound_of(at));
                }
            }

            for (uint64_t cap : { total, total + 1, total ? total - 1 : 0 }) {
                SCOPED_TRACE("release cap " + std::to_string(cap));
                Request request = c.request;
                request.max_contexts = cap;
                Result r;
                std::vector<MCtx> got;
                ASSERT_NO_THROW(got = released_by(*graph, c.text, request, &r));
                ASSERT_TRUE(r.extraction);
                if (total <= cap) {
                    EXPECT_EQ(expected, got);
                } else {
                    EXPECT_TRUE(got.empty());
                }
            }
        }
    }
    for (const char *kind : { "basic odd masked", "basic even unmasked", "canonical odd masked",
                              "canonical even unmasked", "primary odd masked",
                              "primary odd unmasked", "primary even masked",
                              "primary even unmasked" }) {
        EXPECT_LT(5u, cases[kind]) << kind;
    }
    EXPECT_LT(1000u, stops);
    // the stops the palindrome scans' ranges decide at the end of a pass
    EXPECT_LT(10u, pass_end_stops["primary even masked"]);
    EXPECT_LT(10u, pass_end_stops["primary even unmasked"]);
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
    // above the contract's 3, and the paths of a pattern longer than k (its extension, 5m-4):
    // refused as arguments before anything is charged (the route checks them before). A
    // peptide, and a long pattern's anchors, are served
    auto graph = build(5, { "ACGTACGTAC" }, DeBruijnGraph::BASIC);
    const PatternSearch engine(*graph);
    Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 4);
    Budget budget = budget_of();
    EXPECT_THROW(engine.count(parse("ACG"), request, budget), std::invalid_argument);
    request.max_mismatches = 1;
    request.extend_paths = true;
    EXPECT_THROW(engine.count(parse("ACGTAC"), request, budget), std::invalid_argument);
    EXPECT_THROW(engine.enumerate(parse("ACGTAC"), request, budget, [](const Context&) {}),
                 std::invalid_argument);
    EXPECT_THROW(engine.count(Pattern::parse(PatternKind::PROTEIN, "MK"), request, budget),
                 std::invalid_argument);
    EXPECT_EQ(0u, budget.steps_used());
    request.extend_paths = false;
    const Result anchors = engine.count(parse("ACGTAC"), request, budget);
    ASSERT_TRUE(anchors.anchors);
    EXPECT_EQ(2u, anchors.anchors->by_mismatches.size());
    EXPECT_EQ(Relation::UNKNOWN, anchors.anchors->paths.relation);
    const Result peptide = engine.count(Pattern::parse(PatternKind::PROTEIN, "T"), request,
                                        budget);
    ASSERT_TRUE(peptide.contexts);
    EXPECT_EQ(2u, peptide.contexts->by_mismatches.size());
    // without mismatches nothing of the new answer is set
    request.max_mismatches = 0;
    const Result r = engine.count(parse("ACG"), request, budget);
    EXPECT_TRUE(r.contexts->by_mismatches.empty());
    EXPECT_TRUE(r.work.steps_by_mismatches.empty());
    EXPECT_FALSE(r.mismatch_information_bits);
    EXPECT_FALSE(r.max_mismatches_within_floor);
}

TEST(PatternMismatch, StopBeforeThePatternStatesEveryLevel) {
    // a request-wide stop before the pattern (an earlier pattern spent the budget): nothing is
    // searched or charged, every level UNKNOWN, and the work still has an entry per level
    auto graph = build(5, { "ACGTACGTTGCAAGGCTTACG" }, DeBruijnGraph::BASIC);
    const PatternSearch engine(*graph);
    Budget budget = budget_of(1);
    Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 0);
    engine.count(parse("ACG"), request, budget);
    ASSERT_TRUE(budget.stopped());
    const uint64_t spent = budget.steps_used();
    request.max_mismatches = 2;
    for (const std::string &text : { std::string("CGT"), std::string("ACGTACG") }) {
        SCOPED_TRACE(text);
        const Result r = engine.count(parse(text), request, budget);
        ASSERT_TRUE(r.stop);
        EXPECT_EQ(StopPhase::DISCOVERY, r.stop->phase);
        EXPECT_EQ(0, r.stop->mismatch_level.value_or(kMaxMismatches + 1));
        const std::vector<Count> &by_level = r.contexts ? r.contexts->by_mismatches
                                                        : r.anchors->by_mismatches;
        ASSERT_EQ(3u, by_level.size());
        for (const Count &level : by_level) {
            EXPECT_EQ(Relation::UNKNOWN, level.relation);
        }
        EXPECT_EQ(std::vector<uint64_t>(3, 0), r.work.steps_by_mismatches);
        EXPECT_EQ(0u, r.work.steps);
    }
    EXPECT_EQ(spent, budget.steps_used());
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


// ---------------------------------------------------------------- 5m-2: peptides

Pattern peptide_of(const std::string &peptide, int table) {
    return Pattern::parse(PatternKind::PROTEIN, peptide, GeneticCode::get(table));
}

/**
 * check_against() for a peptide against O4's graph walk, and the walk against O4's six-frame
 * record scan with |records|. Returns the number of contexts.
 */
size_t check_peptide(const DeBruijnGraph &graph, const std::string &peptide, int table,
                     const Request &request, const std::vector<std::string> *records = nullptr,
                     DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    const uint8_t m = request.max_mismatches;
    const std::vector<MCtx> expected = peptide_walk(graph, peptide, table, request.scope,
                                                    request.strands, m);
    if (records) {
        EXPECT_EQ(peptide_records(*records, graph.get_k(), mode, peptide, table, request.scope,
                                  request.strands, m),
                  spelled(graph, expected)) << peptide << " table " << table;
    }
    return check_against(graph, peptide_of(peptide, table), request, expected,
                         peptide_orientations(peptide, oracle_code(table), request.strands));
}

// a residue letter of a random peptide: the 20 residues, and X, B, Z, J and '*' less often
char random_letter(std::mt19937 &rng) {
    static const std::string kResidues = "ACDEFGHIKLMNPQRSTVWY";
    static const std::string kOthers = "XBZJ*";
    return rng() % 6 ? kResidues[rng() % kResidues.size()] : kOthers[rng() % kOthers.size()];
}

// a codon of |code| whose residue |accept| admits, drawn at random; "" when none does
std::string random_codon(std::mt19937 &rng, const OracleCode &code,
                         const std::function<bool(char)> &accept) {
    std::vector<std::string> codons;
    for (const std::string &codon : all_codons()) {
        if (accept(code.translate(codon)))
            codons.push_back(codon);
    }
    return codons.empty() ? "" : codons[rng() % codons.size()];
}

/**
 * Records around a peptide's neighbourhood: |copies| back-translations whose residues are each,
 * at random, one of its own codons, a sense codon of another residue (a substitution), a stop
 * codon or any codon, between random flanks, a third of them reverse-complemented, and a
 * random record.
 */
std::vector<std::string> neighbourhood_records(std::mt19937 &rng, const std::string &peptide,
                                               int table, size_t copies) {
    const OracleCode &code = oracle_code(table);
    std::vector<std::string> records;
    for (size_t c = 0; c < copies; ++c) {
        std::string record = random_records(rng, 1, 0, 6)[0];
        for (char letter : peptide) {
            const uint32_t roll = rng() % 10;
            std::string codon;
            if (roll < 6) {
                codon = random_codon(rng, code, [&](char r) { return letter_admits(letter, r); });
            } else if (roll < 8) {
                codon = random_codon(rng, code, [&](char r) {
                    return r != '*' && !letter_admits(letter, r);
                });
            } else if (roll < 9) {
                codon = random_codon(rng, code, [](char r) { return r == '*'; });
            }
            record += codon.empty() ? all_codons()[rng() % 64] : codon;
        }
        record += random_records(rng, 1, 0, 6)[0];
        records.push_back(rng() % 3 ? record : reverse_complement_string(record));
    }
    records.push_back(random_records(rng, 1, 20, 60)[0]);
    return records;
}

// one random peptide case: a table, k, a graph mode, a peptide of L = 3r <= k, records around
// its neighbourhood, m in {1, 2} (3 for k <= 9), a scope (any_offset only on PRIMARY), strands
struct PeptideCase {
    size_t k;
    DeBruijnGraph::Mode mode;
    int table;
    std::string peptide;
    std::vector<std::string> records;
    std::shared_ptr<DeBruijnGraph> graph;
    Request request;
};

PeptideCase random_peptide_case(uint32_t seed,
                                const std::vector<size_t> &ks = { 6, 9, 10, 12 }) {
    std::mt19937 rng(seed);
    PeptideCase c;
    c.k = ks[rng() % ks.size()];
    c.mode = std::vector<DeBruijnGraph::Mode> {
        DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[rng() % 3];
    const std::vector<int> ids = oracle_code_ids();
    c.table = ids[rng() % ids.size()];
    const size_t r = 1 + rng() % (c.k / 3);
    for (size_t i = 0; i < r; ++i) {
        c.peptide.push_back(random_letter(rng));
    }
    c.records = neighbourhood_records(rng, c.peptide, c.table, 3 + rng() % 5);
    c.graph = build(c.k, c.records, c.mode);
    const Scope scope = c.mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0 ? Scope::SUFFIX
                                                                           : Scope::ANY_OFFSET;
    const Strands strands = std::vector<Strands> {
        Strands::BOTH, Strands::BOTH, Strands::FORWARD, Strands::REVERSE }[rng() % 4];
    c.request = make_request(scope, strands, 1 + rng() % (c.k <= 9 ? 3 : 2));
    return c;
}

std::string describe(const PeptideCase &c) {
    return c.peptide + " table " + std::to_string(c.table) + " k " + std::to_string(c.k) + " "
         + mode_name(c.mode) + " m " + std::to_string(c.request.max_mismatches) + " "
         + to_string(c.request.scope) + " " + to_string(c.request.strands);
}

TEST(PatternMismatchPeptide, OraclesAgree) {
    // O4 on the graph (codon by codon) = O4 on the records (six frames translated); its level 0
    // = the search without mismatches
    size_t substituted = 0;
    size_t contexts = 0;
    for (uint32_t seed = 100; seed < 250; ++seed) {
        const PeptideCase c = random_peptide_case(seed);
        SCOPED_TRACE(describe(c));
        const uint8_t m = c.request.max_mismatches;
        const std::vector<MCtx> walk = peptide_walk(*c.graph, c.peptide, c.table,
                                                    c.request.scope, c.request.strands, m);
        EXPECT_EQ(peptide_records(c.records, c.k, c.mode, c.peptide, c.table, c.request.scope,
                                  c.request.strands, m),
                  spelled(*c.graph, walk));
        Request exact = c.request;
        exact.max_mismatches = 0;
        Result r;
        const std::vector<MCtx> v1 = released_by(*c.graph, peptide_of(c.peptide, c.table),
                                                 exact, &r);
        std::vector<MCtx> level0;
        for (const MCtx &x : walk) {
            if (!x.level)
                level0.push_back(x);
            substituted += x.level > 0;
        }
        EXPECT_EQ(level0, v1);
        contexts += walk.size();
    }
    EXPECT_LT(300u, substituted);
    EXPECT_LT(500u, contexts);
}

TEST(PatternMismatchPeptide, RandomPeptidesAgainstOracles) {
    // 300 random cases over every genetic code: every count per level, offset and orientation,
    // the release with each context's level, substituted residues and their codons, the
    // thresholds, against O4 (graph and records)
    std::vector<size_t> at_level(kMaxMismatches + 1, 0);
    std::set<int> tables;
    for (uint32_t seed = 300; seed < 600; ++seed) {
        const PeptideCase c = random_peptide_case(seed);
        SCOPED_TRACE(describe(c));
        check_peptide(*c.graph, c.peptide, c.table, c.request, &c.records, c.mode);
        for (const MCtx &x : peptide_walk(*c.graph, c.peptide, c.table, c.request.scope,
                                          c.request.strands, c.request.max_mismatches)) {
            ++at_level[x.level];
            tables.insert(c.table);
        }
    }
    EXPECT_EQ(oracle_code_ids().size(), tables.size());
    EXPECT_LT(500u, at_level[0]);
    EXPECT_LT(500u, at_level[1]);
    EXPECT_LT(200u, at_level[2]);
    EXPECT_LT(20u, at_level[3]);
}

TEST(PatternMismatchPeptide, RandomPeptidesAtK31) {
    // the served k: peptides of 4 to 10 residues in records of a few hundred bases holding
    // their neighbourhood, m 1 or 2, every mode
    size_t substituted = 0;
    for (uint32_t seed = 700; seed < 720; ++seed) {
        std::mt19937 rng(seed);
        const auto mode = std::vector<DeBruijnGraph::Mode> {
            DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[seed % 3];
        const int table = oracle_code_ids()[rng() % oracle_code_ids().size()];
        std::string peptide;
        for (size_t i = 0, r = 4 + rng() % 7; i < r; ++i) {
            peptide.push_back(random_letter(rng));
        }
        std::vector<std::string> records = neighbourhood_records(rng, peptide, table, 8);
        for (const std::string &extra : random_records(rng, 2, 150, 300)) {
            records.push_back(extra);
        }
        auto graph = build(31, records, mode);
        const Request request = make_request(
                mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0 ? Scope::SUFFIX
                                                                 : Scope::ANY_OFFSET,
                Strands::BOTH, 1 + rng() % 2);
        SCOPED_TRACE(peptide + " table " + std::to_string(table) + " " + mode_name(mode)
                     + " m " + std::to_string(request.max_mismatches));
        check_peptide(*graph, peptide, table, request, &records, mode);
        for (const MCtx &x : peptide_walk(*graph, peptide, table, request.scope,
                                          request.strands, request.max_mismatches)) {
            substituted += x.level > 0;
        }
    }
    EXPECT_LT(20u, substituted);
}

TEST(PatternMismatchPeptide, PassZeroIsTheSearchWithoutMismatches) {
    // as for DNA: pass 0 charges the steps of the peptide search without mismatches in the same
    // order, under every max_steps that stops it, with and without the mask
    size_t stops = 0;
    for (uint32_t seed = 800; seed < 840; ++seed) {
        const PeptideCase c = random_peptide_case(seed);
        const bool unmasked = seed % 3 == 0;
        auto graph = c.graph;
        if (unmasked) {
            graph = build(c.k, c.records, c.mode);
            const_cast<DBGSuccinct&>(stored_graph(*graph)).reset_mask();
        }
        SCOPED_TRACE(describe(c) + (unmasked ? " unmasked" : ""));
        const Pattern pattern = peptide_of(c.peptide, c.table);
        Request with = c.request;
        with.max_checked_entries = unmasked ? 3 : 0;
        Request without = with;
        without.max_mismatches = 0;
        const Result v1 = count_of(*graph, pattern, without);
        const Result full = count_of(*graph, pattern, with);
        ASSERT_FALSE(v1.stop);
        ASSERT_FALSE(full.stop);
        const uint64_t steps = v1.work.steps;
        EXPECT_EQ(steps, full.work.steps_by_mismatches.at(0));
        EXPECT_EQ(v1.contexts->total.value, full.contexts->by_mismatches.at(0).value);
        EXPECT_EQ(v1.contexts->total.relation, full.contexts->by_mismatches.at(0).relation);
        for (uint64_t s = 0; s < steps; ++s) {
            const Result a = count_of(*graph, pattern, without, s);
            const Result b = count_of(*graph, pattern, with, s);
            ASSERT_TRUE(a.stop && b.stop) << s;
            ++stops;
            EXPECT_EQ(0, b.stop->mismatch_level.value_or(9)) << s;
            EXPECT_EQ(a.stop->phase, b.stop->phase) << s;
            EXPECT_EQ(a.work.ranges_visited, b.work.ranges_visited) << s;
            EXPECT_EQ(a.work.steps, b.work.steps) << s;
            EXPECT_EQ(a.contexts->total.relation, b.contexts->by_mismatches.at(0).relation);
            EXPECT_EQ(a.contexts->total.value, b.contexts->by_mismatches.at(0).value) << s;
            Request partial_a = without;
            Request partial_b = with;
            partial_a.mode = partial_b.mode = Mode::PARTIAL;
            Result ra;
            Result rb;
            EXPECT_EQ(released_by(*graph, pattern, partial_a, &ra, s),
                      released_by(*graph, pattern, partial_b, &rb, s)) << s;
        }
    }
    EXPECT_LT(1000u, stops);
}

TEST(PatternMismatchPeptide, StoppedAtEveryStep) {
    // every max_steps below a peptide case's total: the levels below the stop's complete, its
    // own AT_LEAST or BOUNDS, the ones above UNKNOWN; PARTIAL a prefix of the closest levels
    std::set<uint8_t> levels;
    size_t cases = 0;
    for (uint32_t seed = 900; cases < 16 && seed < 2000; ++seed) {
        const PeptideCase c = random_peptide_case(seed, { 9, 12 });
        const Pattern pattern = peptide_of(c.peptide, c.table);
        const Result full = count_of(*c.graph, pattern, c.request);
        if (full.work.steps < 50 || full.work.steps > 20'000
                || !full.contexts->by_mismatches.back().value)
            continue;
        ++cases;
        SCOPED_TRACE(describe(c));
        const std::vector<MCtx> expected = peptide_walk(*c.graph, c.peptide, c.table,
                                                        c.request.scope, c.request.strands,
                                                        c.request.max_mismatches);
        const Truth truth = truth_of(expected, c.request.max_mismatches,
                                     scope_offsets(pattern.length(), c.k, c.request.scope),
                                     full.searched);
        const uint64_t stride = std::max<uint64_t>(1, full.work.steps / 300);
        for (uint64_t s = 0; s < full.work.steps; s += (s < 200 ? 1 : stride)) {
            const Result r = count_of(*c.graph, pattern, c.request, s);
            ASSERT_TRUE(r.stop) << s;
            const uint8_t t = check_stopped(r, truth);
            EXPECT_EQ(pass_of_step(full.work.steps_by_mismatches, s), t) << s;
            levels.insert(t);
            Request partial = c.request;
            partial.mode = Mode::PARTIAL;
            Result p;
            check_partial_after_stop(released_by(*c.graph, pattern, partial, &p, s), expected,
                                     t);
        }
    }
    EXPECT_EQ(16u, cases);
    EXPECT_LE(3u, levels.size());
}

TEST(PatternMismatchPeptide, CheaperOrientationFirstInEveryPass) {
    // the leucines that start LLLLMW branch widely where its reverse complement, which starts
    // with rc(W) and rc(M), is narrow: in pass 1 the reverse orientation runs first (the
    // residue-level estimate), so that a step stop there never leaves the forward one complete
    // before it
    std::mt19937 rng(73);
    auto graph = build(19, random_records(rng, 6, 800, 800), DeBruijnGraph::BASIC);
    const Pattern pattern = peptide_of("LLLLMW", 1);
    const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, 1, Mode::COUNT);
    const Result full = count_of(*graph, pattern, request);
    ASSERT_FALSE(full.stop);
    const uint64_t pass0 = full.work.steps_by_mismatches.at(0);
    const uint64_t stride = std::max<uint64_t>(1, (full.work.steps - pass0) / 2000);
    size_t reverse_only = 0;
    for (uint64_t s = pass0; s < full.work.steps; s += stride) {
        const Result r = count_of(*graph, pattern, request, s);
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

// the level of the context at |kmer| (offset 0, FORWARD) of |peptide| within 1 substitution,
// -1 when it is none
int level_at(const DeBruijnGraph &graph, const std::string &peptide, int table,
             const std::string &kmer, uint8_t m = 1) {
    Result r;
    for (const MCtx &c : released_by(graph, peptide_of(peptide, table),
                                     make_request(Scope::ANY_OFFSET, Strands::FORWARD, m), &r)) {
        if (!c.offset && graph.get_node_sequence(c.node) == kmer)
            return c.level;
    }
    return -1;
}

TEST(PatternMismatchPeptide, SynonymousAndSubstitutedCodons) {
    // a synonymous codon keeps the level, a codon of another residue adds one: Leu TTA -> CTA
    // 0, -> TTT (Phe) 1; Ser TCA -> AGT 0; Met ATG -> ATA (Ile) 1 in table 1, 0 in table 2
    // (where ATA is Met); one substitution per residue however many bases differ
    const std::vector<std::string> records { "TTA", "CTA", "TTT", "TCA", "AGT", "ATG", "ATA",
                                             "GGG" };
    auto graph = build(3, records, DeBruijnGraph::BASIC);
    EXPECT_EQ(0, level_at(*graph, "L", 1, "TTA"));
    EXPECT_EQ(0, level_at(*graph, "L", 1, "CTA"));
    EXPECT_EQ(1, level_at(*graph, "L", 1, "TTT"));
    EXPECT_EQ(0, level_at(*graph, "S", 1, "TCA"));
    EXPECT_EQ(0, level_at(*graph, "S", 1, "AGT"));
    EXPECT_EQ(1, level_at(*graph, "S", 1, "GGG"));
    EXPECT_EQ(0, level_at(*graph, "M", 1, "ATG"));
    EXPECT_EQ(1, level_at(*graph, "M", 1, "ATA"));
    EXPECT_EQ(0, level_at(*graph, "M", 2, "ATA"));
    EXPECT_EQ(-1, level_at(*graph, "M", 1, "ATA", 0));
    // three bases differ, one substitution
    EXPECT_EQ(1, level_at(*graph, "M", 1, "GGG"));
    const auto subs = residue_substitutions(peptide_of("LM", 1), Orientation::FORWARD, "TTTATA");
    ASSERT_EQ(2u, subs.size());
    EXPECT_EQ(0u, subs[0].residue);
    EXPECT_EQ('L', subs[0].expected);
    EXPECT_EQ('F', subs[0].found);
    EXPECT_EQ("TTT", subs[0].codon);
    EXPECT_EQ(1u, subs[1].residue);
    EXPECT_EQ('I', subs[1].found);
    EXPECT_EQ("ATA", subs[1].codon);
    for (const std::string peptide : { "L", "S", "M", "J", "X", "*" }) {
        for (int table : { 1, 2 }) {
            for (uint8_t m = 1; m <= 3; ++m) {
                SCOPED_TRACE(peptide + " " + std::to_string(table) + " " + std::to_string(m));
                check_peptide(*graph, peptide, table,
                              make_request(Scope::ANY_OFFSET, Strands::BOTH, m), &records);
            }
        }
    }
}

TEST(PatternMismatchPeptide, StopsAndXAreNeverSubstituted) {
    // '*' is never substituted (a sense codon there is no instance) and a stop codon never
    // substitutes a residue; X has no substitute; in tables 27, 28 and 31 the codons that end
    // translation only in context are sense codons and substitute other residues
    const std::vector<std::string> records { "AAATAA", "AAATGG", "AAGTGA", "CCCTAG", "AAATGT",
                                             "TGGAAA", "TAAAAA", "CCCGGG", "AAATAG" };
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    // K*: AAATAA and AAGTGA are level 0, CCCTAG level 1 (K -> P), AAATGG none ('*' -> W)
    EXPECT_EQ(0, level_at(*graph, "K*", 1, "AAATAA"));
    EXPECT_EQ(0, level_at(*graph, "K*", 1, "AAGTGA"));
    EXPECT_EQ(1, level_at(*graph, "K*", 1, "CCCTAG"));
    EXPECT_EQ(-1, level_at(*graph, "K*", 1, "AAATGG", 3));
    // KW: TGA is a stop in table 1 (no instance at any level), W in table 4; TGT is C
    EXPECT_EQ(0, level_at(*graph, "KW", 1, "AAATGG"));
    EXPECT_EQ(-1, level_at(*graph, "KW", 1, "AAGTGA", 3));
    EXPECT_EQ(0, level_at(*graph, "KW", 4, "AAGTGA"));
    EXPECT_EQ(1, level_at(*graph, "KW", 1, "AAATGT"));
    // XK: X admits every sense codon and nothing else
    EXPECT_EQ(0, level_at(*graph, "XK", 1, "TGGAAA"));
    EXPECT_EQ(-1, level_at(*graph, "XK", 1, "TAAAAA", 3));
    EXPECT_EQ(1, level_at(*graph, "XK", 1, "CCCGGG"));
    // context stops: KC in table 27 (TGA = W) one substitution at AAGTGA, none in table 1;
    // KE in table 28 (TAA = Q) and KQ in table 31 (TAA = E) one at AAATAA
    EXPECT_EQ(1, level_at(*graph, "KC", 27, "AAGTGA"));
    EXPECT_EQ(-1, level_at(*graph, "KC", 1, "AAGTGA", 3));
    EXPECT_EQ(1, level_at(*graph, "KE", 28, "AAATAA"));
    EXPECT_EQ(1, level_at(*graph, "KQ", 31, "AAATAA"));
    EXPECT_EQ(0, level_at(*graph, "KQ", 28, "AAATAA"));
    for (const std::string peptide : { "K*", "KW", "XK", "KC", "KE", "KQ", "*X", "BZ" }) {
        for (int table : { 1, 4, 27, 28, 31 }) {
            for (uint8_t m = 1; m <= 2; ++m) {
                SCOPED_TRACE(peptide + " " + std::to_string(table) + " " + std::to_string(m));
                check_peptide(*graph, peptide, table,
                              make_request(Scope::ANY_OFFSET, Strands::BOTH, m), &records);
            }
        }
    }
    // '*' read in a table without an unconditional stop has no codon: no instance within any
    // number of substitutions, every level EXACT 0 with nothing searched
    for (int table : { 27, 28, 31 }) {
        for (uint8_t m = 1; m <= 3; ++m) {
            const Pattern stop = peptide_of("K*", table);
            ASSERT_FALSE(stop.has_instances());
            const Result r = count_of(*graph, stop, make_request(Scope::ANY_OFFSET,
                                                                 Strands::BOTH, m));
            ASSERT_TRUE(r.contexts);
            ASSERT_EQ(m + 1u, r.contexts->by_mismatches.size());
            for (const Count &level : r.contexts->by_mismatches) {
                EXPECT_EQ(Relation::EXACT, level.relation);
                EXPECT_EQ(0u, level.value);
            }
            EXPECT_EQ(std::vector<uint64_t>(m + 1, 0), r.work.steps_by_mismatches);
            EXPECT_EQ(0u, r.work.steps);
            EXPECT_TRUE(std::count(r.notes.begin(), r.notes.end(), kNoteNoStopCodon));
            ASSERT_TRUE(r.mismatch_information_bits);
            Result e;
            EXPECT_TRUE(released_by(*graph, stop, make_request(Scope::ANY_OFFSET, Strands::BOTH,
                                                               m), &e).empty());
            EXPECT_TRUE(e.extraction->complete);
        }
    }
}

TEST(PatternMismatchPeptide, ReverseComplementAndBothStrands) {
    // the substitutes of rc(P) are P's reverse-complemented, residues in reverse order; a
    // record holding the reverse complement of a variant is a REVERSE context whose substituted
    // residue and codon read in P's coordinates and reading direction
    std::mt19937 rng(61);
    for (int table : oracle_code_ids()) {
        const OracleCode &code = oracle_code(table);
        std::string peptide;
        for (int i = 0; i < 6; ++i) {
            peptide.push_back(random_letter(rng));
        }
        const Pattern p = peptide_of(peptide, table);
        const Pattern rc = p.reverse_complement();
        ASSERT_EQ(peptide.size(), rc.substitute_sets().size());
        for (size_t i = 0; i < peptide.size(); ++i) {
            // the oracle's substitutes of residue m - 1 - i, reverse-complemented
            std::set<std::string> expected;
            const char letter = peptide[peptide.size() - 1 - i];
            for (const std::string &codon : all_codons()) {
                if (codon_level(letter, code.translate(codon)).value_or(0) == 1)
                    expected.insert(reverse_complement_string(codon));
            }
            std::set<std::string> got;
            for (int c = 0; c < 64; ++c) {
                if (rc.substitute_sets()[i] >> c & 1)
                    got.insert(mtg::test::codon_string(c));
            }
            EXPECT_EQ(expected, got) << peptide << " " << table << " " << i;
        }
    }
    // KDVK (table 1) with D -> E (GAA) on the reverse strand only
    const std::string variant = "AAAGAAGTTAAG";
    const std::vector<std::string> records { "CC" + reverse_complement_string(variant) + "G" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(13, records, mode);
        Result r;
        const std::vector<MCtx> contexts = released_by(
                *graph, peptide_of("KDVK", 1), make_request(Scope::ANY_OFFSET, Strands::BOTH, 1),
                &r);
        size_t reverse = 0;
        for (const MCtx &c : contexts) {
            EXPECT_EQ(1, c.level);
            EXPECT_EQ(std::vector<uint32_t> { 1 }, c.positions);
            reverse += c.orientation == Orientation::REVERSE;
        }
        EXPECT_LT(0u, reverse) << mode_name(mode);
        if (mode == DeBruijnGraph::BASIC) {
            EXPECT_EQ(contexts.size(), reverse);
        }
        check_peptide(*graph, "KDVK", 1, make_request(Scope::ANY_OFFSET, Strands::BOTH, 1),
                      &records, mode);
        const std::string instance = reverse_complement_string(variant);
        const auto subs = residue_substitutions(peptide_of("KDVK", 1), Orientation::REVERSE,
                                                instance);
        ASSERT_EQ(1u, subs.size());
        EXPECT_EQ("GAA", subs[0].codon);
        EXPECT_EQ('E', subs[0].found);
        EXPECT_EQ('D', subs[0].expected);
    }
}

// the bits left within |m| by a peptide's brute-force ball over |n| positions
double peptide_ball_bits(const std::vector<uint64_t> &histogram, size_t n, uint8_t m) {
    uint64_t ball = 0;
    for (size_t u = 0; u <= m && u < histogram.size(); ++u) {
        ball += histogram[u];
    }
    return 2.0 * static_cast<double>(n) - std::log2(static_cast<double>(ball));
}

TEST(PatternMismatchPeptide, BallBitsAgainstTheBruteForceBall) {
    // mismatch_information_bits of a peptide against every string enumerated (O5 by O4's
    // rule), for every table and every letter, on whole codons and on windows that cut them
    std::mt19937 rng(67);
    const std::string kLetters = "ACDEFGHIKLMNPQRSTVWYXBZJ*";
    size_t compared = 0;
    for (int table : oracle_code_ids()) {
        auto compare = [&](const std::string &peptide, size_t begin, size_t end) {
            const Pattern pattern = peptide_of(peptide, table);
            const std::vector<uint64_t> histogram = peptide_ball_histogram(peptide, table,
                                                                          begin, end);
            if (!std::accumulate(histogram.begin(), histogram.end(), uint64_t(0))) {
                // a window without any instance ('*' without a stop codon): the engine's
                // convention, one exact codon with no substitute, keeps it finite (no ball to
                // compare with; alone, the '*' has the bits of every string being exact)
                if (peptide == "*") {
                    EXPECT_DOUBLE_EQ(2.0 * (end - begin),
                                     pattern.mismatch_information_bits(begin, end, 2));
                }
                return;
            }
            for (uint8_t m = 0; m <= kMaxMismatches; ++m) {
                const double expected = peptide_ball_bits(histogram, end - begin, m);
                EXPECT_NEAR(expected, pattern.mismatch_information_bits(begin, end, m), 1e-9)
                    << peptide << " table " << table << " [" << begin << ", " << end
                    << ") m " << int(m);
            }
            ++compared;
        };
        for (char letter : kLetters) {
            for (size_t begin = 0; begin < 3; ++begin) {
                for (size_t end = begin + 1; end <= 3; ++end) {
                    compare(std::string(1, letter), begin, end);
                }
            }
        }
        for (int i = 0; i < 3; ++i) {
            std::string peptide;
            for (int j = 0; j < 3; ++j) {
                peptide.push_back(kLetters[rng() % kLetters.size()]);
            }
            // windows cutting the first and the last codon, and two whole codons
            compare(peptide, 1, 8);
            compare(peptide, 2, 7);
            compare(peptide, 0, 6);
        }
    }
    EXPECT_LT(4000u, compared);
    // the decapeptide of GroEL in table 1 loses 7.8, 14.5 and 20.4 bits (SPEC §21.5); a Leu,
    // Ser, Arg-rich one 6.5, 11.9 and 16.6
    const Pattern groel = peptide_of("KDVKFGNDAR", 1);
    EXPECT_NEAR(45.415, groel.mismatch_information_bits(0), 0.001);
    EXPECT_NEAR(37.570, groel.mismatch_information_bits(1), 0.001);
    EXPECT_NEAR(30.896, groel.mismatch_information_bits(2), 0.001);
    EXPECT_NEAR(24.998, groel.mismatch_information_bits(3), 0.001);
    const Pattern lsr = peptide_of("LSRLSRLSRL", 1);
    EXPECT_NEAR(27.616, lsr.mismatch_information_bits(1), 0.001);
    EXPECT_NEAR(22.231, lsr.mismatch_information_bits(2), 0.001);
    EXPECT_NEAR(17.595, lsr.mismatch_information_bits(3), 0.001);
    // the reverse complement has the same ball, window by mirrored window
    const Pattern rc = groel.reverse_complement();
    for (uint8_t m = 0; m <= 3; ++m) {
        EXPECT_NEAR(groel.mismatch_information_bits(4, 25, m),
                    rc.mismatch_information_bits(5, 26, m), 1e-9);
    }
}

TEST(PatternMismatchPeptide, FloorOnTheResidueBall) {
    // the floor (24 bits) gates the residue-level ball's bits, the refusal naming the largest
    // m' the floor admits; the ball counted here from the oracle's code (per residue: its
    // codons and the sense codons of other residues)
    std::mt19937 rng(71);
    auto graph = build(31, random_records(rng, 3, 200, 200), DeBruijnGraph::BASIC);
    auto bits = [](const std::string &peptide, uint8_t m) {
        const OracleCode &code = oracle_code(1);
        std::vector<double> f(m + 1, 0);
        f[0] = 1;
        for (char letter : peptide) {
            double own = 0;
            double other = 0;
            for (const std::string &codon : all_codons()) {
                const auto level = codon_level(letter, code.translate(codon));
                own += level == std::optional<uint8_t>(0);
                other += level == std::optional<uint8_t>(1);
            }
            for (size_t u = f.size(); u-- > 0; ) {
                f[u] = f[u] * own + (u ? f[u - 1] * other : 0);
            }
        }
        return 6.0 * peptide.size() - std::log2(std::accumulate(f.begin(), f.end(), 0.0));
    };
    for (const std::string peptide : { "KDVKFGNDAR", "LSRLSRLSRL", "MWMWMWMW", "XXXXXMK" }) {
        for (uint8_t m = 1; m <= 3; ++m) {
            Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, m);
            request.min_information_bits = 24;
            const Result r = count_of(*graph, peptide_of(peptide, 1), request);
            ASSERT_TRUE(r.mismatch_information_bits);
            EXPECT_NEAR(bits(peptide, m), *r.mismatch_information_bits, 1e-9) << peptide;
            const bool served = bits(peptide, m) >= 24;
            EXPECT_EQ(served, !r.refusal) << peptide << " m " << int(m);
            if (served)
                continue;
            EXPECT_EQ(0u, r.work.steps);
            std::optional<uint8_t> within;
            for (uint8_t w = m; w-- > 0; ) {
                if (bits(peptide, w) >= 24) {
                    within = w;
                    break;
                }
            }
            EXPECT_EQ(within, r.max_mismatches_within_floor) << peptide << " m " << int(m);
        }
    }
    // a peptide longer than k: the floor reads the least anchor window (reverse: P's last k
    // positions), each window's bits within m
    auto small = build(10, random_records(rng, 3, 100, 100), DeBruijnGraph::BASIC);
    const Pattern longer = peptide_of("WWWWLSRLSR", 1);
    for (uint8_t m = 1; m <= 3; ++m) {
        Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, m);
        request.min_information_bits = 10;
        const Result r = count_of(*small, longer, request);
        const double least = std::min(longer.mismatch_information_bits(0, 10, m),
                                      longer.mismatch_information_bits(20, 30, m));
        ASSERT_TRUE(r.min_anchor_mismatch_information_bits);
        EXPECT_DOUBLE_EQ(least, *r.min_anchor_mismatch_information_bits);
        EXPECT_DOUBLE_EQ(longer.mismatch_information_bits(m), *r.mismatch_information_bits);
        EXPECT_EQ(least < 10, bool(r.refusal)) << int(m);
        if (r.refusal) {
            EXPECT_NE(std::string::npos, r.refusal->message.find("reverse orientation"))
                << r.refusal->message;
        }
    }
}

TEST(PatternMismatchPeptide, CutCodonsInAnchorWindows) {
    // a peptide longer than k: its anchors per level, the anchor window cutting its last codon
    // (k not a multiple of 3) and, on a wrapped PRIMARY graph, the reverse complement's window
    // starting inside one; a cut codon counts as substituted only when its bases inside the
    // window agree with no codon of its own (O4's least level over the agreeing codons)
    std::set<uint8_t> levels;
    size_t anchors = 0;
    for (uint32_t seed = 1100; seed < 1220; ++seed) {
        std::mt19937 rng(seed);
        const size_t k = std::vector<size_t> { 5, 7, 8, 10, 11 }[rng() % 5];
        const auto mode = std::vector<DeBruijnGraph::Mode> {
            DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[rng() % 3];
        const int table = oracle_code_ids()[rng() % oracle_code_ids().size()];
        std::string peptide;
        for (size_t i = 0, r = k / 3 + 1 + rng() % 2; i < r; ++i) {
            peptide.push_back(random_letter(rng));
        }
        const std::vector<std::string> records = neighbourhood_records(rng, peptide, table, 6);
        auto graph = build(k, records, mode);
        const uint8_t m = 1 + rng() % 2;
        const Strands strands = std::vector<Strands> {
            Strands::BOTH, Strands::FORWARD, Strands::REVERSE }[rng() % 3];
        const Pattern pattern = peptide_of(peptide, table);
        Request request = make_request(Scope::ANY_OFFSET, strands, m);
        SCOPED_TRACE(peptide + " table " + std::to_string(table) + " k " + std::to_string(k)
                     + " " + mode_name(mode) + " m " + std::to_string(m));
        std::vector<MCtx> expected = peptide_walk(*graph, peptide, table, request.scope,
                                                  strands, m);
        const Result counted = count_of(*graph, pattern, request);
        ASSERT_TRUE(counted.anchors);
        EXPECT_FALSE(counted.stop);
        EXPECT_EQ(Extension::NOT_REQUESTED, counted.anchors->extension);
        EXPECT_EQ(expected.empty() ? Relation::EXACT : Relation::UNKNOWN,
                  counted.anchors->paths.relation);
        EXPECT_TRUE(std::count(counted.notes.begin(), counted.notes.end(), kNotePathsLater));
        const Truth truth = truth_of(expected, m, { 0 }, counted.searched);
        EXPECT_EQ(Relation::EXACT, counted.anchors->total.relation);
        EXPECT_EQ(truth.total, counted.anchors->total.value);
        ASSERT_EQ(m + 1u, counted.anchors->by_mismatches.size());
        for (uint8_t t = 0; t <= m; ++t) {
            EXPECT_EQ(Relation::EXACT, counted.anchors->by_mismatches[t].relation);
            EXPECT_EQ(truth.by_level[t], counted.anchors->by_mismatches[t].value) << int(t);
        }
        for (const auto &[o, n] : truth.by_orientation) {
            EXPECT_EQ(n, counted.anchors->by_orientation.at(o).value) << orientation_key(o);
        }
        // the anchors released on request: the oracle's list with its levels
        request.release_anchors = true;
        Result released;
        std::vector<MCtx> got = released_by(*graph, pattern, request, &released);
        for (MCtx &c : expected) {
            levels.insert(c.level);
            c.positions.clear();
        }
        EXPECT_EQ(expected, got);
        EXPECT_TRUE(released.extraction->complete);
        anchors += expected.size();
    }
    EXPECT_EQ((std::set<uint8_t> { 0, 1, 2 }), levels);
    EXPECT_LT(300u, anchors);
}

TEST(PatternMismatchPeptide, PalindromicKmerAtALevel) {
    // on an even-k wrapped PRIMARY graph AAATTT is its own reverse complement and holds KK with
    // the second K substituted by F: a substituted residue changes the bases of three mirrored
    // pairs, so the palindrome must be found at level 1 and counted once
    const std::vector<std::string> records { "AAATTT", "CCAAATTTGG", "AAGAAA" };
    for (size_t k : { 6, 8 }) {
        auto graph = build(k, records, DeBruijnGraph::PRIMARY);
        for (const std::string peptide : { "KK", "KF", "FK", "NF" }) {
            for (uint8_t m = 1; m <= 2; ++m) {
                SCOPED_TRACE(peptide + " k " + std::to_string(k) + " m " + std::to_string(m));
                check_peptide(*graph, peptide, 1,
                              make_request(Scope::ANY_OFFSET, Strands::BOTH, m), &records,
                              DeBruijnGraph::PRIMARY);
            }
        }
        if (k == 6) {
            const std::vector<MCtx> kk = peptide_walk(*graph, "KK", 1, Scope::ANY_OFFSET,
                                                      Strands::BOTH, 1);
            EXPECT_TRUE(std::any_of(kk.begin(), kk.end(), [&](const MCtx &c) {
                return graph->get_node_sequence(c.node) == "AAATTT" && c.level == 1;
            }));
        }
    }
}

TEST(PatternMismatchPeptide, UnmaskedTwins) {
    // without the mask every level is a true relation of the masked count, and the release
    // enumerates every candidate, the source dummies dropped: the oracle's list, every level
    // EXACT
    size_t bounded = 0;
    for (uint32_t seed = 1300; seed < 1360; ++seed) {
        const PeptideCase c = random_peptide_case(seed);
        auto unmasked = build(c.k, c.records, c.mode);
        const_cast<DBGSuccinct&>(stored_graph(*unmasked)).reset_mask();
        SCOPED_TRACE(describe(c));
        const Pattern pattern = peptide_of(c.peptide, c.table);
        const std::vector<MCtx> expected = peptide_walk(*c.graph, c.peptide, c.table,
                                                        c.request.scope, c.request.strands,
                                                        c.request.max_mismatches);
        const Truth truth = truth_of(expected, c.request.max_mismatches,
                                     scope_offsets(pattern.length(), c.k, c.request.scope),
                                     peptide_orientations(c.peptide, oracle_code(c.table),
                                                          c.request.strands));
        const Result counted = count_of(*unmasked, pattern, c.request);
        ASSERT_FALSE(counted.stop);
        for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
            const Count &level = counted.contexts->by_mismatches.at(t);
            EXPECT_TRUE(true_relation(level, truth.by_level[t])) << int(t);
            bounded += level.relation == Relation::BOUNDS;
        }
        Result all;
        EXPECT_EQ(expected, released_by(*unmasked, pattern, c.request, &all));
        for (uint8_t t = 0; t <= c.request.max_mismatches; ++t) {
            EXPECT_EQ(Relation::EXACT, all.contexts->by_mismatches.at(t).relation);
            EXPECT_EQ(truth.by_level[t], all.contexts->by_mismatches.at(t).value);
        }
    }
    EXPECT_LT(20u, bounded);
}

TEST(PatternMismatchPeptide, ResidueSubstitutionsInTheSubmittedPeptidesCoordinates) {
    const Pattern p = peptide_of("MKB", 1);
    EXPECT_TRUE(residue_substitutions(p, Orientation::FORWARD, "ATGAAAGAT").empty());
    auto subs = residue_substitutions(p, Orientation::FORWARD, "ATAAAAGAA");
    ASSERT_EQ(2u, subs.size());
    EXPECT_EQ(0u, subs[0].residue);
    EXPECT_EQ('M', subs[0].expected);
    EXPECT_EQ('I', subs[0].found);
    EXPECT_EQ("ATA", subs[0].codon);
    EXPECT_EQ(2u, subs[1].residue);
    EXPECT_EQ('B', subs[1].expected);
    EXPECT_EQ('E', subs[1].found);
    EXPECT_EQ("GAA", subs[1].codon);
    // reverse: the instance instantiates rc(P); residue 1 read as the reverse complement of the
    // instance's middle codon
    subs = residue_substitutions(p, Orientation::REVERSE,
                                 reverse_complement_string("ATGCGTAAT"));
    ASSERT_EQ(1u, subs.size());
    EXPECT_EQ(1u, subs[0].residue);
    EXPECT_EQ('R', subs[0].found);
    EXPECT_EQ("CGT", subs[0].codon);
    // ATA is Met in table 2
    EXPECT_TRUE(residue_substitutions(peptide_of("MKB", 2), Orientation::FORWARD,
                                      "ATAAAAGAT").empty());
    // no instance at any level: a stop at K, a sense codon at '*'; and the arguments
    EXPECT_THROW(residue_substitutions(p, Orientation::FORWARD, "ATGTAAGAT"),
                 std::invalid_argument);
    EXPECT_THROW(residue_substitutions(peptide_of("M*", 1), Orientation::FORWARD, "ATGTGG"),
                 std::invalid_argument);
    EXPECT_THROW(residue_substitutions(p, Orientation::FORWARD, "ATGAAAGA"),
                 std::invalid_argument);
    EXPECT_THROW(residue_substitutions(p, Orientation::FORWARD, "ATGAANGAT"),
                 std::invalid_argument);
    EXPECT_THROW(residue_substitutions(parse("ACG"), Orientation::FORWARD, "ACG"),
                 std::invalid_argument);
}

TEST(PatternMismatch, LongPatternAnchorsPerLevel) {
    // a DNA or IUPAC pattern longer than k: its anchors per level, each the anchor window of its
    // orientation (FORWARD P[0, k), REVERSE rc(P)[0, k)) within m, against O1 on the windows;
    // the paths UNKNOWN (not extended)
    size_t anchors = 0;
    for (uint32_t seed = 1400; seed < 1480; ++seed) {
        std::mt19937 rng(seed);
        const size_t k = std::vector<size_t> { 5, 6, 7, 9 }[rng() % 4];
        const auto mode = std::vector<DeBruijnGraph::Mode> {
            DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }[rng() % 3];
        const std::vector<std::string> records = random_records(rng, 4, 20, 60);
        auto graph = build(k, records, mode);
        const std::string text = random_pattern(rng, k + 1 + rng() % k, rng() % 2);
        if (palindromic(read_pattern(text)))
            continue;
        const uint8_t m = 1 + rng() % 3;
        const Strands strands = std::vector<Strands> {
            Strands::BOTH, Strands::FORWARD, Strands::REVERSE }[rng() % 3];
        SCOPED_TRACE(text + " k " + std::to_string(k) + " " + mode_name(mode) + " m "
                     + std::to_string(m));
        const size_t L = text.size();
        std::vector<MCtx> expected;
        // (a window equal to its reverse complement is O1's PALINDROMIC: relabelled |as|)
        auto add = [&](const std::string &window, Strands window_strands, Orientation as) {
            for (MCtx c : hamming_walk(*graph, window, Scope::ANY_OFFSET, window_strands, m)) {
                c.orientation = as;
                c.positions.clear();
                expected.push_back(c);
            }
        };
        if (strands != Strands::REVERSE)
            add(text.substr(0, k), Strands::FORWARD, Orientation::FORWARD);
        if (strands != Strands::FORWARD)
            add(text.substr(L - k), Strands::REVERSE, Orientation::REVERSE);
        std::sort(expected.begin(), expected.end());
        Request request = make_request(Scope::ANY_OFFSET, strands, m);
        const Result counted = count_of(*graph, text, request);
        ASSERT_TRUE(counted.anchors);
        ASSERT_FALSE(counted.stop);
        const Truth truth = truth_of(expected, m, { 0 }, counted.searched);
        EXPECT_EQ(truth.total, counted.anchors->total.value);
        for (uint8_t t = 0; t <= m; ++t) {
            EXPECT_EQ(truth.by_level[t], counted.anchors->by_mismatches.at(t).value) << int(t);
        }
        request.release_anchors = true;
        Result released;
        EXPECT_EQ(expected, released_by(*graph, text, request, &released));
        anchors += expected.size();
    }
    EXPECT_LT(200u, anchors);
}

#endif // _DNA_GRAPH

} // namespace
