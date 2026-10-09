#include "gtest/gtest.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <type_traits>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/test_trie_oracle.hpp"
#include "tests/graph/traversal/test_trie_checks.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/seq_tools/reverse_complement.hpp"


/*
 * Edge and non-edge cases of the trie contract (spec §6.9), each checked THREE ways:
 *
 *   records  --(test_trie_reference.hpp)-->  what the walk rule admits over the strings
 *   annotate --(structural trie T)--------->  what the walker records with no label logic
 *   constrain --(walker A over P)---------->  what the label-constrained walker claims
 *
 * check_case_on() demands: leaves(T) == structural walks of the records with the same
 * per-node label sets and path reasons; claims(A over P) == the oracle E read off T
 * (test_trie_oracle.hpp) == the per-label claims computed from the records, with the
 * same end reasons; the multi-label run equals the union of the single-label runs (the
 * "complete list of all single-label walks"); a permitted set derived from the seed
 * equals the explicit list; and the §6.3 recurrence at budget 0 reproduces the leaves.
 * Every case adds literal expectations derived from its construction on top, so that
 * the three parties cannot agree on something the records do not say.
 *
 * The cases: linear walks and seeds at record ends, the radius at every boundary, a
 * fork, an unequal bubble, nested bubbles, a one-base tip at the seed boundary, a
 * homopolymer self-loop, a cycle junction, a closed circle, a repeat shared by two
 * records plus a chimera, a seed occurring twice in one record, a seed spanning a
 * bubble, an unlabeled region, a reverse-complement record (basic vs canonical), a
 * hairpin (skip and follow), the even-k palindromic node rule, 100 labels against the
 * per-node cap, one and two label switches against the recurrence, trace support
 * against the records (contrast with k-mer support, record ends, a repeated seed, the
 * homopolymer limitation), a sweep of size caps against the completeness guarantee,
 * dummy nodes, and random graphs at k = 7.
 */
namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;
using namespace mtg::test::trie;

const size_t kDefaultK = 11;

std::string rc(std::string s) { ::reverse_complement(s); return s; }
std::string reversed(std::string s) { std::reverse(s.begin(), s.end()); return s; }

std::string random_seq(size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(len, 'A');
    for (char &c : s) c = "ACGT"[gen() % 4];
    return s;
}

// no repeated k-mer or (k-1)-mer in either orientation, no RC palindromes of length
// k-1, k or k+1 (only even lengths can be)
bool is_clean(const std::string &s, size_t k) {
    for (size_t len : { k - 1, k }) {
        std::set<std::string> seen;
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string km = s.substr(i, len);
            if (seen.count(km) || seen.count(rc(km)))
                return false;
            seen.insert(km);
        }
    }
    for (size_t len : { k - 1, k, k + 1 }) {
        if (len % 2)
            continue;
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string x = s.substr(i, len);
            if (rc(x) == x)
                return false;
        }
    }
    return true;
}

std::vector<std::string> clean_blocks(const std::vector<size_t> &lengths, uint32_t seed,
                                      size_t k = kDefaultK) {
    size_t total = 0;
    for (size_t l : lengths) total += l;
    std::string master;
    for (uint32_t s = seed; ; ++s) {
        master = random_seq(total, s);
        if (is_clean(master, k))
            break;
    }
    std::vector<std::string> blocks;
    size_t pos = 0;
    for (size_t l : lengths) {
        blocks.push_back(master.substr(pos, l));
        pos += l;
    }
    return blocks;
}

// clean blocks satisfying an extra structural condition (distinct fork bases etc.)
template <class Cond>
std::vector<std::string> blocks_where(const std::vector<size_t> &lengths, uint32_t seed,
                                      const Cond &cond, size_t k = kDefaultK) {
    for (uint32_t s = seed; ; ++s) {
        auto b = clean_blocks(lengths, s, k);
        if (cond(b))
            return b;
    }
}

std::vector<DeBruijnGraph::Mode> all_modes() {
    return {
        DeBruijnGraph::BASIC,
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

std::vector<DeBruijnGraph::Mode> canonical_modes() {
    return {
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

template <class Graph, class Annotation>
std::map<DeBruijnGraph::Mode, CaseResult>
check_case(const CaseSpec &c, std::vector<DeBruijnGraph::Mode> modes = all_modes()) {
    std::map<DeBruijnGraph::Mode, CaseResult> out;
    for (auto mode : modes) {
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        trie::StringIndex ix(c.k, mode, c.sequences, c.labels);
        out.emplace(mode, check_case_on(*anno, ix, mode, c));
    }
    return out;
}

} // namespace


template <typename Pair>
class TrieCases : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::RowFlatAnnotator>
> TrieCaseTypes;
TYPED_TEST_SUITE(TrieCases, TrieCaseTypes);

#define CASE_TYPES \
    using Graph = typename TypeParam::first_type; \
    using Annotation = typename TypeParam::second_type


/*********************************** linear ***********************************/

// The non-edge case: one record, the seed inside it, both arms spell the record.
TYPED_TEST(TrieCases, Linear) {
    CASE_TYPES;
    auto b = clean_blocks({ 60, 30, 80 }, 101);
    const std::string &L = b[0], &S = b[1], &R = b[2];
    CaseSpec c { "Linear", kDefaultK, { L + S + R }, { "A" }, S };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A" }, kRight, { { R, { "A" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kLeft, { { reversed(L), { "A" } } }, c.radius, where);
        expect_path_reason(r, kRight, R, EndReason::DEAD_END, where);
        expect_path_reason(r, kLeft, reversed(L), EndReason::DEAD_END, where);
        EXPECT_EQ(1u, r.T.arms[kRight].paths.size());
        EXPECT_TRUE(r.T.arms[kRight].splits.empty());
    }
}

// Seeds at a record's start, its end, and a seed that IS the record: an arm with
// nothing to add is one empty walk ended dead_end, not a missing arm.
TYPED_TEST(TrieCases, SeedAtRecordEnds) {
    CASE_TYPES;
    auto b = clean_blocks({ 50, 30, 50 }, 102);
    const std::string &L = b[0], &S = b[1], &R = b[2];
    struct Sub { std::string name, record; trie::Leaves left, right; };
    std::vector<Sub> subs {
        { "SeedAtRecordStart", S + R, { { "", { "A" } } }, { { R, { "A" } } } },
        { "SeedAtRecordEnd", L + S, { { reversed(L), { "A" } } }, { { "", { "A" } } } },
        { "SeedIsRecord", S, { { "", { "A" } } }, { { "", { "A" } } } },
    };
    for (const Sub &s : subs) {
        CaseSpec c { s.name, kDefaultK, { s.record }, { "A" }, S };
        for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
            const std::string where = c.name + " " + mode_name(mode);
            expect_leaves(r, { "A" }, kLeft, s.left, c.radius, where);
            expect_leaves(r, { "A" }, kRight, s.right, c.radius, where);
            for (size_t a : { kLeft, kRight }) {
                EXPECT_EQ(1u, r.T.arms[a].paths.size()) << where;
                EXPECT_EQ(ArmResult::COMPLETE, r.T.arms[a].status) << where;
            }
        }
    }
}

// The radius at every boundary of a linear walk: 0, 1, 2, one short of the record's
// end, exactly it, one past it. A walk cut by the radius ends max_extension_bp and the
// arm is complete; past the record's end the radius changes nothing.
TYPED_TEST(TrieCases, RadiusAtEveryBoundary) {
    CASE_TYPES;
    auto b = clean_blocks({ 40, 30, 50 }, 103);
    const std::string &L = b[0], &S = b[1], &R = b[2];
    for (uint64_t radius : { 0ul, 1ul, 2ul, R.size() - 1, R.size(), R.size() + 1 }) {
        CaseSpec c { "Radius" + std::to_string(radius), kDefaultK, { L + S + R }, { "A" }, S, radius };
        for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
            const std::string where = c.name + " " + mode_name(mode);
            const std::string right = R.substr(0, std::min<size_t>(radius, R.size()));
            const std::string left = reversed(L).substr(0, std::min<size_t>(radius, L.size()));
            expect_leaves(r, { "A" }, kRight, { { right, { "A" } } }, radius, where);
            expect_leaves(r, { "A" }, kLeft, { { left, { "A" } } }, radius, where);
            // a head AT the radius is ended max_extension_bp before its successors are
            // looked at, so a radius equal to the flank's length reports the radius
            expect_path_reason(r, kRight, right,
                               radius <= R.size() ? EndReason::MAX_EXTENSION : EndReason::DEAD_END, where);
            expect_path_reason(r, kLeft, left,
                               radius <= L.size() ? EndReason::MAX_EXTENSION : EndReason::DEAD_END, where);
            for (size_t a : { kLeft, kRight }) {
                EXPECT_EQ(ArmResult::COMPLETE, r.T.arms[a].status) << where;
                EXPECT_EQ(radius, r.T.arms[a].complete_to_bp) << where;
            }
        }
    }
}


/********************************* branching *********************************/

// Two labels fork right after the seed: a divergence, two leaves, each under one label.
TYPED_TEST(TrieCases, Fork) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 40, 40 }, 104, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    CaseSpec c { "Fork", kDefaultK, { X + P, X + Q }, { "A", "B" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B" }, kRight, { { P, { "A" } }, { Q, { "B" } } }, c.radius, where);
        expect_leaves(r, { "A", "B" }, kLeft, { { "", { "A", "B" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kRight, { { P, { "A" } } }, c.radius, where);
        const ArmResult &t = r.T.arms[kRight];
        ASSERT_EQ(1u, t.splits.size()) << where;
        EXPECT_EQ(0u, t.splits[0].at_bp) << where;
        EXPECT_FALSE(t.splits[0].ambiguous) << where;
        // the trie's view of the fork: the boundary carried by both labels, one per branch,
        // nothing cut; the per-label summary read off the recorded sets
        EXPECT_EQ(0u, t.nodes_labels_truncated) << where;
        EXPECT_EQ(2u, t.max_labels_at_node) << where;
        EXPECT_EQ(2u, t.splits[0].labels_before) << where;
        ASSERT_EQ(2u, t.splits[0].branches.size()) << where;
        for (const auto &br : t.splits[0].branches) {
            EXPECT_EQ(1u, br.labels_distinct) << where;
            EXPECT_EQ(br.labels, t.segments[br.segment].labels_start) << where;
        }
        for (LabelId l = 0; l < 2; ++l) {
            EXPECT_EQ(P.size(), r.T.label_summary[l][kRight].direct_bp) << where;
            EXPECT_EQ(P.size(), r.T.label_summary[l][kRight].reach_bp) << where;
        }
        // a divergence in constrain mode: one lineage per child, not ambiguous
        const ArmResult &ab = r.A.at({ "A", "B" }).arms[kRight];
        ASSERT_EQ(1u, ab.splits.size()) << where;
        EXPECT_FALSE(ab.splits[0].ambiguous) << where;
    }
}

// A bubble whose branches differ in length rejoins at two different depths: both
// walks are complete in the trie, and with merging on the routes stay sound.
TYPED_TEST(TrieCases, UnequalBubble) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 15, 22, 30, 30 }, 105, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Y = b[3], &Z = b[4];
    CaseSpec c { "UnequalBubble", kDefaultK, { X + P + Y + Z, X + Q + Y + Z }, { "A", "B" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B" }, kRight,
                      { { P + Y + Z, { "A" } }, { Q + Y + Z, { "B" } } }, c.radius, where);
        // merging on, otherwise unlimited: whatever joins, every label's route is a
        // prefix of an exhaustive leaf under that label
        Strategy st;
        st.max_extension_bp = c.radius;
        st.max_label_branches = Strategy::kUnlimited;
        st.max_splits_per_path = Strategy::kUnlimited;
        st.merge_reconverge = true;
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        SeedResult merged = run(*anno, X, { "A", "B" }, st);
        EXPECT_EQ(2u, merged.arms[kRight].paths.size()) << where;
        EXPECT_GE(trie::check_routes_subset(r.A.at({ "A", "B" }), merged, kRight, where), 2u);
    }
}

// A bubble inside a bubble branch: three leaves over three labels, splits at 0 and |P|.
TYPED_TEST(TrieCases, NestedBubbles) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 20, 20, 20, 25, 30 }, 106, [](const auto &b) {
        return b[1][0] != b[4][0] && b[2][0] != b[3][0];
    });
    const std::string &X = b[0], &P = b[1], &Y1 = b[2], &Y2 = b[3], &Q = b[4], &Z = b[5];
    CaseSpec c { "NestedBubbles", kDefaultK,
                 { X + P + Y1 + Z, X + P + Y2 + Z, X + Q + Z }, { "A", "B", "C" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B", "C" }, kRight,
                      { { P + Y1 + Z, { "A" } }, { P + Y2 + Z, { "B" } }, { Q + Z, { "C" } } },
                      c.radius, where);
        expect_leaves(r, { "A", "B" }, kRight,
                      { { P + Y1 + Z, { "A" } }, { P + Y2 + Z, { "B" } } }, c.radius, where);
        std::set<uint64_t> at;
        for (const Split &s : r.T.arms[kRight].splits) at.insert(s.at_bp);
        EXPECT_EQ((std::set<uint64_t>{ 0, P.size() }), at) << where;
    }
}

// A one-base tip hanging off the seed boundary node: the shortest possible walk, a
// leaf of its own under its own label, and a divergence (not an ambiguity).
TYPED_TEST(TrieCases, TipAtTheSeedBoundary) {
    CASE_TYPES;
    auto b = clean_blocks({ 40, 40 }, 107);
    const std::string &X = b[0], &R = b[1];
    char t = 'A';
    while (t == R[0]) t = "ACGT"[(std::string("ACGT").find(t) + 1) % 4];
    CaseSpec c { "TipAtTheSeedBoundary", kDefaultK, { X + R, X + t }, { "A", "B" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B" }, kRight,
                      { { R, { "A" } }, { std::string(1, t), { "B" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kRight, { { R, { "A" } } }, c.radius, where);
        expect_path_reason(r, kRight, std::string(1, t), EndReason::DEAD_END, where);
        const ArmResult &ab = r.A.at({ "A", "B" }).arms[kRight];
        ASSERT_EQ(1u, ab.splits.size()) << where;
        EXPECT_FALSE(ab.splits[0].ambiguous) << where;
    }
}


/****************************** cycles and repeats ******************************/

// A homopolymer longer than k + 1 is a self-loop on the node A^k, and every node whose
// (k-1)-suffix is A^(k-1) has the exit A^(k-1)·Z[0] as a successor. The loop edge may
// be used once per walk, so the trie holds A^(k-1)·Z, A^k·Z and A^(k+1)·Z and NOT the
// record's own A^(k+3)·Z. A known consequence of the walk rule, pinned here.
TYPED_TEST(TrieCases, HomopolymerSelfLoop) {
    CASE_TYPES;
    const size_t k = kDefaultK;
    auto b = blocks_where({ 30, 30 }, 108, [](const auto &b) {
        return b[0].back() != 'A' && b[1][0] != 'A' && b[0].back() != 'T' && b[1][0] != 'T';
    });
    const std::string &X = b[0], &Z = b[1];
    const std::string H(k + 3, 'A');
    CaseSpec c { "HomopolymerSelfLoop", k, { X + H + Z }, { "A" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A" }, kRight,
                      { { std::string(k - 1, 'A') + Z, { "A" } }, { std::string(k, 'A') + Z, { "A" } },
                        { std::string(k + 1, 'A') + Z, { "A" } } },
                      c.radius, where);
        EXPECT_EQ(0u, leaf_walks(r.T.arms[kRight]).count(H + Z)) << where << ": the record's walk";
        EXPECT_EQ(1u, count_blocked(r.T.arms[kRight], EndReason::EDGE_REUSE)) << where;
    }
}

// L·J·C·J·E: the junction J is entered twice; the second lap through C is blocked by
// edge reuse, so the trie holds J·E and J·C·J·E.
TYPED_TEST(TrieCases, CycleJunction) {
    CASE_TYPES;
    auto b = blocks_where({ 30, kDefaultK, 40, 30 }, 109, [](const auto &b) { return b[2][0] != b[3][0]; });
    const std::string &L = b[0], &J = b[1], &C = b[2], &E = b[3];
    CaseSpec c { "CycleJunction", kDefaultK, { L + J + C + J + E }, { "A" }, L };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A" }, kRight, { { J + E, { "A" } }, { J + C + J + E, { "A" } } }, c.radius, where);
        EXPECT_EQ(1u, count_blocked(r.T.arms[kRight], EndReason::EDGE_REUSE)) << where;
        // ambiguous in constrain mode: one label on two children
        const ArmResult &a = r.A.at({ "A" }).arms[kRight];
        ASSERT_EQ(1u, a.splits.size()) << where;
        EXPECT_TRUE(a.splits[0].ambiguous) << where;
        EXPECT_EQ(J.size(), a.splits[0].at_bp) << where;
    }
}

// A closed circle: both arms walk round until they would re-enter the seed, and the
// two arms cover the same edges independently.
TYPED_TEST(TrieCases, Circle) {
    CASE_TYPES;
    std::string S = clean_blocks({ 150 }, 110)[0];
    std::string record = S + S.substr(0, kDefaultK - 1 + 10);
    std::string seed = S.substr(50, 40);
    CaseSpec c { "Circle", kDefaultK, { record }, { "A" }, seed };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        const std::string right = S.substr(90) + S.substr(0, 60);
        const std::string left = reversed(S.substr(80) + S.substr(0, 50));
        expect_leaves(r, { "A" }, kRight, { { right, { "A" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kLeft, { { left, { "A" } } }, c.radius, where);
        expect_path_reason(r, kRight, right, EndReason::REACHED_SEED, where);
        expect_path_reason(r, kLeft, left, EndReason::REACHED_SEED, where);
        EXPECT_EQ(S.size() - seed.size() + kDefaultK - 1, right.size());
    }
}

// A repeat R in two records plus a chimeric third: the seed's lineages diverge on
// both sides, and the chimera C shares A's left with B's right.
TYPED_TEST(TrieCases, RepeatInTwoRecordsAndAChimera) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 30, 30, 30, 30 }, 111, [](const auto &b) {
        return b[2][0] != b[4][0] && b[0].back() != b[3].back();
    });
    const std::string &X = b[0], &R = b[1], &Y = b[2], &Z = b[3], &W = b[4];
    CaseSpec c { "RepeatInTwoRecordsAndAChimera", kDefaultK,
                 { X + R + Y, Z + R + W, X + R + W }, { "A", "B", "C" }, R };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B", "C" }, kRight, { { Y, { "A" } }, { W, { "B", "C" } } }, c.radius, where);
        expect_leaves(r, { "A", "B", "C" }, kLeft,
                      { { reversed(X), { "A", "C" } }, { reversed(Z), { "B" } } }, c.radius, where);
        expect_leaves(r, { "B" }, kLeft, { { reversed(Z), { "B" } } }, c.radius, where);
    }
}

// The seed occurs twice in one record (X·P·X·Q): the right arm follows both
// continuations under one label; the P walk goes on through the k-1 junction k-mers
// and ends where it would enter the seed's first node again; the left arm comes back
// to the seed through P the same way.
TYPED_TEST(TrieCases, RepeatedSeedInOneRecord) {
    CASE_TYPES;
    const size_t k = kDefaultK;
    auto b = blocks_where({ 30, 40, 40 }, 112, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    CaseSpec c { "RepeatedSeedInOneRecord", k, { X + P + X + Q }, { "A" }, X };
    const std::string right_p = P + X.substr(0, k - 1);
    const std::string left_p = reversed(X.substr(X.size() - (k - 1)) + P);
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A" }, kRight, { { right_p, { "A" } }, { Q, { "A" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kLeft, { { left_p, { "A" } } }, c.radius, where);
        expect_path_reason(r, kRight, right_p, EndReason::REACHED_SEED, where);
        expect_path_reason(r, kRight, Q, EndReason::DEAD_END, where);
        expect_path_reason(r, kLeft, left_p, EndReason::REACHED_SEED, where);
        EXPECT_EQ(1u, count_label_ends(r.A.at({ "A" }).arms[kRight], EndReason::REACHED_SEED)) << where;
    }
}

// The seed itself spans a bubble: only the record spelling the whole seed carries it.
// Naming the other label drops it (seed not supported) and changes nothing else.
TYPED_TEST(TrieCases, SeedSpansABubble) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 20, 20, 30, 30, 30 }, 113, [](const auto &b) {
        return b[1][0] != b[2][0] && b[4][0] != b[5][0];
    });
    const std::string &X1 = b[0], &P = b[1], &Q = b[2], &X2 = b[3], &Y = b[4], &Z = b[5];
    CaseSpec c { "SeedSpansABubble", kDefaultK, { X1 + P + X2 + Y, X1 + Q + X2 + Z }, { "A", "B" }, X1 + P + X2 };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        EXPECT_EQ((std::set<std::string>{ "A" }), r.carriers) << where;
        expect_leaves(r, { "A" }, kRight, { { Y, { "A" } } }, c.radius, where);
        expect_leaves(r, { "A" }, kLeft, { { "", { "A" } } }, c.radius, where);
        EXPECT_EQ((std::set<std::string>{ Y, Z }), leaf_walks(r.T.arms[kRight])) << where;
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        SeedResult both = run(*anno, c.seed, { "A", "B" }, exhaustive(LabelMode::CONSTRAIN, c.radius));
        ASSERT_EQ(1u, both.dropped_labels.size()) << where;
        EXPECT_EQ("B", both.dropped_labels[0].name) << where;
        for (size_t a : { kLeft, kRight }) {
            trie::expect_equal(trie::constrained_claims(r.A.at({ "A" }), a, c.radius),
                               trie::constrained_claims(both, a, c.radius), where + " [B dropped]");
        }
    }
}


/****************************** labels and strands ******************************/

namespace {

// Part of the graph is carried by no label: the structural trie walks through it
// recording empty sets, the constrained walker stops at its edge with label_lost.
// (The graph holds both records, the annotation only the labelled one.)
template <class Graph>
void unlabeled_region() {
    auto b = clean_blocks({ 30, 30, 30 }, 114);
    const std::string &X = b[0], &P = b[1], &Z = b[2];
    const std::vector<std::string> seqs { X + P, X + P + Z }, labels { "A", "" };
    CaseSpec c { "UnlabeledRegion", kDefaultK, seqs, labels, X };
    for (auto mode : all_modes()) {
        std::shared_ptr<DeBruijnGraph> graph = build_graph_batch<Graph>(c.k, seqs, mode);
        auto canonical = std::dynamic_pointer_cast<const CanonicalDBG>(graph);
        std::shared_ptr<const DeBruijnGraph> base = canonical ? canonical->get_graph_ptr() : graph;
        AnnotatedDBG anno(graph, std::make_unique<annot::ColumnCompressed<>>(base->max_index()));
        anno.annotate_sequence(std::string(seqs[0]), { "A" });
        trie::StringIndex ix(c.k, mode, seqs, labels);
        CaseResult r = check_case_on(anno, ix, mode, c);
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A" }, kRight, { { P, { "A" } } }, c.radius, where);
        EXPECT_EQ((std::set<std::string>{ P + Z }), leaf_walks(r.T.arms[kRight])) << where;
        EXPECT_EQ(1u, count_label_ends(r.A.at({ "A" }).arms[kRight], EndReason::LABEL_LOST)) << where;
        ASSERT_EQ(1u, r.T.arms[kRight].paths.size()) << where;
        bool cut = false;
        const auto at = trie::recorded_labels(r.T, r.T.arms[kRight], r.T.arms[kRight].paths[0], &cut);
        EXPECT_TRUE(at.back().empty()) << where << ": the last node of P·Z is unlabeled";
        EXPECT_EQ((std::set<std::string>{ "A" }), at[P.size()]) << where;
    }
}

} // namespace

TEST(TrieCasesLabels, UnlabeledRegionSuccinct) { unlabeled_region<DBGSuccinct>(); }
TEST(TrieCasesLabels, UnlabeledRegionHashFast) { unlabeled_region<DBGHashFast>(); }

// A record stored as the reverse complement of the continuation: in a basic graph it
// is not even in the same strand as the seed; in canonical and primary graphs it is
// the second branch of the fork.
TYPED_TEST(TrieCases, ReverseComplementRecord) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 40, 40 }, 115, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    CaseSpec c { "ReverseComplementRecord", kDefaultK, { X + P, rc(X + Q) }, { "A", "B" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        if (mode == DeBruijnGraph::BASIC) {
            EXPECT_EQ((std::set<std::string>{ "A" }), r.carriers) << where;
            EXPECT_EQ((std::set<std::string>{ P }), leaf_walks(r.T.arms[kRight])) << where;
            expect_leaves(r, { "A" }, kRight, { { P, { "A" } } }, c.radius, where);
        } else {
            EXPECT_EQ((std::set<std::string>{ "A", "B" }), r.carriers) << where;
            expect_leaves(r, { "A", "B" }, kRight, { { P, { "A" } }, { Q, { "B" } } }, c.radius, where);
            expect_leaves(r, { "A", "B" }, kLeft, { { "", { "A", "B" } } }, c.radius, where);
        }
    }
}

// A self-reverse-complementary (k-1)-mer makes one step a hairpin: skipped, the record
// is reconstructed with one hairpin event; followed, the retrace is blocked at once.
TYPED_TEST(TrieCases, Hairpin) {
    CASE_TYPES;
    const size_t k = kDefaultK;
    const std::string X = "AACGTACGTT";
    ASSERT_EQ(X, rc(X));
    ASSERT_EQ(k - 1, X.size());
    std::string L, R;
    for (uint32_t seed = 116; ; ++seed) {
        L = clean_blocks({ 40 }, seed)[0];
        R = clean_blocks({ 40 }, seed + 100)[0];
        if (R[0] != rc(std::string(1, L.back()))[0] && is_clean(L + X.substr(0, 3), k)
                && is_clean(X.substr(7) + R, k))
            break;
    }
    const std::string S = L + X + R, seed = L.substr(0, 30);
    for (bool skip : { true, false }) {
        CaseSpec c { skip ? "HairpinSkip" : "HairpinFollow", k, { S }, { "A" }, seed, 200, {}, skip };
        for (auto &[mode, r] : check_case<Graph, Annotation>(c, canonical_modes())) {
            const std::string where = c.name + " " + mode_name(mode);
            const ArmResult &right = r.T.arms[kRight];
            EXPECT_EQ(1u, count_events(right, EventType::HAIRPIN)) << where;
            const std::string real = S.substr(30);
            EXPECT_EQ(1u, leaf_walks(right).count(real)) << where;
            expect_path_reason(r, kRight, real, EndReason::DEAD_END, where);
            if (skip) {
                expect_leaves(r, { "A" }, kRight, { { real, { "A" } } }, c.radius, where);
                EXPECT_TRUE(right.splits.empty()) << where;
            } else {
                // the hairpin step plus one base of retrace, blocked as edge_reuse_rc
                // (or as seed re-entry when the retraced edge is the seed's)
                ASSERT_EQ(2u, right.paths.size()) << where;
                const uint64_t hp = L.size() - 30 + k - 1;
                for (const auto &[w, reason] : path_reasons(right)) {
                    if (w == real)
                        continue;
                    EXPECT_EQ(hp + 1, w.size()) << where;
                    EXPECT_TRUE(reason == EndReason::EDGE_REUSE_RC || reason == EndReason::REACHED_SEED)
                        << where << " " << to_string(reason);
                }
            }
        }
    }
}

// For even k a node equal to its own reverse complement is a hairpin in itself: the
// step into it is skipped, so the walk stops one base short of the palindrome.
TYPED_TEST(TrieCases, EvenKPalindromicNode) {
    CASE_TYPES;
    const size_t k = 10;
    const std::string X = "ACGTTAACGT";
    ASSERT_EQ(X, rc(X));
    ASSERT_EQ(k, X.size());
    std::string L, R;
    for (uint32_t seed = 117; ; ++seed) {
        L = clean_blocks({ 40 }, seed, k)[0];
        R = clean_blocks({ 40 }, seed + 100, k)[0];
        if (is_clean(L + X.substr(0, k - 1), k) && is_clean(X.substr(1) + R, k))
            break;
    }
    const std::string S = L + X + R, seed = L.substr(0, 30);
    CaseSpec c { "EvenKPalindromicNode", k, { S }, { "A" }, seed };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c, canonical_modes())) {
        const std::string where = c.name + " " + mode_name(mode);
        const std::string walk = L.substr(30) + X.substr(0, k - 1);
        expect_leaves(r, { "A" }, kRight, { { walk, { "A" } } }, c.radius, where);
        expect_path_reason(r, kRight, walk, EndReason::DEAD_END, where);
        EXPECT_EQ(1u, count_events(r.T.arms[kRight], EventType::HAIRPIN)) << where;
    }
}

// One hundred labels share the seed, fan out four ways and then each goes its own
// way. The default per-node cap (64) cuts the recorded lists and the oracle says so;
// uncapped, the contract holds for a hundred labels and a hundred leaves.
TYPED_TEST(TrieCases, HundredLabels) {
    CASE_TYPES;
    // one graph and annotation type: what costs here is the labels, which every type keeps
    // alike (the other cases run on both)
    if (!std::is_same_v<Graph, DBGSuccinct>)
        GTEST_SKIP() << "run on DBGSuccinct with ColumnCompressed only";
    const size_t n = 100;
    std::vector<size_t> lengths { 30, 20, 20, 20, 20 };
    for (size_t i = 0; i < n; ++i) lengths.push_back(12);
    auto b = blocks_where(lengths, 118, [](const auto &b) {
        std::set<char> first { b[1][0], b[2][0], b[3][0], b[4][0] };
        return first.size() == 4;
    });
    const std::string &X = b[0];
    std::vector<std::string> seqs, labels;
    std::set<std::string> all;
    for (size_t i = 0; i < n; ++i) {
        std::string name = "L" + std::to_string(100 + i);
        seqs.push_back(X + b[1 + i % 4] + b[5 + i]);
        labels.push_back(name);
        all.insert(name);
    }
    CaseSpec c { "HundredLabels", kDefaultK, seqs, labels, X, 200, { all, { "L100" }, { "L137" }, { "L199" } } };
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(c.k, seqs, labels, mode);
        trie::StringIndex ix(c.k, mode, seqs, labels);
        const std::string where = c.name + " " + mode_name(mode);
        // the default cap: the boundary list and the shared stretch are cut, reported
        // as such, and the oracle refuses to prove anything with them
        SeedResult capped = run(*anno, X, {}, exhaustive(LabelMode::ANNOTATE, c.radius));
        EXPECT_GT(capped.arms[kRight].nodes_labels_truncated, 0u) << where;
        EXPECT_EQ(n, capped.arms[kRight].max_labels_at_node) << where;
        EXPECT_EQ(n, trie::root_of(capped.arms[kRight]).labels_start_total) << where;
        EXPECT_EQ(64u, trie::root_of(capped.arms[kRight]).labels_start.size()) << where;
        SeedResult A = run(*anno, X, as_list(all), exhaustive(LabelMode::CONSTRAIN, c.radius));
        EXPECT_TRUE(trie::verify(capped, A, kRight, all).cut) << where;
        // uncapped: the full contract. (The structural trie may hold a few walks more
        // than there are labels: junction k-mers of two records can coincide by chance
        // and open a branch that no label carries in full.)
        CaseResult r = check_case_on(*anno, ix, mode, c);
        EXPECT_GE(r.T.arms[kRight].paths.size(), n) << where;
        EXPECT_EQ(n, r.A.at(all).arms[kRight].paths.size()) << where;
        EXPECT_EQ(n, r.T.arms[kRight].max_labels_at_node) << where;
        expect_leaves(r, { "L137" }, kRight, { { b[1 + 37 % 4] + b[5 + 37], { "L137" } } }, c.radius, where);
        expect_leaves(r, all, kLeft, { { "", all } }, c.radius, where);
    }
}


/********************************* switching *********************************/


// A: X·P·Y, B: X·Q·Y·Z. Under forbid A's walk ends in Y; one switch continues it into
// Z under B — a walk no record spells, which is what one switch means.
TYPED_TEST(TrieCases, OneSwitch) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 30, 30, 30, 30 }, 119, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Y = b[3], &Z = b[4];
    CaseSpec c { "OneSwitch", kDefaultK, { X + P + Y, X + Q + Y + Z }, { "A", "B" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B" }, kRight, { { P + Y, { "A" } }, { Q + Y + Z, { "B" } } }, c.radius, where);
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        trie::StringIndex ix(c.k, mode, c.sequences, c.labels);
        std::map<double, SeedResult> runs;
        check_switching(*anno, ix, c, { "A", "B" }, { 0, 1, 2 }, &runs, where);
        // literally: budget 1 spells P·Y·Z under B at loss 1, switched at |P·Y|
        const trie::SwitchTrie one = switch_leaves(runs.at(1), kRight);
        ASSERT_EQ(1u, one.count(P + Y + Z)) << where;
        EXPECT_EQ((std::map<std::string, double>{ { "B", 1.0 } }), one.at(P + Y + Z).state) << where;
        EXPECT_EQ((std::vector<std::pair<uint64_t, std::string>>{ { P.size() + Y.size(), "B" } }),
                  one.at(P + Y + Z).switches) << where;
        EXPECT_EQ((std::map<std::string, double>{ { "B", 0.0 } }), one.at(Q + Y + Z).state) << where;
        EXPECT_EQ(1u, count_label_ends(runs.at(1).arms[kRight], EndReason::LABEL_LOST)) << where;
        // budget 0 with a finite cost: the switch is reported as the budget it needed
        EXPECT_EQ(1u, count_label_ends(runs.at(0).arms[kRight], EndReason::LOSS_BUDGET)) << where;
        EXPECT_EQ((std::vector<double>{ 1.0 }), runs.at(0).arms[kRight].needed_budgets) << where;
        EXPECT_EQ(0u, count_events(runs.at(0).arms[kRight], EventType::SWITCH)) << where;
    }
}

// A: X·P·Y, B: X·Q·Y·Z, C: X·Q'·Z'·W where Z' is the end of Z: A's walk needs two
// switches to reach W, B's one. Budget 1 cuts A's at |P·Y·Z| with the needed budget
// reported; budget 2 spells P·Y·Z·W under C at loss 2.
TYPED_TEST(TrieCases, TwoSwitches) {
    CASE_TYPES;
    const size_t k = kDefaultK;
    auto b = blocks_where({ 30, 30, 30, 30, 30, 30, 30 }, 120, [](const auto &b) {
        return std::set<char>{ b[1][0], b[2][0], b[3][0] }.size() == 3;
    });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Q2 = b[3], &Y = b[4], &Z = b[5], &W = b[6];
    const std::string Z2 = Z.substr(Z.size() - (k - 1));
    CaseSpec c { "TwoSwitches", k, { X + P + Y, X + Q + Y + Z, X + Q2 + Z2 + W }, { "A", "B", "C" }, X };
    for (auto &[mode, r] : check_case<Graph, Annotation>(c)) {
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B", "C" }, kRight,
                      { { P + Y, { "A" } }, { Q + Y + Z, { "B" } }, { Q2 + Z2 + W, { "C" } } }, c.radius, where);
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        trie::StringIndex ix(c.k, mode, c.sequences, c.labels);
        std::map<double, SeedResult> runs;
        check_switching(*anno, ix, c, { "A", "B", "C" }, { 0, 1, 2 }, &runs, where);
        const trie::SwitchTrie one = switch_leaves(runs.at(1), kRight);
        const trie::SwitchTrie two = switch_leaves(runs.at(2), kRight);
        EXPECT_EQ((std::set<std::string>{ P + Y + Z, Q + Y + Z + W, Q2 + Z2 + W }), keys(one)) << where;
        EXPECT_EQ((std::set<std::string>{ P + Y + Z + W, Q + Y + Z + W, Q2 + Z2 + W }), keys(two)) << where;
        EXPECT_EQ((std::map<std::string, double>{ { "C", 2.0 } }), two.at(P + Y + Z + W).state) << where;
        EXPECT_EQ((std::map<std::string, double>{ { "C", 1.0 } }), one.at(Q + Y + Z + W).state) << where;
        EXPECT_EQ(1u, count_label_ends(runs.at(1).arms[kRight], EndReason::LOSS_BUDGET)) << where;
        EXPECT_EQ((std::vector<double>{ 2.0 }), runs.at(1).arms[kRight].needed_budgets) << where;
    }
}


/******************************** trace support ********************************/


// One label, two records that share a middle stretch and fork on both sides of it
// (X·P·Y·Z[:20]·Z[20:] and X·Q·Y·Z[:20]·W): k-mer support follows the four
// combinations, trace support the two records.
TEST(TrieCasesTrace, KmerVersusTrace) {
    auto b = blocks_where({ 30, 30, 30, 30, 30, 30 }, 121, [](const auto &b) {
        return b[1][0] != b[2][0] && b[4][20] != b[5][0];
    });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Y = b[3], &Z = b[4], &W = b[5];
    const std::string shared = Z.substr(0, 20);
    const std::vector<std::string> seqs { X + P + Y + Z, X + Q + Y + shared + W };
    const std::vector<std::string> labels { "F", "F" };
    const uint64_t n1 = seqs[0].size() - kDefaultK + 1;
    const std::vector<uint64_t> starts { 0, n1 + 100 };
    CaseSpec c { "KmerVersusTrace", kDefaultK, seqs, labels, X };
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        c.k, seqs, labels, DeBruijnGraph::BASIC, true, starts);
    trie::StringIndex ix(c.k, DeBruijnGraph::BASIC, seqs, labels, starts);
    CaseResult r = check_case_on(*anno, ix, DeBruijnGraph::BASIC, c);
    const std::string where = c.name;
    expect_leaves(r, { "F" }, kRight,
                  { { P + Y + Z, { "F" } }, { P + Y + shared + W, { "F" } },
                    { Q + Y + Z, { "F" } }, { Q + Y + shared + W, { "F" } } }, c.radius, where);
    SeedResult trace;
    check_trace(*anno, ix, c, { "F" }, r.A.at({ "F" }), &trace, where);
    EXPECT_EQ((trie::Leaves{ { P + Y + Z, { "F" } }, { Q + Y + shared + W, { "F" } } }),
              trie::constrained_leaves(trace, kRight, c.radius)) << where;
    EXPECT_EQ(2u, count_label_ends(trace.arms[kRight], EndReason::DEAD_END)) << where;
}

// A record ends where another record of the same label goes on: the k-mer walk
// continues, the trace ends with record_end (the label is present, its coordinates
// do not continue).
TEST(TrieCasesTrace, RecordEnd) {
    auto b = clean_blocks({ 30, 30, 30, 30, 30 }, 122);
    const std::string &X = b[0], &P = b[1], &Y = b[2], &Z = b[3], &M = b[4];
    const std::vector<std::string> seqs { X + P + Y + Z, Y + Z + M };
    const std::vector<std::string> labels { "F", "F" };
    const uint64_t n1 = seqs[0].size() - kDefaultK + 1;
    const std::vector<uint64_t> starts { 0, n1 + 100 };
    CaseSpec c { "RecordEnd", kDefaultK, seqs, labels, X };
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        c.k, seqs, labels, DeBruijnGraph::BASIC, true, starts);
    trie::StringIndex ix(c.k, DeBruijnGraph::BASIC, seqs, labels, starts);
    CaseResult r = check_case_on(*anno, ix, DeBruijnGraph::BASIC, c);
    expect_leaves(r, { "F" }, kRight, { { P + Y + Z + M, { "F" } } }, c.radius, c.name);
    SeedResult trace;
    check_trace(*anno, ix, c, { "F" }, r.A.at({ "F" }), &trace, c.name);
    EXPECT_EQ((trie::Leaves{ { P + Y + Z, { "F" } } }), trie::constrained_leaves(trace, kRight, c.radius));
    EXPECT_EQ(1u, count_label_ends(trace.arms[kRight], EndReason::RECORD_END));
    EXPECT_EQ(0u, count_label_ends(trace.arms[kRight], EndReason::DEAD_END));
}

// The seed occurs twice in one record: two traces (two live coordinates), the first
// ending (after the k-1 junction k-mers) where it would enter the seed again, the
// second at the record's end; on the left the empty first occurrence is covered by
// the second.
TEST(TrieCasesTrace, RepeatedSeed) {
    const size_t k = kDefaultK;
    auto b = blocks_where({ 30, 40, 40 }, 123, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    const std::vector<std::string> seqs { X + P + X + Q }, labels { "F" };
    CaseSpec c { "TraceRepeatedSeed", k, seqs, labels, X };
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        c.k, seqs, labels, DeBruijnGraph::BASIC, true, { 0 });
    trie::StringIndex ix(c.k, DeBruijnGraph::BASIC, seqs, labels, { 0 });
    CaseResult r = check_case_on(*anno, ix, DeBruijnGraph::BASIC, c);
    SeedResult trace;
    check_trace(*anno, ix, c, { "F" }, r.A.at({ "F" }), &trace, c.name);
    const std::string right_p = P + X.substr(0, k - 1);
    const std::string left_p = reversed(X.substr(X.size() - (k - 1)) + P);
    EXPECT_EQ((trie::Leaves{ { right_p, { "F" } }, { Q, { "F" } } }), trie::constrained_leaves(trace, kRight, c.radius));
    EXPECT_EQ((trie::Leaves{ { left_p, { "F" } } }), trie::constrained_leaves(trace, kLeft, c.radius));
    EXPECT_EQ(1u, count_label_ends(trace.arms[kRight], EndReason::REACHED_SEED));
    EXPECT_EQ(1u, count_label_ends(trace.arms[kLeft], EndReason::REACHED_SEED));
}

// The structural rule binds trace support too: the record X·A^(k+3)·Z walks the
// self-loop three times, the second use is a reuse, so the trace ends inside the
// homopolymer and never reaches Z, while k-mer support takes the exit at the first
// visit. The limitation is pinned so that it is a documented one.
TEST(TrieCasesTrace, HomopolymerLimitation) {
    const size_t k = kDefaultK;
    auto b = blocks_where({ 30, 30 }, 124, [](const auto &b) {
        return b[0].back() != 'A' && b[1][0] != 'A';
    });
    const std::string &X = b[0], &Z = b[1];
    const std::string H(k + 3, 'A');
    const std::vector<std::string> seqs { X + H + Z }, labels { "F" };
    CaseSpec c { "TraceHomopolymer", k, seqs, labels, X };
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        c.k, seqs, labels, DeBruijnGraph::BASIC, true, { 0 });
    trie::StringIndex ix(c.k, DeBruijnGraph::BASIC, seqs, labels, { 0 });
    CaseResult r = check_case_on(*anno, ix, DeBruijnGraph::BASIC, c);
    SeedResult trace;
    check_trace(*anno, ix, c, { "F" }, r.A.at({ "F" }), &trace, c.name);
    EXPECT_EQ((trie::Leaves{ { std::string(k + 1, 'A'), { "F" } } }),
              trie::constrained_leaves(trace, kRight, c.radius));
    EXPECT_EQ(1u, count_label_ends(trace.arms[kRight], EndReason::EDGE_REUSE));
    EXPECT_EQ(0u, leaf_walks(trace.arms[kRight]).count(H + Z));
}


/************************* caps and the completeness guarantee *************************/

// On the nested-bubble fixture every size cap, swept over its range, leaves a run that
// is a prefix-subset of the exhaustive trie with a reason for every omission, and whose
// walks agree with the exhaustive trie at every depth up to complete_to_bp (§6.10).
// complete_to_bp never decreases as a cap is raised.
TEST(TrieCasesCaps, SweepAgainstTheCompletenessGuarantee) {
    auto b = blocks_where({ 30, 20, 20, 20, 25, 30 }, 125, [](const auto &b) {
        return b[1][0] != b[4][0] && b[2][0] != b[3][0];
    });
    const std::string &X = b[0], &P = b[1], &Y1 = b[2], &Y2 = b[3], &Q = b[4], &Z = b[5];
    const std::vector<std::string> seqs { X + P + Y1 + Z, X + P + Y2 + Z, X + Q + Z };
    const std::vector<std::string> labels { "A", "B", "C" };
    const uint64_t radius = 200;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kDefaultK, seqs, labels, mode);
        const std::string where = "CapSweep " + mode_name(mode);
        SeedResult A = run(*anno, X, labels, exhaustive(LabelMode::CONSTRAIN, radius));
        ASSERT_EQ(ArmResult::COMPLETE, A.arms[kRight].status) << where;

        auto tuned = [&](const Strategy &st, const std::string &what) {
            const trie::SeedContext ctx { X, mode != DeBruijnGraph::BASIC, st };
            SeedResult t = run(*anno, X, labels, st);
            for (size_t a : { kLeft, kRight }) {
                const ArmResult &arm = t.arms[a];
                trie::check_tuned_subset(A, t, a, what, ctx);
                // complete up to complete_to_bp: the same walks at every depth
                for (uint64_t n = 0; n <= arm.complete_to_bp; ++n) {
                    std::set<std::string> want, got;
                    for (const auto &p : A.arms[a].paths) {
                        const std::string w = trie::walk_of(A.arms[a], p);
                        want.insert(w.substr(0, std::min<size_t>(n, w.size())));
                    }
                    for (const auto &p : arm.paths) {
                        const std::string w = trie::walk_of(arm, p);
                        got.insert(w.substr(0, std::min<size_t>(n, w.size())));
                    }
                    EXPECT_EQ(want, got) << what << " arm " << to_string(arm.arm)
                                         << ": walks differ at depth " << n << " <= complete_to_bp "
                                         << arm.complete_to_bp;
                    if (want != got)
                        break;
                }
                if (arm.status != ArmResult::COMPLETE) {
                    EXPECT_TRUE(arm.cap_trigger.has_value()) << what;
                }
            }
            return t;
        };
        Strategy base;
        base.max_extension_bp = radius;
        base.max_label_branches = Strategy::kUnlimited;
        base.max_splits_per_path = Strategy::kUnlimited;
        base.merge_reconverge = false;

        uint64_t last = 0;
        for (uint64_t steps = 1; steps <= 120; steps += steps < 20 ? 1 : 7) {
            Strategy st = base;
            st.max_steps = steps;
            SeedResult t = tuned(st, where + " max_steps " + std::to_string(steps));
            const uint64_t c = std::min(t.arms[kLeft].complete_to_bp, t.arms[kRight].complete_to_bp);
            EXPECT_GE(c, last) << where << " max_steps " << steps;
            last = c;
        }
        for (size_t live : { 1, 2, 3 }) {
            Strategy st = base;
            st.max_live_paths = live;
            tuned(st, where + " max_live_paths " + std::to_string(live));
            st.on_overflow = Strategy::BEAM;
            tuned(st, where + " beam " + std::to_string(live));
        }
        for (uint64_t out : { 1, 5, 17, 40, 61, 90 }) {
            Strategy st = base;
            st.max_output_bp = out;
            tuned(st, where + " max_output_bp " + std::to_string(out));
        }
        for (size_t paths : { 1, 2 }) {
            Strategy st = base;
            st.max_paths = paths;
            tuned(st, where + " max_paths " + std::to_string(paths));
        }
    }
}


/********************************* dummy nodes *********************************/

// An unmasked DBGSuccinct has '$' dummy k-mers around every record end: they are never
// successors, never recorded, never a branch — the nested-bubble contract holds on it.
TEST(TrieCasesDummies, UnmaskedSuccinct) {
    auto b = blocks_where({ 30, 20, 20, 20, 25, 30 }, 126, [](const auto &b) {
        return b[1][0] != b[4][0] && b[2][0] != b[3][0];
    });
    const std::string &X = b[0], &P = b[1], &Y1 = b[2], &Y2 = b[3], &Q = b[4], &Z = b[5];
    const std::vector<std::string> seqs { X + P + Y1 + Z, X + P + Y2 + Z, X + Q + Z };
    const std::vector<std::string> labels { "A", "B", "C" };
    CaseSpec c { "UnmaskedSuccinct", kDefaultK, seqs, labels, X };
    for (auto mode : all_modes()) {
        std::shared_ptr<DeBruijnGraph> graph;
        std::shared_ptr<DBGSuccinct> base;
        if (mode == DeBruijnGraph::PRIMARY) {
            auto canonical = build_graph<DBGSuccinct>(c.k, seqs, DeBruijnGraph::CANONICAL);
            std::vector<std::string> contigs;
            canonical->call_sequences([&](const std::string &s, const auto &) { contigs.push_back(s); }, 1, true);
            base = std::make_shared<DBGSuccinct>(c.k, DeBruijnGraph::PRIMARY);
            for (const auto &s : contigs) base->add_sequence(s);
            graph = std::make_shared<CanonicalDBG>(base);
        } else {
            base = std::make_shared<DBGSuccinct>(c.k, mode);
            for (const auto &s : seqs) base->add_sequence(s);
            graph = base;
        }
        ASSERT_EQ(nullptr, base->get_mask());
        AnnotatedDBG anno(graph, std::make_unique<annot::ColumnCompressed<>>(base->max_index()));
        for (size_t i = 0; i < seqs.size(); ++i)
            anno.annotate_sequence(std::string(seqs[i]), { labels[i] });
        trie::StringIndex ix(c.k, mode, seqs, labels);
        CaseResult r = check_case_on(anno, ix, mode, c);
        const std::string where = c.name + " " + mode_name(mode);
        expect_leaves(r, { "A", "B", "C" }, kRight,
                      { { P + Y1 + Z, { "A" } }, { P + Y2 + Z, { "B" } }, { Q + Z, { "C" } } }, c.radius, where);
        for (const auto &seg : r.T.arms[kRight].segments)
            EXPECT_EQ(std::string::npos, seg.sequence.find('$')) << where;
    }
}


/********************************** random graphs **********************************/

// Random records at k = 7 from a small pool of blocks, so that repeats, cycles and
// bubbles arise on their own: the full contract on every fixture, every mode, both arms.
TEST(TrieCasesRandom, SmallKRandomRecords) {
    const size_t k = 7;
    std::mt19937 gen(127);
    size_t fixtures_with_cycles = 0, leaves_total = 0;
    for (size_t fixture = 0; fixture < 16; ++fixture) {
        std::vector<size_t> lengths { 12 };
        for (size_t i = 0; i < 6; ++i) lengths.push_back(8 + gen() % 10);
        auto blocks = clean_blocks(lengths, 2000 + fixture, k);
        std::vector<std::string> seqs, labels { "A", "B", "C", "D" };
        for (size_t l = 0; l < labels.size(); ++l) {
            size_t n = 3 + gen() % 3;
            std::vector<size_t> order;
            for (size_t i = 0; i < n; ++i) order.push_back(gen() % blocks.size());
            order[gen() % n] = 0;   // every label contains the seed block
            std::string s;
            for (size_t i : order) s += blocks[i];
            seqs.push_back(s);
        }
        CaseSpec c { "Random" + std::to_string(fixture), k, seqs, labels, blocks[0], 24 };
        for (auto mode : all_modes()) {
            auto anno = build_anno_graph<DBGHashFast, annot::ColumnCompressed<>>(k, seqs, labels, mode);
            trie::StringIndex ix(k, mode, seqs, labels);
            CaseResult r = check_case_on(*anno, ix, mode, c);
            for (size_t a : { kLeft, kRight }) {
                leaves_total += r.T.arms[a].paths.size();
                fixtures_with_cycles += count_blocked(r.T.arms[a], EndReason::EDGE_REUSE) > 0;
            }
        }
    }
    // the fixtures are not trivial
    EXPECT_GT(fixtures_with_cycles, 0u);
    EXPECT_GT(leaves_total, 16u * 3 * 2);
}
