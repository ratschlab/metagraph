#include "gtest/gtest.h"

#include <algorithm>
#include <random>

#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/walker_test_helpers.hpp"
#include "tests/graph/traversal/test_trie_reference.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

/**
 * The walker against an independent model (test_trie_reference.hpp): for every fixture,
 * mode and permitted set, the labeled walks the walker claims with an unlimited branch
 * allowance — every (walk, label) where a label ends, with its reason — are exactly the
 * maximal label-consistent walks the records admit under the walk rule. And on top:
 *   - a run over several labels is the union of the single-label runs;
 *   - a permitted set derived from the seed is its carriers and runs like the explicit list;
 *   - with a branch allowance of 0 every claim is a claim of the reference cut at the
 *     label's first fork, where it ends with BRANCH.
 */

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

std::vector<DeBruijnGraph::Mode> all_modes() {
    return {
        DeBruijnGraph::BASIC,
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

std::string reversed(std::string s) {
    std::reverse(s.begin(), s.end());
    return s;
}

struct Case {
    std::string name;
    size_t k;
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    std::string seed;
    uint64_t radius = 200;
};

// the walk (outward from the seed boundary) of |segment|'s chain up to |at_bp|
std::string outward_walk(const ArmResult &arm, size_t segment, uint64_t at_bp) {
    std::vector<size_t> chain;
    for (size_t s = segment; ; s = arm.segments[s].parents[0]) {
        chain.push_back(s);
        if (arm.segments[s].parents.empty())
            break;
    }
    std::string walk;
    for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
        const std::string &seq = arm.segments[*it].sequence;
        walk += arm.arm == Arm::RIGHT ? seq : reversed(seq);
    }
    return walk.substr(0, at_bp);
}

trie::RefClaims walker_claims(const SeedResult &res, Arm arm) {
    trie::RefClaims out;
    const ArmResult &a = res.arms[static_cast<size_t>(arm)];
    for (const Segment &seg : a.segments) {
        for (const LabelEnd &e : seg.label_ends) {
            auto [it, inserted] = out.emplace(
                    std::make_pair(outward_walk(a, seg.id, e.at_bp), res.label_dict[e.label].name),
                    e.reason);
            EXPECT_TRUE(inserted) << "a (walk, label) claimed twice";
        }
    }
    return out;
}

std::string show(const trie::RefClaims &claims) {
    std::string out;
    for (const auto &[key, reason] : claims) {
        out += "\n  " + key.second + " " + to_string(reason) + " [" + key.first + "]";
    }
    return out;
}

SeedResult walk(const AnnotatedDBG &anno, const Case &c, const std::vector<std::string> &labels,
                size_t branches = Strategy::kUnlimited) {
    LabelOracle oracle(anno);
    Seed seed;
    seed.sequence = c.seed;
    seed.labels = labels;
    Strategy st;
    st.max_label_branches = branches;
    st.max_extension_bp = c.radius;
    // no cap may trip: the comparison is with the whole reference trie
    st.max_splits_per_path = std::numeric_limits<size_t>::max();
    st.max_live_paths = 1'000'000;
    st.max_paths = 1'000'000;
    return traverse_seed(oracle, seed, st);
}

template <class Graph, class Annotation>
void check_case(const Case &c, std::vector<DeBruijnGraph::Mode> modes = all_modes()) {
    for (auto mode : modes) {
        const std::string where = c.name + " mode " + std::to_string(mode);
        auto anno = build_anno_graph<Graph, Annotation>(c.k, c.sequences, c.labels, mode);
        trie::StringIndex ix(c.k, mode, c.sequences, c.labels);
        const std::set<std::string> carriers = ix.carriers(c.seed);
        ASSERT_FALSE(carriers.empty()) << where;
        const std::vector<std::string> permitted(carriers.begin(), carriers.end());

        const SeedResult all = walk(*anno, c, permitted);
        for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
            const trie::RefClaims expected = trie::claims(ix, carriers, c.seed, arm, c.radius);
            const trie::RefClaims got = walker_claims(all, arm);
            EXPECT_EQ(expected, got) << where << " arm " << to_string(arm)
                                     << "\nexpected:" << show(expected) << "\ngot:" << show(got);
            EXPECT_EQ(ArmResult::COMPLETE, all.arms[static_cast<size_t>(arm)].status) << where;

            // the union of the single-label runs
            trie::RefClaims single;
            for (const std::string &l : permitted) {
                for (const auto &claim : walker_claims(walk(*anno, c, { l }), arm)) {
                    single.insert(claim);
                }
            }
            EXPECT_EQ(got, single) << where << " arm " << to_string(arm);

            // limit 0: a label's claim ends at its first fork (BRANCH) or as the reference
            for (const auto &[key, reason] : walker_claims(walk(*anno, c, permitted, 0), arm)) {
                const auto &[w, l] = key;
                if (reason == EndReason::BRANCH) {
                    size_t continuing = 0;
                    for (const auto &[ref_key, ref_reason] : expected) {
                        continuing += ref_key.second == l && ref_key.first.size() > w.size()
                                && ref_key.first.compare(0, w.size(), w) == 0;
                    }
                    EXPECT_GE(continuing, 2u) << where << ": " << l << " [" << w << "]";
                } else {
                    auto it = expected.find(key);
                    ASSERT_NE(expected.end(), it) << where << ": " << l << " [" << w << "]";
                    EXPECT_EQ(it->second, reason) << where;
                }
            }
        }

        // a derived set is the carriers and runs as the explicit list
        const SeedResult derived = walk(*anno, c, {});
        std::set<std::string> names;
        for (const LabelRef &l : derived.label_dict) {
            names.insert(l.name);
        }
        EXPECT_EQ(carriers, names) << where;
        for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
            EXPECT_EQ(walker_claims(all, arm), walker_claims(derived, arm)) << where;
        }
    }
}

template <typename Pair>
class WalkerReference : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::RowFlatAnnotator>
> WalkerReferenceTypes;
TYPED_TEST_SUITE(WalkerReference, WalkerReferenceTypes);

#define CASE_TYPES \
    using Graph = typename TypeParam::first_type; \
    using Annotation = typename TypeParam::second_type

template <class Cond>
std::vector<std::string> blocks_where(const std::vector<size_t> &lengths, uint32_t seed,
                                      const Cond &cond) {
    for (uint32_t s = seed; ; ++s) {
        auto b = clean_blocks(lengths, s);
        if (cond(b))
            return b;
    }
}

TYPED_TEST(WalkerReference, Linear) {
    CASE_TYPES;
    auto b = clean_blocks({ 60, 30, 80 }, 101);
    check_case<Graph, Annotation>({ "Linear", kK, { b[0] + b[1] + b[2] }, { "A" }, b[1] });
}

TYPED_TEST(WalkerReference, SeedAtRecordEnds) {
    CASE_TYPES;
    auto b = clean_blocks({ 50, 30, 50 }, 102);
    const std::string &L = b[0], &S = b[1], &R = b[2];
    check_case<Graph, Annotation>({ "SeedAtRecordStart", kK, { S + R }, { "A" }, S });
    check_case<Graph, Annotation>({ "SeedAtRecordEnd", kK, { L + S }, { "A" }, S });
    check_case<Graph, Annotation>({ "SeedIsRecord", kK, { S }, { "A" }, S });
}

TYPED_TEST(WalkerReference, RadiusAtEveryBoundary) {
    CASE_TYPES;
    auto b = clean_blocks({ 40, 30, 50 }, 103);
    const std::string &L = b[0], &S = b[1], &R = b[2];
    for (uint64_t radius : { 0ul, 1ul, 2ul, R.size() - 1, R.size(), R.size() + 1 }) {
        check_case<Graph, Annotation>({ "Radius" + std::to_string(radius), kK,
                                        { L + S + R }, { "A" }, S, radius });
    }
}

TYPED_TEST(WalkerReference, Fork) {
    CASE_TYPES;
    auto b = fork_blocks(104);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    check_case<Graph, Annotation>({ "Fork", kK, { X + P, X + Q }, { "A", "B" }, X });
    check_case<Graph, Annotation>({ "ForkOneLabel", kK, { X + P, X + Q }, { "A", "A" }, X });
}

TYPED_TEST(WalkerReference, Bubbles) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 15, 22, 30, 30 }, 105, [](const auto &b) { return b[1][0] != b[2][0]; });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Y = b[3], &Z = b[4];
    check_case<Graph, Annotation>({ "UnequalBubble", kK, { X + P + Y + Z, X + Q + Y + Z },
                                    { "A", "B" }, X });
    check_case<Graph, Annotation>({ "UnequalBubbleOneLabel", kK, { X + P + Y + Z, X + Q + Y + Z },
                                    { "A", "A" }, X });
    auto n = blocks_where({ 30, 20, 20, 20, 25, 30 }, 106, [](const auto &b) {
        return b[1][0] != b[4][0] && b[2][0] != b[3][0];
    });
    check_case<Graph, Annotation>({ "NestedBubbles", kK,
                                    { n[0] + n[1] + n[2] + n[5], n[0] + n[1] + n[3] + n[5],
                                      n[0] + n[4] + n[5] },
                                    { "A", "B", "C" }, n[0] });
}

TYPED_TEST(WalkerReference, TipAtTheSeedBoundary) {
    CASE_TYPES;
    auto b = clean_blocks({ 40, 40 }, 107);
    const std::string &X = b[0], &R = b[1];
    char t = 'A';
    while (t == R[0]) {
        t = "ACGT"[(std::string("ACGT").find(t) + 1) % 4];
    }
    check_case<Graph, Annotation>({ "Tip", kK, { X + R, X + t }, { "A", "B" }, X });
}

TYPED_TEST(WalkerReference, HomopolymerSelfLoop) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 30 }, 108, [](const auto &b) {
        return b[0].back() != 'A' && b[1][0] != 'A' && b[0].back() != 'T' && b[1][0] != 'T';
    });
    check_case<Graph, Annotation>({ "HomopolymerSelfLoop", kK,
                                    { b[0] + std::string(kK + 3, 'A') + b[1] }, { "A" }, b[0] });
}

TYPED_TEST(WalkerReference, CycleJunction) {
    CASE_TYPES;
    auto b = blocks_where({ 30, kK, 40, 30 }, 109, [](const auto &b) { return b[2][0] != b[3][0]; });
    const std::string &L = b[0], &J = b[1], &C = b[2], &E = b[3];
    check_case<Graph, Annotation>({ "CycleJunction", kK, { L + J + C + J + E }, { "A" }, L });
}

TYPED_TEST(WalkerReference, Circle) {
    CASE_TYPES;
    std::string S = clean_block(150, 110);
    check_case<Graph, Annotation>({ "Circle", kK, { S + S.substr(0, kK - 1 + 10) }, { "A" },
                                    S.substr(50, 40) });
}

TYPED_TEST(WalkerReference, Repeats) {
    CASE_TYPES;
    auto b = blocks_where({ 30, 30, 30, 30, 30 }, 111, [](const auto &b) {
        return b[2][0] != b[4][0] && b[0].back() != b[3].back();
    });
    const std::string &X = b[0], &R = b[1], &Y = b[2], &Z = b[3], &W = b[4];
    check_case<Graph, Annotation>({ "RepeatInTwoRecordsAndAChimera", kK,
                                    { X + R + Y, Z + R + W, X + R + W }, { "A", "B", "C" }, R });
    auto f = fork_blocks(112);
    check_case<Graph, Annotation>({ "RepeatedSeedInOneRecord", kK,
                                    { f[0] + f[1] + f[0] + f[2] }, { "A" }, f[0] });
}

TYPED_TEST(WalkerReference, ReverseComplementRecord) {
    CASE_TYPES;
    auto b = fork_blocks(115);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    check_case<Graph, Annotation>({ "ReverseComplementRecord", kK, { X + P, rc(X + Q) },
                                    { "A", "B" }, X });
}

TYPED_TEST(WalkerReference, Hairpin) {
    CASE_TYPES;
    const std::string X = "AACGTACGTT";
    ASSERT_EQ(X, rc(X));
    std::string L, R;
    for (uint32_t seed = 116; ; ++seed) {
        L = clean_block(40, seed);
        R = clean_block(40, seed + 100);
        if (R[0] != rc(std::string(1, L.back()))[0] && is_clean(L + X.substr(0, 3))
                && is_clean(X.substr(7) + R))
            break;
    }
    check_case<Graph, Annotation>({ "Hairpin", kK, { L + X + R }, { "A" }, L.substr(0, 30) });
}

// random records at k = 7 over a small alphabet of shared pieces: repeats, cycles, both
// strands and forks everywhere
TYPED_TEST(WalkerReference, RandomGraphs) {
    CASE_TYPES;
    std::mt19937 gen(118);
    for (size_t trial = 0; trial < 12; ++trial) {
        std::vector<std::string> pieces;
        for (size_t i = 0; i < 6; ++i) {
            pieces.push_back(random_seq(6 + gen() % 10, gen()));
        }
        std::vector<std::string> seqs, labels;
        for (size_t r = 0; r < 4; ++r) {
            std::string s;
            for (size_t j = 0; j < 5; ++j) {
                s += pieces[gen() % pieces.size()];
            }
            seqs.push_back(s);
            labels.push_back(std::string(1, 'A' + gen() % 3));
        }
        const std::string &host = seqs[gen() % seqs.size()];
        const size_t at = gen() % (host.size() - 10);
        check_case<Graph, Annotation>({ "Random" + std::to_string(trial), 7, seqs, labels,
                                        host.substr(at, 10), 25 });
    }
}

} // namespace
