#include "gtest/gtest.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <random>
#include <set>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/test_trie_oracle.hpp"
#include "tests/graph/traversal/test_trie_reference.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/seq_tools/reverse_complement.hpp"


/*
 * The exhaustive trie oracle (spec §6.9 / §6.10). First the smoke tests: the annotate
 * label mode, the exhaustive preset, the BFS completeness boundary and the per-node
 * label cap. Then the verification contract itself (suite TrieOracle, helpers in
 * test_trie_oracle.hpp): leaves(A) == leaves(E) on a fixture with a bubble, a tandem
 * repeat and a tip over basic / canonical / primary graphs and two annotation
 * representations, the tuned-subset property with its every-omission-has-a-reason
 * clause, the partial-level exclusion of the completeness guarantee, and the cut
 * recorded lists that make the oracle unusable (and are reported as such). Last the
 * checkers' own regressions: results corrupted the way a reviewer did must fail them.
 */
namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

const size_t kK = 11;
const size_t kLeft = static_cast<size_t>(Arm::LEFT);
const size_t kRight = static_cast<size_t>(Arm::RIGHT);

std::string rc(std::string s) { ::reverse_complement(s); return s; }

std::string random_seq(size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(len, 'A');
    for (char &c : s) c = "ACGT"[gen() % 4];
    return s;
}

// no repeated k-mer or (k-1)-mer in either orientation, no RC palindromes
bool is_clean(const std::string &s) {
    for (size_t len : { kK - 1, kK }) {
        std::set<std::string> seen;
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string km = s.substr(i, len);
            if (seen.count(km) || seen.count(rc(km)))
                return false;
            seen.insert(km);
        }
    }
    for (size_t len : { kK - 1, kK + 1 }) {
        for (size_t i = 0; i + len <= s.size(); ++i) {
            std::string x = s.substr(i, len);
            if (rc(x) == x)
                return false;
        }
    }
    return true;
}

std::vector<std::string> clean_blocks(const std::vector<size_t> &lengths, uint32_t seed) {
    size_t total = 0;
    for (size_t l : lengths) total += l;
    std::string master;
    for (uint32_t s = seed; ; ++s) {
        master = random_seq(total, s);
        if (is_clean(master))
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

// blocks {X, P, Q} with P[0] != Q[0]: X+P and X+Q fork right after the seed X
std::vector<std::string> fork_blocks(uint32_t seed, std::vector<size_t> lengths = { 30, 40, 40 }) {
    for (uint32_t s = seed; ; ++s) {
        auto b = clean_blocks(lengths, s);
        if (b[1][0] != b[2][0])
            return b;
    }
}

Strategy exhaustive(LabelMode mode) {
    Strategy st;
    st.exhaustive = true;
    st.label_mode = mode;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.merge_reconverge = false;
    return st;
}

SeedResult run(const AnnotatedDBG &anno, const std::string &seq,
               const std::vector<std::string> &labels, const Strategy &st) {
    LabelOracle oracle(anno);
    Seed seed;
    seed.sequence = seq;
    seed.labels = labels;
    return traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
}

std::vector<std::string> names(const SeedResult &r, const std::vector<LabelId> &ids) {
    std::vector<std::string> out;
    for (LabelId l : ids) out.push_back(r.label_dict[l].name);
    std::sort(out.begin(), out.end());
    return out;
}

// leaf flank -> sorted names of the labels present on every node of that flank, read
// from the RECORDED sets of an annotate run (no walker label machinery involved)
std::map<std::string, std::vector<std::string>> structural_leaves(const SeedResult &r, size_t a) {
    const ArmResult &arm = r.arms[a];
    std::map<std::string, std::vector<std::string>> out;
    for (const auto &path : arm.paths) {
        std::vector<LabelId> alive = arm.segments[path.segments.front()].labels_start;
        for (size_t s : path.segments) {
            for (const auto &run : arm.segments[s].label_sets) {
                EXPECT_FALSE(run.truncated());
                std::vector<LabelId> still;
                std::set_intersection(alive.begin(), alive.end(), run.labels.begin(),
                                      run.labels.end(), std::back_inserter(still));
                alive.swap(still);
            }
        }
        out[spell_path(arm, path)] = names(r, alive);
    }
    return out;
}

} // namespace


// Annotate mode follows every structural successor and records what is there: on a
// fork carried by two labels, the structural trie has both branches, each branch's
// recorded sets name exactly the label that carries it, and the constrained walk over
// the same seed yields the same leaves with the same labels (the §6.9 contract in
// miniature).
TEST(Trie, AnnotateRecordsWhatConstrainFilters) {
    auto b = fork_blocks(3);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    std::vector<DeBruijnGraph::Mode> modes { DeBruijnGraph::BASIC };
#if ! _PROTEIN_GRAPH
    modes.push_back(DeBruijnGraph::CANONICAL);
    modes.push_back(DeBruijnGraph::PRIMARY);
#endif
    for (auto mode : modes) {
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, { X + P, X + Q }, { "A", "B" }, mode);

        Strategy st = exhaustive(LabelMode::ANNOTATE);
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 100;
        auto t = run(*anno, X, {}, st);
        EXPECT_EQ(0u, t.num_seed_labels);
        EXPECT_FALSE(t.labels_from_seed);
        ASSERT_EQ(2u, t.label_dict.size());
        const ArmResult &arm = t.arms[kRight];
        EXPECT_EQ(ArmResult::COMPLETE, arm.status);
        EXPECT_EQ(100u, arm.complete_to_bp);
        EXPECT_EQ(0u, arm.nodes_labels_truncated);
        EXPECT_EQ(2u, arm.max_labels_at_node);
        // the root's entry node is the seed boundary, carried by both labels
        ASSERT_FALSE(arm.segments.empty());
        EXPECT_EQ((std::vector<std::string>{ "A", "B" }), names(t, arm.segments[0].labels_start));
        EXPECT_EQ(0u, arm.segments[0].length_bp);
        EXPECT_TRUE(arm.segments[0].label_sets.empty());
        ASSERT_EQ(1u, arm.splits.size());
        const Split &split = arm.splits[0];
        EXPECT_EQ(0u, split.at_bp);
        EXPECT_FALSE(split.ambiguous);
        EXPECT_EQ(2u, split.labels_before);
        ASSERT_EQ(2u, split.branches.size());
        std::set<char> chars;
        for (const auto &br : split.branches) {
            chars.insert(br.ch);
            EXPECT_EQ(1u, br.labels_distinct);
            EXPECT_EQ(1u, br.labels.size());
            EXPECT_EQ(br.labels, arm.segments[br.segment].labels_start);
        }
        EXPECT_EQ((std::set<char>{ P[0], Q[0] }), chars);
        // every leaf is a structural end with a path reason and no label ends
        ASSERT_EQ(2u, arm.paths.size());
        for (const auto &path : arm.paths) {
            ASSERT_TRUE(path.path_reason.has_value());
            EXPECT_EQ(EndReason::DEAD_END, *path.path_reason);
            EXPECT_TRUE(path.end_labels.empty());
            EXPECT_FALSE(path.continuation.has_value());
            // one run per segment: the label set never changes along a branch
            for (size_t s : path.segments) {
                const Segment &seg = arm.segments[s];
                if (!seg.length_bp) continue;
                ASSERT_EQ(1u, seg.label_sets.size());
                EXPECT_EQ(seg.from_bp, seg.label_sets[0].from_bp);
                EXPECT_EQ(seg.from_bp + seg.length_bp, seg.label_sets[0].to_bp);
                EXPECT_EQ(seg.label_sets[0].labels, seg.labels_start);
                EXPECT_EQ(seg.label_sets[0].labels, seg.labels_end);
            }
        }
        auto structural = structural_leaves(t, kRight);
        ASSERT_EQ(2u, structural.size());
        EXPECT_EQ((std::vector<std::string>{ "A" }), structural.at(P));
        EXPECT_EQ((std::vector<std::string>{ "B" }), structural.at(Q));
        // the per-label summary reads off the recorded sets
        for (LabelId l = 0; l < 2; ++l) {
            EXPECT_EQ(P.size(), t.label_summary[l][kRight].direct_bp);
            EXPECT_EQ(P.size(), t.label_summary[l][kRight].reach_bp);
            EXPECT_TRUE(t.label_summary[l][kRight].runs.empty());
        }

        // the constrained exhaustive walk agrees leaf by leaf
        Strategy sc = exhaustive(LabelMode::CONSTRAIN);
        sc.direction = Strategy::RIGHT;
        sc.max_extension_bp = 100;
        auto a = run(*anno, X, { "A", "B" }, sc);
        const ArmResult &carm = a.arms[kRight];
        EXPECT_EQ(ArmResult::COMPLETE, carm.status);
        EXPECT_EQ(100u, carm.complete_to_bp);
        std::map<std::string, std::vector<std::string>> constrained;
        for (const auto &path : carm.paths) {
            std::vector<LabelId> ids;
            for (const auto &e : path.end_labels) ids.push_back(e.label);
            constrained[spell_path(carm, path)] = names(a, ids);
        }
        EXPECT_EQ(structural, constrained) << "mode " << mode;
        // and its trie view reports the same branches
        ASSERT_EQ(1u, carm.splits.size());
        EXPECT_EQ(2u, carm.splits[0].labels_before);
        ASSERT_EQ(2u, carm.splits[0].branches.size());
        for (const auto &br : carm.splits[0].branches) {
            EXPECT_EQ(1u, br.labels_distinct);
            EXPECT_EQ(br.labels, carm.segments[br.segment].labels_start);
        }
    }
}

// The preset refuses what would silently prune it, in both modes; annotate mode
// refuses the label machinery and a seed label list.
TEST(Trie, ExhaustiveRejectsConflictingKnobs) {
    auto b = fork_blocks(4);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { b[0] + b[1] }, { "A" }, DeBruijnGraph::BASIC);
    auto expect_reject = [&](const Strategy &st, const std::vector<std::string> &labels,
                             const std::string &knob) {
        try {
            run(*anno, b[0], labels, st);
            FAIL() << "accepted a strategy conflicting on " << knob;
        } catch (const std::invalid_argument &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find(knob)) << e.what();
        }
    };
    for (LabelMode mode : { LabelMode::CONSTRAIN, LabelMode::ANNOTATE }) {
        const std::vector<std::string> labels = mode == LabelMode::CONSTRAIN
            ? std::vector<std::string>{ "A" } : std::vector<std::string>{};
        // the preset itself is fine
        EXPECT_NO_THROW(run(*anno, b[0], labels, exhaustive(mode)));

        Strategy st = exhaustive(mode);
        st.max_label_branches = 3;
        expect_reject(st, labels, "max_label_branches");
        st = exhaustive(mode);
        st.merge_reconverge = true;
        expect_reject(st, labels, "on_reconverge");
        st = exhaustive(mode);
        st.max_splits_per_path = 5;
        expect_reject(st, labels, "max_splits_per_path");
        st = exhaustive(mode);
        st.min_successor_labels = 2;
        expect_reject(st, labels, "min_successor_labels");
        st = exhaustive(mode);
        st.min_successor_fraction = 0.5;
        expect_reject(st, labels, "min_successor_fraction");
        st = exhaustive(mode);
        st.min_live_labels = 2;
        expect_reject(st, labels, "min_live_labels");
        st = exhaustive(mode);
        st.on_overflow = Strategy::BEAM;
        expect_reject(st, labels, "on_overflow");
    }
    // annotate mode: no permitted set, so none of the label machinery
    Strategy st = exhaustive(LabelMode::ANNOTATE);
    st.extra = { "A" };
    expect_reject(st, {}, "labels.extra");
    st = exhaustive(LabelMode::ANNOTATE);
    st.loss_budget = 1;
    expect_reject(st, {}, "loss_budget");
    st = exhaustive(LabelMode::ANNOTATE);
    st.support = Support::TRACE;
    expect_reject(st, {}, "support");
    st = exhaustive(LabelMode::ANNOTATE);
    expect_reject(st, { "A" }, "seeds[].labels");
    {
        LabelOracle oracle(*anno);
        Seed seed;
        seed.sequence = b[0];
        EXPECT_THROW(traverse_seed(oracle, seed, exhaustive(LabelMode::ANNOTATE),
                                   LabelChangeCost::constant(1)),
                     std::invalid_argument);
    }
    // a pairwise (table) cost under the preset: derive() keeps only the cheapest
    // max_switch_sources sources, and a target reachable only from a cut source is
    // not entered — a walk pruned with no path-level reason. The default (64) is
    // refused, "unlimited" is required; CONSTANT and FORBID never cut, so the knob
    // does not bind there and is accepted
    {
        LabelOracle oracle(*anno);
        Seed seed;
        seed.sequence = b[0];
        seed.labels = { "A" };
        const auto table = LabelChangeCost::table({}, kInfiniteLoss);
        Strategy st = exhaustive(LabelMode::CONSTRAIN);
        EXPECT_EQ(64u, st.max_switch_sources);
        try {
            traverse_seed(oracle, seed, st, table);
            FAIL() << "accepted a bounded max_switch_sources under exhaustive with a table cost";
        } catch (const std::invalid_argument &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find("max_switch_sources"))
                << e.what();
        }
        st.max_switch_sources = Strategy::kUnlimited;
        EXPECT_NO_THROW(traverse_seed(oracle, seed, st, table));
        st.max_switch_sources = 64;
        EXPECT_NO_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::constant(1)));
        EXPECT_NO_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()));
    }
    // the plain (non-exhaustive) constrain default still accepts its own defaults
    EXPECT_NO_THROW(run(*anno, b[0], { "A" }, Strategy()));
}

// A size cap trips between two heads of a level: the level is partial, complete_to_bp
// is the last complete depth, the status is never "complete", and every walk up to
// that depth is present.
TEST(Trie, TrippedCapReportsTheCompleteDepth) {
    auto b = fork_blocks(5);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P, X + Q }, { "A", "B" }, DeBruijnGraph::BASIC);
    for (LabelMode mode : { LabelMode::ANNOTATE, LabelMode::CONSTRAIN }) {
        const std::vector<std::string> labels = mode == LabelMode::CONSTRAIN
            ? std::vector<std::string>{ "A", "B" } : std::vector<std::string>{};
        // two walks of 40 bases; a budget of 15 steps ends at a partial level
        Strategy st = exhaustive(mode);
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 100;
        st.max_steps = 15;
        auto r = run(*anno, X, labels, st);
        const ArmResult &arm = r.arms[kRight];
        EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
        ASSERT_TRUE(arm.cap_trigger.has_value());
        EXPECT_EQ(EndReason::MAX_STEPS, arm.cap_trigger->reason);
        EXPECT_LT(arm.complete_to_bp, 100u);
        // 15 steps over two heads: 7 full levels (14 steps), the 8th level expands one
        // head and trips on the other, so walks of length 7 are all present
        EXPECT_EQ(7u, arm.complete_to_bp);
        std::set<std::string> walks;
        for (const auto &path : arm.paths) {
            std::string flank = spell_path(arm, path);
            EXPECT_GE(flank.size(), arm.complete_to_bp);
            walks.insert(flank.substr(0, arm.complete_to_bp));
            ASSERT_TRUE(path.path_reason.has_value());
            EXPECT_EQ(EndReason::MAX_STEPS, *path.path_reason);
            ASSERT_TRUE(path.continuation.has_value());
        }
        EXPECT_EQ((std::set<std::string>{ P.substr(0, 7), Q.substr(0, 7) }), walks);

        // the other arm has nothing to do and is complete to the radius
        auto both = run(*anno, X, labels, exhaustive(mode));
        EXPECT_EQ(ArmResult::COMPLETE, both.arms[kLeft].status);
        EXPECT_EQ(both.arms[kLeft].complete_to_bp, Strategy().max_extension_bp);

        // a per-arm cap: output capped at 10 bases per arm, 5 complete levels
        st = exhaustive(mode);
        st.direction = Strategy::RIGHT;
        st.max_output_bp = 10;
        r = run(*anno, X, labels, st);
        EXPECT_EQ(ArmResult::TRUNCATED, r.arms[kRight].status);
        EXPECT_EQ(5u, r.arms[kRight].complete_to_bp);
        EXPECT_EQ(EndReason::MAX_OUTPUT, r.arms[kRight].cap_trigger->reason);

        // reaching the radius is complete
        st = exhaustive(mode);
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 12;
        r = run(*anno, X, labels, st);
        EXPECT_EQ(ArmResult::COMPLETE, r.arms[kRight].status);
        EXPECT_EQ(12u, r.arms[kRight].complete_to_bp);
        for (const auto &path : r.arms[kRight].paths) {
            EXPECT_EQ(12u, path.length_bp);
            EXPECT_EQ(EndReason::MAX_EXTENSION, *path.path_reason);
        }
    }
}

// The per-node label cap never hides that it cut a list: the true count and the cut
// are reported on the run, the branch and the arm.
TEST(Trie, AnnotateReportsLabelListTruncation) {
    auto b = fork_blocks(6);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P, X + Q, X + P }, { "A", "B", "C" }, DeBruijnGraph::BASIC);
    Strategy st = exhaustive(LabelMode::ANNOTATE);
    st.direction = Strategy::RIGHT;
    st.max_labels_per_node = 1;
    auto r = run(*anno, X, {}, st);
    const ArmResult &arm = r.arms[kRight];
    EXPECT_EQ(ArmResult::COMPLETE, arm.status);
    EXPECT_EQ(3u, arm.max_labels_at_node);
    EXPECT_GT(arm.nodes_labels_truncated, 0u);
    // the boundary (3 labels) and the P branch (2 labels) are cut, the Q branch is not
    EXPECT_EQ(1u, arm.segments[0].labels_start.size());
    ASSERT_EQ(1u, arm.splits.size());
    EXPECT_EQ(3u, arm.splits[0].labels_before);
    for (const auto &br : arm.splits[0].branches) {
        EXPECT_EQ(1u, br.labels.size());
        EXPECT_EQ(br.ch == P[0] ? 2u : 1u, br.labels_distinct);
        const Segment &seg = arm.segments[br.segment];
        ASSERT_EQ(1u, seg.label_sets.size());
        EXPECT_EQ(br.labels_distinct, seg.label_sets[0].labels_total);
        EXPECT_EQ(br.ch == P[0], seg.label_sets[0].truncated());
    }
    // with a cap that fits, nothing is cut and the counts agree with the lists
    st.max_labels_per_node = 3;
    r = run(*anno, X, {}, st);
    EXPECT_EQ(0u, r.arms[kRight].nodes_labels_truncated);
    EXPECT_EQ(3u, r.arms[kRight].segments[0].labels_start.size());
    for (const auto &seg : r.arms[kRight].segments) {
        for (const auto &run : seg.label_sets) {
            EXPECT_FALSE(run.truncated());
            EXPECT_EQ(run.labels.size(), run.labels_total);
        }
    }
}

// The dictionary and the recorded sets do not depend on how far ahead rows are
// prefetched, and the same request gives the same result twice.
TEST(Trie, AnnotateIsDeterministicAndBatchInvariant) {
    auto b = fork_blocks(7, { 30, 60, 60 });
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    auto anno = build_anno_graph<DBGHashFast, annot::ColumnCompressed<>>(
            kK, { X + P, X + Q, X + P.substr(0, 20) + Q.substr(20) }, { "A", "B", "C" },
            DeBruijnGraph::BASIC);
    auto reduce = [](const SeedResult &r) {
        std::string s;
        for (const auto &l : r.label_dict) s += l.name + ",";
        for (const ArmResult &arm : r.arms) {
            for (const auto &seg : arm.segments) {
                s += "\n" + std::to_string(seg.id) + " " + seg.sequence;
                for (const auto &run : seg.label_sets) {
                    s += " [" + std::to_string(run.from_bp) + "," + std::to_string(run.to_bp) + ")";
                    for (LabelId l : run.labels) s += std::to_string(l) + ";";
                    s += "/" + std::to_string(run.labels_total);
                }
            }
            s += "\n" + std::to_string(arm.complete_to_bp) + " " + std::to_string(arm.status);
        }
        return s;
    };
    std::string first;
    for (size_t batch : { 1u, 2u, 7u, 64u }) {
        Strategy st = exhaustive(LabelMode::ANNOTATE);
        st.batch_kmers = batch;
        auto r = run(*anno, X, {}, st);
        EXPECT_EQ(3u, r.label_dict.size());
        if (first.empty()) {
            first = reduce(r);
        } else {
            EXPECT_EQ(first, reduce(r)) << "batch_kmers " << batch;
        }
    }
}


/*****************************************************************************
 *                       The verification contract (§6.9)                    *
 *****************************************************************************/

namespace {

std::vector<DeBruijnGraph::Mode> all_modes() {
    return {
        DeBruijnGraph::BASIC,
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

const uint32_t kFixtureSeed = 11;
const uint64_t kRadius = 200;
const std::set<std::string> kAllLabels { "A", "B", "C", "D", "S" };

std::vector<std::string> as_list(const std::set<std::string> &s) {
    return { s.begin(), s.end() };
}

std::string reversed(std::string s) { std::reverse(s.begin(), s.end()); return s; }

Strategy exhaustive_at(LabelMode mode, uint64_t radius) {
    Strategy st = exhaustive(mode);
    st.max_extension_bp = radius;
    return st;
}

/*
 * The oracle fixture (brief: a bubble, a tandem repeat and a short tip), k = 11. The
 * blocks come from one clean master (no k-mer or (k-1)-mer repeated in either
 * orientation, no RC palindromes), so the only structure is the one built in:
 *
 *   left of the seed:   V · R · R · X            R is exactly k bases: a tandem repeat,
 *                                                i.e. a cycle of length k through node R
 *   right of the seed:  X · P · Y · Z            A, C      P | Q: a bubble closing in Y
 *                       X · Q · Y · Z            B, C      (C follows BOTH branches)
 *                       X · P · Y · Z[:m] · T    C, D      T: a 6 bp tip off Z at m
 *                       X · P · Y[:10]           S         a label that stops early
 *
 * Every label but S also carries the left context V · R · R · X. With k-mer support a
 * label carries a walk iff it holds every k-mer of it, so C (which holds Q·Y, Y·Z[:m]
 * and Z[:m]·T) carries the walk Q·Y·Z[:m]·T that no record spells, and every label
 * carries both the single and the double lap of the repeat.
 */
struct OracleFixture {
    static constexpr size_t m = 20;
    std::string V, R, X, P, Q, Y, Z, T;
    std::vector<std::string> sequences, labels;

    explicit OracleFixture(uint32_t seed) {
        std::vector<std::string> b;
        for (uint32_t s = seed; ; ++s) {
            b = clean_blocks({ 25, kK, 30, 15, 15, 25, 40, 6 }, s);
            // real forks at the bubble and the tip; and the repeat's cycle nodes
            // (rotations of R, absent from the master) must not coincide with the
            // context k-mers V[-1]·R[:k-1] and R[1:]·X[0], which they do exactly when
            // V ends like R or X starts like R
            if (b[3][0] != b[4][0] && b[7][0] != b[6][m]
                    && b[0].back() != b[1].back() && b[2][0] != b[1][0])
                break;
        }
        V = b[0]; R = b[1]; X = b[2]; P = b[3]; Q = b[4]; Y = b[5]; Z = b[6]; T = b[7];
        const std::string left = V + R + R + X;
        auto add = [&](const std::string &label, const std::string &seq) {
            sequences.push_back(seq);
            labels.push_back(label);
        };
        add("A", left + P + Y + Z);
        add("B", left + Q + Y + Z);
        add("C", left + P + Y + Z);
        add("C", left + Q + Y + Z);
        add("C", left + P + Y + Z.substr(0, m) + T);
        add("D", left + P + Y + Z.substr(0, m) + T);
        add("S", X + P + Y.substr(0, 10));
    }

    // the k-mers with more than one successor or predecessor: the bubble's open and
    // close, the tip's fork and the repeat node
    std::set<std::string> branching() const {
        return { X.substr(X.size() - kK), Y.substr(0, kK), Z.substr(m - kK, kK), R };
    }
    // the structural walks, outward from the seed boundary
    std::string pyz() const { return P + Y + Z; }
    std::string pyt() const { return P + Y + Z.substr(0, m) + T; }
    std::string qyz() const { return Q + Y + Z; }
    std::string qyt() const { return Q + Y + Z.substr(0, m) + T; }
    std::set<std::string> right_walks() const { return { pyz(), pyt(), qyz(), qyt() }; }
    std::string one_lap() const { return reversed(V + R); }
    std::string two_laps() const { return reversed(V + R + R); }
    std::set<std::string> left_walks() const { return { one_lap(), two_laps() }; }
};

// the graph has exactly the branching the fixture builds in (both orientations fold
// onto one k-mer in canonical and primary mode)
bool topology_is(const AnnotatedDBG &anno, DeBruijnGraph::Mode mode,
                 const std::set<std::string> &expected) {
    const DeBruijnGraph &g = anno.get_graph();
    auto norm = [&](std::string s) {
        if (mode != DeBruijnGraph::BASIC)
            s = std::min(s, rc(s));
        return s;
    };
    std::set<std::string> found, want;
    g.call_nodes([&](node_index n) {
        if (g.outdegree(n) > 1 || g.indegree(n) > 1) {
            std::string s = g.get_node_sequence(n);
            if (s.find('$') == std::string::npos)
                found.insert(norm(s));
        }
    });
    for (const auto &s : expected) want.insert(norm(s));
    EXPECT_EQ(want, found) << "unexpected graph structure in mode " << mode;
    return want == found;
}

// the walks of exactly |n| bases that |arm| reached
std::set<std::string> walks_at(const ArmResult &arm, uint64_t n) {
    std::set<std::string> out;
    for (const auto &p : arm.paths) {
        std::string w = trie::walk_of(arm, p);
        if (w.size() >= n)
            out.insert(w.substr(0, n));
    }
    return out;
}

const Split* split_at(const ArmResult &arm, uint64_t at) {
    for (const Split &s : arm.splits) {
        if (s.at_bp == at)
            return &s;
    }
    return nullptr;
}

size_t count_blocked(const ArmResult &arm, EndReason reason) {
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &ev : seg.events)
            n += ev.type == EventType::BLOCKED && ev.reason == reason;
    }
    return n;
}

// runs of the label called |name| that ended with |reason| on arm |a|
size_t count_ends_named(const SeedResult &r, size_t a, const std::string &name, EndReason reason) {
    size_t n = 0;
    for (const auto &run : r.arms[a].runs) {
        n += run.ended && run.end_reason == reason && r.label_dict[run.label].name == name;
    }
    return n;
}

std::set<std::string> leaf_walks(const ArmResult &arm) {
    std::set<std::string> out;
    for (const auto &p : arm.paths) out.insert(trie::walk_of(arm, p));
    return out;
}

} // namespace


template <typename Pair>
class TrieOracle : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::RowFlatAnnotator>
> TrieOracleTypes;
TYPED_TEST_SUITE(TrieOracle, TrieOracleTypes);


// Brief test 1. T = the structural trie with labels recorded; E = its walks filtered
// by the recorded sets with test-side code; A = the label-constrained exhaustive walk.
// leaves(A) == leaves(E) for several permitted sets, on both arms, in every graph mode.
TYPED_TEST(TrieOracle, ConstrainedTrieEqualsTheFilteredStructuralTrie) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const OracleFixture f(kFixtureSeed);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.sequences, f.labels, mode);
        ASSERT_TRUE(topology_is(*anno, mode, f.branching()));
        const std::string where = "mode " + std::to_string(mode);

        auto T = run(*anno, f.X, {}, exhaustive_at(LabelMode::ANNOTATE, kRadius));
        ASSERT_EQ(5u, T.label_dict.size());
        for (size_t a : { kLeft, kRight }) {
            EXPECT_EQ(ArmResult::COMPLETE, T.arms[a].status) << where;
            EXPECT_EQ(kRadius, T.arms[a].complete_to_bp);
            EXPECT_EQ(0u, T.arms[a].nodes_labels_truncated);
            EXPECT_EQ(5u, T.arms[a].max_labels_at_node);
            EXPECT_TRUE(trie::is_trie(T.arms[a]));
            for (const auto &p : T.arms[a].paths) {
                ASSERT_TRUE(p.path_reason.has_value());
                EXPECT_EQ(EndReason::DEAD_END, *p.path_reason) << where;
                EXPECT_TRUE(p.end_labels.empty());
            }
        }
        // the structural trie: bubble x tip on the right, one or two laps of the
        // repeat on the left, the second lap's re-entry blocked by edge reuse
        EXPECT_EQ(f.right_walks(), leaf_walks(T.arms[kRight])) << where;
        EXPECT_EQ(f.left_walks(), leaf_walks(T.arms[kLeft])) << where;
        EXPECT_EQ(1u, count_blocked(T.arms[kLeft], EndReason::EDGE_REUSE)) << where;
        EXPECT_EQ(3u, T.arms[kRight].splits.size());
        EXPECT_EQ(1u, T.arms[kLeft].splits.size());

        // the trie view of the bubble: C follows both branches, so the branch counts
        // (4 + 2) exceed labels_before (5) and the split is not "ambiguous" (no lineage)
        const Split *bubble = split_at(T.arms[kRight], 0);
        ASSERT_NE(nullptr, bubble);
        EXPECT_FALSE(bubble->ambiguous);
        EXPECT_EQ(5u, bubble->labels_before);
        ASSERT_EQ(2u, bubble->branches.size());
        std::map<char, size_t> counts;
        for (const auto &br : bubble->branches) {
            counts[br.ch] = br.labels_distinct;
            EXPECT_EQ(br.labels_distinct, br.labels.size());
        }
        EXPECT_EQ(4u, counts[f.P[0]]);
        EXPECT_EQ(2u, counts[f.Q[0]]);

        // the contract
        const std::vector<std::set<std::string>> permitted_sets {
            kAllLabels, { "A", "D" }, { "S" }, { "C" }, { "B", "D" },
        };
        for (const auto &P : permitted_sets) {
            std::string what = where + " P={";
            for (const auto &l : P) what += " " + l;
            what += " }";
            auto A = run(*anno, f.X, as_list(P), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
            EXPECT_EQ(P.size(), A.num_seed_labels);
            for (size_t a : { kLeft, kRight }) {
                EXPECT_EQ(ArmResult::COMPLETE, A.arms[a].status) << what;
                auto v = trie::verify(T, A, a, P);
                EXPECT_EQ(kRadius, v.depth);
                EXPECT_FALSE(v.cut);
                trie::expect_equal(v.expected, v.actual, what + " arm " + std::to_string(a));
                EXPECT_FALSE(v.expected.empty()) << what;
            }
            if (P == kAllLabels) {
                // the leaves, literally
                trie::Leaves right {
                    { f.pyz(), { "A", "C" } }, { f.pyt(), { "C", "D" } },
                    { f.qyz(), { "B", "C" } }, { f.qyt(), { "C" } },
                };
                trie::Leaves left {
                    { f.one_lap(), { "A", "B", "C", "D" } },
                    { f.two_laps(), { "A", "B", "C", "D" } },
                };
                EXPECT_EQ(right, trie::constrained_leaves(A, kRight, kRadius)) << what;
                EXPECT_EQ(left, trie::constrained_leaves(A, kLeft, kRadius)) << what;
                // S stops inside Y, D at the tip: label ends on paths that go on
                EXPECT_EQ(1u, count_ends_named(A, kRight, "S", EndReason::LABEL_LOST));
                EXPECT_EQ(1u, count_ends_named(A, kLeft, "S", EndReason::LABEL_LOST));
                // the CLAIMS, literally: the leaves plus S's own, which ends inside the
                // P stretch that A, C and D continue (right) and at the boundary (left).
                // This is what the contract compares, on both sides: the oracle E keeps
                // S's claim although it is a prefix of the PYZ and PYT claims, and the
                // walker's label end for S is at exactly that depth.
                trie::Leaves right_claims = right, left_claims = left;
                right_claims[f.P + f.Y.substr(0, 10)] = { "S" };
                left_claims[""] = { "S" };
                EXPECT_EQ(right_claims, trie::constrained_claims(A, kRight, kRadius)) << what;
                EXPECT_EQ(left_claims, trie::constrained_claims(A, kLeft, kRadius)) << what;
                EXPECT_EQ(right_claims, trie::verify(T, A, kRight, P).expected) << what;
                EXPECT_EQ(left_claims, trie::verify(T, A, kLeft, P).expected) << what;
                // the constrained bubble IS ambiguous (C), with the same branch counts
                const Split *cb = split_at(A.arms[kRight], 0);
                ASSERT_NE(nullptr, cb);
                EXPECT_TRUE(cb->ambiguous);
                EXPECT_EQ(5u, cb->labels_before);
                std::map<char, size_t> cc;
                for (const auto &br : cb->branches) cc[br.ch] = br.labels_distinct;
                EXPECT_EQ(counts, cc);
            } else if (P == std::set<std::string>{ "S" }) {
                trie::Leaves right { { f.P + f.Y.substr(0, 10), { "S" } } };
                trie::Leaves left { { "", { "S" } } };
                EXPECT_EQ(right, trie::constrained_leaves(A, kRight, kRadius)) << what;
                EXPECT_EQ(left, trie::constrained_leaves(A, kLeft, kRadius)) << what;
            }
        }
    }
}


// Brief test 2. Tuned runs (branch limits, quorums, a split limit, a beam, caps) are
// prefix-subsets of the exhaustive trie and record a reason for every omission; with
// merging on, the per-label routes are.
TYPED_TEST(TrieOracle, TunedRunsArePrefixSubsetsAndEveryOmissionHasAReason) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    const OracleFixture f(kFixtureSeed);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.sequences, f.labels, mode);
        ASSERT_TRUE(topology_is(*anno, mode, f.branching()));
        const std::string where = "mode " + std::to_string(mode) + " ";
        auto A = run(*anno, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
        ASSERT_EQ(ArmResult::COMPLETE, A.arms[kRight].status);
        ASSERT_EQ(ArmResult::COMPLETE, A.arms[kLeft].status);
        // the claims of A: 7 (leaf, label) pairs on the right plus S stopping inside Y
        // on the shared P stretch; 2 leaves x 4 labels on the left plus S at the boundary
        const size_t claims_right = 8, claims_left = 9;

        auto tuned = [](size_t branches) {
            Strategy st;
            st.max_extension_bp = kRadius;
            st.max_label_branches = branches;
            st.merge_reconverge = false;
            return st;
        };
        struct Case {
            const char *name;
            Strategy st;
            std::string omitted_walk;      // a walk of A the case must drop entirely
            std::string omitted_label;     // ... or a label it must drop from every leaf
        };
        std::vector<Case> cases;
        cases.push_back({ "max_label_branches 0", tuned(0), f.qyt(), "C" });
        cases.push_back({ "max_label_branches 1", tuned(1), f.qyt(), "C" });
        Strategy st = tuned(Strategy::kUnlimited);
        st.min_successor_labels = 2;
        cases.push_back({ "min_successor_labels 2", st, f.qyt(), "" });
        st = tuned(Strategy::kUnlimited);
        st.min_successor_fraction = 0.5;
        cases.push_back({ "min_successor_fraction 0.5", st, f.qyz(), "B" });
        st = tuned(Strategy::kUnlimited);
        st.max_splits_per_path = 1;
        cases.push_back({ "max_splits_per_path 1", st, f.pyt(), "" });
        st = tuned(Strategy::kUnlimited);
        st.on_overflow = Strategy::BEAM;
        st.max_live_paths = 1;
        cases.push_back({ "beam 1", st, f.qyt(), "" });
        st = tuned(Strategy::kUnlimited);
        st.max_live_paths = 1;
        cases.push_back({ "max_live_paths 1", st, f.pyz(), "" });
        st = tuned(Strategy::kUnlimited);
        st.max_steps = 50;
        cases.push_back({ "max_steps 50", st, f.pyz(), "" });
        st = tuned(Strategy::kUnlimited);
        st.max_output_bp = 130;
        cases.push_back({ "max_output_bp 130", st, f.pyt(), "" });

        for (const Case &c : cases) {
            const std::string what = where + c.name;
            const trie::SeedContext ctx { f.X, mode != DeBruijnGraph::BASIC, c.st };
            auto t = run(*anno, f.X, as_list(kAllLabels), c.st);
            trie::SubsetReport right = trie::check_tuned_subset(A, t, kRight, what, ctx);
            trie::SubsetReport left = trie::check_tuned_subset(A, t, kLeft, what, ctx);
            EXPECT_EQ(claims_right, right.present + right.omitted) << what;
            EXPECT_EQ(claims_left, left.present + left.omitted) << what;
            // the case is not vacuous: it does prune
            EXPECT_GT(right.omitted, 0u) << what;
            const trie::Leaves leaves = trie::constrained_leaves(t, kRight, kRadius);
            EXPECT_EQ(0u, leaves.count(c.omitted_walk)) << what << ": " << c.omitted_walk;
            if (!c.omitted_label.empty()) {
                for (const auto &[w, ls] : leaves)
                    EXPECT_EQ(0u, ls.count(c.omitted_label)) << what << ": " << w;
            }
        }

        // merging on: the bubble's walks join in Y and the leaf through the first parent
        // carries the other branch's label with route_bp > 0 — not evidence for the
        // spelled bases before the join, which the route check honours
        st = tuned(Strategy::kUnlimited);
        st.merge_reconverge = true;
        auto merged = run(*anno, f.X, as_list(kAllLabels), st);
        EXPECT_FALSE(trie::is_trie(merged.arms[kRight])) << where << "no join?";
        EXPECT_EQ(2u, merged.arms[kRight].paths.size());
        size_t merged_in = 0;
        for (const auto &p : merged.arms[kRight].paths) {
            for (const auto &e : p.end_labels) merged_in += e.route_bp > 0;
        }
        // the Z leaf carries the other bubble branch's label through the second
        // parent (and the T leaf carries D that way when the first parent is Q)
        EXPECT_GE(merged_in, 1u) << where;
        EXPECT_LE(merged_in, 2u) << where;
        // every label end is checked on its own route: the 5 leaf labels on the right
        // plus S ending inside the joined Y stretch, the 8 on the left plus S at the
        // boundary
        EXPECT_EQ(6u, trie::check_routes_subset(A, merged, kRight, where + "merge"));
        EXPECT_EQ(9u, trie::check_routes_subset(A, merged, kLeft, where + "merge"));
    }
}


/*
 * Checker regressions (review round 2). The checkers of the tuned-run property had
 * accepted three deliberate corruptions. Each test below corrupts a REAL result the way
 * the reviewer did and asserts that the checker's report is NOT empty — the reports are
 * read directly (tuned_subset_report, routes_subset_report), so the checks that fail
 * on the corrupted result do not fail these tests — after asserting that the genuine
 * result passes with an empty report.
 */
namespace {

// a LABEL_END for |label| at |at| inside the segment of |arm| that starts at depth 0
// with base |first|; returns how many segments were hit (one, on the fixture)
size_t insert_label_end(const SeedResult &r, ArmResult &arm, const std::string &label,
                        char first, uint64_t at, EndReason reason) {
    const LabelId id = *trie::label_id(r, label);
    size_t hit = 0;
    for (Segment &seg : arm.segments) {
        if (seg.from_bp != 0 || seg.sequence.empty()
                || trie::outward(arm, seg.sequence)[0] != first)
            continue;
        Event ev;
        ev.type = EventType::LABEL_END;
        ev.label = id;
        ev.at_bp = at;
        ev.reason = reason;
        seg.events.push_back(ev);
        hit++;
    }
    return hit;
}

bool any_mentions(const trie::Problems &problems, const std::string &needle) {
    for (const std::string &p : problems) {
        if (p.find(needle) != std::string::npos)
            return true;
    }
    return false;
}

} // namespace

// Finding 1 (round 2). The Q subtree is deleted from an (exhaustive, hence
// tuned-shaped) result and the bubble's ordinary ambiguity event is kept: it says C
// followed BOTH branches and nothing was dropped. The old branch_recorded() took any
// branch event listing the missing base as the recorded reason for every omission
// through it — B's claims among them, which nothing dropped: 5 present, 3 omitted, no
// failure. Discard evidence is about the label and the successor, and (round 3,
// finding 1) only what the walker states explicitly: the event has no refusal for Q,
// so neither B's claim nor C's two claims through Q have a record — the round-2
// checker still excused C's, inferring "refused" from the deleted child.
TEST(Trie, TunedCheckerRejectsASilentlyDeletedBranch) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    const trie::SeedContext ctx { f.X, false, exhaustive_at(LabelMode::CONSTRAIN, kRadius) };
    auto A = run(*anno, f.X, as_list(kAllLabels), ctx.strategy);
    const trie::SubsetReport genuine = trie::tuned_subset_report(A, A, kRight, "genuine", ctx);
    EXPECT_TRUE(genuine.problems.empty()) << trie::listed(genuine.problems);
    EXPECT_EQ(8u, genuine.present);
    EXPECT_EQ(0u, genuine.omitted);

    SeedResult tuned = A;
    ArmResult &arm = tuned.arms[kRight];
    const size_t root = trie::root_of(arm).id;
    ASSERT_EQ(2u, arm.segments[root].children.size());
    size_t removed = arm.segments.size();
    for (size_t s : arm.segments[root].children) {
        if (trie::outward(arm, arm.segments[s].sequence)[0] == f.Q[0])
            removed = s;
    }
    ASSERT_LT(removed, arm.segments.size());
    auto &children = arm.segments[root].children;
    children.erase(std::remove(children.begin(), children.end(), removed), children.end());
    arm.paths.erase(std::remove_if(arm.paths.begin(), arm.paths.end(), [&](const PathResult &p) {
        return std::find(p.segments.begin(), p.segments.end(), removed) != p.segments.end();
    }), arm.paths.end());
    ASSERT_EQ(2u, arm.paths.size());
    // the bubble's event is untouched: C ambiguous, nothing dropped, Q's base listed
    const BranchEvent *bubble = nullptr;
    for (const BranchEvent &be : arm.branch_events) {
        if (be.segment == root && be.at_bp == 0)
            bubble = &be;
    }
    ASSERT_NE(nullptr, bubble);
    EXPECT_TRUE(bubble->dropped.empty());
    EXPECT_TRUE(bubble->refused.empty());
    EXPECT_EQ((std::set<std::string>{ "C" }), trie::name_set(tuned, bubble->ambiguous));
    EXPECT_NE(bubble->chars.end(), std::find(bubble->chars.begin(), bubble->chars.end(), f.Q[0]));
    // ... and is no evidence that anybody was refused Q, on the corrupted result as on
    // the genuine one
    for (const char *l : { "B", "C" }) {
        EXPECT_FALSE(trie::branch_recorded(tuned, arm, ctx, "", root, 0, f.Q[0], l)) << l;
        EXPECT_FALSE(trie::branch_recorded(A, A.arms[kRight], ctx, "", root, 0, f.Q[0], l)) << l;
    }

    const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight, "deleted Q", ctx);
    EXPECT_EQ(5u, rep.present);
    EXPECT_EQ(3u, rep.omitted);
    ASSERT_FALSE(rep.problems.empty()) << "the checker accepted a silently deleted branch";
    // all three claims through Q are unexplained: (Q·Y·Z, B), (Q·Y·Z, C), (Q·Y·T, C)
    EXPECT_EQ(3u, rep.problems.size()) << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "with B alive and no recorded reason"))
        << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "with C alive and no recorded reason"))
        << trie::listed(rep.problems);
}

// Round 3, finding 1: the same deletion where the deleted branch is carried by the
// ambiguous label ONLY. k = 3, records AAA·C and AAA·G both labelled C, seed AAA: C is
// ambiguous at the boundary and followed on both, so the event is {ambiguous: C,
// dropped: none}. Deleting the G child left that event unchanged, and the round-2
// checker read "C ambiguous, G not followed" as "G refused to C by a quorum": 1
// present, 1 omitted, no failure. A refusal is now something the walker states
// (BranchEvent::refused), never something inferred from a missing child.
TEST(Trie, TunedCheckerRejectsADeletedBranchUnderASharedLabel) {
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAAC", "AAAG" }, { "C", "C" }, DeBruijnGraph::BASIC);
    const trie::SeedContext ctx { "AAA", false, exhaustive_at(LabelMode::CONSTRAIN, 100) };
    auto A = run(*anno, "AAA", { "C" }, ctx.strategy);
    ASSERT_EQ(2u, A.arms[kRight].paths.size());
    EXPECT_TRUE(trie::tuned_subset_report(A, A, kRight, "genuine", ctx).problems.empty());

    SeedResult tuned = A;
    ArmResult &arm = tuned.arms[kRight];
    const size_t root = trie::root_of(arm).id;
    ASSERT_EQ(2u, arm.segments[root].children.size());
    const size_t removed = arm.segments[root].children.back();
    const char removed_ch = arm.segments[removed].sequence[0];
    arm.segments[root].children.pop_back();
    arm.paths.erase(std::remove_if(arm.paths.begin(), arm.paths.end(), [&](const PathResult &p) {
        return std::find(p.segments.begin(), p.segments.end(), removed) != p.segments.end();
    }), arm.paths.end());
    ASSERT_EQ(1u, arm.branch_events.size());
    EXPECT_EQ((std::set<std::string>{ "C" }), trie::name_set(tuned, arm.branch_events[0].ambiguous));
    EXPECT_TRUE(arm.branch_events[0].dropped.empty());
    EXPECT_TRUE(arm.branch_events[0].refused.empty());

    const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight, "deleted branch", ctx);
    EXPECT_EQ(1u, rep.present);
    EXPECT_EQ(1u, rep.omitted);
    ASSERT_EQ(1u, rep.problems.size()) << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, std::string("drops ") + removed_ch + " at 0"))
        << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "with C alive and no recorded reason"))
        << trie::listed(rep.problems);
}

// Round 3, finding 1, second half: a deleted branch excused by a FORGED structural
// block. The round-2 checker accepted a BLOCKED event on its reason being structural;
// an edge_reuse at depth 0 of a clean fixture is no reuse (no (k+1)-mer of the seed
// is the step), and a rejoined_seed there re-enters nothing (no k-mer of the seed is
// the one stepped into). Both are verified on the seed and the walk now and both
// forgeries fail; the walker's genuine blocks on the same fixture still verify.
TEST(Trie, TunedCheckerRejectsAForgedStructuralBlock) {
    auto b = fork_blocks(3);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P, X + Q }, { "A", "B" }, DeBruijnGraph::BASIC);
    const trie::SeedContext ctx { X, false, exhaustive_at(LabelMode::CONSTRAIN, 100) };
    auto A = run(*anno, X, { "A", "B" }, ctx.strategy);
    ASSERT_EQ(2u, A.arms[kRight].paths.size());
    for (EndReason forged : { EndReason::EDGE_REUSE, EndReason::EDGE_REUSE_RC, EndReason::REACHED_SEED }) {
        SeedResult tuned = A;
        ArmResult &arm = tuned.arms[kRight];
        const size_t root = trie::root_of(arm).id;
        const size_t removed = arm.segments[root].children.back();
        Event ev;
        ev.type = EventType::BLOCKED;
        ev.at_bp = 0;
        ev.ch = arm.segments[removed].sequence[0];
        ev.labels = arm.segments[removed].labels_start;
        ev.labels_total = ev.labels.size();
        ev.reason = forged;
        arm.segments[root].events.push_back(ev);
        arm.segments[root].children.pop_back();
        arm.paths.erase(std::remove_if(arm.paths.begin(), arm.paths.end(), [&](const PathResult &p) {
            return std::find(p.segments.begin(), p.segments.end(), removed) != p.segments.end();
        }), arm.paths.end());
        const std::string name = tuned.label_dict.at(ev.labels[0]).name;
        EXPECT_FALSE(trie::branch_recorded(tuned, arm, ctx, "", root, 0, ev.ch, name)) << to_string(forged);
        const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight,
                                                                 std::string("forged ") + to_string(forged), ctx);
        EXPECT_EQ(1u, rep.present) << to_string(forged);
        EXPECT_EQ(1u, rep.omitted) << to_string(forged);
        ASSERT_EQ(1u, rep.problems.size()) << to_string(forged) << trie::listed(rep.problems);
        EXPECT_TRUE(any_mentions(rep.problems, "no recorded reason")) << trie::listed(rep.problems);
    }
    // the verification accepts what is true: on the oracle fixture's tandem repeat the
    // second lap's loop edge IS reused, and the walker's BLOCKED edge_reuse there verifies
    const OracleFixture f(kFixtureSeed);
    auto fixture = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    auto E = run(*fixture, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
    const ArmResult &left = E.arms[kLeft];
    const trie::SeedContext fixture_ctx { f.X, false, exhaustive_at(LabelMode::CONSTRAIN, kRadius) };
    size_t verified = 0;
    for (const Segment &seg : left.segments) {
        for (const Event &ev : seg.events) {
            if (ev.type != EventType::BLOCKED || ev.reason != EndReason::EDGE_REUSE)
                continue;
            const std::string walk = f.two_laps().substr(0, ev.at_bp);
            verified += trie::branch_recorded(E, left, fixture_ctx, walk,
                                              seg.id, ev.at_bp, ev.ch, E.label_dict.at(ev.labels[0]).name);
        }
    }
    EXPECT_EQ(1u, verified) << "the loop edge's genuine reuse should verify exactly once";
}

// Finding 2. A label end invented INSIDE a segment of a tuned result: LABEL_END(S,
// dead_end) one base into the Q branch, where S is not even on the first node. The old
// clause 1 looked at the leaves' end_labels only and clause 2 walks the exhaustive
// claims, so the extra claim escaped both (8 present, 0 omitted, no failure). Every
// claim of the tuned run is checked now: (Q[0], S) is a prefix of no exhaustive walk
// carrying S. The same end one base into the P branch — where S IS alive, going on to
// |P| + 10 — is an invented END: semantic, where the exhaustive run continues the
// label.
TEST(Trie, TunedCheckerRejectsAnInventedInteriorClaim) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    const trie::SeedContext ctx { f.X, false, exhaustive_at(LabelMode::CONSTRAIN, kRadius) };
    auto A = run(*anno, f.X, as_list(kAllLabels), ctx.strategy);

    SeedResult invented = A;
    ASSERT_EQ(1u, insert_label_end(invented, invented.arms[kRight], "S", f.Q[0], 1,
                                   EndReason::DEAD_END));
    ASSERT_TRUE(trie::constrained_claims(invented, kRight, kRadius).at(f.Q.substr(0, 1)).count("S"));
    const trie::SubsetReport rep = trie::tuned_subset_report(A, invented, kRight, "S invented on Q", ctx);
    // the exhaustive claims know nothing of it, so the forward walk still finds all 8
    EXPECT_EQ(8u, rep.present);
    EXPECT_EQ(0u, rep.omitted);
    ASSERT_FALSE(rep.problems.empty()) << "the checker accepted an invented interior claim";
    EXPECT_EQ(1u, rep.problems.size()) << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "under S is not a prefix")) << trie::listed(rep.problems);

    SeedResult early = A;
    ASSERT_EQ(1u, insert_label_end(early, early.arms[kRight], "S", f.P[0], 1, EndReason::DEAD_END));
    const trie::SubsetReport rep2 = trie::tuned_subset_report(A, early, kRight, "S ended early on P", ctx);
    ASSERT_FALSE(rep2.problems.empty()) << "the checker accepted an invented semantic end";
    // clause 1 sees an invented end, clause 2 (following S's exhaustive claim) the
    // same end as a disagreement between the two walkers
    EXPECT_EQ(2u, rep2.problems.size()) << trie::listed(rep2.problems);
    EXPECT_TRUE(any_mentions(rep2.problems, "an INVENTED end")) << trie::listed(rep2.problems);
    EXPECT_TRUE(any_mentions(rep2.problems, "the two walkers DISAGREE")) << trie::listed(rep2.problems);
}

// Finding 2, merging on: the same invented end in the merged run's Q segment. The old
// check_routes_subset() reconstructed routes for the leaves' labels only, so an end
// inside a segment was never looked at (all 5 leaf labels checked, no failure). Every
// label end is checked on its own route now, against the exhaustive CLAIMS.
TEST(Trie, MergedCheckerRejectsAnInventedInteriorClaim) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    auto A = run(*anno, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
    Strategy st;
    st.max_extension_bp = kRadius;
    st.max_label_branches = Strategy::kUnlimited;
    st.merge_reconverge = true;
    auto merged = run(*anno, f.X, as_list(kAllLabels), st);
    ASSERT_FALSE(trie::is_trie(merged.arms[kRight]));
    const trie::RouteReport genuine = trie::routes_subset_report(A, merged, kRight, "genuine");
    EXPECT_TRUE(genuine.problems.empty()) << trie::listed(genuine.problems);
    // the 5 leaf labels plus S, which ends inside the joined Y stretch
    EXPECT_EQ(6u, genuine.checked);

    SeedResult bad = merged;
    ASSERT_EQ(1u, insert_label_end(bad, bad.arms[kRight], "S", f.Q[0], 1, EndReason::DEAD_END));
    const trie::RouteReport rep = trie::routes_subset_report(A, bad, kRight, "S invented on Q");
    EXPECT_EQ(7u, rep.checked);
    ASSERT_FALSE(rep.problems.empty()) << "the merged checker accepted an invented interior claim";
    EXPECT_EQ(1u, rep.problems.size()) << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "the route of S to segment")) << trie::listed(rep.problems);
    EXPECT_TRUE(any_mentions(rep.problems, "INVENTED")) << trie::listed(rep.problems);

    // Round 3, the coverage boundary: the checkers read the label end EVENTS, so a
    // result whose leaves' end_labels were cleared passed while the lists a client
    // reads were empty. The leaf records are cross-checked against the events now, in
    // both checkers.
    SeedResult cleared = merged;
    for (PathResult &p : cleared.arms[kRight].paths) {
        p.end_labels.clear();
    }
    const trie::RouteReport rep2 = trie::routes_subset_report(A, cleared, kRight, "cleared end_labels");
    ASSERT_FALSE(rep2.problems.empty()) << "the merged checker accepted cleared leaf labels";
    EXPECT_TRUE(any_mentions(rep2.problems, "lists 0 end label(s)")) << trie::listed(rep2.problems);
    SeedResult rewritten = A;
    for (PathResult &p : rewritten.arms[kRight].paths) {
        p.end_reasons.fill(0);
        p.end_reasons[static_cast<size_t>(EndReason::RECORD_END)] = p.end_labels.size();
    }
    const trie::SubsetReport rep3 = trie::tuned_subset_report(
            A, rewritten, kRight, "rewritten end_reasons",
            trie::SeedContext{ f.X, false, exhaustive_at(LabelMode::CONSTRAIN, kRadius) });
    ASSERT_FALSE(rep3.problems.empty()) << "the tuned checker accepted rewritten leaf reasons";
    EXPECT_TRUE(any_mentions(rep3.problems, "end_reasons disagree")) << trie::listed(rep3.problems);
}

// Round 3, finding 2: a VALID merged result the round-2 checker rejected. k = 3, the
// one record AAAGTAATAA under C, seed AAA, radius 8 (the graph is node-centric: every
// pair of k-mers overlapping by k - 1 is an edge, so AAA leads to AAG and to AAT). The
// walks G·T·A·A·T·A·A and T·A·A·G·T·A·A both reach node TAA at depth 7 and merge there.
// Each on its own could go on (TAA→AAG for the first, TAA→AAT for the second), but
// under the united edge history (§6.10) the first has used TAA→AAT and the second
// TAA→AAG, and TAA→AAA re-enters the seed: every continuation is barred and C ends
// with rejoined_seed by block precedence. The keep trie continues both walks to the
// radius, so the old rule — a rejoined_seed end must be where the trie ends the label
// — called it an invented end. Termination on a joined route is judged under the
// united history now: the structural block passes on prefix support, and the seed
// re-entry itself, a fact about the node, is confirmed by the trie's own BLOCKED
// event there.
TEST(Trie, MergedCheckerAcceptsAUnitedHistoryTermination) {
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAAGTAATAA" }, { "C" }, DeBruijnGraph::BASIC);
    Strategy st = exhaustive_at(LabelMode::CONSTRAIN, 8);
    st.direction = Strategy::RIGHT;
    auto A = run(*anno, "AAA", { "C" }, st);
    ASSERT_EQ(ArmResult::COMPLETE, A.arms[kRight].status);
    EXPECT_EQ(8u, A.arms[kRight].complete_to_bp);
    EXPECT_STREQ("per_path", A.arms[kRight].completeness_scope);
    st.exhaustive = false;
    st.merge_reconverge = true;
    auto merged = run(*anno, "AAA", { "C" }, st);
    const ArmResult &arm = merged.arms[kRight];
    EXPECT_STREQ("united_history", arm.completeness_scope);
    ASSERT_FALSE(trie::is_trie(arm)) << "no join?";
    // the join at depth 7 ends C with rejoined_seed, where the trie goes on
    size_t joined_ends = 0;
    for (const Segment &seg : arm.segments) {
        if (seg.parents.size() < 2)
            continue;
        for (const Event &ev : seg.events) {
            if (ev.type != EventType::LABEL_END)
                continue;
            joined_ends++;
            EXPECT_EQ(7u, ev.at_bp);
            EXPECT_STREQ("rejoined_seed", to_string(ev.reason));
        }
    }
    EXPECT_GE(joined_ends, 1u);
    EXPECT_TRUE(trie::constrained_claims(A, kRight, 8).count("GTAATAAG"));
    const trie::RouteReport rep = trie::routes_subset_report(A, merged, kRight, "AAAGTAATAA merged");
    EXPECT_TRUE(rep.problems.empty()) << trie::listed(rep.problems);
    EXPECT_GE(rep.checked, 2u);

    // ... while a rejoined_seed the trie does not confirm is still invented: forged in
    // place of a genuine dead_end at a joined leaf of the oracle fixture (the end of Z,
    // where no successor leads into the seed)
    const OracleFixture f(kFixtureSeed);
    auto fixture = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    auto E = run(*fixture, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
    Strategy mst;
    mst.max_extension_bp = kRadius;
    mst.max_label_branches = Strategy::kUnlimited;
    mst.merge_reconverge = true;
    SeedResult forged = run(*fixture, f.X, as_list(kAllLabels), mst);
    ASSERT_TRUE(trie::routes_subset_report(E, forged, kRight, "genuine").problems.empty());
    size_t forged_ends = 0;
    for (PathResult &p : forged.arms[kRight].paths) {
        // the Z leaf, spelled through whichever bubble branch is the first parent
        const std::string w = trie::walk_of(forged.arms[kRight], p);
        if (w != f.pyz() && w != f.qyz())
            continue;
        Segment &last = forged.arms[kRight].segments[p.segments.back()];
        for (Event &ev : last.events) {
            if (ev.type == EventType::LABEL_END && ev.at_bp == p.length_bp) {
                ev.reason = EndReason::REACHED_SEED;
                forged_ends++;
            }
        }
        p.end_reasons.fill(0);
        p.end_reasons[static_cast<size_t>(EndReason::REACHED_SEED)] = forged_ends;
    }
    ASSERT_EQ(3u, forged_ends) << "A, B and C end at the end of Z";
    const trie::RouteReport rep2 = trie::routes_subset_report(E, forged, kRight, "forged rejoined_seed");
    EXPECT_EQ(3u, rep2.problems.size()) << trie::listed(rep2.problems);
    EXPECT_TRUE(any_mentions(rep2.problems, "records no seed-node successor")) << trie::listed(rep2.problems);
}

// Round 3, finding 2, at scale: on dense random k = 3 graphs (three records of 30
// random bases after the seed AAA, one label) every genuine merged run passes the
// route checker against its exhaustive trie. Reconvergences are everywhere at k = 3,
// so the united-history terminations the checker has to accept and the per-path ends
// it still compares exactly both occur in numbers — a wrong rule either way would
// show as false rejections here (the reviewer's valid-merge stress, 250 seeds).
TEST(Trie, MergedCheckerAcceptsGenuineDenseMerges) {
    size_t joined_ends = 0, checked = 0;
    for (uint32_t seed = 1; seed <= 80; ++seed) {
        const std::vector<std::string> seqs { "AAA" + random_seq(30, seed * 3),
                                              "AAA" + random_seq(30, seed * 3 + 1),
                                              "AAA" + random_seq(30, seed * 3 + 2) };
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                3, seqs, { "C", "C", "C" }, DeBruijnGraph::BASIC);
        Strategy st = exhaustive_at(LabelMode::CONSTRAIN, 8);
        st.direction = Strategy::RIGHT;
        st.max_live_paths = 100000;
        st.max_paths = 100000;
        auto A = run(*anno, "AAA", { "C" }, st);
        ASSERT_EQ(ArmResult::COMPLETE, A.arms[kRight].status) << "seed " << seed;
        st.exhaustive = false;
        st.merge_reconverge = true;
        auto merged = run(*anno, "AAA", { "C" }, st);
        const std::string what = "seed " + std::to_string(seed);
        const trie::RouteReport rep = trie::routes_subset_report(A, merged, kRight, what);
        EXPECT_TRUE(rep.problems.empty()) << what << trie::listed(rep.problems);
        checked += rep.checked;
        for (const Segment &seg : merged.arms[kRight].segments) {
            if (seg.parents.size() < 2)
                continue;
            for (const Event &ev : seg.events)
                joined_ends += ev.type == EventType::LABEL_END;
        }
    }
    // the sweep is not vacuous: ends under the united history did occur
    EXPECT_GT(joined_ends, 0u);
    EXPECT_GT(checked, 80u);
}

// Round 3, minor: a REFERENCE whose own claim ended with a cap establishes prefix
// support, not where the label stops. k = 3, AAA·CG and AAA·G under C, seed AAA; the
// reference runs with max_steps 2 and is complete through depth 1, where both of its
// claims are censored (max_steps). Uncapped tuned and merged runs continue C·G to
// depth 2 and end G at depth 1 with dead_end; the old checkers rejected the first for
// outliving the reference and the second for ending at its boundary for another
// reason. A censored reference end is unknown: prefix support only.
TEST(Trie, CensoredReferenceBoundaryIsUnknown) {
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAACG", "AAAG" }, { "C", "C" }, DeBruijnGraph::BASIC);
    Strategy st = exhaustive_at(LabelMode::CONSTRAIN, 100);
    st.direction = Strategy::RIGHT;
    st.max_steps = 2;
    auto A = run(*anno, "AAA", { "C" }, st);
    ASSERT_EQ(ArmResult::TRUNCATED, A.arms[kRight].status);
    ASSERT_EQ(1u, A.arms[kRight].complete_to_bp);
    for (const PathResult &p : A.arms[kRight].paths) {
        EXPECT_EQ(1u, p.length_bp);
        EXPECT_EQ(1u, p.end_reasons[static_cast<size_t>(EndReason::MAX_STEPS)]);
    }
    st.exhaustive = false;
    st.max_steps = 100;
    auto keep = run(*anno, "AAA", { "C" }, st);
    EXPECT_EQ((std::set<std::string>{ "CG", "G" }), leaf_walks(keep.arms[kRight]));
    const trie::SubsetReport tuned = trie::tuned_subset_report(A, keep, kRight, "uncapped keep",
                                                               trie::SeedContext{ "AAA", false, st });
    EXPECT_TRUE(tuned.problems.empty()) << trie::listed(tuned.problems);
    EXPECT_EQ(2u, tuned.present + tuned.omitted);
    st.merge_reconverge = true;
    auto merged = run(*anno, "AAA", { "C" }, st);
    const trie::RouteReport routes = trie::routes_subset_report(A, merged, kRight, "uncapped merged");
    EXPECT_TRUE(routes.problems.empty()) << trie::listed(routes.problems);
    EXPECT_EQ(2u, routes.checked);
}

#if ! _PROTEIN_GRAPH
// Review round 4, finding 1: a FOLLOWED hairpin excused the deletion of its child.
// CANONICAL, k = 3, records AAATT and AAATG under C, seed AAA, exhaustive (keep) with
// hairpins followed: at depth 1 the step AAT → ATT is its own reverse complement, so the
// walker follows it and flags it HAIRPIN(T, "followed", {C}) on the parent. With that
// child and its paths deleted, the round-3 checker verified the step to be a hairpin and
// took the event for the reason the child is missing (1 present, 1 omitted, no problem):
// it checked the geometry, not the policy. A hairpin is discard evidence only where the
// strategy skips hairpins, and a followed one never. The genuine skip still verifies:
// walked with hairpins skipped, the fixture omits the T child for a reason the checker
// accepts — but not under a context whose strategy follows hairpins.
TEST(Trie, TunedCheckerRejectsADeletedFollowedHairpin) {
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAATT", "AAATG" }, { "C", "C" }, DeBruijnGraph::CANONICAL);
    Strategy follow = exhaustive_at(LabelMode::CONSTRAIN, 100);
    follow.direction = Strategy::RIGHT;
    follow.skip_hairpins = false;
    const trie::SeedContext ctx { "AAA", true, follow };
    auto A = run(*anno, "AAA", { "C" }, follow);
    const trie::SubsetReport genuine = trie::tuned_subset_report(A, A, kRight, "genuine", ctx);
    EXPECT_TRUE(genuine.problems.empty()) << trie::listed(genuine.problems);

    SeedResult tuned = A;
    ArmResult &arm = tuned.arms[kRight];
    size_t parent = arm.segments.size(), removed = arm.segments.size();
    char hairpin = 0;
    for (const Segment &s : arm.segments) {
        for (const Event &ev : s.events) {
            if (ev.type != EventType::HAIRPIN || ev.text != "followed")
                continue;
            for (size_t c : s.children) {
                if (trie::outward(arm, arm.segments[c].sequence)[0] == ev.ch) {
                    parent = s.id;
                    removed = c;
                    hairpin = ev.ch;
                }
            }
        }
    }
    ASSERT_LT(removed, arm.segments.size()) << "no followed hairpin with a child";
    EXPECT_EQ('T', hairpin);
    auto &children = arm.segments[parent].children;
    children.erase(std::remove(children.begin(), children.end(), removed), children.end());
    arm.paths.erase(std::remove_if(arm.paths.begin(), arm.paths.end(), [&](const PathResult &p) {
        return std::find(p.segments.begin(), p.segments.end(), removed) != p.segments.end();
    }), arm.paths.end());
    // the split is at the parent's end: the walk to it is the parent's route
    const uint64_t at = arm.segments[parent].from_bp + arm.segments[parent].length_bp;
    const std::string walk = trie::label_route(arm, parent, 0);
    ASSERT_EQ(at, walk.size());
    EXPECT_FALSE(trie::branch_recorded(tuned, arm, ctx, walk, parent, at, hairpin, "C"));
    const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight, "deleted followed hairpin", ctx);
    EXPECT_GT(rep.omitted, 0u);
    ASSERT_FALSE(rep.problems.empty()) << "the checker accepted a deleted followed hairpin";
    EXPECT_TRUE(any_mentions(rep.problems, "with C alive and no recorded reason"))
        << trie::listed(rep.problems);

    // hairpins skipped: the T child is not walked, the HAIRPIN event says why
    Strategy skip = follow;
    skip.skip_hairpins = true;
    auto skipped = run(*anno, "AAA", { "C" }, skip);
    const trie::SeedContext skip_ctx { "AAA", true, skip };
    const trie::SubsetReport accepted = trie::tuned_subset_report(A, skipped, kRight, "skipped hairpin", skip_ctx);
    EXPECT_TRUE(accepted.problems.empty()) << trie::listed(accepted.problems);
    EXPECT_GT(accepted.omitted, 0u) << "the skip case is vacuous";
    size_t skipped_events = 0;
    for (const Segment &s : skipped.arms[kRight].segments) {
        for (const Event &ev : s.events) {
            if (ev.type != EventType::HAIRPIN)
                continue;
            EXPECT_TRUE(ev.text.empty()) << "a followed hairpin under skip_hairpins";
            const std::string to = trie::label_route(skipped.arms[kRight], s.id, 0).substr(0, ev.at_bp);
            EXPECT_TRUE(trie::branch_recorded(skipped, skipped.arms[kRight], skip_ctx, to, s.id,
                                              ev.at_bp, ev.ch, "C"));
            // the same event is no evidence where the strategy follows hairpins
            EXPECT_FALSE(trie::branch_recorded(skipped, skipped.arms[kRight], ctx, to, s.id,
                                               ev.at_bp, ev.ch, "C"));
            skipped_events++;
        }
    }
    EXPECT_GE(skipped_events, 1u);
    const trie::SubsetReport mismatched = trie::tuned_subset_report(A, skipped, kRight, "skip read as follow", ctx);
    EXPECT_FALSE(mismatched.problems.empty()) << "a skipped hairpin excused under a follow strategy";
}
#endif

// Review round 4, the trust boundary: a refusal was taken on its cause string. k = 3,
// AAA·C and AAA·G under C, seed AAA, the G child deleted from the exhaustive trie: with
// no refusal the omission is unexplained, but adding {G, "branch", C} — under unlimited
// branching — or {G, "not_a_cause", C} made the checker pass. A refusal is now checked
// against the strategy the tuned run was made with (refusal_problem): the cause must
// name an active knob whose condition holds there. Every cause forged onto the exhaustive
// run is rejected, the forgery under a finite branch limit too (C goes on along the C
// child, which a branch-limit exclusion forbids). Genuine refusals of each cause pass,
// including a quorum stop counted after a branch-limit exclusion on the same successor.
TEST(Trie, TunedCheckerRejectsARefusalItsStrategyDoesNotMake) {
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAAC", "AAAG" }, { "C", "C" }, DeBruijnGraph::BASIC);
    const Strategy full = exhaustive_at(LabelMode::CONSTRAIN, 100);
    const trie::SeedContext ctx { "AAA", false, full };
    auto A = run(*anno, "AAA", { "C" }, full);
    SeedResult tuned = A;
    ArmResult &arm = tuned.arms[kRight];
    const size_t root = trie::root_of(arm).id;
    const size_t removed = arm.segments[root].children.back();
    const char ch = arm.segments[removed].sequence[0];
    arm.segments[root].children.pop_back();
    arm.paths.erase(std::remove_if(arm.paths.begin(), arm.paths.end(), [&](const PathResult &p) {
        return std::find(p.segments.begin(), p.segments.end(), removed) != p.segments.end();
    }), arm.paths.end());
    ASSERT_EQ(1u, arm.branch_events.size());
    ASSERT_EQ(1u, trie::tuned_subset_report(A, tuned, kRight, "no refusal", ctx).problems.size());
    const LabelId c = *trie::label_id(tuned, "C");
    for (const char *cause : { "branch", "not_a_cause", "minority", "below_min_labels",
                               "split_limit", "loss_budget" }) {
        arm.branch_events[0].refused = { { ch, cause, { c } } };
        EXPECT_FALSE(trie::branch_recorded(tuned, arm, ctx, "", root, 0, ch, "C")) << cause;
        const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight, cause, ctx);
        EXPECT_EQ(2u, rep.problems.size()) << cause << trie::listed(rep.problems);
        EXPECT_TRUE(any_mentions(rep.problems, "UNSUPPORTED")) << cause << trie::listed(rep.problems);
        EXPECT_TRUE(any_mentions(rep.problems, "no recorded reason")) << cause << trie::listed(rep.problems);
    }
    // a finite allowance supports the cause, not this result: C was not excluded
    Strategy limited = full;
    limited.exhaustive = false;
    limited.max_label_branches = 0;
    arm.branch_events[0].refused = { { ch, "branch", { c } } };
    const trie::SubsetReport rep = trie::tuned_subset_report(A, tuned, kRight, "branch, limit 0",
                                                             trie::SeedContext{ "AAA", false, limited });
    EXPECT_TRUE(any_mentions(rep.problems, "does not end there with branch")) << trie::listed(rep.problems);

    // genuine refusals of every cause the walker emits under forbid pass
    struct Case {
        const char *cause;
        std::vector<std::string> seqs, labels, permitted;
        Strategy st;
    };
    std::vector<Case> cases;
    Strategy st = full;
    st.exhaustive = false;
    st.max_label_branches = 0;
    cases.push_back({ "branch", { "AAAC", "AAAG" }, { "C", "C" }, { "C" }, st });
    st = full;
    st.exhaustive = false;
    st.min_live_labels = 2;
    cases.push_back({ "below_min_labels", { "AAAC", "AAAC", "AAAG" }, { "C", "D", "C" }, { "C", "D" }, st });
    st = full;
    st.exhaustive = false;
    st.max_splits_per_path = 0;
    cases.push_back({ "split_limit", { "AAAC", "AAAG" }, { "C", "D" }, { "C", "D" }, st });
    // A ambiguous and excluded; the C successor keeps B alone, below the quorum of 2
    // (labels_per_successor says 2: it is counted before the exclusion)
    st = full;
    st.exhaustive = false;
    st.max_label_branches = 0;
    st.min_successor_labels = 2;
    cases.push_back({ "minority", { "AAAC", "AAAC", "AAAG" }, { "A", "B", "A" }, { "A", "B" }, st });
    for (const Case &k : cases) {
        auto g = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, k.seqs, k.labels, DeBruijnGraph::BASIC);
        Strategy reference = full;
        reference.direction = Strategy::RIGHT;
        Strategy tuned_st = k.st;
        tuned_st.direction = Strategy::RIGHT;
        auto E = run(*g, "AAA", k.permitted, reference);
        auto t = run(*g, "AAA", k.permitted, tuned_st);
        size_t with_cause = 0;
        for (const BranchEvent &be : t.arms[kRight].branch_events) {
            for (const auto &rf : be.refused)
                with_cause += std::string(rf.cause) == k.cause;
        }
        EXPECT_GT(with_cause, 0u) << k.cause << ": the case refuses nothing";
        const trie::SubsetReport r = trie::tuned_subset_report(E, t, kRight, k.cause,
                                                               trie::SeedContext{ "AAA", false, tuned_st });
        EXPECT_TRUE(r.problems.empty()) << k.cause << trie::listed(r.problems);
        EXPECT_GT(r.omitted, 0u) << k.cause;
    }
}

// Review round 4: the stricter checkers must not reject genuine output. On dense random
// graphs (k = 3, five records after the seed under three labels, every mode, hairpins
// followed and skipped) every tuned run — each branch, quorum and split knob, all of
// them at once, a step cap — passes against its exhaustive trie: every refusal it states
// is one its strategy makes, and a skipped hairpin excuses its omission (the reviewer's
// dense sweep for finding 1, widened to every refusal cause).
// The refusals ARE the evidence, and on a graph this dense the default 100 branch events
// do not hold them all. The cut is stated (branch_events_complete_to_bp, a level
// boundary), so the sweep runs with the DEFAULT cap: below the boundary every omission
// must be explained exactly as without a cap, at or beyond it an unexplained one counts
// as unexplained_capped; with "unlimited" nothing is cut and nothing is unexplained.
namespace {

struct DenseSweep {
    std::map<std::string, size_t> causes;
    size_t cells = 0;
    size_t cut_arms = 0;       // arms whose branch events were cut by the cap
    size_t capped = 0;         // omissions counted unexplained_capped, over all cells
};

DenseSweep dense_tuned_sweep(size_t max_branch_events) {
    DenseSweep out;
    for (auto mode : all_modes()) {
        for (bool skip : { false, true }) {
            for (uint32_t s = 1; s <= 6; ++s) {
                std::vector<std::string> seqs, labels;
                for (uint32_t i = 0; i < 5; ++i) {
                    seqs.push_back("AAA" + random_seq(12, s * 7 + i));
                    labels.push_back(std::string(1, "CDECD"[i]));
                }
                auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, seqs, labels, mode);
                Strategy full = exhaustive_at(LabelMode::CONSTRAIN, 8);
                full.skip_hairpins = skip;
                full.max_live_paths = 100000;
                full.max_paths = 100000;
                // the checker reads the exhaustive run's claims, never its events
                auto A = run(*anno, "AAA", { "C", "D", "E" }, full);
                for (int knob = 0; knob < 8; ++knob) {
                    Strategy st = full;
                    st.exhaustive = false;
                    st.max_branch_events = max_branch_events;
                    switch (knob) {
                        case 0: st.max_label_branches = 0; break;
                        case 1: st.max_label_branches = 1; break;
                        case 2: st.min_successor_labels = 2; break;
                        case 3: st.min_successor_fraction = 0.5; break;
                        case 4: st.min_live_labels = 2; break;
                        case 5: st.max_splits_per_path = 1; break;
                        case 6:
                            st.max_label_branches = 0;
                            st.min_successor_labels = 2;
                            st.min_successor_fraction = 0.4;
                            st.max_splits_per_path = 2;
                            break;
                        case 7: st.max_steps = 10; break;
                    }
                    auto t = run(*anno, "AAA", { "C", "D", "E" }, st);
                    // the walker's guarantee, against the same run with every event kept:
                    // the events kept are the first ones, and every event below the
                    // boundary is among them
                    std::optional<SeedResult> all;
                    if (max_branch_events != Strategy::kUnlimited) {
                        Strategy st_all = st;
                        st_all.max_branch_events = Strategy::kUnlimited;
                        all = run(*anno, "AAA", { "C", "D", "E" }, st_all);
                    }
                    const trie::SeedContext ctx { "AAA", mode != DeBruijnGraph::BASIC, st };
                    for (size_t side : { kLeft, kRight }) {
                        const ArmResult &ta = t.arms[side];
                        if (all) {
                            const auto &every = all->arms[side].branch_events;
                            EXPECT_EQ(every.size(), ta.branch_events_total);
                            for (size_t i = 0; i < every.size(); ++i) {
                                if (i < ta.branch_events.size()) {
                                    const BranchEvent &x = ta.branch_events[i], &y = every[i];
                                    EXPECT_TRUE(x.at_bp == y.at_bp && x.segment == y.segment
                                                && x.chars == y.chars && x.ambiguous == y.ambiguous
                                                && x.dropped == y.dropped
                                                && x.refused.size() == y.refused.size())
                                        << "event " << i << " differs from the uncapped run's";
                                } else {
                                    EXPECT_GE(every[i].at_bp, ta.branch_events_complete_to_bp)
                                        << "an event below the boundary was not kept";
                                }
                            }
                        }
                        for (const BranchEvent &be : ta.branch_events) {
                            for (const auto &rf : be.refused) out.causes[rf.cause]++;
                        }
                        const std::string what = "mode " + std::to_string(mode) + " skip "
                            + std::to_string(skip) + " random " + std::to_string(s) + " knob "
                            + std::to_string(knob) + " arm " + std::to_string(side);
                        const trie::SubsetReport rep = trie::tuned_subset_report(A, t, side, what, ctx);
                        if (!rep.problems.empty()) {
                            ADD_FAILURE() << trie::listed(rep.problems);
                            return out;
                        }
                        const uint64_t boundary = ta.branch_events_complete_to_bp;
                        if (boundary == std::numeric_limits<uint64_t>::max()) {
                            // nothing cut: every omission carries its reason
                            EXPECT_EQ(ta.branch_events_total, ta.branch_events.size()) << what;
                            EXPECT_EQ(0u, rep.unexplained_capped) << what;
                        } else {
                            out.cut_arms++;
                            EXPECT_LT(ta.branch_events.size(), ta.branch_events_total) << what;
                            EXPECT_GE(rep.min_capped_bp, boundary) << what;
                        }
                        out.capped += rep.unexplained_capped;
                        out.cells++;
                    }
                }
            }
        }
    }
    return out;
}

} // namespace

TEST(Trie, CheckersAcceptGenuineTunedRunsOnDenseGraphs) {
    const DenseSweep sweep = dense_tuned_sweep(Strategy().max_branch_events);
    // not vacuous: every cause the walker emits under forbid occurred ...
    for (const char *cause : { "branch", "minority", "below_min_labels", "split_limit" })
        EXPECT_GT(sweep.causes.count(cause) ? sweep.causes.at(cause) : 0, 0u) << cause;
    EXPECT_EQ(all_modes().size() * 2 * 6 * 8 * 2, sweep.cells);
    // ... the default cap does cut evidence here, and omissions resting on the cut
    // events exist and are accepted only at or beyond the boundary (checked per cell)
    EXPECT_GT(sweep.cut_arms, 0u);
    EXPECT_GT(sweep.capped, 0u);
}

// The same sweep with every branch event kept: no boundary, and every omission of every
// tuned run is explained by a recorded reason.
TEST(Trie, CheckersAcceptGenuineTunedRunsOnDenseGraphsWithAllEvents) {
    const DenseSweep sweep = dense_tuned_sweep(Strategy::kUnlimited);
    for (const char *cause : { "branch", "minority", "below_min_labels", "split_limit" })
        EXPECT_GT(sweep.causes.count(cause) ? sweep.causes.at(cause) : 0, 0u) << cause;
    EXPECT_EQ(all_modes().size() * 2 * 6 * 8 * 2, sweep.cells);
    EXPECT_EQ(0u, sweep.cut_arms);
    EXPECT_EQ(0u, sweep.capped);
}

// Finding 3 (the reference model, support: trace). Two occurrences of the seed under
// one label spell the same walk and stop there for different reasons: the earlier one
// would re-enter the seed (rejoined_seed, a structural block), the later one reaches
// its record's end (record_end). The walker keeps one lineage with the coordinates of
// both and the blocked successor sets the reason; the model used to let the LAST
// occurrence choose, so with two records its answer flipped with their order. k = 3
// keeps the records legible: AAA·C·AAA·C·AA holds the seed AAA at 0 and at 4, both
// continuing C·A·A on the right (and, mirrored, on the left).
TEST(Trie, TraceReferenceAggregatesTheOccurrencesOfTheSeed) {
    const size_t k = 3;
    const uint64_t radius = 100;
    const std::string seed = "AAA";
    struct Case {
        const char *name;
        std::vector<std::string> seqs;
        std::vector<uint64_t> coords;
    };
    const std::vector<Case> cases {
        { "one record", { "AAACAAACAA" }, { 0 } },
        { "blocked occurrence first", { "AAACAAA", "AAACAA" }, { 0, 100 } },
        { "record end first", { "AAACAA", "AAACAAA" }, { 0, 100 } },
    };
    for (const Case &c : cases) {
        const std::vector<std::string> labels(c.seqs.size(), "F");
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                k, c.seqs, labels, DeBruijnGraph::BASIC, true, c.coords);
        const trie::StringIndex ix(k, DeBruijnGraph::BASIC, c.seqs, labels, c.coords);
        Strategy st = exhaustive_at(LabelMode::CONSTRAIN, radius);
        st.support = Support::TRACE;
        const SeedResult actual = run(*anno, seed, { "F" }, st);
        EXPECT_STREQ("tuples", actual.access_path) << c.name;
        for (size_t a : { kLeft, kRight }) {
            const Arm arm = actual.arms[a].arm;
            const std::string what = std::string(c.name) + " arm " + to_string(arm);
            EXPECT_EQ(ArmResult::COMPLETE, actual.arms[a].status) << what;
            const trie::RefClaims ref = trie::trace_claims(ix, { "F" }, seed, arm, radius);
            // literally: one walk C·A·A on either side, ended by the seed re-entry
            EXPECT_EQ((trie::Leaves{ { "CAA", { "F" } } }), ref.leaves) << what;
            ASSERT_EQ(1u, ref.reasons.count({ "CAA", "F" })) << what;
            EXPECT_STREQ(to_string(EndReason::REACHED_SEED), to_string(ref.reasons.at({ "CAA", "F" })))
                << what;
            // and the walker agrees, walk by walk and reason by reason
            trie::expect_equal(ref.leaves, trie::constrained_claims(actual, a, radius), what);
            size_t compared = 0;
            for (const auto &[w, ls] : trie::label_claims(actual, a, radius)) {
                for (const auto &[l, reason] : ls) {
                    auto it = ref.reasons.find({ w, l });
                    ASSERT_NE(ref.reasons.end(), it) << what << ": " << w << " " << l;
                    ASSERT_TRUE(reason.has_value()) << what;
                    EXPECT_STREQ(to_string(it->second), to_string(*reason)) << what << " on " << w;
                    compared++;
                }
            }
            EXPECT_EQ(1u, compared) << what;
        }
    }
    // the two record orders give the same model, walk by walk
    const std::vector<std::string> labels { "F", "F" };
    const trie::StringIndex fwd(k, DeBruijnGraph::BASIC, cases[1].seqs, labels, cases[1].coords);
    const trie::StringIndex rev(k, DeBruijnGraph::BASIC, cases[2].seqs, labels, cases[2].coords);
    for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
        const trie::RefTrie a = trie::trace_trie(fwd, "F", seed, arm, radius);
        const trie::RefTrie b = trie::trace_trie(rev, "F", seed, arm, radius);
        ASSERT_EQ(a.size(), b.size()) << to_string(arm);
        for (const auto &[w, leaf] : a) {
            ASSERT_EQ(1u, b.count(w)) << to_string(arm) << " " << w;
            EXPECT_STREQ(to_string(leaf.reason), to_string(b.at(w).reason)) << to_string(arm) << " " << w;
            EXPECT_EQ(leaf.labels, b.at(w).labels) << to_string(arm) << " " << w;
        }
    }
}


// Brief test 3b. A size cap trips between two heads of one level: complete_to_bp is
// the last complete depth, every walk up to it is present, and at complete_to_bp + 1
// a walk IS missing — the partial level is excluded from the guarantee, in both modes.
TEST(Trie, PartialLevelIsExcludedFromTheGuarantee) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    ASSERT_TRUE(topology_is(*anno, DeBruijnGraph::BASIC, f.branching()));
    for (LabelMode mode : { LabelMode::ANNOTATE, LabelMode::CONSTRAIN }) {
        const std::vector<std::string> labels = mode == LabelMode::CONSTRAIN
            ? as_list(kAllLabels) : std::vector<std::string>{};
        Strategy full = exhaustive_at(mode, kRadius);
        full.direction = Strategy::RIGHT;
        auto ref = run(*anno, f.X, labels, full);
        const ArmResult &rarm = ref.arms[kRight];
        ASSERT_EQ(ArmResult::COMPLETE, rarm.status);
        ASSERT_EQ(4u, rarm.paths.size());

        // right arm alone: two heads per level up to the tip fork at 60, four beyond
        struct Case { const char *name; EndReason reason; uint64_t expect; Strategy st; };
        std::vector<Case> cases;
        Strategy st = full;
        st.max_steps = 125;          // 120 to depth 60, 124 at 61, trips on the 2nd head of 62
        cases.push_back({ "max_steps", EndReason::MAX_STEPS, 61, st });
        st = full;
        st.max_output_bp = 125;
        cases.push_back({ "max_output_bp", EndReason::MAX_OUTPUT, 61, st });
        st = full;
        st.max_live_paths = 3;       // the 2nd tip fork would make 4 live heads
        cases.push_back({ "max_live_paths", EndReason::MAX_LIVE_PATHS, 60, st });
        st = full;
        st.max_paths = 3;
        cases.push_back({ "max_paths", EndReason::MAX_PATHS, 60, st });
        st = full;
        st.time_budget_ms = 0;       // stops after the first level
        cases.push_back({ "time_budget_ms", EndReason::TIME_BUDGET, 1, st });

        for (const Case &c : cases) {
            const std::string what = std::string(to_string(mode)) + " " + c.name;
            auto r = run(*anno, f.X, labels, c.st);
            const ArmResult &arm = r.arms[kRight];
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status) << what;
            ASSERT_TRUE(arm.cap_trigger.has_value()) << what;
            EXPECT_EQ(c.reason, arm.cap_trigger->reason) << what;
            EXPECT_LT(arm.complete_to_bp, kRadius) << what;
            EXPECT_EQ(c.expect, arm.complete_to_bp) << what;
            const uint64_t d = arm.complete_to_bp;
            // every walk of every length up to d is present ...
            for (uint64_t n = 0; n <= d; ++n) {
                EXPECT_EQ(walks_at(rarm, n), walks_at(arm, n)) << what << " at " << n;
            }
            // ... and at d + 1 some walk is missing: a strict subset
            const auto have = walks_at(arm, d + 1), want = walks_at(rarm, d + 1);
            EXPECT_TRUE(std::includes(want.begin(), want.end(), have.begin(), have.end())) << what;
            EXPECT_LT(have.size(), want.size()) << what;
            // the cut walks are honestly ended with the cap, inside the region nothing is
            for (const auto &p : arm.paths) {
                ASSERT_TRUE(p.path_reason.has_value()) << what;
                if (p.length_bp < d) {
                    EXPECT_FALSE(is_resource_stop(*p.path_reason)) << what;
                } else if (*p.path_reason != EndReason::DEAD_END) {
                    EXPECT_EQ(c.reason, *p.path_reason) << what;
                    EXPECT_TRUE(p.continuation.has_value()) << what;
                }
            }
            if (mode == LabelMode::CONSTRAIN) {
                // cut to the boundary, the capped run has the complete run's leaves
                // with the same labels alive
                EXPECT_EQ(trie::constrained_leaves(ref, kRight, d),
                          trie::constrained_leaves(r, kRight, d)) << what;
            }
        }
        // ... and the structural oracle cut by a cap still verifies the complete
        // constrained run up to its boundary (the "usable on hard loci" clause)
        if (mode == LabelMode::ANNOTATE) {
            Strategy capped = full;
            capped.max_output_bp = 125;
            auto T = run(*anno, f.X, {}, capped);
            EXPECT_EQ(61u, T.arms[kRight].complete_to_bp);
            auto A = run(*anno, f.X, as_list(kAllLabels), [&] {
                Strategy s = exhaustive_at(LabelMode::CONSTRAIN, kRadius);
                s.direction = Strategy::RIGHT;
                return s;
            }());
            auto v = trie::verify(T, A, kRight, kAllLabels);
            EXPECT_EQ(61u, v.depth);
            EXPECT_FALSE(v.cut);
            trie::expect_equal(v.expected, v.actual, "structural oracle cut at 61");
            // four walks of 61 bases (the tip fork at 60 has split them) plus S's own
            // claim, which ends at |P| + 10 = 25 inside the P stretch
            EXPECT_EQ(5u, v.expected.size());
            for (const auto &[w, ls] : v.expected) {
                if (ls == std::set<std::string>{ "S" }) {
                    EXPECT_EQ(f.P + f.Y.substr(0, 10), w);
                } else {
                    EXPECT_EQ(61u, w.size()) << w;
                }
            }
        }
    }
}


// Brief test 4b. A cut recorded list is reported, and the report matters: an oracle
// built from cut lists disagrees with a correct walker. Caps 1 and 2 cut C, which
// reaches every leaf; cap 4 cuts only S, which reaches no leaf, so a comparison of
// LEAVES would agree there by luck — the claim comparison does not (S's claim is
// missing from E), and nodes_labels_truncated tells a reader why.
TEST(Trie, CutRecordedListsAreReportedAndBreakTheOracle) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    auto A = run(*anno, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
    for (size_t cap : { 1u, 2u, 4u }) {
        Strategy st = exhaustive_at(LabelMode::ANNOTATE, kRadius);
        st.max_labels_per_node = cap;
        auto T = run(*anno, f.X, {}, st);
        const ArmResult &arm = T.arms[kRight];
        EXPECT_EQ(ArmResult::COMPLETE, arm.status);
        EXPECT_EQ(5u, arm.max_labels_at_node);
        EXPECT_GT(arm.nodes_labels_truncated, 0u) << "cap " << cap;
        EXPECT_LE(arm.segments[0].labels_start.size(), cap);
        const Split *bubble = split_at(arm, 0);
        ASSERT_NE(nullptr, bubble);
        EXPECT_EQ(5u, bubble->labels_before);
        for (const auto &br : bubble->branches) {
            EXPECT_LE(br.labels.size(), cap);
            EXPECT_EQ(br.ch == f.P[0] ? 4u : 2u, br.labels_distinct);
        }
        auto v = trie::verify(T, A, kRight, kAllLabels);
        EXPECT_TRUE(v.cut) << "cap " << cap;
        // the kept labels are the first |cap| by column: {A}, {A, B}, {A, B, C, D}; every
        // one of these caps makes the claim comparison fail ...
        EXPECT_NE(v.expected, v.actual)
            << "cap " << cap << ": oracle" << trie::describe(v.expected) << "\nwalker"
            << trie::describe(v.actual);
        // ... while the leaf view (the claims on leaf walks) agrees at cap 4 by luck
        const auto walks = leaf_walks(T.arms[kRight]);
        auto leaves_of = [&](const trie::Leaves &claims) {
            trie::Leaves out;
            for (const auto &[w, ls] : claims) {
                if (walks.count(w))
                    out.insert({ w, ls });
            }
            return out;
        };
        EXPECT_EQ(cap <= 2, leaves_of(v.expected) != leaves_of(v.actual)) << "cap " << cap;
    }
    // a cap that fits: nothing cut, the oracle holds
    Strategy st = exhaustive_at(LabelMode::ANNOTATE, kRadius);
    st.max_labels_per_node = 5;
    auto T = run(*anno, f.X, {}, st);
    EXPECT_EQ(0u, T.arms[kRight].nodes_labels_truncated);
    auto v = trie::verify(T, A, kRight, kAllLabels);
    EXPECT_FALSE(v.cut);
    trie::expect_equal(v.expected, v.actual, "cap 5");
}


// The contract is label-granular, not leaf-granular. S's maximal label-consistent
// walk ends inside the P stretch that A, C and D continue, so no LEAF lists S: a walker
// that ended S two nodes early (or late) would pass a comparison of leaves. The claim
// (P·Y[:10], S) is in E and in claims(A), and a shifted end is caught.
TEST(Trie, OracleComparesPerLabelClaimsNotLeaves) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    auto T = run(*anno, f.X, {}, exhaustive_at(LabelMode::ANNOTATE, kRadius));
    auto A = run(*anno, f.X, as_list(kAllLabels), exhaustive_at(LabelMode::CONSTRAIN, kRadius));
    const std::string s_walk = f.P + f.Y.substr(0, 10);
    auto v = trie::verify(T, A, kRight, kAllLabels);
    ASSERT_FALSE(v.cut);
    trie::expect_equal(v.expected, v.actual, "right");
    // S's claim is a proper prefix of other labels' claims and is kept as its own
    ASSERT_EQ(1u, v.expected.count(s_walk));
    EXPECT_EQ((std::set<std::string>{ "S" }), v.expected.at(s_walk));
    for (const auto &[w, ls] : v.expected) {
        if (w != s_walk)
            EXPECT_EQ(0u, ls.count("S")) << w;
    }
    // no leaf carries S: a leaf comparison alone says nothing about where S ends
    for (const auto &[w, ls] : trie::constrained_leaves(A, kRight, kRadius)) {
        EXPECT_EQ(0u, ls.count("S")) << w;
    }

    // a walker that ends S two nodes early: the same leaves, a different claim
    SeedResult wrong = A;
    size_t moved = 0;
    for (Segment &seg : wrong.arms[kRight].segments) {
        for (Event &ev : seg.events) {
            if (ev.type == EventType::LABEL_END && wrong.label_dict[ev.label].name == "S") {
                ASSERT_EQ(s_walk.size(), ev.at_bp);
                ev.at_bp -= 2;
                moved++;
            }
        }
    }
    ASSERT_EQ(1u, moved);
    EXPECT_EQ(trie::constrained_leaves(A, kRight, kRadius),
              trie::constrained_leaves(wrong, kRight, kRadius));
    const trie::Leaves wrong_claims = trie::constrained_claims(wrong, kRight, kRadius);
    EXPECT_NE(v.expected, wrong_claims);
    EXPECT_EQ(1u, wrong_claims.count(s_walk.substr(0, s_walk.size() - 2)));
    EXPECT_EQ(0u, wrong_claims.count(s_walk));
}


// Annotate mode: the labels of a successor that is NOT followed appear only on its
// BLOCKED / HAIRPIN event, and the root's boundary list is in no run. Both carry the
// true count, a cut of either is counted in nodes_labels_truncated, and the counter
// is exactly the number of cut lists in the output.
TEST(Trie, CutEventListsAndTheRootTotalAreReported) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    for (size_t cap : { 2u, 5u }) {
        Strategy st = exhaustive_at(LabelMode::ANNOTATE, kRadius);
        st.max_labels_per_node = cap;
        auto r = run(*anno, f.X, {}, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = r.arms[a];
            const std::string what = "cap " + std::to_string(cap) + " arm " + std::to_string(a);
            // the root: all five labels carry X
            const Segment &root = trie::root_of(arm);
            EXPECT_EQ(5u, root.labels_start_total) << what;
            EXPECT_EQ(std::min<size_t>(5, cap), root.labels_start.size()) << what;
            // every child segment's entry total agrees with its first run's (the
            // root's first run is the node AFTER the boundary, which is in no run)
            for (const Segment &seg : arm.segments) {
                if (seg.parents.empty() || seg.label_sets.empty())
                    continue;
                EXPECT_EQ(seg.label_sets.front().labels_total, seg.labels_start_total) << what;
                EXPECT_EQ(seg.label_sets.front().labels, seg.labels_start) << what;
            }
            // the counter is the number of cut lists: the root, every node entered
            // (a run covers to_bp - from_bp of them), every not followed successor
            size_t cut = root.labels_start_total > root.labels_start.size();
            size_t blocked = 0;
            for (const Segment &seg : arm.segments) {
                for (const LabelSetRun &run : seg.label_sets) {
                    if (run.truncated())
                        cut += run.to_bp - run.from_bp;
                }
                for (const Event &ev : seg.events) {
                    if (ev.type == EventType::BLOCKED || ev.type == EventType::HAIRPIN) {
                        blocked++;
                        EXPECT_EQ(std::min(ev.labels_total, cap), ev.labels.size()) << what;
                        cut += ev.truncated();
                    }
                }
            }
            EXPECT_EQ(cut, arm.nodes_labels_truncated) << what;
            EXPECT_EQ(cap < 5, arm.nodes_labels_truncated > 0) << what;
            if (a == kLeft) {
                // the second lap's re-entry is blocked (edge reuse); the labels present
                // on that successor are the four carrying the repeat, not S
                ASSERT_EQ(1u, blocked) << what;
                for (const Segment &seg : arm.segments) {
                    for (const Event &ev : seg.events) {
                        if (ev.type != EventType::BLOCKED)
                            continue;
                        EXPECT_EQ(EndReason::EDGE_REUSE, ev.reason);
                        EXPECT_EQ(4u, ev.labels_total) << what;
                        EXPECT_EQ(std::min<size_t>(4, cap), ev.labels.size()) << what;
                        EXPECT_EQ(cap < 4, ev.truncated()) << what;
                        if (cap >= 4) {
                            EXPECT_EQ((std::vector<std::string>{ "A", "B", "C", "D" }),
                                      names(r, ev.labels));
                        }
                    }
                }
            } else {
                EXPECT_EQ(0u, blocked) << what;
            }
        }
    }
}


// Under `exhaustive` a derived permitted set is never cut: more carriers than
// max_seed_labels is a per-seed refusal naming the knob, not a trie over a subset.
// Without the preset the cap applies and is reported, as before.
TEST(Trie, ExhaustiveRefusesACutDerivedSet) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    Strategy st = exhaustive_at(LabelMode::CONSTRAIN, kRadius);
    st.max_seed_labels = 3;             // five labels carry X
    try {
        run(*anno, f.X, {}, st);
        FAIL() << "a cut derived set was accepted under exhaustive";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("max_seed_labels")) << e.what();
        EXPECT_NE(std::string::npos, std::string(e.what()).find("5 labels")) << e.what();
        // the cause as a code, with the cap and what it met (the response's derivation
        // limitation is built from these, not from the message)
        EXPECT_EQ(SeedDerivationError::OVER_SEED_LABEL_CAP, e.cause());
        EXPECT_EQ(3.0, e.limit());
        EXPECT_EQ(5.0, e.observed());
    }
    st.max_seed_labels = 5;
    auto full = run(*anno, f.X, {}, st);
    EXPECT_TRUE(full.labels_from_seed);
    EXPECT_EQ(5u, full.num_seed_labels);
    EXPECT_EQ(0u, full.labels_dropped);
    EXPECT_EQ(ArmResult::COMPLETE, full.arms[kRight].status);
    // the plain walk still truncates and reports it
    Strategy plain;
    plain.max_extension_bp = kRadius;
    plain.max_seed_labels = 3;
    auto capped = run(*anno, f.X, {}, plain);
    EXPECT_EQ(3u, capped.num_seed_labels);
    EXPECT_EQ(5u, capped.labels_supporting_total);
    EXPECT_EQ(2u, capped.labels_dropped);
}


// Annotate mode: the live-label counts (frontier_remaining, cap_trigger, growth) are
// taken over the CUT lists, so they are lower bounds whenever a counted head was cut,
// and say so; with a cap that fits they are exact. Constrain mode never cuts.
TEST(Trie, CutListsMakeLiveLabelCountsInexact) {
    const OracleFixture f(kFixtureSeed);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.sequences, f.labels, DeBruijnGraph::BASIC);
    for (size_t cap : { 2u, 5u }) {
        Strategy st = exhaustive_at(LabelMode::ANNOTATE, kRadius);
        st.direction = Strategy::RIGHT;
        st.max_labels_per_node = cap;
        // two heads per level inside the bubble: trips on the first head of level 25,
        // where the P head still carries A, C, D and S (S ends at |P| + 10 = 25)
        st.max_steps = 50;
        auto r = run(*anno, f.X, {}, st);
        const ArmResult &arm = r.arms[kRight];
        const std::string what = "cap " + std::to_string(cap);
        ASSERT_EQ(ArmResult::TRUNCATED, arm.status) << what;
        ASSERT_TRUE(arm.cap_trigger.has_value()) << what;
        EXPECT_EQ(25u, arm.complete_to_bp) << what;
        const bool exact = cap >= 5;
        EXPECT_EQ(exact, arm.cap_trigger->live_labels_exact) << what;
        EXPECT_EQ(exact, arm.frontier_live_labels_exact) << what;
        // the P side carries four labels, the Q side two: a cap of 2 cuts the P head
        // at every level, so every bin with a P head is inexact
        ASSERT_FALSE(arm.growth.empty());
        for (const GrowthBin &b : arm.growth) {
            if (b.max_live_paths)
                EXPECT_EQ(exact, b.live_labels_exact) << what << " bin " << b.from_bp;
        }
        if (exact) {
            EXPECT_EQ(5u, arm.growth.front().max_live_labels) << what;
        } else {
            EXPECT_LT(arm.growth.front().max_live_labels, 5u) << what;
            EXPECT_LT(arm.cap_trigger->live_labels, 5u) << what;
        }
    }
    Strategy sc = exhaustive_at(LabelMode::CONSTRAIN, kRadius);
    sc.direction = Strategy::RIGHT;
    sc.max_labels_per_node = 2;        // cuts the trie-view branches, never the state
    sc.max_steps = 50;
    auto r = run(*anno, f.X, as_list(kAllLabels), sc);
    const ArmResult &arm = r.arms[kRight];
    ASSERT_TRUE(arm.cap_trigger.has_value());
    EXPECT_GT(arm.nodes_labels_truncated, 0u);
    EXPECT_TRUE(arm.cap_trigger->live_labels_exact);
    EXPECT_TRUE(arm.frontier_live_labels_exact);
    for (const GrowthBin &b : arm.growth) EXPECT_TRUE(b.live_labels_exact);
    EXPECT_EQ(5u, arm.growth.front().max_live_labels);
}


// The walk rule states what the walker enforces on THIS index: seed nodes of either
// strand and the even-k hairpin rule in the stranded regimes, the hash identity of
// edges when a (k+1)-mer does not pack into 64 bits.
TEST(Trie, WalkRuleStatesWhatTheWalkerEnforces) {
    auto b = fork_blocks(8);
    auto has = [](const std::string &s, const char *needle) {
        return s.find(needle) != std::string::npos;
    };
    Strategy st;
    {
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, { b[0] + b[1] }, { "A" }, DeBruijnGraph::BASIC);
        LabelOracle oracle(*anno);
        const std::string s = walk_rule_statement(st, oracle);
        EXPECT_TRUE(has(s, "no (k+1)-mer edge twice")) << s;
        EXPECT_TRUE(has(s, "enters no seed node, and")) << s;
        EXPECT_FALSE(has(s, "either strand")) << s;
        EXPECT_FALSE(has(s, "canonical")) << s;
        EXPECT_FALSE(has(s, "even k")) << s;
        EXPECT_FALSE(has(s, "128-bit")) << s;
        // a basic graph holds one strand: check_structure() never tests for hairpins
        // there, and the statement says so whatever the knob
        EXPECT_TRUE(has(s, "single-strand graph has no hairpins")) << s;
        EXPECT_FALSE(has(s, "hairpin step")) << s;
        st.skip_hairpins = false;
        EXPECT_TRUE(has(walk_rule_statement(st, oracle), "single-strand graph has no hairpins")) << s;
        st.skip_hairpins = true;
        // merging unites the edge histories of the routes it joins, so the set of walks
        // the certificate quantifies over is qualified whenever merging is on
        EXPECT_TRUE(has(s, "edge histories are united")) << s;
        st.merge_reconverge = false;
        EXPECT_FALSE(has(walk_rule_statement(st, oracle), "edge histories are united")) << s;
        st.merge_reconverge = true;
        st.label_mode = LabelMode::ANNOTATE;
        EXPECT_TRUE(has(walk_rule_statement(st, oracle), "labels present at its nodes are recorded"));
        st.label_mode = LabelMode::CONSTRAIN;
    }
#if ! _PROTEIN_GRAPH
    for (auto mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        // odd k: both strands of the seed count, a self-RC k-mer node does not
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, { b[0] + b[1] }, { "A" }, mode);
        LabelOracle oracle(*anno);
        const std::string s = walk_rule_statement(st, oracle);
        EXPECT_TRUE(has(s, "canonical (k+1)-mers")) << s;
        EXPECT_TRUE(has(s, "of either strand")) << s;
        EXPECT_FALSE(has(s, "even k")) << s;
        // even k: a step into or out of a self-RC k-mer node is a hairpin too
        auto even = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK - 1, { b[0] + b[1] }, { "A" }, mode);
        LabelOracle even_oracle(*even);
        EXPECT_TRUE(has(walk_rule_statement(st, even_oracle), "self-reverse-complementary k-mer node (even k)"));
    }
#endif
    {
        // (k+1) * 2 bits > 64: edges are identified by a hash, a collision blocks
        const std::string seq = random_seq(80, 9);
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                32, { seq }, { "A" }, DeBruijnGraph::BASIC);
        LabelOracle oracle(*anno);
        EXPECT_TRUE(has(walk_rule_statement(st, oracle), "128-bit FNV-1a hash"));
        auto packable = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                31, { seq }, { "A" }, DeBruijnGraph::BASIC);
        LabelOracle packable_oracle(*packable);
        EXPECT_FALSE(has(walk_rule_statement(st, packable_oracle), "128-bit"));
    }
}
