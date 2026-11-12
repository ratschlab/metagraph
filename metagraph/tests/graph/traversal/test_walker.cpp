#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>

#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "tests/graph/traversal/walker_test_helpers.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/traversal/resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/coord_to_header.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using namespace mtg::test;

std::vector<DeBruijnGraph::Mode> canonical_modes() {
    return {
#if ! _PROTEIN_GRAPH
        DeBruijnGraph::CANONICAL,
        DeBruijnGraph::PRIMARY,
#endif
    };
}

Strategy strategy(size_t max_label_branches = 0) {
    Strategy st;
    st.max_label_branches = max_label_branches;
    return st;
}

SeedResult run(const AnnotatedDBG &anno,
               const std::string &seq,
               const std::vector<std::string> &labels,
               const Strategy &st = Strategy()) {
    LabelOracle oracle(anno);
    Seed seed;
    seed.sequence = seq;
    seed.labels = labels;
    return traverse_seed(oracle, seed, st);
}

// every k-mer of seed + flank (or flank + seed) is in the graph
void check_spelled(const AnnotatedDBG &anno, const ArmResult &arm, const PathResult &path,
                   const std::string &seed) {
    std::string flank = spell_path(arm, path);
    EXPECT_EQ(path.length_bp, flank.size());
    std::string full = arm.arm == Arm::RIGHT ? seed + flank : flank + seed;
    for (node_index n : map_to_nodes_sequentially(anno.get_graph(), full)) {
        ASSERT_NE(npos, n) << full;
    }
}

size_t count_ends(const ArmResult &arm, EndReason reason) {
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &e : seg.label_ends) {
            n += e.reason == reason;
        }
    }
    return n;
}

const LabelEnd* label_end(const ArmResult &arm, LabelId label) {
    for (const auto &seg : arm.segments) {
        for (const auto &e : seg.label_ends) {
            if (e.label == label)
                return &e;
        }
    }
    return nullptr;
}

size_t path_ends(const PathResult &path, EndReason reason) {
    return std::count_if(path.end_labels.begin(), path.end_labels.end(),
                         [&](const PathEnd &e) { return e.reason == reason; });
}

// The invariants of every arm: steps are the bases added, every label end of a leaf is on
// its path, and a path's end labels are the labels at its leaf's last node
void check_invariants(const ArmResult &arm) {
    uint64_t bases = 0;
    for (const auto &seg : arm.segments) {
        bases += seg.length_bp;
        EXPECT_LE(seg.parents.size(), 1u);
        for (size_t c : seg.children) {
            EXPECT_EQ(std::vector<size_t>{ seg.id }, arm.segments[c].parents);
        }
    }
    EXPECT_EQ(arm.steps, bases);
    for (const auto &path : arm.paths) {
        const Segment &leaf = arm.segments[path.leaf];
        EXPECT_TRUE(leaf.children.empty());
        EXPECT_EQ(leaf.from_bp + leaf.length_bp, path.length_bp);
        std::vector<LabelId> labels;
        for (const auto &e : path.end_labels) {
            labels.push_back(e.label);
            const LabelEnd *end = nullptr;
            for (const auto &x : leaf.label_ends) {
                if (x.label == e.label && x.at_bp == path.length_bp)
                    end = &x;
            }
            ASSERT_NE(nullptr, end) << "label " << e.label;
            EXPECT_EQ(end->reason, e.reason);
        }
        EXPECT_EQ(leaf.labels_end, labels);
    }
    EXPECT_EQ(arm.status == ArmResult::COMPLETE, !arm.cap_trigger.has_value());
}

std::string serialize(const SeedResult &r) {
    std::ostringstream os;
    // without labels_from_seed and its counters: a derived run is compared with a named one
    os << r.length_bp << '|' << r.num_kmers << '|';
    for (const LabelRef &l : r.label_dict) {
        os << to_string(l.kind) << ':' << l.name << ';';
    }
    os << '|';
    for (const auto &d : r.dropped_labels) {
        os << d.name;
        for (const KmerInterval &run : d.runs) {
            os << ' ' << run.begin << '-' << run.end;
        }
        os << ';';
    }
    for (const ArmResult &a : r.arms) {
        os << "\nARM " << to_string(a.arm) << ' ' << a.requested << ' ' << a.status << ' '
           << a.steps << ' ' << a.frontier_live_paths << ' '
           << a.complete_to_bp;
        if (a.cap_trigger)
            os << " cap " << to_string(a.cap_trigger->reason) << ' ' << a.cap_trigger->at_bp;
        for (const Segment &s : a.segments) {
            os << "\n S" << s.id << " p";
            for (size_t p : s.parents) {
                os << p << ',';
            }
            os << " c";
            for (size_t c : s.children) {
                os << c << ',';
            }
            os << ' ' << s.from_bp << ' ' << s.length_bp << ' ' << s.sequence << " ls";
            for (LabelId l : s.labels_start) {
                os << l << ',';
            }
            os << " le";
            for (LabelId l : s.labels_end) {
                os << l << ',';
            }
            for (const LabelEnd &e : s.label_ends) {
                os << "\n  E" << e.at_bp << ' ' << e.label << ' ' << to_string(e.reason)
                   << ' ' << e.text;
            }
        }
        for (const Split &sp : a.splits) {
            os << "\n X" << sp.at_bp << ' ' << sp.segment << ' ' << sp.ambiguous << ' ';
            for (size_t c : sp.children) {
                os << c << ',';
            }
        }
        for (const PathResult &p : a.paths) {
            os << "\n P" << p.id << ' ' << p.leaf << ' ' << p.length_bp << ' '
               << (p.path_reason ? to_string(*p.path_reason) : "-") << ' ';
            for (const auto &e : p.end_labels) {
                os << e.label << ':' << to_string(e.reason) << ':' << e.branches << ',';
            }
        }
    }
    return os.str();
}


template <typename Pair>
class WalkerTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSSHash, annot::ColumnCompressed<>>
> WalkerTypes;
TYPED_TEST_SUITE(WalkerTest, WalkerTypes);


// one record, a seed in its middle: both arms run to the record's ends
TYPED_TEST(WalkerTest, Linear) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(200, 1);
    std::string seed = S.substr(60, 40);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto res = run(*anno, seed, { "A" });
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_TRUE(res.dropped_labels.empty());
        EXPECT_EQ(40u, res.length_bp);
        EXPECT_EQ(30u, res.num_kmers);

        const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
        for (const ArmResult *arm : { &left, &right }) {
            ASSERT_TRUE(arm->requested);
            EXPECT_EQ(ArmResult::COMPLETE, arm->status);
            EXPECT_EQ(Strategy().max_extension_bp, arm->complete_to_bp);
            ASSERT_EQ(1u, arm->paths.size());
            ASSERT_EQ(1u, arm->segments.size());
            const PathResult &path = arm->paths[0];
            EXPECT_FALSE(path.path_reason.has_value());
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(0u, path.end_labels[0].label);
            EXPECT_EQ(EndReason::DEAD_END, path.end_labels[0].reason);
            EXPECT_EQ(0u, path.end_labels[0].branches);
            EXPECT_EQ((std::vector<LabelId>{ 0 }), arm->segments[0].labels_start);
            EXPECT_EQ((std::vector<LabelId>{ 0 }), arm->segments[0].labels_end);
            EXPECT_TRUE(arm->splits.empty());
            check_spelled(*anno, *arm, path, seed);
            check_invariants(*arm);
            EXPECT_EQ(arm->steps, path.length_bp);
        }
        EXPECT_EQ(S.substr(100), spell_path(right, right.paths[0])) << mode;
        EXPECT_EQ(S.substr(0, 60), spell_path(left, left.paths[0])) << mode;

        // single-direction requests produce one arm only
        Strategy st;
        st.direction = Strategy::LEFT;
        auto only_left = run(*anno, seed, { "A" }, st);
        EXPECT_FALSE(only_left.arms[kRight].requested);
        EXPECT_TRUE(only_left.arms[kRight].paths.empty());
        EXPECT_TRUE(only_left.arms[kRight].segments.empty());
        ASSERT_EQ(1u, only_left.arms[kLeft].paths.size());
        EXPECT_EQ(S.substr(0, 60), spell_path(only_left.arms[kLeft], only_left.arms[kLeft].paths[0]));
        st.direction = Strategy::RIGHT;
        auto only_right = run(*anno, seed, { "A" }, st);
        EXPECT_FALSE(only_right.arms[kLeft].requested);
        ASSERT_EQ(1u, only_right.arms[kRight].paths.size());
        EXPECT_EQ(S.substr(100), spell_path(only_right.arms[kRight], only_right.arms[kRight].paths[0]));

        // a lower-case seed maps the same way and gives the same result
        std::string lower = seed;
        for (char &c : lower) {
            c = std::tolower(c);
        }
        auto lc = run(*anno, lower, { "A" });
        EXPECT_EQ(serialize(res), serialize(lc));
    }
}

// exact label-end positions: A on W·R·X, B on R·X·Y, seed R
TYPED_TEST(WalkerTest, LabelEndCoordinates) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    LabelEndFixture f;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        auto res = run(*anno, f.R, { "A", "B" });
        const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
        check_invariants(right);
        check_invariants(left);
        // right: both carry X, A ends after it, B goes on through Y
        ASSERT_EQ(1u, right.paths.size()) << mode;
        EXPECT_EQ(f.X + f.Y, spell_path(right, right.paths[0]));
        const LabelEnd *a = label_end(right, 0), *b = label_end(right, 1);
        ASSERT_TRUE(a && b);
        EXPECT_EQ(f.X.size(), a->at_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, a->reason);
        EXPECT_EQ(f.X.size() + f.Y.size(), b->at_bp);
        EXPECT_EQ(EndReason::DEAD_END, b->reason);
        ASSERT_EQ(1u, right.paths[0].end_labels.size());
        EXPECT_EQ(1u, right.paths[0].end_labels[0].label);
        // left: only A carries W; B is lost at the seed boundary
        ASSERT_EQ(1u, left.paths.size());
        EXPECT_EQ(f.W, spell_path(left, left.paths[0]));
        EXPECT_EQ(0u, label_end(left, 1)->at_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, label_end(left, 1)->reason);
        EXPECT_EQ(f.W.size(), label_end(left, 0)->at_bp);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), left.segments[0].labels_start);
        EXPECT_EQ((std::vector<LabelId>{ 0 }), left.segments[0].labels_end);
    }
}

// two labels diverging at the seed boundary: a split, one label per child
TYPED_TEST(WalkerTest, Divergence) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(3);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);
        auto res = run(*anno, X, { "A", "B" });
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        ASSERT_EQ(1u, right.splits.size());
        EXPECT_EQ(0u, right.splits[0].at_bp);
        EXPECT_EQ(0u, right.splits[0].segment);
        EXPECT_FALSE(right.splits[0].ambiguous);
        ASSERT_EQ(2u, right.splits[0].children.size());
        ASSERT_EQ(3u, right.segments.size());
        EXPECT_EQ(0u, right.segments[0].length_bp);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), right.segments[0].labels_end);
        ASSERT_EQ(2u, right.paths.size());
        std::map<std::string, LabelId> leaves;
        for (const auto &path : right.paths) {
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(EndReason::DEAD_END, path.end_labels[0].reason);
            EXPECT_EQ(0u, path.end_labels[0].branches);
            leaves[spell_path(right, path)] = path.end_labels[0].label;
            check_spelled(*anno, right, path, X);
        }
        ASSERT_EQ(2u, leaves.size());
        EXPECT_EQ(0u, leaves.at(P));
        EXPECT_EQ(1u, leaves.at(Q));
        EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
        EXPECT_EQ(80u, right.steps);
    }
}

// one label on both sides of a structural branch
TYPED_TEST(WalkerTest, AmbiguousBranch) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t seed = 4; ; ++seed) {
        b = clean_blocks({ 30, 25, 30, 25, 30 }, seed);
        if (b[1][0] != b[3][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Y = b[2], &Q = b[3], &Z = b[4];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P + Y, X + Q + Z },
                                                        { "A", "A" }, mode);
        // limit 0: the label stops at the branch
        auto res = run(*anno, X, { "A" }, strategy(0));
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(0u, right.paths[0].length_bp);
        EXPECT_EQ(1u, path_ends(right.paths[0], EndReason::BRANCH));
        EXPECT_EQ(0u, label_end(right, 0)->at_bp);
        EXPECT_TRUE(right.splits.empty());
        EXPECT_EQ(0u, right.steps);

        // limit 1 and unlimited: both continuations are followed, each charged one branch
        for (size_t limit : { size_t(1), Strategy::kUnlimited }) {
            res = run(*anno, X, { "A" }, strategy(limit));
            const ArmResult &r1 = res.arms[kRight];
            check_invariants(r1);
            ASSERT_EQ(2u, r1.paths.size());
            std::set<std::string> flanks;
            for (const auto &path : r1.paths) {
                flanks.insert(spell_path(r1, path));
                ASSERT_EQ(1u, path.end_labels.size());
                EXPECT_EQ(1u, path.end_labels[0].branches);
                EXPECT_EQ(EndReason::DEAD_END, path.end_labels[0].reason);
                check_spelled(*anno, r1, path, X);
            }
            EXPECT_EQ((std::set<std::string>{ P + Y, Q + Z }), flanks);
            ASSERT_EQ(1u, r1.splits.size());
            EXPECT_TRUE(r1.splits[0].ambiguous);
            EXPECT_EQ(0u, count_ends(r1, EndReason::BRANCH));
        }
    }
}

// a cycle entered through junction J: the walk rule lets it be walked once
TYPED_TEST(WalkerTest, CycleJunction) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t seed = 6; ; ++seed) {
        b = clean_blocks({ 30, kK, 40, 30 }, seed);
        if (b[2][0] != b[3][0])
            break;
    }
    const std::string &L = b[0], &J = b[1], &C = b[2], &E = b[3];
    std::string S = L + J + C + J + E;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto res = run(*anno, L, { "A" }, strategy(0));
        const ArmResult &r0 = res.arms[kRight];
        check_invariants(r0);
        ASSERT_EQ(1u, r0.paths.size());
        EXPECT_EQ(J, spell_path(r0, r0.paths[0]));
        EXPECT_EQ(1u, count_ends(r0, EndReason::BRANCH));

        res = run(*anno, L, { "A" }, strategy(1));
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1);
        ASSERT_EQ(2u, r1.paths.size());
        std::set<std::string> flanks;
        for (const auto &path : r1.paths) {
            flanks.insert(spell_path(r1, path));
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(1u, path.end_labels[0].branches);
            EXPECT_EQ(EndReason::DEAD_END, path.end_labels[0].reason);
            check_spelled(*anno, r1, path, L);
        }
        // the second time at J the cycle's first edge was taken: only the exit is followed
        EXPECT_EQ((std::set<std::string>{ J + E, J + C + J + E }), flanks);
        EXPECT_EQ(0u, count_ends(r1, EndReason::BRANCH));
        EXPECT_EQ(0u, count_ends(r1, EndReason::EDGE_REUSE));
    }
}

// an inverted repeat Y ... rc(Y): on a canonical or primary graph the walk re-enters Y's
// edges in the opposite orientation and ends with EDGE_REUSE_RC; a basic graph has no such edge
TYPED_TEST(WalkerTest, InvertedRepeatReusesEdgesInReverse) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 30, 25, 20, 30 }, 77);
    const std::string &L = b[0], &Y = b[1], &B = b[2], &E = b[3];
    const std::string S = L + Y + B + rc(Y) + E;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto res = run(*anno, L, { "A" }, strategy(Strategy::kUnlimited));
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        EXPECT_EQ(0u, count_ends(right, EndReason::EDGE_REUSE)) << mode;
        if (mode == DeBruijnGraph::BASIC) {
            EXPECT_EQ(0u, count_ends(right, EndReason::EDGE_REUSE_RC)) << mode;
        } else {
            EXPECT_GE(count_ends(right, EndReason::EDGE_REUSE_RC), 1u) << mode;
        }
    }
}

// a homopolymer self-loop: each loop step re-takes the loop edge
TYPED_TEST(WalkerTest, SelfLoop) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t seed = 7; ; ++seed) {
        b = clean_blocks({ 30, 30 }, seed);
        if (b[0].back() != 'A' && b[1][0] != 'A' && b[0].back() != 'T')
            break;
    }
    const std::string &L = b[0], &E = b[1];
    std::string H(kK + 1, 'A');
    std::string S = L + H + E;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto res = run(*anno, L, { "A" }, strategy(0));
        const ArmResult &r0 = res.arms[kRight];
        check_invariants(r0);
        ASSERT_EQ(1u, r0.paths.size());
        EXPECT_EQ(kK - 1, r0.paths[0].length_bp);
        EXPECT_EQ(1u, count_ends(r0, EndReason::BRANCH));
        check_spelled(*anno, r0, r0.paths[0], L);

        res = run(*anno, L, { "A" }, strategy(1));
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1);
        ASSERT_EQ(2u, r1.paths.size());
        std::set<std::string> flanks;
        for (const auto &path : r1.paths) {
            flanks.insert(spell_path(r1, path));
            check_spelled(*anno, r1, path, L);
        }
        EXPECT_EQ((std::set<std::string>{ H.substr(1), H.substr(2) + E }), flanks);
        EXPECT_EQ(1u, count_ends(r1, EndReason::BRANCH));
        EXPECT_EQ(1u, count_ends(r1, EndReason::DEAD_END));

        res = run(*anno, L, { "A" }, strategy(2));
        const ArmResult &r2 = res.arms[kRight];
        check_invariants(r2);
        ASSERT_EQ(3u, r2.paths.size());
        flanks.clear();
        for (const auto &path : r2.paths) {
            flanks.insert(spell_path(r2, path));
            check_spelled(*anno, r2, path, L);
            EXPECT_EQ(1u, path_ends(path, EndReason::DEAD_END));
        }
        EXPECT_EQ((std::set<std::string>{ H + E, H.substr(1) + E, H.substr(2) + E }), flanks);
        EXPECT_EQ(0u, count_ends(r2, EndReason::BRANCH));
    }
}

// a circular record: both arms go round and stop where they would step back into the seed
TYPED_TEST(WalkerTest, Circle) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(150, 8);
    std::string record = S + S.substr(0, kK - 1 + 10);
    std::string seed = S.substr(50, 40);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { record }, { "A" }, mode);
        auto res = run(*anno, seed, { "A" });
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm);
            ASSERT_EQ(1u, arm.paths.size()) << mode << " arm " << a;
            EXPECT_EQ(ArmResult::COMPLETE, arm.status);
            EXPECT_EQ(1u, path_ends(arm.paths[0], EndReason::REACHED_SEED));
            EXPECT_EQ(S.size() - seed.size() + kK - 1, arm.paths[0].length_bp);
            check_spelled(*anno, arm, arm.paths[0], seed);
        }
        EXPECT_EQ(S.substr(90) + S.substr(0, 60), spell_path(res.arms[kRight], res.arms[kRight].paths[0]));
        EXPECT_EQ(S.substr(80) + S.substr(0, 50), spell_path(res.arms[kLeft], res.arms[kLeft].paths[0]));
    }
}

// a self-reverse-complementary stretch: the hairpin step back is never taken
TYPED_TEST(WalkerTest, Hairpin) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string X = "AACGTACGTT";
    ASSERT_EQ(X, rc(X));
    ASSERT_EQ(kK - 1, X.size());
    std::string L, R;
    for (uint32_t seed = 9; ; ++seed) {
        L = clean_block(40, seed);
        R = clean_block(40, seed + 100);
        if (R[0] != rc(std::string(1, L.back()))[0] && is_clean(L + X.substr(0, 3))
                && is_clean(X.substr(7) + R))
            break;
    }
    std::string S = L + X + R;
    std::string seed = L.substr(0, 30);
    for (auto mode : canonical_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto res = run(*anno, seed, { "A" });
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size()) << mode;
        EXPECT_EQ(S.substr(30), spell_path(right, right.paths[0])) << mode;
        EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
        EXPECT_TRUE(right.splits.empty());
    }
}

// CLI-built DBGSuccinct graphs keep their dummy k-mers: the walk is the same
TEST(Walker, UnmaskedSuccinct) {
    LabelEndFixture f;
    for (auto mode : all_modes()) {
        auto masked = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, f.seqs(), f.labels(), mode);
        std::shared_ptr<DeBruijnGraph> graph;
        std::shared_ptr<DBGSuccinct> base;
        if (mode == DeBruijnGraph::PRIMARY) {
            auto canonical = build_graph<DBGSuccinct>(kK, f.seqs(), DeBruijnGraph::CANONICAL);
            std::vector<std::string> contigs;
            canonical->call_sequences([&](const std::string &c, const auto &) { contigs.push_back(c); }, 1, true);
            base = std::make_shared<DBGSuccinct>(kK, DeBruijnGraph::PRIMARY);
            for (const auto &c : contigs) {
                base->add_sequence(c);
            }
            graph = std::make_shared<CanonicalDBG>(base);
        } else {
            base = std::make_shared<DBGSuccinct>(kK, mode);
            for (const auto &s : f.seqs()) {
                base->add_sequence(s);
            }
            graph = base;
        }
        ASSERT_EQ(nullptr, base->get_mask());
        auto unmasked = std::make_unique<AnnotatedDBG>(
            graph, std::make_unique<annot::ColumnCompressed<>>(base->max_index()));
        for (size_t i = 0; i < f.seqs().size(); ++i) {
            unmasked->annotate_sequence(f.seqs()[i], { f.labels()[i] });
        }
        for (const std::string &seed : { f.R, f.X, f.W + f.R }) {
            std::vector<std::string> labels = seed == f.W + f.R ? std::vector<std::string>{ "A" }
                                                                : std::vector<std::string>{ "A", "B" };
            EXPECT_EQ(serialize(run(*masked, seed, labels)), serialize(run(*unmasked, seed, labels)))
                << mode << " seed " << seed;
        }
    }
}

TEST(Walker, InconsistentPrimaryGraph) {
    std::string S = clean_block(80, 17);
    auto base = std::make_shared<DBGSuccinct>(kK, DeBruijnGraph::PRIMARY);
    base->add_sequence(S);
    base->add_sequence(rc(S));
    base->mask_dummy_kmers(1, false);
    auto graph = std::make_shared<CanonicalDBG>(base);
    auto anno = std::make_unique<AnnotatedDBG>(
        graph, std::make_unique<annot::ColumnCompressed<>>(base->max_index()));
    anno->annotate_sequence(S, { "A" });
    LabelOracle oracle(*anno);
    ASSERT_EQ(Regime::PRIMARY, oracle.regime());
    Seed seed;
    seed.sequence = S.substr(20, 30);
    seed.labels = { "A" };
    try {
        traverse_seed(oracle, seed, Strategy());
        FAIL() << "no exception";
    } catch (const std::runtime_error &e) {
        std::string what = e.what();
        EXPECT_NE(std::string::npos, what.find("Inconsistent primary graph at node")) << what;
        EXPECT_NE(std::string::npos, what.find("on arm")) << what;
    }
}

TYPED_TEST(WalkerTest, Caps) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(200, 10);
    std::string seed = S.substr(60, 40);
    auto b = fork_blocks(11);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto linear = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto fork = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);

        // the extension length completes the arm
        Strategy st;
        st.max_extension_bp = 10;
        auto res = run(*linear, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm);
            EXPECT_EQ(ArmResult::COMPLETE, arm.status);
            EXPECT_EQ(10u, arm.complete_to_bp);
            ASSERT_EQ(1u, arm.paths.size());
            EXPECT_EQ(10u, arm.paths[0].length_bp);
            ASSERT_TRUE(arm.paths[0].path_reason.has_value());
            EXPECT_EQ(EndReason::MAX_EXTENSION, *arm.paths[0].path_reason);
            EXPECT_EQ(1u, path_ends(arm.paths[0], EndReason::MAX_EXTENSION));
            EXPECT_EQ(10u, arm.steps);
        }

        // the steps cap stops the seed, both arms
        st = Strategy();
        st.max_steps = 5;
        res = run(*linear, seed, { "A" }, st);
        EXPECT_EQ(5u, res.arms[kLeft].steps + res.arms[kRight].steps);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm);
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
            ASSERT_TRUE(arm.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_STEPS, arm.cap_trigger->reason);
            EXPECT_EQ(arm.cap_trigger->at_bp, arm.complete_to_bp);
            EXPECT_EQ(1u, arm.frontier_live_paths);
            EXPECT_EQ(1u, arm.frontier_live_labels);
            // the stopped head's segment is a leaf
            ASSERT_LT(arm.cap_trigger->segment, arm.segments.size());
            EXPECT_TRUE(arm.segments[arm.cap_trigger->segment].children.empty());
            // a bound ended the label: the node's successors were not looked at
            for (const Segment &s : arm.segments) {
                for (const LabelEnd &e : s.label_ends) {
                    EXPECT_FALSE(e.successors.has_value());
                }
            }
            ASSERT_EQ(1u, arm.paths.size());
            EXPECT_EQ(EndReason::MAX_STEPS, *arm.paths[0].path_reason);
            EXPECT_EQ(1u, count_ends(arm, EndReason::MAX_STEPS));
        }

        // the live-path and path caps stop one arm only
        st = Strategy();
        st.max_live_paths = 1;
        res = run(*fork, X, { "A", "B" }, st);
        {
            const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
            check_invariants(right);
            EXPECT_EQ(ArmResult::TRUNCATED, right.status);
            ASSERT_TRUE(right.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_LIVE_PATHS, right.cap_trigger->reason);
            EXPECT_EQ(0u, right.cap_trigger->at_bp);
            EXPECT_EQ(0u, right.complete_to_bp);
            EXPECT_EQ(2u, count_ends(right, EndReason::MAX_LIVE_PATHS));
            EXPECT_TRUE(right.splits.empty());
            EXPECT_EQ(ArmResult::COMPLETE, left.status);
            EXPECT_EQ(2u, count_ends(left, EndReason::DEAD_END));
        }
        // a fork at the last level needs no further expansion, so it is not capped
        st = Strategy();
        st.max_live_paths = 1;
        st.max_extension_bp = 1;
        st.direction = Strategy::RIGHT;
        res = run(*fork, X, { "A", "B" }, st);
        EXPECT_EQ(ArmResult::COMPLETE, res.arms[kRight].status);
        EXPECT_FALSE(res.arms[kRight].cap_trigger.has_value());
        EXPECT_EQ(2u, res.arms[kRight].paths.size());
        st = Strategy();
        st.max_paths = 1;
        res = run(*fork, X, { "A", "B" }, st);
        EXPECT_EQ(ArmResult::TRUNCATED, res.arms[kRight].status);
        ASSERT_TRUE(res.arms[kRight].cap_trigger.has_value());
        EXPECT_EQ(EndReason::MAX_PATHS, res.arms[kRight].cap_trigger->reason);
        EXPECT_EQ(2u, count_ends(res.arms[kRight], EndReason::MAX_PATHS));
        EXPECT_EQ(ArmResult::COMPLETE, res.arms[kLeft].status);
        st.max_paths = 2;
        res = run(*fork, X, { "A", "B" }, st);
        EXPECT_EQ(ArmResult::COMPLETE, res.arms[kRight].status);
        EXPECT_EQ(2u, res.arms[kRight].paths.size());

        // a zero time budget still delivers the first level
        st = Strategy();
        st.time_budget_ms = 0;
        res = run(*linear, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
            EXPECT_EQ(1u, arm.steps);
            ASSERT_TRUE(arm.cap_trigger.has_value());
            EXPECT_EQ(EndReason::TIME_BUDGET, arm.cap_trigger->reason);
            EXPECT_EQ(1u, arm.cap_trigger->at_bp);
            EXPECT_EQ(1u, arm.complete_to_bp);
            EXPECT_EQ(1u, count_ends(arm, EndReason::TIME_BUDGET));
        }
    }
}

// a head that reached the extension length when a cap trips ends complete, not censored
TEST(Walker, HeadsAtTheRadiusAreCompleteWhenACapTrips) {
    auto b = fork_blocks(62);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { X + P, X + Q }, { "A", "B" }, DeBruijnGraph::BASIC);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.max_extension_bp = 2;
    st.max_steps = 3;
    auto res = run(*anno, X, { "A", "B" }, st);
    const ArmResult &arm = res.arms[kRight];
    EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
    EXPECT_EQ(1u, arm.complete_to_bp);
    ASSERT_EQ(2u, arm.paths.size());
    std::map<uint64_t, EndReason> ends;
    for (const auto &p : arm.paths) {
        ASSERT_EQ(1u, p.end_labels.size());
        ends[p.length_bp] = p.end_labels[0].reason;
    }
    ASSERT_EQ(2u, ends.size());
    EXPECT_STREQ("max_extension_bp", to_string(ends.at(2)));
    EXPECT_STREQ("max_steps", to_string(ends.at(1)));
    EXPECT_EQ(1u, arm.frontier_live_paths);
    ASSERT_TRUE(arm.cap_trigger.has_value());
    EXPECT_EQ(1u, arm.cap_trigger->at_bp);
}

// the split allowance of a path is a cap: the labels it cuts end with MAX_SPLITS, and the arm
// is complete only up to the fork
TYPED_TEST(WalkerTest, SplitLimit) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(3);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);
        Strategy st;
        st.max_splits_per_path = 0;
        auto res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size()) << mode;
        EXPECT_EQ(0u, right.paths[0].length_bp);
        EXPECT_EQ(2u, path_ends(right.paths[0], EndReason::MAX_SPLITS));
        ASSERT_TRUE(right.paths[0].path_reason.has_value());
        EXPECT_EQ(EndReason::MAX_SPLITS, *right.paths[0].path_reason);
        EXPECT_TRUE(right.splits.empty());
        EXPECT_EQ(0u, right.steps);
        EXPECT_EQ(ArmResult::TRUNCATED, right.status);
        EXPECT_EQ(0u, right.complete_to_bp);
        ASSERT_TRUE(right.cap_trigger);
        EXPECT_EQ(EndReason::MAX_SPLITS, right.cap_trigger->reason);
        EXPECT_EQ(0u, right.cap_trigger->at_bp);
        for (LabelId l : { 0u, 1u }) {
            const LabelEnd *e = label_end(right, l);
            ASSERT_NE(nullptr, e);
            EXPECT_EQ(0u, e->at_bp);
            EXPECT_EQ(EndReason::MAX_SPLITS, e->reason);
            // the fork has both successors
            EXPECT_EQ(std::optional<size_t>(2), e->successors);
        }
        st.max_splits_per_path = 1;
        auto one = run(*anno, X, { "A", "B" }, st);
        EXPECT_EQ(2u, one.arms[kRight].paths.size());
        EXPECT_EQ(ArmResult::COMPLETE, one.arms[kRight].status);

        // two forks in a row: the path of A and B splits from C after X and again after Y
        auto c = fork_blocks(7);
        const std::string &Y = c[0], &R = c[1], &T = c[2];
        auto two = build_anno_graph<Graph, Annotation>(
                kK, { X + P + Y + R, X + P + Y + T, X + Q }, { "A", "B", "C" }, mode);
        st.max_splits_per_path = 1;
        auto cut = run(*two, X, { "A", "B", "C" }, st);
        check_invariants(cut.arms[kRight]);
        EXPECT_EQ(ArmResult::TRUNCATED, cut.arms[kRight].status);
        EXPECT_EQ(P.size() + Y.size(), cut.arms[kRight].complete_to_bp);
        for (size_t limit : { size_t(2), Strategy::kUnlimited }) {
            st.max_splits_per_path = limit;
            auto all = run(*two, X, { "A", "B", "C" }, st);
            check_invariants(all.arms[kRight]);
            EXPECT_EQ(ArmResult::COMPLETE, all.arms[kRight].status) << mode;
            EXPECT_EQ(3u, all.arms[kRight].paths.size()) << mode;
        }
    }
}

TYPED_TEST(WalkerTest, SeedValidation) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    LabelEndFixture f;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        // A does not carry the seed's last k-mers: dropped, with the k-mers it carries
        std::string seed = f.R + f.X + f.Y.substr(0, 5);
        auto res = run(*anno, seed, { "A", "B" });
        ASSERT_EQ(1u, res.dropped_labels.size());
        EXPECT_EQ("A", res.dropped_labels[0].name);
        uint64_t supported = f.R.size() + f.X.size() - kK + 1;
        EXPECT_EQ((std::vector<KmerInterval>{ { 0, supported } }),
                  res.dropped_labels[0].runs);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ("B", res.label_dict[0].name);
        EXPECT_EQ(1u, res.labels_supporting_total);
        EXPECT_EQ(f.Y.substr(5), spell_path(res.arms[kRight], res.arms[kRight].paths[0]));

        try {
            run(*anno, f.W + f.R, { "B" });  // unsupported
            ADD_FAILURE() << "a seed no named label carries must fail";
        } catch (const SeedDerivationError &e) {
            EXPECT_EQ(SeedDerivationError::NO_CARRIER, e.cause());
        }
        EXPECT_THROW(run(*anno, random_seq(40, 99), { "A" }), SeedNotInGraph);
        EXPECT_THROW(run(*anno, f.R + "ACGTACGTACGTACGT", { "A", "B" }), SeedNotInGraph);
        EXPECT_THROW(run(*anno, f.R, { "C" }), std::invalid_argument);           // unknown
        EXPECT_THROW(run(*anno, f.R, { "A", "C" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R.substr(0, kK - 1), { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R.substr(0, 20) + "N" + f.R.substr(21, 20), { "A" }),
                     std::invalid_argument);
        EXPECT_THROW(run(*anno, "", { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R, { "A", "A" }), std::invalid_argument);     // duplicate
        EXPECT_EQ(2u, run(*anno, f.R, {}).label_dict.size());
        EXPECT_THROW(run(*anno, f.W + f.R + f.X + f.Y, {}), SeedDerivationError);
    }
}

// labels derived from the seed: its carriers, resubmittable as an explicit list
TYPED_TEST(WalkerTest, DeriveSeedLabels) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    DerivedFixture f;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        Strategy st;
        auto derived = run(*anno, f.R, {}, st);
        EXPECT_TRUE(derived.labels_from_seed);
        EXPECT_EQ(2u, derived.label_dict.size()) << mode;
        EXPECT_EQ(2u, derived.labels_supporting_total);
        EXPECT_EQ(0u, derived.labels_dropped);
        EXPECT_TRUE(derived.dropped_labels.empty());

        std::vector<std::string> names;
        for (const auto &l : derived.label_dict) {
            EXPECT_EQ(LabelKind::COLUMN, l.kind);   // no CoordToHeader: column labels
            names.push_back(l.name);
        }
        EXPECT_EQ((std::set<std::string>{ "A", "B" }), std::set<std::string>(names.begin(), names.end()));
        for (size_t i = 1; i < derived.label_dict.size(); ++i) {   // ascending column id
            EXPECT_LT(derived.label_dict[i - 1].column, derived.label_dict[i].column);
        }

        auto named = run(*anno, f.R, names, st);
        EXPECT_FALSE(named.labels_from_seed);
        EXPECT_EQ(serialize(named), serialize(derived)) << mode;

        auto all_three = run(*anno, f.prefix(), {}, st);
        EXPECT_EQ(3u, all_three.label_dict.size()) << mode;
        EXPECT_EQ(3u, all_three.labels_supporting_total);

        st.max_seed_labels = 2;
        auto capped = run(*anno, f.prefix(), {}, st);
        ASSERT_EQ(2u, capped.label_dict.size()) << mode;
        EXPECT_EQ(3u, capped.labels_supporting_total);
        EXPECT_EQ(1u, capped.labels_dropped);
        EXPECT_TRUE(capped.dropped_labels.empty());   // capped, not unsupported
        EXPECT_EQ(capped.label_dict[0].name, all_three.label_dict[0].name);
        EXPECT_EQ(capped.label_dict[1].name, all_three.label_dict[1].name);
        st.max_seed_labels = 0;   // a cap of zero would permit nothing
        EXPECT_THROW(run(*anno, f.prefix(), {}, st), std::invalid_argument);

        try {
            run(*anno, f.W + f.R + f.X + f.Y, {});
            FAIL() << "no carrier";
        } catch (const SeedDerivationError &e) {
            EXPECT_EQ(SeedDerivationError::NO_CARRIER, e.cause());
        }
    }
}

TYPED_TEST(WalkerTest, Determinism) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(5);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q, X + P + Q },
                                                        { "A", "B", "C" }, mode);
        const std::string first = serialize(run(*anno, X, { "A", "B", "C" }, strategy(1)));
        for (int i = 0; i < 2; ++i) {
            EXPECT_EQ(first, serialize(run(*anno, X, { "A", "B", "C" }, strategy(1))));
        }
        // the label order of the request changes the ids, not the walk
        auto a = run(*anno, X, { "A", "B", "C" }, strategy(1));
        auto c = run(*anno, X, { "C", "B", "A" }, strategy(1));
        EXPECT_EQ(a.arms[kRight].paths.size(), c.arms[kRight].paths.size());
        EXPECT_EQ(a.arms[kRight].steps, c.arms[kRight].steps);
    }
}

TEST(WalkerDerive, ManyDerivedLabelsMatchTheExplicitList) {
    constexpr size_t kN = 24;
    std::vector<size_t> lengths { 30 };
    lengths.insert(lengths.end(), kN, 20);
    auto b = clean_blocks(lengths, 97);
    const std::string S = b[0];

    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < kN; ++i) {
        seqs.push_back(S + b[1 + i]);
        labels.push_back("full" + std::to_string(i));
    }
    seqs.push_back(S.substr(0, S.size() - 10));
    labels.push_back("partial");
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, labels);

    auto derived = run(*anno, S, {}, strategy(Strategy::kUnlimited));
    EXPECT_TRUE(derived.labels_from_seed);
    ASSERT_EQ(kN, derived.label_dict.size());
    EXPECT_EQ(kN, derived.labels_supporting_total);
    EXPECT_EQ(kN, derived.arms[kRight].paths.size());

    std::vector<std::string> names;
    for (const auto &l : derived.label_dict) {
        names.push_back(l.name);
    }
    auto named = run(*anno, S, names, strategy(Strategy::kUnlimited));
    EXPECT_EQ(serialize(named), serialize(derived));

    names.push_back("partial");
    auto mixed = run(*anno, S, names, strategy(Strategy::kUnlimited));
    EXPECT_EQ(kN, mixed.label_dict.size());
    ASSERT_EQ(1u, mixed.dropped_labels.size());
    EXPECT_EQ("partial", mixed.dropped_labels[0].name);
    EXPECT_FALSE(mixed.dropped_labels[0].runs.empty());
}

TEST(WalkerDerive, DerivationHonoursTheTimeBudget) {
    DerivedFixture f;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.seqs(), f.labels(), DeBruijnGraph::BASIC);
    Strategy st;
    st.time_budget_ms = 1e-9;    // spent by the time the first seed k-mer is consumed
    try {
        run(*anno, f.R, {}, st);
        FAIL() << "the derivation ignored its time budget";
    } catch (const SeedDerivationError &e) {
        EXPECT_EQ(SeedDerivationError::TIME_BUDGET, e.cause());
        EXPECT_EQ(1.0, e.observed());
    }
    // explicit labels need no derivation: the walk stops at its first level instead
    SeedResult named = run(*anno, f.R, { "A", "B" }, st);
    EXPECT_EQ(2u, named.label_dict.size());
    EXPECT_EQ(ArmResult::TRUNCATED, named.arms[kRight].status);
}

// a seed starting in a repetitive k-mer: the intersection starts from the narrowest row
TEST(WalkerDerive, DerivationNarrowsAWideFirstKmerToItsCarriers) {
    constexpr size_t kDecoys = 100;
    auto b = clean_blocks({ kK, 30 }, 137);
    const std::string U = b[0], V = b[1];     // U is exactly one k-mer
    const std::string seed_seq = U + V;

    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < kDecoys; ++i) {
        seqs.push_back(U);
        labels.push_back("decoy" + std::to_string(i));
    }
    seqs.push_back(seed_seq);
    labels.push_back("carrierA");
    seqs.push_back(seed_seq);
    labels.push_back("carrierB");
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, seqs, labels, DeBruijnGraph::BASIC, true);

    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(labels.size());
    std::vector<std::vector<uint64_t>> num_kmers(labels.size());
    for (size_t i = 0; i < labels.size(); ++i) {
        Column c = static_cast<Column>(encoder.encode(labels[i]));
        headers[c] = { "h_" + labels[i] };
        num_kmers[c].assign(1, seqs[i].size() - kK + 1);
    }
    anno = with_coord_to_header(std::move(anno), std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno);
    Seed seed;
    seed.sequence = seed_seq;
    auto derived = traverse_seed(oracle, seed, Strategy());
    ASSERT_EQ(2u, derived.label_dict.size());
    std::set<std::string> names;
    for (const auto &l : derived.label_dict) {
        EXPECT_EQ(LabelKind::HEADER, l.kind);
        names.insert(l.name);
    }
    EXPECT_EQ((std::set<std::string>{ "h_carrierA", "h_carrierB" }), names);
}

// a header in two columns cannot be resubmitted by name: refused, with the lever
TEST(WalkerDerive, AmbiguousDerivedHeaderNamesAreRefused) {
    auto b = clean_blocks({ 30, 20 }, 211);
    const std::string S = b[0], T = b[1];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { S + T, S }, { "colA", "colB" }, DeBruijnGraph::BASIC, true);

    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[static_cast<Column>(encoder.encode("colA"))] = { "ACC1" };
    num_kmers[static_cast<Column>(encoder.encode("colA"))].assign(1, (S + T).size() - kK + 1);
    headers[static_cast<Column>(encoder.encode("colB"))] = { "ACC1" };
    num_kmers[static_cast<Column>(encoder.encode("colB"))].assign(1, S.size() - kK + 1);
    anno = with_coord_to_header(std::move(anno), std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno);
    Seed seed;
    seed.sequence = S;                  // carried in full by both columns
    Strategy st;
    try {
        traverse_seed(oracle, seed, st);
        FAIL() << "expected an ambiguous derived header set to be refused";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(nullptr, std::strstr(e.what(), "ACC1")) << e.what();
        EXPECT_NE(nullptr, std::strstr(e.what(), "seed_label_kind")) << e.what();
        EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
        EXPECT_EQ("ACC1", e.subject());
    }
    st.seed_label_kind = LabelKind::COLUMN;
    auto by_column = traverse_seed(oracle, seed, st);
    ASSERT_EQ(2u, by_column.label_dict.size());
    EXPECT_EQ("colA", by_column.label_dict[0].name);
    EXPECT_EQ("colB", by_column.label_dict[1].name);
    // named explicitly, the header of two records is refused as ambiguous
    seed.labels = { "ACC1" };
    EXPECT_THROW(traverse_seed(oracle, seed, Strategy()), std::invalid_argument);
}

/************ HEADER labels: a sequence record is followed in coordinate order ************/

template <typename Pair>
class WalkerRecordTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>
> WalkerRecordTypes;
TYPED_TEST_SUITE(WalkerRecordTest, WalkerRecordTypes);

// the one end of a label on an arm with a single path
const LabelEnd& only_end(const ArmResult &arm, LabelId label) {
    EXPECT_EQ(1u, arm.paths.size());
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &e : seg.label_ends) {
            n += e.label == label;
        }
    }
    EXPECT_EQ(1u, n);
    const LabelEnd *end = label_end(arm, label);
    static const LabelEnd none {};
    return end ? *end : none;
}

// A record whose two ends share a repeat longer than k. By membership the record's first
// and last k-mers have neighbours carrying the record (the other copy's), and both arms
// walked on into them. In coordinate order both arms end where the record does.
TYPED_TEST(WalkerRecordTest, TerminalRepeatEndsAtTheRecordEnds) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 30, 40, 40 }, 301);
    const std::string &T = b[0], &A = b[1], &B = b[2];
    const std::string record = T + A + B + T;
    const std::string seed = A.substr(5, 25);
    auto anno = build_records<Graph, Annotation>(kK, { { "F", "rec", record } });

    for (const std::vector<std::string> &labels : { std::vector<std::string>{ "rec" },
                                                    std::vector<std::string>{} }) {
        auto res = run(*anno, seed, labels);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ(LabelKind::HEADER, res.label_dict[0].kind);
        EXPECT_EQ("rec", res.label_dict[0].name);
        EXPECT_TRUE(res.dropped_labels.empty());
        const ArmResult &left = res.arms[kLeft], &right = res.arms[kRight];
        check_invariants(left);
        check_invariants(right);
        ASSERT_EQ(1u, left.paths.size());
        ASSERT_EQ(1u, right.paths.size());
        // exactly the record's bases on both sides of the seed
        EXPECT_EQ(T + A.substr(0, 5), spell_path(left, left.paths[0]));
        EXPECT_EQ(A.substr(30) + B + T, spell_path(right, right.paths[0]));
        for (const ArmResult *arm : { &left, &right }) {
            EXPECT_EQ(ArmResult::COMPLETE, arm->status);
            EXPECT_TRUE(arm->splits.empty());
            const LabelEnd &end = only_end(*arm, 0);
            EXPECT_EQ(EndReason::RECORD_END, end.reason) << to_string(end.reason);
            EXPECT_EQ(arm->paths[0].length_bp, end.at_bp);
            EXPECT_FALSE(arm->paths[0].path_reason.has_value());
            EXPECT_EQ(1u, path_ends(arm->paths[0], EndReason::RECORD_END));
        }
        EXPECT_STREQ("record_end", to_string(EndReason::RECORD_END));
    }
}

// A repeat inside one record, its two copies followed by different bases. By membership
// both successors after the first copy carry the record: a fork, BRANCH at allowance 0.
// In coordinate order the walk goes through both copies, re-using the repeat's edges.
TYPED_TEST(WalkerRecordTest, RepeatInsideARecordIsWalkedThrough) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t s = 302; ; ++s) {
        b = clean_blocks({ 40, 30, 40, 40 }, s);
        if (b[2][0] != b[3][0])
            break;
    }
    const std::string &L = b[0], &R = b[1], &M = b[2], &N = b[3];
    const std::string record = L + R + M + R + N;
    const std::string seed = L.substr(5, 25);
    auto anno = build_records<Graph, Annotation>(kK, { { "F", "rec", record } });

    for (size_t branches : { size_t(0), Strategy::kUnlimited }) {
        auto res = run(*anno, seed, { "rec" }, strategy(branches));
        const ArmResult &left = res.arms[kLeft], &right = res.arms[kRight];
        check_invariants(left);
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(record.substr(30), spell_path(right, right.paths[0]));
        EXPECT_TRUE(right.splits.empty());
        EXPECT_EQ(0u, right.paths[0].end_labels[0].branches);
        for (EndReason never : { EndReason::BRANCH, EndReason::EDGE_REUSE,
                                 EndReason::EDGE_REUSE_RC, EndReason::REACHED_SEED,
                                 EndReason::LABEL_LOST }) {
            EXPECT_EQ(0u, count_ends(right, never)) << to_string(never);
        }
        // the record ends in a tip of the graph on both sides
        EXPECT_EQ(EndReason::DEAD_END, only_end(right, 0).reason);
        ASSERT_EQ(1u, left.paths.size());
        EXPECT_EQ(L.substr(0, 5), spell_path(left, left.paths[0]));
        EXPECT_EQ(EndReason::DEAD_END, only_end(left, 0).reason);
    }

    // the extension length still cuts the walk, inside the second copy
    Strategy st;
    st.max_extension_bp = 10 + R.size() + M.size() + 15;
    auto cut = run(*anno, seed, { "rec" }, st);
    const ArmResult &right = cut.arms[kRight];
    ASSERT_EQ(1u, right.paths.size());
    EXPECT_EQ(record.substr(30, st.max_extension_bp), spell_path(right, right.paths[0]));
    EXPECT_EQ(EndReason::MAX_EXTENSION, only_end(right, 0).reason);
    EXPECT_EQ(ArmResult::COMPLETE, right.status);
    EXPECT_EQ(st.max_extension_bp, right.complete_to_bp);
}

// A record holding the seed twice is followed from the lower place. The record is the
// second of its column, so its coordinates do not start at 0.
TYPED_TEST(WalkerRecordTest, ASeedHeldTwiceIsPlacedAtTheLowestCoordinate) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 50, 30, 40, 40 }, 303);
    const std::string &X = b[0], &S = b[1], &A = b[2], &B = b[3];
    const std::string record = S + A + S + B;
    auto anno = build_records<Graph, Annotation>(kK, { { "F", "other", X },
                                                      { "F", "rec", record } });

    for (const std::vector<std::string> &labels : { std::vector<std::string>{ "rec" },
                                                    std::vector<std::string>{} }) {
        auto res = run(*anno, S, labels);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ("rec", res.label_dict[0].name);
        EXPECT_EQ(1u, res.label_dict[0].seq_id);
        const ArmResult &left = res.arms[kLeft], &right = res.arms[kRight];
        check_invariants(left);
        check_invariants(right);
        // from the first copy: nothing to the left, the whole rest to the right
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(A + S + B, spell_path(right, right.paths[0]));
        EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
        EXPECT_EQ(0u, count_ends(right, EndReason::REACHED_SEED));
        EXPECT_EQ(EndReason::DEAD_END, only_end(right, 0).reason);
        ASSERT_EQ(1u, left.paths.size());
        EXPECT_EQ(0u, left.paths[0].length_bp);
        // the record's first k-mer has a predecessor in the graph (the second copy's)
        EXPECT_EQ(EndReason::RECORD_END, only_end(left, 0).reason);
        EXPECT_EQ(0u, only_end(left, 0).at_bp);
    }

    // A record with every k-mer of the seed, but not at consecutive coordinates (the
    // seed's two overlapping halves side by side), does not carry it: dropped, with the
    // runs of k-mers it is on
    const std::string halves = S.substr(0, 25) + S.substr(5);
    auto gapped = build_records<Graph, Annotation>(kK, { { "F", "whole", S + B },
                                                        { "F", "gapped", halves } });
    auto res = run(*gapped, S, { "whole", "gapped" });
    ASSERT_EQ(1u, res.label_dict.size());
    EXPECT_EQ("whole", res.label_dict[0].name);
    ASSERT_EQ(1u, res.dropped_labels.size());
    EXPECT_EQ("gapped", res.dropped_labels[0].name);
    EXPECT_EQ((std::vector<KmerInterval>{ { 0, S.size() - kK + 1 } }),
              res.dropped_labels[0].runs);
    auto derived = run(*gapped, S, {});
    ASSERT_EQ(1u, derived.label_dict.size());
    EXPECT_EQ("whole", derived.label_dict[0].name);
    try {
        run(*gapped, S, { "gapped" });
        ADD_FAILURE() << "a record not holding the seed in a row must not carry it";
    } catch (const SeedDerivationError &e) {
        EXPECT_EQ(SeedDerivationError::NO_CARRIER, e.cause());
    }
}

// A column and one of its records on one seed: the column keeps its membership walk and
// the walk rule, the record goes on alone
TYPED_TEST(WalkerRecordTest, ColumnAndHeaderLabelsOnOneSeed) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t s = 304; ; ++s) {
        b = clean_blocks({ 40, 30, 40, 40 }, s);
        if (b[2][0] != b[3][0])
            break;
    }
    const std::string &L = b[0], &R = b[1], &M = b[2], &N = b[3];
    const std::string record = L + R + M + R + N;
    const std::string seed = L.substr(5, 25);
    const uint64_t fork_bp = 10 + R.size();
    const uint64_t end_bp = record.size() - 30;
    auto anno = build_records<Graph, Annotation>(kK, { { "F", "rec", record } });

    // allowance 0: the column ends at the fork after the first copy, as it does alone
    auto res = run(*anno, seed, { "F", "rec" });
    ASSERT_EQ(2u, res.label_dict.size());
    EXPECT_EQ(LabelKind::COLUMN, res.label_dict[0].kind);
    EXPECT_EQ(LabelKind::HEADER, res.label_dict[1].kind);
    {
        const ArmResult &right = res.arms[kRight];
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(record.substr(30), spell_path(right, right.paths[0]));
        EXPECT_EQ(EndReason::BRANCH, only_end(right, 0).reason);
        EXPECT_EQ(fork_bp, only_end(right, 0).at_bp);
        EXPECT_EQ(EndReason::DEAD_END, only_end(right, 1).reason);
        EXPECT_EQ(end_bp, only_end(right, 1).at_bp);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), right.segments[0].labels_start);
        EXPECT_EQ((std::vector<LabelId>{ 1 }), right.segments[0].labels_end);

        auto alone = run(*anno, seed, { "F" });
        const ArmResult &column = alone.arms[kRight];
        ASSERT_EQ(1u, column.paths.size());
        EXPECT_EQ(EndReason::BRANCH, only_end(column, 0).reason);
        EXPECT_EQ(fork_bp, only_end(column, 0).at_bp);
        EXPECT_EQ(record.substr(30, fork_bp), spell_path(column, column.paths[0]));
    }

    // unlimited: the column follows both successors of the first copy; on the way back
    // into the repeat the walk rule stops it, and only it
    res = run(*anno, seed, { "F", "rec" }, strategy(Strategy::kUnlimited));
    const ArmResult &right = res.arms[kRight];
    check_invariants(right);
    ASSERT_EQ(2u, right.paths.size());
    std::set<std::string> flanks;
    for (const auto &path : right.paths) {
        flanks.insert(spell_path(right, path));
    }
    EXPECT_EQ((std::set<std::string>{ record.substr(30), L.substr(30) + R + N }), flanks);
    EXPECT_EQ(1u, count_ends(right, EndReason::EDGE_REUSE));
    EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
    size_t header_ends = 0;
    for (const auto &seg : right.segments) {
        for (const auto &e : seg.label_ends) {
            if (e.label == 1) {
                header_ends++;
                EXPECT_EQ(EndReason::DEAD_END, e.reason);
                EXPECT_EQ(end_bp, e.at_bp);
            } else {
                EXPECT_TRUE(e.reason == EndReason::EDGE_REUSE || e.reason == EndReason::DEAD_END)
                    << to_string(e.reason);
            }
        }
    }
    EXPECT_EQ(1u, header_ends);
    // the column's walks are those of the column alone
    auto alone = run(*anno, seed, { "F" }, strategy(Strategy::kUnlimited));
    const ArmResult &column = alone.arms[kRight];
    EXPECT_EQ(1u, count_ends(column, EndReason::EDGE_REUSE));
    EXPECT_EQ(1u, count_ends(column, EndReason::DEAD_END));
    ASSERT_EQ(2u, column.paths.size());
}

// A record with a base outside the alphabet has no k-mers over it, and so a gap in its
// coordinates. Here the gap is annotated directly: the record's two parts are given
// coordinates kK apart, as around one such base. The label is lost at the gap (the last
// k-mer before it has a successor, of another record), not at a record end.
TYPED_TEST(WalkerRecordTest, AGapInTheRecordsCoordinatesLosesTheLabel) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 60, 40, 30 }, 306);
    const std::string &X = b[0], &Y = b[1], &W = b[2];
    const std::string other = X.substr(40) + W;
    const uint64_t nx = X.size() - kK + 1, ny = Y.size() - kK + 1;
    auto anno = build_anno_graph<Graph, Annotation>(
            kK, { X, Y, other }, { "F", "F", "G" }, DeBruijnGraph::BASIC, true,
            { 0, nx + kK, 0 });
    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[encoder.encode("F")] = { "rec" };
    num_kmers[encoder.encode("F")] = { nx + kK + ny };
    headers[encoder.encode("G")] = { "other" };
    num_kmers[encoder.encode("G")] = { other.size() - kK + 1 };
    anno = with_coord_to_header(std::move(anno), std::move(headers), std::move(num_kmers));

    const std::string seed = X.substr(10, 25);
    for (const std::vector<std::string> &labels : { std::vector<std::string>{ "rec" },
                                                    std::vector<std::string>{} }) {
        auto res = run(*anno, seed, labels);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ("rec", res.label_dict[0].name);
        const ArmResult &left = res.arms[kLeft], &right = res.arms[kRight];
        check_invariants(left);
        check_invariants(right);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(X.substr(35), spell_path(right, right.paths[0]));
        EXPECT_EQ(EndReason::LABEL_LOST, only_end(right, 0).reason)
            << to_string(only_end(right, 0).reason);
        EXPECT_EQ(X.size() - 35, only_end(right, 0).at_bp);
        EXPECT_EQ(0u, count_ends(right, EndReason::RECORD_END));
        ASSERT_EQ(1u, left.paths.size());
        EXPECT_EQ(X.substr(0, 10), spell_path(left, left.paths[0]));
        EXPECT_EQ(EndReason::DEAD_END, only_end(left, 0).reason);
    }
    // after the gap the record goes on at its own coordinates
    auto after = run(*anno, Y.substr(5, 25), { "rec" });
    ASSERT_EQ(1u, after.arms[kLeft].paths.size());
    EXPECT_EQ(Y.substr(0, 5), spell_path(after.arms[kLeft], after.arms[kLeft].paths[0]));
    // the record's first k-mer after the gap has no predecessor in the graph
    EXPECT_EQ(EndReason::DEAD_END, only_end(after.arms[kLeft], 0).reason);
}

#if ! _PROTEIN_GRAPH
// In canonical and primary graphs a seed that is the reverse complement of a part of the
// record is held at descending coordinates: the arms follow the record backwards
template <class Annotation>
void check_reverse_complement_seed() {
    auto b = clean_blocks({ 30, 40, 40 }, 305);
    const std::string &T = b[0], &A = b[1], &B = b[2];
    const std::string record = T + A + B + T;
    const std::string forward = A.substr(5, 25);
    for (auto mode : canonical_modes()) {
        auto anno = build_records<DBGSuccinct, Annotation>(kK, { { "F", "rec", record } }, mode);
        for (const std::vector<std::string> &labels : { std::vector<std::string>{ "rec" },
                                                        std::vector<std::string>{} }) {
            auto fwd = run(*anno, forward, labels);
            ASSERT_EQ(1u, fwd.label_dict.size()) << mode;
            ASSERT_EQ(1u, fwd.arms[kRight].paths.size()) << mode;
            ASSERT_EQ(1u, fwd.arms[kLeft].paths.size()) << mode;
            EXPECT_EQ(A.substr(30) + B + T,
                      spell_path(fwd.arms[kRight], fwd.arms[kRight].paths[0])) << mode;
            EXPECT_EQ(T + A.substr(0, 5),
                      spell_path(fwd.arms[kLeft], fwd.arms[kLeft].paths[0])) << mode;

            auto bwd = run(*anno, rc(forward), labels);
            ASSERT_EQ(1u, bwd.label_dict.size()) << mode;
            EXPECT_EQ(LabelKind::HEADER, bwd.label_dict[0].kind);
            ASSERT_EQ(1u, bwd.arms[kRight].paths.size()) << mode;
            ASSERT_EQ(1u, bwd.arms[kLeft].paths.size()) << mode;
            // the right arm of the reverse complement is the record's left side
            EXPECT_EQ(rc(T + A.substr(0, 5)),
                      spell_path(bwd.arms[kRight], bwd.arms[kRight].paths[0])) << mode;
            EXPECT_EQ(rc(A.substr(30) + B + T),
                      spell_path(bwd.arms[kLeft], bwd.arms[kLeft].paths[0])) << mode;
            for (const SeedResult *res : { &fwd, &bwd }) {
                for (const ArmResult &arm : res->arms) {
                    check_invariants(arm);
                    EXPECT_EQ(EndReason::RECORD_END, only_end(arm, 0).reason) << mode;
                }
            }
        }
    }
}

TEST(WalkerRecord, ReverseComplementSeedFollowsTheRecordBackwards) {
    check_reverse_complement_seed<annot::ColumnCompressed<>>();
}

// the same on a row-diff coordinate annotation, whose rows are decoded through the cache
TEST(WalkerRecord, ReverseComplementSeedOnRowDiff) {
    check_reverse_complement_seed<annot::RowDiffColumnAnnotator>();
}

// A seed and its reverse complement are one seed in canonical and primary graphs (one
// seed id), so they take the same place in a record and their walks are mirror images.
// Here the record holds rc(S) at its coordinates c..c+19 and S, overlapping it, at
// c+10..c+29: both seeds are placed on the lower occurrence.
TEST(WalkerRecord, ASeedAndItsReverseComplementTakeTheSamePlace) {
    auto b = clean_blocks({ 40, 40, 10, 10 }, 307);
    const std::string &L = b[0], &R = b[1], &H = b[2], &T = b[3];
    const std::string S = H + rc(H) + T;
    // rc(S) = rc(T) + H + rc(H): the two occurrences share H + rc(H)
    const std::string record = L + rc(T) + H + rc(H) + T + R;
    ASSERT_EQ(L.size(), record.find(rc(S)));
    ASSERT_EQ(L.size() + 10, record.find(S));
    for (auto mode : canonical_modes()) {
        auto anno = build_records<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, { { "F", "rec", record } }, mode);
        for (const std::vector<std::string> &labels : { std::vector<std::string>{ "rec" },
                                                        std::vector<std::string>{} }) {
            auto fwd = run(*anno, S, labels);
            auto bwd = run(*anno, rc(S), labels);
            ASSERT_EQ(1u, fwd.label_dict.size()) << mode;
            ASSERT_EQ(1u, bwd.label_dict.size()) << mode;
            for (const SeedResult *res : { &fwd, &bwd }) {
                for (const ArmResult &arm : res->arms) {
                    check_invariants(arm);
                    ASSERT_EQ(1u, arm.paths.size()) << mode;
                }
            }
            // rc(S) lies forwards at the lower occurrence: the record's own flanks
            EXPECT_EQ(L, spell_path(bwd.arms[kLeft], bwd.arms[kLeft].paths[0])) << mode;
            EXPECT_EQ(T + R, spell_path(bwd.arms[kRight], bwd.arms[kRight].paths[0])) << mode;
            // S is placed there too, backwards: the mirror image
            EXPECT_EQ(rc(spell_path(bwd.arms[kLeft], bwd.arms[kLeft].paths[0])),
                      spell_path(fwd.arms[kRight], fwd.arms[kRight].paths[0])) << mode;
            EXPECT_EQ(rc(spell_path(bwd.arms[kRight], bwd.arms[kRight].paths[0])),
                      spell_path(fwd.arms[kLeft], fwd.arms[kLeft].paths[0])) << mode;
            EXPECT_EQ(only_end(bwd.arms[kLeft], 0).reason, only_end(fwd.arms[kRight], 0).reason);
            EXPECT_EQ(only_end(bwd.arms[kRight], 0).reason, only_end(fwd.arms[kLeft], 0).reason);
            EXPECT_EQ(only_end(bwd.arms[kLeft], 0).at_bp, only_end(fwd.arms[kRight], 0).at_bp);
            EXPECT_EQ(only_end(bwd.arms[kRight], 0).at_bp, only_end(fwd.arms[kLeft], 0).at_bp);
        }
    }
}
#endif

// zero caps, a radius of 0 and a non-finite or negative time budget are refused
TEST(Walker, RefusesInvalidCaps) {
    const std::string s = random_seq(60, 17);
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { s }, { "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Seed seed;
    seed.sequence = s.substr(10, 30);
    EXPECT_NO_THROW(traverse_seed(oracle, seed, Strategy()));
    const std::vector<std::function<void(Strategy &)>> invalid {
        [](Strategy &st) { st.max_seed_labels = 0; },
        [](Strategy &st) { st.max_steps = 0; },
        [](Strategy &st) { st.max_live_paths = 0; },
        [](Strategy &st) { st.max_paths = 0; },
        [](Strategy &st) { st.max_extension_bp = 0; },
        [](Strategy &st) { st.time_budget_ms = std::numeric_limits<double>::infinity(); },
        [](Strategy &st) { st.time_budget_ms = std::nan(""); },
        [](Strategy &st) { st.time_budget_ms = -1; },
    };
    for (size_t i = 0; i < invalid.size(); ++i) {
        Strategy st;
        invalid[i](st);
        EXPECT_THROW(traverse_seed(oracle, seed, st), std::invalid_argument) << "case " << i;
    }
}

// a derived header that an explicit list would resolve to a column is refused
TEST(WalkerDerive, DerivedHeaderMustNotBeAColumnName) {
    auto b = clean_blocks({ 40, 40, 30 }, 941);
    const std::string &D = b[0], &M = b[1], &T = b[2];
    const std::string in_acc = D, in_b = M + T;   // the seed M is only in column B's record
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { in_acc, in_b }, { "ACC", "B" }, DeBruijnGraph::BASIC, true, { 0, 0 });
    const auto &enc = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[enc.encode("ACC")] = { "decoy" };
    num_kmers[enc.encode("ACC")].assign(1, in_acc.size() - kK + 1);
    headers[enc.encode("B")] = { "ACC" };
    num_kmers[enc.encode("B")].assign(1, in_b.size() - kK + 1);
    anno = with_coord_to_header(std::move(anno), std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno);

    Seed seed;
    seed.sequence = M;
    Strategy st;
    st.direction = Strategy::RIGHT;
    try {
        traverse_seed(oracle, seed, st);
        FAIL() << "a derived header that an explicit list resolves to a column was accepted";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("'ACC'")) << e.what();
        EXPECT_NE(std::string::npos, std::string(e.what()).find("name of an annotation column")) << e.what();
        EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
    }
    st.seed_label_kind = LabelKind::COLUMN;
    auto by_column = traverse_seed(oracle, seed, st);
    ASSERT_EQ(1u, by_column.label_dict.size());
    EXPECT_EQ("B", by_column.label_dict[0].name);
    EXPECT_EQ(T.size(), by_column.arms[kRight].paths[0].length_bp);
}

} // namespace
