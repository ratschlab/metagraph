#include "gtest/gtest.h"

#include <algorithm>
#include <cstring>
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
#include "annotation/coord_to_header.hpp"


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
    EXPECT_EQ(arm.output_bp, bases);
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
    os << r.validated_seed_id << '|' << r.label_dict.size() << '|' << r.length_bp << '|';
    for (const auto &d : r.dropped_labels) {
        os << d.name << ';';
    }
    for (const ArmResult &a : r.arms) {
        os << "\nARM " << to_string(a.arm) << ' ' << a.requested << ' ' << a.status << ' '
           << a.steps << ' ' << a.output_bp << ' ' << a.frontier_live_paths << ' '
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
        EXPECT_EQ(16u, res.validated_seed_id.size());
        EXPECT_FALSE(res.seed_id_mismatch);
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

        // a lower-case seed maps the same way and gets the same seed id
        std::string lower = seed;
        for (char &c : lower) {
            c = std::tolower(c);
        }
        auto lc = run(*anno, lower, { "A" });
        EXPECT_EQ(serialize(res), serialize(lc));
        EXPECT_EQ(make_seed_id(seed, mode != DeBruijnGraph::BASIC, { "A" }), lc.validated_seed_id);
        LabelOracle oracle(*anno);
        Seed s;
        s.sequence = lower;
        s.labels = { "A" };
        s.seed_id = res.validated_seed_id;
        EXPECT_FALSE(traverse_seed(oracle, s, Strategy()).seed_id_mismatch);
        s.seed_id = "0000000000000000";
        auto mismatch = traverse_seed(oracle, s, Strategy());
        EXPECT_TRUE(mismatch.seed_id_mismatch);
        EXPECT_EQ("0000000000000000", mismatch.seed_id);
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
            ASSERT_EQ(1u, arm.paths.size());
            EXPECT_EQ(EndReason::MAX_STEPS, *arm.paths[0].path_reason);
            EXPECT_EQ(1u, count_ends(arm, EndReason::MAX_STEPS));
        }

        st = Strategy();
        st.max_output_bp = 10;
        res = run(*linear, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm);
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
            EXPECT_EQ(10u, arm.output_bp);
            ASSERT_TRUE(arm.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_OUTPUT, arm.cap_trigger->reason);
            EXPECT_EQ(10u, arm.cap_trigger->at_bp);
            EXPECT_EQ(1u, count_ends(arm, EndReason::MAX_OUTPUT));
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

// the split allowance of a path: past it the labels end with BRANCH, "split_limit"
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
        EXPECT_EQ(2u, path_ends(right.paths[0], EndReason::BRANCH));
        EXPECT_TRUE(right.splits.empty());
        EXPECT_EQ(0u, right.steps);
        for (LabelId l : { 0u, 1u }) {
            const LabelEnd *e = label_end(right, l);
            ASSERT_NE(nullptr, e);
            EXPECT_EQ(0u, e->at_bp);
            EXPECT_EQ(EndReason::BRANCH, e->reason);
            EXPECT_STREQ("split_limit", e->text);
        }
        st.max_splits_per_path = 1;
        EXPECT_EQ(2u, run(*anno, X, { "A", "B" }, st).arms[kRight].paths.size());
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
        EXPECT_EQ((std::vector<std::pair<uint64_t, uint64_t>>{ { 0, supported } }),
                  res.dropped_labels[0].runs);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ("B", res.label_dict[0].name);
        EXPECT_EQ(1u, res.labels_supporting_total);
        EXPECT_EQ(make_seed_id(seed, mode != DeBruijnGraph::BASIC, { "B" }), res.validated_seed_id);
        EXPECT_EQ(f.Y.substr(5), spell_path(res.arms[kRight], res.arms[kRight].paths[0]));

        EXPECT_THROW(run(*anno, f.W + f.R, { "B" }), std::invalid_argument);     // unsupported
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
        EXPECT_TRUE(derived.labels_dropped_digest.empty());
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
        const std::string &cut = all_three.label_dict[2].name;
        EXPECT_EQ(hex64(fnv1a64(cut + "\n")), capped.labels_dropped_digest) << mode;
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
        EXPECT_EQ(a.validated_seed_id, c.validated_seed_id);
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
TEST(WalkerDerive, IntersectionStartsFromTheCheapestRowOfTheFirstWindow) {
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
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno, &cth);
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
    // the decoys' coordinates of U were never mapped to their sequences
    EXPECT_LT(derived.annotation_counters.coords_mapped, kDecoys);
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
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno, &cth);
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
    seed.labels = { "ACC1" };
    auto one = traverse_seed(oracle, seed, Strategy());
    EXPECT_EQ(1u, one.label_dict.size());
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
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);

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
