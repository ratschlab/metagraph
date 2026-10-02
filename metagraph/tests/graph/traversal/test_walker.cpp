#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

#include "tests/test_helpers.hpp"
#include "tests/graph/all/test_dbg_helpers.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/traversal/resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "cli/traverse.hpp"


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

// no repeated k-mer or (k-1)-mer (in either orientation; the graph is node-centric,
// so a repeated (k-1)-mer implies an edge), no RC-palindromic (k-1)-mer or (k+1)-mer
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

std::string clean_block(size_t len, uint32_t seed) {
    for (uint32_t s = seed; ; ++s) {
        std::string b = random_seq(len, s);
        if (is_clean(b))
            return b;
    }
}

// split a clean master sequence into consecutive blocks of the given lengths
std::vector<std::string> clean_blocks(const std::vector<size_t> &lengths, uint32_t seed) {
    size_t total = 0;
    for (size_t l : lengths) total += l;
    std::string master = clean_block(total, seed);
    std::vector<std::string> blocks;
    size_t pos = 0;
    for (size_t l : lengths) {
        blocks.push_back(master.substr(pos, l));
        pos += l;
    }
    return blocks;
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

Strategy strategy(size_t max_label_branches = 0, bool merge = true) {
    Strategy st;
    st.max_label_branches = max_label_branches;
    st.merge_reconverge = merge;
    return st;
}

SeedResult run(const AnnotatedDBG &anno,
               const std::string &seq,
               const std::vector<std::string> &labels,
               const Strategy &st = Strategy(),
               const LabelChangeCost &cost = LabelChangeCost::forbid()) {
    LabelOracle oracle(anno);
    Seed seed;
    seed.sequence = seq;
    seed.labels = labels;
    return traverse_seed(oracle, seed, st, cost);
}

// the flank plus the seed must map to the graph without gaps
void check_spelled(const AnnotatedDBG &anno, const ArmResult &arm, const PathResult &path,
                   const std::string &seed) {
    std::string flank = spell_path(arm, path);
    EXPECT_EQ(path.length_bp, flank.size());
    std::string full = arm.arm == Arm::RIGHT ? seed + flank : flank + seed;
    for (node_index n : map_to_nodes_sequentially(anno.get_graph(), full)) {
        ASSERT_NE(npos, n) << full;
    }
}

uint64_t sum_steps(const ArmResult &arm) {
    uint64_t s = 0;
    for (const auto &b : arm.growth) s += b.steps;
    return s;
}

size_t count_events(const ArmResult &arm, EventType type) {
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &ev : seg.events) n += ev.type == type;
    }
    return n;
}

size_t count_ends(const ArmResult &arm, EndReason reason) {
    size_t n = 0;
    for (const auto &run : arm.runs) n += run.ended && run.end_reason == reason;
    return n;
}

std::vector<const Event*> events_of(const ArmResult &arm, EventType type) {
    std::vector<const Event*> out;
    for (const auto &seg : arm.segments) {
        for (const auto &ev : seg.events) {
            if (ev.type == type)
                out.push_back(&ev);
        }
    }
    return out;
}

// the LABEL_END event of |label| (nullptr if the label never ended)
const Event* label_end_event(const ArmResult &arm, LabelId label) {
    for (const Event *ev : events_of(arm, EventType::LABEL_END)) {
        if (ev->label == label)
            return ev;
    }
    return nullptr;
}

// does |label_seq| contain the k-mer (in either orientation for canonical modes)?
bool contains_kmer(const std::string &label_seq, const std::string &kmer, bool canonical) {
    return label_seq.find(kmer) != std::string::npos
        || (canonical && label_seq.find(rc(kmer)) != std::string::npos);
}

// does |label_seq| contain every k-mer of |seq|?
bool supports_sequence(const std::string &label_seq, const std::string &seq, bool canonical) {
    for (size_t i = 0; i + kK <= seq.size(); ++i) {
        if (!contains_kmer(label_seq, seq.substr(i, kK), canonical))
            return false;
    }
    return true;
}

// The flank of |label| ending at segment |leaf|, in natural orientation: at a join
// the parent through which the label entered (Segment::labels_via_parent) is taken.
std::string route_of(const ArmResult &arm, size_t leaf, LabelId label) {
    std::vector<size_t> route;
    size_t s = leaf;
    while (true) {
        route.push_back(s);
        const Segment &seg = arm.segments[s];
        if (seg.parents.empty())
            break;
        size_t next = seg.parents[0];
        if (seg.parents.size() > 1) {
            EXPECT_EQ(seg.parents.size(), seg.labels_via_parent.size());
            bool found = false;
            for (size_t p = 0; p < seg.labels_via_parent.size(); ++p) {
                const auto &via = seg.labels_via_parent[p];
                if (std::find(via.begin(), via.end(), label) != via.end()) {
                    next = seg.parents[p];
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "label " << label << " enters segment " << s << " through no parent";
        }
        s = next;
    }
    std::string out;
    if (arm.arm == Arm::RIGHT) {
        for (auto it = route.rbegin(); it != route.rend(); ++it) out += arm.segments[*it].sequence;
    } else {
        for (size_t x : route) out += arm.segments[x].sequence;
    }
    return out;
}

// blocks {X, P, Q} with P[0] != Q[0]: X+P and X+Q fork right after the seed X
std::vector<std::string> fork_blocks(uint32_t seed, std::vector<size_t> lengths = { 30, 40, 40 }) {
    for (uint32_t s = seed; ; ++s) {
        auto b = clean_blocks(lengths, s);
        if (b[1][0] != b[2][0])
            return b;
    }
}

// §7.2 invariants
void check_invariants(const ArmResult &arm, const Strategy &st) {
    EXPECT_EQ(arm.steps, sum_steps(arm));
    std::array<uint32_t, kNumEndReasons> by_bins {};
    size_t ambiguous = 0;
    for (const auto &b : arm.growth) {
        for (size_t r = 0; r < kNumEndReasons; ++r) by_bins[r] += b.label_ends[r];
        ambiguous += b.ambiguous_branches;
        EXPECT_LE(b.max_live_paths, st.max_live_paths);
    }
    for (size_t r = 0; r < kNumEndReasons; ++r) {
        EXPECT_EQ(count_ends(arm, static_cast<EndReason>(r)), by_bins[r])
            << to_string(static_cast<EndReason>(r));
    }
    size_t ambiguous_splits = 0;
    for (const auto &sp : arm.splits) ambiguous_splits += sp.ambiguous;
    EXPECT_EQ(ambiguous_splits, ambiguous);
    // every leaf reports an end for each of its final labels
    for (const auto &path : arm.paths) {
        uint32_t ends = 0;
        for (uint32_t x : path.end_reasons) ends += x;
        EXPECT_EQ(path.end_labels.size(), ends);
        for (const auto &e : path.end_labels) {
            EXPECT_TRUE(arm.runs[e.run].ended);
            EXPECT_EQ(path.length_bp, arm.runs[e.run].to_bp);
        }
    }
}

// the essential result fields, as a string (no timing, no cache counters)
std::string serialize(const SeedResult &r) {
    std::ostringstream os;
    os << r.validated_seed_id << '|' << r.num_seed_labels << '|' << r.length_bp << '|';
    for (const auto &d : r.dropped_labels) os << d.name << ':' << d.reason << ';';
    for (const ArmResult &a : r.arms) {
        os << "\nARM " << to_string(a.arm) << ' ' << a.requested << ' ' << a.status << ' '
           << a.steps << ' ' << a.output_bp << ' ' << a.successor_enumerations << ' '
           << a.branch_events_total << ' ' << a.frontier_live_paths << ' ' << a.frontier_live_labels;
        if (a.cap_trigger) {
            os << " cap " << to_string(a.cap_trigger->reason) << ' ' << a.cap_trigger->at_bp
               << ' ' << a.cap_trigger->segment << ' ' << a.cap_trigger->live_paths
               << ' ' << a.cap_trigger->live_labels;
        }
        for (const Segment &s : a.segments) {
            os << "\n S" << s.id << " p";
            for (size_t p : s.parents) os << p << ',';
            os << " c";
            for (size_t c : s.children) os << c << ',';
            os << ' ' << s.from_bp << ' ' << s.length_bp << ' ' << s.sequence << " ls";
            for (LabelId l : s.labels_start) os << l << ',';
            os << " le";
            for (LabelId l : s.labels_end) os << l << ',';
            for (const Event &ev : s.events) {
                os << "\n  E" << ev.at_bp << ' ' << to_string(ev.type) << ' ' << ev.label << ' '
                   << ev.to << ' ' << ev.cost << ' ' << ev.needed_budget << ' '
                   << to_string(ev.reason) << ' ' << ev.ch << ' ' << ev.length_bp << ' '
                   << ev.structural_successors << ' ';
                for (LabelId l : ev.labels) os << l << ',';
                os << ' ' << ev.text;
            }
        }
        for (const Split &sp : a.splits) {
            os << "\n X" << sp.at_bp << ' ' << sp.segment << ' ' << sp.ambiguous << ' ';
            for (size_t c : sp.children) os << c << ',';
        }
        for (const PathResult &p : a.paths) {
            os << "\n P" << p.id << ' ' << p.length_bp << ' ';
            for (size_t s : path_segments(a, p)) os << s << ',';
            os << ' ';
            for (uint32_t x : p.end_reasons) os << x << ',';
            os << ' ' << (p.path_reason ? to_string(*p.path_reason) : "-") << ' ';
            for (const auto &e : p.end_labels) os << e.label << ':' << e.loss << ':' << e.branches << ':' << e.run << ',';
            if (p.continuation) {
                os << " cont " << p.continuation->sequence << ' ' << p.continuation->loss_used
                   << ' ' << p.continuation->branches_used << ' ';
                for (LabelId l : p.continuation->labels) os << l << ',';
            }
        }
        for (const LabelRun &run : a.runs) {
            os << "\n R" << run.label << ' ' << run.from_bp << ' ' << run.to_bp << ' '
               << run.entered_by_switch << ' ' << run.from_label << ' ' << run.switch_cost
               << ' ' << run.prev_run << ' ' << run.ended << ' ' << to_string(run.end_reason);
        }
        for (const GrowthBin &b : a.growth) {
            os << "\n G" << b.from_bp << ' ' << b.max_live_paths << ' ' << b.max_live_labels << ' '
               << b.max_live_pairs << ' ' << b.steps << ' ' << b.divergences << ' '
               << b.ambiguous_branches << ' ' << b.splits << ' ' << b.reconvergences << ' '
               << b.blocked_repeat << ' ';
            for (uint32_t x : b.label_ends) os << x << ',';
        }
        for (const BranchEvent &be : a.branch_events) {
            os << "\n B" << be.at_bp << ' ' << be.segment << ' ';
            for (char c : be.chars) os << c;
            os << ' ';
            for (size_t n : be.labels_per_successor) os << n << ',';
            os << ' ';
            for (LabelId l : be.ambiguous) os << l << ',';
            os << ' ';
            for (LabelId l : be.dropped) os << l << ',';
            os << ' ' << be.labels_affected;
        }
        for (double d : a.needed_budgets) os << " nb" << d;
    }
    for (size_t l = 0; l < r.label_summary.size(); ++l) {
        for (size_t a = 0; a < 2; ++a) {
            const auto &s = r.label_summary[l][a];
            os << "\n L" << l << ' ' << a << ' ' << s.direct_bp << ' ' << s.reach_bp << ' '
               << s.reentries << ' ';
            for (uint32_t x : s.runs) os << x << ',';
        }
    }
    return os.str();
}

// the per-seed outcome (spec §7.0) as "walks/branch_diagnostics/label_evidence/delivery",
// so that an assertion names all four axes at once
std::string outcome_of(const Json::Value &result) {
    const Json::Value &o = result["outcome"];
    EXPECT_TRUE(o.isObject());
    EXPECT_EQ(4u, o.size());
    return o["walks"].asString() + "/" + o["branch_diagnostics"].asString() + "/"
         + o["label_evidence"].asString() + "/" + o["delivery"].asString();
}


template <typename Pair>
class WalkerTest : public ::testing::Test {};

typedef ::testing::Types<
    std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
    std::pair<DBGHashFast, annot::ColumnCompressed<>>,
    std::pair<DBGSSHash, annot::ColumnCompressed<>>
> WalkerTypes;
TYPED_TEST_SUITE(WalkerTest, WalkerTypes);


// T5: a linear labeled sequence is reconstructed exactly on both sides of the seed
TYPED_TEST(WalkerTest, Linear) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(200, 1);
    std::string seed = S.substr(60, 40);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        Strategy st;
        auto res = run(*anno, seed, { "A" }, st);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ(1u, res.num_seed_labels);
        EXPECT_TRUE(res.dropped_labels.empty());
        EXPECT_EQ(16u, res.validated_seed_id.size());
        EXPECT_FALSE(res.seed_id_mismatch);
        EXPECT_EQ(40u, res.length_bp);
        EXPECT_EQ(30u, res.num_kmers);

        const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
        for (const ArmResult *arm : { &left, &right }) {
            ASSERT_TRUE(arm->requested);
            EXPECT_EQ(ArmResult::COMPLETE, arm->status);
            EXPECT_FALSE(arm->cap_trigger.has_value());
            ASSERT_EQ(1u, arm->paths.size());
            ASSERT_EQ(1u, arm->segments.size());
            ASSERT_EQ(1u, arm->runs.size());
            const PathResult &path = arm->paths[0];
            EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
            EXPECT_FALSE(path.path_reason.has_value());
            EXPECT_FALSE(path.continuation.has_value());
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(0u, path.end_labels[0].label);
            EXPECT_EQ(0.0, path.end_labels[0].loss);
            EXPECT_EQ(0u, path.end_labels[0].branches);
            EXPECT_TRUE(arm->runs[0].ended);
            EXPECT_EQ(EndReason::DEAD_END, arm->runs[0].end_reason);
            EXPECT_EQ(0u, arm->runs[0].from_bp);
            EXPECT_FALSE(arm->runs[0].entered_by_switch);
            EXPECT_EQ((std::vector<LabelId>{ 0 }), arm->segments[0].labels_start);
            EXPECT_EQ((std::vector<LabelId>{ 0 }), arm->segments[0].labels_end);
            EXPECT_TRUE(arm->splits.empty());
            EXPECT_EQ(0u, arm->branch_events_total);
            check_spelled(*anno, *arm, path, seed);
            check_invariants(*arm, st);
            EXPECT_EQ(arm->steps, path.length_bp);
            EXPECT_EQ(arm->output_bp, path.length_bp);
            EXPECT_GE(arm->successor_enumerations, path.length_bp);
        }
        std::string rflank = spell_path(right, right.paths[0]);
        std::string lflank = spell_path(left, left.paths[0]);
        EXPECT_EQ(S.substr(100), rflank) << mode;
        EXPECT_EQ(S.substr(0, 60), lflank) << mode;
        EXPECT_EQ(S, lflank + seed + rflank);
        EXPECT_EQ(100u, right.runs[0].to_bp);
        EXPECT_EQ(60u, left.runs[0].to_bp);
        EXPECT_EQ(100u, res.label_summary[0][kRight].direct_bp);
        EXPECT_EQ(100u, res.label_summary[0][kRight].reach_bp);
        EXPECT_EQ(60u, res.label_summary[0][kLeft].direct_bp);
        EXPECT_EQ(0u, res.label_summary[0][kLeft].reentries);
        EXPECT_STRNE("", res.access_path);

        // single-direction requests produce one arm only
        st.direction = Strategy::LEFT;
        auto only_left = run(*anno, seed, { "A" }, st);
        EXPECT_FALSE(only_left.arms[kRight].requested);
        EXPECT_TRUE(only_left.arms[kRight].paths.empty());
        EXPECT_TRUE(only_left.arms[kRight].segments.empty());
        EXPECT_EQ(0u, only_left.arms[kRight].steps);
        ASSERT_EQ(1u, only_left.arms[kLeft].paths.size());
        EXPECT_EQ(S.substr(0, 60), spell_path(only_left.arms[kLeft], only_left.arms[kLeft].paths[0]));

        st.direction = Strategy::RIGHT;
        auto only_right = run(*anno, seed, { "A" }, st);
        EXPECT_FALSE(only_right.arms[kLeft].requested);
        EXPECT_TRUE(only_right.arms[kLeft].paths.empty());
        ASSERT_EQ(1u, only_right.arms[kRight].paths.size());
        EXPECT_EQ(S.substr(100), spell_path(only_right.arms[kRight], only_right.arms[kRight].paths[0]));

        // a lower-case seed maps the same way, is not rewritten, and gets the same
        // seed id (ids are computed over the case-mapped sequence, spec §4.3/§6.1)
        std::string lower = seed;
        for (char &c : lower) c = std::tolower(c);
        auto lc = run(*anno, lower, { "A" }, Strategy());
        ASSERT_EQ(1u, lc.arms[kRight].paths.size());
        ASSERT_EQ(1u, lc.arms[kLeft].paths.size());
        check_invariants(lc.arms[kRight], Strategy());
        check_invariants(lc.arms[kLeft], Strategy());
        EXPECT_EQ(S.substr(100), spell_path(lc.arms[kRight], lc.arms[kRight].paths[0]));
        EXPECT_EQ(S.substr(0, 60), spell_path(lc.arms[kLeft], lc.arms[kLeft].paths[0]));
        EXPECT_EQ(res.validated_seed_id, lc.validated_seed_id);
        EXPECT_EQ(make_seed_id("", seed, mode != DeBruijnGraph::BASIC, { "A" }), lc.validated_seed_id);
        {
            LabelOracle oracle(*anno);
            Seed s;
            s.sequence = lower;
            s.labels = { "A" };
            s.seed_id = res.validated_seed_id;
            EXPECT_FALSE(traverse_seed(oracle, s, Strategy(), LabelChangeCost::forbid()).seed_id_mismatch);
        }
    }
}


// T6: exact label-end coordinates. A on W·R·X, B on R·X·Y, seed R.
struct LabelEndFixture {
    std::string W, R, X, Y;
    LabelEndFixture() {
        auto b = clean_blocks({ 20, 60, 25, 30 }, 2);
        W = b[0]; R = b[1]; X = b[2]; Y = b[3];
    }
    std::vector<std::string> seqs() const { return { W + R + X, R + X + Y }; }
    std::vector<std::string> labels() const { return { "A", "B" }; }
};

TYPED_TEST(WalkerTest, LabelEndCoordinates) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    LabelEndFixture f;
    for (auto mode : all_modes()) {
        bool canonical = mode != DeBruijnGraph::BASIC;
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        Strategy st;
        auto res = run(*anno, f.R, { "A", "B" }, st);
        ASSERT_EQ(2u, res.label_dict.size());
        const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
        check_invariants(right, st);
        check_invariants(left, st);

        // right arm: A ends label_lost at 25, B dead_end at 55
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(55u, right.paths[0].length_bp);
        EXPECT_EQ(f.X + f.Y, spell_path(right, right.paths[0]));
        ASSERT_EQ(2u, right.runs.size());
        EXPECT_EQ(0u, right.runs[0].label);
        EXPECT_EQ(25u, right.runs[0].to_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, right.runs[0].end_reason);
        EXPECT_EQ(1u, right.runs[1].label);
        EXPECT_EQ(55u, right.runs[1].to_bp);
        EXPECT_EQ(EndReason::DEAD_END, right.runs[1].end_reason);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), right.segments[0].labels_start);
        EXPECT_EQ((std::vector<LabelId>{ 1 }), right.segments[0].labels_end);
        ASSERT_EQ(2u, right.segments[0].events.size());
        const Event &lost = right.segments[0].events[0];
        EXPECT_EQ(25u, lost.at_bp);
        EXPECT_EQ(EventType::LABEL_END, lost.type);
        EXPECT_EQ(0u, lost.label);
        EXPECT_EQ(EndReason::LABEL_LOST, lost.reason);
        EXPECT_EQ(1u, lost.structural_successors);
        EXPECT_EQ(55u, right.segments[0].events[1].at_bp);
        EXPECT_EQ(1u, right.paths[0].end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
        EXPECT_EQ(0u, right.paths[0].end_reasons[static_cast<size_t>(EndReason::LABEL_LOST)]);

        // left arm: B label_lost at 0, A dead_end at 20
        ASSERT_EQ(1u, left.paths.size());
        EXPECT_EQ(20u, left.paths[0].length_bp);
        EXPECT_EQ(f.W, spell_path(left, left.paths[0]));
        ASSERT_EQ(2u, left.runs.size());
        EXPECT_EQ(20u, left.runs[0].to_bp);
        EXPECT_EQ(EndReason::DEAD_END, left.runs[0].end_reason);
        EXPECT_EQ(0u, left.runs[1].to_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, left.runs[1].end_reason);
        EXPECT_EQ(25u, res.label_summary[0][kRight].direct_bp);
        EXPECT_EQ(55u, res.label_summary[1][kRight].direct_bp);
        EXPECT_EQ(20u, res.label_summary[0][kLeft].direct_bp);
        EXPECT_EQ(0u, res.label_summary[1][kLeft].direct_bp);

        // brute force: every run end equals the first outward base whose k-mer the
        // label's sequence does not contain
        std::string rflank = f.X + f.Y, lflank = f.W;
        const auto label_seqs = f.seqs();
        for (LabelId l = 0; l < 2; ++l) {
            const std::string &lseq = label_seqs[l];
            uint64_t to = rflank.size();
            for (uint64_t i = 0; i < rflank.size(); ++i) {
                std::string kmer = (f.R + rflank).substr(f.R.size() + i + 1 - kK, kK);
                if (!contains_kmer(lseq, kmer, canonical)) { to = i; break; }
            }
            EXPECT_EQ(to, right.runs[l].to_bp) << "label " << l;
            to = lflank.size();
            for (uint64_t i = 0; i < lflank.size(); ++i) {
                std::string kmer = (lflank + f.R).substr(lflank.size() - 1 - i, kK);
                if (!contains_kmer(lseq, kmer, canonical)) { to = i; break; }
            }
            EXPECT_EQ(to, left.runs[l].to_bp) << "label " << l;
        }
        check_spelled(*anno, right, right.paths[0], f.R);
        check_spelled(*anno, left, left.paths[0], f.R);
    }
}


// T7: a divergence of label sets is a split but not a branch
TYPED_TEST(WalkerTest, Divergence) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t seed = 3; ; ++seed) {
        b = clean_blocks({ 30, 40, 40 }, seed);
        if (b[1][0] != b[2][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);
        Strategy st;
        auto res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right, st);
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
            EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
            EXPECT_EQ(0u, path.end_labels[0].branches);
            leaves[spell_path(right, path)] = path.end_labels[0].label;
            check_spelled(*anno, right, path, X);
        }
        ASSERT_EQ(2u, leaves.size());
        EXPECT_EQ(0u, leaves.at(P));
        EXPECT_EQ(1u, leaves.at(Q));
        EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
        EXPECT_EQ(0u, right.branch_events_total);
        EXPECT_TRUE(right.branch_events.empty());
        size_t divergences = 0, ambiguous = 0;
        for (const auto &g : right.growth) { divergences += g.divergences; ambiguous += g.ambiguous_branches; }
        EXPECT_EQ(1u, divergences);
        EXPECT_EQ(0u, ambiguous);
        EXPECT_EQ(80u, right.steps);
    }
}


// T8: one label on both sides of a structural branch
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
        Strategy st = strategy(0);
        auto res = run(*anno, X, { "A" }, st);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right, st);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(0u, right.paths[0].length_bp);
        EXPECT_EQ(1u, right.paths[0].end_reasons[static_cast<size_t>(EndReason::BRANCH)]);
        ASSERT_EQ(1u, right.runs.size());
        EXPECT_EQ(0u, right.runs[0].to_bp);
        EXPECT_EQ(EndReason::BRANCH, right.runs[0].end_reason);
        EXPECT_EQ(1u, right.branch_events_total);
        ASSERT_EQ(1u, right.branch_events.size());
        const BranchEvent &be = right.branch_events[0];
        EXPECT_EQ(0u, be.at_bp);
        std::vector<char> chars { P[0], Q[0] };
        std::sort(chars.begin(), chars.end());
        EXPECT_EQ(chars, be.chars);
        EXPECT_EQ((std::vector<size_t>{ 1, 1 }), be.labels_per_successor);
        EXPECT_EQ((std::vector<LabelId>{ 0 }), be.ambiguous);
        EXPECT_EQ((std::vector<LabelId>{ 0 }), be.dropped);
        EXPECT_TRUE(right.splits.empty());
        EXPECT_EQ(0u, right.steps);

        // limit 1: both continuations are followed, each charged one branch
        st = strategy(1);
        res = run(*anno, X, { "A" }, st);
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1, st);
        ASSERT_EQ(2u, r1.paths.size());
        std::vector<uint64_t> lengths;
        for (const auto &path : r1.paths) {
            lengths.push_back(path.length_bp);
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(1u, path.end_labels[0].branches);
            EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
            check_spelled(*anno, r1, path, X);
        }
        std::sort(lengths.begin(), lengths.end());
        EXPECT_EQ((std::vector<uint64_t>{ P.size() + Y.size(), Q.size() + Z.size() }), lengths);
        ASSERT_EQ(1u, r1.splits.size());
        EXPECT_TRUE(r1.splits[0].ambiguous);
        size_t ambiguous = 0;
        for (const auto &g : r1.growth) ambiguous += g.ambiguous_branches;
        EXPECT_EQ(1u, ambiguous);
        EXPECT_EQ(0u, count_ends(r1, EndReason::BRANCH));
        EXPECT_EQ(1u, r1.branch_events_total);
        // the growth profile of the first bin: both paths live, one label
        ASSERT_FALSE(r1.growth.empty());
        EXPECT_EQ(2u, r1.growth[0].max_live_paths);
        EXPECT_EQ(1u, r1.growth[0].max_live_labels);
        EXPECT_EQ(2u, r1.growth[0].max_live_pairs);
    }
}


// The fixpoint over excluded sources at an ambiguous node (§6.3) derives every
// successor's state again after an exclusion; the rounds after the first are counted
// (ArmResult::reminimisation_rounds, and the largest at one node).
TYPED_TEST(WalkerTest, ReminimisationRounds) {
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
        // A on both branches, B on the P branch only
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P + Y, X + Q + Z, X + P + Y },
                                                        { "A", "A", "B" }, mode);
        // limit 0, no switching: A is excluded in the first round, B goes on along P
        // alone, so one re-derivation settles the node
        Strategy st = strategy(0);
        auto res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &r0 = res.arms[kRight];
        check_invariants(r0, st);
        EXPECT_EQ(1u, r0.reminimisation_rounds);
        EXPECT_EQ(1u, r0.max_reminimisation_rounds);
        ASSERT_EQ(1u, r0.branch_events.size());
        EXPECT_EQ((std::vector<LabelId>{ 0 }), r0.branch_events[0].ambiguous);
        ASSERT_EQ(1u, r0.paths.size());
        EXPECT_EQ(P.size() + Y.size(), r0.paths[0].length_bp);
        // the seed starts the sequences: nothing to re-minimise on the left
        EXPECT_EQ(0u, res.arms[kLeft].reminimisation_rounds);

        // With a switch within budget the exclusion cascades: once A is excluded, B is
        // the only source left for A's node on Q and switches into it, so B now goes
        // on along both branches and is excluded in the second round. Two
        // re-derivations, and every label ends at the node.
        st.loss_budget = 1;
        res = run(*anno, X, { "A", "B" }, st, LabelChangeCost::constant(1));
        const ArmResult &r2 = res.arms[kRight];
        check_invariants(r2, st);
        EXPECT_EQ(2u, r2.reminimisation_rounds);
        EXPECT_EQ(2u, r2.max_reminimisation_rounds);
        ASSERT_EQ(1u, r2.branch_events.size());
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), r2.branch_events[0].ambiguous);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), r2.branch_events[0].dropped);
        EXPECT_EQ(2u, count_ends(r2, EndReason::BRANCH));
        ASSERT_EQ(1u, r2.paths.size());
        EXPECT_EQ(0u, r2.paths[0].length_bp);

        // an ambiguity within the allowance excludes nothing: no re-derivation
        st = strategy(1);
        res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1, st);
        EXPECT_EQ(0u, r1.reminimisation_rounds);
        EXPECT_EQ(0u, r1.max_reminimisation_rounds);
        EXPECT_EQ(2u, r1.paths.size());
    }
}


// T10: diamond. A on X·P·Y·Q, B on X·P'·Y·Q'
TYPED_TEST(WalkerTest, Diamond) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::vector<std::string> b;
    for (uint32_t seed = 5; ; ++seed) {
        b = clean_blocks({ 30, 25, 25, 30, 25, 25 }, seed);
        // P and P' must differ at both ends: the paths diverge at the first base and
        // reconverge exactly at the first k-mer inside Y
        if (b[1][0] != b[2][0] && b[1].back() != b[2].back() && b[4][0] != b[5][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &P2 = b[2], &Y = b[3], &Q = b[4], &Q2 = b[5];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P + Y + Q, X + P2 + Y + Q2 },
                                                        { "A", "B" }, mode);
        // no merging: two independent leaves, each label on its own sequence
        Strategy st = strategy(0, false);
        auto res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right, st);
        ASSERT_EQ(2u, right.paths.size());
        std::map<std::string, std::vector<LabelId>> leaves;
        for (const auto &path : right.paths) {
            std::vector<LabelId> labels;
            for (const auto &e : path.end_labels) labels.push_back(e.label);
            leaves[spell_path(right, path)] = labels;
            check_spelled(*anno, right, path, X);
        }
        ASSERT_EQ(2u, leaves.size());
        EXPECT_EQ((std::vector<LabelId>{ 0 }), leaves.at(P + Y + Q));
        EXPECT_EQ((std::vector<LabelId>{ 1 }), leaves.at(P2 + Y + Q2));
        EXPECT_EQ(1u, right.splits.size());
        EXPECT_EQ(0u, count_events(right, EventType::RECONVERGE));
        // the second path shares all |Y|-k+1 nodes of Y with the first one at the
        // same distance: reported once per stretch, not once per node
        ASSERT_EQ(1u, count_events(right, EventType::REVISIT));
        EXPECT_EQ(P.size() + kK, events_of(right, EventType::REVISIT)[0]->at_bp);
        EXPECT_EQ("same_distance", events_of(right, EventType::REVISIT)[0]->text);
        size_t reconv = 0;
        for (const auto &g : right.growth) reconv += g.reconvergences;
        EXPECT_EQ(0u, reconv);

        // merging: the Y segment is shared with two parents and splits again
        st = strategy(0, true);
        res = run(*anno, X, { "A", "B" }, st);
        const ArmResult &rm = res.arms[kRight];
        check_invariants(rm, st);
        const Segment *y = nullptr;
        for (const auto &seg : rm.segments) {
            if (seg.parents.size() == 2) {
                ASSERT_EQ(nullptr, y);
                y = &seg;
            }
        }
        ASSERT_NE(nullptr, y);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), y->labels_start);
        EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), y->labels_end);
        EXPECT_EQ(P.size() + kK, y->from_bp);
        EXPECT_EQ(Y.size() - kK, y->length_bp);
        EXPECT_EQ(Y.substr(kK), y->sequence);
        ASSERT_EQ(2u, y->children.size());
        std::map<std::string, std::vector<LabelId>> children;
        for (size_t c : y->children) {
            children[rm.segments[c].sequence] = rm.segments[c].labels_start;
            EXPECT_EQ((std::vector<size_t>{ y->id }), rm.segments[c].parents);
        }
        EXPECT_EQ((std::vector<LabelId>{ 0 }), children.at(Q));
        EXPECT_EQ((std::vector<LabelId>{ 1 }), children.at(Q2));
        ASSERT_EQ(1u, count_events(rm, EventType::RECONVERGE));
        const Event &ev = y->events.front();
        EXPECT_EQ(EventType::RECONVERGE, ev.type);
        EXPECT_EQ(y->from_bp, ev.at_bp);
        EXPECT_EQ(2u, ev.labels.size());
        reconv = 0;
        for (const auto &g : rm.growth) reconv += g.reconvergences;
        EXPECT_EQ(1u, reconv);
        EXPECT_EQ(2u, rm.splits.size());
        EXPECT_EQ(2u, rm.paths.size());
        EXPECT_EQ(0u, count_ends(rm, EndReason::BRANCH));
        // the parents' runs closed by the merge are not label ends
        EXPECT_EQ(2u, count_ends(rm, EndReason::DEAD_END));
        EXPECT_EQ(res.label_summary[0][kRight].direct_bp, P.size() + Y.size() + Q.size());
        EXPECT_EQ(res.label_summary[1][kRight].direct_bp, P2.size() + Y.size() + Q2.size());

        // per-label routes stay reconstructible (§6.5): A entered Y through P, B
        // through P', and each end label supports the seed plus its own route
        const bool canonical = mode != DeBruijnGraph::BASIC;
        const std::string seqA = X + P + Y + Q, seqB = X + P2 + Y + Q2;
        ASSERT_EQ(2u, y->labels_via_parent.size());
        // the parents spell P / P' plus the first k bases of Y (the merge node)
        std::map<std::string, std::vector<LabelId>> via;
        for (size_t p = 0; p < 2; ++p) via[rm.segments[y->parents[p]].sequence] = y->labels_via_parent[p];
        ASSERT_EQ(1u, via.count(P + Y.substr(0, kK)));
        ASSERT_EQ(1u, via.count(P2 + Y.substr(0, kK)));
        EXPECT_EQ((std::vector<LabelId>{ 0 }), via.at(P + Y.substr(0, kK)));
        EXPECT_EQ((std::vector<LabelId>{ 1 }), via.at(P2 + Y.substr(0, kK)));
        for (const auto &path : rm.paths) {
            ASSERT_EQ(1u, path.end_labels.size());
            LabelId l = path.end_labels[0].label;
            std::string route = route_of(rm, path.leaf, l);
            EXPECT_EQ(l == 0 ? P + Y + Q : P2 + Y + Q2, route);
            EXPECT_TRUE(supports_sequence(l == 0 ? seqA : seqB, X + route, canonical));
        }

        // continuations after a merge list only labels covering the spelled tail
        // (§7.1: valid /traverse input): stopping inside Y (merged leaf) and right
        // after the second split (one leaf per label, spelled through the first parent)
        for (uint64_t extra : { 10u, 20u }) {
            Strategy sc = strategy(0, true);
            sc.max_extension_bp = P.size() + kK + extra;
            auto rc_ = run(*anno, X, { "A", "B" }, sc);
            const ArmResult &arm = rc_.arms[kRight];
            check_invariants(arm, sc);
            ASSERT_EQ(extra == 10u ? 1u : 2u, arm.paths.size()) << mode << " extra " << extra;
            for (const auto &path : arm.paths) {
                ASSERT_TRUE(path.continuation.has_value());
                const Continuation &c = *path.continuation;
                ASSERT_FALSE(c.labels.empty());
                EXPECT_GE(c.sequence.size(), kK);
                std::vector<std::string> names;
                for (LabelId l : c.labels) {
                    names.push_back(rc_.label_dict[l].name);
                    EXPECT_TRUE(supports_sequence(l == 0 ? seqA : seqB, c.sequence, canonical))
                        << mode << " extra " << extra << " label " << l << " " << c.sequence;
                }
                // re-submitting the continuation as a seed keeps every listed label
                auto again = run(*anno, c.sequence, names, Strategy());
                EXPECT_TRUE(again.dropped_labels.empty()) << mode << " extra " << extra;
                EXPECT_EQ(names.size(), again.num_seed_labels);
                // every end label supports its own route
                for (const auto &e : path.end_labels) {
                    std::string route = route_of(arm, path.leaf, e.label);
                    EXPECT_TRUE(supports_sequence(e.label == 0 ? seqA : seqB, X + route, canonical))
                        << mode << " extra " << extra << " label " << e.label;
                }
            }
            if (extra == 20u) {
                // both labels are live inside Y but the Q' leaf is spelled through P
                std::set<std::vector<LabelId>> label_sets;
                for (const auto &path : arm.paths) label_sets.insert(path.continuation->labels);
                EXPECT_EQ((std::set<std::vector<LabelId>>{ { 0 }, { 1 } }), label_sets);
            }
        }
    }
}


// T11a: simple cycle L·J·C·J·E with J a single k-mer
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
        Strategy st = strategy(0);
        auto res = run(*anno, L, { "A" }, st);
        const ArmResult &r0 = res.arms[kRight];
        check_invariants(r0, st);
        ASSERT_EQ(1u, r0.paths.size());
        EXPECT_EQ(kK, r0.paths[0].length_bp);
        EXPECT_EQ(J, spell_path(r0, r0.paths[0]));
        EXPECT_EQ(1u, count_ends(r0, EndReason::BRANCH));
        EXPECT_EQ(kK, r0.runs[0].to_bp);

        st = strategy(1);
        res = run(*anno, L, { "A" }, st);
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1, st);
        ASSERT_EQ(2u, r1.paths.size());
        std::set<std::string> flanks;
        for (const auto &path : r1.paths) {
            flanks.insert(spell_path(r1, path));
            ASSERT_EQ(1u, path.end_labels.size());
            EXPECT_EQ(1u, path.end_labels[0].branches);
            EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
            check_spelled(*anno, r1, path, L);
        }
        EXPECT_EQ((std::set<std::string>{ J + E, J + C + J + E }), flanks);
        EXPECT_EQ(1u, count_events(r1, EventType::BLOCKED));
        EXPECT_EQ(1u, count_events(r1, EventType::REVISIT));
        for (const auto &seg : r1.segments) {
            for (const auto &ev : seg.events) {
                if (ev.type == EventType::BLOCKED) {
                    EXPECT_EQ(EndReason::EDGE_REUSE, ev.reason);
                    EXPECT_EQ(C[0], ev.ch);
                    EXPECT_EQ(J.size() + C.size() + J.size(), ev.at_bp);
                    EXPECT_EQ((std::vector<LabelId>{ 0 }), ev.labels);
                }
            }
        }
        size_t blocked = 0;
        for (const auto &g : r1.growth) blocked += g.blocked_repeat;
        EXPECT_EQ(1u, blocked);
        EXPECT_EQ(0u, count_ends(r1, EndReason::BRANCH));
        EXPECT_EQ(0u, count_ends(r1, EndReason::EDGE_REUSE));
    }
}

// 2^n labels through a cascade of n bubbles (A_i | B_i of equal length, then a shared
// block S_i) and on through the |tail| blocks: label c takes B_i where bit i of c is
// set, so every label has its own path and, with merging off, every one of the 2^n
// paths walks the shared blocks as a head of its own.
struct CascadeFixture {
    size_t n;
    std::vector<std::string> blocks;     // X, A_1, B_1, S_1, ..., A_n, B_n, S_n, tail...

    CascadeFixture(size_t n, const std::vector<size_t> &tail_lengths, uint32_t seed,
                   std::function<bool(const CascadeFixture&)> accept = nullptr) : n(n) {
        std::vector<size_t> lengths { 30 };
        for (size_t i = 0; i < n; ++i) {
            lengths.insert(lengths.end(), { 20, 20, 15 });
        }
        lengths.insert(lengths.end(), tail_lengths.begin(), tail_lengths.end());
        for (uint32_t s = seed; ; ++s) {
            blocks = clean_blocks(lengths, s);
            bool ok = true;
            for (size_t i = 0; i < n; ++i) {
                ok &= blocks[1 + 3 * i][0] != blocks[2 + 3 * i][0];
            }
            if (ok && (!accept || accept(*this)))
                break;
        }
    }
    const std::string& X() const { return blocks[0]; }
    const std::string& tail(size_t j) const { return blocks[1 + 3 * n + j]; }
    std::string bubbles(size_t c) const {
        std::string s;
        for (size_t i = 0; i < n; ++i) {
            s += blocks[1 + 3 * i + ((c >> i) & 1)] + blocks[3 + 3 * i];
        }
        return s;
    }
    std::vector<std::string> labels() const {
        std::vector<std::string> out;
        for (size_t c = 0; c < (size_t(1) << n); ++c) out.push_back("L" + std::to_string(c));
        return out;
    }
};

// The per-path edge-reuse check costs a step at most the path's depth in segments,
// not the number of live paths that took the same edge (ArmResult::edge_reuse_probes).
// With merging off (keep mode, the exhaustive preset) the 2^5 paths of a cascade all
// walk the shared tail and all record its edges: scanning those uses costs ~P/2 per
// step there, probing the path's own segments at most n + 1.
TYPED_TEST(WalkerTest, EdgeReuseProbes) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    {
        const size_t n = 5;
        CascadeFixture f(n, { 100 }, 11);
        std::vector<std::string> seqs;
        for (size_t c = 0; c < (size_t(1) << n); ++c) seqs.push_back(f.X() + f.bubbles(c) + f.tail(0));
        for (auto mode : all_modes()) {
            auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, f.labels(), mode);
            Strategy st = strategy(0, false);
            st.direction = Strategy::RIGHT;
            auto res = run(*anno, f.X(), f.labels(), st);
            const ArmResult &right = res.arms[kRight];
            check_invariants(right, st);
            ASSERT_EQ(seqs.size(), right.paths.size());
            for (const auto &path : right.paths) {
                ASSERT_EQ(1u, path.end_labels.size());
                EXPECT_EQ(seqs[path.end_labels[0].label].substr(f.X().size()),
                          spell_path(right, path));
            }
            EXPECT_EQ(0u, count_events(right, EventType::BLOCKED));
            // every path crosses the n splits, so a head has at most n + 1 segments
            EXPECT_GT(right.edge_reuse_probes, 0u);
            EXPECT_LE(right.edge_reuse_probes, (n + 1) * right.steps);
        }
    }
    // The probe has to find a use by the path's own segments and ignore the sibling
    // paths' uses of the same edge: 2^3 paths into the cycle L·J·C·J·E of T11a. Every
    // path takes J→C once (eight uses, more than its segments) and is blocked at it
    // the second time; J→E, taken by the siblings only, stays open.
    {
        const size_t n = 3;
        CascadeFixture f(n, { 15, kK, 40, 30 }, 12,
                         [](const CascadeFixture &x) { return x.tail(2)[0] != x.tail(3)[0]; });
        const std::string &L = f.tail(0), &J = f.tail(1), &C = f.tail(2), &E = f.tail(3);
        std::vector<std::string> seqs;
        std::set<std::string> expected;
        for (size_t c = 0; c < (size_t(1) << n); ++c) {
            seqs.push_back(f.X() + f.bubbles(c) + L + J + C + J + E);
            expected.insert(f.bubbles(c) + L + J + E);
            expected.insert(f.bubbles(c) + L + J + C + J + E);
        }
        for (auto mode : all_modes()) {
            auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, f.labels(), mode);
            Strategy st = strategy(1, false);
            st.direction = Strategy::RIGHT;
            auto res = run(*anno, f.X(), f.labels(), st);
            const ArmResult &right = res.arms[kRight];
            check_invariants(right, st);
            std::set<std::string> flanks;
            for (const auto &path : right.paths) {
                flanks.insert(spell_path(right, path));
                ASSERT_EQ(1u, path.end_labels.size());
                EXPECT_EQ(1u, path.end_labels[0].branches);
                EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
                check_spelled(*anno, right, path, f.X());
            }
            EXPECT_EQ(expected, flanks);
            EXPECT_EQ(2 * seqs.size(), right.paths.size());
            auto blocked = events_of(right, EventType::BLOCKED);
            EXPECT_EQ(seqs.size(), blocked.size());
            for (const Event *ev : blocked) {
                EXPECT_EQ(EndReason::EDGE_REUSE, ev->reason);
                EXPECT_EQ(C[0], ev->ch);
                EXPECT_EQ(1u, ev->labels.size());
            }
            EXPECT_EQ(0u, count_ends(right, EndReason::EDGE_REUSE));
            EXPECT_EQ(0u, count_ends(right, EndReason::EDGE_REUSE_RC));
        }
    }
}

// T11b: a homopolymer self-loop is taken once, then blocked
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
    // The graph is node-centric: the node L[-1]·A^(k-1) already has the two successors
    // A^k and A^(k-1)·E[0], and A^k has the self-loop plus the same exit. Unrolling
    // the loop once therefore needs two branches on the lineage.
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        Strategy st = strategy(0);
        auto res = run(*anno, L, { "A" }, st);
        const ArmResult &r0 = res.arms[kRight];
        check_invariants(r0, st);
        ASSERT_EQ(1u, r0.paths.size());
        EXPECT_EQ(kK - 1, r0.paths[0].length_bp);
        EXPECT_EQ(1u, count_ends(r0, EndReason::BRANCH));
        check_spelled(*anno, r0, r0.paths[0], L);

        st = strategy(1);
        res = run(*anno, L, { "A" }, st);
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1, st);
        ASSERT_EQ(2u, r1.paths.size());
        std::set<std::string> flanks;
        for (const auto &path : r1.paths) {
            flanks.insert(spell_path(r1, path));
            check_spelled(*anno, r1, path, L);
        }
        EXPECT_EQ((std::set<std::string>{ H.substr(1), H.substr(2) + E }), flanks);
        EXPECT_EQ(1u, count_ends(r1, EndReason::BRANCH));
        EXPECT_EQ(1u, count_ends(r1, EndReason::DEAD_END));
        EXPECT_EQ(0u, count_events(r1, EventType::BLOCKED));

        // limit 2: the self-loop is taken once and then blocked as an edge reuse
        st = strategy(2);
        res = run(*anno, L, { "A" }, st);
        const ArmResult &r2 = res.arms[kRight];
        check_invariants(r2, st);
        ASSERT_EQ(3u, r2.paths.size());
        flanks.clear();
        for (const auto &path : r2.paths) {
            flanks.insert(spell_path(r2, path));
            check_spelled(*anno, r2, path, L);
            EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
        }
        EXPECT_EQ((std::set<std::string>{ H + E, H.substr(1) + E, H.substr(2) + E }), flanks);
        EXPECT_EQ(1u, count_events(r2, EventType::BLOCKED));
        EXPECT_GE(count_events(r2, EventType::REVISIT), 1u);
        for (const auto &seg : r2.segments) {
            for (const auto &ev : seg.events) {
                if (ev.type == EventType::BLOCKED) {
                    EXPECT_EQ(EndReason::EDGE_REUSE, ev.reason);
                    EXPECT_EQ('A', ev.ch);
                    EXPECT_EQ(kK + 1, ev.at_bp);
                }
            }
        }
        EXPECT_EQ(0u, count_ends(r2, EndReason::BRANCH));
    }
}

// T11c: a circular molecule: each arm walks around the circle back into the seed
TYPED_TEST(WalkerTest, Circle) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(150, 8);
    std::string record = S + S.substr(0, kK - 1 + 10);
    std::string seed = S.substr(50, 40);
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { record }, { "A" }, mode);
        Strategy st;
        auto res = run(*anno, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm, st);
            ASSERT_EQ(1u, arm.paths.size()) << mode << " arm " << a;
            EXPECT_EQ(ArmResult::COMPLETE, arm.status);
            EXPECT_EQ(1u, arm.paths[0].end_reasons[static_cast<size_t>(EndReason::REACHED_SEED)]);
            EXPECT_EQ(S.size() - seed.size() + kK - 1, arm.paths[0].length_bp);
            EXPECT_EQ(1u, count_events(arm, EventType::BLOCKED));
            check_spelled(*anno, arm, arm.paths[0], seed);
        }
        EXPECT_EQ(S.substr(90) + S.substr(0, 60), spell_path(res.arms[kRight], res.arms[kRight].paths[0]));
        EXPECT_EQ(S.substr(80) + S.substr(0, 50), spell_path(res.arms[kLeft], res.arms[kLeft].paths[0]));
    }
}


// T12: an RC-palindromic (k-1)-mer creates a hairpin successor in canonical graphs
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
        // the real continuation must not itself be the hairpin step
        if (R[0] != rc(std::string(1, L.back()))[0] && is_clean(L + X.substr(0, 3))
                && is_clean(X.substr(7) + R))
            break;
    }
    std::string S = L + X + R;
    std::string seed = L.substr(0, 30);
    for (auto mode : canonical_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        Strategy st;
        auto res = run(*anno, seed, { "A" }, st);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right, st);
        ASSERT_EQ(1u, right.paths.size()) << mode;
        EXPECT_EQ(S.substr(30), spell_path(right, right.paths[0])) << mode;
        EXPECT_EQ(0u, count_ends(right, EndReason::BRANCH));
        EXPECT_EQ(0u, right.branch_events_total);
        EXPECT_EQ(1u, count_events(right, EventType::HAIRPIN));
        EXPECT_EQ(0u, count_events(right, EventType::BLOCKED));
        for (const auto &ev : right.segments[0].events) {
            if (ev.type == EventType::HAIRPIN) {
                // the hairpin successor hangs off the node a·X, reached after the
                // remaining 10 bases of L plus the k-1 bases of X
                EXPECT_EQ(L.size() - 30 + kK - 1, ev.at_bp);
                EXPECT_EQ(rc(std::string(1, L.back()))[0], ev.ch);
                EXPECT_EQ((std::vector<LabelId>{ 0 }), ev.labels);
            }
        }
        EXPECT_EQ(0u, count_events(res.arms[kLeft], EventType::HAIRPIN));
        EXPECT_EQ(L.substr(0, 0), spell_path(res.arms[kLeft], res.arms[kLeft].paths[0]));

        // following hairpins instead of skipping them (§6.5 "follow"): the hairpin
        // step is followed and flagged, never counts toward ambiguity, and the RC
        // retrace that follows it is blocked as edge_reuse_rc
        st.skip_hairpins = false;
        st.max_label_branches = 0;
        res = run(*anno, seed, { "A" }, st);
        const ArmResult &rf = res.arms[kRight];
        check_invariants(rf, st);
        check_invariants(res.arms[kLeft], st);
        EXPECT_EQ(0u, count_ends(rf, EndReason::BRANCH));
        EXPECT_EQ(0u, rf.branch_events_total);
        ASSERT_EQ(1u, count_events(rf, EventType::HAIRPIN)) << mode;
        const Event *hp = events_of(rf, EventType::HAIRPIN)[0];
        EXPECT_EQ(L.size() - 30 + kK - 1, hp->at_bp);
        EXPECT_EQ(rc(std::string(1, L.back()))[0], hp->ch);
        EXPECT_EQ((std::vector<LabelId>{ 0 }), hp->labels);
        EXPECT_EQ("followed", hp->text);
        ASSERT_EQ(1u, rf.splits.size());
        EXPECT_FALSE(rf.splits[0].ambiguous);
        EXPECT_EQ(hp->at_bp, rf.splits[0].at_bp);
        ASSERT_EQ(2u, rf.paths.size());
        bool retrace = false, real = false;
        for (const auto &path : rf.paths) {
            std::string flank = spell_path(rf, path);
            check_spelled(*anno, rf, path, seed);
            if (flank == S.substr(30)) {
                real = true;
                EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
            } else {
                retrace = true;
                EXPECT_EQ(hp->at_bp + 1, path.length_bp);
                EXPECT_EQ(rc(std::string(1, L.back()))[0], flank[hp->at_bp]);
                EXPECT_EQ(1u, path.end_reasons[static_cast<size_t>(EndReason::EDGE_REUSE_RC)]
                              + path.end_reasons[static_cast<size_t>(EndReason::REACHED_SEED)]);
            }
        }
        EXPECT_TRUE(real && retrace) << mode;
        EXPECT_EQ(1u, count_events(rf, EventType::BLOCKED));
    }
}


// T14: unmasked DBGSuccinct (dummy k-mers present): identical results, '$' is never a branch
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
            for (const auto &c : contigs) base->add_sequence(c);
            graph = std::make_shared<CanonicalDBG>(base);
        } else {
            base = std::make_shared<DBGSuccinct>(kK, mode);
            for (const auto &s : f.seqs()) base->add_sequence(s);
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
            Strategy st;
            auto a = run(*masked, seed, labels, st);
            auto b = run(*unmasked, seed, labels, st);
            EXPECT_EQ(serialize(a), serialize(b)) << mode << " seed " << seed;
            // tips end with dead_end ('$' sinks and sources are never successors),
            // and the dummy k-mers never count as a branch
            size_t dead_ends = 0;
            for (size_t arm : { kLeft, kRight }) {
                EXPECT_EQ(0u, count_ends(b.arms[arm], EndReason::BRANCH));
                EXPECT_EQ(0u, b.arms[arm].branch_events_total);
                EXPECT_TRUE(b.arms[arm].splits.empty());
                dead_ends += count_ends(b.arms[arm], EndReason::DEAD_END);
                for (const auto &run : b.arms[arm].runs) {
                    EXPECT_TRUE(run.end_reason == EndReason::DEAD_END
                                || run.end_reason == EndReason::LABEL_LOST);
                }
            }
            EXPECT_GE(dead_ends, 1u);
        }
    }
}


// §8.1: CanonicalDBG's empty runtime_error on an inconsistent primary graph (both
// strands stored) is rethrown naming the node, the k-mer and the arm
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
        traverse_seed(oracle, seed, Strategy(), LabelChangeCost::forbid());
        FAIL() << "no exception";
    } catch (const std::runtime_error &e) {
        std::string what = e.what();
        EXPECT_NE(std::string::npos, what.find("Inconsistent primary graph at node")) << what;
        EXPECT_NE(std::string::npos, what.find("on arm")) << what;
    }
}


// T15: caps
TYPED_TEST(WalkerTest, Caps) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::string S = clean_block(200, 10);
    std::string seed = S.substr(60, 40);
    std::vector<std::string> b;
    for (uint32_t s = 11; ; ++s) {
        b = clean_blocks({ 30, 40, 40 }, s);
        if (b[1][0] != b[2][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto linear = build_anno_graph<Graph, Annotation>(kK, { S }, { "A" }, mode);
        auto fork = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);

        // max_extension_bp: the requested domain is complete
        Strategy st;
        st.max_extension_bp = 10;
        auto res = run(*linear, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm, st);
            EXPECT_EQ(ArmResult::COMPLETE, arm.status);
            EXPECT_FALSE(arm.cap_trigger.has_value());
            ASSERT_EQ(1u, arm.paths.size());
            EXPECT_EQ(10u, arm.paths[0].length_bp);
            ASSERT_TRUE(arm.paths[0].path_reason.has_value());
            EXPECT_EQ(EndReason::MAX_EXTENSION, *arm.paths[0].path_reason);
            EXPECT_EQ(1u, arm.paths[0].end_reasons[static_cast<size_t>(EndReason::MAX_EXTENSION)]);
            ASSERT_TRUE(arm.paths[0].continuation.has_value());
            const Continuation &c = *arm.paths[0].continuation;
            EXPECT_EQ(kK, c.sequence.size());
            EXPECT_EQ((std::vector<LabelId>{ 0 }), c.labels);
            EXPECT_EQ(0.0, c.loss_used);
            EXPECT_EQ(0u, c.branches_used);
            EXPECT_EQ(10u, arm.steps);
        }
        EXPECT_EQ(S.substr(99, kK), res.arms[kRight].paths[0].continuation->sequence);
        EXPECT_EQ(S.substr(50, kK), res.arms[kLeft].paths[0].continuation->sequence);
        EXPECT_EQ(10u, res.label_summary[0][kRight].reach_bp);

        // max_steps: per seed, both arms end
        st = Strategy();
        st.max_steps = 5;
        res = run(*linear, seed, { "A" }, st);
        EXPECT_EQ(5u, res.arms[kLeft].steps + res.arms[kRight].steps);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm, st);
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
            ASSERT_TRUE(arm.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_STEPS, arm.cap_trigger->reason);
            EXPECT_EQ(1u, arm.cap_trigger->live_paths);
            EXPECT_EQ(1u, arm.cap_trigger->live_labels);
            ASSERT_EQ(1u, arm.paths.size());
            EXPECT_EQ(1u, arm.paths[0].end_reasons[static_cast<size_t>(EndReason::MAX_STEPS)]);
            EXPECT_EQ(EndReason::MAX_STEPS, *arm.paths[0].path_reason);
            EXPECT_TRUE(arm.paths[0].continuation.has_value());
            EXPECT_EQ(1u, count_ends(arm, EndReason::MAX_STEPS));
        }

        // max_output_bp: per arm
        st = Strategy();
        st.max_output_bp = 10;
        res = run(*linear, seed, { "A" }, st);
        for (size_t a : { kLeft, kRight }) {
            const ArmResult &arm = res.arms[a];
            check_invariants(arm, st);
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
            EXPECT_EQ(10u, arm.output_bp);
            EXPECT_EQ(10u, arm.steps);
            ASSERT_TRUE(arm.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_OUTPUT, arm.cap_trigger->reason);
            EXPECT_EQ(10u, arm.cap_trigger->at_bp);
            EXPECT_EQ(1u, count_ends(arm, EndReason::MAX_OUTPUT));
        }

        // max_live_paths (stop): the split trips it on the right arm only
        st = Strategy();
        st.max_live_paths = 1;
        res = run(*fork, X, { "A", "B" }, st);
        {
            const ArmResult &right = res.arms[kRight], &left = res.arms[kLeft];
            check_invariants(right, st);
            EXPECT_EQ(ArmResult::TRUNCATED, right.status);
            ASSERT_TRUE(right.cap_trigger.has_value());
            EXPECT_EQ(EndReason::MAX_LIVE_PATHS, right.cap_trigger->reason);
            EXPECT_EQ(0u, right.cap_trigger->at_bp);
            EXPECT_EQ(2u, right.cap_trigger->live_labels);
            EXPECT_EQ(2u, count_ends(right, EndReason::MAX_LIVE_PATHS));
            EXPECT_TRUE(right.splits.empty());
            EXPECT_EQ(ArmResult::COMPLETE, left.status);
            EXPECT_EQ(2u, count_ends(left, EndReason::DEAD_END));
        }
        // beam: one path survives, the other is pruned
        st.on_overflow = Strategy::BEAM;
        res = run(*fork, X, { "A", "B" }, st);
        {
            const ArmResult &right = res.arms[kRight];
            check_invariants(right, st);
            EXPECT_EQ(ArmResult::PRUNED, right.status);
            ASSERT_TRUE(right.cap_trigger.has_value());
            EXPECT_EQ(EndReason::BEAM_PRUNED, right.cap_trigger->reason);
            ASSERT_EQ(2u, right.paths.size());
            EXPECT_EQ(1u, count_ends(right, EndReason::BEAM_PRUNED));
            EXPECT_EQ(1u, count_ends(right, EndReason::DEAD_END));
            size_t pruned = 0;
            for (const auto &path : right.paths) {
                if (path.path_reason && *path.path_reason == EndReason::BEAM_PRUNED) {
                    ++pruned;
                    EXPECT_EQ(1u, path.length_bp);
                    EXPECT_TRUE(path.continuation.has_value());
                } else {
                    EXPECT_EQ(40u, path.length_bp);
                }
            }
            EXPECT_EQ(1u, pruned);
        }

        // max_paths
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

        // time budget 0: stops after the first level of each arm
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
            EXPECT_EQ(1u, count_ends(arm, EndReason::TIME_BUDGET));
        }
    }
}


// T9: switching chain {A,B} -> {A} -> {A,B} -> {B}
struct ChainFixture {
    std::string b1, b2, b3, b4;
    ChainFixture() {
        auto b = clean_blocks({ 40, 30, 30, 30 }, 12);
        b1 = b[0]; b2 = b[1]; b3 = b[2]; b4 = b[3];
    }
    std::vector<std::string> seqs() const { return { b1 + b2 + b3, b1, b3 + b4 }; }
    std::vector<std::string> labels() const { return { "A", "B", "B" }; }
};

TYPED_TEST(WalkerTest, SwitchingChain) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    ChainFixture f;
    const uint64_t l2 = f.b2.size(), l3 = f.b3.size(), l4 = f.b4.size();
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);

        // forbid
        Strategy st;
        auto res = run(*anno, f.b1, { "A", "B" }, st);
        const ArmResult &rf = res.arms[kRight];
        check_invariants(rf, st);
        ASSERT_EQ(1u, rf.paths.size());
        EXPECT_EQ(l2 + l3, rf.paths[0].length_bp);
        ASSERT_EQ(2u, rf.runs.size());
        EXPECT_EQ(0u, rf.runs[0].label);
        EXPECT_EQ(l2 + l3, rf.runs[0].to_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, rf.runs[0].end_reason);
        EXPECT_EQ(1u, rf.runs[1].label);
        EXPECT_EQ(0u, rf.runs[1].to_bp);
        EXPECT_EQ(EndReason::LABEL_LOST, rf.runs[1].end_reason);
        EXPECT_EQ(0u, count_events(rf, EventType::SWITCH));
        EXPECT_EQ(2u, count_ends(res.arms[kLeft], EndReason::DEAD_END));
        EXPECT_EQ(0u, res.arms[kLeft].paths[0].length_bp);

        // constant cost 1 within budget 1: A -> B at the b3/b4 boundary
        st.loss_budget = 1;
        res = run(*anno, f.b1, { "A", "B" }, st, LabelChangeCost::constant(1));
        const ArmResult &rc1 = res.arms[kRight];
        check_invariants(rc1, st);
        ASSERT_EQ(1u, rc1.paths.size());
        EXPECT_EQ(l2 + l3 + l4, rc1.paths[0].length_bp);
        EXPECT_EQ(f.b2 + f.b3 + f.b4, spell_path(rc1, rc1.paths[0]));
        ASSERT_EQ(1u, rc1.paths[0].end_labels.size());
        EXPECT_EQ(1u, rc1.paths[0].end_labels[0].label);
        EXPECT_EQ(1.0, rc1.paths[0].end_labels[0].loss);
        EXPECT_EQ(1u, rc1.paths[0].end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
        ASSERT_EQ(1u, count_events(rc1, EventType::SWITCH));
        const Event *sw = nullptr;
        for (const auto &ev : rc1.segments[0].events) {
            if (ev.type == EventType::SWITCH) sw = &ev;
        }
        ASSERT_NE(nullptr, sw);
        EXPECT_EQ(l2 + l3, sw->at_bp);
        EXPECT_EQ(0u, sw->label);
        EXPECT_EQ(1u, sw->to);
        EXPECT_EQ(1.0, sw->cost);
        ASSERT_EQ(3u, rc1.runs.size());
        EXPECT_EQ(l2 + l3, rc1.runs[0].to_bp);
        EXPECT_TRUE(rc1.runs[0].ended);
        const LabelRun &switched = rc1.runs[2];
        EXPECT_EQ(1u, switched.label);
        EXPECT_TRUE(switched.entered_by_switch);
        EXPECT_EQ(0u, switched.from_label);
        EXPECT_EQ(1.0, switched.switch_cost);
        EXPECT_EQ(l2 + l3, switched.from_bp);
        EXPECT_EQ(l2 + l3 + l4, switched.to_bp);
        EXPECT_EQ(0u, switched.prev_run);
        EXPECT_EQ(EndReason::DEAD_END, switched.end_reason);
        EXPECT_EQ(l2 + l3, res.label_summary[0][kRight].direct_bp);
        EXPECT_EQ(l2 + l3 + l4, res.label_summary[0][kRight].reach_bp);
        EXPECT_EQ(0u, res.label_summary[1][kRight].direct_bp);
        EXPECT_EQ(1u, res.label_summary[1][kRight].reentries);
        EXPECT_EQ((std::vector<uint32_t>{ 1, 2 }), res.label_summary[1][kRight].runs);

        // cost 1 above budget 0: the switch is reported as needed budget
        st.loss_budget = 0;
        res = run(*anno, f.b1, { "A", "B" }, st, LabelChangeCost::constant(1));
        EXPECT_EQ(l2 + l3, res.arms[kRight].paths[0].length_bp);
        EXPECT_EQ(1u, count_ends(res.arms[kRight], EndReason::LOSS_BUDGET));
        EXPECT_EQ((std::vector<double>{ 1.0 }), res.arms[kRight].needed_budgets);

        // cost 0: same path, loss 0
        res = run(*anno, f.b1, { "A", "B" }, st, LabelChangeCost::constant(0));
        const ArmResult &rc0 = res.arms[kRight];
        check_invariants(rc0, st);
        ASSERT_EQ(1u, rc0.paths.size());
        EXPECT_EQ(f.b2 + f.b3 + f.b4, spell_path(rc0, rc0.paths[0]));
        ASSERT_EQ(1u, rc0.paths[0].end_labels.size());
        EXPECT_EQ(0.0, rc0.paths[0].end_labels[0].loss);
        EXPECT_EQ(1u, count_events(rc0, EventType::SWITCH));
        EXPECT_EQ(l2 + l3, rc0.runs[0].to_bp);

        // switch_on any: B enters as soon as it is present (the first k-mer inside
        // b3, decided at the node before it), A continues
        st.switch_on_loss_only = false;
        res = run(*anno, f.b1, { "A", "B" }, st, LabelChangeCost::constant(0));
        const ArmResult &ra = res.arms[kRight];
        check_invariants(ra, st);
        ASSERT_EQ(1u, count_events(ra, EventType::SWITCH));
        for (const auto &ev : ra.segments[0].events) {
            if (ev.type == EventType::SWITCH) EXPECT_EQ(l2 + kK - 1, ev.at_bp);
        }
        ASSERT_EQ(3u, ra.runs.size());
        EXPECT_EQ(l2 + kK - 1, ra.runs[2].from_bp);
        EXPECT_FALSE(ra.runs[0].entered_by_switch);
        EXPECT_EQ(EndReason::LABEL_LOST, ra.runs[0].end_reason);
        EXPECT_EQ(f.b2 + f.b3 + f.b4, spell_path(ra, ra.paths[0]));
        EXPECT_EQ(l2 + l3, ra.runs[0].to_bp);

        // a table cost: the pair (A -> B) priced below the default
        st.switch_on_loss_only = true;
        st.loss_budget = 0.5;
        auto table = LabelChangeCost::table({ { { 0, 1 }, 0.5 } }, kInfiniteLoss);
        res = run(*anno, f.b1, { "A", "B" }, st, table);
        EXPECT_EQ(l2 + l3 + l4, res.arms[kRight].paths[0].length_bp);
        EXPECT_EQ(0.5, res.arms[kRight].paths[0].end_labels[0].loss);
    }
}


// T17b: seed validation
TYPED_TEST(WalkerTest, SeedValidation) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    LabelEndFixture f;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        // A does not support the last k-mers of R·X·Y[0,5): dropped
        std::string seed = f.R + f.X + f.Y.substr(0, 5);
        auto res = run(*anno, seed, { "A", "B" });
        ASSERT_EQ(1u, res.dropped_labels.size());
        EXPECT_EQ("A", res.dropped_labels[0].name);
        EXPECT_EQ("seed_unsupported", res.dropped_labels[0].reason);
        uint64_t supported = f.R.size() + f.X.size() - kK + 1;
        EXPECT_EQ((std::vector<std::pair<uint64_t, uint64_t>>{ { 0, supported } }),
                  res.dropped_labels[0].runs);
        ASSERT_EQ(1u, res.label_dict.size());
        EXPECT_EQ("B", res.label_dict[0].name);
        EXPECT_EQ(1u, res.num_seed_labels);
        EXPECT_EQ(make_seed_id("", seed, mode != DeBruijnGraph::BASIC, { "B" }), res.validated_seed_id);
        EXPECT_EQ(f.Y.substr(5), spell_path(res.arms[kRight], res.arms[kRight].paths[0]));

        // a supplied seed id is checked
        {
            LabelOracle oracle(*anno);
            Seed s;
            s.sequence = seed;
            s.labels = { "B" };
            s.seed_id = "0000000000000000";
            auto r = traverse_seed(oracle, s, Strategy(), LabelChangeCost::forbid());
            EXPECT_TRUE(r.seed_id_mismatch);
            EXPECT_EQ("0000000000000000", r.seed_id);
            s.seed_id = r.validated_seed_id;
            EXPECT_FALSE(traverse_seed(oracle, s, Strategy(), LabelChangeCost::forbid()).seed_id_mismatch);
            // the release id is folded into the seed id
            EXPECT_NE(r.validated_seed_id,
                      traverse_seed(oracle, s, Strategy(), LabelChangeCost::forbid(), "rel").validated_seed_id);
        }

        // all labels dropped
        EXPECT_THROW(run(*anno, f.W + f.R, { "B" }), std::invalid_argument);
        // a k-mer absent from the graph
        EXPECT_THROW(run(*anno, random_seq(40, 99), { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R + "ACGTACGTACGTACGT", { "A", "B" }), std::invalid_argument);
        // unknown label
        EXPECT_THROW(run(*anno, f.R, { "C" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R, { "A", "C" }), std::invalid_argument);
        // short or invalid seeds
        EXPECT_THROW(run(*anno, f.R.substr(0, kK - 1), { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R.substr(0, 20) + "N" + f.R.substr(21, 20), { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, "", { "A" }), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R, { "A", "A" }), std::invalid_argument);
        // an omitted label list is NOT an error: it derives the set from the seed
        // (DeriveSeedLabels below). A seed no label carries in full still is one.
        EXPECT_EQ(2u, run(*anno, f.R, {}).num_seed_labels);
        EXPECT_THROW(run(*anno, f.W + f.R + f.X + f.Y, {}), std::invalid_argument);
        // extra labels: unreachable under forbid, or above the budget
        Strategy st;
        st.extra = { "B" };
        EXPECT_THROW(run(*anno, f.R, { "A" }, st), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R, { "A" }, st, LabelChangeCost::constant(1)), std::invalid_argument);
        EXPECT_THROW(run(*anno, f.R, { "A", "B" }, st, LabelChangeCost::constant(0)), std::invalid_argument);
        st.extra = { "C" };
        EXPECT_THROW(run(*anno, f.R, { "A" }, st, LabelChangeCost::constant(0)), std::invalid_argument);
        st.extra = { "B" };
        st.loss_budget = 1;
        auto ok = run(*anno, f.R, { "A" }, st, LabelChangeCost::constant(1));
        ASSERT_EQ(2u, ok.label_dict.size());
        EXPECT_EQ(1u, ok.num_seed_labels);
        EXPECT_EQ("B", ok.label_dict[1].name);
        // the extra label is switched into where A ends
        EXPECT_EQ(f.X.size() + f.Y.size(), ok.arms[kRight].paths[0].length_bp);
        EXPECT_EQ(1u, count_events(ok.arms[kRight], EventType::SWITCH));
        // trace support needs coordinates
        st = Strategy();
        st.support = Support::TRACE;
        EXPECT_THROW(run(*anno, f.R, { "A" }, st), std::invalid_argument);
        st = Strategy();
        st.max_live_paths = 0;
        EXPECT_THROW(run(*anno, f.R, { "A" }, st), std::invalid_argument);
    }
}


// T17c: the permitted set derived from the seed itself (the design note's
// `permit: per_hit` default). A on W·R·X, B on R·X·Y, C on the first 40 bp of R, so
// the seed R is carried in full by A and B and only partially by C.
struct DerivedFixture : public LabelEndFixture {
    std::string prefix() const { return R.substr(0, 40); }
    std::vector<std::string> seqs() const { return { W + R + X, R + X + Y, prefix() }; }
    std::vector<std::string> labels() const { return { "A", "B", "C" }; }
};

TYPED_TEST(WalkerTest, DeriveSeedLabels) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    DerivedFixture f;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, f.seqs(), f.labels(), mode);
        Strategy st;
        auto derived = run(*anno, f.R, {}, st);
        EXPECT_TRUE(derived.labels_from_seed);
        EXPECT_EQ(2u, derived.num_seed_labels) << mode;
        EXPECT_EQ(2u, derived.labels_supporting_total);
        EXPECT_EQ(0u, derived.labels_dropped);
        EXPECT_TRUE(derived.labels_dropped_digest.empty());
        // C supports only a prefix of the seed, so it is not derived at all: it was
        // never permitted, which is not the same as a named label that got dropped
        EXPECT_TRUE(derived.dropped_labels.empty());

        std::vector<std::string> names;
        for (const auto &l : derived.label_dict) {
            EXPECT_EQ(LabelKind::COLUMN, l.kind);   // no CoordToHeader: column labels
            names.push_back(l.name);
        }
        EXPECT_EQ((std::set<std::string>{ "A", "B" }),
                  std::set<std::string>(names.begin(), names.end())) << mode;
        for (size_t i = 1; i < derived.label_dict.size(); ++i) {   // ascending column id
            EXPECT_LT(derived.label_dict[i - 1].column, derived.label_dict[i].column);
        }

        // ... and the run is the one naming exactly those labels in that order
        auto named = run(*anno, f.R, names, st);
        EXPECT_FALSE(named.labels_from_seed);
        EXPECT_EQ(serialize(named), serialize(derived)) << mode;
        EXPECT_EQ(named.validated_seed_id, derived.validated_seed_id);

        // all three labels carry the 40 bp prefix of R in full
        auto all_three = run(*anno, f.prefix(), {}, st);
        EXPECT_EQ(3u, all_three.num_seed_labels) << mode;
        EXPECT_EQ(3u, all_three.labels_supporting_total);
        EXPECT_EQ(0u, all_three.labels_dropped);

        // the cap keeps the lowest two columns and reports what it cut
        st.max_seed_labels = 2;
        auto capped = run(*anno, f.prefix(), {}, st);
        ASSERT_EQ(2u, capped.num_seed_labels) << mode;
        EXPECT_EQ(3u, capped.labels_supporting_total);
        EXPECT_EQ(1u, capped.labels_dropped);
        EXPECT_TRUE(capped.dropped_labels.empty());   // capped, not unsupported
        EXPECT_EQ(capped.label_dict[0].name, all_three.label_dict[0].name);
        EXPECT_EQ(capped.label_dict[1].name, all_three.label_dict[1].name);
        const std::string &cut = all_three.label_dict[2].name;
        EXPECT_EQ(hex64(fnv1a64(cut + "\n")), capped.labels_dropped_digest) << mode;
        st.max_seed_labels = 0;   // a cap of zero would permit nothing
        EXPECT_THROW(run(*anno, f.prefix(), {}, st), std::invalid_argument);
        st.max_seed_labels = 1000;

        // no label carries W·R·X·Y in full: the derived set is empty and the seed is
        // rejected, exactly as when every named label is dropped
        EXPECT_THROW(run(*anno, f.W + f.R + f.X + f.Y, {}, st), std::invalid_argument);

        // `extra` applies on top of the derived set
        Strategy sx;
        sx.extra = { "C" };
        sx.loss_budget = 1;
        auto with_extra = run(*anno, f.R, {}, sx, LabelChangeCost::constant(1));
        ASSERT_EQ(3u, with_extra.label_dict.size()) << mode;
        EXPECT_EQ(2u, with_extra.num_seed_labels);
        EXPECT_EQ("C", with_extra.label_dict[2].name);
        EXPECT_EQ(serialize(run(*anno, f.R, names, sx, LabelChangeCost::constant(1))),
                  serialize(with_extra)) << mode;
        // an extra label the derivation also produced is a duplicate, as when named
        sx.extra = { names[0] };
        EXPECT_THROW(run(*anno, f.R, {}, sx, LabelChangeCost::constant(1)), std::invalid_argument);
    }
}


// T16: determinism
TYPED_TEST(WalkerTest, Determinism) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    ChainFixture chain;
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, chain.seqs(), chain.labels(), mode);
        Strategy st;
        st.loss_budget = 1;
        std::string first;
        for (int i = 0; i < 3; ++i) {
            auto res = run(*anno, chain.b1, { "A", "B" }, st, LabelChangeCost::constant(1));
            std::string s = serialize(res);
            if (i == 0) {
                first = s;
            } else {
                EXPECT_EQ(first, s);
            }
        }
        // the label order of the seed changes ids but not the structure
        auto a = run(*anno, chain.b1, { "A", "B" }, st, LabelChangeCost::constant(1));
        auto b = run(*anno, chain.b1, { "B", "A" }, st, LabelChangeCost::constant(1));
        EXPECT_EQ(a.validated_seed_id, b.validated_seed_id);
        EXPECT_EQ(a.arms[kRight].paths[0].length_bp, b.arms[kRight].paths[0].length_bp);
        // batch_kmers never changes which steps are taken, nor the contract counters
        // (successor enumerations, keys mapped, rows requested; §6.8, T23)
        auto base = run(*anno, chain.b1, { "A", "B" }, st, LabelChangeCost::constant(1));
        for (size_t batch : { 0, 1, 2, 7, 64 }) {
            Strategy sb = st;
            sb.batch_kmers = batch;
            auto res = run(*anno, chain.b1, { "A", "B" }, sb, LabelChangeCost::constant(1));
            EXPECT_EQ(first, serialize(res)) << "batch " << batch;
            EXPECT_EQ(base.annotation_counters.keys_mapped, res.annotation_counters.keys_mapped)
                << "batch " << batch;
            EXPECT_EQ(base.annotation_counters.rows_requested, res.annotation_counters.rows_requested)
                << "batch " << batch;
        }
    }
}


// §6.4 quorum and split limits on the fork X+P (A), X+Q (B), seed X {A, B}
TYPED_TEST(WalkerTest, QuorumAndSplitLimits) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(3);
    const std::string &X = b[0], &P = b[1], &Q = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q }, { "A", "B" }, mode);
        struct Case { const char *name; Strategy st; const char *text; };
        std::vector<Case> cases;
        Strategy st;
        st.min_successor_labels = 2;
        cases.push_back({ "min_successor_labels", st, "minority" });
        st = Strategy();
        st.min_successor_fraction = 0.6;
        cases.push_back({ "min_successor_fraction", st, "minority" });
        st = Strategy();
        st.max_splits_per_path = 0;
        cases.push_back({ "max_splits_per_path", st, "split_limit" });
        st = Strategy();
        st.min_live_labels = 2;
        cases.push_back({ "min_live_labels", st, "below_min_labels" });
        for (const Case &c : cases) {
            auto res = run(*anno, X, { "A", "B" }, c.st);
            const ArmResult &right = res.arms[kRight];
            check_invariants(right, c.st);
            ASSERT_EQ(1u, right.paths.size()) << c.name << " " << mode;
            EXPECT_EQ(0u, right.paths[0].length_bp) << c.name;
            EXPECT_EQ(2u, right.paths[0].end_reasons[static_cast<size_t>(EndReason::BRANCH)]) << c.name;
            EXPECT_TRUE(right.splits.empty()) << c.name;
            EXPECT_EQ(0u, right.steps) << c.name;
            EXPECT_EQ(2u, count_ends(right, EndReason::BRANCH)) << c.name;
            for (LabelId l : { 0u, 1u }) {
                const Event *ev = label_end_event(right, l);
                ASSERT_NE(nullptr, ev) << c.name;
                EXPECT_EQ(0u, ev->at_bp) << c.name;
                EXPECT_EQ(EndReason::BRANCH, ev->reason) << c.name;
                EXPECT_EQ(c.text, ev->text) << c.name;
                EXPECT_EQ(2u, ev->structural_successors) << c.name;
            }
            EXPECT_EQ(1u, right.branch_events_total) << c.name;
            ASSERT_EQ(1u, right.branch_events.size()) << c.name;
            const BranchEvent &be = right.branch_events[0];
            EXPECT_EQ(0u, be.at_bp);
            std::vector<char> chars { P[0], Q[0] };
            std::sort(chars.begin(), chars.end());
            EXPECT_EQ(chars, be.chars) << c.name;
            EXPECT_EQ((std::vector<size_t>{ 1, 1 }), be.labels_per_successor) << c.name;
            EXPECT_TRUE(be.ambiguous.empty()) << c.name;
            EXPECT_EQ((std::vector<LabelId>{ 0, 1 }), be.dropped) << c.name;
            EXPECT_EQ(2u, be.labels_affected) << c.name;
            // the left arm is a plain dead end for both
            EXPECT_EQ(2u, count_ends(res.arms[kLeft], EndReason::DEAD_END)) << c.name;
        }
        // thresholds that are met change nothing
        st = Strategy();
        st.min_successor_labels = 1;
        st.min_successor_fraction = 0.5;
        st.max_splits_per_path = 1;
        st.min_live_labels = 1;
        auto ok = run(*anno, X, { "A", "B" }, st);
        EXPECT_EQ(2u, ok.arms[kRight].paths.size());
        EXPECT_EQ(0u, count_ends(ok.arms[kRight], EndReason::BRANCH));
        EXPECT_EQ(0u, ok.arms[kRight].branch_events_total);
    }
}


// T15b: the frontier order decides which head advances when a per-seed cap trips.
// Seed X {A}, extra {B, C}: A on X+P (loss 0, one label), B and C on X+Q (switched
// in at cost 1, two labels). With max_steps = 3 the split uses two steps and only
// the head processed first in the next level takes the third.
TYPED_TEST(WalkerTest, FrontierOrder) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(13);
    const std::string &X = b[0];
    for (int swap = 0; swap < 2; ++swap) {
        // both successor orders: the loss-0 branch first, then the 2-label branch first
        const std::string &P = swap ? b[2] : b[1], &Q = swap ? b[1] : b[2];
        for (auto mode : all_modes()) {
            auto anno = build_anno_graph<Graph, Annotation>(kK, { X + P, X + Q, X + Q },
                                                            { "A", "B", "C" }, mode);
            Strategy st;
            st.direction = Strategy::RIGHT;
            st.extra = { "B", "C" };
            st.loss_budget = 1;
            st.switch_on_loss_only = false;
            st.max_label_branches = 1;
            st.max_steps = 3;
            for (auto order : { Strategy::BREADTH_FIRST, Strategy::LOWEST_LOSS_FIRST,
                                Strategy::MOST_SUPPORTED_FIRST }) {
                st.order = order;
                auto res = run(*anno, X, { "A" }, st, LabelChangeCost::constant(1));
                const ArmResult &right = res.arms[kRight];
                check_invariants(right, st);
                EXPECT_EQ(ArmResult::TRUNCATED, right.status);
                ASSERT_TRUE(right.cap_trigger.has_value());
                EXPECT_EQ(EndReason::MAX_STEPS, right.cap_trigger->reason);
                EXPECT_EQ(3u, right.steps);
                ASSERT_EQ(1u, right.splits.size());
                EXPECT_TRUE(right.splits[0].ambiguous);
                ASSERT_EQ(2u, right.paths.size());
                const PathResult *advanced = nullptr, *stalled = nullptr;
                for (const auto &path : right.paths) {
                    EXPECT_EQ(EndReason::MAX_STEPS, path.path_reason.value_or(EndReason::DEAD_END));
                    (path.length_bp == 2 ? advanced : stalled) = &path;
                }
                ASSERT_NE(nullptr, advanced) << mode << " order " << order;
                ASSERT_NE(nullptr, stalled) << mode << " order " << order;
                EXPECT_EQ(1u, stalled->length_bp);
                std::vector<LabelId> labels;
                for (const auto &e : advanced->end_labels) labels.push_back(e.label);
                const std::vector<LabelId> loss0 { 0 }, supported { 1, 2 };
                // BREADTH_FIRST follows the path ids, i.e. the successor character order
                std::vector<LabelId> expected;
                switch (order) {
                    case Strategy::BREADTH_FIRST: expected = P[0] < Q[0] ? loss0 : supported; break;
                    case Strategy::LOWEST_LOSS_FIRST: expected = loss0; break;
                    case Strategy::MOST_SUPPORTED_FIRST: expected = supported; break;
                }
                EXPECT_EQ(expected, labels) << mode << " order " << order << " swap " << swap;
            }
        }
    }
}


// TABLE costs: the switch sources are the max_switch_sources cheapest by (loss,
// branches, column), not the cheapest pairs; labels cut from that set are never
// reported as loss_budget, and an over-budget switch into a target another source
// carries is a plain label_lost (not "superseded")
TYPED_TEST(WalkerTest, SwitchSources) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 60, 40 }, 34);
    const std::string &S = b[0], &Z = b[1];
    // A, B, C on S; E continues S into Z
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { S, S, S, S.substr(S.size() - kK) + Z },
                                                        { "A", "B", "C", "E" }, mode);
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.extra = { "E" };
        st.loss_budget = 1;
        // request order: A 0, B 1, C 2, E 3
        auto table = LabelChangeCost::table({ { { 0, 3 }, 2.0 }, { { 1, 3 }, 0.5 } }, kInfiniteLoss);

        // one source: A (column 0) is the only candidate and its switch is over budget
        st.max_switch_sources = 1;
        auto res = run(*anno, S, { "A", "B", "C" }, st, table);
        const ArmResult &r1 = res.arms[kRight];
        check_invariants(r1, st);
        ASSERT_EQ(4u, res.label_dict.size());
        ASSERT_EQ(1u, r1.paths.size());
        EXPECT_EQ(0u, r1.paths[0].length_bp);
        EXPECT_EQ(0u, count_events(r1, EventType::SWITCH));
        EXPECT_EQ(1u, count_ends(r1, EndReason::LOSS_BUDGET));
        EXPECT_EQ(2u, count_ends(r1, EndReason::LABEL_LOST));
        EXPECT_EQ((std::vector<double>{ 2.0 }), r1.needed_budgets);
        {
            const Event *a = label_end_event(r1, 0), *bb = label_end_event(r1, 1), *c = label_end_event(r1, 2);
            ASSERT_TRUE(a && bb && c);
            EXPECT_EQ(EndReason::LOSS_BUDGET, a->reason);
            EXPECT_EQ(2.0, a->needed_budget);
            EXPECT_EQ(EndReason::LABEL_LOST, bb->reason);
            EXPECT_EQ("switch_sources", bb->text) << mode;
            EXPECT_EQ(EndReason::LABEL_LOST, c->reason);
            EXPECT_EQ("", c->text);
        }
        // the cut of B (whose B -> E is priced) is counted as a derivation it may have
        // changed, as well as by B's label end
        EXPECT_EQ(1u, r1.switch_sources_cut) << mode;

        // two sources: B's switch is taken, A's over-budget switch is a plain loss
        st.max_switch_sources = 2;
        res = run(*anno, S, { "A", "B", "C" }, st, table);
        const ArmResult &r2 = res.arms[kRight];
        check_invariants(r2, st);
        ASSERT_EQ(1u, r2.paths.size());
        EXPECT_EQ(Z.size(), r2.paths[0].length_bp);
        EXPECT_EQ(S.substr(S.size() - kK) + Z, S.substr(S.size() - kK) + spell_path(r2, r2.paths[0]));
        ASSERT_EQ(1u, r2.paths[0].end_labels.size());
        EXPECT_EQ(3u, r2.paths[0].end_labels[0].label);
        EXPECT_EQ(0.5, r2.paths[0].end_labels[0].loss);
        ASSERT_EQ(1u, count_events(r2, EventType::SWITCH));
        EXPECT_EQ(1u, events_of(r2, EventType::SWITCH)[0]->label);
        EXPECT_EQ(3u, events_of(r2, EventType::SWITCH)[0]->to);
        EXPECT_EQ(0.5, events_of(r2, EventType::SWITCH)[0]->cost);
        EXPECT_EQ(0u, count_ends(r2, EndReason::LOSS_BUDGET));
        EXPECT_TRUE(r2.needed_budgets.empty());
        {
            const Event *a = label_end_event(r2, 0), *c = label_end_event(r2, 2);
            ASSERT_TRUE(a && c);
            EXPECT_EQ(EndReason::LABEL_LOST, a->reason);
            EXPECT_EQ("", a->text) << mode;
            EXPECT_EQ(EndReason::LABEL_LOST, c->reason);
            EXPECT_EQ("", c->text);
        }
        // the default (64) behaves like 2 here
        st.max_switch_sources = 64;
        auto r64 = run(*anno, S, { "A", "B", "C" }, st, table);
        EXPECT_EQ(serialize(res), serialize(r64));
        // with two sources the only cut one is C, which has no finite switch into E, so
        // nothing that could matter was cut (the one-source run is counted above)
        EXPECT_EQ(0u, r2.switch_sources_cut) << mode;
        EXPECT_EQ(0u, r64.arms[kRight].switch_sources_cut) << mode;
    }
}

// A cut switch-source list that leaves NO label end behind (spec §6.3, §7.0): the cut
// source goes on along another successor, so nothing ends with `switch_sources`, yet the
// target it was the cheapest way into is entered at a higher loss. A, B, C on S; past S
// the graph forks into Z1, where B goes on, and Z2, carried only by E (extra). B -> E
// costs 0.5, A -> E 0.8. With one source A (the lower column) is the only one priced on
// Z2, so E is entered at 0.8, and no label ends because of the cut (B goes on along Z1,
// C has no finite switch). The arm's switch_sources_cut counts the derivation and the
// response states it; "unlimited" prices B, enters E at 0.5 (B then follows both
// successors, so it needs one branch) and states nothing.
TYPED_TEST(WalkerTest, SwitchSourcesCutWithoutALabelEnd) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = fork_blocks(34, { 60, 40, 40 });
    const std::string &S = b[0], &Z1 = b[1], &Z2 = b[2];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(
                kK, { S, S + Z1, S, S.substr(S.size() - kK) + Z2 }, { "A", "B", "C", "E" }, mode);
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.extra = { "E" };
        st.loss_budget = 1;
        st.max_label_branches = 1;
        st.merge_reconverge = false;     // no `scope` limitation beside the one tested
        // request order: A 0, B 1, C 2, E 3
        auto table = LabelChangeCost::table({ { { 0, 3 }, 0.8 }, { { 1, 3 }, 0.5 } }, kInfiniteLoss);
        auto loss_of_e = [&](const ArmResult &arm) {
            for (const PathResult &p : arm.paths) {
                for (const LabelEnd &e : p.end_labels) {
                    if (e.label == 3)
                        return e.loss;
                }
            }
            return -1.0;
        };

        st.max_switch_sources = 1;
        auto cut = run(*anno, S, { "A", "B", "C" }, st, table);
        const ArmResult &rc1 = cut.arms[kRight];
        check_invariants(rc1, st);
        ASSERT_EQ(2u, rc1.paths.size()) << mode;
        for (const PathResult &p : rc1.paths) {
            EXPECT_EQ(Z1.size(), p.length_bp) << mode;
        }
        EXPECT_EQ(0.8, loss_of_e(rc1)) << mode;   // overestimated: B's 0.5 was cut
        for (const Event *ev : events_of(rc1, EventType::LABEL_END)) {
            EXPECT_NE("switch_sources", ev->text) << mode << ": the cut source did not end";
        }
        EXPECT_EQ(1u, rc1.switch_sources_cut) << mode;
        // every walk is there; a loss in it may be too high
        EXPECT_EQ("complete/complete/lower_bound/inline",
                  outcome_of(cli::seed_result_to_json(cut, st, "summary", false))) << mode;
        Json::Value arm = cli::seed_result_to_json(cut, st, "summary", false)["arms"]["right"];
        ASSERT_EQ(1u, arm["limitations"].size()) << mode;
        const Json::Value &l = arm["limitations"][0];
        EXPECT_EQ("switch_sources", l["kind"].asString());
        EXPECT_EQ("labels.max_switch_sources", l["knob"].asString());
        EXPECT_EQ(1u, l["limit"].asUInt64());
        EXPECT_EQ(1u, l["observed"].asUInt64());       // the derivations
        EXPECT_EQ(0u, l["label_ends"].asUInt64());     // ... none of which ended a label
        EXPECT_NE(std::string::npos, l["effect"].asString().find("overestimated"));
        EXPECT_EQ(1u, arm["counters"]["switch_sources_cut"].asUInt64());

        st.max_switch_sources = Strategy::kUnlimited;
        auto full = run(*anno, S, { "A", "B", "C" }, st, table);
        const ArmResult &ru = full.arms[kRight];
        check_invariants(ru, st);
        EXPECT_EQ(2u, ru.paths.size()) << mode;
        EXPECT_EQ(0.5, loss_of_e(ru)) << mode;
        EXPECT_EQ(0u, ru.switch_sources_cut) << mode;
        EXPECT_EQ("complete/complete/complete/inline",
                  outcome_of(cli::seed_result_to_json(full, st, "summary", false))) << mode;
        arm = cli::seed_result_to_json(full, st, "summary", false)["arms"]["right"];
        EXPECT_EQ(0u, arm["limitations"].size()) << mode;
        EXPECT_EQ(0u, arm["counters"]["switch_sources_cut"].asUInt64());
    }
}


// TABLE indices are request indices: a dropped seed label does not renumber the
// pairs (A: W·R·X is dropped on the seed R·X·Y[0:5]; B: R·X·Y; extra C: Y·Z)
TYPED_TEST(WalkerTest, TableRequestOrder) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    auto b = clean_blocks({ 20, 60, 25, 30, 30 }, 35);
    const std::string &W = b[0], &R = b[1], &X = b[2], &Y = b[3], &Z = b[4];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<Graph, Annotation>(kK, { W + R + X, R + X + Y, Y + Z },
                                                        { "A", "B", "C" }, mode);
        std::string seed = R + X + Y.substr(0, 5);
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.extra = { "C" };
        st.loss_budget = 1;
        // request order: A 0, B 1, C 2; B -> C priced, everything else forbidden
        auto table = LabelChangeCost::table({ { { 1, 2 }, 0.5 } }, kInfiniteLoss);
        auto res = run(*anno, seed, { "A", "B" }, st, table);
        ASSERT_EQ(1u, res.dropped_labels.size());
        EXPECT_EQ("A", res.dropped_labels[0].name);
        ASSERT_EQ(2u, res.label_dict.size());
        EXPECT_EQ("B", res.label_dict[0].name);
        EXPECT_EQ("C", res.label_dict[1].name);
        const ArmResult &right = res.arms[kRight];
        check_invariants(right, st);
        ASSERT_EQ(1u, right.paths.size());
        EXPECT_EQ(Y.size() - 5 + Z.size(), right.paths[0].length_bp) << mode;
        EXPECT_EQ(Y.substr(5) + Z, spell_path(right, right.paths[0]));
        ASSERT_EQ(1u, count_events(right, EventType::SWITCH));
        const Event *sw = events_of(right, EventType::SWITCH)[0];
        EXPECT_EQ(Y.size() - 5, sw->at_bp);
        EXPECT_EQ(0u, sw->label);     // dictionary ids in the output
        EXPECT_EQ(1u, sw->to);
        EXPECT_EQ(0.5, sw->cost);
        ASSERT_EQ(1u, right.paths[0].end_labels.size());
        EXPECT_EQ(0.5, right.paths[0].end_labels[0].loss);

        // the same request without the dropped label, priced in its own index space
        auto same = run(*anno, seed, { "B" }, st, LabelChangeCost::table({ { { 0, 1 }, 0.5 } }, kInfiniteLoss));
        EXPECT_EQ(right.paths[0].length_bp, same.arms[kRight].paths[0].length_bp);
        EXPECT_EQ(res.validated_seed_id, same.validated_seed_id);
        EXPECT_EQ(serialize(res).substr(serialize(res).find('\n')),
                  serialize(same).substr(serialize(same).find('\n')));

        // an entry naming the dropped label is ignored: C becomes unreachable
        auto dropped_pair = LabelChangeCost::table({ { { 0, 2 }, 0.5 } }, kInfiniteLoss);
        EXPECT_THROW(run(*anno, seed, { "A", "B" }, st, dropped_pair), std::invalid_argument);
    }
}


// T22: reference walker. A naive per-label walk applying the forbid / limit-0 rules.
struct RefEnd {
    std::string label;
    Arm arm;
    uint64_t at_bp;
    EndReason reason;
    bool operator<(const RefEnd &o) const {
        return std::tie(label, arm, at_bp, reason) < std::tie(o.label, o.arm, o.at_bp, o.reason);
    }
    bool operator==(const RefEnd &o) const {
        return std::tie(label, arm, at_bp, reason) == std::tie(o.label, o.arm, o.at_bp, o.reason);
    }
};

RefEnd reference_walk(const LabelOracle &oracle, const std::string &label, Arm arm,
                      const std::string &seed) {
    const DeBruijnGraph &graph = oracle.graph();
    const size_t k = oracle.get_k();
    const bool canonical = oracle.regime() != Regime::BASIC;
    LabelQuery query(oracle, { oracle.resolve_label(label) }, false);
    std::set<node_index> seed_nodes;
    for (node_index n : map_to_nodes_sequentially(graph, seed)) seed_nodes.insert(n);
    if (canonical) {
        for (node_index n : map_to_nodes_sequentially(graph, rc(seed))) {
            if (n != npos) seed_nodes.insert(n);
        }
    }
    auto nodes = map_to_nodes_sequentially(graph, seed);
    node_index node = arm == Arm::RIGHT ? nodes.back() : nodes.front();
    std::string kmer = arm == Arm::RIGHT ? seed.substr(seed.size() - k) : seed.substr(0, k);
    std::map<std::string, std::set<bool>> used;
    uint64_t at = 0;
    while (true) {
        std::vector<std::pair<node_index, char>> succs;
        auto cb = [&](node_index n, char c) { if (c != '$') succs.emplace_back(n, c); };
        if (arm == Arm::RIGHT) {
            graph.call_outgoing_kmers(node, cb);
        } else {
            graph.call_incoming_kmers(node, cb);
        }
        std::sort(succs.begin(), succs.end(), [](const auto &a, const auto &b) { return a.second < b.second; });
        if (succs.empty())
            return { label, arm, at, EndReason::DEAD_END };
        std::vector<std::pair<node_index, std::string>> admissible;
        std::vector<std::string> admissible_steps;
        uint8_t blocked = 0;
        bool hairpin = false;
        for (const auto &[v, c] : succs) {
            std::string vk = arm == Arm::RIGHT ? kmer.substr(1) + c : c + kmer.substr(0, k - 1);
            if (query.fetch(oracle.key_of(v, vk)).empty())
                continue;
            std::string step = arm == Arm::RIGHT ? kmer + c : c + kmer;
            if (canonical && rc(step) == step) {
                hairpin = true;
                continue;
            }
            if (seed_nodes.count(v)) {
                blocked = std::max<uint8_t>(blocked, 3);
                continue;
            }
            std::string key = canonical ? std::min(step, rc(step)) : step;
            bool orientation = key != step;
            auto it = used.find(key);
            if (it != used.end()) {
                if (it->second.count(!orientation)) {
                    blocked = std::max<uint8_t>(blocked, 2);
                } else {
                    blocked = std::max<uint8_t>(blocked, 1);
                }
                continue;
            }
            admissible.emplace_back(v, vk);
            admissible_steps.push_back(step);
        }
        if (admissible.size() >= 2)
            return { label, arm, at, EndReason::BRANCH };
        if (admissible.empty()) {
            if (blocked == 3) return { label, arm, at, EndReason::REACHED_SEED };
            if (blocked == 2) return { label, arm, at, EndReason::EDGE_REUSE_RC };
            if (blocked == 1) return { label, arm, at, EndReason::EDGE_REUSE };
            if (hairpin) return { label, arm, at, EndReason::DEAD_END };
            return { label, arm, at, EndReason::LABEL_LOST };
        }
        const std::string &step = admissible_steps[0];
        std::string key = canonical ? std::min(step, rc(step)) : step;
        used[key].insert(key != step);
        node = admissible[0].first;
        kmer = admissible[0].second;
        ++at;
    }
}

TYPED_TEST(WalkerTest, ReferenceWalker) {
    using Graph = typename TypeParam::first_type;
    using Annotation = typename TypeParam::second_type;
    std::mt19937 gen(22);
    for (size_t fixture = 0; fixture < 30; ++fixture) {
        // a pool of blocks; every label contains the seed block B0 at least once
        std::vector<size_t> lengths;
        for (size_t i = 0; i < 6; ++i) lengths.push_back(15 + gen() % 26);
        lengths[0] = 30;
        auto blocks = clean_blocks(lengths, 1000 + fixture);
        std::vector<std::string> seqs, labels { "A", "B", "C" };
        for (size_t l = 0; l < 3; ++l) {
            size_t n = 3 + gen() % 3;
            std::vector<size_t> order;
            for (size_t i = 0; i < n; ++i) order.push_back(gen() % blocks.size());
            order[gen() % n] = 0;
            std::string s;
            for (size_t i : order) s += blocks[i];
            seqs.push_back(s);
        }
        for (auto mode : all_modes()) {
            auto anno = build_anno_graph<Graph, Annotation>(kK, seqs, labels, mode);
            Strategy st = strategy(0, false);
            auto res = run(*anno, blocks[0], labels, st);
            ASSERT_EQ(3u, res.label_dict.size()) << "fixture " << fixture;
            std::set<RefEnd> got, expected;
            for (size_t a : { kLeft, kRight }) {
                const ArmResult &arm = res.arms[a];
                check_invariants(arm, st);
                ASSERT_EQ(ArmResult::COMPLETE, arm.status);
                for (const auto &run : arm.runs) {
                    ASSERT_TRUE(run.ended);
                    got.insert({ res.label_dict[run.label].name, arm.arm, run.to_bp, run.end_reason });
                }
                for (const auto &path : arm.paths) {
                    check_spelled(*anno, arm, path, blocks[0]);
                }
            }
            LabelOracle oracle(*anno);
            for (const auto &label : labels) {
                for (Arm arm : { Arm::LEFT, Arm::RIGHT }) {
                    expected.insert(reference_walk(oracle, label, arm, blocks[0]));
                }
            }
            EXPECT_EQ(expected.size(), got.size()) << "fixture " << fixture << " mode " << mode;
            for (const auto &e : expected) {
                EXPECT_TRUE(got.count(e)) << "fixture " << fixture << " mode " << mode << ": "
                    << e.label << " " << to_string(e.arm) << " " << e.at_bp << " " << to_string(e.reason);
            }
            for (const auto &g : got) {
                EXPECT_TRUE(expected.count(g)) << "fixture " << fixture << " mode " << mode << ": "
                    << g.label << " " << to_string(g.arm) << " " << g.at_bp << " " << to_string(g.reason);
            }
        }
    }
}


// Trace support on a coordinate index (the refseq33m shape): a repeat inside one
// accession is followed under k-mer support and kept apart under trace support.
TEST(WalkerCoord, TraceSupport) {
    // acc1 = L·R·M·R (R twice), acc2 = R, one column "F". L and M must end with
    // different bases so that the two entries into R are distinct nodes.
    std::vector<std::string> b;
    for (uint32_t seed = 29; ; ++seed) {
        b = clean_blocks({ 40, 40, 40 }, seed);
        if (b[0].back() != b[2].back())
            break;
    }
    const std::string &L = b[0], &R = b[1], &M = b[2];
    std::string acc1 = L + R + M + R, acc2 = R;
    uint64_t n1 = acc1.size() - kK + 1, n2 = acc2.size() - kK + 1;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { acc1, acc2 }, { "F", "F" }, DeBruijnGraph::BASIC, true, { 0, n1 });
    std::vector<std::vector<std::string>> headers { { "acc1", "acc2" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n1, n2 } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);
    ASSERT_TRUE(oracle.has_coordinates());

    for (const std::string &label : { "acc1", "F" }) {
        Seed seed;
        seed.sequence = M;
        seed.labels = { label };

        // k-mer support: the right arm walks R and re-enters the seed through M,
        // the left arm reaches the ambiguous junction before R
        Strategy st;
        auto res = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        const ArmResult &rk = res.arms[kRight], &lk = res.arms[kLeft];
        check_invariants(rk, st);
        ASSERT_EQ(1u, rk.paths.size());
        EXPECT_EQ(R.size() + kK - 1, rk.paths[0].length_bp);
        EXPECT_EQ(1u, rk.paths[0].end_reasons[static_cast<size_t>(EndReason::REACHED_SEED)]);
        ASSERT_EQ(1u, lk.paths.size());
        EXPECT_EQ(R, spell_path(lk, lk.paths[0]));
        EXPECT_EQ(R.size(), lk.paths[0].length_bp);
        EXPECT_EQ(1u, lk.paths[0].end_reasons[static_cast<size_t>(EndReason::BRANCH)]);

        // trace support: the second R ends the record; the left arm follows the
        // coordinates back through the first R into L
        st.support = Support::TRACE;
        st.merge_reconverge = false;   // trace evidence cannot survive a merge
        res = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        EXPECT_STREQ("tuples", res.access_path);
        const ArmResult &rt = res.arms[kRight], &lt = res.arms[kLeft];
        check_invariants(rt, st);
        check_invariants(lt, st);
        ASSERT_EQ(1u, rt.paths.size()) << label;
        EXPECT_EQ(R.size(), rt.paths[0].length_bp) << label;
        EXPECT_EQ(R, spell_path(rt, rt.paths[0]));
        EXPECT_EQ(1u, rt.paths[0].end_reasons[static_cast<size_t>(EndReason::RECORD_END)]);
        ASSERT_EQ(1u, lt.paths.size());
        EXPECT_EQ(L.size() + R.size(), lt.paths[0].length_bp) << label;
        EXPECT_EQ(L + R, spell_path(lt, lt.paths[0]));
        EXPECT_EQ(1u, lt.paths[0].end_reasons[static_cast<size_t>(EndReason::DEAD_END)]);
        EXPECT_EQ(0u, count_ends(lt, EndReason::BRANCH));
        EXPECT_TRUE(lt.splits.empty());
    }

    // a seed whose k-mers are not trace-consistent for the label is rejected
    Seed jump;
    jump.sequence = R + M.substr(0, 20);
    jump.labels = { "acc2" };
    Strategy st;
    st.support = Support::TRACE;
    st.merge_reconverge = false;
    EXPECT_THROW(traverse_seed(oracle, jump, st, LabelChangeCost::forbid()), std::invalid_argument);
    jump.labels = { "acc1", "acc2" };
    auto res = traverse_seed(oracle, jump, st, LabelChangeCost::forbid());
    ASSERT_EQ(1u, res.dropped_labels.size());
    EXPECT_EQ("acc2", res.dropped_labels[0].name);
    EXPECT_EQ(1u, res.label_dict.size());
}

// Deriving the permitted set on a coordinate index (the refseq33m shape): with a
// CoordToHeader the derived labels are the indexed sequences (accessions), without one
// (or on request) the column; and under `support: trace` a derived label must also be
// coordinate-consecutive over the seed.
TEST(WalkerCoord, DeriveSeedLabels) {
    // acc1 = L·R·M·R, acc2 = R, acc3 = M·R·M, one column "F". The seed M·R·M occurs
    // verbatim in acc3; acc1 carries every one of its k-mers but not consecutively
    // (its M·R is at the end of the record), and acc2 carries no M at all.
    std::vector<std::string> b;
    for (uint32_t seed = 29; ; ++seed) {
        b = clean_blocks({ 40, 40, 40 }, seed);
        if (b[0].back() != b[2].back())
            break;
    }
    const std::string &L = b[0], &R = b[1], &M = b[2];
    std::string acc1 = L + R + M + R, acc2 = R, acc3 = M + R + M;
    uint64_t n1 = acc1.size() - kK + 1, n2 = acc2.size() - kK + 1, n3 = acc3.size() - kK + 1;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { acc1, acc2, acc3 }, { "F", "F", "F" }, DeBruijnGraph::BASIC, true,
        { 0, n1, n1 + n2 });
    std::vector<std::vector<std::string>> headers { { "acc1", "acc2", "acc3" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n1, n2, n3 } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    // ---- header kind (the default when the index has a CoordToHeader)
    {
        LabelOracle oracle(*anno, &cth);
        Seed seed;
        seed.sequence = M;
        Strategy st;
        auto derived = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        EXPECT_TRUE(derived.labels_from_seed);
        ASSERT_EQ(2u, derived.num_seed_labels);
        EXPECT_EQ(LabelKind::HEADER, derived.label_dict[0].kind);
        EXPECT_EQ("acc1", derived.label_dict[0].name);   // ascending seq_id
        EXPECT_EQ("acc3", derived.label_dict[1].name);
        EXPECT_EQ(2u, derived.labels_supporting_total);

        // the seed R is carried by all three records
        seed.sequence = R;
        auto all_three = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        ASSERT_EQ(3u, all_three.num_seed_labels);
        EXPECT_EQ("acc2", all_three.label_dict[1].name);
    }
    // ---- column kind on request: the derived label is the column
    {
        LabelOracle oracle(*anno, &cth);
        Seed seed;
        seed.sequence = M;
        Strategy st;
        st.seed_label_kind = LabelKind::COLUMN;
        auto derived = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        EXPECT_TRUE(derived.labels_from_seed);
        ASSERT_EQ(1u, derived.num_seed_labels);
        EXPECT_EQ(LabelKind::COLUMN, derived.label_dict[0].kind);
        EXPECT_EQ("F", derived.label_dict[0].name);
        // identical to naming the column
        Seed named = seed;
        named.labels = { "F" };
        LabelOracle other(*anno, &cth);
        EXPECT_EQ(serialize(traverse_seed(other, named, st, LabelChangeCost::forbid())),
                  serialize(derived));
    }
    // ---- without a CoordToHeader the derived kind falls back to the column
    {
        LabelOracle oracle(*anno);
        Seed seed;
        seed.sequence = M;
        auto derived = traverse_seed(oracle, seed, Strategy(), LabelChangeCost::forbid());
        ASSERT_EQ(1u, derived.num_seed_labels);
        EXPECT_EQ(LabelKind::COLUMN, derived.label_dict[0].kind);
        // header labels cannot be derived without the sidecar: rejected, not downgraded
        Strategy st;
        st.seed_label_kind = LabelKind::HEADER;
        EXPECT_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()),
                     std::invalid_argument);
    }
    // ---- trace: derived from k-mer presence, then held to coordinate continuity
    {
        LabelOracle oracle(*anno, &cth);
        Seed seed;
        seed.sequence = M + R + M;
        Strategy st;
        auto by_kmer = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        ASSERT_EQ(2u, by_kmer.num_seed_labels);     // both records carry every k-mer
        EXPECT_EQ("acc1", by_kmer.label_dict[0].name);
        EXPECT_EQ("acc3", by_kmer.label_dict[1].name);

        st.support = Support::TRACE;
        st.merge_reconverge = false;                // trace evidence cannot survive a merge
        auto by_trace = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        EXPECT_TRUE(by_trace.labels_from_seed);
        // the derived set is still the two presence carriers, but acc1 has no
        // consecutive occurrence of the seed and is dropped by the seed validation
        EXPECT_EQ(2u, by_trace.labels_supporting_total);
        ASSERT_EQ(1u, by_trace.num_seed_labels);
        EXPECT_EQ("acc3", by_trace.label_dict[0].name);
        ASSERT_EQ(1u, by_trace.dropped_labels.size());
        EXPECT_EQ("acc1", by_trace.dropped_labels[0].name);
        EXPECT_EQ("seed_unsupported", by_trace.dropped_labels[0].reason);
        // the cap applies before the trace check: with one label taken, acc1 (the first
        // by seq_id) is checked and fails while acc3, which passes, was cut unchecked.
        // The failure says so: two presence carriers, one of them cut by the cap
        Strategy capped = st;
        capped.max_seed_labels = 1;
        try {
            traverse_seed(oracle, seed, capped, LabelChangeCost::forbid());
            FAIL() << "expected no trace carrier among the labels taken";
        } catch (const SeedDerivationError &e) {
            EXPECT_EQ(SeedDerivationError::NO_TRACE_CARRIER, e.cause()) << e.what();
            EXPECT_EQ(2.0, e.observed());
            EXPECT_EQ(1u, e.labels_cut());
        }
        // the same under column labels, whose coordinates live in the column frame
        Strategy sc = st;
        sc.seed_label_kind = LabelKind::COLUMN;
        Seed one;
        one.sequence = M;
        auto column_trace = traverse_seed(oracle, one, sc, LabelChangeCost::forbid());
        ASSERT_EQ(1u, column_trace.num_seed_labels);
        EXPECT_EQ("F", column_trace.label_dict[0].name);
        Seed column_named = one;
        column_named.labels = { "F" };
        LabelOracle fresh(*anno, &cth);
        EXPECT_EQ(serialize(traverse_seed(fresh, column_named, sc, LabelChangeCost::forbid())),
                  serialize(column_trace));

        // identical to naming the trace-consistent record
        Seed named = seed;
        named.labels = { "acc3" };
        LabelOracle other(*anno, &cth);
        auto explicitly = traverse_seed(other, named, st, LabelChangeCost::forbid());
        EXPECT_EQ(serialize(explicitly).substr(serialize(explicitly).find('\n')),
                  serialize(by_trace).substr(serialize(by_trace).find('\n')));
        EXPECT_EQ(explicitly.validated_seed_id, by_trace.validated_seed_id);
    }
}

// Trace support across a reconvergence merge: one column F with S·X·R·T and
// S·Y·R·V (|X| = |Y|), seed S. The two traces of F reconverge at R; the merged
// state keeps the live coordinates of both, so both exits of R stay visible.
TEST(WalkerCoord, TraceMerge) {
    std::vector<std::string> b;
    for (uint32_t seed = 41; ; ++seed) {
        b = clean_blocks({ 30, 25, 25, 40, 30, 30 }, seed);
        if (b[1][0] != b[2][0] && b[1].back() != b[2].back() && b[4][0] != b[5][0])
            break;
    }
    const std::string &S = b[0], &X = b[1], &Y = b[2], &R = b[3], &T = b[4], &V = b[5];
    std::string acc1 = S + X + R + T, acc2 = S + Y + R + V;
    uint64_t n1 = acc1.size() - kK + 1;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { acc1, acc2 }, { "F", "F" }, DeBruijnGraph::BASIC, true, { 0, n1 });
    LabelOracle oracle(*anno);
    ASSERT_TRUE(oracle.has_coordinates());
    Seed seed;
    seed.sequence = S;
    seed.labels = { "F" };
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.support = Support::TRACE;
    st.max_label_branches = 1;

    // no merging: each trace is followed to its own end
    st.merge_reconverge = false;
    auto res = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    const ArmResult &rk = res.arms[kRight];
    check_invariants(rk, st);
    ASSERT_EQ(2u, rk.paths.size());
    std::set<std::string> flanks;
    for (const auto &path : rk.paths) flanks.insert(spell_path(rk, path));
    EXPECT_EQ((std::set<std::string>{ X + R + T, Y + R + V }), flanks);
    EXPECT_EQ(0u, count_ends(rk, EndReason::BRANCH));
    EXPECT_EQ(2u, count_ends(rk, EndReason::DEAD_END));

    // Merging is refused under trace support. A merge unites the two parents'
    // coordinate sets but keeps only one parent's (loss, branches, run, route), so a
    // later step can continue on the discarded occurrence's coordinates while carrying
    // the retained one's evidence — reporting direct support, a loss or a run that no
    // single occurrence justifies. This fixture is exactly that shape: acc1 = S·X·R·T
    // and acc2 = S·Y·R·V reconverge on R.
    st.merge_reconverge = true;
    EXPECT_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()),
                 std::invalid_argument);
    try {
        traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        FAIL() << "expected the trace+merge combination to be rejected";
    } catch (const std::invalid_argument &e) {
        EXPECT_NE(nullptr, std::strstr(e.what(), "on_reconverge")) << e.what();
    }
    // k-mer support is unaffected: merging stays available there
    st.support = Support::KMER;
    EXPECT_NO_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()));
}


// A permitted set DERIVED from the seed is chosen by the machine, not typed by the
// caller, so its size and its cost are the server's problem: the per-label validation
// must not be quadratic in it, the derivation must watch the clock it runs before, the
// derived names must stay usable as an explicit list, and which seed k-mer seeds the
// intersection must not be an accident of where the caller cut the seed.

// Many derived labels: the validation (a binary search over the hits of each k-mer, or
// nothing at all when the set is derived without `trace`) must give exactly what naming
// the same labels gives, including what it drops.
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
    // one more label carries only a prefix of the seed
    seqs.push_back(S.substr(0, S.size() - 10));
    labels.push_back("partial");
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, labels);

    Strategy st;
    auto derived = run(*anno, S, {}, st);
    EXPECT_TRUE(derived.labels_from_seed);
    ASSERT_EQ(kN, derived.num_seed_labels);
    EXPECT_EQ(kN, derived.labels_supporting_total);
    EXPECT_EQ(0u, derived.labels_dropped);
    EXPECT_TRUE(derived.dropped_labels.empty());

    std::vector<std::string> names;
    for (size_t i = 0; i < derived.num_seed_labels; ++i) {
        names.push_back(derived.label_dict[i].name);
    }
    auto named = run(*anno, S, names, st);
    EXPECT_FALSE(named.labels_from_seed);
    EXPECT_EQ(serialize(named), serialize(derived));
    // an explicit list is its own support total, never 0 (a client reads
    // labels_supporting_total > |labels| as "truncated")
    EXPECT_EQ(kN, named.labels_supporting_total);

    // the prefix carrier is dropped when NAMED, and never derived in the first place
    std::vector<std::string> with_partial = names;
    with_partial.push_back("partial");
    auto mixed = run(*anno, S, with_partial, st);
    EXPECT_EQ(kN, mixed.num_seed_labels);
    EXPECT_EQ(kN, mixed.labels_supporting_total);
    ASSERT_EQ(1u, mixed.dropped_labels.size());
    EXPECT_EQ("partial", mixed.dropped_labels[0].name);
    EXPECT_EQ("seed_unsupported", mixed.dropped_labels[0].reason);
    EXPECT_FALSE(mixed.dropped_labels[0].runs.empty());
}

// The derivation reads one FULL annotation row per seed k-mer and runs entirely before
// the walk loop, i.e. before the only other place that looks at the clock. So it has to
// watch `bounds.time_budget_ms` itself.
TEST(WalkerDerive, DerivationHonoursTheTimeBudget) {
    DerivedFixture f;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, f.seqs(), f.labels(), DeBruijnGraph::BASIC);
    Strategy st;
    st.time_budget_ms = 1e-9;    // spent by the time the first seed k-mer is consumed
    EXPECT_THROW(run(*anno, f.R, {}, st), SeedDerivationError);
    try {
        run(*anno, f.R, {}, st);
        FAIL() << "expected the derivation to stop at the time budget";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(nullptr, std::strstr(e.what(), "time budget")) << e.what();
        EXPECT_EQ(SeedDerivationError::TIME_BUDGET, e.cause());
        EXPECT_EQ(1e-9, e.limit());
        EXPECT_GT(e.observed(), e.limit());   // the elapsed milliseconds
    }

    // An explicit list reads only its own columns, so the budget stops the WALK instead
    // of rejecting the seed: the same budget must not turn a named request into an error.
    SeedResult named = run(*anno, f.R, { "A", "B" }, st);
    EXPECT_EQ(2u, named.num_seed_labels);
    EXPECT_EQ(ArmResult::TRUNCATED, named.arms[kRight].status);
}

// Which k-mer seeds the intersection decides the derivation's peak work: that row is the
// only one no intersection has narrowed yet. Fixture: 100 decoy records are exactly the
// seed's FIRST k-mer and two records carry the whole seed, so row 0 has 102 candidates
// while every other row of the seed has two. Starting from k-mer 0 costs one map_coord
// per candidate of row 0; starting from the cheapest row of the first batch costs two
// per k-mer, and the derived set is the same either way.
TEST(WalkerDerive, IntersectionStartsFromTheCheapestRowOfTheFirstBatch) {
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

    // one indexed sequence per column, named after its label; the column ids are the
    // label encoder's, not the input order
    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(labels.size());
    std::vector<std::vector<uint64_t>> num_kmers(labels.size());
    for (size_t i = 0; i < labels.size(); ++i) {
        Column c = static_cast<Column>(encoder.encode(labels[i]));
        headers[c] = { "h_" + labels[i] };
        num_kmers[c] = { seqs[i].size() - kK + 1 };
    }
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno, &cth);
    Seed seed;
    seed.sequence = seed_seq;
    auto derived = traverse_seed(oracle, seed, Strategy(), LabelChangeCost::forbid());
    ASSERT_EQ(2u, derived.num_seed_labels);
    EXPECT_EQ(2u, derived.labels_supporting_total);
    std::set<std::string> names;
    for (const auto &l : derived.label_dict) {
        EXPECT_EQ(LabelKind::HEADER, l.kind);
        names.insert(l.name);
    }
    EXPECT_EQ((std::set<std::string>{ "h_carrierA", "h_carrierB" }), names);

    // every seed k-mer costs one map_coord per live candidate (two), plus whatever the
    // (immediately dead-ending) walk needs. Starting from k-mer 0 would alone cost
    // kDecoys + 2 for that one row.
    EXPECT_LT(derived.annotation_counters.coords_mapped, kDecoys)
            << "the intersection was seeded from a row it did not have to pay for";
}

// A FASTA header is unique inside a column, not inside the index: the same accession can
// sit in two columns, and the derived set — deduplicated by (column, seq_id) — keeps
// both. That list is not resubmittable (the explicit path rejects duplicate names, and a
// name resolves to the first column holding it), so it is refused instead of echoed.
TEST(WalkerDerive, AmbiguousDerivedHeaderNamesAreRefused) {
    auto b = clean_blocks({ 30, 20 }, 211);
    const std::string S = b[0], T = b[1];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { S + T, S }, { "colA", "colB" }, DeBruijnGraph::BASIC, true);

    const auto &encoder = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[static_cast<Column>(encoder.encode("colA"))] = { "ACC1" };
    num_kmers[static_cast<Column>(encoder.encode("colA"))] = { (S + T).size() - kK + 1 };
    headers[static_cast<Column>(encoder.encode("colB"))] = { "ACC1" };
    num_kmers[static_cast<Column>(encoder.encode("colB"))] = { S.size() - kK + 1 };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    LabelOracle oracle(*anno, &cth);
    Seed seed;
    seed.sequence = S;                  // carried in full by both columns
    Strategy st;
    EXPECT_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()),
                 SeedDerivationError);
    try {
        traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        FAIL() << "expected an ambiguous derived header set to be refused";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(nullptr, std::strstr(e.what(), "ACC1")) << e.what();
        EXPECT_NE(nullptr, std::strstr(e.what(), "seed_label_kind")) << e.what();
        EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
        EXPECT_EQ("ACC1", e.subject());
    }
    // the column kind is unambiguous here: two columns, two distinct names
    st.seed_label_kind = LabelKind::COLUMN;
    auto by_column = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_EQ(2u, by_column.num_seed_labels);
    EXPECT_EQ("colA", by_column.label_dict[0].name);
    EXPECT_EQ("colB", by_column.label_dict[1].name);
    // naming the ambiguous accession explicitly still works: it resolves to ONE column
    Seed explicitly;
    explicitly.sequence = S;
    explicitly.labels = { "ACC1" };
    LabelOracle other(*anno, &cth);
    auto one = traverse_seed(other, explicitly, Strategy(), LabelChangeCost::forbid());
    EXPECT_EQ(1u, one.num_seed_labels);
    EXPECT_EQ(1u, one.labels_supporting_total);
}

// Under `support: trace` the derivation keeps the live coordinates of every candidate at
// every seed k-mer, so that the per-label validation needs no second read. Those entries
// are only ever read for the candidates that survive to the end, and the prefix is
// compacted once the live set halves. Fixture: eight records carry the seed's first 100
// k-mers (so the first batch cannot narrow it) and only two carry the rest, which halves
// the live set in the second batch — exactly when the prefix is compacted.
TEST(WalkerDerive, TraceCoordinatesSurviveCompaction) {
    auto b = clean_blocks({ 110, 100, 40, 40, 40, 40, 40, 40 }, 53);
    const std::string U = b[0], V = b[1];
    const std::string seed_seq = U + V;
    ASSERT_GT(U.size() - kK + 1, 64u) << "the narrowing must fall outside the first batch";

    std::vector<std::string> seqs { seed_seq, seed_seq };
    std::vector<std::string> labels { "cA", "cB" };
    for (size_t i = 0; i < 6; ++i) {
        seqs.push_back(U + b[2 + i]);
        labels.push_back("d" + std::to_string(i));
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, seqs, labels, DeBruijnGraph::BASIC, true);
    LabelOracle oracle(*anno);
    ASSERT_TRUE(oracle.has_coordinates());

    Seed seed;
    seed.sequence = seed_seq;
    Strategy st;
    st.support = Support::TRACE;
    st.merge_reconverge = false;         // trace evidence cannot survive a merge
    auto derived = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    EXPECT_TRUE(derived.labels_from_seed);
    ASSERT_EQ(2u, derived.num_seed_labels);
    EXPECT_EQ(2u, derived.labels_supporting_total);
    EXPECT_TRUE(derived.dropped_labels.empty());
    EXPECT_EQ("cA", derived.label_dict[0].name);
    EXPECT_EQ("cB", derived.label_dict[1].name);

    // the coordinates carried through the compaction are the ones the explicit run finds
    LabelOracle other(*anno);
    Seed named = seed;
    named.labels = { "cA", "cB" };
    auto explicitly = traverse_seed(other, named, st, LabelChangeCost::forbid());
    EXPECT_EQ(serialize(explicitly), serialize(derived));
}


// A derived header must resolve back to its own column. Two columns hold a record
// named ACC; only column B's carries the seed. The derivation finds B's, but an
// explicit list resolves "ACC" to the FIRST column holding it, so the derived list
// would not be resubmittable: it is refused, naming the header. The column kind and
// an explicit column name still work.
TEST(WalkerCoord, DerivedHeaderMustResolveBackToItsColumn) {
    auto b = clean_blocks({ 60, 40, 60 }, 61);
    const std::string &L = b[0], &M = b[1], &T = b[2];
    const std::string in_a = L, in_b = M + T;   // the seed M is only in column B's record
    const uint64_t na = in_a.size() - kK + 1, nb = in_b.size() - kK + 1;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { in_a, in_b }, { "A", "B" }, DeBruijnGraph::BASIC, true, { 0, 0 });
    std::vector<std::vector<std::string>> headers { { "ACC" }, { "ACC" } };
    std::vector<std::vector<uint64_t>> num_kmers { { na }, { nb } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);

    Seed seed;
    seed.sequence = M;
    Strategy st;
    try {
        traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        FAIL() << "a derived header that resolves to another column was accepted";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("'ACC'")) << e.what();
        EXPECT_NE(std::string::npos, std::string(e.what()).find("another annotation column")) << e.what();
        EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
        EXPECT_EQ("ACC", e.subject());
    }
    // the column kind derives B, and B explicitly still works
    st.seed_label_kind = LabelKind::COLUMN;
    auto by_column = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_EQ(1u, by_column.label_dict.size());
    EXPECT_EQ("B", by_column.label_dict[0].name);
    EXPECT_EQ(T.size(), by_column.arms[kRight].paths[0].length_bp);
    seed.labels = { "B" };
    auto explicitly = traverse_seed(oracle, seed, Strategy(), LabelChangeCost::forbid());
    EXPECT_EQ(T.size(), explicitly.arms[kRight].paths[0].length_bp);
}


// Round 3, finding 3: a derived header that is ALSO the name of an annotation column.
// Column "ACC" holds an unrelated record (header "decoy"); column B's record carries
// the seed under the header "ACC". find_header("ACC") gives B's record, so the old
// round-trip check passed — but an explicit label list resolves "ACC" to the COLUMN
// first, a label that does not carry the seed. The derivation now round-trips through
// the resolver the explicit path uses and refuses; the column kind still derives B,
// and B explicitly still works.
TEST(WalkerCoord, DerivedHeaderMustNotBeAColumnName) {
    auto b = clean_blocks({ 40, 40, 30 }, 941);
    const std::string &D = b[0], &M = b[1], &T = b[2];
    const std::string in_acc = D, in_b = M + T;   // the seed M is only in column B's record
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { in_acc, in_b }, { "ACC", "B" }, DeBruijnGraph::BASIC, true, { 0, 0 });
    const auto &enc = anno->get_annotator().get_label_encoder();
    std::vector<std::vector<std::string>> headers(2);
    std::vector<std::vector<uint64_t>> num_kmers(2);
    headers[enc.encode("ACC")] = { "decoy" };
    num_kmers[enc.encode("ACC")] = { in_acc.size() - kK + 1 };
    headers[enc.encode("B")] = { "ACC" };
    num_kmers[enc.encode("B")] = { in_b.size() - kK + 1 };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*anno, &cth);
    // what the explicit path makes of the name: the column, and the header is B's
    EXPECT_EQ(LabelKind::COLUMN, oracle.resolve_label("ACC").kind);
    ASSERT_TRUE(oracle.find_header("ACC").has_value());
    EXPECT_EQ(static_cast<uint64_t>(enc.encode("B")),
              static_cast<uint64_t>(oracle.find_header("ACC")->column));

    Seed seed;
    seed.sequence = M;
    Strategy st;
    st.direction = Strategy::RIGHT;
    try {
        traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        FAIL() << "a derived header that an explicit list resolves to a column was accepted";
    } catch (const SeedDerivationError &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("'ACC'")) << e.what();
        EXPECT_NE(std::string::npos, std::string(e.what()).find("name of an annotation column")) << e.what();
        EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
        EXPECT_EQ("ACC", e.subject());
    }
    // the column kind derives B, and B explicitly still works
    st.seed_label_kind = LabelKind::COLUMN;
    auto by_column = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_EQ(1u, by_column.label_dict.size());
    EXPECT_EQ("B", by_column.label_dict[0].name);
    EXPECT_EQ(LabelKind::COLUMN, by_column.label_dict[0].kind);
    EXPECT_EQ(T.size(), by_column.arms[kRight].paths[0].length_bp);
    seed.labels = { "B" };
    auto explicitly = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    EXPECT_EQ(T.size(), explicitly.arms[kRight].paths[0].length_bp);
    // and the name the old check would have echoed does not carry the seed explicitly
    seed.labels = { "ACC" };
    EXPECT_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()), std::invalid_argument);
}


// Whether a seed's labels can be derived must not depend on annotation.batch_kmers. A
// header record A^66000·tail: the k-mer A^11 alone carries ~66k coordinates, over the
// guard for max_seed_labels = 1; the cheapest-row choice and the guard look at a window
// of 64 k-mers whatever the fetch batch is, so a seed starting in the homopolymer is
// accepted with batch 1 as with batch 64, and the two runs are identical.
TEST(WalkerCoord, DerivationDoesNotDependOnBatchKmers) {
    const std::string tail = "CGTCGACTGCTACGTACGATCGATGC";
    const std::string record = std::string(66000, 'A') + tail;
    const uint64_t n = record.size() - kK + 1;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { record }, { "F" }, DeBruijnGraph::BASIC, true, { 0 });
    std::vector<std::vector<std::string>> headers { { "H" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));

    Seed seed;
    seed.sequence = std::string(kK, 'A') + tail.substr(0, 8);
    std::string first;
    for (size_t batch : { 1, 64 }) {
        LabelOracle oracle(*anno, &cth);
        Strategy st;
        st.max_seed_labels = 1;
        st.batch_kmers = batch;
        st.max_extension_bp = 200;
        SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        ASSERT_EQ(1u, r.label_dict.size()) << "batch " << batch;
        EXPECT_EQ("H", r.label_dict[0].name);
        EXPECT_EQ(tail.size() - 8, r.arms[kRight].paths[0].length_bp) << "batch " << batch;
        if (first.empty()) {
            first = serialize(r);
        } else {
            EXPECT_EQ(first, serialize(r)) << "batch " << batch;
        }
    }
}


// A head that has already reached the radius when a seed-level cap trips is complete,
// not out of budget: two linear branches, radius 2, max_steps 3 — the first head
// reaches depth 2 (max_extension_bp), the second is cut at depth 1 (max_steps), and
// the boundary is 1.
TEST(Walker, HeadsAtTheRadiusAreCompleteWhenACapTrips) {
    auto b = clean_blocks({ 30, 40, 40 }, 62);
    for (uint32_t s = 63; b[1][0] == b[2][0]; ++s)
        b = clean_blocks({ 30, 40, 40 }, s);
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
        for (size_t r = 0; r < kNumEndReasons; ++r) {
            if (p.end_reasons[r])
                ends[p.length_bp] = static_cast<EndReason>(r);
        }
    }
    ASSERT_EQ(2u, ends.size());
    EXPECT_STREQ("max_extension_bp", to_string(ends.at(2)));
    EXPECT_STREQ("max_steps", to_string(ends.at(1)));
    // the head at the radius is complete, not remaining: one path and its one label
    // were cut, and that is what the frontier and the cap trigger report
    EXPECT_EQ(1u, arm.frontier_live_paths);
    EXPECT_EQ(1u, arm.frontier_live_labels);
    ASSERT_TRUE(arm.cap_trigger.has_value());
    EXPECT_EQ(1u, arm.cap_trigger->live_paths);
    EXPECT_EQ(1u, arm.cap_trigger->live_labels);
    EXPECT_EQ(1u, arm.cap_trigger->at_bp);
}


namespace {

std::string random_dna(std::mt19937 &rng, size_t len) {
    std::string s;
    for (size_t i = 0; i < len; ++i) s += "ACGT"[rng() % 4];
    return s;
}

// does some root-to-leaf route of |arm| spell |wanted| as a prefix?
bool spells_prefix(const ArmResult &arm, size_t s, const std::string &wanted, size_t at = 0) {
    const Segment &seg = arm.segments[s];
    const std::string x = seg.sequence.substr(0, wanted.size() - at);
    if (wanted.compare(at, x.size(), x) != 0)
        return false;
    at += x.size();
    if (at == wanted.size())
        return true;
    for (size_t c : seg.children) {
        if (spells_prefix(arm, c, wanted, at))
            return true;
    }
    return false;
}

} // namespace

// Annotate mode over a MERGED DAG: two records of one label that disagree in 26 places
// and rejoin after each — O(N) segments but 2^N routes. The summary must be computed per
// segment (union of the parents' surviving sets), not per route; the per-route walk it
// replaces took seconds here and grew exponentially.
TEST(Walker, AnnotateSummaryIsLinearOnAMergedDag) {
    const size_t k = 21, n = 26;
    std::mt19937 rng(18881);
    const std::string seed_seq = random_dna(rng, 35);
    std::string a = seed_seq, b = seed_seq;
    for (size_t i = 0; i < n; ++i) {
        std::string x = random_dna(rng, 35), y = random_dna(rng, 35), join = random_dna(rng, 35);
        while (x[0] == y[0]) y = random_dna(rng, 35);
        a += x + join;
        b += y + join;
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, { a, b }, { "A", "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Seed seed;
    seed.sequence = seed_seq;
    Strategy st;
    st.label_mode = LabelMode::ANNOTATE;
    st.direction = Strategy::RIGHT;
    st.merge_reconverge = true;
    st.max_extension_bp = 10000;
    const auto start = std::chrono::steady_clock::now();
    auto r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const ArmResult &arm = r.arms[kRight];
    EXPECT_EQ(ArmResult::COMPLETE, arm.status);
    EXPECT_LT(arm.segments.size(), 4 * n);   // a compact DAG, not a trie of 2^n leaves
    EXPECT_EQ(1u, arm.paths.size());
    ASSERT_EQ(1u, r.label_dict.size());
    // the label is on every node of both records, so its direct support is the whole flank
    EXPECT_EQ(a.size() - seed_seq.size(), r.label_summary[0][kRight].direct_bp);
    EXPECT_EQ(a.size() - seed_seq.size(), r.label_summary[0][kRight].reach_bp);
    EXPECT_LT(elapsed, 5.0) << "the summary is not linear in the DAG";
}


// The completeness certificate under merging: a merge unites the edge histories of the
// routes it joins, so a walk admissible on its own path (A's P·R·Q·E, 190 bp) is present
// under `keep` and lost under `merge` where B's Q edges are imported into A's history.
// The walk rule the response carries must say so whenever merging is on.
TEST(Walker, MergedWalkRuleIsQualified) {
    const size_t k = 21;
    std::mt19937 rng(8042);
    const std::string S = random_dna(rng, 30), P = random_dna(rng, 70), Y = random_dna(rng, 40),
                      R = random_dna(rng, 40), Q = random_dna(rng, 40), E = random_dna(rng, 40);
    const std::string A = S + P + R + Q + E, B = S + Y.substr(0, 30) + Q + R;
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, { A, B }, { "A", "B" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    for (bool merge : { false, true }) {
        Seed seed;
        seed.sequence = S;
        Strategy st;
        st.label_mode = LabelMode::ANNOTATE;
        st.direction = Strategy::RIGHT;
        st.merge_reconverge = merge;
        st.max_extension_bp = 300;
        auto r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        const ArmResult &arm = r.arms[kRight];
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << merge;
        // the scope of the certificate is a field, not only prose
        EXPECT_STREQ(merge ? "united_history" : "per_path", arm.completeness_scope);
        const std::string rule = walk_rule_statement(st, oracle);
        if (!merge) {
            EXPECT_TRUE(spells_prefix(arm, 0, P + R + Q + E)) << "keep lost A's own walk";
            EXPECT_EQ(std::string::npos, rule.find("edge histories are united")) << rule;
        } else {
            // the conservative merge: pinned as the documented behaviour, not endorsed
            EXPECT_FALSE(spells_prefix(arm, 0, P + R + Q + E));
            EXPECT_NE(std::string::npos, rule.find("edge histories are united")) << rule;
            EXPECT_NE(std::string::npos, rule.find("on_reconverge: keep gives the per-path set")) << rule;
        }
    }
}


// Label-free exploration (spec §6.11): annotate mode with a beam of width 1 and
// most_supported_first follows, at every fork, the branch carried by more labels — here
// P (3 labels) over Q (1) at the seed boundary and then P1 (2) over P2 (1) — and reports
// itself as pruned with the per-node support recorded; breadth_first would keep whichever
// head was created first.
TEST(Walker, LabelFreeBeamFollowsTheMostSupportedBranch) {
    auto b = clean_blocks({ 30, 30, 30, 30, 30 }, 64);
    for (uint32_t s = 65; b[1][0] == b[2][0] || b[3][0] == b[4][0]; ++s)
        b = clean_blocks({ 30, 30, 30, 30, 30 }, s);
    const std::string &X = b[0], &P = b[1], &Q = b[2], &P1 = b[3], &P2 = b[4];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { X + P + P1, X + P + P1, X + P + P2, X + Q }, { "A", "B", "C", "D" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Seed seed;
    seed.sequence = X;
    Strategy st;
    st.label_mode = LabelMode::ANNOTATE;
    st.direction = Strategy::RIGHT;
    st.on_overflow = Strategy::BEAM;
    st.max_live_paths = 1;
    st.order = Strategy::MOST_SUPPORTED_FIRST;
    st.max_extension_bp = 200;
    auto r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    const ArmResult &arm = r.arms[kRight];
    EXPECT_EQ(ArmResult::PRUNED, arm.status);
    ASSERT_TRUE(arm.cap_trigger.has_value());
    EXPECT_STREQ("beam_pruned", to_string(arm.cap_trigger->reason));
    // the pruned heads stay on record as stubs ended beam_pruned; exactly one walk goes on
    const PathResult *kept = nullptr;
    for (const PathResult &p : arm.paths) {
        ASSERT_TRUE(p.path_reason.has_value());
        if (*p.path_reason == EndReason::BEAM_PRUNED) {
            EXPECT_LE(p.length_bp, P.size() + 1);
        } else {
            EXPECT_EQ(nullptr, kept) << "two walks survived a beam of width 1";
            kept = &p;
        }
    }
    ASSERT_NE(nullptr, kept);
    EXPECT_EQ(P + P1, spell_path(arm, *kept));
    EXPECT_STREQ("dead_end", to_string(*kept->path_reason));
    // pruned right after the first level: the 1-base walks are all present, nothing longer
    EXPECT_EQ(1u, arm.complete_to_bp);
    // the per-node support is on record: 3 labels along P, 2 along P1
    auto labels_at_depth = [&](uint64_t d) -> size_t {
        for (size_t s : path_segments(arm, *kept)) {
            for (const LabelSetRun &run : arm.segments[s].label_sets) {
                if (run.from_bp < d && d <= run.to_bp) {
                    EXPECT_FALSE(run.truncated());
                    return run.labels.size();
                }
            }
        }
        ADD_FAILURE() << "no recorded set at depth " << d;
        return 0;
    };
    EXPECT_EQ(3u, labels_at_depth(P.size()));
    EXPECT_EQ(2u, labels_at_depth(P.size() + P1.size()));
    // and per sample, how far some recorded route carries it: A and B the whole walk, C
    // to the P1 fork plus the one-base P2 stub the beam left there, D only the one-base Q
    // stub — the stubs are recorded nodes and count, the spelled walk is not what
    // direct_bp measures
    std::map<std::string, uint64_t> direct;
    for (size_t i = 0; i < r.label_dict.size(); ++i)
        direct[r.label_dict[i].name] = r.label_summary[i][kRight].direct_bp;
    EXPECT_EQ(P.size() + P1.size(), direct["A"]);
    EXPECT_EQ(P.size() + P1.size(), direct["B"]);
    EXPECT_EQ(P.size() + 1, direct["C"]);
    EXPECT_EQ(1u, direct["D"]);
}


// Round 3, finding 4: the beam ranked heads by the RECORDED label list, which
// labels.max_labels_per_node cuts — with a cap of 1 a branch carried by one label and
// one carried by three tied, and the tie went to the head created first. The minority
// branch here is the one enumerated first (its first base sorts lower); the TRUE count
// ranks the three-label branch above it, under a cap of 1 as under none, and under
// either frontier order, since the beam keeps the most supported heads whatever order
// expands a level (§6.11).
TEST(Walker, LabelFreeBeamRanksByTheTrueCount) {
    auto b = clean_blocks({ 30, 40, 40 }, 664);
    for (uint32_t s = 665; b[1][0] >= b[2][0]; ++s)
        b = clean_blocks({ 30, 40, 40 }, s);
    const std::string &X = b[0], &minor = b[1], &major = b[2];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        kK, { X + minor, X + major, X + major, X + major }, { "A", "B", "C", "D" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    for (auto order : { Strategy::BREADTH_FIRST, Strategy::MOST_SUPPORTED_FIRST }) {
        for (size_t cap : { size_t(1), size_t(64) }) {
            const std::string what = "order " + std::to_string(order) + " cap " + std::to_string(cap);
            Seed seed;
            seed.sequence = X;
            Strategy st;
            st.label_mode = LabelMode::ANNOTATE;
            st.direction = Strategy::RIGHT;
            st.merge_reconverge = false;
            st.on_overflow = Strategy::BEAM;
            st.max_live_paths = 1;
            st.order = order;
            st.max_extension_bp = 100;
            st.max_labels_per_node = cap;
            auto r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
            const ArmResult &arm = r.arms[kRight];
            EXPECT_EQ(ArmResult::PRUNED, arm.status) << what;
            // the cut is reported; it must not have decided the ranking
            EXPECT_EQ(cap == 1, arm.nodes_labels_truncated > 0) << what;
            const PathResult *kept = nullptr;
            for (const PathResult &p : arm.paths) {
                ASSERT_TRUE(p.path_reason.has_value()) << what;
                if (*p.path_reason == EndReason::BEAM_PRUNED)
                    continue;
                EXPECT_EQ(nullptr, kept) << what << ": two walks survived a beam of width 1";
                kept = &p;
            }
            ASSERT_NE(nullptr, kept) << what;
            EXPECT_EQ(major, spell_path(arm, *kept)) << what << ": the beam kept the minority branch";
        }
    }
}


// Review round 4, finding 2: recording the refusals of an ambiguous node is linear in its
// labels. k = 3, n labels each carrying AAA·C and AAA·G, seed AAA, branch allowance 0:
// every label is ambiguous and excluded in the first round. Each successor used to be
// scanned once PER excluded source for that source's entries — n(n + 1) tests for the 2n
// refused label ids (10 000 labels: 7.5 ms before refusals were recorded, 104 ms after).
// One scan of each successor's state per round now: 2n entries, so the work doubles with
// n. Pinned through ArmResult::refusal_scans, not wall time; the refusals are unchanged.
TEST(Walker, RefusalRecordingIsLinearInTheLabels) {
    for (size_t n : { 250, 500 }) {
        std::vector<std::string> seqs, labels, names;
        for (size_t i = 0; i < n; ++i) {
            names.push_back("L" + std::to_string(i));
            for (const char *s : { "AAAC", "AAAG" }) {
                seqs.push_back(s);
                labels.push_back(names.back());
            }
        }
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                3, seqs, labels, DeBruijnGraph::BASIC);
        Strategy st = strategy(0, false);
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 10;
        auto res = run(*anno, "AAA", names, st);
        const ArmResult &arm = res.arms[kRight];
        check_invariants(arm, st);
        EXPECT_EQ(n, count_ends(arm, EndReason::BRANCH));
        ASSERT_EQ(1u, arm.branch_events.size());
        const BranchEvent &be = arm.branch_events[0];
        EXPECT_EQ(n, be.ambiguous.size());
        EXPECT_EQ(n, be.dropped.size());
        std::vector<LabelId> all(n);
        std::iota(all.begin(), all.end(), 0);
        ASSERT_EQ(2u, be.refused.size());
        std::set<char> chars;
        for (const auto &rf : be.refused) {
            chars.insert(rf.ch);
            EXPECT_STREQ("branch", rf.cause);
            EXPECT_EQ(all, rf.labels);
        }
        EXPECT_EQ((std::set<char>{ 'C', 'G' }), chars);
        // two successors of n entries, each scanned once: linear, not n(n + 1)
        EXPECT_EQ(2 * n, arm.refusal_scans) << n;
        EXPECT_LE(arm.refusal_scans, 4 * n) << n;
    }
}

// Review round 4, finding 3: a successor the loss budget refuses to a source is stated
// whether or not the source goes on along another one. Blocks X, U, P, V, Q: A carries
// X·U·P, B carries X·V and U's last k - 1 bases + Q, the seed is X under {A, B}. B leaves
// along V; at the end of U, A goes on along P and could enter Q only by switching to B
// (constant cost 1). With budget 0 Q is refused to A by the budget alone — recorded only
// for a source that ENDED at the node until now (has_cont skipped A), so the fork had no
// branch event and Q's absence no reason. With budget 1, Q is followed through the switch.
TEST(Walker, LossBudgetRefusalWhileTheSourceGoesOn) {
    std::vector<std::string> b;
    for (uint32_t seed = 582; ; ++seed) {
        b = clean_blocks({ 30, 30, 40, 40, 40 }, seed);
        if (b[1][0] != b[3][0] && b[2][0] != b[4][0])
            break;
    }
    const std::string &X = b[0], &U = b[1], &P = b[2], &V = b[3], &Q = b[4];
    for (auto mode : all_modes()) {
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, { X + U + P, X + V, U.substr(U.size() - kK + 1) + Q }, { "A", "B", "B" }, mode);
        const std::string where = "mode " + std::to_string(mode);
        Strategy st = strategy(Strategy::kUnlimited, false);
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 100;
        st.loss_budget = 0;
        auto r0 = run(*anno, X, { "A", "B" }, st, LabelChangeCost::constant(1));
        const ArmResult &arm = r0.arms[kRight];
        check_invariants(arm, st);
        std::set<std::string> walks;
        for (const auto &p : arm.paths) walks.insert(spell_path(arm, p));
        EXPECT_EQ((std::set<std::string>{ U + P, V }), walks) << where;
        EXPECT_EQ(0u, count_events(arm, EventType::SWITCH)) << where;
        const BranchEvent *fork = nullptr;
        for (const BranchEvent &be : arm.branch_events) {
            EXPECT_EQ(U.size(), be.at_bp) << where << ": a branch event off the fork";
            if (be.at_bp == U.size())
                fork = &be;
        }
        ASSERT_NE(nullptr, fork) << where << ": the fork has no branch event";
        ASSERT_EQ(1u, fork->refused.size()) << where;
        EXPECT_EQ(Q[0], fork->refused[0].ch) << where;
        EXPECT_STREQ("loss_budget", fork->refused[0].cause) << where;
        EXPECT_EQ((std::vector<LabelId>{ 0 }), fork->refused[0].labels) << where;
        // A goes on: nothing ambiguous, no label end at the fork, only the refusal
        EXPECT_TRUE(fork->ambiguous.empty()) << where;
        EXPECT_TRUE(fork->dropped.empty()) << where;
        EXPECT_EQ(0u, count_ends(arm, EndReason::LOSS_BUDGET)) << where;

        st.loss_budget = 1;
        auto r1 = run(*anno, X, { "A", "B" }, st, LabelChangeCost::constant(1));
        const ArmResult &arm1 = r1.arms[kRight];
        check_invariants(arm1, st);
        walks.clear();
        for (const auto &p : arm1.paths) walks.insert(spell_path(arm1, p));
        EXPECT_EQ((std::set<std::string>{ U + P, U + Q, V }), walks) << where;
        EXPECT_EQ(1u, count_events(arm1, EventType::SWITCH)) << where;
        for (const BranchEvent &be : arm1.branch_events) {
            EXPECT_TRUE(be.refused.empty()) << where << ": a refusal within the budget";
        }
    }
}


Json::Value parse_json(const std::string &text) {
    Json::Value v;
    Json::CharReaderBuilder builder;
    std::string errs;
    std::istringstream in(text);
    EXPECT_TRUE(Json::parseFromStream(builder, in, &v, &errs)) << errs;
    return v;
}

std::vector<std::string> kinds_of(const Json::Value &limitations) {
    EXPECT_TRUE(limitations.isArray());
    std::vector<std::string> out;
    for (const Json::Value &l : limitations) {
        out.push_back(l["kind"].asString());
    }
    return out;
}

// Branch events are kept in level order up to output.max_branch_events, and the cut is
// a depth boundary (spec §7.2): one label on a main walk X·Y with side branches off Y at
// three depths is ambiguous at each fork, so the right arm produces one event per
// depth. Keeping n of them puts the boundary at the depth of the (n+1)-th; every event
// kept lies below it, and the response says so — `evidence` and a `branch_events`
// limitation naming the knob — in every detail level. Keeping exactly as many as were
// produced cuts nothing.
TEST(Walker, BranchEventsCompleteToBp) {
    const std::vector<uint64_t> depths { 10, 25, 40 };
    std::vector<std::string> b;
    for (uint32_t seed = 50; ; ++seed) {
        b = clean_blocks({ 30, 60, 20, 20, 20 }, seed);
        bool forks = true;
        for (size_t i = 0; i < depths.size(); ++i) {
            forks &= b[2 + i][0] != b[1][depths[i]];
        }
        if (forks)
            break;
    }
    const std::string &X = b[0], &Y = b[1];
    std::vector<std::string> seqs { X + Y };
    for (size_t i = 0; i < depths.size(); ++i) {
        seqs.push_back(X + Y.substr(0, depths[i]) + b[2 + i]);
    }
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, seqs, std::vector<std::string>(seqs.size(), "A"), DeBruijnGraph::BASIC);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.max_label_branches = Strategy::kUnlimited;
    st.merge_reconverge = false;
    for (size_t cap : { size_t(0), size_t(1), size_t(2), size_t(3), Strategy::kUnlimited }) {
        st.max_branch_events = cap;
        const std::string where = "max_branch_events "
            + (cap == Strategy::kUnlimited ? std::string("unlimited") : std::to_string(cap));
        const SeedResult res = run(*anno, X, { "A" }, st);
        const ArmResult &arm = res.arms[kRight];
        check_invariants(arm, st);
        // the walk does not depend on the cap: the main walk and one leaf per side branch
        EXPECT_EQ(depths.size() + 1, arm.paths.size()) << where;
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << where;
        ASSERT_EQ(depths.size(), arm.branch_events_total) << where;
        const size_t kept = std::min(cap, depths.size());
        ASSERT_EQ(kept, arm.branch_events.size()) << where;
        for (size_t i = 0; i < kept; ++i) {
            EXPECT_EQ(depths[i], arm.branch_events[i].at_bp) << where;
            EXPECT_EQ((std::vector<LabelId>{ 0 }), arm.branch_events[i].ambiguous) << where;
        }
        // the boundary is the depth of the first event not kept, and every event kept
        // lies below it
        const bool complete = kept == depths.size();
        const uint64_t boundary = complete ? std::numeric_limits<uint64_t>::max() : depths[kept];
        EXPECT_EQ(boundary, arm.branch_events_complete_to_bp) << where;
        for (const BranchEvent &be : arm.branch_events) {
            EXPECT_LT(be.at_bp, arm.branch_events_complete_to_bp) << where;
        }

        for (const char *detail : { "summary", "tree", "full" }) {
            const Json::Value seed_json = cli::seed_result_to_json(res, st, detail, false);
            const Json::Value &j = seed_json["arms"]["right"];
            const std::string at = where + ", detail " + detail;
            // the walks are complete either way; only the branch diagnostics are cut
            EXPECT_EQ("complete", j["status"].asString()) << at;
            EXPECT_EQ(complete ? "complete/complete/complete/inline" : "complete/cut/complete/inline",
                      outcome_of(seed_json)) << at;
            ASSERT_TRUE(j.isMember("evidence")) << at;
            EXPECT_EQ(complete, j["evidence"]["complete"].asBool()) << at;
            if (complete) {
                EXPECT_TRUE(j["evidence"]["complete_to_bp"].isNull()) << at;
                EXPECT_EQ(std::vector<std::string>{}, kinds_of(j["limitations"])) << at;
            } else {
                EXPECT_EQ(boundary, j["evidence"]["complete_to_bp"].asUInt64()) << at;
                ASSERT_EQ(std::vector<std::string>{ "branch_events" }, kinds_of(j["limitations"])) << at;
                const Json::Value &l = j["limitations"][0];
                EXPECT_EQ("output.max_branch_events", l["knob"].asString()) << at;
                EXPECT_EQ(cap, l["limit"].asUInt64()) << at;
                // more were produced than the knob kept
                EXPECT_EQ(depths.size(), l["observed"].asUInt64()) << at;
                EXPECT_EQ(boundary, l["complete_to_bp"].asUInt64()) << at;
                EXPECT_NE(std::string::npos, l["effect"].asString().find("unlimited")) << at;
            }
            EXPECT_EQ(depths.size() - kept, j["branch_events_truncated"].asUInt64()) << at;
        }
    }

    // the knob takes "unlimited" and echoes it so; the default stays 100
    Json::Value request = parse_json(R"({"seeds": [{"sequence": "ACGTACGTACGTACGT"}],
                                         "strategy": {"output": {"max_branch_events": "unlimited"}}})");
    cli::TraverseRequest parsed = cli::parse_traverse_request(request);
    EXPECT_EQ(Strategy::kUnlimited, parsed.strategy.max_branch_events);
    EXPECT_EQ("unlimited",
              cli::strategy_to_json(parsed.strategy, parsed.cost)["output"]["max_branch_events"].asString());
    request["strategy"]["output"]["max_branch_events"] = 7;
    EXPECT_EQ(7u, cli::parse_traverse_request(request).strategy.max_branch_events);
    request["strategy"]["output"]["max_branch_events"] = "lots";
    EXPECT_THROW(cli::parse_traverse_request(request), cli::InvalidRequest);
    request["strategy"].removeMember("output");
    parsed = cli::parse_traverse_request(request);
    EXPECT_EQ(100u, parsed.strategy.max_branch_events);
    EXPECT_EQ(100u, cli::strategy_to_json(parsed.strategy, parsed.cost)["output"]["max_branch_events"].asUInt64());
}

// Every cap that limited a result is stated in `limitations` (spec §7.0) — and only
// those: each kind is produced here by the one knob that causes it and is absent from
// the unconstrained run. A bubble on the right of the seed X (A and B on P, C on Q, both
// closing into Y) gives a fork, a reconvergence and three labels at every node.
TEST(Walker, LimitationsStateExactlyWhatLimitedTheResult) {
    const auto b = fork_blocks(60, { 30, 20, 20, 30 });
    const std::string &X = b[0], &P = b[1], &Q = b[2], &Y = b[3];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P + Y, X + P + Y, X + Q + Y }, { "A", "B", "C" }, DeBruijnGraph::BASIC);
    auto request = [&](const std::vector<std::string> &labels, const std::string &strategy) {
        Json::Value r;
        Json::Value seed;
        seed["sequence"] = X;
        for (const std::string &l : labels) {
            seed["labels"].append(l);
        }
        r["seeds"].append(seed);
        r["strategy"] = parse_json(strategy);
        return r;
    };
    auto traverse = [&](const Json::Value &req, const cli::TraverseLimits &limits) {
        return cli::process_traverse_request(req, *anno, "", limits);
    };
    const cli::TraverseLimits none;
    const std::vector<std::string> nothing;
    const std::vector<std::string> ABC { "A", "B", "C" };

    // nothing limits a complete per-path run with every label and every event
    Json::Value out = traverse(request(ABC, R"({"direction": "right",
                                                "branching": {"on_reconverge": "keep"}})"), none);
    Json::Value result = out["results"][0];
    Json::Value right = result["arms"]["right"];
    EXPECT_EQ("complete", right["status"].asString());
    EXPECT_EQ(2u, right["paths"].size());
    EXPECT_EQ(nothing, kinds_of(right["limitations"]));
    EXPECT_EQ(nothing, kinds_of(result["limitations"]));
    EXPECT_TRUE(right["evidence"]["complete"].asBool());
    EXPECT_TRUE(right["evidence"]["complete_to_bp"].isNull());
    EXPECT_EQ("complete/complete/complete/inline", outcome_of(result));

    // walk_domain: the fork needs two live paths
    out = traverse(request(ABC, R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                    "bounds": {"max_live_paths": 1}})"), none);
    right = out["results"][0]["arms"]["right"];
    EXPECT_EQ("partial/complete/complete/inline", outcome_of(out["results"][0]));
    EXPECT_EQ("truncated", right["status"].asString());
    ASSERT_EQ(std::vector<std::string>{ "walk_domain" }, kinds_of(right["limitations"]));
    const Json::Value &walk = right["limitations"][0];
    EXPECT_EQ("bounds.max_live_paths", walk["knob"].asString());
    EXPECT_EQ(1u, walk["limit"].asUInt64());
    EXPECT_EQ(2u, walk["observed"].asUInt64());   // what the fork needed
    EXPECT_EQ(right["complete_to_bp"], walk["complete_to_bp"]);
    EXPECT_EQ(0u, walk["complete_to_bp"].asUInt64());
    EXPECT_EQ(nothing, kinds_of(out["results"][0]["limitations"]));

    // ... and a beam of width one prunes at the fork: pruned, so partial as well
    out = traverse(request(ABC, R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                    "frontier": {"on_overflow": "beam"},
                                    "bounds": {"max_live_paths": 1}})"), none);
    right = out["results"][0]["arms"]["right"];
    EXPECT_EQ("pruned", right["status"].asString());
    EXPECT_EQ("partial/complete/complete/inline", outcome_of(out["results"][0]));
    ASSERT_EQ(std::vector<std::string>{ "walk_domain" }, kinds_of(right["limitations"]));
    EXPECT_EQ("bounds.max_live_paths", right["limitations"][0]["knob"].asString());

    // scope: merging (constrain's default) closes the bubble once. Every walk of the
    // united-history rule is present, but not every walk per path: a merge that united a
    // history makes the walks partial (conservative outcome, DESIGN-traverse-graphlet.md
    // §14 v5.2), stated by the limitation and by completeness_scope
    out = traverse(request(ABC, R"({"direction": "right"})"), none);
    EXPECT_EQ("partial/complete/complete/inline", outcome_of(out["results"][0]));
    right = out["results"][0]["arms"]["right"];
    EXPECT_EQ("united_history", right["completeness_scope"].asString());
    ASSERT_EQ(std::vector<std::string>{ "scope" }, kinds_of(right["limitations"]));
    EXPECT_EQ("branching.on_reconverge", right["limitations"][0]["knob"].asString());
    EXPECT_EQ("merge", right["limitations"][0]["limit"].asString());
    EXPECT_EQ(1u, right["limitations"][0]["observed"].asUInt64());

    // label_lists and inexact_counts: annotate mode, one label per recorded list ...
    out = traverse(request(nothing, R"({"direction": "right",
                                        "labels": {"mode": "annotate", "max_labels_per_node": 1}})"), none);
    right = out["results"][0]["arms"]["right"];
    EXPECT_EQ("complete", right["status"].asString());
    // every walk is there, its recorded lists are not
    EXPECT_EQ("complete/complete/lower_bound/inline", outcome_of(out["results"][0]));
    ASSERT_EQ((std::vector<std::string>{ "label_lists", "inexact_counts" }), kinds_of(right["limitations"]));
    for (const Json::Value &l : right["limitations"]) {
        EXPECT_EQ("labels.max_labels_per_node", l["knob"].asString());
        EXPECT_EQ(1u, l["limit"].asUInt64());
    }
    EXPECT_EQ(3u, right["limitations"][0]["observed"].asUInt64());   // the largest true count
    // ... and nothing with the default cap
    out = traverse(request(nothing, R"({"direction": "right", "labels": {"mode": "annotate"}})"), none);
    EXPECT_EQ(nothing, kinds_of(out["results"][0]["arms"]["right"]["limitations"]));

    // seed_labels: the derived carrier set {A, B, C} cut to two
    const std::string derived_two = R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                        "labels": {"max_seed_labels": 2}})";
    out = traverse(request(nothing, derived_two), none);
    result = out["results"][0];
    ASSERT_EQ(std::vector<std::string>{ "seed_labels" }, kinds_of(result["limitations"]));
    EXPECT_EQ("labels.max_seed_labels", result["limitations"][0]["knob"].asString());
    EXPECT_EQ(2u, result["limitations"][0]["limit"].asUInt64());
    EXPECT_EQ(3u, result["limitations"][0]["observed"].asUInt64());
    EXPECT_FALSE(result["limitations"][0].isMember("server_limit"));
    EXPECT_EQ(nothing, kinds_of(result["arms"]["right"]["limitations"]));
    // a carrier the cap left out takes its walks and its evidence with it: both axes say
    // so (conservative outcome, DESIGN-traverse-graphlet.md §14 v5.2), and the seed_labels
    // entry says which knob brings it back
    EXPECT_EQ("partial/complete/lower_bound/inline", outcome_of(result));

    // server_clamp: the server lowers the derived-set cap to 2. The derived seed runs
    // into it (and its seed_labels entry names the server's maximum); a seed with an
    // explicit list in the same request does not.
    Json::Value req = request(nothing, R"({"direction": "right", "branching": {"on_reconverge": "keep"}})");
    Json::Value named;
    named["sequence"] = X;
    for (const std::string &l : ABC) {
        named["labels"].append(l);
    }
    req["seeds"].append(named);
    cli::TraverseLimits limits;
    limits.max_seed_labels = 2;
    out = traverse(req, limits);
    ASSERT_EQ(1u, out["strategy"]["clamped"].size());
    // an integer knob is echoed as an integer (1000, not 1000.0), in the clamp and in the
    // server_clamp limitation built from it
    EXPECT_EQ(Json::uintValue, out["strategy"]["clamped"][0]["requested"].type());
    EXPECT_EQ(Json::uintValue, out["strategy"]["clamped"][0]["effective"].type());
    EXPECT_EQ("1000", Json::writeString(Json::StreamWriterBuilder(),
                                        out["strategy"]["clamped"][0]["requested"]));
    result = out["results"][0];
    ASSERT_EQ((std::vector<std::string>{ "seed_labels", "server_clamp" }), kinds_of(result["limitations"]));
    EXPECT_EQ(2u, result["limitations"][0]["server_limit"].asUInt64());
    const Json::Value &clamp = result["limitations"][1];
    EXPECT_EQ("labels.max_seed_labels", clamp["knob"].asString());
    EXPECT_EQ(2u, clamp["limit"].asUInt64());
    EXPECT_EQ(1000u, clamp["observed"].asUInt64());   // what the request asked for
    EXPECT_EQ(Json::uintValue, clamp["limit"].type());
    EXPECT_EQ(Json::uintValue, clamp["observed"].type());
    EXPECT_EQ(nothing, kinds_of(out["results"][1]["limitations"]));

    // a time budget raised from zero (it also bounds the derivation) binds every seed's
    // walk; a lowered one that never tripped binds none
    limits = cli::TraverseLimits();
    limits.max_time_ms = 60000;
    out = traverse(request(nothing, R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                        "bounds": {"time_budget_ms": 0}})"), limits);
    result = out["results"][0];
    ASSERT_EQ(std::vector<std::string>{ "server_clamp" }, kinds_of(result["limitations"]));
    EXPECT_EQ("bounds.time_budget_ms", result["limitations"][0]["knob"].asString());
    EXPECT_EQ(60000.0, result["limitations"][0]["limit"].asDouble());
    EXPECT_EQ(0.0, result["limitations"][0]["observed"].asDouble());
    out = traverse(request(nothing, R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                        "bounds": {"time_budget_ms": 120000}})"), limits);
    EXPECT_EQ(1u, out["strategy"]["clamped"].size());
    EXPECT_EQ(nothing, kinds_of(out["results"][0]["limitations"]));

    // switch_sources (the SwitchSources fixture): with one priced source B's cheap
    // switch into E is never considered; with two it is
    auto s = clean_blocks({ 60, 40 }, 34);
    const std::string &S = s[0], &Z = s[1];
    auto switching = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { S, S, S, S.substr(S.size() - kK) + Z }, { "A", "B", "C", "E" }, DeBruijnGraph::BASIC);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.extra = { "E" };
    st.loss_budget = 1;
    st.merge_reconverge = false;
    auto table = LabelChangeCost::table({ { { 0, 3 }, 2.0 }, { { 1, 3 }, 0.5 } }, kInfiniteLoss);
    st.max_switch_sources = 1;
    Json::Value seed_json = cli::seed_result_to_json(run(*switching, S, { "A", "B", "C" }, st, table),
                                                     st, "summary", false);
    EXPECT_EQ("complete/complete/lower_bound/inline", outcome_of(seed_json));
    Json::Value arm = seed_json["arms"]["right"];
    ASSERT_EQ(std::vector<std::string>{ "switch_sources" }, kinds_of(arm["limitations"]));
    EXPECT_EQ("labels.max_switch_sources", arm["limitations"][0]["knob"].asString());
    EXPECT_EQ(1u, arm["limitations"][0]["limit"].asUInt64());
    // one derivation was cut while B could switch into E, and it ended B
    EXPECT_EQ(1u, arm["limitations"][0]["observed"].asUInt64());
    EXPECT_EQ(1u, arm["limitations"][0]["label_ends"].asUInt64());
    st.max_switch_sources = 2;
    seed_json = cli::seed_result_to_json(run(*switching, S, { "A", "B", "C" }, st, table),
                                         st, "summary", false);
    EXPECT_EQ("complete/complete/complete/inline", outcome_of(seed_json));
    EXPECT_EQ(nothing, kinds_of(seed_json["arms"]["right"]["limitations"]));
}

// Every seed result carries an `outcome` (spec §7.0). A seed whose permitted set cannot
// be derived has `walks: failed` and states its cause as a `derivation` limitation naming
// the request field that would get past it, while the other seeds of the batch are
// traversed: `walks: complete` when every requested arm reached the radius, `partial`
// when one was truncated or pruned. Records A = Y1·X·P, B = Y2·X·Q (a fork past X), C = W·Z (linear); the seed
// Y1·X·Q is in the graph but no record carries it.
TEST(Walker, FailedDerivationIsAStatedOutcome) {
    std::vector<std::string> b;
    for (uint32_t seed = 71; ; ++seed) {
        b = clean_blocks({ 30, 30, 30, 30, 30, 30, 30 }, seed);
        if (b[3][0] != b[4][0])
            break;
    }
    const std::string &Y1 = b[0], &Y2 = b[1], &X = b[2], &P = b[3], &Q = b[4], &W = b[5], &Z = b[6];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { Y1 + X + P, Y2 + X + Q, W + Z }, { "A", "B", "C" }, DeBruijnGraph::BASIC);
    auto request = [&](const std::vector<std::string> &seeds, const std::string &strategy) {
        Json::Value r;
        for (const std::string &s : seeds) {
            Json::Value seed;
            seed["sequence"] = s;
            r["seeds"].append(seed);
        }
        r["strategy"] = parse_json(strategy);
        return r;
    };
    const std::string orphan = Y1 + X + Q;

    // the walk of W is linear (complete), the fork past X needs two live paths
    // (truncated: partial), and the orphan fails alone
    Json::Value out = cli::process_traverse_request(
            request({ W, orphan, X }, R"({"direction": "right", "branching": {"on_reconverge": "keep"},
                                          "bounds": {"max_extension_bp": 100, "max_live_paths": 1}})"),
            *anno, "", cli::TraverseLimits());
    ASSERT_EQ(3u, out["results"].size());
    const Json::Value &linear = out["results"][0], &failed = out["results"][1], &fork = out["results"][2];
    EXPECT_EQ("complete/complete/complete/inline", outcome_of(linear));
    EXPECT_EQ("complete", linear["arms"]["right"]["status"].asString());
    EXPECT_EQ("partial/complete/complete/inline", outcome_of(fork));
    EXPECT_EQ("truncated", fork["arms"]["right"]["status"].asString());
    for (const Json::Value *ok : { &linear, &fork }) {
        EXPECT_FALSE(ok->isMember("error"));
        EXPECT_TRUE((*ok)["arms"].isMember("right"));
    }

    EXPECT_EQ("failed/complete/complete/inline", outcome_of(failed));
    EXPECT_TRUE(failed.isMember("error"));
    EXPECT_FALSE(failed.isMember("arms"));
    EXPECT_TRUE(failed["seed"]["labels_from_seed"].asBool());
    EXPECT_EQ(orphan.size(), failed["seed"]["length_bp"].asUInt64());
    ASSERT_EQ(std::vector<std::string>{ "derivation" }, kinds_of(failed["limitations"]));
    const Json::Value &d = failed["limitations"][0];
    EXPECT_EQ("no_carrier", d["cause"].asString());
    EXPECT_EQ("seeds[].sequence", d["knob"].asString());
    EXPECT_EQ(orphan.size() - kK + 1, d["limit"].asUInt64());   // the seed's k-mers
    EXPECT_GE(d["observed"].asUInt64(), 1u);                    // read until none was left
    EXPECT_LE(d["observed"].asUInt64(), d["limit"].asUInt64());
    EXPECT_FALSE(d["effect"].asString().empty());
    EXPECT_FALSE(d.isMember("server_limit"));

    // `exhaustive` refuses to cut the derived set {A, B} of X: the cap is the knob, and
    // the server's clamp of it is its server_limit (raising it further does nothing)
    cli::TraverseLimits limits;
    limits.max_seed_labels = 1;
    out = cli::process_traverse_request(
            request({ X }, R"({"exhaustive": true, "direction": "right",
                               "labels": {"max_seed_labels": 5}})"), *anno, "", limits);
    const Json::Value &refused = out["results"][0];
    EXPECT_EQ("failed", refused["outcome"]["walks"].asString());
    ASSERT_EQ(std::vector<std::string>{ "derivation" }, kinds_of(refused["limitations"]));
    EXPECT_EQ("over_seed_label_cap", refused["limitations"][0]["cause"].asString());
    EXPECT_EQ("labels.max_seed_labels", refused["limitations"][0]["knob"].asString());
    EXPECT_EQ(1u, refused["limitations"][0]["limit"].asUInt64());
    EXPECT_EQ(2u, refused["limitations"][0]["observed"].asUInt64());
    EXPECT_EQ(1u, refused["limitations"][0]["server_limit"].asUInt64());

    // the time budget runs out during the derivation; the server lowered it, so the
    // entry carries the server's value as well
    limits = cli::TraverseLimits();
    limits.max_time_ms = 1e-9;
    out = cli::process_traverse_request(
            request({ X }, R"({"direction": "right", "bounds": {"time_budget_ms": 60000}})"),
            *anno, "", limits);
    const Json::Value &late = out["results"][0];
    EXPECT_EQ("failed", late["outcome"]["walks"].asString());
    ASSERT_EQ(std::vector<std::string>{ "derivation" }, kinds_of(late["limitations"]));
    EXPECT_EQ("time_budget", late["limitations"][0]["cause"].asString());
    EXPECT_EQ("bounds.time_budget_ms", late["limitations"][0]["knob"].asString());
    EXPECT_EQ(1e-9, late["limitations"][0]["limit"].asDouble());
    EXPECT_GT(late["limitations"][0]["observed"].asDouble(), 1e-9);
    EXPECT_EQ(1e-9, late["limitations"][0]["server_limit"].asDouble());
}

// Knobs that would be accepted and then do nothing are refused (spec §7.0): the tip and
// bubble windows are not implemented, so a non-zero window is a 400 naming the field (0,
// the echoed value, stays accepted). A continuation is resubmittable only when it is at
// least k long, so output.continuation_bp 1 .. k - 1 is refused naming k; 0 means "no
// continuation sequence" and still reports the labels and the loss.
TEST(Walker, UnimplementedWindowsAndShortContinuationsAreRefused) {
    for (const char *window : { "tip_window_bp", "bubble_window_bp" }) {
        Json::Value request = parse_json(R"({"seeds": [{"sequence": "ACGTACGTACGTACGT"}], "strategy": {}})");
        request["strategy"]["branching"][window] = 20;
        try {
            cli::parse_traverse_request(request);
            FAIL() << "accepted an unimplemented " << window;
        } catch (const cli::InvalidRequest &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find(std::string("strategy.branching.") + window))
                    << e.what();
            EXPECT_NE(std::string::npos, std::string(e.what()).find("not implemented")) << e.what();
        }
        request["strategy"]["branching"][window] = 0;
        cli::TraverseRequest parsed = cli::parse_traverse_request(request);
        EXPECT_EQ(0u, cli::strategy_to_json(parsed.strategy, parsed.cost)["branching"][window].asUInt64());
        // the C++ API refuses it as well
        Strategy st;
        (std::string(window) == "tip_window_bp" ? st.tip_window_bp : st.bubble_window_bp) = 20;
        EXPECT_THROW(validate_strategy(st, LabelChangeCost::forbid()), std::invalid_argument);
    }

    auto b = clean_blocks({ 30, 40 }, 19);
    const std::string &X = b[0], &P = b[1];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P }, { "A" }, DeBruijnGraph::BASIC);
    auto traverse = [&](const std::string &seed, uint64_t continuation_bp) {
        Json::Value request;
        request["seeds"][0]["sequence"] = seed;
        request["seeds"][0]["labels"].append("A");
        request["strategy"] = parse_json(R"({"direction": "right", "bounds": {"max_extension_bp": 10}})");
        request["strategy"]["output"]["continuation_bp"] = static_cast<Json::UInt64>(continuation_bp);
        return cli::process_traverse_request(request, *anno, "", cli::TraverseLimits());
    };
    for (uint64_t short_bp : { uint64_t(1), uint64_t(5), uint64_t(kK - 1) }) {
        try {
            traverse(X, short_bp);
            FAIL() << "accepted continuation_bp " << short_bp << " below k";
        } catch (const cli::InvalidRequest &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find("continuation_bp")) << e.what();
            EXPECT_NE(std::string::npos, std::string(e.what()).find("k = " + std::to_string(kK))) << e.what();
        }
        Strategy st;
        st.continuation_bp = short_bp;
        EXPECT_THROW(run(*anno, X, { "A" }, st), std::invalid_argument);
    }
    // k: the continuation is a valid seed, and resubmitting it works
    Json::Value out = traverse(X, kK);
    const Json::Value &cont = out["results"][0]["arms"]["right"]["paths"][0]["continuation"];
    ASSERT_EQ(kK, cont["sequence"].asString().size());
    EXPECT_EQ(1u, cont["labels"].size());
    EXPECT_NO_THROW(traverse(cont["sequence"].asString(), kK));
    // 0: no continuation sequence, the labels and the loss are still reported
    out = traverse(X, 0);
    const Json::Value &none = out["results"][0]["arms"]["right"]["paths"][0]["continuation"];
    EXPECT_EQ("", none["sequence"].asString());
    EXPECT_EQ(1u, none["labels"].size());
    EXPECT_EQ(0.0, none["loss_used"].asDouble());
    EXPECT_EQ(0u, out["strategy"]["output"]["continuation_bp"].asUInt64());
}


// The two gaps that depend on the request and the data rather than on a cap are stated
// too (spec §7.0): greedy re-minimisation under a finite cost and a finite branch limit,
// and column-label traces that cannot see a record boundary.
TEST(Walker, GreedyLossesAndColumnTracesAreStated) {
    // greedy losses: the ReminimisationRounds fixture (A on both branches, B on P only)
    // with a switch within budget makes the exclusion cascade through a second round
    std::vector<std::string> b;
    for (uint32_t seed = 4; ; ++seed) {
        b = clean_blocks({ 30, 25, 30, 25, 30 }, seed);
        if (b[1][0] != b[3][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Y = b[2], &Q = b[3], &Z = b[4];
    auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { X + P + Y, X + Q + Z, X + P + Y }, { "A", "A", "B" }, DeBruijnGraph::BASIC);
    auto traverse = [&](const std::string &strategy) {
        Json::Value r;
        Json::Value seed;
        seed["sequence"] = X;
        seed["labels"].append("A");
        seed["labels"].append("B");
        r["seeds"].append(seed);
        r["strategy"] = parse_json(strategy);
        return cli::process_traverse_request(r, *anno, "", cli::TraverseLimits())["results"][0];
    };
    Json::Value greedy = traverse(R"({"direction": "right",
        "labels": {"change_cost": {"model": "constant", "value": 1}, "loss_budget": 1}})");
    const Json::Value &lims = greedy["arms"]["right"]["limitations"];
    const auto greedy_kinds = kinds_of(lims);
    EXPECT_EQ(1, std::count(greedy_kinds.begin(), greedy_kinds.end(), "greedy_losses"))
        << lims.toStyledString();
    for (const auto &l : lims) {
        if (l["kind"].asString() == "greedy_losses") {
            EXPECT_EQ("branching.max_label_branches", l["knob"].asString());
            EXPECT_GE(l["observed"].asUInt64(), 1u);
        }
    }
    EXPECT_EQ("lower_bound", greedy["outcome"]["label_evidence"].asString());
    // forbid prices no switch, so its zero losses are exact; an unlimited branch limit
    // never excludes; either way nothing is stated
    for (const char *strategy : { R"({"direction": "right"})",
                                  R"({"direction": "right", "branching": {"max_label_branches": "unlimited"},
                                      "labels": {"change_cost": {"model": "constant", "value": 1}, "loss_budget": 1}})" }) {
        Json::Value r = traverse(strategy);
        const auto kinds = kinds_of(r["arms"]["right"]["limitations"]);
        EXPECT_EQ(0, std::count(kinds.begin(), kinds.end(), "greedy_losses")) << strategy;
    }

    // column traces: a coordinate index with one column F and a CoordToHeader
    std::vector<std::string> c;
    for (uint32_t seed = 29; ; ++seed) {
        c = clean_blocks({ 40, 40 }, seed);
        if (c[0].back() != c[1].back())
            break;
    }
    const std::string acc = c[0] + c[1];
    const uint64_t n = acc.size() - kK + 1;
    auto coord = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            kK, { acc }, { "F" }, DeBruijnGraph::BASIC, true, { 0 });
    std::vector<std::vector<std::string>> headers { { "acc" } };
    std::vector<std::vector<uint64_t>> num_kmers { { n } };
    annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
    LabelOracle oracle(*coord, &cth);
    for (const char *kind : { "column", "header" }) {
        Seed seed;
        seed.sequence = c[0];
        Strategy st;
        st.support = Support::TRACE;
        st.merge_reconverge = false;
        st.seed_label_kind = std::string(kind) == "column" ? LabelKind::COLUMN : LabelKind::HEADER;
        SeedResult res = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        Json::Value j = cli::seed_result_to_json(res, st, "summary", false);
        const auto kinds = kinds_of(j["limitations"]);
        EXPECT_EQ(std::string(kind) == "column" ? 1 : 0,
                  std::count(kinds.begin(), kinds.end(), "trace_record_boundaries")) << kind;
        for (const auto &d : j["seed"]["dropped_labels"]) {
            EXPECT_EQ("presence", d["runs_kind"].asString());
        }
    }
}

// Every counter belongs to the arm whose head was processed (review: Walker::run used to
// assign ALL pair evaluations to the right arm, so a left-arm greedy re-minimisation read
// pair_evaluations 0, its greedy_losses was not stated and its label_evidence read
// complete). The ReminimisationRounds fixture (two rounds under a constant cost and
// loss_budget 1) on the right, its reverse complement walked to the left, and both at
// once with a plain block on the right: each re-minimising arm states greedy_losses and
// lower_bound, and each arm's counters equal its own single-arm run's.
TEST(Walker, CountersAreAttributedToTheirArm) {
    std::vector<std::string> b;
    for (uint32_t seed = 4; ; ++seed) {
        b = clean_blocks({ 30, 25, 30, 25, 30, 30 }, seed);
        if (b[1][0] != b[3][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Y = b[2], &Q = b[3], &Z = b[4], &R = b[5];
    struct Case {
        const char *name;
        std::vector<std::string> records;
        std::string seed;
        const char *direction;
        bool left_greedy, right_greedy;
    };
    const std::vector<Case> cases {
        { "right", { X + P + Y, X + Q + Z, X + P + Y }, X, "right", false, true },
        // the mirror: walking left from rc(X) meets rc(P) / rc(Q), whose last bases differ
        { "left", { rc(Y) + rc(P) + rc(X), rc(Z) + rc(Q) + rc(X), rc(Y) + rc(P) + rc(X) },
          rc(X), "left", true, false },
        // only the left arm re-minimises; the right arm prices its own switches along R
        { "both", { rc(Y) + rc(P) + rc(X) + R, rc(Z) + rc(Q) + rc(X) + R,
                    rc(Y) + rc(P) + rc(X) + R }, rc(X), "both", true, false },
    };
    for (const Case &c : cases) {
        auto anno = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, c.records, { "A", "A", "B" }, DeBruijnGraph::BASIC);
        auto traverse = [&](const std::string &direction) {
            Json::Value r;
            Json::Value seed;
            seed["sequence"] = c.seed;
            seed["labels"].append("A");
            seed["labels"].append("B");
            r["seeds"].append(seed);
            r["strategy"] = parse_json(R"({"labels": {"change_cost": {"model": "constant",
                "value": 1}, "loss_budget": 1}})");
            r["strategy"]["direction"] = direction;
            return cli::process_traverse_request(r, *anno, "", cli::TraverseLimits())["results"][0];
        };
        const Json::Value res = traverse(c.direction);
        for (const char *side : { "left", "right" }) {
            if (!res["arms"].isMember(side))
                continue;
            const Json::Value &arm = res["arms"][side];
            const bool expected = std::string(side) == "left" ? c.left_greedy : c.right_greedy;
            const auto kinds = kinds_of(arm["limitations"]);
            EXPECT_EQ(expected ? 1 : 0, std::count(kinds.begin(), kinds.end(), "greedy_losses"))
                << c.name << ' ' << side << arm["limitations"].toStyledString();
            if (expected) {
                EXPECT_EQ(2u, arm["counters"]["max_reminimisation_rounds"].asUInt64())
                    << c.name << ' ' << side;
                EXPECT_GT(arm["counters"]["pair_evaluations"].asUInt64(), 0u) << c.name << ' ' << side;
            }
            // the arm's counters are its own: the same as when it is walked alone
            const Json::Value alone = traverse(side);
            EXPECT_EQ(alone["arms"][side]["counters"], arm["counters"]) << c.name << ' ' << side;
        }
        EXPECT_EQ("lower_bound", res["outcome"]["label_evidence"].asString()) << c.name;
        if (std::string(c.direction) == "both") {
            // the right arm priced switches of its own along R (and only those)
            EXPECT_GT(res["arms"]["right"]["counters"]["pair_evaluations"].asUInt64(), 0u);
            EXPECT_EQ(0u, res["arms"]["right"]["counters"]["max_reminimisation_rounds"].asUInt64());
        }
    }
}


} // namespace
