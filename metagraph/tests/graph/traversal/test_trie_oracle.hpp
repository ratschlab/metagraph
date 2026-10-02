#ifndef __TEST_TRIE_ORACLE_HPP__
#define __TEST_TRIE_ORACLE_HPP__

#include "gtest/gtest.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "graph/traversal/walker.hpp"


/*
 * Test-side half of the verification contract (spec §6.9):
 *
 *     T = trie(S, R, mode=annotate)                    structural, labels recorded
 *     E = { (w, l) : l in P, w a walk of T on EVERY node of which l is present,
 *                    and no longer such walk of l extends w }
 *     A = trie(S, R, mode=constrain, P, cost=forbid)   the label-constrained walker
 *     claims(A) == E, both cut to min(complete_to_bp)
 *
 * The unit of comparison is a (walk, label) CLAIM — "label l carries exactly these
 * bases from the seed boundary" — not a leaf: a label whose maximal label-consistent
 * walk ends INSIDE a walk that other labels continue (S stopping in Y while A and C
 * go on) is a claim of its own, and where it ends is exactly what the walker has to
 * get right per label. Comparing leaves only would never check it (the leaf lists
 * the labels that reach it, and S does not). The claims of A are its label ends:
 * every LABEL_END event at its depth, which includes every leaf's end_labels.
 *
 * E is computed from the RECORDED label sets with plain code that shares nothing with
 * the walker's label state: the only walker outputs it reads are segment sequences, the
 * recorded sets (Segment::label_sets, the root's labels_start and its total) and the
 * path / segment structure. Walks are compared in WALKING order, outward from the seed
 * boundary, so that a prefix is a prefix on both arms.
 *
 * The second half is the tuned-run property: every claim of a tuned run (every label
 * end, leaf or not) is a prefix of an exhaustive walk under its label and a semantic
 * end is one the exhaustive run makes too, and every (walk, label) of the exhaustive
 * trie that the tuned run does not reach has a recorded reason — about THAT label and
 * THAT successor — where the tuned run left it. Nothing is dropped silently, nothing
 * is invented.
 *
 * The checkers of the second half collect every violation as a line of text: the core
 * of each (tuned_subset_report, routes_subset_report) returns the report, the EXPECT
 * wrapper (check_tuned_subset, check_routes_subset) asserts that it is empty, and a
 * test can corrupt a result on purpose and assert that the report is NOT empty — the
 * checkers are themselves under test (review round 2: three accepted corruptions).
 */
namespace mtg {
namespace test {
namespace trie {

using namespace mtg::graph::traversal;

// a set of (walk, label) claims, grouped by walk: walk (outward from the seed
// boundary) -> names of the labels claiming exactly it. The leaves of a run are the
// claims whose walk is a leaf's.
using Leaves = std::map<std::string, std::set<std::string>>;

// the violations a checker found, one line each; empty iff the result passes
using Problems = std::vector<std::string>;

// a segment's or path's bases in walking order (the left arm spells them in natural
// orientation, i.e. reversed)
inline std::string outward(const ArmResult &arm, std::string s) {
    if (arm.arm == Arm::LEFT)
        std::reverse(s.begin(), s.end());
    return s;
}

inline std::string walk_of(const ArmResult &arm, const PathResult &path) {
    return outward(arm, spell_path(arm, path));
}

inline std::set<std::string> name_set(const SeedResult &r, const std::vector<LabelId> &ids) {
    std::set<std::string> out;
    for (LabelId l : ids)
        out.insert(r.label_dict.at(l).name);
    return out;
}

inline std::optional<LabelId> label_id(const SeedResult &r, const std::string &name) {
    for (LabelId l = 0; l < r.label_dict.size(); ++l) {
        if (r.label_dict[l].name == name)
            return l;
    }
    return std::nullopt;
}

inline const Segment& root_of(const ArmResult &arm) {
    static const Segment none;
    const Segment *root = nullptr;
    for (const Segment &seg : arm.segments) {
        if (seg.parents.empty()) {
            EXPECT_EQ(nullptr, root) << "two root segments";
            root = &seg;
        }
    }
    EXPECT_NE(nullptr, root) << "no root segment";
    return root ? *root : none;
}

// no joins: the segment DAG is a trie (on_reconverge: keep)
inline bool is_trie(const ArmResult &arm) {
    for (const Segment &seg : arm.segments) {
        if (seg.parents.size() > 1)
            return false;
    }
    return true;
}

// a recorded non-semantic end: the walk was cut, not exhausted
inline bool is_cut(EndReason reason) {
    return reason == EndReason::BRANCH || reason == EndReason::LOSS_BUDGET
        || reason == EndReason::BEAM_PRUNED || is_resource_stop(reason);
}

// a censoring end: a size cap or the beam, which the completeness guarantee rules out
// inside the complete region (a branch-limit or quorum end is a reported decision,
// not a censoring)
inline bool is_censored(EndReason reason) {
    return reason == EndReason::BEAM_PRUNED || is_resource_stop(reason);
}

// a structural block: the walk rule itself refused the step
inline bool is_block(EndReason reason) {
    return reason == EndReason::EDGE_REUSE || reason == EndReason::EDGE_REUSE_RC
        || reason == EndReason::REACHED_SEED;
}

inline std::string shown(const std::string &w) { return w.empty() ? "<empty>" : w; }

inline std::string listed(const Problems &problems) {
    std::string s;
    for (const std::string &p : problems)
        s += "\n  - " + p;
    return s;
}

inline std::string describe(const Leaves &leaves) {
    std::ostringstream os;
    for (const auto &[w, ls] : leaves) {
        os << "\n  " << shown(w) << " (" << w.size() << " bp) {";
        for (const auto &l : ls) os << ' ' << l;
        os << " }";
    }
    return os.str();
}


/********************************* the oracle E *********************************/

// The labels present at every node of |path| (index 0 = the seed boundary node), read
// off the recorded sets of an ANNOTATE run. Sets |cut| when any recorded list was cut
// by max_labels_per_node: the oracle is then unusable for exact label claims.
inline std::vector<std::set<std::string>>
recorded_labels(const SeedResult &T, const ArmResult &arm, const PathResult &path, bool *cut) {
    std::vector<std::set<std::string>> at(path.length_bp + 1);
    std::vector<bool> filled(path.length_bp + 1, false);
    at[0] = name_set(T, root_of(arm).labels_start);
    filled[0] = true;
    for (size_t s : path.segments) {
        const Segment &seg = arm.segments[s];
        for (const LabelSetRun &run : seg.label_sets) {
            if (run.truncated())
                *cut = true;
            const std::set<std::string> names = name_set(T, run.labels);
            for (uint64_t d = run.from_bp + 1; d <= run.to_bp; ++d) {
                EXPECT_LE(d, path.length_bp) << "a run past the end of path " << path.id;
                if (d > path.length_bp)
                    break;
                EXPECT_FALSE(filled[d]) << "overlapping runs at depth " << d;
                at[d] = names;
                filled[d] = true;
            }
        }
    }
    for (size_t d = 0; d <= path.length_bp; ++d) {
        EXPECT_TRUE(filled[d]) << "no recorded set at depth " << d << " of path " << path.id;
    }
    return at;
}

struct Oracle {
    Leaves leaves;      // E: the per-label maximal (walk, label) claims
    bool cut = false;   // some recorded list was truncated: the labels are unreliable
};

// E for the permitted set |permitted|, cut to |depth|: for every walk of T and every
// permitted label, the longest prefix on which the label is present at every node;
// PER LABEL, the maximal ones among those are its claims. (Maximal per label, not per
// walk: a label's longest prefix may be a proper prefix of another label's, and it is
// a claim of its own — the walker has to end that label exactly there.)
inline Oracle expected_leaves(const SeedResult &T, size_t a,
                              const std::set<std::string> &permitted, uint64_t depth) {
    const ArmResult &arm = T.arms[a];
    Oracle out;
    EXPECT_TRUE(is_trie(arm)) << "the structural run is not a trie (merging on?)";
    // the root's boundary list is in no run; its own total tells whether it was cut
    // (the contract: a cut list must not be read as "these labels and no other")
    const Segment &root = root_of(arm);
    if (root.labels_start_total > root.labels_start.size())
        out.cut = true;
    std::map<std::string, std::set<std::string>> prefixes;   // label -> its prefixes
    for (const PathResult &path : arm.paths) {
        EXPECT_TRUE(path.end_labels.empty()) << "not an annotate run: path " << path.id
                                             << " has label ends";
        const std::string w = walk_of(arm, path);
        EXPECT_EQ(path.length_bp, w.size());
        const auto at = recorded_labels(T, arm, path, &out.cut);
        for (const std::string &l : permitted) {
            if (!at[0].count(l))
                continue;   // not even on the boundary node: cannot carry the seed
            uint64_t j = 0;
            while (j < w.size() && j < depth && at[j + 1].count(l))
                ++j;
            prefixes[l].insert(w.substr(0, j));
        }
    }
    // keep each label's maximal prefixes: in lexicographic order, every extension of
    // s directly follows s, so s is maximal iff its successor does not start with it
    for (const auto &[l, ps] : prefixes) {
        for (auto it = ps.begin(); it != ps.end(); ++it) {
            auto next = std::next(it);
            const bool extended = next != ps.end() && next->size() > it->size()
                && next->compare(0, it->size(), *it) == 0;
            if (!extended)
                out.leaves[*it].insert(l);
        }
    }
    return out;
}


/**************************** the walker's leaves A ****************************/

// Labels alive at |depth| along |path| of a CONSTRAIN run under cost forbid. A label
// can leave a path without a label end (at a split it may follow only the other
// child), so this reads the segment that holds the depth: its labels at entry
// (Segment::labels_start, the seed boundary for the root) minus the lineages that
// ended inside it before |depth|, or its labels at the last node when |depth| is it.
inline std::set<std::string> alive_at(const SeedResult &r, const ArmResult &arm,
                                      const PathResult &path, uint64_t depth) {
    for (size_t s : path.segments) {
        const Segment &seg = arm.segments[s];
        const uint64_t end = seg.from_bp + seg.length_bp;
        if (depth > end)
            continue;
        if (depth == end)
            return name_set(r, seg.labels_end);
        if (depth == seg.from_bp)   // the root's own boundary node
            return name_set(r, seg.labels_start);
        std::set<std::string> alive = name_set(r, seg.labels_start);
        for (const Event &ev : seg.events) {
            EXPECT_NE(EventType::SWITCH, ev.type) << "cost forbid expected";
            if (ev.type == EventType::LABEL_END && ev.at_bp < depth)
                alive.erase(r.label_dict.at(ev.label).name);
        }
        return alive;
    }
    ADD_FAILURE() << "depth " << depth << " is beyond path " << path.id;
    return {};
}

// The leaves of a constrain run, cut to |depth|: a leaf longer than that is cut and
// carries the labels alive at the cut. Inside the arm's guaranteed region every end
// must be semantic.
inline Leaves constrained_leaves(const SeedResult &A, size_t a, uint64_t depth) {
    const ArmResult &arm = A.arms[a];
    EXPECT_TRUE(is_trie(arm)) << "the constrained run is not a trie (merging on?)";
    Leaves out;
    for (const PathResult &path : arm.paths) {
        std::string w = walk_of(arm, path);
        EXPECT_EQ(path.length_bp, w.size());
        std::set<std::string> labels;
        if (w.size() <= depth) {
            for (const LabelEnd &e : path.end_labels) {
                EXPECT_EQ(0u, e.route_bp) << "a trie has no merged-in labels";
                labels.insert(A.label_dict.at(e.label).name);
            }
            if (w.size() < arm.complete_to_bp) {
                if (path.path_reason) {
                    EXPECT_FALSE(is_censored(*path.path_reason))
                        << "path " << path.id << " (" << w.size() << " bp) ended with "
                        << to_string(*path.path_reason) << " inside the guaranteed region ("
                        << arm.complete_to_bp << " bp)";
                }
                for (size_t r = 0; r < kNumEndReasons; ++r) {
                    if (path.end_reasons[r]) {
                        EXPECT_FALSE(is_censored(static_cast<EndReason>(r)))
                            << "path " << path.id << " has a label ending with "
                            << to_string(static_cast<EndReason>(r))
                            << " inside the guaranteed region";
                    }
                }
            }
        } else {
            labels = alive_at(A, arm, path, depth);
            w.resize(depth);
        }
        auto [it, inserted] = out.emplace(w, labels);
        if (!inserted) {
            EXPECT_EQ(it->second, labels)
                << "two paths through the walk " << w
                << " disagree on the labels alive at depth " << depth;
        }
    }
    return out;
}


// forward: the (walk, label) ends of a constrain run, defined with the tuned-run
// helpers below
inline std::map<std::string, std::map<std::string, std::optional<EndReason>>>
label_claims(const SeedResult &r, size_t a, uint64_t depth);

// The claims of a constrain run cut to |depth|: every (walk prefix, label) at which
// it ends a label — the leaves' labels (constrained_leaves, whose checks run too) plus
// every label ended INSIDE a walk that goes on under other labels.
inline Leaves constrained_claims(const SeedResult &A, size_t a, uint64_t depth) {
    Leaves out = constrained_leaves(A, a, depth);
    for (const auto &[w, labels] : label_claims(A, a, depth)) {
        for (const auto &[l, reason] : labels) {
            out[w].insert(l);
        }
    }
    return out;
}


/********************************* the contract *********************************/

struct Verdict {
    uint64_t depth = 0;   // min(complete_to_bp) of the two runs
    bool cut = false;     // the oracle's lists were truncated: it proves nothing
    Leaves expected;      // E, the per-label maximal (walk, label) claims
    Leaves actual;        // claims(A)
};

inline Verdict verify(const SeedResult &T, const SeedResult &A, size_t a,
                      const std::set<std::string> &permitted) {
    Verdict v;
    v.depth = std::min(T.arms[a].complete_to_bp, A.arms[a].complete_to_bp);
    // a derived permitted set cut by max_seed_labels would make A a trie over FEWER
    // labels than E filters by, and the two would disagree for a reason that is not
    // the walker's; the exhaustive preset refuses that cut, so it never shows here
    EXPECT_EQ(0u, A.labels_dropped) << "the constrained run's permitted set was cut";
    EXPECT_TRUE(A.dropped_labels.empty()) << "the constrained run dropped seed labels";
    Oracle o = expected_leaves(T, a, permitted, v.depth);
    v.cut = o.cut;
    v.expected = std::move(o.leaves);
    v.actual = constrained_claims(A, a, v.depth);
    return v;
}

inline void expect_equal(const Leaves &expected, const Leaves &actual, const std::string &what) {
    std::ostringstream diff;
    auto list = [](const std::set<std::string> &ls) {
        std::string s = "{";
        for (const auto &l : ls) s += " " + l;
        return s + " }";
    };
    for (const auto &[w, ls] : expected) {
        auto it = actual.find(w);
        if (it == actual.end()) {
            diff << "\n  only the structural oracle has " << shown(w)
                 << " (" << w.size() << " bp) " << list(ls);
        } else if (it->second != ls) {
            diff << "\n  labels differ on " << shown(w) << " ("
                 << w.size() << " bp): oracle " << list(ls) << ", walker "
                 << list(it->second);
        }
    }
    for (const auto &[w, ls] : actual) {
        if (!expected.count(w)) {
            diff << "\n  only the constrained walker has " << shown(w)
                 << " (" << w.size() << " bp) " << list(ls);
        }
    }
    EXPECT_TRUE(diff.str().empty())
        << what << ": the structural oracle and the label-constrained walker disagree:"
        << diff.str();
}


/****************************** the tuned-run property ******************************/

inline const PathResult* leaf_path(const ArmResult &arm, size_t segment) {
    for (const PathResult &p : arm.paths) {
        if (!p.segments.empty() && p.segments.back() == segment)
            return &p;
    }
    return nullptr;
}

// the bases the trie continues with at depth |at| of |segment|: inside the segment its
// own next base, at its end the first base of every child
inline std::set<char> continuations(const ArmResult &arm, size_t segment, uint64_t at) {
    const Segment &seg = arm.segments[segment];
    const std::string seq = outward(arm, seg.sequence);
    std::set<char> out;
    if (at < seg.from_bp)
        return out;
    const uint64_t j = at - seg.from_bp;
    if (j < seq.size()) {
        out.insert(seq[j]);
    } else if (j == seq.size()) {
        for (size_t c : seg.children) {
            const std::string cs = outward(arm, arm.segments[c].sequence);
            if (!cs.empty())
                out.insert(cs[0]);
        }
    }
    return out;
}

// Discard evidence for label |name| not taking base |ch| at depth |at| of |segment| of
// a constrain run |r| — evidence about THAT label and THAT successor:
//  - a branch event there listing |ch| among the labelled admissible successors whose
//    |dropped| names the label: it was ended there by the branch limit, a quorum, the
//    split limit or the loss budget; or whose |ambiguous| names the label while the
//    trie does not follow |ch| there: the label went on along another successor and
//    |ch|, labelled and admissible, was refused by a quorum or the split limit (a
//    labelled admissible successor is followed otherwise, and a label that lost every
//    successor is in |dropped|);
//  - a BLOCKED event there for |ch| naming the label among the labels that would have
//    continued on it, with a structural reason (edge_reuse, edge_reuse_rc,
//    rejoined_seed);
//  - a HAIRPIN event there for |ch| naming the label: the step was skipped.
// A branch event alone proves nothing: an ambiguity that was FOLLOWED on every
// successor emits one too (|dropped| empty), and a result from which a subtree was
// deleted still carries it (review round 2, finding 1).
inline bool branch_recorded(const SeedResult &r, const ArmResult &arm, size_t segment,
                            uint64_t at, char ch, const std::string &name) {
    const std::optional<LabelId> id = label_id(r, name);
    if (!id)
        return false;
    auto names = [&](const std::vector<LabelId> &ids) {
        return std::find(ids.begin(), ids.end(), *id) != ids.end();
    };
    const bool followed = continuations(arm, segment, at).count(ch) > 0;
    for (const BranchEvent &be : arm.branch_events) {
        if (be.segment != segment || be.at_bp != at
                || std::find(be.chars.begin(), be.chars.end(), ch) == be.chars.end())
            continue;
        if (names(be.dropped))
            return true;
        if (!followed && names(be.ambiguous))
            return true;
    }
    for (const Event &ev : arm.segments[segment].events) {
        if (ev.at_bp != at || ev.ch != ch || !names(ev.labels))
            continue;
        if (ev.type == EventType::BLOCKED && is_block(ev.reason))
            return true;
        if (ev.type == EventType::HAIRPIN && !followed)
            return true;
    }
    return false;
}

inline const Event* label_end_at(const SeedResult &r, const Segment &seg,
                                 const std::string &name, uint64_t at) {
    for (const Event &ev : seg.events) {
        if (ev.type == EventType::LABEL_END && ev.at_bp == at
                && r.label_dict.at(ev.label).name == name)
            return &ev;
    }
    return nullptr;
}

struct SubsetReport {
    size_t present = 0;   // (walk, label) claims of A the tuned run reaches in full
    size_t omitted = 0;   // ... it leaves early, each with a recorded reason
    Problems problems;    // every violation found; empty iff the tuned run passes
};

// The claims of a constrain run: every (walk prefix, label) at which the run ends a
// label, with the reason — i.e. every stretch the run says the label carries, leaf or
// not (a label that stops inside a walk that goes on under other labels is a claim
// too). Cut to |depth|: a label alive at the cut is claimed up to it, reason unset.
// (A label that leaves a path at a split by following only the other child has no
// end there and makes no claim on this path; its claim is on the child it took.)
inline std::map<std::string, std::map<std::string, std::optional<EndReason>>>
label_claims(const SeedResult &r, size_t a, uint64_t depth) {
    const ArmResult &arm = r.arms[a];
    std::map<std::string, std::map<std::string, std::optional<EndReason>>> out;
    for (const PathResult &p : arm.paths) {
        const std::string w = walk_of(arm, p);
        std::set<std::string> ended;
        for (size_t s : p.segments) {
            for (const Event &ev : arm.segments[s].events) {
                if (ev.type != EventType::LABEL_END)
                    continue;
                const std::string &l = r.label_dict.at(ev.label).name;
                ended.insert(l);
                if (ev.at_bp <= depth) {
                    out[w.substr(0, ev.at_bp)][l] = ev.reason;
                } else {
                    out[w.substr(0, depth)].emplace(l, std::nullopt);
                }
            }
        }
        if (w.size() > depth) {
            for (const std::string &l : alive_at(r, arm, p, depth))
                out[w.substr(0, depth)].emplace(l, std::nullopt);
        }
    }
    return out;
}

// |tuned| must be a trie (merging off). Checks, against the exhaustive run |A|:
//  1. nothing is invented: every CLAIM of the tuned run — every (walk prefix, label) at
//     which it ends a label, leaf or not, and every label alive at the cut — is a
//     prefix of an exhaustive walk on which that label is alive at that depth; and a
//     claim ended for a SEMANTIC reason (not a cut) is one the exhaustive run makes
//     too, i.e. it ends the label at exactly that position (the leaves alone would
//     miss a label end inserted inside a segment: review round 2, finding 2);
//  2. every omission has a reason: for every claim (walk prefix, label) of A, following
//     the prefix through the tuned trie either reaches its end with the label alive
//     and ended there too (present, for the same reason unless the tuned run ran into
//     a cap exactly there), or stops earlier at a RECORDED reason: a label end by the
//     branch limit, a quorum, the split limit, the loss budget, a cap or the beam, or
//     discard evidence for the label and the base not taken (branch_recorded). A label
//     ended with a SEMANTIC reason (label_lost, dead_end, ...) where A continues it,
//     or a label that vanishes without any record, is a disagreement between the two
//     runs and fails.
// The report lists every violation; check_tuned_subset() asserts that it is empty.
inline SubsetReport tuned_subset_report(const SeedResult &A, const SeedResult &tuned,
                                        size_t a, const std::string &what) {
    const ArmResult &ta = tuned.arms[a];
    const ArmResult &aa = A.arms[a];
    SubsetReport rep;
    Problems &problems = rep.problems;
    if (!is_trie(ta))
        problems.push_back(what + ": the tuned run is not a trie (merging on)");
    const uint64_t depth = aa.complete_to_bp;
    const auto a_claims = label_claims(A, a, depth);
    const Leaves a_all = constrained_claims(A, a, depth);
    auto bp = [](uint64_t n) { return std::to_string(n); };

    // 1. nothing is invented
    const auto t_claims = label_claims(tuned, a, depth);
    for (const auto &[w, labels] : constrained_claims(tuned, a, depth)) {
        for (const std::string &l : labels) {
            bool found = false;
            for (const PathResult &ap : aa.paths) {
                const std::string aw = walk_of(aa, ap);
                if (aw.size() < w.size() || aw.compare(0, w.size(), w) != 0)
                    continue;
                if (alive_at(A, aa, ap, w.size()).count(l)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                problems.push_back(what + ": tuned claim " + shown(w) + " (" + bp(w.size())
                                   + " bp) under " + l + " is not a prefix of an exhaustive "
                                     "walk carrying that label to that depth: INVENTED");
                continue;
            }
            std::optional<EndReason> reason;
            if (auto it = t_claims.find(w); it != t_claims.end()) {
                if (auto jt = it->second.find(l); jt != it->second.end())
                    reason = jt->second;
            }
            if (reason && !is_cut(*reason)) {
                auto it = a_all.find(w);
                if (it == a_all.end() || !it->second.count(l)) {
                    problems.push_back(what + ": the tuned run ends " + l + " at " + bp(w.size())
                                       + " on " + shown(w) + " with " + to_string(*reason)
                                       + " while the exhaustive trie continues it: an INVENTED end");
                }
            }
        }
    }

    // 2. every omission has a reason
    const Segment &root = root_of(ta);
    const std::set<std::string> permitted = name_set(tuned, root.labels_start);
    for (const auto &[w, labels] : a_claims) {
        for (const auto &[l, a_reason] : labels) {
            if (!permitted.count(l)) {
                bool dropped = false;
                for (const DroppedLabel &d : tuned.dropped_labels)
                    dropped |= d.name == l;
                if (!dropped) {
                    problems.push_back(what + ": label " + l
                                       + " is neither permitted nor reported dropped");
                }
                rep.omitted++;
                continue;
            }
            // follow w through the tuned trie with l alive
            const Segment *seg = &root;
            bool present = false, decided = false;
            while (!decided) {
                const std::string seq = outward(ta, seg->sequence);
                const Segment *next = nullptr;
                for (size_t j = 0; !decided; ++j) {
                    const uint64_t d = seg->from_bp + j;
                    if (const Event *end = label_end_at(tuned, *seg, l, d)) {
                        if (d == w.size()) {
                            // both runs end the label here: the reasons agree unless
                            // the tuned run ran into a cap exactly there
                            if (!is_cut(end->reason) && a_reason
                                    && std::string(to_string(*a_reason)) != to_string(end->reason)) {
                                problems.push_back(what + ": " + l + " ends at " + bp(d) + " on "
                                                   + shown(w) + " for different reasons: exhaustive "
                                                   + to_string(*a_reason) + ", tuned "
                                                   + to_string(end->reason));
                            }
                            present = true;
                        } else if (!is_cut(end->reason)) {
                            problems.push_back(what + ": the tuned run ends " + l + " at " + bp(d)
                                               + " on " + shown(w) + " with " + to_string(end->reason)
                                               + (end->text.empty() ? "" : " (" + end->text + ")")
                                               + " while the exhaustive trie continues it: the two "
                                                 "walkers DISAGREE");
                        }
                        decided = true;
                        break;
                    }
                    if (d == w.size()) {
                        // the exhaustive claim was itself cut at the boundary: nothing
                        // is claimed beyond it. Otherwise the tuned run outlives it.
                        if (a_reason) {
                            problems.push_back(what + ": label " + l + " outlives the exhaustive walk "
                                               + shown(w) + " (" + bp(w.size()) + " bp, ended there with "
                                               + to_string(*a_reason) + ") in the tuned run");
                        }
                        present = true;
                        decided = true;
                        break;
                    }
                    if (j < seq.size()) {
                        if (seq[j] == w[d])
                            continue;   // the tuned run follows w
                        if (!branch_recorded(tuned, ta, seg->id, d, w[d], l)) {
                            problems.push_back(what + ": the tuned run leaves " + shown(w) + " at "
                                               + bp(d) + " (takes " + seq[j] + ", not " + w[d]
                                               + ") with " + l + " alive and no recorded reason");
                        }
                        decided = true;
                        break;
                    }
                    // the segment ends at depth d: the child continuing w, if any
                    for (size_t c : seg->children) {
                        const std::string cs = outward(ta, ta.segments[c].sequence);
                        if (!cs.empty() && cs[0] == w[d])
                            next = &ta.segments[c];
                    }
                    if (next) {
                        // the label must enter the child it is said to follow: a label
                        // that vanishes at a split without an end was dropped silently
                        if (!name_set(tuned, next->labels_start).count(l)) {
                            problems.push_back(what + ": label " + l + " is on " + shown(w) + " at "
                                               + bp(d) + " and the tuned run takes " + w[d]
                                               + ", but the label is not on that child and nothing "
                                                 "recorded its end: dropped SILENTLY");
                            decided = true;
                        }
                        break;
                    }
                    if (seg->children.empty()) {
                        const PathResult *p = leaf_path(ta, seg->id);
                        const bool ok = p && p->path_reason && is_cut(*p->path_reason);
                        if (!ok) {
                            problems.push_back(what + ": tuned leaf " + shown(w.substr(0, d)) + " ("
                                               + bp(d) + " bp) ends with " + l
                                               + " alive and no recorded reason while the "
                                                 "exhaustive trie continues to " + shown(w));
                        }
                    } else if (!branch_recorded(tuned, ta, seg->id, d, w[d], l)) {
                        problems.push_back(what + ": the tuned run drops " + w[d] + " at " + bp(d)
                                           + " on " + shown(w) + " with " + l
                                           + " alive and no recorded reason");
                    }
                    decided = true;
                }
                if (!decided)
                    seg = next;
            }
            if (present)
                rep.present++;
            else
                rep.omitted++;
        }
    }
    return rep;
}

inline SubsetReport check_tuned_subset(const SeedResult &A, const SeedResult &tuned,
                                       size_t a, const std::string &what) {
    SubsetReport rep = tuned_subset_report(A, tuned, a, what);
    EXPECT_TRUE(rep.problems.empty())
        << what << ": the tuned run fails the prefix-subset property ("
        << rep.problems.size() << " problem(s)):" << listed(rep.problems);
    return rep;
}

// The flank of |label| ending at segment |leaf| in walking order, following at each
// join the parent the label entered through (Segment::labels_via_parent). A label that
// enters a join through no listed parent goes into |problems| when given, else fails.
inline std::string label_route(const ArmResult &arm, size_t leaf, LabelId label,
                               Problems *problems = nullptr) {
    std::string rev;   // leaf -> root
    size_t s = leaf;
    while (true) {
        const Segment &seg = arm.segments[s];
        const std::string piece = outward(arm, seg.sequence);
        rev.append(piece.rbegin(), piece.rend());
        if (seg.parents.empty())
            break;
        size_t next = seg.parents[0];
        if (seg.parents.size() > 1) {
            bool found = false;
            for (size_t p = 0; p < seg.labels_via_parent.size() && !found; ++p) {
                const auto &via = seg.labels_via_parent[p];
                if (std::find(via.begin(), via.end(), label) != via.end()) {
                    next = seg.parents[p];
                    found = true;
                }
            }
            if (!found) {
                const std::string msg = "label " + std::to_string(label) + " enters segment "
                    + std::to_string(s) + " through no parent";
                if (problems) {
                    problems->push_back(msg);
                } else {
                    ADD_FAILURE() << msg;
                }
            }
        }
        s = next;
    }
    return std::string(rev.rbegin(), rev.rend());
}

struct RouteReport {
    size_t checked = 0;   // label ends checked: every leaf's labels and every interior end
    Problems problems;    // every violation found; empty iff the merged run passes
};

// Merging on: a leaf's spelled walk is not evidence for a label merged in at a join
// (LabelEnd::route_bp > 0) — the label's own route is. EVERY label end of the merged
// run is checked on its own route: the leaves' labels and every LABEL_END inside a
// segment (a leaf-only check would miss an end inserted inside a segment: review round
// 2, finding 2). The route must be a prefix of an exhaustive CLAIM under that label (a
// claim, not a leaf: the exhaustive run may end the label inside a walk other labels go
// on with), and an end for a reason that does not depend on the edge history
// (dead_end, label_lost, rejoined_seed, the radius, a trace break) must be where the
// exhaustive run ends the label too, for the same reason. An edge-reuse end may come
// earlier under merging (§6.10: a merge unites the edge histories of the routes it
// joins), so for it the prefix condition is all that is required. A merge closes a
// duplicate lineage without a label end, so every LABEL_END event is a real end.
// The report lists every violation; check_routes_subset() asserts that it is empty.
inline RouteReport routes_subset_report(const SeedResult &A, const SeedResult &merged,
                                        size_t a, const std::string &what) {
    const ArmResult &ma = merged.arms[a];
    const uint64_t depth = A.arms[a].complete_to_bp;
    const Leaves full = constrained_claims(A, a, depth);
    const auto a_claims = label_claims(A, a, depth);
    RouteReport rep;
    Problems &problems = rep.problems;
    auto bp = [](uint64_t n) { return std::to_string(n); };
    auto supported = [&](const std::string &route, const std::string &l) {
        for (auto it = full.lower_bound(route);
                it != full.end() && it->first.compare(0, route.size(), route) == 0; ++it) {
            if (it->second.count(l))
                return true;
        }
        return false;
    };
    // one label end: |route| the label's own route to it, |reason| how it ended there
    auto check_end = [&](std::string route, const std::string &l,
                         std::optional<EndReason> reason, const std::string &where) {
        rep.checked++;
        if (route.size() > depth) {
            route.resize(depth);
            reason.reset();
        }
        if (!supported(route, l)) {
            problems.push_back(what + ": the route of " + l + " to " + where + " (" + shown(route)
                               + ", " + bp(route.size()) + " bp) is not a prefix of an exhaustive "
                                 "claim under that label: INVENTED");
            return;
        }
        if (!reason || is_cut(*reason) || *reason == EndReason::EDGE_REUSE
                || *reason == EndReason::EDGE_REUSE_RC)
            return;
        auto it = full.find(route);
        if (it == full.end() || !it->second.count(l)) {
            problems.push_back(what + ": the merged run ends " + l + " at " + where + " ("
                               + shown(route) + ", " + bp(route.size()) + " bp) with "
                               + to_string(*reason) + " while the exhaustive trie continues it: "
                                 "an INVENTED end");
            return;
        }
        if (auto jt = a_claims.find(route); jt != a_claims.end()) {
            auto kt = jt->second.find(l);
            if (kt != jt->second.end() && kt->second
                    && std::string(to_string(*kt->second)) != to_string(*reason)) {
                problems.push_back(what + ": " + l + " ends at " + where + " for different reasons: "
                                   "exhaustive " + to_string(*kt->second) + ", merged "
                                   + to_string(*reason));
            }
        }
    };

    // the leaves: a leaf label's route and the spelled walk share everything from
    // route_bp on, and the leaf lists exactly the labels ended at its depth
    for (const PathResult &p : ma.paths) {
        const std::string spelled = walk_of(ma, p);
        if (p.segments.empty()) {
            problems.push_back(what + ": path " + bp(p.id) + " has no segments");
            continue;
        }
        const Segment &last = ma.segments[p.segments.back()];
        for (const LabelEnd &e : p.end_labels) {
            const std::string &l = merged.label_dict.at(e.label).name;
            const std::string route = label_route(ma, p.segments.back(), e.label, &problems);
            if (spelled.size() != route.size()) {
                problems.push_back(what + ": leaf " + bp(p.id) + ": the route of " + l + " spells "
                                   + bp(route.size()) + " bases, the walk " + bp(spelled.size()));
            }
            const size_t from = std::min<size_t>(e.route_bp, std::min(spelled.size(), route.size()));
            if (spelled.substr(from) != route.substr(from)) {
                problems.push_back(what + ": leaf " + bp(p.id) + ": the route of " + l
                                   + " and the spelled walk differ after route_bp " + bp(e.route_bp));
            }
            if (e.route_bp == 0 && spelled != route) {
                problems.push_back(what + ": leaf " + bp(p.id) + ": " + l + " has route_bp 0 but its "
                                     "route is not the spelled walk");
            } else if (e.route_bp > 0 && spelled == route) {
                problems.push_back(what + ": " + l + " is reported merged in at " + bp(e.route_bp)
                                   + " but its route is the spelled walk");
            }
            if (!label_end_at(merged, last, l, spelled.size())) {
                problems.push_back(what + ": leaf " + bp(p.id) + " lists " + l
                                   + " with no label end event at its depth");
            }
        }
    }
    // every label end, leaf or interior, on the label's own route to it
    for (const Segment &seg : ma.segments) {
        const uint64_t end = seg.from_bp + seg.length_bp;
        for (const Event &ev : seg.events) {
            if (ev.type != EventType::LABEL_END)
                continue;
            const std::string &l = merged.label_dict.at(ev.label).name;
            const std::string where = "segment " + bp(seg.id) + " at " + bp(ev.at_bp);
            std::string route = label_route(ma, seg.id, ev.label, &problems);
            if (route.size() != end || ev.at_bp > end || ev.at_bp < seg.from_bp) {
                problems.push_back(what + ": the route of " + l + " to " + where + " spells "
                                   + bp(route.size()) + " bases for a segment over ["
                                   + bp(seg.from_bp) + ", " + bp(end) + ")");
                continue;
            }
            route.resize(ev.at_bp);
            check_end(route, l, ev.reason, where);
        }
    }
    return rep;
}

inline size_t check_routes_subset(const SeedResult &A, const SeedResult &merged,
                                  size_t a, const std::string &what) {
    RouteReport rep = routes_subset_report(A, merged, a, what);
    EXPECT_TRUE(rep.problems.empty())
        << what << ": the merged run's label ends are not all supported by the exhaustive trie ("
        << rep.problems.size() << " problem(s)):" << listed(rep.problems);
    return rep.checked;
}

} // namespace trie
} // namespace test
} // namespace mtg

#endif // __TEST_TRIE_ORACLE_HPP__
