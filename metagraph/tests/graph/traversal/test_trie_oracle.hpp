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
 * The second half is the tuned-run property: every leaf of a tuned run is a prefix of
 * an exhaustive leaf under each of its labels, and every (walk, label) of the
 * exhaustive trie that the tuned run does not reach has a recorded reason where the
 * tuned run left it. Nothing is dropped silently.
 */
namespace mtg {
namespace test {
namespace trie {

using namespace mtg::graph::traversal;

// a set of (walk, label) claims, grouped by walk: walk (outward from the seed
// boundary) -> names of the labels claiming exactly it. The leaves of a run are the
// claims whose walk is a leaf's.
using Leaves = std::map<std::string, std::set<std::string>>;

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

inline std::string describe(const Leaves &leaves) {
    std::ostringstream os;
    for (const auto &[w, ls] : leaves) {
        os << "\n  " << (w.empty() ? "<empty>" : w) << " (" << w.size() << " bp) {";
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
            diff << "\n  only the structural oracle has " << (w.empty() ? "<empty>" : w)
                 << " (" << w.size() << " bp) " << list(ls);
        } else if (it->second != ls) {
            diff << "\n  labels differ on " << (w.empty() ? "<empty>" : w) << " ("
                 << w.size() << " bp): oracle " << list(ls) << ", walker "
                 << list(it->second);
        }
    }
    for (const auto &[w, ls] : actual) {
        if (!expected.count(w)) {
            diff << "\n  only the constrained walker has " << (w.empty() ? "<empty>" : w)
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

// a dropped continuation |ch| at depth |at| of |segment| is on record: a branch event
// listing it (quorum, split limit, ambiguity) or a blocked-successor event
inline bool branch_recorded(const ArmResult &arm, size_t segment, uint64_t at, char ch) {
    for (const BranchEvent &be : arm.branch_events) {
        if (be.segment == segment && be.at_bp == at
                && std::find(be.chars.begin(), be.chars.end(), ch) != be.chars.end())
            return true;
    }
    for (const Event &ev : arm.segments[segment].events) {
        if (ev.type == EventType::BLOCKED && ev.at_bp == at && ev.ch == ch)
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
//  1. every tuned leaf, under each of its labels, is a prefix of an exhaustive walk on
//     which that label is alive at the leaf's depth (nothing is invented);
//  2. for every claim (walk prefix, label) of A, following the prefix through the
//     tuned trie either reaches its end with the label alive and ended there too
//     (present), or stops earlier at a RECORDED reason: a label end by the branch
//     limit, a quorum, the split limit, the loss budget, a cap or the beam, or a
//     branch event naming the base not taken. A label ended with a SEMANTIC reason
//     (label_lost, dead_end, ...) where A continues it, or a label that vanishes
//     without any record, is a disagreement between the two runs and fails.
inline SubsetReport check_tuned_subset(const SeedResult &A, const SeedResult &tuned,
                                       size_t a, const std::string &what) {
    const ArmResult &ta = tuned.arms[a];
    const ArmResult &aa = A.arms[a];
    SubsetReport rep;
    EXPECT_TRUE(is_trie(ta)) << what << ": the tuned run is not a trie (merging on)";
    const uint64_t depth = aa.complete_to_bp;
    auto name = [&](const SeedResult &r, LabelId l) { return r.label_dict.at(l).name; };

    // 1. prefix-subset
    for (const PathResult &p : ta.paths) {
        std::string w = walk_of(ta, p);
        std::set<std::string> labels;
        if (w.size() <= depth) {
            for (const LabelEnd &e : p.end_labels) {
                EXPECT_EQ(0u, e.route_bp) << what;
                labels.insert(name(tuned, e.label));
            }
        } else {
            labels = alive_at(tuned, ta, p, depth);
            w.resize(depth);
        }
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
            EXPECT_TRUE(found) << what << ": tuned leaf " << w << " (" << w.size()
                               << " bp) under " << l << " is not a prefix of an exhaustive "
                                  "walk carrying that label to that depth";
        }
    }

    // 2. every omission has a reason
    const Segment &root = root_of(ta);
    const std::set<std::string> permitted = name_set(tuned, root.labels_start);
    for (const auto &[w, labels] : label_claims(A, a, depth)) {
        for (const auto &[l, a_reason] : labels) {
            if (!permitted.count(l)) {
                bool dropped = false;
                for (const DroppedLabel &d : tuned.dropped_labels)
                    dropped |= d.name == l;
                EXPECT_TRUE(dropped) << what << ": label " << l
                                     << " is neither permitted nor reported dropped";
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
                            if (!is_cut(end->reason) && a_reason) {
                                EXPECT_EQ(to_string(*a_reason), to_string(end->reason))
                                    << what << ": " << l << " ends at " << d << " on " << w
                                    << " for different reasons";
                            }
                            present = true;
                        } else {
                            EXPECT_TRUE(is_cut(end->reason))
                                << what << ": the tuned run ends " << l << " at " << d
                                << " on " << w << " with " << to_string(end->reason)
                                << (end->text.empty() ? "" : " (" + end->text + ")")
                                << " while the exhaustive trie continues it: the two "
                                   "walkers DISAGREE";
                        }
                        decided = true;
                        break;
                    }
                    if (d == w.size()) {
                        // the exhaustive claim was itself cut at the boundary: nothing
                        // is claimed beyond it. Otherwise the tuned run outlives it.
                        EXPECT_FALSE(a_reason.has_value())
                            << what << ": label " << l << " outlives the exhaustive walk "
                            << w << " (" << w.size() << " bp, ended there with "
                            << to_string(a_reason.value_or(EndReason::DEAD_END))
                            << ") in the tuned run";
                        present = true;
                        decided = true;
                        break;
                    }
                    if (j < seq.size()) {
                        if (seq[j] == w[d])
                            continue;   // the tuned run follows w
                        EXPECT_TRUE(branch_recorded(ta, seg->id, d, w[d]))
                            << what << ": the tuned run leaves " << w << " at " << d
                            << " (takes " << seq[j] << ", not " << w[d] << ") with " << l
                            << " alive and no recorded reason";
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
                            ADD_FAILURE() << what << ": label " << l << " is on " << w
                                          << " at " << d << " and the tuned run takes "
                                          << w[d] << ", but the label is not on that child "
                                             "and nothing recorded its end: dropped SILENTLY";
                            decided = true;
                        }
                        break;
                    }
                    if (seg->children.empty()) {
                        const PathResult *p = leaf_path(ta, seg->id);
                        const bool ok = p && p->path_reason && is_cut(*p->path_reason);
                        EXPECT_TRUE(ok) << what << ": tuned leaf " << w.substr(0, d) << " ("
                                        << d << " bp) ends with " << l
                                        << " alive and no recorded reason while the "
                                           "exhaustive trie continues to " << w;
                    } else {
                        EXPECT_TRUE(branch_recorded(ta, seg->id, d, w[d]))
                            << what << ": the tuned run drops " << w[d] << " at " << d
                            << " on " << w << " with " << l
                            << " alive and no recorded reason";
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

// The flank of |label| ending at segment |leaf| in walking order, following at each
// join the parent the label entered through (Segment::labels_via_parent).
inline std::string label_route(const ArmResult &arm, size_t leaf, LabelId label) {
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
            EXPECT_TRUE(found) << "label " << label << " enters segment " << s
                               << " through no parent";
        }
        s = next;
    }
    return std::string(rev.rbegin(), rev.rend());
}

// Merging on: a leaf's spelled walk is not evidence for a label merged in at a join
// (LabelEnd::route_bp > 0) — the label's own route is. Every (route, label) of the
// merged run must be a prefix of an exhaustive leaf under that label. Returns the
// number of (leaf, label) pairs checked.
inline size_t check_routes_subset(const SeedResult &A, const SeedResult &merged,
                                  size_t a, const std::string &what) {
    const ArmResult &ma = merged.arms[a];
    const uint64_t depth = A.arms[a].complete_to_bp;
    const Leaves full = constrained_leaves(A, a, depth);
    size_t checked = 0;
    for (const PathResult &p : ma.paths) {
        const std::string spelled = walk_of(ma, p);
        for (const LabelEnd &e : p.end_labels) {
            std::string route = label_route(ma, p.segments.back(), e.label);
            const std::string &l = merged.label_dict.at(e.label).name;
            // the route and the spelled walk share everything from route_bp on
            EXPECT_EQ(spelled.size(), route.size()) << what;
            EXPECT_EQ(spelled.substr(std::min<size_t>(e.route_bp, spelled.size())),
                      route.substr(std::min<size_t>(e.route_bp, route.size()))) << what;
            if (e.route_bp == 0) {
                EXPECT_EQ(spelled, route) << what;
            } else {
                EXPECT_NE(spelled, route) << what << ": " << l << " is reported merged in at "
                                          << e.route_bp << " but its route is the spelled walk";
            }
            if (route.size() > depth)
                route.resize(depth);
            bool found = false;
            for (auto it = full.lower_bound(route);
                    it != full.end() && it->first.compare(0, route.size(), route) == 0; ++it) {
                if (it->second.count(l)) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << what << ": the route of " << l << " to leaf " << p.id
                               << " (" << route << ") is not a prefix of an exhaustive leaf "
                                  "under that label";
            checked++;
        }
    }
    return checked;
}

} // namespace trie
} // namespace test
} // namespace mtg

#endif // __TEST_TRIE_ORACLE_HPP__
