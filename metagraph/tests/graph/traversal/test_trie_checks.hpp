#ifndef __TEST_TRIE_CHECKS_HPP__
#define __TEST_TRIE_CHECKS_HPP__

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "tests/graph/traversal/test_trie_oracle.hpp"
#include "tests/graph/traversal/test_trie_reference.hpp"

#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"


/*
 * The checks of the trie contract (spec §6.9) that every fixture runs, shared by the
 * unit cases (test_trie_cases.cpp) and the real-index tests (test_mini_refseq.cpp):
 * check_case_on() demands, for one annotated graph and one seed,
 *
 *   1. leaves(T) == the structural walks of the records, with the same labels at every
 *      node and the same path reasons (check_structural);
 *   2. for every permitted set P: claims(A over P) == the oracle E read off T
 *      (test_trie_oracle.hpp) == the per-label claims of the records, with the same
 *      end reasons, and the §6.3 recurrence at budget 0 gives the same leaves
 *      (check_constrained);
 *   3. the run over P is the union of the single-label runs (check_single_label_union);
 *   4. a permitted set derived from the seed is the carriers and runs like the explicit
 *      list (check_derived).
 *
 * check_switching() compares constant-cost runs with the recurrence over the records,
 * check_trace() compares support: trace with the records themselves.
 */
namespace mtg {
namespace test {
namespace trie {

using namespace mtg::graph;
using namespace mtg::graph::traversal;

inline constexpr size_t kLeft = static_cast<size_t>(Arm::LEFT);
inline constexpr size_t kRight = static_cast<size_t>(Arm::RIGHT);

inline std::string mode_name(DeBruijnGraph::Mode mode) {
    return mode == DeBruijnGraph::BASIC ? "basic"
         : mode == DeBruijnGraph::CANONICAL ? "canonical" : "primary";
}

inline Strategy exhaustive(LabelMode mode, uint64_t radius) {
    Strategy st;
    st.exhaustive = true;
    st.label_mode = mode;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.merge_reconverge = false;
    st.max_extension_bp = radius;
    // the oracle must never be cut by a size cap on these fixtures
    st.max_live_paths = 100'000;
    st.max_paths = 100'000;
    return st;
}

inline SeedResult run(const AnnotatedDBG &anno, const std::string &seq,
               const std::vector<std::string> &labels, const Strategy &st,
               const LabelChangeCost &cost = LabelChangeCost::forbid()) {
    LabelOracle oracle(anno);
    Seed seed;
    seed.sequence = seq;
    seed.labels = labels;
    return traverse_seed(oracle, seed, st, cost);
}

inline std::vector<std::string> as_list(const std::set<std::string> &s) {
    return { s.begin(), s.end() };
}

inline std::string show(const std::set<std::string> &s) {
    std::string out = "{";
    for (const auto &x : s) out += " " + x;
    return out + " }";
}

inline std::set<std::string> leaf_walks(const ArmResult &arm) {
    std::set<std::string> out;
    for (const auto &p : arm.paths) out.insert(trie::walk_of(arm, p));
    return out;
}

template <class Map>
inline std::set<std::string> keys(const Map &m) {
    std::set<std::string> out;
    for (const auto &kv : m) out.insert(kv.first);
    return out;
}

inline std::set<std::string> dict_names(const SeedResult &r) {
    std::set<std::string> out;
    for (const auto &ref : r.label_dict) out.insert(ref.name);
    return out;
}

inline size_t count_events(const ArmResult &arm, EventType type) {
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &ev : seg.events) n += ev.type == type;
    }
    return n;
}

inline size_t count_blocked(const ArmResult &arm, EndReason reason) {
    size_t n = 0;
    for (const auto &seg : arm.segments) {
        for (const auto &ev : seg.events)
            n += ev.type == EventType::BLOCKED && ev.reason == reason;
    }
    return n;
}

inline size_t count_label_ends(const ArmResult &arm, EndReason reason) {
    size_t n = 0;
    for (const auto &run : arm.runs) n += run.ended && run.end_reason == reason;
    return n;
}

// walk -> the leaf's path reason (annotate mode: always set)
inline std::map<std::string, EndReason> path_reasons(const ArmResult &arm) {
    std::map<std::string, EndReason> out;
    for (const auto &p : arm.paths) {
        if (p.path_reason)
            out[trie::walk_of(arm, p)] = *p.path_reason;
    }
    return out;
}

// the leaves of a switching run in the reference's shape: walk -> (labels alive at the
// leaf with their losses, the switch events on the way as (at_bp, label entered))
inline trie::SwitchTrie switch_leaves(const SeedResult &r, size_t a) {
    const ArmResult &arm = r.arms[a];
    trie::SwitchTrie out;
    for (const auto &p : arm.paths) {
        trie::SwitchLeaf leaf;
        for (const LabelEnd &e : p.end_labels)
            leaf.state[r.label_dict.at(e.label).name] = e.loss;
        for (size_t s : path_segments(arm, p)) {
            for (const Event &ev : arm.segments[s].events) {
                if (ev.type == EventType::SWITCH)
                    leaf.switches.emplace_back(ev.at_bp, r.label_dict.at(ev.to).name);
            }
        }
        std::sort(leaf.switches.begin(), leaf.switches.end());
        out[trie::walk_of(arm, p)] = std::move(leaf);
    }
    return out;
}

inline std::string show(const trie::SwitchTrie &t) {
    std::ostringstream os;
    for (const auto &[w, leaf] : t) {
        os << "\n  " << (w.empty() ? "<empty>" : w) << " (" << w.size() << " bp) {";
        for (const auto &[l, loss] : leaf.state) os << ' ' << l << ':' << loss;
        os << " } switches [";
        for (const auto &[at, l] : leaf.switches) os << ' ' << at << "->" << l;
        os << " ]";
    }
    return os.str();
}

inline void expect_switch_equal(trie::SwitchTrie expected, trie::SwitchTrie actual, const std::string &what) {
    for (auto &[w, leaf] : expected) std::sort(leaf.switches.begin(), leaf.switches.end());
    for (auto &[w, leaf] : actual) std::sort(leaf.switches.begin(), leaf.switches.end());
    bool same = expected.size() == actual.size();
    for (auto it = expected.begin(); same && it != expected.end(); ++it) {
        auto jt = actual.find(it->first);
        same = jt != actual.end() && jt->second.state == it->second.state
            && jt->second.switches == it->second.switches;
    }
    EXPECT_TRUE(same) << what << ": the §6.3 recurrence over the records and the walker disagree."
                      << "\nrecords:" << show(expected) << "\nwalker:" << show(actual);
}


/********************************* the checks *********************************/

struct CaseSpec {
    std::string name;
    size_t k;
    std::vector<std::string> sequences, labels;
    std::string seed;
    uint64_t radius;
    // the permitted sets to check; empty = the carriers, every single carrier and
    // (for up to five carriers) every pair
    std::vector<std::set<std::string>> permitted;
    bool skip_hairpins;

    CaseSpec(std::string name_, size_t k_, std::vector<std::string> sequences_,
             std::vector<std::string> labels_, std::string seed_, uint64_t radius_ = 200,
             std::vector<std::set<std::string>> permitted_ = {}, bool skip_hairpins_ = true)
          : name(std::move(name_)), k(k_), sequences(std::move(sequences_)),
            labels(std::move(labels_)), seed(std::move(seed_)), radius(radius_),
            permitted(std::move(permitted_)), skip_hairpins(skip_hairpins_) {}
};

inline std::vector<std::set<std::string>> default_permitted(const std::set<std::string> &carriers) {
    std::vector<std::set<std::string>> out { carriers };
    if (carriers.size() > 1) {
        for (const auto &l : carriers) out.push_back({ l });
    }
    if (carriers.size() > 2 && carriers.size() <= 5) {
        for (auto it = carriers.begin(); it != carriers.end(); ++it) {
            for (auto jt = std::next(it); jt != carriers.end(); ++jt)
                out.push_back({ *it, *jt });
        }
    }
    return out;
}

// 1. The structural trie T against the records: the same walks, the same labels at
// every node, the same path reasons. Returns T for the constrained checks.
inline SeedResult check_structural(const AnnotatedDBG &anno, const trie::StringIndex &ix,
                            const CaseSpec &c, const std::string &where) {
    Strategy st = exhaustive(LabelMode::ANNOTATE, c.radius);
    st.max_labels_per_node = Strategy::kUnlimited;
    st.skip_hairpins = c.skip_hairpins;
    SeedResult T = run(anno, c.seed, {}, st);
    for (size_t a : { kLeft, kRight }) {
        const ArmResult &arm = T.arms[a];
        const std::string what = where + " arm " + to_string(arm.arm) + " [structural]";
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << what;
        EXPECT_EQ(c.radius, arm.complete_to_bp) << what;
        EXPECT_EQ(0u, arm.nodes_labels_truncated) << what;
        EXPECT_TRUE(trie::is_trie(arm)) << what;
        const trie::RefTrie ref = trie::structural_trie(ix, c.seed, arm.arm, c.radius, c.skip_hairpins);
        EXPECT_EQ(keys(ref), leaf_walks(arm))
            << what << ": the structural walks of the records and the annotate run differ";
        for (const PathResult &p : arm.paths) {
            const std::string w = trie::walk_of(arm, p);
            auto it = ref.find(w);
            if (it == ref.end())
                continue;
            EXPECT_TRUE(p.path_reason.has_value()) << what << " " << w;
            if (p.path_reason) {
                EXPECT_STREQ(to_string(it->second.reason), to_string(*p.path_reason))
                    << what << ": path " << (w.empty() ? "<empty>" : w) << " ends for another reason";
            }
            bool cut = false;
            const auto recorded = trie::recorded_labels(T, arm, p, &cut);
            EXPECT_FALSE(cut) << what;
            EXPECT_EQ(it->second.labels.size(), recorded.size()) << what << " " << w;
            if (it->second.labels.size() != recorded.size())
                continue;
            for (size_t d = 0; d < recorded.size(); ++d) {
                EXPECT_EQ(it->second.labels[d], recorded[d])
                    << what << ": labels recorded at depth " << d << " of "
                    << (w.empty() ? "<empty>" : w) << " differ from the records";
            }
        }
    }
    for (const std::string &name : dict_names(T))
        EXPECT_TRUE(ix.by_label.count(name)) << where << ": recorded label " << name << " is in no record";
    return T;
}

// 2. The constrained walker over P against the oracle E read off T and against the
// per-label claims of the records; the end reasons; the recurrence at budget 0.
inline SeedResult check_constrained(const AnnotatedDBG &anno, const trie::StringIndex &ix,
                             const CaseSpec &c, const SeedResult &T,
                             const std::set<std::string> &P, const std::string &where) {
    const std::string what = where + " P=" + show(P);
    Strategy st = exhaustive(LabelMode::CONSTRAIN, c.radius);
    st.skip_hairpins = c.skip_hairpins;
    SeedResult A = run(anno, c.seed, as_list(P), st);
    EXPECT_EQ(P.size(), A.num_seed_labels) << what;
    EXPECT_TRUE(A.dropped_labels.empty()) << what;
    for (size_t a : { kLeft, kRight }) {
        const ArmResult &arm = A.arms[a];
        const std::string here = what + " arm " + to_string(arm.arm);
        EXPECT_EQ(ArmResult::COMPLETE, arm.status) << here;
        EXPECT_EQ(c.radius, arm.complete_to_bp) << here;
        trie::Verdict v = trie::verify(T, A, a, P);
        EXPECT_EQ(c.radius, v.depth) << here;
        EXPECT_FALSE(v.cut) << here;
        trie::expect_equal(v.expected, v.actual, here + " [oracle E]");
        const trie::RefClaims ref = trie::claims(ix, P, c.seed, arm.arm, c.radius, c.skip_hairpins);
        trie::expect_equal(ref.leaves, v.actual, here + " [records]");
        for (const auto &[w, labels] : trie::label_claims(A, a, c.radius)) {
            for (const auto &[l, reason] : labels) {
                auto it = ref.reasons.find({ w, l });
                if (!reason || it == ref.reasons.end())
                    continue;
                EXPECT_STREQ(to_string(it->second), to_string(*reason))
                    << here << ": " << l << " ends on " << (w.empty() ? "<empty>" : w)
                    << " (" << w.size() << " bp) for another reason than the records say";
            }
        }
        // the recurrence with nothing to spend is cost forbid
        const trie::SwitchTrie s0 = trie::switch_trie(ix, P, c.seed, arm.arm, c.radius, 1, 0, c.skip_hairpins);
        trie::Leaves from_recurrence;
        for (const auto &[w, leaf] : s0) from_recurrence[w] = keys(leaf.state);
        trie::expect_equal(from_recurrence, trie::constrained_leaves(A, a, c.radius),
                           here + " [recurrence at budget 0]");
    }
    return A;
}

// 3. The multi-label run is the union of the single-label runs: every claim of the
// run over P is a claim of exactly the single-label runs of the labels claiming it.
inline void check_single_label_union(const AnnotatedDBG &anno, const CaseSpec &c, const SeedResult &A,
                              const std::set<std::string> &P, const std::string &where) {
    if (P.size() < 2)
        return;
    Strategy st = exhaustive(LabelMode::CONSTRAIN, c.radius);
    st.skip_hairpins = c.skip_hairpins;
    std::array<trie::Leaves, 2> unions;
    for (const std::string &l : P) {
        SeedResult single = run(anno, c.seed, { l }, st);
        for (size_t a : { kLeft, kRight }) {
            for (const auto &[w, labels] : trie::constrained_claims(single, a, c.radius)) {
                EXPECT_EQ((std::set<std::string>{ l }), labels) << where << " single-label run of " << l;
                unions[a][w].insert(l);
            }
        }
    }
    for (size_t a : { kLeft, kRight }) {
        trie::expect_equal(unions[a], trie::constrained_claims(A, a, c.radius),
                           where + " arm " + to_string(A.arms[a].arm)
                               + " [union of single-label runs vs P=" + show(P) + "]");
    }
}

// 4. A permitted set derived from the seed (seeds[].labels omitted) is the carriers,
// and the run over it equals the run over the explicit list.
inline void check_derived(const AnnotatedDBG &anno, const CaseSpec &c, const SeedResult &A,
                   const std::set<std::string> &carriers, const std::string &where) {
    Strategy st = exhaustive(LabelMode::CONSTRAIN, c.radius);
    st.skip_hairpins = c.skip_hairpins;
    SeedResult D = run(anno, c.seed, {}, st);
    EXPECT_TRUE(D.labels_from_seed) << where;
    EXPECT_EQ(carriers, dict_names(D)) << where << " [derived set]";
    EXPECT_EQ(carriers.size(), D.num_seed_labels) << where;
    for (size_t a : { kLeft, kRight }) {
        trie::expect_equal(trie::constrained_claims(A, a, c.radius),
                           trie::constrained_claims(D, a, c.radius),
                           where + " arm " + to_string(A.arms[a].arm) + " [derived vs explicit]");
    }
}

struct CaseResult {
    SeedResult T;                                  // the structural trie
    std::map<std::set<std::string>, SeedResult> A; // per permitted set
    std::set<std::string> carriers;
};

// phase timings on stderr when TRIE_TIMINGS is set (the real-index tests are slow
// enough to want to know where the time goes)
struct PhaseTimer {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    void lap(const std::string &what) {
        static const bool on = std::getenv("TRIE_TIMINGS") != nullptr;
        if (!on)
            return;
        const auto now = std::chrono::steady_clock::now();
        std::cerr << "[trie] " << what << ": "
                  << std::chrono::duration<double>(now - start).count() << " s\n";
        start = now;
    }
};

inline CaseResult check_case_on(const AnnotatedDBG &anno, const trie::StringIndex &ix,
                         DeBruijnGraph::Mode mode, const CaseSpec &c) {
    const std::string where = c.name + " (" + mode_name(mode) + ", k=" + std::to_string(c.k) + ")";
    PhaseTimer timer;
    CaseResult out;
    out.carriers = ix.carriers(c.seed);
    EXPECT_FALSE(out.carriers.empty()) << where << ": no label carries the seed";
    out.T = check_structural(anno, ix, c, where);
    timer.lap(where + " structural");
    const auto sets = c.permitted.empty() ? default_permitted(out.carriers) : c.permitted;
    for (const auto &P : sets) {
        SeedResult A = check_constrained(anno, ix, c, out.T, P, where);
        timer.lap(where + " constrained " + show(P));
        check_single_label_union(anno, c, A, P, where);
        timer.lap(where + " single-label union " + show(P));
        if (P == out.carriers) {
            check_derived(anno, c, A, out.carriers, where);
            timer.lap(where + " derived");
        }
        out.A.emplace(P, std::move(A));
    }
    return out;
}


// literal expectations: the leaves of the run over P, and a path reason of T
inline void expect_leaves(const CaseResult &r, const std::set<std::string> &P, size_t a,
                   const trie::Leaves &expected, uint64_t radius, const std::string &where) {
    auto it = r.A.find(P);
    ASSERT_NE(r.A.end(), it) << where << ": no run over " << show(P);
    EXPECT_EQ(expected, trie::constrained_leaves(it->second, a, radius)) << where << " P=" << show(P);
}

inline void expect_path_reason(const CaseResult &r, size_t a, const std::string &walk,
                        EndReason reason, const std::string &where) {
    const auto reasons = path_reasons(r.T.arms[a]);
    auto it = reasons.find(walk);
    ASSERT_NE(reasons.end(), it) << where << ": no structural leaf " << walk;
    EXPECT_STREQ(to_string(reason), to_string(it->second)) << where << " leaf " << walk;
}


// forbid, then constant(1) at the given budgets, against the §6.3 recurrence over the
// records (switch_on: loss); the budget-0 run is also the forbid run's leaves
inline void check_switching(const AnnotatedDBG &anno, const trie::StringIndex &ix, const CaseSpec &c,
                     const std::set<std::string> &P, const std::vector<double> &budgets,
                     std::map<double, SeedResult> *runs, const std::string &where) {
    for (double budget : budgets) {
        Strategy st = exhaustive(LabelMode::CONSTRAIN, c.radius);
        st.loss_budget = budget;
        SeedResult A = run(anno, c.seed, as_list(P), st, LabelChangeCost::constant(1));
        for (size_t a : { kLeft, kRight }) {
            EXPECT_EQ(ArmResult::COMPLETE, A.arms[a].status) << where;
            expect_switch_equal(trie::switch_trie(ix, P, c.seed, A.arms[a].arm, c.radius, 1, budget),
                                switch_leaves(A, a),
                                where + " arm " + to_string(A.arms[a].arm) + " budget " + std::to_string(budget));
        }
        runs->emplace(budget, std::move(A));
    }
}


// support: trace against the records (BASIC graphs with k-mer coordinates): every
// label's claims are its records' own continuations, the k-mer run is a superset
inline void check_trace(const AnnotatedDBG &anno, const trie::StringIndex &ix, const CaseSpec &c,
                 const std::set<std::string> &P, const SeedResult &kmer, SeedResult *out,
                 const std::string &where) {
    Strategy st = exhaustive(LabelMode::CONSTRAIN, c.radius);
    st.support = Support::TRACE;
    *out = run(anno, c.seed, as_list(P), st);
    EXPECT_STREQ("tuples", out->access_path) << where;
    for (size_t a : { kLeft, kRight }) {
        const Arm arm = out->arms[a].arm;
        const std::string here = where + " arm " + to_string(arm) + " [trace]";
        EXPECT_EQ(ArmResult::COMPLETE, out->arms[a].status) << here;
        const trie::RefClaims ref = trie::trace_claims(ix, P, c.seed, arm, c.radius);
        const trie::Leaves actual = trie::constrained_claims(*out, a, c.radius);
        trie::expect_equal(ref.leaves, actual, here);
        for (const auto &[w, labels] : trie::label_claims(*out, a, c.radius)) {
            for (const auto &[l, reason] : labels) {
                auto it = ref.reasons.find({ w, l });
                if (!reason || it == ref.reasons.end())
                    continue;
                EXPECT_STREQ(to_string(it->second), to_string(*reason))
                    << here << ": " << l << " ends on " << (w.empty() ? "<empty>" : w) << " for another reason";
            }
        }
        // every trace claim is a prefix of a k-mer claim under the same label (a claim,
        // not a leaf: the label's k-mer walk may end inside a walk other labels go on
        // with, e.g. where its record ends)
        const trie::Leaves full = trie::constrained_claims(kmer, a, c.radius);
        for (const auto &[w, labels] : actual) {
            for (const std::string &l : labels) {
                bool found = false;
                for (const auto &[fw, fl] : full)
                    found |= fw.size() >= w.size() && fw.compare(0, w.size(), w) == 0 && fl.count(l);
                EXPECT_TRUE(found) << here << ": trace walk " << w << " under " << l
                                   << " is not a prefix of a k-mer walk under it";
            }
        }
    }
}


} // namespace trie
} // namespace test
} // namespace mtg

#endif // __TEST_TRIE_CHECKS_HPP__
