#ifndef __TEST_TRIE_REFERENCE_HPP__
#define __TEST_TRIE_REFERENCE_HPP__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <tsl/hopscotch_set.h>

#include "graph/traversal/traversal_types.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"


/*
 * A string-level model of an annotated index and of the walk rule (spec §6.10), built
 * from the records a fixture is made of and nothing else: no graph object, no
 * annotation matrix, no walker code. It is the third party in the verification
 * contract of §6.9 — the structural trie T and the constrained walker A are compared
 * with each other there, and here both are compared with what the records say:
 *
 *   structural_trie()  the walks the rule admits over ALL k-mers, with the labels
 *                      present at every node — what an annotate run must record;
 *   label_trie()       the walks the rule admits over ONE label's k-mers — a label's
 *                      maximal label-consistent walks, i.e. its claims under cost forbid;
 *   switch_trie()      the §6.3 loss recurrence over a permitted set with a constant
 *                      switch cost and a budget (switch_on: loss) — what a constrain
 *                      run with 0, 1, 2 ... switches must yield;
 *   trace_trie()       the records' own continuations (support: trace, BASIC graphs).
 *
 * Orientation: a k-mer and its reverse complement are one node in canonical and
 * primary graphs, so there both orientations of every record k-mer are in the model;
 * in a basic graph only the record's orientation is. Walks are spelled OUTWARD from
 * the seed boundary (the left arm's bases in order of distance from the seed), the
 * convention of test_trie_oracle.hpp.
 */
namespace mtg {
namespace test {
namespace trie {

using mtg::graph::traversal::Arm;
using mtg::graph::traversal::EndReason;

inline std::string rc_of(std::string s) { ::reverse_complement(s); return s; }

struct StringIndex {
    struct Record {
        std::string label;      // "" = in the graph, carried by no label
        std::string seq;
        uint64_t coord_start;   // coordinate of the record's first k-mer (trace model)
    };

    // nodes are kept 2-bit packed (k <= 31 on ACGT) so that a real index of a few
    // million k-mers fits the model too
    using Kmer = uint64_t;
    static constexpr Kmer kNone = std::numeric_limits<Kmer>::max();
    // packed k-mers are far from uniform in their low bits (std::hash is the identity),
    // which an open-addressing table with power-of-two buckets must not see
    struct KmerHash {
        size_t operator()(Kmer x) const {
            x ^= x >> 33;
            x *= 0xff51afd7ed558ccdULL;
            x ^= x >> 33;
            x *= 0xc4ceb9fe1a85ec53ULL;
            x ^= x >> 33;
            return static_cast<size_t>(x);
        }
    };
    using KmerSet = tsl::hopscotch_set<Kmer, KmerHash>;

    size_t k;
    bool canonical;
    std::vector<Record> records;
    KmerSet unlabeled;                            // nodes of records with label ""
    std::map<std::string, KmerSet> by_label;     // label -> its nodes (both orientations if canonical)

    StringIndex(size_t k_, mtg::graph::DeBruijnGraph::Mode mode,
                const std::vector<std::string> &seqs,
                const std::vector<std::string> &labels,
                const std::vector<uint64_t> &coord_starts = {})
          : k(k_), canonical(mode != mtg::graph::DeBruijnGraph::BASIC) {
        assert(k >= 2 && k <= 31);
        for (size_t i = 0; i < seqs.size(); ++i) {
            records.push_back({ labels[i], seqs[i],
                                coord_starts.empty() ? 0 : coord_starts[i] });
            for (size_t j = 0; j + k <= seqs[i].size(); ++j) {
                add(labels[i], seqs[i].substr(j, k));
            }
        }
    }

    static Kmer encode(const std::string &kmer) {
        Kmer code = 0;
        for (char c : kmer) {
            const int v = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : -1;
            if (v < 0)
                return kNone;
            code = (code << 2) | static_cast<Kmer>(v);
        }
        return code;
    }

    void add(const std::string &label, const std::string &kmer) {
        auto &set = label.empty() ? unlabeled : by_label[label];
        set.insert(encode(kmer));
        if (canonical)
            set.insert(encode(rc_of(kmer)));
    }

    bool has(const std::string &kmer) const {
        const Kmer code = encode(kmer);
        if (code == kNone)
            return false;
        if (unlabeled.count(code))
            return true;
        for (const auto &[l, kmers] : by_label) {
            if (kmers.count(code))
                return true;
        }
        return false;
    }

    bool carries(const std::string &label, const std::string &kmer) const {
        auto it = by_label.find(label);
        const Kmer code = encode(kmer);
        return code != kNone && it != by_label.end() && it->second.count(code);
    }

    std::set<std::string> labels_at(const std::string &kmer) const {
        std::set<std::string> out;
        const Kmer code = encode(kmer);
        if (code == kNone)
            return out;
        for (const auto &[l, kmers] : by_label) {
            if (kmers.count(code))
                out.insert(l);
        }
        return out;
    }

    // the labels carrying every k-mer of |seed| (the walker's derived permitted set)
    std::set<std::string> carriers(const std::string &seed) const {
        std::set<std::string> out;
        for (const auto &[l, kmers] : by_label) {
            bool all = seed.size() >= k;
            for (size_t j = 0; all && j + k <= seed.size(); ++j) {
                const Kmer code = encode(seed.substr(j, k));
                all = code != kNone && kmers.count(code);
            }
            if (all)
                out.insert(l);
        }
        return out;
    }

    // the seed's nodes (both orientations if canonical): re-entering one blocks a step
    std::set<std::string> seed_nodes(const std::string &seed) const {
        std::set<std::string> out;
        for (size_t j = 0; j + k <= seed.size(); ++j) {
            const std::string km = seed.substr(j, k);
            out.insert(km);
            if (canonical)
                out.insert(rc_of(km));
        }
        return out;
    }

    std::string boundary(const std::string &seed, Arm arm) const {
        return arm == Arm::RIGHT ? seed.substr(seed.size() - k) : seed.substr(0, k);
    }

    // one step outward from node |u| by base |c|: the node entered and the (k+1)-mer
    // the step spells
    void step(Arm arm, const std::string &u, char c, std::string *v, std::string *edge) const {
        if (arm == Arm::RIGHT) {
            *v = u.substr(1) + c;
            *edge = u + c;
        } else {
            *v = c + u.substr(0, k - 1);
            *edge = c + u;
        }
    }

    // §6.5: a step spelling a self-reverse-complementary (k+1)-mer; for even k also a
    // step into or out of a self-reverse-complementary node
    bool is_hairpin(const std::string &u, const std::string &v, const std::string &edge) const {
        if (!canonical)
            return false;
        if (rc_of(edge) == edge)
            return true;
        if (k % 2 == 0 && (rc_of(v) == v || rc_of(u) == u))
            return true;
        return false;
    }

    // §6.6: the edge identity (canonical (k+1)-mer in canonical regimes) and whether the
    // step uses it in the flipped orientation
    std::string edge_key(const std::string &edge, bool *flipped) const {
        *flipped = false;
        if (!canonical)
            return edge;
        const std::string r = rc_of(edge);
        if (r < edge) {
            *flipped = true;
            return r;
        }
        return edge;
    }
};


/******************************* the walk rule *******************************/

inline uint8_t block_rank(EndReason r) {
    switch (r) {
        case EndReason::REACHED_SEED: return 3;
        case EndReason::EDGE_REUSE_RC: return 2;
        case EndReason::EDGE_REUSE: return 1;
        default: return 0;
    }
}

inline EndReason block_of_rank(uint8_t rank) {
    return rank == 3 ? EndReason::REACHED_SEED
         : rank == 2 ? EndReason::EDGE_REUSE_RC : EndReason::EDGE_REUSE;
}

// a structural successor of the current node and its admissibility under the rule
struct Step {
    char ch;
    std::string v;        // the node entered
    std::string key;      // edge identity
    bool flipped;
    bool hairpin;         // skipped (hairpins: skip)
    bool blocked;
    EndReason block = EndReason::DEAD_END;
    bool admissible() const { return !hairpin && !blocked; }
};

struct WalkRule {
    const StringIndex &ix;
    std::set<std::string> seed;    // the seed's nodes
    Arm arm;
    bool skip_hairpins = true;

    WalkRule(const StringIndex &ix_, const std::string &seed_seq, Arm arm_)
          : ix(ix_), seed(ix_.seed_nodes(seed_seq)), arm(arm_) {}

    // every node reachable by one step, in base order, with the rule applied against
    // the edges |used| so far on this path (key -> orientation of the use)
    std::vector<Step> steps(const std::string &u, const std::map<std::string, bool> &used) const {
        std::vector<Step> out;
        for (char c : std::string("ACGT")) {
            Step s;
            s.ch = c;
            std::string edge;
            ix.step(arm, u, c, &s.v, &edge);
            if (!ix.has(s.v))
                continue;
            s.hairpin = skip_hairpins && ix.is_hairpin(u, s.v, edge);
            s.blocked = false;
            s.key = ix.edge_key(edge, &s.flipped);
            if (!s.hairpin) {
                if (seed.count(s.v)) {
                    s.blocked = true;
                    s.block = EndReason::REACHED_SEED;
                } else if (auto it = used.find(s.key); it != used.end()) {
                    s.blocked = true;
                    s.block = it->second != s.flipped ? EndReason::EDGE_REUSE_RC
                                                      : EndReason::EDGE_REUSE;
                }
            }
            out.push_back(std::move(s));
        }
        return out;
    }
};


/**************************** structural / per-label ****************************/

struct RefLeaf {
    EndReason reason = EndReason::DEAD_END;
    // the labels present at every node of the walk, index 0 = the seed boundary
    std::vector<std::set<std::string>> labels;
};
using RefTrie = std::map<std::string, RefLeaf>;   // walk (outward) -> leaf

namespace detail {

// The end reason of a walk none of whose candidate steps is followed. |passing| are
// the steps whose node passes the filter (every step for the structural trie):
// none -> the label is on no successor (label_lost); some blocked -> the strongest
// block; only skipped hairpins -> dead_end (the walker's "hairpin" text).
inline EndReason stop_reason(const std::vector<Step> &steps, const std::vector<bool> &passing) {
    bool any_pass = false;
    uint8_t rank = 0;
    for (size_t i = 0; i < steps.size(); ++i) {
        if (!passing[i])
            continue;
        any_pass = true;
        if (steps[i].blocked)
            rank = std::max(rank, block_rank(steps[i].block));
    }
    if (!any_pass)
        return EndReason::LABEL_LOST;
    return rank ? block_of_rank(rank) : EndReason::DEAD_END;
}

template <class Filter>
void dfs(const WalkRule &rule, const std::string &u, uint64_t radius,
         std::map<std::string, bool> &used, std::string &walk,
         std::vector<std::set<std::string>> &at, const Filter &passes, RefTrie *out) {
    if (walk.size() >= radius) {
        (*out)[walk] = RefLeaf{ EndReason::MAX_EXTENSION, at };
        return;
    }
    const std::vector<Step> steps = rule.steps(u, used);
    if (steps.empty()) {
        (*out)[walk] = RefLeaf{ EndReason::DEAD_END, at };
        return;
    }
    std::vector<bool> passing(steps.size());
    bool followed = false;
    for (size_t i = 0; i < steps.size(); ++i) {
        const Step &s = steps[i];
        passing[i] = passes(s.v);
        if (!passing[i] || !s.admissible())
            continue;
        followed = true;
        used[s.key] = s.flipped;
        walk.push_back(s.ch);
        at.push_back(rule.ix.labels_at(s.v));
        dfs(rule, s.v, radius, used, walk, at, passes, out);
        at.pop_back();
        walk.pop_back();
        used.erase(s.key);
    }
    if (!followed)
        (*out)[walk] = RefLeaf{ stop_reason(steps, passing), at };
}

} // namespace detail

// The structural trie from the seed boundary: every walk the rule admits over all
// nodes, with the labels present at each node. A leaf is a walk no step extends.
inline RefTrie structural_trie(const StringIndex &ix, const std::string &seed, Arm arm,
                               uint64_t radius, bool skip_hairpins = true) {
    WalkRule rule(ix, seed, arm);
    rule.skip_hairpins = skip_hairpins;
    RefTrie out;
    std::map<std::string, bool> used;
    std::string walk;
    const std::string b = ix.boundary(seed, arm);
    std::vector<std::set<std::string>> at { ix.labels_at(b) };
    detail::dfs(rule, b, radius, used, walk, at, [](const std::string &) { return true; }, &out);
    return out;
}

// One label's maximal label-consistent walks (its claims under cost forbid), with
// the reason each ends. Empty when the label does not carry the seed boundary.
inline RefTrie label_trie(const StringIndex &ix, const std::string &label,
                          const std::string &seed, Arm arm, uint64_t radius,
                          bool skip_hairpins = true) {
    RefTrie out;
    const std::string b = ix.boundary(seed, arm);
    if (!ix.carries(label, b))
        return out;
    WalkRule rule(ix, seed, arm);
    rule.skip_hairpins = skip_hairpins;
    std::map<std::string, bool> used;
    std::string walk;
    std::vector<std::set<std::string>> at { ix.labels_at(b) };
    detail::dfs(rule, b, radius, used, walk, at,
                [&](const std::string &v) { return ix.carries(label, v); }, &out);
    return out;
}

// The claims of a permitted set under cost forbid: walk -> labels claiming exactly it,
// and the reason of every (walk, label).
struct RefClaims {
    std::map<std::string, std::set<std::string>> leaves;
    std::map<std::pair<std::string, std::string>, EndReason> reasons;   // (walk, label)
};

inline RefClaims claims(const StringIndex &ix, const std::set<std::string> &permitted,
                        const std::string &seed, Arm arm, uint64_t radius,
                        bool skip_hairpins = true) {
    RefClaims out;
    for (const std::string &l : permitted) {
        for (const auto &[w, leaf] : label_trie(ix, l, seed, arm, radius, skip_hairpins)) {
            out.leaves[w].insert(l);
            out.reasons[{ w, l }] = leaf.reason;
        }
    }
    return out;
}


/**************************** the §6.3 loss recurrence ****************************/

struct SwitchLeaf {
    std::map<std::string, double> state;                    // label -> loss at the leaf
    std::vector<std::pair<uint64_t, std::string>> switches; // (at_bp, label entered)
};
using SwitchTrie = std::map<std::string, SwitchLeaf>;

namespace detail {

// |state| the (label, loss) entries at node |u|. One step: targets are the permitted
// labels present at the node entered; the eligible sources are the entries lost there
// (switch_on: loss); a target stays at its own loss when it has an entry, otherwise
// (or when cheaper) it enters from the cheapest source at source loss + cost, within
// the budget. The step is followed iff the derived state is non-empty.
inline void switch_dfs(const WalkRule &rule, const std::set<std::string> &permitted,
                       const std::string &u, uint64_t radius, double cost, double budget,
                       std::map<std::string, bool> &used, std::string &walk,
                       const std::map<std::string, double> &state,
                       std::vector<std::pair<uint64_t, std::string>> &switches,
                       SwitchTrie *out) {
    if (walk.size() >= radius) {
        (*out)[walk] = SwitchLeaf{ state, switches };
        return;
    }
    bool followed = false;
    for (const Step &s : rule.steps(u, used)) {
        if (!s.admissible())
            continue;
        std::set<std::string> targets;
        for (const std::string &l : rule.ix.labels_at(s.v)) {
            if (permitted.count(l))
                targets.insert(l);
        }
        double cheapest_source = std::numeric_limits<double>::infinity();
        for (const auto &[l, loss] : state) {
            if (!targets.count(l))
                cheapest_source = std::min(cheapest_source, loss);
        }
        std::map<std::string, double> next;
        std::vector<std::string> entered;
        for (const std::string &t : targets) {
            auto it = state.find(t);
            const double stay = it != state.end() ? it->second
                                                  : std::numeric_limits<double>::infinity();
            const double sw = cheapest_source + cost;
            if (it != state.end() && stay <= sw) {
                next[t] = stay;
            } else if (sw <= budget) {
                next[t] = sw;
                entered.push_back(t);
            }
        }
        if (next.empty())
            continue;
        followed = true;
        used[s.key] = s.flipped;
        walk.push_back(s.ch);
        const size_t n_sw = switches.size();
        for (const std::string &t : entered)
            switches.emplace_back(walk.size() - 1, t);
        switch_dfs(rule, permitted, s.v, radius, cost, budget, used, walk, next, switches, out);
        switches.resize(n_sw);
        walk.pop_back();
        used.erase(s.key);
    }
    if (!followed)
        (*out)[walk] = SwitchLeaf{ state, switches };
}

} // namespace detail

// The leaves of a constrain run over |permitted| with LabelChangeCost::constant(cost),
// loss_budget |budget| and switch_on: loss: walk -> (state at the leaf, switches on
// the way). Budget 0 (or an infinite cost) is cost forbid: the leaves then carry the
// labels alive at them.
inline SwitchTrie switch_trie(const StringIndex &ix, const std::set<std::string> &permitted,
                              const std::string &seed, Arm arm, uint64_t radius,
                              double cost, double budget, bool skip_hairpins = true) {
    WalkRule rule(ix, seed, arm);
    rule.skip_hairpins = skip_hairpins;
    SwitchTrie out;
    std::map<std::string, bool> used;
    std::string walk;
    std::map<std::string, double> state;
    const std::string b = ix.boundary(seed, arm);
    for (const std::string &l : permitted) {
        if (ix.carries(l, b))
            state[l] = 0;
    }
    std::vector<std::pair<uint64_t, std::string>> switches;
    detail::switch_dfs(rule, permitted, b, radius, cost, budget, used, walk, state, switches, &out);
    return out;
}


/********************************* trace support *********************************/

// support: trace on a BASIC graph: a label's walks are its records' own continuations
// from every occurrence of the seed, subject to the structural rule (seed re-entry and
// edge reuse still block), maximal per label. The reason at a record's end: dead_end
// when the node has no successor at all, record_end when a successor carries the label
// (its coordinates do not continue), label_lost otherwise.
//
// Several occurrences of the seed (in one record or in several records of the label)
// may spell the SAME walk and stop there for different reasons: one is blocked where it
// would re-enter the seed, another reaches its record's end. The walker keeps ONE
// lineage per label with the live coordinates of every occurrence, and a blocked
// successor whose coordinates continue sets the block reason, which outranks a trace
// break — so the reason of a walk is derived from the AGGREGATE of its occurrences
// with the walker's precedence: any structural block wins (the highest block rank:
// rejoined_seed > edge_reuse_rc > edge_reuse), otherwise the occurrences agree. (Taking
// whichever occurrence came last would make the reason depend on the order of the
// records.)
inline RefTrie trace_trie(const StringIndex &ix, const std::string &label,
                          const std::string &seed, Arm arm, uint64_t radius) {
    struct Outcomes {
        std::vector<EndReason> reasons;              // one per occurrence spelling the walk
        std::vector<std::set<std::string>> labels;   // the same nodes for all of them
    };
    std::map<std::string, Outcomes> seen;
    WalkRule rule(ix, seed, arm);
    for (const StringIndex::Record &r : ix.records) {
        if (r.label != label)
            continue;
        for (size_t p = r.seq.find(seed); p != std::string::npos; p = r.seq.find(seed, p + 1)) {
            std::map<std::string, bool> used;
            std::string walk, u = ix.boundary(seed, arm);
            std::vector<std::set<std::string>> at { ix.labels_at(u) };
            EndReason reason = EndReason::DEAD_END;
            while (true) {
                if (walk.size() >= radius) {
                    reason = EndReason::MAX_EXTENSION;
                    break;
                }
                const std::vector<Step> steps = rule.steps(u, used);
                const int64_t idx = arm == Arm::RIGHT
                        ? static_cast<int64_t>(p + seed.size() + walk.size())
                        : static_cast<int64_t>(p) - 1 - static_cast<int64_t>(walk.size());
                if (idx < 0 || idx >= static_cast<int64_t>(r.seq.size())) {
                    // the record ends here
                    bool present = false;
                    for (const Step &s : steps)
                        present |= ix.carries(label, s.v);
                    reason = steps.empty() ? EndReason::DEAD_END
                           : present ? EndReason::RECORD_END : EndReason::LABEL_LOST;
                    break;
                }
                const char c = r.seq[idx];
                const Step *next = nullptr;
                for (const Step &s : steps) {
                    if (s.ch == c)
                        next = &s;
                }
                if (!next || next->blocked) {
                    reason = next ? next->block : EndReason::DEAD_END;
                    break;
                }
                used[next->key] = next->flipped;
                walk.push_back(c);
                at.push_back(ix.labels_at(next->v));
                u = next->v;
            }
            Outcomes &o = seen[walk];
            o.reasons.push_back(reason);
            o.labels = std::move(at);
        }
    }
    // one reason per walk from the aggregate of its occurrences (see above)
    RefTrie all;
    for (auto &[walk, o] : seen) {
        uint8_t rank = 0;
        std::set<EndReason> plain;
        for (EndReason r : o.reasons) {
            if (block_rank(r)) {
                rank = std::max(rank, block_rank(r));
            } else {
                plain.insert(r);
            }
        }
        EndReason reason;
        if (rank) {
            reason = block_of_rank(rank);
        } else if (plain.size() == 1) {
            reason = *plain.begin();
        } else {
            // the same node at the same depth with the same edge history: a stop that
            // is not a block is the radius or the record's end, and both depend on
            // nothing but the node — the model would be wrong, not the walker
            std::string what = "trace_trie: occurrences of the seed spelling the walk "
                + (walk.empty() ? std::string("<empty>") : walk) + " end for different reasons:";
            for (EndReason r : plain)
                what += std::string(" ") + mtg::graph::traversal::to_string(r);
            throw std::logic_error(what);
        }
        all[walk] = RefLeaf{ reason, std::move(o.labels) };
    }
    // maximal walks only: an occurrence whose record ends inside another occurrence's
    // walk makes no claim of its own (the label still continues on the other one)
    RefTrie out;
    for (auto it = all.begin(); it != all.end(); ++it) {
        auto next = std::next(it);
        const bool extended = next != all.end() && next->first.size() > it->first.size()
            && next->first.compare(0, it->first.size(), it->first) == 0;
        if (!extended)
            out.insert(*it);
    }
    return out;
}

inline RefClaims trace_claims(const StringIndex &ix, const std::set<std::string> &permitted,
                              const std::string &seed, Arm arm, uint64_t radius) {
    RefClaims out;
    for (const std::string &l : permitted) {
        for (const auto &[w, leaf] : trace_trie(ix, l, seed, arm, radius)) {
            out.leaves[w].insert(l);
            out.reasons[{ w, l }] = leaf.reason;
        }
    }
    return out;
}

} // namespace trie
} // namespace test
} // namespace mtg

#endif // __TEST_TRIE_REFERENCE_HPP__
