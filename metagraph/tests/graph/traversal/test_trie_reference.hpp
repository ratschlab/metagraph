#ifndef __TEST_TRIE_REFERENCE_HPP__
#define __TEST_TRIE_REFERENCE_HPP__

#include <cassert>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <tsl/hopscotch_set.h>

#include "graph/representation/base/sequence_graph.hpp"
#include "graph/traversal/traversal_types.hpp"
#include "common/seq_tools/reverse_complement.hpp"

/**
 * A string-level model of an annotated index and of the walk rule, built from the records
 * a fixture is made of and nothing else: no graph object, no annotation matrix, no walker
 * code. label_trie() gives the walks the rule admits over ONE label's k-mers — the label's
 * maximal label-consistent walks, which is what the walker claims for that label with an
 * unlimited branch allowance.
 *
 * Orientation: a k-mer and its reverse complement are one node in canonical and primary
 * graphs, so there both orientations of every record k-mer are in the model; in a basic
 * graph only the record's orientation is. Walks are spelled OUTWARD from the seed boundary
 * (the left arm's bases in order of distance from the seed).
 */

namespace mtg {
namespace test {
namespace trie {

using mtg::graph::traversal::Arm;
using mtg::graph::traversal::EndReason;

inline std::string rc_of(std::string s) { ::reverse_complement(s); return s; }

struct StringIndex {
    // 2-bit packed k-mers (k <= 31 on ACGT)
    using Kmer = uint64_t;
    static constexpr Kmer kNone = std::numeric_limits<Kmer>::max();
    // packed k-mers are far from uniform in their low bits
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
    std::map<std::string, KmerSet> by_label;     // label -> its nodes (both orientations if canonical)

    StringIndex(size_t k_, mtg::graph::DeBruijnGraph::Mode mode,
                const std::vector<std::string> &seqs,
                const std::vector<std::string> &labels)
          : k(k_), canonical(mode != mtg::graph::DeBruijnGraph::BASIC) {
        assert(k >= 2 && k <= 31);
        for (size_t i = 0; i < seqs.size(); ++i) {
            for (size_t j = 0; j + k <= seqs[i].size(); ++j) {
                const std::string kmer = seqs[i].substr(j, k);
                by_label[labels[i]].insert(encode(kmer));
                if (canonical)
                    by_label[labels[i]].insert(encode(rc_of(kmer)));
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

    bool has(const std::string &kmer) const {
        const Kmer code = encode(kmer);
        for (const auto &[l, kmers] : by_label) {
            if (code != kNone && kmers.count(code))
                return true;
        }
        return false;
    }

    bool carries(const std::string &label, const std::string &kmer) const {
        auto it = by_label.find(label);
        const Kmer code = encode(kmer);
        return code != kNone && it != by_label.end() && it->second.count(code);
    }

    // the labels carrying every k-mer of |seed| (the walker's derived permitted set)
    std::set<std::string> carriers(const std::string &seed) const {
        std::set<std::string> out;
        for (const auto &[l, kmers] : by_label) {
            bool all = seed.size() >= k;
            for (size_t j = 0; all && j + k <= seed.size(); ++j) {
                all = carries(l, seed.substr(j, k));
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

    // one step outward from node |u| by base |c|: the node entered and the (k+1)-mer the
    // step spells
    void step(Arm arm, const std::string &u, char c, std::string *v, std::string *edge) const {
        if (arm == Arm::RIGHT) {
            *v = u.substr(1) + c;
            *edge = u + c;
        } else {
            *v = c + u.substr(0, k - 1);
            *edge = c + u;
        }
    }

    // a step spelling a self-reverse-complementary (k+1)-mer; for even k also a step into
    // or out of a self-reverse-complementary node
    bool is_hairpin(const std::string &u, const std::string &v, const std::string &edge) const {
        if (!canonical)
            return false;
        if (rc_of(edge) == edge)
            return true;
        return k % 2 == 0 && (rc_of(v) == v || rc_of(u) == u);
    }

    // the edge identity (canonical (k+1)-mer in canonical regimes) and whether the step
    // uses it in the flipped orientation
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


inline uint8_t block_rank(EndReason r) {
    switch (r) {
        case EndReason::REACHED_SEED: return 3;
        case EndReason::EDGE_REUSE_RC: return 2;
        case EndReason::EDGE_REUSE: return 1;
        default: return 0;
    }
}

// a structural successor of the current node and its admissibility under the rule
struct Step {
    char ch;
    std::string v;        // the node entered
    std::string key;      // edge identity
    bool flipped;
    bool hairpin;         // skipped
    bool blocked;
    EndReason block = EndReason::DEAD_END;
    bool admissible() const { return !hairpin && !blocked; }
};

struct WalkRule {
    const StringIndex &ix;
    std::set<std::string> seed;    // the seed's nodes
    Arm arm;

    WalkRule(const StringIndex &ix_, const std::string &seed_seq, Arm arm_)
          : ix(ix_), seed(ix_.seed_nodes(seed_seq)), arm(arm_) {}

    // every node reachable by one step, in base order, with the rule applied against the
    // edges |used| so far on this path (key -> orientation of the use)
    std::vector<Step> steps(const std::string &u, const std::map<std::string, bool> &used) const {
        std::vector<Step> out;
        for (char c : std::string("ACGT")) {
            Step s;
            s.ch = c;
            std::string edge;
            ix.step(arm, u, c, &s.v, &edge);
            if (!ix.has(s.v))
                continue;
            s.hairpin = ix.is_hairpin(u, s.v, edge);
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

// walk (outward) -> the reason it ends
using RefTrie = std::map<std::string, EndReason>;

namespace detail {

// The end reason of a walk none of whose candidate steps is followed. |passing| are the
// steps whose node carries the label: none -> label_lost; some blocked -> the strongest
// block; only skipped hairpins -> dead_end.
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
    return rank == 3 ? EndReason::REACHED_SEED
         : rank == 2 ? EndReason::EDGE_REUSE_RC
         : rank == 1 ? EndReason::EDGE_REUSE : EndReason::DEAD_END;
}

inline void dfs(const WalkRule &rule, const std::string &label, const std::string &u,
                uint64_t radius, std::map<std::string, bool> &used, std::string &walk,
                RefTrie *out) {
    if (walk.size() >= radius) {
        (*out)[walk] = EndReason::MAX_EXTENSION;
        return;
    }
    const std::vector<Step> steps = rule.steps(u, used);
    if (steps.empty()) {
        (*out)[walk] = EndReason::DEAD_END;
        return;
    }
    std::vector<bool> passing(steps.size());
    bool followed = false;
    for (size_t i = 0; i < steps.size(); ++i) {
        const Step &s = steps[i];
        passing[i] = rule.ix.carries(label, s.v);
        if (!passing[i] || !s.admissible())
            continue;
        followed = true;
        used[s.key] = s.flipped;
        walk.push_back(s.ch);
        dfs(rule, label, s.v, radius, used, walk, out);
        walk.pop_back();
        used.erase(s.key);
    }
    if (!followed)
        (*out)[walk] = stop_reason(steps, passing);
}

} // namespace detail

// One label's maximal label-consistent walks, with the reason each ends. Empty when the
// label does not carry the seed boundary.
inline RefTrie label_trie(const StringIndex &ix, const std::string &label,
                          const std::string &seed, Arm arm, uint64_t radius) {
    RefTrie out;
    const std::string b = ix.boundary(seed, arm);
    if (!ix.carries(label, b))
        return out;
    WalkRule rule(ix, seed, arm);
    std::map<std::string, bool> used;
    std::string walk;
    detail::dfs(rule, label, b, radius, used, walk, &out);
    return out;
}

// The claims of a permitted set: (walk, label) -> the reason the label's walk ends
using RefClaims = std::map<std::pair<std::string, std::string>, EndReason>;

inline RefClaims claims(const StringIndex &ix, const std::set<std::string> &permitted,
                        const std::string &seed, Arm arm, uint64_t radius) {
    RefClaims out;
    for (const std::string &l : permitted) {
        for (const auto &[w, reason] : label_trie(ix, l, seed, arm, radius)) {
            out[{ w, l }] = reason;
        }
    }
    return out;
}

} // namespace trie
} // namespace test
} // namespace mtg

#endif // __TEST_TRIE_REFERENCE_HPP__
