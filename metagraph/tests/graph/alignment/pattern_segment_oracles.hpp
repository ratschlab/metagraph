#ifndef __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__
#define __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__

/**
 * Oracles of the segment search at the graph level (the windows, their anchors and admission,
 * ownership, the roots' extents and the work), independent of the engine:
 *  - a pattern is read through this file's own tables: the IUPAC codes as the bases each
 *    admits, a genetic code as NCBI's 64 residues in TCAG order (the data of
 *    GeneticCode::amino_acids, decoded here), X, B, Z, J and '*' derived here. An oriented
 *    pattern Q is a predicate on strings placed at an offset (OracleQ::admits): for a peptide,
 *    every residue the string overlaps has a codon agreeing with it on the positions inside it;
 *  - the graph is its set of k-mers, spelled (the node ids only name them), its edges the
 *    (k - 1)-overlaps;
 *  - the windows are listed by their own formulation, the anchors found by scanning every
 *    k-mer, the arms computed by breadth-first layers over (k-mer, offset) states, a step
 *    allowed when the (k + 1)-mer it forms is admitted at its offset.
 * They share with the engine only the node ids and the definitions of the contract.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "graph/alignment/genetic_code.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace test {
namespace segments {

using graph::DeBruijnGraph;
using graph::pattern::Orientation;
using graph::pattern::PatternKind;
using graph::pattern::Strands;
using graph::pattern::WindowChoice;
typedef DeBruijnGraph::node_index node_index;

constexpr char kOracleBases[] = "ACGT";

inline std::string oracle_iupac(char code) {
    switch (code) {
        case 'A': return "A";
        case 'C': return "C";
        case 'G': return "G";
        case 'T': return "T";
        case 'R': return "AG";
        case 'Y': return "CT";
        case 'S': return "CG";
        case 'W': return "AT";
        case 'K': return "GT";
        case 'M': return "AC";
        case 'B': return "CGT";
        case 'D': return "AGT";
        case 'H': return "ACT";
        case 'V': return "ACG";
        case 'N': return "ACGT";
        default: return "";
    }
}

inline char oracle_complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return '?';
    }
}

inline std::string oracle_rev_comp(std::string_view s) {
    std::string rc(s.rbegin(), s.rend());
    for (char &c : rc) {
        c = oracle_complement(c);
    }
    return rc;
}

// the 64 codons as strings, in A C G T order
inline const std::vector<std::string>& oracle_codons() {
    static const std::vector<std::string> codons = []() {
        std::vector<std::string> result;
        for (char a : std::string(kOracleBases)) {
            for (char b : std::string(kOracleBases)) {
                for (char c : std::string(kOracleBases)) {
                    result.push_back({ a, b, c });
                }
            }
        }
        return result;
    }();
    return codons;
}

// the residue a codon codes in NCBI table |code_id|: its ncbieaa string read in TCAG order
inline char oracle_translate(int code_id, std::string_view codon) {
    static const std::string order = "TCAG";
    const std::string_view table = graph::pattern::GeneticCode::get(code_id).amino_acids();
    size_t index = 0;
    for (char c : codon) {
        index = 4 * index + order.find(c);
    }
    return table[index];
}

// a peptide letter as the amino acids it admits ('*' the stops)
inline bool oracle_residue_admits(char residue, char amino) {
    switch (residue) {
        case 'X': return amino != '*';
        case 'B': return amino == 'D' || amino == 'N';
        case 'Z': return amino == 'E' || amino == 'Q';
        case 'J': return amino == 'I' || amino == 'L';
        default: return amino == residue;
    }
}

// the codon strings of |residue| in table |code_id|
inline std::vector<std::string> oracle_residue_codons(char residue, int code_id) {
    std::vector<std::string> result;
    for (const std::string &codon : oracle_codons()) {
        if (oracle_residue_admits(residue, oracle_translate(code_id, codon)))
            result.push_back(codon);
    }
    return result;
}

/**
 * One oriented pattern Q as a predicate on strings placed at an offset. DNA and IUPAC: the
 * bases of each position. PROTEIN: per residue of Q (for rc(P): P's residues reversed), the
 * codon strings it admits in Q's reading direction (for rc(P) each reverse-complemented).
 */
struct OracleQ {
    Orientation orientation = Orientation::FORWARD;
    bool protein = false;
    size_t length = 0;
    std::vector<std::string> sets;
    std::vector<std::vector<std::string>> codons;

    // |x| placed at Q-offset |t| is admitted (by the cylinder of each residue it overlaps)
    bool admits(std::string_view x, size_t t) const {
        if (t + x.size() > length)
            return false;
        if (!protein) {
            for (size_t i = 0; i < x.size(); ++i) {
                if (sets[t + i].find(x[i]) == std::string::npos)
                    return false;
            }
            return true;
        }
        for (size_t r = t / 3; 3 * r < t + x.size(); ++r) {
            bool any = false;
            for (const std::string &c : codons[r]) {
                bool ok = true;
                for (size_t q = 0; q < 3 && ok; ++q) {
                    const size_t p = 3 * r + q;
                    if (p >= t && p < t + x.size())
                        ok = c[q] == x[p - t];
                }
                if (ok) {
                    any = true;
                    break;
                }
            }
            if (!any)
                return false;
        }
        return true;
    }

    // the information bits of Q[s, s + k): 2k less log2 of the distinct strings each position
    // (a peptide: each residue, on its positions inside the window) admits
    double bits(size_t s, size_t k) const {
        double bits = 2.0 * k;
        if (!protein) {
            for (size_t i = s; i < s + k; ++i) {
                bits -= std::log2(static_cast<double>(sets[i].size()));
            }
            return bits;
        }
        for (size_t r = s / 3; 3 * r < s + k; ++r) {
            std::set<std::string> projections;
            for (const std::string &c : codons[r]) {
                std::string piece;
                for (size_t q = 0; q < 3; ++q) {
                    if (3 * r + q >= s && 3 * r + q < s + k)
                        piece.push_back(c[q]);
                }
                projections.insert(piece);
            }
            bits -= std::log2(static_cast<double>(std::max<size_t>(1, projections.size())));
        }
        return bits;
    }
};

// the oriented patterns a request searches, in plan order (a pattern equal to its reverse
// complement once, as PALINDROMIC)
inline std::vector<OracleQ> oracle_orientations(PatternKind kind, const std::string &text,
                                                Strands strands, int code_id = 1) {
    OracleQ forward;
    OracleQ reverse;
    forward.orientation = Orientation::FORWARD;
    reverse.orientation = Orientation::REVERSE;
    if (kind == PatternKind::PROTEIN) {
        forward.protein = reverse.protein = true;
        forward.length = reverse.length = 3 * text.size();
        for (char residue : text) {
            forward.codons.push_back(oracle_residue_codons(residue, code_id));
        }
        for (auto it = forward.codons.rbegin(); it != forward.codons.rend(); ++it) {
            std::vector<std::string> rc;
            for (const std::string &c : *it) {
                rc.push_back(oracle_rev_comp(c));
            }
            reverse.codons.push_back(rc);
        }
    } else {
        forward.length = reverse.length = text.size();
        for (char c : text) {
            forward.sets.push_back(oracle_iupac(c));
        }
        for (auto it = forward.sets.rbegin(); it != forward.sets.rend(); ++it) {
            std::string set;
            for (char b : *it) {
                set.push_back(oracle_complement(b));
            }
            std::sort(set.begin(), set.end());
            reverse.sets.push_back(set);
        }
    }
    // palindromic: the two predicates are the same, position by position (residue by residue)
    bool palindromic = forward.sets == reverse.sets;
    if (kind == PatternKind::PROTEIN) {
        palindromic = true;
        for (size_t r = 0; r < forward.codons.size(); ++r) {
            std::vector<std::string> a = forward.codons[r];
            std::vector<std::string> b = reverse.codons[r];
            std::sort(a.begin(), a.end());
            std::sort(b.begin(), b.end());
            palindromic &= a == b;
        }
    }
    if (palindromic) {
        forward.orientation = Orientation::PALINDROMIC;
        return { forward };
    }
    std::vector<OracleQ> result;
    if (strands != Strands::REVERSE)
        result.push_back(forward);
    if (strands != Strands::FORWARD)
        result.push_back(reverse);
    return result;
}

// the served graph as its k-mers (spelled) with their node ids: the valid edges of the
// DBGSuccinct, and on a wrapped PRIMARY graph also the wrapper's reverse complements
struct OracleGraph {
    size_t k = 0;
    std::map<std::string, node_index> kmers;

    bool has(const std::string &kmer) const { return kmers.count(kmer); }

    // the k-mers following |u| (by overlap), whatever their base
    std::vector<std::string> outgoing(const std::string &u) const {
        std::vector<std::string> result;
        for (char c : std::string(kOracleBases)) {
            const std::string v = u.substr(1) + c;
            if (has(v))
                result.push_back(v);
        }
        return result;
    }

    std::vector<std::string> incoming(const std::string &u) const {
        std::vector<std::string> result;
        for (char c : std::string(kOracleBases)) {
            const std::string w = c + u.substr(0, u.size() - 1);
            if (has(w))
                result.push_back(w);
        }
        return result;
    }
};

inline const graph::DBGSuccinct& oracle_base_dbg(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const graph::CanonicalDBG*>(&graph))
        return dynamic_cast<const graph::DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const graph::DBGSuccinct&>(graph);
}

inline OracleGraph oracle_graph(const DeBruijnGraph &graph) {
    OracleGraph g;
    g.k = graph.get_k();
    const graph::DBGSuccinct &dbg_succ = oracle_base_dbg(graph);
    const auto *canonical = dynamic_cast<const graph::CanonicalDBG*>(&graph);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (!dbg_succ.in_graph(y))
            continue;
        const std::string kmer = graph.get_node_sequence(y);
        if (kmer.find('$') != std::string::npos)
            continue;
        g.kmers.emplace(kmer, y);
        if (canonical) {
            const node_index z = canonical->reverse_complement(y);
            if (z != y)
                g.kmers.emplace(graph.get_node_sequence(z), z);
        }
    }
    return g;
}

// the window starts of one orientation, by their own formulation: {min(i S, L - W')}, ascending
inline std::vector<uint32_t> oracle_window_grid(uint32_t L, uint32_t w, uint32_t step,
                                                WindowChoice choice) {
    std::set<uint32_t> grid;
    const uint32_t last = L - w;
    switch (choice) {
        case WindowChoice::FIRST:
            grid.insert(0);
            break;
        case WindowChoice::ENDS:
            grid.insert(0);
            grid.insert(last);
            break;
        case WindowChoice::STEPS:
            for (uint64_t i = 0;; ++i) {
                const uint64_t s = std::min<uint64_t>(i * step, last);
                grid.insert(static_cast<uint32_t>(s));
                if (s == last)
                    break;
            }
            break;
    }
    return std::vector<uint32_t>(grid.begin(), grid.end());
}

struct OraclePlan {
    // W (0: k), S, the choice
    uint32_t window = 0;
    uint32_t step = 15;
    WindowChoice choice = WindowChoice::STEPS;
    uint64_t max_anchors = 1000;
    double min_information_bits = 0;
};

struct OracleWindow {
    uint32_t offset = 0;
    double bits = 0;
    // the k-mers admitted at the offset, in node order
    std::vector<std::pair<node_index, std::string>> anchors;
    bool below_floor = false;
    bool admitted = false;
    uint64_t roots = 0;
    uint64_t owned = 0;
};

// one root's extent, in the search's order
struct OracleExtent {
    Orientation orientation;
    uint32_t window;
    node_index root;
    uint32_t from;
    uint32_t to;
    bool reported;
};

typedef std::tuple<Orientation, uint32_t, uint32_t> OracleSegmentKey;

struct OracleRun {
    uint32_t window = 0;
    std::map<Orientation, std::vector<OracleWindow>> windows;
    std::vector<OracleExtent> extents;
    // reported extents grouped by (orientation, P's from, P's to): their roots
    std::map<OracleSegmentKey, uint64_t> segments;
    // the work of the memoised arms: states memoised (expanded or reached beside the roots),
    // the edges examined, the branchings of the right arm, the horizon cuts
    uint64_t states = 0;
    uint64_t right_edges = 0;
    uint64_t right_branches = 0;
    uint64_t left_edges = 0;
    uint64_t left_cuts = 0;
    // the k-mers the left arm expanded (once per window): the bwd cost of their incoming edges
    std::vector<std::string> left_expanded;
    // the roots at an offset > 0 (each pays the bwd steps of its own pair on a DBGSuccinct)
    uint64_t left_roots = 0;
};

/**
 * The segment search at the graph level by its definitions: window-major over the orientations
 * in plan order, a window admitted when it is above the floor and has at most max_anchors
 * anchors; per anchor in node order, owned when an earlier root's right walks reach its k-mer at
 * the window's start, else a root whose extent is [least start of its left walks (not below
 * the horizon: the last admitted window's start + 1), greatest end of its right walks); reported
 * when its end is >= its window's start + W'. The work as the memoised arms do it: each distinct
 * state expanded once (the right arm's per orientation, the left arm's per window).
 */
inline OracleRun oracle_segments(const OracleGraph &g, const std::vector<OracleQ> &searched,
                                 const OraclePlan &plan) {
    typedef std::pair<std::string, uint32_t> State;
    const size_t k = g.k;
    const uint32_t L = static_cast<uint32_t>(searched.front().length);
    OracleRun run;
    run.window = std::min<uint32_t>(plan.window ? plan.window : k, L);
    const std::vector<uint32_t> grid = oracle_window_grid(L, run.window, plan.step, plan.choice);

    for (const OracleQ &q : searched) {
        std::vector<OracleWindow> &windows = run.windows[q.orientation];
        for (uint32_t s : grid) {
            OracleWindow w;
            w.offset = s;
            w.bits = q.bits(s, k);
            w.below_floor = w.bits < plan.min_information_bits;
            for (const auto &[kmer, node] : g.kmers) {
                if (q.admits(kmer, s))
                    w.anchors.emplace_back(node, kmer);
            }
            std::sort(w.anchors.begin(), w.anchors.end());
            w.admitted = !w.below_floor && w.anchors.size() <= plan.max_anchors;
            windows.push_back(w);
        }
    }

    std::map<Orientation, std::set<State>> right_seen;
    std::map<Orientation, std::set<State>> right_expanded;
    std::map<Orientation, uint32_t> horizon;
    for (size_t i = 0; i < grid.size(); ++i) {
        for (const OracleQ &q : searched) {
            OracleWindow &w = run.windows[q.orientation][i];
            if (!w.admitted)
                continue;
            const uint32_t s = w.offset;
            const uint32_t h = horizon[q.orientation];
            std::set<State> left_seen;
            std::set<State> left_expanded;
            for (const auto &[node, kmer] : w.anchors) {
                if (right_seen[q.orientation].count({ kmer, s })) {
                    ++w.owned;
                    continue;
                }
                ++w.roots;
                // left: layers of k-mers by offset, down to the horizon
                uint32_t from = s;
                std::set<std::string> layer { kmer };
                for (uint32_t t = s; !layer.empty(); --t) {
                    from = t;
                    std::set<std::string> next;
                    for (const std::string &u : layer) {
                        if (!t)
                            continue;
                        const bool fresh = left_expanded.insert({ u, t }).second;
                        if (fresh) {
                            run.left_expanded.push_back(u);
                            run.left_edges += g.incoming(u).size();
                        }
                        for (const std::string &x : g.incoming(u)) {
                            if (!q.admits(x[0] + u, t - 1))
                                continue;
                            if (t - 1 < h) {
                                run.left_cuts += fresh;
                                continue;
                            }
                            next.insert(x);
                            left_seen.insert({ x, t - 1 });
                        }
                    }
                    layer.swap(next);
                    if (!t)
                        break;
                }
                // right: layers up to L
                uint32_t to = s + static_cast<uint32_t>(k);
                layer = { kmer };
                for (uint32_t t = s; !layer.empty(); ++t) {
                    to = t + static_cast<uint32_t>(k);
                    std::set<std::string> next;
                    for (const std::string &u : layer) {
                        if (t + k == L)
                            continue;
                        const bool fresh = right_expanded[q.orientation].insert({ u, t }).second;
                        uint64_t allowed = 0;
                        for (const std::string &v : g.outgoing(u)) {
                            if (!q.admits(u + v.back(), t))
                                continue;
                            ++allowed;
                            next.insert(v);
                            right_seen[q.orientation].insert({ v, t + 1 });
                        }
                        if (fresh) {
                            run.right_edges += g.outgoing(u).size();
                            run.right_branches += allowed > 1;
                        }
                    }
                    layer.swap(next);
                    if (t + k == L)
                        break;
                }
                const bool reported = to >= s + run.window;
                run.extents.push_back({ q.orientation, s, node, from, to, reported });
                if (reported) {
                    const auto [pf, pt] = graph::pattern::pattern_interval(q.orientation, L, from,
                                                                           to);
                    ++run.segments[{ q.orientation, pf, pt }];
                }
                if (s)
                    ++run.left_roots;
            }
            run.states += left_seen.size();
            horizon[q.orientation] = s + 1;
        }
    }
    for (const auto &[o, seen] : right_seen) {
        run.states += seen.size();
    }
    return run;
}

} // namespace segments
} // namespace test
} // namespace mtg

#endif // __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__
