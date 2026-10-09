#ifndef __METAGRAPH_TESTS_PATTERN_MISMATCH_ORACLES_HPP__
#define __METAGRAPH_TESTS_PATTERN_MISMATCH_ORACLES_HPP__

/**
 * The oracles of bounded mismatches (SPEC §21; TESTS.md of the mismatch plan, O1, O2, O3, O5):
 * brute-force Hamming neighbourhoods of a DNA or IUPAC pattern that share nothing with the
 * engine but the graph.
 *  - O1, the Hamming graph-walk: every k-mer of the served graph, every offset of the scope,
 *    every orientation searched, the substitutions counted position by position with this
 *    file's own IUPAC table: the contexts within m, each with its level and its substituted
 *    positions in the submitted pattern's coordinates, in answer order (level, node, offset,
 *    orientation).
 *  - O2, the Hamming record-scan: the same over the records' k-mers (both orientations unless
 *    BASIC), as (spelled k-mer, offset, orientation, level): O1 independently of the graph's
 *    ids.
 *  - O3, the candidate enumeration (the probe baseline of docs/NOTE-mismatch-search.md): every
 *    string of the ball searched as an exact DNA pattern by the engine without mismatches, the
 *    contexts united, each at the level of the string that found it: the level DFS against a
 *    different algorithm.
 *  - O5, the brute-force ball: every string over the pattern's positions enumerated and its
 *    substitutions counted: |ball| and the bits left, against the engine's log-space DP.
 * The IUPAC table, the reverse complement and palindromy are written out here; Pattern is read
 * only by O3, which hands the engine exact DNA strings.
 */

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace test {
namespace mismatch {

using graph::DeBruijnGraph;
using graph::DBGSuccinct;
using graph::CanonicalDBG;
using graph::pattern::Orientation;
using graph::pattern::Scope;
using graph::pattern::Strands;

typedef DeBruijnGraph::node_index node_index;


// ---------------------------------------------------------------- the oracles' IUPAC table

// the bases each IUPAC code admits, sorted; empty for anything else
inline std::string iupac_bases(char code) {
    switch (std::toupper(static_cast<unsigned char>(code))) {
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

inline char complement_base(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return '?';
    }
}

inline std::string reverse_complement_string(std::string_view s) {
    std::string rc(s.rbegin(), s.rend());
    for (char &c : rc) {
        c = complement_base(c);
    }
    return rc;
}

// an oriented pattern as the oracles read it: per position the bases it admits (sorted)
typedef std::vector<std::string> Positions;

inline Positions read_pattern(std::string_view text) {
    Positions q;
    for (char c : text) {
        q.push_back(iupac_bases(c));
        assert(!q.back().empty());
    }
    return q;
}

inline Positions reverse_complement(const Positions &q) {
    Positions rc(q.rbegin(), q.rend());
    for (std::string &bases : rc) {
        for (char &base : bases) {
            base = complement_base(base);
        }
        std::sort(bases.begin(), bases.end());
    }
    return rc;
}

inline bool palindromic(const Positions &q) {
    return q == reverse_complement(q);
}

// the orientations a request searches, each with its oriented pattern: a palindromic P once
inline std::vector<std::pair<Orientation, Positions>> orientations(std::string_view text,
                                                                  Strands strands) {
    const Positions q = read_pattern(text);
    if (palindromic(q))
        return { { Orientation::PALINDROMIC, q } };
    std::vector<std::pair<Orientation, Positions>> result;
    if (strands != Strands::REVERSE)
        result.emplace_back(Orientation::FORWARD, q);
    if (strands != Strands::FORWARD)
        result.emplace_back(Orientation::REVERSE, reverse_complement(q));
    return result;
}

inline std::vector<uint32_t> scope_offsets(size_t L, size_t k, Scope scope) {
    if (scope == Scope::SUFFIX)
        return { static_cast<uint32_t>(k - L) };
    std::vector<uint32_t> offsets;
    for (uint32_t p = 0; p + L <= k; ++p) {
        offsets.push_back(p);
    }
    return offsets;
}

/**
 * The positions of |q| where |s| holds a base outside the position's set (an N position never
 * does), ascending; nullopt when |s| is not an instance at all: another length, or a symbol
 * other than A, C, G, T at a pattern position (the sentinel '$', or a record's N, never fill
 * one: SPEC §21.3).
 */
inline std::optional<std::vector<uint32_t>> substituted(const Positions &q, std::string_view s) {
    if (s.size() != q.size())
        return std::nullopt;
    std::vector<uint32_t> positions;
    for (size_t i = 0; i < q.size(); ++i) {
        if (std::string("ACGT").find(s[i]) == std::string::npos)
            return std::nullopt;
        if (q[i].find(s[i]) == std::string::npos)
            positions.push_back(static_cast<uint32_t>(i));
    }
    return positions;
}

// positions of the oriented pattern of |orientation| in the submitted pattern's coordinates
// (a REVERSE context instantiates rc(P): its position j is P's L - 1 - j), ascending
inline std::vector<uint32_t> in_submitted(std::vector<uint32_t> positions,
                                          Orientation orientation, size_t L) {
    if (orientation == Orientation::REVERSE) {
        for (uint32_t &i : positions) {
            i = static_cast<uint32_t>(L - 1 - i);
        }
        // (ascending positions mirrored are descending)
        std::reverse(positions.begin(), positions.end());
    }
    return positions;
}


// ---------------------------------------------------------------- O1, the Hamming graph-walk

// one context within m: (orientation, node, offset), its level and its substituted positions
// in P's coordinates
struct MCtx {
    node_index node;
    uint32_t offset;
    Orientation orientation;
    uint8_t level;
    std::vector<uint32_t> positions = {};

    // the answer order of a request with mismatches: (level, node, offset, orientation)
    bool operator<(const MCtx &o) const {
        return std::tie(level, node, offset, orientation)
                < std::tie(o.level, o.node, o.offset, o.orientation);
    }
    bool operator==(const MCtx &o) const {
        return std::tie(level, node, offset, orientation, positions)
                == std::tie(o.level, o.node, o.offset, o.orientation, o.positions);
    }
};

inline std::ostream& operator<<(std::ostream &out, const MCtx &c) {
    out << "(L" << int(c.level) << ", " << c.node << ", " << c.offset << ", "
        << graph::pattern::orientation_key(c.orientation) << ", [";
    for (uint32_t i : c.positions) {
        out << " " << i;
    }
    return out << " ])";
}

inline const DBGSuccinct& stored_graph(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// every k-mer of the served graph with its node id: the valid edges of the DBGSuccinct and,
// on a wrapped PRIMARY graph, the wrapper's reverse complements (a palindrome once)
inline std::vector<std::pair<node_index, std::string>> graph_kmers(const DeBruijnGraph &graph) {
    std::vector<std::pair<node_index, std::string>> kmers;
    const DBGSuccinct &dbg_succ = stored_graph(graph);
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
    for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
        if (!dbg_succ.in_graph(y))
            continue;
        kmers.emplace_back(y, graph.get_node_sequence(y));
        if (canonical) {
            node_index z = canonical->reverse_complement(y);
            if (z != y)
                kmers.emplace_back(z, graph.get_node_sequence(z));
        }
    }
    return kmers;
}

/**
 * O1: every context of |text| within |m| substitutions on |graph| (on a graph without its mask,
 * pass the masked twin: its k-mers are the true ones, its ids those of the twin), in answer
 * order.
 */
inline std::vector<MCtx> hamming_walk(const DeBruijnGraph &graph, std::string_view text,
                                      Scope scope, Strands strands, uint8_t m) {
    const size_t k = graph.get_k();
    const size_t L = text.size();
    std::vector<MCtx> result;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        for (const auto &[orientation, q] : orientations(text, strands)) {
            for (uint32_t p : scope_offsets(L, k, scope)) {
                auto positions = substituted(q, std::string_view(kmer).substr(p, L));
                if (positions && positions->size() <= m) {
                    result.push_back(MCtx { node, p, orientation,
                                            static_cast<uint8_t>(positions->size()),
                                            in_submitted(*positions, orientation, L) });
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}


// ---------------------------------------------------------------- O2, the Hamming record-scan

typedef std::set<std::tuple<std::string, uint32_t, Orientation, uint8_t>> SpelledMCtx;

/**
 * O2: the contexts within |m| that the records' k-mers spell (both orientations of each unless
 * |mode| is BASIC; a k-mer holding a symbol other than A, C, G, T is not indexed), never reading
 * the graph.
 */
inline SpelledMCtx hamming_records(const std::vector<std::string> &records, size_t k,
                                   DeBruijnGraph::Mode mode, std::string_view text, Scope scope,
                                   Strands strands, uint8_t m) {
    std::set<std::string> kmers;
    for (const std::string &record : records) {
        for (size_t i = 0; i + k <= record.size(); ++i) {
            const std::string kmer = record.substr(i, k);
            if (kmer.find_first_not_of("ACGT") != std::string::npos)
                continue;
            kmers.insert(kmer);
            if (mode != DeBruijnGraph::BASIC)
                kmers.insert(reverse_complement_string(kmer));
        }
    }
    SpelledMCtx result;
    for (const std::string &kmer : kmers) {
        for (const auto &[orientation, q] : orientations(text, strands)) {
            for (uint32_t p : scope_offsets(text.size(), k, scope)) {
                auto positions = substituted(q, std::string_view(kmer).substr(p, text.size()));
                if (positions && positions->size() <= m)
                    result.emplace(kmer, p, orientation, positions->size());
            }
        }
    }
    return result;
}

inline SpelledMCtx spelled(const DeBruijnGraph &graph, const std::vector<MCtx> &contexts) {
    SpelledMCtx result;
    for (const MCtx &c : contexts) {
        result.emplace(graph.get_node_sequence(c.node), c.offset, c.orientation, c.level);
    }
    return result;
}


// ---------------------------------------------------------------- O3, the candidate enumeration

// every string over A, C, G, T within |m| substitutions of |q|, with its number of them
inline void for_each_ball_string(const Positions &q, uint8_t m,
                                 const std::function<void(const std::string&, uint8_t)> &f) {
    std::string s(q.size(), '?');
    std::function<void(size_t, uint8_t)> fill = [&](size_t i, uint8_t used) {
        if (i == q.size()) {
            f(s, used);
            return;
        }
        for (char base : std::string("ACGT")) {
            const bool same = q[i].find(base) != std::string::npos;
            if (!same && used == m)
                continue;
            s[i] = base;
            fill(i + 1, used + !same);
        }
    };
    fill(0, 0);
}

/**
 * O3: every string of the ball of each oriented pattern searched as an exact DNA pattern, on
 * its own strand only, by |engine| at max_mismatches 0 (the search of v1); the contexts it
 * releases are those of the oriented pattern at the string's level. Each (node, offset) of an
 * orientation spells one string, so the union has no duplicate within an orientation.
 */
inline std::vector<MCtx> candidate_enumeration(const graph::pattern::PatternSearch &engine,
                                               std::string_view text, Scope scope,
                                               Strands strands, uint8_t m) {
    using namespace graph::pattern;
    const size_t L = text.size();
    std::vector<MCtx> result;
    for (const auto &oriented : orientations(text, strands)) {
        const Orientation orientation = oriented.first;
        const Positions &q = oriented.second;
        for_each_ball_string(q, m, [&](const std::string &s, uint8_t level) {
            Request request;
            request.scope = scope;
            request.strands = Strands::FORWARD;
            request.min_information_bits = 0;
            request.max_contexts = 1'000'000'000;
            request.mode = Mode::ALL_OR_COUNT;
            Budget budget(1'000'000'000, Deadline(Deadline::Clock::now(),
                                                  std::numeric_limits<double>::infinity(), 0));
            const Pattern exact = Pattern::parse(PatternKind::DNA, s);
            Result r = engine.enumerate(exact, request, budget, [&](const Context &c) {
                result.push_back(MCtx { c.node, c.offset, orientation, level,
                                        in_submitted(*substituted(q, s), orientation, L) });
            });
            if (r.refusal || r.stop || !r.extraction || !r.extraction->complete)
                throw std::runtime_error("O3: the exact search of " + s + " is incomplete");
        });
    }
    std::sort(result.begin(), result.end());
    return result;
}


// ---------------------------------------------------------------- O5, the brute-force ball

/**
 * O5: per number of substitutions u, the strings over A, C, G, T of |text|'s length holding
 * exactly u substitutions against it (every string enumerated: 4^L of them, L <= 11).
 */
inline std::vector<uint64_t> ball_histogram(std::string_view text) {
    const Positions q = read_pattern(text);
    const size_t L = q.size();
    assert(L <= 11);
    std::vector<uint64_t> histogram(L + 1, 0);
    std::string s(L, 'A');
    for (uint64_t code = 0; code < (uint64_t(1) << (2 * L)); ++code) {
        uint64_t c = code;
        for (size_t i = 0; i < L; ++i, c >>= 2) {
            s[i] = "ACGT"[c & 3];
        }
        ++histogram[substituted(q, s)->size()];
    }
    return histogram;
}

// the bits left within |m| substitutions: 2L - log2 |ball|
inline double ball_bits(const std::vector<uint64_t> &histogram, uint8_t m) {
    uint64_t ball = 0;
    for (size_t u = 0; u <= m && u < histogram.size(); ++u) {
        ball += histogram[u];
    }
    return 2.0 * static_cast<double>(histogram.size() - 1) - std::log2(static_cast<double>(ball));
}

// the ball of an exact DNA pattern of |L| bases within |m|: sum over u <= m of C(L, u) 3^u
inline double exact_ball(size_t L, uint8_t m) {
    double ball = 0;
    double choose = 1;
    double three = 1;
    for (size_t u = 0; u <= m && u <= L; ++u) {
        ball += choose * three;
        choose = choose * static_cast<double>(L - u) / static_cast<double>(u + 1);
        three *= 3;
    }
    return ball;
}

} // namespace mismatch
} // namespace test
} // namespace mtg

#endif // __METAGRAPH_TESTS_PATTERN_MISMATCH_ORACLES_HPP__
