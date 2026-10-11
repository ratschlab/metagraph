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
 *  - O4, the peptide residue oracle: the records translated in their six frames with this
 *    file's own genetic codes and each window of residues compared letter by letter with the
 *    peptide (a residue's codon of another residue's sense is one substitution; X and '*' are
 *    never substituted, a stop never substitutes), and the same rule on the graph's k-mers
 *    codon by codon (a codon cut by a window's edge: the least level over the codons agreeing
 *    with its bases inside), with each context's level and substituted residues;
 *  - O5, the brute-force ball: every string over the pattern's positions enumerated and its
 *    substitutions counted: |ball| and the bits left, against the engine's log-space DP (a
 *    peptide's by O4's rule, cut windows included).
 * The IUPAC table, the genetic codes, the reverse complement and palindromy are written out
 * here; Pattern is read only by O3, which hands the engine exact DNA strings.
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


// ---------------------------------------------------------------- O4, the peptide residue oracle

/**
 * The standard code (NCBI table 1) by residue, '*' the stops, as textbooks list it: written
 * apart from the engine's 64-letter strings.
 */
inline const std::map<char, std::vector<std::string>>& standard_code() {
    static const std::map<char, std::vector<std::string>> code {
        { 'A', { "GCT", "GCC", "GCA", "GCG" } },
        { 'C', { "TGT", "TGC" } },
        { 'D', { "GAT", "GAC" } },
        { 'E', { "GAA", "GAG" } },
        { 'F', { "TTT", "TTC" } },
        { 'G', { "GGT", "GGC", "GGA", "GGG" } },
        { 'H', { "CAT", "CAC" } },
        { 'I', { "ATT", "ATC", "ATA" } },
        { 'K', { "AAA", "AAG" } },
        { 'L', { "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" } },
        { 'M', { "ATG" } },
        { 'N', { "AAT", "AAC" } },
        { 'P', { "CCT", "CCC", "CCA", "CCG" } },
        { 'Q', { "CAA", "CAG" } },
        { 'R', { "CGT", "CGC", "CGA", "CGG", "AGA", "AGG" } },
        { 'S', { "TCT", "TCC", "TCA", "TCG", "AGT", "AGC" } },
        { 'T', { "ACT", "ACC", "ACA", "ACG" } },
        { 'V', { "GTT", "GTC", "GTA", "GTG" } },
        { 'W', { "TGG" } },
        { 'Y', { "TAT", "TAC" } },
        { '*', { "TAA", "TAG", "TGA" } },
    };
    return code;
}

/**
 * Every other NCBI table as the codons it codes differently from the standard code (the NCBI
 * Taxonomy's "Differences from the Standard Code"; for 27, 28 and 31 the residue NCBI's ncbieaa
 * gives a codon that can also end translation in context, which codes that residue).
 */
inline const std::map<int, std::vector<std::pair<std::string, char>>>& code_differences() {
    static const std::map<int, std::vector<std::pair<std::string, char>>> diffs {
        { 1, {} },
        { 2, { { "AGA", '*' }, { "AGG", '*' }, { "ATA", 'M' }, { "TGA", 'W' } } },
        { 3, { { "ATA", 'M' }, { "CTT", 'T' }, { "CTC", 'T' }, { "CTA", 'T' }, { "CTG", 'T' },
               { "TGA", 'W' } } },
        { 4, { { "TGA", 'W' } } },
        { 5, { { "AGA", 'S' }, { "AGG", 'S' }, { "ATA", 'M' }, { "TGA", 'W' } } },
        { 6, { { "TAA", 'Q' }, { "TAG", 'Q' } } },
        { 9, { { "AAA", 'N' }, { "AGA", 'S' }, { "AGG", 'S' }, { "TGA", 'W' } } },
        { 10, { { "TGA", 'C' } } },
        { 11, {} },
        { 12, { { "CTG", 'S' } } },
        { 13, { { "AGA", 'G' }, { "AGG", 'G' }, { "ATA", 'M' }, { "TGA", 'W' } } },
        { 14, { { "AAA", 'N' }, { "AGA", 'S' }, { "AGG", 'S' }, { "TAA", 'Y' },
                { "TGA", 'W' } } },
        { 15, { { "TAG", 'Q' } } },
        { 16, { { "TAG", 'L' } } },
        { 21, { { "TGA", 'W' }, { "ATA", 'M' }, { "AGA", 'S' }, { "AGG", 'S' },
                { "AAA", 'N' } } },
        { 22, { { "TCA", '*' }, { "TAG", 'L' } } },
        { 23, { { "TTA", '*' } } },
        { 24, { { "AGA", 'S' }, { "AGG", 'K' }, { "TGA", 'W' } } },
        { 25, { { "TGA", 'G' } } },
        { 26, { { "CTG", 'A' } } },
        { 27, { { "TAG", 'Q' }, { "TAA", 'Q' }, { "TGA", 'W' } } },
        { 28, { { "TAA", 'Q' }, { "TAG", 'Q' }, { "TGA", 'W' } } },
        { 29, { { "TAA", 'Y' }, { "TAG", 'Y' } } },
        { 30, { { "TAA", 'E' }, { "TAG", 'E' } } },
        { 31, { { "TGA", 'W' }, { "TAG", 'E' }, { "TAA", 'E' } } },
        { 32, { { "TAG", 'W' } } },
        { 33, { { "TAA", 'Y' }, { "TGA", 'W' }, { "AGA", 'S' }, { "AGG", 'K' } } },
    };
    return diffs;
}

// one genetic code as the oracles read it: codon -> residue ('*' a stop)
struct OracleCode {
    int id;
    std::map<std::string, char> residue;

    // '?' for anything but a codon of A, C, G, T
    char translate(std::string_view codon) const {
        auto it = residue.find(std::string(codon));
        return it == residue.end() ? '?' : it->second;
    }
};

inline const OracleCode& oracle_code(int id) {
    static std::map<int, OracleCode> cache;
    auto it = cache.find(id);
    if (it != cache.end())
        return it->second;
    OracleCode code { id, {} };
    for (const auto &[residue, codons] : standard_code()) {
        for (const std::string &codon : codons) {
            code.residue[codon] = residue;
        }
    }
    assert(code.residue.size() == 64);
    for (const auto &[codon, residue] : code_differences().at(id)) {
        code.residue[codon] = residue;
    }
    return cache.emplace(id, std::move(code)).first->second;
}

inline std::vector<int> oracle_code_ids() {
    std::vector<int> ids;
    for (const auto &[id, diffs] : code_differences()) {
        ids.push_back(id);
    }
    return ids;
}

// every codon, in A C G T order
inline const std::vector<std::string>& all_codons() {
    static const std::vector<std::string> codons = []() {
        std::vector<std::string> result;
        for (char a : std::string("ACGT")) {
            for (char b : std::string("ACGT")) {
                for (char c : std::string("ACGT")) {
                    result.push_back({ a, b, c });
                }
            }
        }
        return result;
    }();
    return codons;
}

// a peptide letter admits a coded residue without a substitution: X any residue, B D or N, Z E
// or Q, J I or L, a residue itself, '*' a stop; no letter but '*' admits a stop
inline bool letter_admits(char letter, char coded) {
    if (letter == '*')
        return coded == '*';
    if (coded == '*' || coded == '?')
        return false;
    switch (letter) {
        case 'X': return true;
        case 'B': return coded == 'D' || coded == 'N';
        case 'Z': return coded == 'E' || coded == 'Q';
        case 'J': return coded == 'I' || coded == 'L';
        default: return coded == letter;
    }
}

/**
 * The level of a codon coding |coded| at a residue |letter| (SPEC §21.3): 0 when the letter
 * admits it, 1 when it is a sense codon of another residue and the letter can be substituted
 * (every letter but X and '*'), nullopt otherwise (a stop where the letter is no '*', a sense
 * codon at '*': no instance at any level).
 */
inline std::optional<uint8_t> codon_level(char letter, char coded) {
    if (letter_admits(letter, coded))
        return 0;
    if (letter != 'X' && letter != '*' && coded != '*' && coded != '?')
        return 1;
    return std::nullopt;
}

/**
 * The level of |piece|, the bases of a codon of residue |letter| at the codon positions
 * [first, first + |piece|): a whole codon by codon_level, a cut one by the least level over
 * the codons agreeing with it there (the least count any completion gives); nullopt when no
 * codon of either kind agrees.
 */
inline std::optional<uint8_t> piece_level(char letter, const OracleCode &code, size_t first,
                                          std::string_view piece) {
    if (piece.size() == 3)
        return codon_level(letter, code.translate(piece));
    // (a cut piece's level is looked up again for every string of a brute-force ball)
    static std::map<std::tuple<int, char, size_t, std::string>, std::optional<uint8_t>> cache;
    const auto key = std::make_tuple(code.id, letter, first, std::string(piece));
    if (auto it = cache.find(key); it != cache.end())
        return it->second;
    std::optional<uint8_t> least;
    for (const std::string &codon : all_codons()) {
        if (std::string_view(codon).substr(first, piece.size()) != piece)
            continue;
        const auto level = codon_level(letter, code.translate(codon));
        if (level && (!least || *level < *least))
            least = level;
    }
    cache.emplace(key, least);
    return least;
}

// the substituted residues of an instance, in P's coordinates (ascending); its level is their
// number
typedef std::vector<uint32_t> Residues;

/**
 * |s| as the positions [begin, begin + |s|) of the instances of |peptide| (P's letters): the
 * residues whose piece is substituted, nullopt when a piece is no instance at any level.
 */
inline std::optional<Residues> peptide_residues(const std::string &peptide,
                                                const OracleCode &code, size_t begin,
                                                std::string_view s) {
    const size_t end = begin + s.size();
    if (end > 3 * peptide.size())
        return std::nullopt;
    Residues residues;
    for (size_t i = begin / 3; 3 * i < end; ++i) {
        const size_t lo = std::max(begin, 3 * i);
        const size_t hi = std::min(end, 3 * i + 3);
        const auto level = piece_level(peptide[i], code, lo - 3 * i,
                                       s.substr(lo - begin, hi - lo));
        if (!level)
            return std::nullopt;
        if (*level)
            residues.push_back(static_cast<uint32_t>(i));
    }
    return residues;
}

/**
 * The same for the oriented peptide of |orientation|: REVERSE reads rc(P), whose positions
 * [b, e) are the reverse complement of P's [L - e, L - b), so its residues are read on P.
 */
inline std::optional<Residues> oriented_residues(const std::string &peptide,
                                                 const OracleCode &code,
                                                 Orientation orientation, size_t begin,
                                                 std::string_view s) {
    if (orientation != Orientation::REVERSE)
        return peptide_residues(peptide, code, begin, s);
    const size_t L = 3 * peptide.size();
    if (begin + s.size() > L)
        return std::nullopt;
    return peptide_residues(peptide, code, L - begin - s.size(), reverse_complement_string(s));
}

// the codons a letter admits without a substitution
inline std::set<std::string> letter_codons(char letter, const OracleCode &code) {
    std::set<std::string> codons;
    for (const std::string &codon : all_codons()) {
        if (letter_admits(letter, code.translate(codon)))
            codons.insert(codon);
    }
    return codons;
}

// P's instances equal rc(P)'s: residue i's codons are the reverse complements of residue
// m - 1 - i's, for every i (the instances are a product of codon sets)
inline bool peptide_palindromic(const std::string &peptide, const OracleCode &code) {
    for (size_t i = 0, j = peptide.size(); i < peptide.size(); ++i) {
        std::set<std::string> mirrored;
        for (const std::string &codon : letter_codons(peptide[--j], code)) {
            mirrored.insert(reverse_complement_string(codon));
        }
        if (letter_codons(peptide[i], code) != mirrored)
            return false;
    }
    return true;
}

inline std::vector<Orientation> peptide_orientations(const std::string &peptide,
                                                     const OracleCode &code, Strands strands) {
    if (peptide_palindromic(peptide, code))
        return { Orientation::PALINDROMIC };
    std::vector<Orientation> result;
    if (strands != Strands::REVERSE)
        result.push_back(Orientation::FORWARD);
    if (strands != Strands::FORWARD)
        result.push_back(Orientation::REVERSE);
    return result;
}

/**
 * O4 on the graph: every context of |peptide| (table |table|) within |m| substituted residues
 * on |graph|, in answer order, its positions the substituted residues. L <= k: every k-mer,
 * offset of the scope and orientation, the k-mer's bases at the offset read codon by codon; L >
 * k: the anchors, every k-mer (offset 0) holding the anchor window of its orientation (FORWARD
 * P[0, k), REVERSE rc(P)[0, k)) within |m|, a codon cut by the window's edge by its bases
 * inside.
 */
inline std::vector<MCtx> peptide_walk(const DeBruijnGraph &graph, const std::string &peptide,
                                      int table, Scope scope, Strands strands, uint8_t m) {
    const OracleCode &code = oracle_code(table);
    const size_t k = graph.get_k();
    const size_t L = 3 * peptide.size();
    const size_t window = std::min(L, k);
    const std::vector<uint32_t> offsets = L > k ? std::vector<uint32_t> { 0 }
                                                : scope_offsets(L, k, scope);
    std::vector<MCtx> result;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        for (Orientation orientation : peptide_orientations(peptide, code, strands)) {
            for (uint32_t p : offsets) {
                const auto residues = oriented_residues(peptide, code, orientation, 0,
                                                        std::string_view(kmer).substr(p, window));
                if (residues && residues->size() <= m) {
                    result.push_back(MCtx { node, p, orientation,
                                            static_cast<uint8_t>(residues->size()), *residues });
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

// the translation of |s| from |frame|: one letter per whole codon ('?' if it holds a symbol
// other than A, C, G, T)
inline std::string translate_frame(std::string_view s, size_t frame, const OracleCode &code) {
    std::string protein;
    for (size_t i = frame; i + 3 <= s.size(); i += 3) {
        protein.push_back(code.translate(s.substr(i, 3)));
    }
    return protein;
}

/**
 * O4 on the records (L <= k): every strand sequence (the records, and on a CANONICAL or
 * PRIMARY graph their reverse complements too) translated in its three frames, and those of
 * its reverse complement for the reverse orientation; every window of r residues compared with
 * the peptide letter by letter (codon_level on the translated residue); each instance within
 * |m| at its start, and every indexed k-mer holding it at an offset of the scope, as (spelled
 * k-mer, offset, orientation, level). Never reads the graph.
 */
inline SpelledMCtx peptide_records(const std::vector<std::string> &records, size_t k,
                                   DeBruijnGraph::Mode mode, const std::string &peptide,
                                   int table, Scope scope, Strands strands, uint8_t m) {
    const OracleCode &code = oracle_code(table);
    const size_t r = peptide.size();
    const size_t L = 3 * r;
    assert(L <= k);
    std::vector<std::string> sequences = records;
    if (mode != DeBruijnGraph::BASIC) {
        for (const std::string &record : records) {
            sequences.push_back(reverse_complement_string(record));
        }
    }
    // the starts in |s| of the peptide's instances within m, with their levels
    auto instances = [&](const std::string &s) {
        std::vector<std::pair<size_t, uint8_t>> found;
        for (size_t frame = 0; frame < 3; ++frame) {
            const std::string protein = translate_frame(s, frame, code);
            for (size_t start = 0; start + r <= protein.size(); ++start) {
                uint32_t level = 0;
                bool instance = true;
                for (size_t i = 0; i < r && instance; ++i) {
                    const auto l = codon_level(peptide[i], protein[start + i]);
                    instance = l.has_value();
                    level += l.value_or(0);
                }
                if (instance && level <= m)
                    found.emplace_back(frame + 3 * start, static_cast<uint8_t>(level));
            }
        }
        return found;
    };
    SpelledMCtx result;
    for (const std::string &s : sequences) {
        for (Orientation o : peptide_orientations(peptide, code, strands)) {
            std::vector<std::pair<size_t, uint8_t>> at;
            if (o != Orientation::REVERSE) {
                at = instances(s);
            } else {
                // an instance of P in rc(s) at j is one of rc(P) in s at |s| - j - L
                for (const auto &[j, level] : instances(reverse_complement_string(s))) {
                    at.emplace_back(s.size() - j - L, level);
                }
            }
            for (const auto &[i, level] : at) {
                for (uint32_t p : scope_offsets(L, k, scope)) {
                    if (i < p || i - p + k > s.size())
                        continue;
                    const std::string kmer = s.substr(i - p, k);
                    if (kmer.find_first_not_of("ACGT") == std::string::npos)
                        result.emplace(kmer, p, o, level);
                }
            }
        }
    }
    return result;
}

/**
 * O5 for a peptide: per number of substituted residues u, the strings over the positions
 * [begin, end) of |peptide|'s instances (every string of A, C, G, T enumerated, 4^(end - begin)
 * of them, end - begin <= 10) whose level by O4's rule is u; a string that is no instance is in
 * no entry.
 */
inline std::vector<uint64_t> peptide_ball_histogram(const std::string &peptide, int table,
                                                    size_t begin, size_t end) {
    const OracleCode &code = oracle_code(table);
    const size_t n = end - begin;
    assert(n <= 10);
    std::vector<uint64_t> histogram(peptide.size() + 1, 0);
    std::string s(n, 'A');
    for (uint64_t c = 0; c < (uint64_t(1) << (2 * n)); ++c) {
        uint64_t x = c;
        for (size_t i = 0; i < n; ++i, x >>= 2) {
            s[i] = "ACGT"[x & 3];
        }
        if (const auto residues = peptide_residues(peptide, code, begin, s))
            ++histogram[residues->size()];
    }
    return histogram;
}

} // namespace mismatch
} // namespace test
} // namespace mtg

#endif // __METAGRAPH_TESTS_PATTERN_MISMATCH_ORACLES_HPP__
