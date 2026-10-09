#ifndef __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__
#define __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__

/**
 * Oracles of the segment search, independent of the engine:
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
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "graph/alignment/genetic_code.hpp"
#include "graph/alignment/pattern_search.hpp"


namespace mtg {
namespace test {
namespace segments {

using graph::DeBruijnGraph;
using graph::pattern::Orientation;
using graph::pattern::PatternKind;
using graph::pattern::Strands;
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

} // namespace segments
} // namespace test
} // namespace mtg

#endif // __METAGRAPH_TESTS_PATTERN_SEGMENT_ORACLES_HPP__
