/**
 * The oracle suite of peptide patterns (docs/DESIGN-pattern-search.md §6, §13 increment 5;
 * owner decision #15 of 2026-10-08): the genetic codes, the codon automaton, and the engine
 * on tiny graphs against two brute-force oracles that never read the engine's automaton:
 *  - the record-scan oracle translates every record in all six frames (the three frames of
 *    the record and of its reverse complement) with this file's own genetic code and lists
 *    the peptide's occurrences on both strands: an occurrence of P in the reverse complement
 *    is an occurrence of rc(P) in the record. The graph contexts (L <= k: every indexed
 *    k-mer of the record holding an occurrence, with its offset) and the paths (L > k: the
 *    occurrences whose k-mers were all retained) follow from them;
 *  - the graph-walk oracle enumerates the built graph's k-mers (and for L > k its walks along
 *    the k-mers present) and keeps every string whose bases spell the peptide codon by codon
 *    in the same table (a codon cut by a window's edge: the prefix or suffix of a codon of
 *    its residue, compared as strings).
 * This file's genetic code is its own: the standard code written out by residue, as
 * textbooks list it, and every other NCBI table as the differences from it that NCBI
 * documents; GeneticCodeTables checks the engine's 64-letter tables against it, codon by
 * codon. The residue classes X, B, Z, J are re-derived here from the table. What the
 * graph-walk oracle shares with the engine: the graph (in_graph, node ids); the record-scan
 * oracle shares neither, and the two must agree wherever no k-mer was pruned and no path
 * joins records.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"
#include "pattern_test_support.hpp"

#include "common/vectors/bit_vector_dyn.hpp"
#include "graph/alignment/genetic_code.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

// the nucleotide builds the engine serves, as in test_pattern_search.cpp
#if _DNA_GRAPH || _DNA5_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;
using mtg::test::build_graph_batch;
using mtg::test::codon_string;
using mtg::test::unbounded_deadline;

typedef DeBruijnGraph::node_index node_index;

constexpr uint64_t kManySteps = 1'000'000'000;


// ---------------------------------------------------------------- the oracles' sequences

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return base;
    }
}

std::string rev_comp(const std::string &s) {
    std::string r(s.rbegin(), s.rend());
    for (char &c : r) {
        c = complement(c);
    }
    return r;
}

// the record symbols the build indexes (§3): A, C, G, T, and N on a DNA5 build
bool indexed(std::string_view s) {
    return std::all_of(s.begin(), s.end(), [](char c) {
#if _DNA5_GRAPH
        return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
#else
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
#endif
    });
}

// every codon, in A C G T order
const std::vector<std::string>& all_codons() {
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


// ---------------------------------------------------------------- the oracles' genetic codes

/**
 * The standard code (NCBI table 1) by residue, '*' the stops: this file's reading of it,
 * written independently of the engine's 64-letter strings.
 */
const std::map<char, std::vector<std::string>>& standard_code() {
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
 * Every other NCBI table as the codons it codes differently from the standard code (the
 * NCBI Taxonomy's "Differences from the Standard Code"; for 27, 28 and 31 the residue NCBI's
 * ncbieaa gives a codon that can also end translation in context).
 */
const std::map<int, std::vector<std::pair<std::string, char>>>& differences() {
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

// one table as the oracles read it: codon -> residue ('*' a stop)
struct OracleCode {
    int id;
    std::map<std::string, char> residue;

    // '?' for anything but a codon of A, C, G, T
    char translate(std::string_view codon) const {
        auto it = residue.find(std::string(codon));
        return it == residue.end() ? '?' : it->second;
    }
};

const OracleCode& oracle_code(int id) {
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
    EXPECT_EQ(64u, code.residue.size());
    for (const auto &[codon, residue] : differences().at(id)) {
        EXPECT_NE(code.residue[codon], residue) << id << " " << codon;
        code.residue[codon] = residue;
    }
    return cache.emplace(id, std::move(code)).first->second;
}

// a peptide letter admits a coded residue: X any, B D or N, Z E or Q, J I or L; never a
// stop, which only the letter '*' admits (owner decision #19 of 2026-10-08: a stop of the
// table, the codons this file's table codes '*')
bool admits(char letter, char coded) {
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

// the codons a peptide letter admits in |code|
const std::vector<std::string>& codons_of(char letter, const OracleCode &code) {
    static std::map<std::pair<int, char>, std::vector<std::string>> cache;
    auto [it, inserted] = cache.emplace(std::make_pair(code.id, letter),
                                        std::vector<std::string>());
    if (inserted) {
        for (const std::string &codon : all_codons()) {
            if (admits(letter, code.translate(codon)))
                it->second.push_back(codon);
        }
    }
    return it->second;
}

/**
 * |s| instantiates the positions [begin, begin + |s|) of the peptide's DNA instances: the
 * piece of |s| over each residue's positions equals the same positions of one of its codons
 * (a whole codon: a codon the residue admits; one cut by the edge: its prefix or suffix).
 */
bool instantiates(const std::string &peptide, const OracleCode &code, size_t begin,
                  std::string_view s) {
    const size_t end = begin + s.size();
    if (end > 3 * peptide.size())
        return false;
    for (size_t i = begin / 3; 3 * i < end; ++i) {
        const size_t lo = std::max(begin, 3 * i);
        const size_t hi = std::min(end, 3 * i + 3);
        const std::string_view piece = s.substr(lo - begin, hi - lo);
        bool any = false;
        for (const std::string &codon : codons_of(peptide[i], code)) {
            if (std::string_view(codon).substr(lo - 3 * i, hi - lo) == piece) {
                any = true;
                break;
            }
        }
        if (!any)
            return false;
    }
    return true;
}

// |s| instantiates the positions [begin, begin + |s|) of the oriented peptide: P itself, or
// rc(P), whose positions [b, e) are the reverse complement of P's [L - e, L - b)
bool oriented_instantiates(const std::string &peptide, const OracleCode &code,
                           Orientation orientation, size_t begin, const std::string &s) {
    if (orientation != Orientation::REVERSE)
        return instantiates(peptide, code, begin, s);
    const size_t L = 3 * peptide.size();
    if (begin + s.size() > L)
        return false;
    return instantiates(peptide, code, L - begin - s.size(), rev_comp(s));
}

// P's instances equal rc(P)'s: residue i's codons are the reverse complements of those of
// residue m - 1 - i, for every i (the instances are a product of codon sets)
bool oracle_palindromic(const std::string &peptide, const OracleCode &code) {
    for (size_t i = 0, j = peptide.size(); i < peptide.size(); ++i) {
        std::set<std::string> left(codons_of(peptide[i], code).begin(),
                                   codons_of(peptide[i], code).end());
        std::set<std::string> right;
        for (const std::string &codon : codons_of(peptide[--j], code)) {
            right.insert(rev_comp(codon));
        }
        if (left != right)
            return false;
    }
    return true;
}

std::vector<Orientation> oracle_orientations(const std::string &peptide,
                                             const OracleCode &code, Strands strands) {
    if (oracle_palindromic(peptide, code))
        return { Orientation::PALINDROMIC };
    std::vector<Orientation> result;
    if (strands != Strands::REVERSE)
        result.push_back(Orientation::FORWARD);
    if (strands != Strands::FORWARD)
        result.push_back(Orientation::REVERSE);
    return result;
}

// the translation of |s| from |frame|: one letter per whole codon ('?' if it holds a symbol
// other than A, C, G, T)
std::string translate_frame(const std::string &s, size_t frame, const OracleCode &code) {
    std::string protein;
    for (size_t i = frame; i + 3 <= s.size(); i += 3) {
        protein.push_back(code.translate(std::string_view(s).substr(i, 3)));
    }
    return protein;
}

// the starts in |s| of the peptide's instances, from the three frames' translations
std::vector<size_t> occurrences(const std::string &s, const std::string &peptide,
                                const OracleCode &code) {
    std::vector<size_t> result;
    for (size_t frame = 0; frame < 3; ++frame) {
        const std::string protein = translate_frame(s, frame, code);
        for (size_t r = 0; r + peptide.size() <= protein.size(); ++r) {
            bool all = true;
            for (size_t i = 0; i < peptide.size() && all; ++i) {
                all = admits(peptide[i], protein[r + i]);
            }
            if (all)
                result.push_back(frame + 3 * r);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

/**
 * The record-scan oracle: the six-frame translation of every record (and, on a CANONICAL
 * or PRIMARY graph, of its reverse complement as a record of its own: both orientations of
 * every k-mer are k-mers there), each instance of the oriented peptide at its start i in the
 * strand sequence s: FORWARD (and PALINDROMIC) the instances of P in s's frames, REVERSE
 * those of P in rc(s)'s frames (an instance of rc(P) in s at |s| - j - L).
 */
typedef std::vector<std::tuple<std::string, size_t, Orientation>> Occurrences;

Occurrences record_occurrences(const std::vector<std::string> &records,
                               DeBruijnGraph::Mode mode, const std::string &peptide,
                               const OracleCode &code, Strands strands) {
    std::vector<std::string> sequences = records;
    if (mode != DeBruijnGraph::BASIC) {
        for (const std::string &record : records) {
            sequences.push_back(rev_comp(record));
        }
    }
    const size_t L = 3 * peptide.size();
    Occurrences result;
    for (const std::string &s : sequences) {
        for (Orientation o : oracle_orientations(peptide, code, strands)) {
            if (o != Orientation::REVERSE) {
                for (size_t i : occurrences(s, peptide, code)) {
                    result.emplace_back(s, i, o);
                }
            } else {
                for (size_t j : occurrences(rev_comp(s), peptide, code)) {
                    result.emplace_back(s, s.size() - j - L, o);
                }
            }
        }
    }
    return result;
}

std::vector<uint32_t> scope_offsets(size_t L, size_t k, Scope scope) {
    if (L > k)
        return { 0 };
    if (scope == Scope::SUFFIX)
        return { static_cast<uint32_t>(k - L) };
    std::vector<uint32_t> offsets;
    for (uint32_t p = 0; p + L <= k; ++p) {
        offsets.push_back(p);
    }
    return offsets;
}

typedef std::set<std::tuple<std::string, uint32_t, Orientation>> SpelledContexts;

// L <= k: the indexed k-mers of the strand sequences holding an occurrence, at its offset
SpelledContexts record_contexts(const std::vector<std::string> &records, size_t k,
                                DeBruijnGraph::Mode mode, const std::string &peptide,
                                const OracleCode &code, const Request &request) {
    const size_t L = 3 * peptide.size();
    SpelledContexts result;
    for (const auto &[s, i, o] : record_occurrences(records, mode, peptide, code,
                                                     request.strands)) {
        for (uint32_t p : scope_offsets(L, k, request.scope)) {
            if (i < p || i - p + k > s.size())
                continue;
            std::string kmer = s.substr(i - p, k);
            if (indexed(kmer))
                result.emplace(kmer, p, o);
        }
    }
    return result;
}

typedef std::set<std::pair<std::string, Orientation>> SpelledPaths;

// L > k: the occurrences whose k-mers were all retained (§3, "Covered sequence")
SpelledPaths record_paths(const std::vector<std::string> &records, size_t k,
                          DeBruijnGraph::Mode mode, const std::string &peptide,
                          const OracleCode &code, Strands strands,
                          const std::set<std::string> &pruned = {}) {
    const size_t L = 3 * peptide.size();
    SpelledPaths result;
    for (const auto &[s, i, o] : record_occurrences(records, mode, peptide, code, strands)) {
        const std::string occurrence = s.substr(i, L);
        bool retained = indexed(occurrence);
        for (size_t j = 0; retained && j + k <= L; ++j) {
            retained = !pruned.count(occurrence.substr(j, k));
        }
        if (retained)
            result.emplace(occurrence, o);
    }
    return result;
}


// ---------------------------------------------------------------- graphs

std::shared_ptr<DeBruijnGraph> build(size_t k, const std::vector<std::string> &records,
                                     DeBruijnGraph::Mode mode, bool batch = false) {
    return batch ? build_graph_batch<DBGSuccinct>(k, records, mode)
                 : build_graph<DBGSuccinct>(k, records, mode);
}

const DBGSuccinct& base_dbg(const DeBruijnGraph &graph) {
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        return dynamic_cast<const DBGSuccinct&>(canonical->get_graph());
    return dynamic_cast<const DBGSuccinct&>(graph);
}

// every k-mer of the served graph with its node id (on a wrapped PRIMARY graph also the
// wrapper's reverse complements)
std::vector<std::pair<node_index, std::string>> graph_kmers(const DeBruijnGraph &graph) {
    std::vector<std::pair<node_index, std::string>> kmers;
    const DBGSuccinct &dbg_succ = base_dbg(graph);
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

// the stored node of a node of the served graph, by spelling (never the engine's arithmetic)
class StoredNodes {
  public:
    explicit StoredNodes(const DeBruijnGraph &graph)
          : graph_(graph), wrapped_(dynamic_cast<const CanonicalDBG*>(&graph)) {
        if (!wrapped_)
            return;
        const DBGSuccinct &dbg_succ = dynamic_cast<const DBGSuccinct&>(wrapped_->get_graph());
        for (node_index y = 1; y <= dbg_succ.max_index(); ++y) {
            if (dbg_succ.in_graph(y))
                stored_.emplace(dbg_succ.get_node_sequence(y), y);
        }
    }

    node_index of(node_index node) const {
        if (!wrapped_)
            return node;
        const std::string kmer = graph_.get_node_sequence(node);
        auto it = stored_.find(kmer);
        if (it == stored_.end())
            it = stored_.find(rev_comp(kmer));
        return it == stored_.end() ? DeBruijnGraph::npos : it->second;
    }

  private:
    const DeBruijnGraph &graph_;
    const CanonicalDBG *wrapped_;
    std::map<std::string, node_index> stored_;
};

Request make_request(Scope scope = Scope::ANY_OFFSET, Strands strands = Strands::BOTH,
                     bool extend_paths = false) {
    Request request;
    request.scope = scope;
    request.strands = strands;
    request.min_information_bits = 0;
    request.max_contexts = 1'000'000'000;
    request.max_anchors = 1'000'000'000;
    request.max_paths = 1'000'000'000;
    request.extend_paths = extend_paths;
    return request;
}

Pattern peptide_pattern(const std::string &peptide, int table) {
    return Pattern::parse(PatternKind::PROTEIN, peptide, GeneticCode::get(table));
}

Result count_of(const DeBruijnGraph &graph, const Pattern &pattern, const Request &request,
                uint64_t max_steps = kManySteps) {
    Budget budget(max_steps, unbounded_deadline());
    return PatternSearch(graph).count(pattern, request, budget);
}


// ---------------------------------------------------------------- the graph-walk oracle

struct Ctx {
    node_index node;
    uint32_t offset;
    Orientation orientation;

    bool operator<(const Ctx &o) const {
        return std::tie(node, offset, orientation) < std::tie(o.node, o.offset, o.orientation);
    }
    bool operator==(const Ctx &o) const {
        return std::tie(node, offset, orientation) == std::tie(o.node, o.offset, o.orientation);
    }
};

std::ostream& operator<<(std::ostream &out, const Ctx &c) {
    return out << "(" << c.node << ", " << c.offset << ", " << orientation_key(c.orientation)
               << ")";
}

// L <= k: every (node, offset, orientation) whose k-mer holds an instance at that offset;
// L > k: every anchor (offset 0) whose k-mer instantiates the anchor window. Answer order.
std::vector<Ctx> walk_contexts(const DeBruijnGraph &graph, const std::string &peptide,
                               const OracleCode &code, const Request &request) {
    const size_t k = graph.get_k();
    const size_t L = 3 * peptide.size();
    const size_t m = std::min(L, k);
    std::vector<Ctx> result;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        for (Orientation o : oracle_orientations(peptide, code, request.strands)) {
            for (uint32_t p : scope_offsets(L, k, request.scope)) {
                if (oriented_instantiates(peptide, code, o, 0, kmer.substr(p, m)))
                    result.push_back(Ctx { node, p, o });
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

struct PathCtx {
    node_index anchor;
    Orientation orientation;
    std::string sequence;
    std::vector<node_index> path;

    bool operator<(const PathCtx &o) const {
        return std::tie(anchor, orientation, sequence)
                < std::tie(o.anchor, o.orientation, o.sequence);
    }
    bool operator==(const PathCtx &o) const {
        return std::tie(anchor, orientation, sequence, path)
                == std::tie(o.anchor, o.orientation, o.sequence, o.path);
    }
};

std::ostream& operator<<(std::ostream &out, const PathCtx &p) {
    out << "(" << p.anchor << ", " << orientation_key(p.orientation) << ", " << p.sequence
        << ", [";
    for (node_index n : p.path) {
        out << " " << n;
    }
    return out << " ])";
}

struct PathOracle {
    std::vector<PathCtx> paths;
    std::map<Orientation, uint64_t> anchors;
    // the walks of k + 1 .. L bases whose spelling is a prefix of an instance
    uint64_t candidates = 0;
};

// L > k: every walk of L - k + 1 k-mers of the served graph spelling an instance of the
// oriented peptide, grown base by base along the k-mers present and kept while its spelling
// is a prefix of an instance (whole codons translated, the last one a codon prefix)
PathOracle walk_paths(const DeBruijnGraph &graph, const std::string &peptide,
                      const OracleCode &code, Strands strands) {
    const size_t k = graph.get_k();
    const size_t L = 3 * peptide.size();
    std::map<std::string, node_index> node_of;
    for (const auto &[node, kmer] : graph_kmers(graph)) {
        node_of.emplace(kmer, node);
    }
    PathOracle result;
    for (Orientation o : oracle_orientations(peptide, code, strands)) {
        std::function<void(std::string&, std::vector<node_index>&)> grow
                = [&](std::string &s, std::vector<node_index> &p) {
            if (s.size() == L) {
                result.paths.push_back(PathCtx { p.front(), o, s, p });
                return;
            }
            for (char b : std::string("ACGT")) {
                auto it = node_of.find(s.substr(s.size() - k + 1) + b);
                if (it == node_of.end())
                    continue;
                s.push_back(b);
                if (oriented_instantiates(peptide, code, o, 0, s)) {
                    ++result.candidates;
                    p.push_back(it->second);
                    grow(s, p);
                    p.pop_back();
                }
                s.pop_back();
            }
        };
        for (const auto &[kmer, node] : node_of) {
            if (!oriented_instantiates(peptide, code, o, 0, kmer))
                continue;
            ++result.anchors[o];
            std::string s = kmer;
            std::vector<node_index> p { node };
            grow(s, p);
        }
    }
    std::sort(result.paths.begin(), result.paths.end());
    return result;
}


// ---------------------------------------------------------------- the checks

std::vector<Ctx> run_contexts(const PatternSearch &engine, const StoredNodes &stored,
                              const Pattern &pattern, const Request &request, Result *result) {
    Budget budget(kManySteps, unbounded_deadline());
    std::vector<Ctx> contexts;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        EXPECT_EQ(stored.of(c.node), c.base_node) << c.node;
        EXPECT_TRUE(c.path.empty());
        contexts.push_back(Ctx { c.node, c.offset, c.orientation });
    });
    return contexts;
}

std::vector<PathCtx> run_paths(const PatternSearch &engine, const StoredNodes &stored,
                               const Pattern &pattern, const Request &request, Result *result) {
    Budget budget(kManySteps, unbounded_deadline());
    std::vector<PathCtx> paths;
    *result = engine.enumerate(pattern, request, budget, [&](const Context &c) {
        EXPECT_EQ(0u, c.offset);
        EXPECT_EQ(stored.of(c.node), c.base_node);
        for (node_index n : c.path) {
            EXPECT_EQ(stored.of(n), engine.base_node(n)) << n;
        }
        paths.push_back(PathCtx { c.node, c.orientation, c.sequence, c.path });
    });
    return paths;
}

/**
 * Every check of one (graph, peptide, table, request) against both oracles. L <= k: the
 * contexts (counts per orientation and offset, the released list in answer order, the
 * record-scan's spelled contexts, PARTIAL's prefixes, ALL_OR_COUNT's threshold). L > k: the
 * anchors, and with request.extend_paths the paths (counts, candidates_examined, the released
 * list, each path's nodes and sequence, the record-scan's occurrences: included in the
 * graph's paths, equal with |records_exact|). Returns the number of contexts or paths.
 */
size_t check_peptide(const DeBruijnGraph &graph, const std::vector<std::string> &records,
                     DeBruijnGraph::Mode mode, const std::string &peptide, int table,
                     const Request &request, bool records_exact = false,
                     const std::set<std::string> &pruned = {}) {
    const size_t k = graph.get_k();
    const size_t L = 3 * peptide.size();
    const OracleCode &code = oracle_code(table);
    const Pattern pattern = peptide_pattern(peptide, table);
    SCOPED_TRACE("peptide " + peptide + " table " + std::to_string(table) + " k "
                 + std::to_string(k) + " mode " + std::to_string(mode));
    EXPECT_EQ(L, pattern.length());
    EXPECT_EQ(oracle_palindromic(peptide, code), pattern.is_palindromic());

    PatternSearch engine(graph);
    const StoredNodes stored(graph);
    const std::vector<Orientation> orientations
        = oracle_orientations(peptide, code, request.strands);

    Budget budget(kManySteps, unbounded_deadline());
    Result counted = engine.count(pattern, request, budget);
    EXPECT_FALSE(counted.refusal);
    if (counted.refusal)
        return 0;
    EXPECT_FALSE(counted.stop);
    EXPECT_EQ(orientations, counted.searched);

    if (L <= k) {
        const std::vector<Ctx> expected = walk_contexts(graph, peptide, code, request);
        if (pruned.empty()) {
            SpelledContexts in_graph;
            for (const Ctx &c : expected) {
                in_graph.emplace(graph.get_node_sequence(c.node), c.offset, c.orientation);
            }
            EXPECT_EQ(record_contexts(records, k, mode, peptide, code, request), in_graph);
        }
        EXPECT_FALSE(counted.anchors);
        EXPECT_TRUE(counted.contexts);
        if (!counted.contexts)
            return 0;
        EXPECT_EQ(Relation::EXACT, counted.contexts->total.relation);
        EXPECT_EQ(expected.size(), counted.contexts->total.value);
        std::map<Orientation, uint64_t> by_orientation;
        std::map<uint32_t, uint64_t> by_offset;
        for (const Ctx &c : expected) {
            ++by_orientation[c.orientation];
            ++by_offset[c.offset];
        }
        EXPECT_EQ(orientations.size(), counted.contexts->by_orientation.size());
        for (const auto &[o, count] : counted.contexts->by_orientation) {
            EXPECT_EQ(Relation::EXACT, count.relation);
            EXPECT_EQ(by_orientation[o], count.value) << orientation_key(o);
        }
        for (const auto &[p, count] : counted.contexts->by_offset) {
            EXPECT_EQ(by_offset[p], count.value) << "offset " << p;
        }

        Request all = request;
        all.mode = Mode::ALL_OR_COUNT;
        Result enumerated;
        std::vector<Ctx> released = run_contexts(engine, stored, pattern, all, &enumerated);
        EXPECT_EQ(counted.work.steps, enumerated.work.steps);
        EXPECT_TRUE(enumerated.extraction && enumerated.extraction->complete);
        EXPECT_EQ(expected, released);
        for (const Ctx &c : released) {
            // the instance spells the oriented peptide, codon by codon
            const std::string instance = graph.get_node_sequence(c.node).substr(c.offset, L);
            EXPECT_TRUE(oriented_instantiates(peptide, code, c.orientation, 0, instance))
                << instance;
        }
        if (!expected.empty()) {
            Request tight = all;
            tight.max_contexts = expected.size() - 1;
            Result withheld;
            EXPECT_TRUE(run_contexts(engine, stored, pattern, tight, &withheld).empty());
            EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, withheld.extraction->withheld);
        }
        Request partial = request;
        partial.mode = Mode::PARTIAL;
        for (uint64_t cap : { uint64_t(1), uint64_t(expected.size() / 2) }) {
            partial.max_contexts = cap;
            Result cut;
            std::vector<Ctx> prefix = run_contexts(engine, stored, pattern, partial, &cut);
            EXPECT_EQ(std::min<uint64_t>(cap, expected.size()), prefix.size());
            EXPECT_TRUE(prefix.size() <= expected.size()
                            && std::equal(prefix.begin(), prefix.end(), expected.begin()));
        }
        return expected.size();
    }

    // L > k: the anchors, then the paths
    EXPECT_TRUE(counted.anchors);
    if (!counted.anchors)
        return 0;
    const std::vector<Ctx> anchors = walk_contexts(graph, peptide, code, request);
    EXPECT_EQ(Relation::EXACT, counted.anchors->total.relation);
    EXPECT_EQ(anchors.size(), counted.anchors->total.value);
    // the floor's window bits are the exact ones (§3)
    EXPECT_DOUBLE_EQ(pattern.information_bits(0, k), *counted.anchor_information_bits);

    if (!request.extend_paths) {
        EXPECT_EQ(Extension::NOT_REQUESTED, counted.anchors->extension);
        EXPECT_EQ(anchors.empty() ? Relation::EXACT : Relation::UNKNOWN,
                  counted.anchors->paths.relation);
        // the anchors themselves, released on request, in answer order
        Request release = request;
        release.release_anchors = true;
        Result r;
        EXPECT_EQ(anchors, run_contexts(engine, stored, pattern, release, &r));
        return anchors.size();
    }

    const PathOracle oracle = walk_paths(graph, peptide, code, request.strands);
    uint64_t num_anchors = 0;
    for (const auto &[o, n] : oracle.anchors) {
        num_anchors += n;
    }
    EXPECT_EQ(num_anchors, anchors.size());
    const AnchorCounts &a = *counted.anchors;
    EXPECT_EQ(num_anchors ? Extension::COMPLETED : Extension::NO_ANCHORS, a.extension);
    EXPECT_EQ(Relation::EXACT, a.paths.relation);
    EXPECT_EQ(oracle.paths.size(), a.paths.value);
    EXPECT_EQ(oracle.candidates, a.candidates_examined);
    std::map<Orientation, uint64_t> by_orientation;
    for (const PathCtx &p : oracle.paths) {
        ++by_orientation[p.orientation];
    }
    for (const auto &[o, count] : a.paths_by_orientation) {
        EXPECT_EQ(Relation::EXACT, count.relation);
        EXPECT_EQ(by_orientation[o], count.value) << orientation_key(o);
    }

    SpelledPaths in_graph;
    for (const PathCtx &p : oracle.paths) {
        in_graph.emplace(p.sequence, p.orientation);
    }
    const SpelledPaths in_records
        = record_paths(records, k, mode, peptide, code, request.strands, pruned);
    for (const auto &occurrence : in_records) {
        EXPECT_TRUE(in_graph.count(occurrence))
            << occurrence.first << " " << orientation_key(occurrence.second);
    }
    if (records_exact) {
        EXPECT_EQ(in_records, in_graph);
    }

    Request all = request;
    all.mode = Mode::ALL_OR_COUNT;
    Result enumerated;
    std::vector<PathCtx> released = run_paths(engine, stored, pattern, all, &enumerated);
    EXPECT_EQ(counted.work.steps, enumerated.work.steps);
    EXPECT_TRUE(enumerated.extraction && enumerated.extraction->complete);
    EXPECT_EQ(oracle.paths, released);
    for (const PathCtx &p : released) {
        EXPECT_EQ(L, p.sequence.size());
        EXPECT_EQ(L - k + 1, p.path.size());
        for (size_t i = 0; i < p.path.size(); ++i) {
            EXPECT_EQ(p.sequence.substr(i, k), graph.get_node_sequence(p.path[i]));
        }
        EXPECT_TRUE(oriented_instantiates(peptide, code, p.orientation, 0, p.sequence)) << p;
    }
    Request partial = request;
    partial.mode = Mode::PARTIAL;
    for (uint64_t cap : { uint64_t(1), uint64_t(oracle.paths.size() / 2) }) {
        partial.max_paths = cap;
        Result cut;
        std::vector<PathCtx> prefix = run_paths(engine, stored, pattern, partial, &cut);
        EXPECT_EQ(std::min<uint64_t>(cap, oracle.paths.size()), prefix.size());
        EXPECT_TRUE(prefix.size() <= oracle.paths.size()
                        && std::equal(prefix.begin(), prefix.end(), oracle.paths.begin()));
    }
    return oracle.paths.size();
}


// ---------------------------------------------------------------- the genetic codes

TEST(PatternPeptide, GeneticCodeTables) {
    const std::vector<int> expected_ids {
        1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33
    };
    EXPECT_EQ(expected_ids, GeneticCode::ids());
    EXPECT_EQ("1-6, 9-16, 21-33", GeneticCode::ids_text());
    for (int id : { -1, 0, 7, 8, 17, 18, 19, 20, 34, 1000 }) {
        EXPECT_EQ(nullptr, GeneticCode::find(id)) << id;
        EXPECT_THROW(GeneticCode::get(id), std::invalid_argument) << id;
    }
    EXPECT_EQ(1, GeneticCode::standard().id());
    EXPECT_STREQ("Standard", GeneticCode::standard().name());

    for (int id : expected_ids) {
        SCOPED_TRACE(id);
        const GeneticCode &table = GeneticCode::get(id);
        ASSERT_EQ(&table, GeneticCode::find(id));
        EXPECT_EQ(id, table.id());
        const OracleCode &oracle = oracle_code(id);

        // NCBI's ncbieaa string, TCAG order, against this file's table
        std::string ncbieaa;
        for (char a : std::string("TCAG")) {
            for (char b : std::string("TCAG")) {
                for (char c : std::string("TCAG")) {
                    ncbieaa.push_back(oracle.translate(std::string { a, b, c }));
                }
            }
        }
        EXPECT_EQ(ncbieaa, table.amino_acids());

        CodonSet stops = 0;
        std::map<char, CodonSet> by_residue;
        for (const std::string &codon : all_codons()) {
            const char residue = oracle.translate(codon);
            EXPECT_EQ(residue, table.translate(codon)) << codon;
            EXPECT_EQ(residue, table.translate(codon_index(codon))) << codon;
            // case-insensitive
            std::string lower = codon;
            for (char &c : lower) {
                c = std::tolower(c);
            }
            EXPECT_EQ(residue, table.translate(lower));
            const CodonSet bit = CodonSet(1) << codon_index(codon);
            if (residue == '*') {
                stops |= bit;
            } else {
                by_residue[residue] |= bit;
            }
        }
        EXPECT_EQ(stops, table.stops());
        for (char residue : std::string("ACDEFGHIKLMNPQRSTVWY")) {
            EXPECT_EQ(by_residue[residue], table.codons(residue)) << residue;
            EXPECT_EQ(by_residue[residue], table.codons(std::tolower(residue))) << residue;
            // every residue has a codon in every table
            EXPECT_NE(0u, table.codons(residue)) << residue;
        }
        for (char other : std::string("*XBZJUO-1 ")) {
            EXPECT_EQ(0u, table.codons(other)) << other;
        }
        EXPECT_EQ(0, table.translate("ACN"));
        EXPECT_EQ(0, table.translate("AC"));
        EXPECT_EQ(0, table.translate(64));
    }

    // documented differences, by name
    const GeneticCode &standard = GeneticCode::get(1);
    EXPECT_EQ('*', standard.translate("TGA"));
    EXPECT_EQ('I', standard.translate("ATA"));
    EXPECT_EQ('R', standard.translate("AGA"));
    const GeneticCode &vertebrate_mito = GeneticCode::get(2);
    EXPECT_EQ('*', vertebrate_mito.translate("AGA"));
    EXPECT_EQ('*', vertebrate_mito.translate("AGG"));
    EXPECT_EQ('M', vertebrate_mito.translate("ATA"));
    EXPECT_EQ('W', vertebrate_mito.translate("TGA"));
    EXPECT_EQ('W', GeneticCode::get(4).translate("TGA"));
    EXPECT_EQ('Q', GeneticCode::get(6).translate("TAA"));
    EXPECT_EQ('S', GeneticCode::get(12).translate("CTG"));
    EXPECT_EQ('*', GeneticCode::get(23).translate("TTA"));
    // table 11 (bacteria, archaea, plastids) codes as the standard code
    EXPECT_EQ(standard.amino_acids(), GeneticCode::get(11).amino_acids());
    // the codons that are a residue and also a terminator in context
    for (int id : expected_ids) {
        CodonSet expected = 0;
        if (id == 27)
            expected = CodonSet(1) << codon_index("TGA");
        if (id == 28) {
            expected = CodonSet(1) << codon_index("TAA") | CodonSet(1) << codon_index("TAG")
                        | CodonSet(1) << codon_index("TGA");
        }
        if (id == 31)
            expected = CodonSet(1) << codon_index("TAA") | CodonSet(1) << codon_index("TAG");
        EXPECT_EQ(expected, GeneticCode::get(id).context_stops()) << id;
        EXPECT_EQ(0u, GeneticCode::get(id).context_stops() & GeneticCode::get(id).stops());
    }
}

TEST(PatternPeptide, CodonHelpers) {
    for (const std::string &codon : all_codons()) {
        const int c = codon_index(codon);
        ASSERT_LE(0, c);
        ASSERT_GT(64, c);
        EXPECT_EQ(codon, codon_string(c));
        EXPECT_EQ(CodonSet(1) << codon_index(rev_comp(codon)),
                  reverse_complement_codons(CodonSet(1) << c));
    }
    EXPECT_EQ(-1, codon_index("ACN"));
    EXPECT_EQ(-1, codon_index("ACGT"));
    EXPECT_EQ(codon_index("ACG"), codon_index("acg"));
    EXPECT_EQ(kAllCodons, reverse_complement_codons(kAllCodons));

    std::mt19937 rng(7);
    for (int t = 0; t < 200; ++t) {
        CodonSet set = (CodonSet(rng()) << 32 | rng()) & (CodonSet(rng()) << 32 | rng());
        // next_bases against the strings
        for (size_t j = 0; j < 3; ++j) {
            for (int b0 = -1; b0 < 4; ++b0) {
                for (int b1 = -1; b1 < 4; ++b1) {
                    uint8_t expected = 0;
                    for (const std::string &codon : all_codons()) {
                        if (!(set >> codon_index(codon) & 1))
                            continue;
                        if (j > 0 && b0 >= 0 && codon[0] != "ACGT"[b0])
                            continue;
                        if (j > 1 && b1 >= 0 && codon[1] != "ACGT"[b1])
                            continue;
                        expected |= 1 << std::string("ACGT").find(codon[j]);
                    }
                    EXPECT_EQ(expected, next_bases(set, j, b0, b1));
                }
            }
        }
        // cylinder and distinct_projections against the strings
        for (size_t begin = 0; begin <= 3; ++begin) {
            for (size_t end = begin; end <= 3; ++end) {
                std::set<std::string> pieces;
                for (const std::string &codon : all_codons()) {
                    if (set >> codon_index(codon) & 1)
                        pieces.insert(codon.substr(begin, end - begin));
                }
                EXPECT_EQ(pieces.size(), distinct_projections(set, begin, end));
                CodonSet expected = 0;
                for (const std::string &codon : all_codons()) {
                    if (pieces.count(codon.substr(begin, end - begin)))
                        expected |= CodonSet(1) << codon_index(codon);
                }
                EXPECT_EQ(expected, cylinder(set, begin, end));
            }
        }
    }
}


// ---------------------------------------------------------------- parsing and the automaton

TEST(PatternPeptide, ParsePeptides) {
    Pattern p = Pattern::parse(PatternKind::PROTEIN, "mKvLaXbzj");
    EXPECT_EQ(PatternKind::PROTEIN, p.kind());
    EXPECT_STREQ("protein", to_string(p.kind()));
    EXPECT_EQ("MKVLAXBZJ", p.text());
    EXPECT_EQ(27u, p.length());
    EXPECT_EQ(9u, p.codon_sets().size());
    EXPECT_EQ(1, p.genetic_code());
    EXPECT_FALSE(p.is_exact());
    // the default code is the standard one
    const Pattern standard = Pattern::parse(PatternKind::PROTEIN, "MKVLAXBZJ",
                                            GeneticCode::get(1));
    EXPECT_EQ(standard.codon_sets(), p.codon_sets());
    EXPECT_EQ(standard.positions(), p.positions());
    // the code is read: ATA is Met in table 2
    const Pattern mito = Pattern::parse(PatternKind::PROTEIN, "M", GeneticCode::get(2));
    EXPECT_EQ(2, mito.genetic_code());
    EXPECT_EQ(CodonSet(1) << codon_index("ATG") | CodonSet(1) << codon_index("ATA"),
              mito.codon_sets()[0]);
    // the code is not read for DNA and IUPAC
    const Pattern dna = Pattern::parse(PatternKind::DNA, "ACGT", GeneticCode::get(2));
    EXPECT_EQ(0, dna.genetic_code());
    EXPECT_TRUE(dna.codon_sets().empty());
    EXPECT_EQ(Pattern::parse(PatternKind::DNA, "ACGT").positions(), dna.positions());
    // M and W have one codon each in the standard code: exact
    EXPECT_TRUE(Pattern::parse(PatternKind::PROTEIN, "MWMW").is_exact());
    EXPECT_FALSE(Pattern::parse(PatternKind::PROTEIN, "MWMW", GeneticCode::get(2)).is_exact());

    // the stop '*' is a residue (owner decision #19 of 2026-10-08; the refusal
    // stop_unsupported of 4596bb3b is gone): the table's stop codons at that position
    for (const char *text : { "M*K", "*", "MK**w*" }) {
        const Pattern p = Pattern::parse(PatternKind::PROTEIN, text);
        EXPECT_EQ(3 * std::string(text).size(), p.length()) << text;
        EXPECT_TRUE(p.has_instances()) << text;
    }
    const Pattern stop = Pattern::parse(PatternKind::PROTEIN, "M*K");
    EXPECT_EQ("M*K", stop.text());
    EXPECT_EQ(CodonSet(1) << codon_index("TAA") | CodonSet(1) << codon_index("TAG")
                  | CodonSet(1) << codon_index("TGA"),
              stop.codon_sets()[1]);
    EXPECT_EQ(GeneticCode::standard().stops(), stop.codon_sets()[1]);
    // a table without an unconditional stop (27, 28, 31): '*' admits no codon, the peptide
    // has no instance (answered EXACT 0 with a note), its bits stay finite
    for (int table : { 27, 28, 31 }) {
        const Pattern none = Pattern::parse(PatternKind::PROTEIN, "M*K", GeneticCode::get(table));
        EXPECT_FALSE(none.has_instances()) << table;
        EXPECT_EQ(0u, none.codon_sets()[1]) << table;
        EXPECT_TRUE(std::isfinite(none.information_bits())) << table;
        EXPECT_FALSE(none.reverse_complement().has_instances()) << table;
    }
    EXPECT_TRUE(Pattern::parse(PatternKind::PROTEIN, "MXK", GeneticCode::get(27)).has_instances());
    EXPECT_TRUE(Pattern::parse(PatternKind::DNA, "ACGT").has_instances());

    // refusals: the residues not served, every other character (bad_alphabet, named
    // wherever a stop is), the empty text
    struct Bad {
        std::string text;
        size_t position;
        std::string says;
        std::string code = "bad_alphabet";
    };
    for (const Bad &bad : std::vector<Bad> {
            { "*MU", 2, "protein alphabet" }, { "M*K-", 3, "protein alphabet" },
            { "MUK", 1, "protein alphabet" },
            { "MKO", 2, "protein alphabet" }, { "M K", 1, "protein alphabet" },
            { "M-K", 1, "protein alphabet" }, { "M1", 1, "protein alphabet" },
            { "MK.", 2, "protein alphabet" } }) {
        try {
            Pattern::parse(PatternKind::PROTEIN, bad.text);
            ADD_FAILURE() << bad.text << " parsed";
        } catch (const PatternError &e) {
            EXPECT_EQ(bad.code, e.code()) << bad.text;
            const std::string message = e.what();
            EXPECT_NE(std::string::npos,
                      message.find("at position " + std::to_string(bad.position)))
                << message;
            EXPECT_NE(std::string::npos, message.find(bad.says)) << message;
        }
    }
    EXPECT_THROW(Pattern::parse(PatternKind::PROTEIN, ""), PatternError);
}

// the bases the oracle's codons of |letter| allow after |prefix| (the codon's bases so far)
uint8_t oracle_next(char letter, const OracleCode &code, const std::string &prefix) {
    uint8_t result = 0;
    for (const std::string &codon : codons_of(letter, code)) {
        if (codon.compare(0, prefix.size(), prefix) == 0)
            result |= 1 << std::string("ACGT").find(codon[prefix.size()]);
    }
    return result;
}

TEST(PatternPeptide, AutomatonIsExact) {
    // every table, every letter, every codon prefix: allowed() is exactly the bases that
    // extend the prefix to a codon of the residue (no superset), never into a stop
    const std::string letters = "ACDEFGHIKLMNPQRSTVWYXBZJ";
    for (int id : GeneticCode::ids()) {
        const OracleCode &code = oracle_code(id);
        for (char letter : letters) {
            const Pattern p = peptide_pattern(std::string("M") + letter, id);
            for (const std::string &codon : all_codons()) {
                for (size_t j = 0; j < 3; ++j) {
                    const std::string prefix = codon.substr(0, j);
                    // what precedes the codon is M's ATG, as the extension would pass it
                    EXPECT_EQ(oracle_next(letter, code, prefix),
                              p.allowed(3 + j, "ATG" + prefix))
                        << id << " " << letter << " " << prefix;
                }
            }
            // a symbol other than A, C, G, T in the codon admits nothing
            EXPECT_EQ(0, p.allowed(4, "ATGN"));
            EXPECT_EQ(0, p.allowed(5, "ATGAN"));
        }
    }

    // the cases a per-position IUPAC reading gets wrong, in the standard code
    const Pattern lsr = peptide_pattern("LSR", 1);
    EXPECT_EQ(kBaseC | kBaseT, lsr.allowed(0, ""));
    EXPECT_EQ(kBaseT, lsr.allowed(1, "C"));
    EXPECT_EQ(kBaseT, lsr.allowed(1, "T"));
    EXPECT_EQ(kAllBases, lsr.allowed(2, "CT"));
    EXPECT_EQ(kBaseA | kBaseG, lsr.allowed(2, "TT"));  // TTA, TTG; TTC and TTT are Phe
    EXPECT_EQ(kBaseA | kBaseT, lsr.allowed(3, "CTG"));
    EXPECT_EQ(kBaseG, lsr.allowed(4, "CTGA"));       // AGT, AGC; ACN is Thr
    EXPECT_EQ(kBaseC, lsr.allowed(4, "CTGT"));       // TCN; TGN is Cys, Trp, stop
    EXPECT_EQ(kBaseC | kBaseT, lsr.allowed(5, "CTGAG"));
    EXPECT_EQ(kAllBases, lsr.allowed(5, "CTGTC"));
    EXPECT_EQ(kBaseA | kBaseC, lsr.allowed(6, "CTGAGC"));
    EXPECT_EQ(kBaseG, lsr.allowed(7, "CTGAGCA"));    // AGA, AGG; AGT and AGC are Ser
    EXPECT_EQ(kBaseA | kBaseG, lsr.allowed(8, "CTGAGCAG"));
    EXPECT_EQ(kAllBases, lsr.allowed(8, "CTGAGCCG"));
    // X never enters a stop: TAA, TAG (standard) and AGA, AGG (vertebrate mitochondrial)
    EXPECT_EQ(kBaseC | kBaseT, peptide_pattern("X", 1).allowed(2, "TA"));
    EXPECT_EQ(kBaseC | kBaseG | kBaseT, peptide_pattern("X", 1).allowed(2, "TG"));
    EXPECT_EQ(kBaseC | kBaseT, peptide_pattern("X", 2).allowed(2, "AG"));
    EXPECT_EQ(kAllBases, peptide_pattern("X", 2).allowed(2, "TG"));
    // fewer spelled bases than the codon position: the missing ones are any base
    EXPECT_EQ(kBaseT, lsr.allowed(1, ""));
    EXPECT_EQ(kAllBases, lsr.allowed(2, "T") | lsr.allowed(2, "C"));
}

// the instances of a pattern as allowed() admits them, by walking the automaton
std::set<std::string> automaton_language(const Pattern &p) {
    std::set<std::string> result;
    std::function<void(std::string&)> grow = [&](std::string &s) {
        if (s.size() == p.length()) {
            result.insert(s);
            return;
        }
        const BaseSet allowed = p.allowed(s.size(), s);
        for (size_t b = 0; b < 4; ++b) {
            if (allowed & (1 << b)) {
                s.push_back("ACGT"[b]);
                grow(s);
                s.pop_back();
            }
        }
    };
    std::string s;
    grow(s);
    return result;
}

TEST(PatternPeptide, ReverseComplementAutomaton) {
    // rc(P)'s instances are exactly the reverse complements of P's; P's are the codon strings
    // of its residues
    std::mt19937 rng(11);
    const std::string letters = "ACDEFGHIKLMNPQRSTVWYXBZJ";
    for (int t = 0; t < 60; ++t) {
        const int id = GeneticCode::ids()[rng() % GeneticCode::ids().size()];
        const OracleCode &code = oracle_code(id);
        std::string peptide;
        for (size_t i = 0, m = 1 + rng() % 2; i < m; ++i) {
            peptide.push_back(letters[rng() % letters.size()]);
        }
        const Pattern p = peptide_pattern(peptide, id);
        std::set<std::string> expected { "" };
        for (char letter : peptide) {
            std::set<std::string> longer;
            for (const std::string &s : expected) {
                for (const std::string &codon : codons_of(letter, code)) {
                    longer.insert(s + codon);
                }
            }
            expected = std::move(longer);
        }
        EXPECT_EQ(expected, automaton_language(p)) << peptide << " " << id;
        std::set<std::string> expected_rc;
        for (const std::string &s : expected) {
            expected_rc.insert(rev_comp(s));
        }
        const Pattern rc = p.reverse_complement();
        EXPECT_EQ(PatternKind::PROTEIN, rc.kind());
        EXPECT_EQ(p.length(), rc.length());
        EXPECT_EQ(expected_rc, automaton_language(rc)) << peptide << " " << id;
        EXPECT_EQ(p.codon_sets(), rc.reverse_complement().codon_sets());
        // positions(): the per-position union, a superset of every instance's base
        for (const std::string &s : expected) {
            for (size_t i = 0; i < s.size(); ++i) {
                EXPECT_TRUE(p.positions()[i] & (1 << std::string("ACGT").find(s[i])));
            }
        }
        EXPECT_EQ(expected == expected_rc, p.is_palindromic()) << peptide << " " << id;
    }
    // the palindromic peptides: X runs in the tables without a stop
    for (int id : GeneticCode::ids()) {
        const bool no_stop = id == 27 || id == 28 || id == 31;
        EXPECT_EQ(no_stop, peptide_pattern("XXX", id).is_palindromic()) << id;
        EXPECT_FALSE(peptide_pattern("XMX", id).is_palindromic()) << id;
    }
}

TEST(PatternPeptide, InformationBitsAreExact) {
    // whole residues: log2(64 / codons)
    EXPECT_DOUBLE_EQ(6, peptide_pattern("M", 1).information_bits());
    EXPECT_DOUBLE_EQ(6 - std::log2(6.0), peptide_pattern("L", 1).information_bits());
    EXPECT_DOUBLE_EQ(6 - std::log2(61.0), peptide_pattern("X", 1).information_bits());
    // four stops in the vertebrate mitochondrial code: TAA, TAG, AGA, AGG
    EXPECT_DOUBLE_EQ(6 - std::log2(60.0), peptide_pattern("X", 2).information_bits());
    EXPECT_DOUBLE_EQ(0, peptide_pattern("XXX", 27).information_bits());
    EXPECT_DOUBLE_EQ(6 - std::log2(4.0), peptide_pattern("B", 1).information_bits());
    // Ser is 6 codons (TCN, AGY): 3.42 bits, where a per-position reading ({A,T}, {C,G}, N)
    // would say 2
    EXPECT_DOUBLE_EQ(6 - std::log2(6.0), peptide_pattern("S", 1).information_bits());

    // every window [begin, end) of random peptides: 2 (end - begin) - log2(the distinct
    // strings the instances spell there), the strings enumerated
    std::mt19937 rng(13);
    const std::string letters = "ACDEFGHIKLMNPQRSTVWYXBZJ";
    for (int t = 0; t < 40; ++t) {
        const int id = GeneticCode::ids()[rng() % GeneticCode::ids().size()];
        std::string peptide;
        for (size_t i = 0, m = 1 + rng() % 3; i < m; ++i) {
            peptide.push_back(letters[rng() % letters.size()]);
        }
        const Pattern p = peptide_pattern(peptide, id);
        const std::set<std::string> language = automaton_language(p);
        for (size_t begin = 0; begin <= p.length(); ++begin) {
            for (size_t end = begin; end <= p.length(); ++end) {
                std::set<std::string> pieces;
                for (const std::string &s : language) {
                    pieces.insert(s.substr(begin, end - begin));
                }
                EXPECT_NEAR(2.0 * (end - begin) - std::log2(double(pieces.size())),
                            p.information_bits(begin, end), 1e-9)
                    << peptide << " " << id << " [" << begin << ", " << end << ")";
            }
        }
    }
}

TEST(PatternPeptide, InformationFloor) {
    std::vector<std::string> records { "ATGCTGAGCAGATGGATGCCATGGCCAAGT" };
    auto graph = build(21, records, DeBruijnGraph::BASIC);
    Request request = make_request();
    request.min_information_bits = 24;
    // LLLL: 4 * (6 - log2 6) = 13.7 bits, below the floor
    Result low = count_of(*graph, peptide_pattern("LLLL", 1), request);
    ASSERT_TRUE(low.refusal);
    EXPECT_EQ("information_below_floor", low.refusal->code);
    EXPECT_EQ(0u, low.work.steps);
    // MLSRW: 6 + 3 * 3.42 + 6 = 22.2 bits; with one more Met, 28.2
    EXPECT_TRUE(count_of(*graph, peptide_pattern("MLSRW", 1), request).refusal);
    EXPECT_FALSE(count_of(*graph, peptide_pattern("MLSRWM", 1), request).refusal);
    // an exact peptide in suffix scope is admitted however short (one range)
    Request suffix = request;
    suffix.scope = Scope::SUFFIX;
    EXPECT_FALSE(count_of(*graph, peptide_pattern("MW", 1), suffix).refusal);
    EXPECT_TRUE(count_of(*graph, peptide_pattern("MW", 1), request).refusal);
    EXPECT_TRUE(count_of(*graph, peptide_pattern("ML", 1), suffix).refusal);

    // a long peptide: the floor reads the exact bits of every searched anchor window
    auto small = build(5, records, DeBruijnGraph::BASIC);
    Request floor = make_request(Scope::ANY_OFFSET, Strands::BOTH, true);
    floor.min_information_bits = 9;
    // MLSRW at k = 5: P[0, 5) = ATG + Leu's first two (CT, TT): 6 + 4 - 1 = 9 bits; the
    // reverse window rc(P[10, 15)) = Arg's last two (GT GC GA GG) and Trp: 4 - 2 + 6 = 8
    Result r = count_of(*small, peptide_pattern("MLSRW", 1), floor);
    ASSERT_TRUE(r.refusal);
    EXPECT_NE(std::string::npos, r.refusal->message.find("reverse orientation"));
    floor.strands = Strands::FORWARD;
    r = count_of(*small, peptide_pattern("MLSRW", 1), floor);
    EXPECT_FALSE(r.refusal);
    EXPECT_DOUBLE_EQ(9, *r.anchor_information_bits);
    EXPECT_DOUBLE_EQ(9, *r.min_anchor_information_bits);
}


// ---------------------------------------------------------------- the engine against the oracles

TEST(PatternPeptide, LeuSerArgNoSuperset) {
    // LSR: Leu TTR|CTN, Ser TCN|AGY, Arg CGN|AGR. Each decoy record changes one codon to one
    // a per-position reading of the codon sets would admit (TTC Phe, ACC Thr, AGC Ser)
    const std::vector<std::string> records {
        "GG" "CTGAGCAGA" "GG",  // L S R
        "GG" "TTCAGCAGA" "GG",  // F S R
        "GG" "CTGACCAGA" "GG",  // L T R
        "GG" "CTGAGCAGC" "GG",  // L S S
        "GG" "TTAAGTCGG" "GG",  // L S R
        "GG" "CTCTCAAGG" "GG",  // L S R
    };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(11, records, mode);
        for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
            check_peptide(*graph, records, mode, "LSR", 1, make_request(Scope::ANY_OFFSET,
                                                                         strands));
            if (mode != DeBruijnGraph::PRIMARY)
                check_peptide(*graph, records, mode, "LSR", 1, make_request(Scope::SUFFIX,
                                                                             strands));
        }
    }
    // on the + strand: three records with LSR, three k-mers each
    auto graph = build(11, records, DeBruijnGraph::BASIC);
    Result r = count_of(*graph, peptide_pattern("LSR", 1),
                        make_request(Scope::ANY_OFFSET, Strands::FORWARD));
    EXPECT_EQ(9u, r.contexts->total.value);
    // in suffix scope one per record
    r = count_of(*graph, peptide_pattern("LSR", 1), make_request(Scope::SUFFIX, Strands::FORWARD));
    EXPECT_EQ(3u, r.contexts->total.value);
    // the decoys of the IUPAC reading YTN WSN MGN are found by it, not by the peptide
    r = count_of(*graph, Pattern::parse(PatternKind::IUPAC, "YTNWSNMGN"),
                 make_request(Scope::SUFFIX, Strands::FORWARD));
    EXPECT_EQ(6u, r.contexts->total.value);
}

TEST(PatternPeptide, AmbiguityCodes) {
    // X any residue, B Asp or Asn, Z Glu or Gln, J Ile or Leu
    const std::vector<std::string> records {
        "C" "ATGGATGAAATTTGG" "CC",   // M D E I W
        "C" "ATGAACCAGCTGTGG" "CC",   // M N Q L W
        "C" "ATGGAGCAAATATGG" "AA",   // M E Q I W
        "C" "ATGCATAAACTTTGG" "TT",   // M H K L W
        "C" "ATGTAAGAAATTTGG" "CC",   // M * E I W
    };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
        auto graph = build(16, records, mode);
        for (std::string peptide : { "MBZJW", "MXZJW", "MBXJW", "MBZXW", "MXXXW", "BZJ", "MXQ" }) {
            check_peptide(*graph, records, mode, peptide, 1, make_request());
        }
    }
    // MBZJW: records 1 and 2 (D E I, N Q L); MXZJW adds record 3 (E Q I); MXXXW adds record
    // 4 (H K L) and never the stop of record 5. One context per record in suffix scope, two
    // (offsets 0 and 1) in any_offset
    auto graph = build(16, records, DeBruijnGraph::BASIC);
    auto forward = make_request(Scope::SUFFIX, Strands::FORWARD);
    EXPECT_EQ(2u, count_of(*graph, peptide_pattern("MBZJW", 1), forward).contexts->total.value);
    EXPECT_EQ(3u, count_of(*graph, peptide_pattern("MXZJW", 1), forward).contexts->total.value);
    EXPECT_EQ(4u, count_of(*graph, peptide_pattern("MXXXW", 1), forward).contexts->total.value);
    forward.scope = Scope::ANY_OFFSET;
    EXPECT_EQ(4u, count_of(*graph, peptide_pattern("MBZJW", 1), forward).contexts->total.value);
    EXPECT_EQ(6u, count_of(*graph, peptide_pattern("MXZJW", 1), forward).contexts->total.value);
    EXPECT_EQ(8u, count_of(*graph, peptide_pattern("MXXXW", 1), forward).contexts->total.value);
}

TEST(PatternPeptide, OnlyInstanceCrossesAStop) {
    // WXW: the record's only candidate is TGG TAG TGG, whose middle codon is a stop in the
    // standard code: no hit, although every base of it is one X's codons admit at its
    // position. In tables 6 and 15 TAG is Gln, in 16 and 22 Leu: one hit
    const std::vector<std::string> records { "AA" "TGGTAGTGG" "AA" };
    auto graph = build(11, records, DeBruijnGraph::BASIC);
    auto forward = make_request(Scope::ANY_OFFSET, Strands::FORWARD);
    EXPECT_EQ(0u, count_of(*graph, peptide_pattern("WXW", 1), forward).contexts->total.value);
    EXPECT_EQ(Relation::EXACT,
              count_of(*graph, peptide_pattern("WXW", 1), forward).contexts->total.relation);
    for (int id : { 6, 15, 16, 22 }) {
        EXPECT_EQ(3u, count_of(*graph, peptide_pattern("WXW", id), forward)
                          .contexts->total.value) << id;
    }
    // the IUPAC reading of X (NNN) finds it
    EXPECT_EQ(3u, count_of(*graph, Pattern::parse(PatternKind::IUPAC, "TGGNNNTGG"), forward)
                      .contexts->total.value);
    for (int id : { 1, 2, 6, 15, 16, 22 }) {
        check_peptide(*graph, records, DeBruijnGraph::BASIC, "WXW", id, make_request());
    }

    // longer than k: the anchor TGGTA exists (TA starts Tyr's TAC and TAT) and its only path
    // runs through the stop (the other anchor, TGGAA, ends the record): anchors 2, paths 0 in
    // the standard code, 1 in table 6
    auto small = build(5, records, DeBruijnGraph::BASIC);
    auto paths = make_request(Scope::ANY_OFFSET, Strands::FORWARD, true);
    Result r = count_of(*small, peptide_pattern("WXW", 1), paths);
    EXPECT_EQ(2u, r.anchors->total.value);
    EXPECT_EQ(Extension::COMPLETED, r.anchors->extension);
    EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
    EXPECT_EQ(0u, r.anchors->paths.value);
    r = count_of(*small, peptide_pattern("WXW", 6), paths);
    EXPECT_EQ(1u, r.anchors->paths.value);
    for (int id : { 1, 6 }) {
        for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                           DeBruijnGraph::PRIMARY }) {
            auto g = build(5, records, mode);
            check_peptide(*g, records, mode, "WXW", id, make_request(Scope::ANY_OFFSET,
                                                                     Strands::BOTH, true),
                          true);
        }
    }
}

TEST(PatternPeptide, NonStandardTableChangesAHit) {
    // vertebrate mitochondrial (2): ATA is Met, AGA and AGG are stops
    const std::vector<std::string> records { "CC" "ATATGG" "CC", "GG" "CGTAGA" "GG" };
    auto graph = build(8, records, DeBruijnGraph::BASIC);
    auto forward = make_request(Scope::SUFFIX, Strands::FORWARD);
    EXPECT_EQ(0u, count_of(*graph, peptide_pattern("MW", 1), forward).contexts->total.value);
    EXPECT_EQ(1u, count_of(*graph, peptide_pattern("MW", 2), forward).contexts->total.value);
    EXPECT_EQ(1u, count_of(*graph, peptide_pattern("RR", 1), forward).contexts->total.value);
    EXPECT_EQ(0u, count_of(*graph, peptide_pattern("RR", 2), forward).contexts->total.value);
    // ... and in table 5 AGA is Ser
    EXPECT_EQ(1u, count_of(*graph, peptide_pattern("RS", 5), forward).contexts->total.value);
    for (int id : { 1, 2, 5 }) {
        for (std::string peptide : { "MW", "RR", "RS", "JW", "XX" }) {
            check_peptide(*graph, records, DeBruijnGraph::BASIC, peptide, id, make_request());
            check_peptide(*graph, records, DeBruijnGraph::BASIC, peptide, id,
                          make_request(Scope::SUFFIX));
        }
    }
}

TEST(PatternPeptide, BothStrands) {
    // MKW = ATG AAR TGG, deposited on the - strand only
    const std::string coding = "ATGAAATGG";
    const std::vector<std::string> records { "CC" + rev_comp(coding) + "CC" };
    auto graph = build(10, records, DeBruijnGraph::BASIC);
    const Pattern mkw = peptide_pattern("MKW", 1);
    EXPECT_EQ(0u, count_of(*graph, mkw, make_request(Scope::ANY_OFFSET, Strands::FORWARD))
                      .contexts->total.value);
    Result r = count_of(*graph, mkw, make_request(Scope::ANY_OFFSET, Strands::REVERSE));
    EXPECT_EQ(2u, r.contexts->total.value);
    EXPECT_EQ(2u, r.contexts->by_orientation.at(Orientation::REVERSE).value);
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto g = build(10, records, mode);
        for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
            check_peptide(*g, records, mode, "MKW", 1, make_request(Scope::ANY_OFFSET, strands));
        }
        // longer than k, mid-codon cuts on both orientations' windows: on BASIC the path
        // of rc(P) only; where both orientations of every k-mer are k-mers, one of each
        const bool basic = mode == DeBruijnGraph::BASIC;
        for (size_t k : { 4, 5, 7, 8 }) {
            auto small = build(k, records, mode);
            for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
                const size_t expected = strands == Strands::BOTH ? (basic ? 1 : 2)
                                      : strands == Strands::FORWARD ? (basic ? 0 : 1)
                                      : 1;
                EXPECT_EQ(expected,
                          check_peptide(*small, records, mode, "MKW", 1,
                                        make_request(Scope::ANY_OFFSET, strands, true), basic));
                check_peptide(*small, records, mode, "MKW", 1,
                              make_request(Scope::ANY_OFFSET, strands, false));
            }
        }
    }
}

TEST(PatternPeptide, AutomatonStateAcrossTheKBoundary) {
    // MLW at k = 5: the anchor window ATG + Leu's first two bases ends inside the Leu codon,
    // so the extension's next base depends on the anchor's spelled bases: after TT only A or
    // G (TTC is Phe), after CT any. Records: TTA (Leu), TTC (Phe), CTC (Leu)
    const std::vector<std::string> records { "ATGTTATGG", "ATGTTCTGG", "ATGCTCTGG" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    auto forward = make_request(Scope::ANY_OFFSET, Strands::FORWARD, true);
    Result r = count_of(*graph, peptide_pattern("MLW", 1), forward);
    // the anchors ATGTT (records 1 and 2 share it) and ATGCT
    EXPECT_EQ(2u, r.anchors->total.value);
    EXPECT_EQ(2u, r.anchors->paths.value);
    std::vector<std::string> sequences;
    {
        Budget budget(kManySteps, unbounded_deadline());
        Request all = forward;
        PatternSearch(*graph).enumerate(peptide_pattern("MLW", 1), all, budget,
                                        [&](const Context &c) {
            sequences.push_back(c.sequence);
        });
    }
    std::sort(sequences.begin(), sequences.end());
    EXPECT_EQ(std::vector<std::string>({ "ATGCTCTGG", "ATGTTATGG" }), sequences);
    // the per-position reading ATG YTN TGG extends ATGTT through C too
    r = count_of(*graph, Pattern::parse(PatternKind::IUPAC, "ATGYTNTGG"), forward);
    EXPECT_EQ(3u, r.anchors->paths.value);

    for (size_t k : { 3, 4, 5, 6, 7, 8 }) {
        for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                           DeBruijnGraph::PRIMARY }) {
            auto g = build(k, records, mode);
            for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
                check_peptide(*g, records, mode, "MLW", 1,
                              make_request(Scope::ANY_OFFSET, strands, true),
                              k >= 5 && mode == DeBruijnGraph::BASIC);
            }
        }
    }
}

TEST(PatternPeptide, ShortAndLongPeptides) {
    // the same records, k from below one codon to above the peptide: one k-mer (3m <= k) or
    // the extension (3m > k), with the window cut at every codon position
    const std::vector<std::string> records {
        "TTGATGCTGAGCAGATGGAAACATTAGGG",
        "CCCATGTTAAGTCGGTGGAAGCAC",
        "ATGCTGAGCAGGTGTAAACATGGC",
    };
    for (size_t k : { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19 }) {
        for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                           DeBruijnGraph::PRIMARY }) {
            auto graph = build(k, records, mode);
            for (std::string peptide : { "MLSRWKH", "MLSR", "SRWK", "MJBR", "MLXR", "LSRW" }) {
                const bool long_peptide = 3 * peptide.size() > k;
                check_peptide(*graph, records, mode, peptide, 1,
                              make_request(Scope::ANY_OFFSET, Strands::BOTH, long_peptide));
                if (long_peptide) {
                    check_peptide(*graph, records, mode, peptide, 1,
                                  make_request(Scope::ANY_OFFSET, Strands::BOTH, false));
                }
            }
        }
    }
    // the 7-residue peptide is in the first two records only with M L S R W K H; in the
    // first with its codons ATG CTG AGC AGA TGG AAA CAT
    auto graph = build(31, records, DeBruijnGraph::BASIC);
    Result r = count_of(*graph, peptide_pattern("MLSRWKH", 1),
                        make_request(Scope::SUFFIX, Strands::FORWARD));
    EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
}

TEST(PatternPeptide, PrimaryReverseWindowStartsMidCodon) {
    // on a wrapped PRIMARY graph a long peptide's anchors are found as Q[0, k) and, mapped,
    // as rc(Q[0, k)) = rc(Q)[L - k, L): with L = 9 and k = 7 that window starts at the third
    // base of a codon of the reverse automaton
    const std::vector<std::string> records {
        "ATGCTGTGG", rev_comp("ATGTTGTGGAA"), "GGATGCTATGGC", rev_comp("ATGTTCTGG")
    };
    for (size_t k : { 4, 5, 7, 8 }) {
        auto graph = build(k, records, DeBruijnGraph::PRIMARY);
        for (Strands strands : { Strands::BOTH, Strands::FORWARD, Strands::REVERSE }) {
            EXPECT_EQ(strands == Strands::BOTH ? 6u : 3u,
                      check_peptide(*graph, records, DeBruijnGraph::PRIMARY, "MLW", 1,
                                    make_request(Scope::ANY_OFFSET, strands, true), true));
        }
    }
}

TEST(PatternPeptide, PalindromicPeptide) {
    // in the Karyorelict code (27) no codon is a stop, so XXX admits every 9-mer, as does its
    // reverse complement: searched once, its contexts "palindromic" (strand "=")
    const std::vector<std::string> records { "ACGTTGCAAGGCTTACGG", "TTTACGGATCAAC" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto graph = build(10, records, mode);
        Result r = count_of(*graph, peptide_pattern("XXX", 27), make_request());
        EXPECT_TRUE(r.palindromic);
        EXPECT_EQ(std::vector<Orientation>({ Orientation::PALINDROMIC }), r.searched);
        check_peptide(*graph, records, mode, "XXX", 27, make_request());
        check_peptide(*graph, records, mode, "XXX", 1, make_request());
        auto small = build(5, records, mode);
        check_peptide(*small, records, mode, "XXX", 27,
                      make_request(Scope::ANY_OFFSET, Strands::BOTH, true));
    }
}

TEST(PatternPeptide, LowComplexityNoteOnTheCodons) {
    // an exact peptide is noted low-complexity as the DNA pattern of its codons is (sdust on
    // ATG ATG ..., not on the letters); a peptide that is not exact never is
    std::vector<std::string> records { "CCATGATGATGATGATGATGATGCCATGTGGATGTGGTGGACC" };
    auto graph = build(31, records, DeBruijnGraph::BASIC);
    auto noted = [](const Result &r) {
        return std::find(r.notes.begin(), r.notes.end(), std::string(kNoteLowComplexity))
                != r.notes.end();
    };
    std::string mmm;
    for (int i = 0; i < 30; ++i) {
        mmm += "ATG";
    }
    const bool dna_noted = noted(count_of(*graph, Pattern::parse(PatternKind::DNA, mmm),
                                          make_request()));
    EXPECT_TRUE(dna_noted);
    EXPECT_EQ(dna_noted, noted(count_of(*graph, peptide_pattern(std::string(30, 'M'), 1),
                                        make_request())));
    for (std::string peptide : { "MWMWMWMWMW", "MWWMWMMWMMMW", "WWWWWWWWWWWWWWWWWWWW" }) {
        std::string bases;
        for (char residue : peptide) {
            bases += residue == 'M' ? "ATG" : "TGG";
        }
        EXPECT_EQ(noted(count_of(*graph, Pattern::parse(PatternKind::DNA, bases),
                                 make_request())),
                  noted(count_of(*graph, peptide_pattern(peptide, 1), make_request())))
            << peptide;
    }
    EXPECT_FALSE(noted(count_of(*graph, peptide_pattern(std::string(30, 'K'), 1),
                                make_request())));
}

TEST(PatternPeptide, PrunedMiddleKmer) {
    // a path exists only where all its k-mers were retained (§3, "Covered sequence")
    const std::vector<std::string> records { "ATGCTGTGGAAA" };
    auto graph = build(5, records, DeBruijnGraph::BASIC);
    auto &dbg_succ = const_cast<DBGSuccinct&>(base_dbg(*graph));
    const std::string middle = "GCTGT";
    node_index node = dbg_succ.kmer_to_node(middle);
    ASSERT_NE(DeBruijnGraph::npos, node);
    auto *mask = dynamic_cast<bit_vector_dyn*>(const_cast<bit_vector*>(dbg_succ.get_mask()));
    ASSERT_TRUE(mask);
    mask->set(node, false);
    auto forward = make_request(Scope::ANY_OFFSET, Strands::FORWARD, true);
    Result r = count_of(*graph, peptide_pattern("MLW", 1), forward);
    EXPECT_EQ(1u, r.anchors->total.value);
    EXPECT_EQ(0u, r.anchors->paths.value);
    check_peptide(*graph, records, DeBruijnGraph::BASIC, "MLW", 1, make_request(
                      Scope::ANY_OFFSET, Strands::BOTH, true), true, { middle });
}


// ---------------------------------------------------------------- the stop '*'

TEST(PatternPeptide, StopResidue) {
    // owner decisions #19 and #21 of 2026-10-08: '*' is a stop codon of the table at that
    // position, X never one; in tables 27, 28 and 31 the codons that stop only in context
    // code their residue (Q, W, E) and '*' matches nothing
    const std::vector<std::string> records {
        "CC" "ATGTAA" "GG",       // M then TAA: a stop in 1, 2, 11; Q in 27 and 28, E in 31
        "AA" "ATGAGA" "TGGCC",    // M then AGA: R in 1 and 11, a stop in 2
        "TT" "TGGTGA" "AAA",      // W then TGA: a stop in 1 and 11, W in 2, 27, 28
    };
    auto graph = build(8, records, DeBruijnGraph::BASIC);
    auto forward = make_request(Scope::ANY_OFFSET, Strands::FORWARD);
    auto count = [&](const char *peptide, int table) {
        const Result r = count_of(*graph, peptide_pattern(peptide, table), forward);
        EXPECT_EQ(Relation::EXACT, r.contexts->total.relation) << peptide << " " << table;
        return r.contexts->total.value;
    };
    // ATGTAA sits in 3 k-mers (offsets 0, 1, 2): "M*" in 1 and 11 (and 2, below); not in 27
    // (TAA is Q)
    for (int table : { 1, 2, 11 }) {
        if (table != 2) {
            EXPECT_EQ(3u, count("M*", table)) << table;
        }
        EXPECT_EQ(0u, count("MQ", table)) << table;
    }
    EXPECT_EQ(0u, count("M*", 27));
    EXPECT_EQ(3u, count("MQ", 27));
    EXPECT_EQ(3u, count("MQ", 28));
    EXPECT_EQ(3u, count("ME", 31));
    // ATGAGA: "M*" in table 2 only (AGA a stop there), with ATGTAA's: 6; "MR" in 1 and 11
    EXPECT_EQ(3u + 3u, count("M*", 2));
    EXPECT_EQ(3u, count("MR", 1));
    EXPECT_EQ(0u, count("MR", 2));
    // TGGTGA: "W*" in 1 and 11 (TGA), "WW" in 2 and 27
    EXPECT_EQ(3u, count("W*", 1));
    EXPECT_EQ(3u, count("W*", 11));
    EXPECT_EQ(0u, count("W*", 2));
    EXPECT_EQ(3u, count("WW", 2));
    EXPECT_EQ(3u, count("WW", 27));
    // X never matches a stop: "MX" misses ATGTAA in 1, finds it in 27 (TAA is Q there);
    // ATGAGA (3 k-mers) and ATGGCC (1, at the end of the second record) in both
    EXPECT_EQ(4u, count("MX", 1));
    EXPECT_EQ(7u, count("MX", 27));
    // the note says why a 0 is a 0 in a table without a stop codon, and only there
    for (int table : { 1, 2, 11, 27, 28, 31 }) {
        const Result r = count_of(*graph, peptide_pattern("M*", table), forward);
        const bool none = !GeneticCode::get(table).stops();
        EXPECT_EQ(none, table == 27 || table == 28 || table == 31);
        EXPECT_EQ(none, std::count(r.notes.begin(), r.notes.end(),
                                   std::string(kNoteNoStopCodon)) > 0) << table;
        const Result plain = count_of(*graph, peptide_pattern("MW", table), forward);
        EXPECT_EQ(0, std::count(plain.notes.begin(), plain.notes.end(),
                                std::string(kNoteNoStopCodon))) << table;
    }

    // against both oracles, every mode, short and long (k 4: the extension through a stop)
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL,
                       DeBruijnGraph::PRIMARY }) {
        for (size_t k : { size_t(4), size_t(5), size_t(8), size_t(11) }) {
            auto g = build(k, records, mode);
            for (int table : { 1, 2, 11, 27, 28, 31 }) {
                for (std::string peptide : { "M*", "*", "W*", "MR", "M*W", "X*", "*X", "R*",
                                             "**", "MX", "*G", "Z*", "M*G" }) {
                    check_peptide(*g, records, mode, peptide, table,
                                  make_request(Scope::ANY_OFFSET, Strands::BOTH,
                                               3 * peptide.size() > k));
                    if (3 * peptide.size() <= k && mode != DeBruijnGraph::PRIMARY) {
                        check_peptide(*g, records, mode, peptide, table,
                                      make_request(Scope::SUFFIX, Strands::FORWARD));
                    }
                }
            }
        }
    }
}

TEST(PatternPeptide, StopIsNeverAMirror) {
    // a peptide is palindromic iff each residue's codons are the reverse complements of its
    // mirror's: in a table with stop codons the stops are no residue's mirror, nor their own
    // (so '*' never makes a peptide palindromic); without stops '*' is the empty set
    for (int id : GeneticCode::ids()) {
        const GeneticCode &code = GeneticCode::get(id);
        const CodonSet stops = code.stops();
        if (!stops) {
            EXPECT_TRUE(id == 27 || id == 28 || id == 31) << id;
            continue;
        }
        const CodonSet mirror = reverse_complement_codons(stops);
        EXPECT_NE(stops, mirror) << id;
        for (char letter : std::string("ACDEFGHIKLMNPQRSTVWYXBZJ")) {
            const Pattern p = Pattern::parse(PatternKind::PROTEIN, std::string(1, letter), code);
            EXPECT_NE(p.codon_sets()[0], mirror) << id << " " << letter;
        }
        EXPECT_EQ(stops, Pattern::parse(PatternKind::PROTEIN, "*", code).codon_sets()[0]);
        EXPECT_FALSE(Pattern::parse(PatternKind::PROTEIN, "*", code).is_palindromic()) << id;
    }
}

TEST(PatternPeptide, NoStopCodonAnswered) {
    // a peptide holding '*' in a table without a stop codon has no instance. L <= k: answered
    // EXACT 0 in every count, nothing searched or charged, the information floor not
    // consulted (nothing to gate). L > k: searched as any other (the anchors are the anchor
    // window's, which may not reach the '*'), no path. Both with the note; in a table with
    // stops the floor applies as ever
    const std::vector<std::string> records { "ATGTAAGGCTGGTGAATGCCC", "ATGCAATGGTAG" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
        auto graph = build(6, records, mode);
        for (int table : { 27, 28, 31 }) {
            for (const char *peptide : { "*", "M*", "**", "*W" }) {
                SCOPED_TRACE(std::string(peptide) + " " + std::to_string(table));
                const Pattern p = peptide_pattern(peptide, table);
                EXPECT_FALSE(p.has_instances());
                Request request = make_request();
                request.min_information_bits = 24;
                const Result r = count_of(*graph, p, request, 0);
                EXPECT_FALSE(r.refusal);
                EXPECT_FALSE(r.stop);
                EXPECT_EQ(0u, r.work.steps);
                EXPECT_TRUE(std::isfinite(r.information_bits));
                EXPECT_EQ(1, std::count(r.notes.begin(), r.notes.end(),
                                        std::string(kNoteNoStopCodon)));
                EXPECT_EQ(Relation::EXACT, r.contexts->total.relation);
                EXPECT_EQ(0u, r.contexts->total.value);
                EXPECT_EQ(7 - p.length(), r.contexts->by_offset.size());
                for (const auto &[offset, c] : r.contexts->by_offset) {
                    EXPECT_EQ(Relation::EXACT, c.relation);
                    EXPECT_EQ(0u, c.value);
                }
                EXPECT_EQ(Relation::EXACT, r.contexts->suffix.relation);
                // enumerate(): nothing, and that is complete
                request.mode = Mode::PARTIAL;
                Budget budget(kManySteps, unbounded_deadline());
                size_t called = 0;
                const Result e = PatternSearch(*graph).enumerate(
                        p, request, budget, [&](const Context &) { ++called; });
                EXPECT_EQ(0u, called);
                EXPECT_TRUE(e.extraction && e.extraction->complete);
                EXPECT_EQ(0u, e.work.steps);
                check_peptide(*graph, records, mode, peptide, table, make_request());
            }
            // longer than k: MQW's anchors exist in table 27/28 (CAA, TAA, TAG are Q), "MQ*"'s
            // anchor window [0, 6) is MQ's: anchors, and no path through the '*'
            for (const char *peptide : { "MQ*", "M*W", "*MQ", "MQW*" }) {
                SCOPED_TRACE(std::string(peptide) + " " + std::to_string(table));
                const Pattern p = peptide_pattern(peptide, table);
                const Result r = count_of(*graph, p, make_request(Scope::ANY_OFFSET,
                                                                  Strands::BOTH, true));
                ASSERT_FALSE(r.refusal);
                EXPECT_EQ(Relation::EXACT, r.anchors->paths.relation);
                EXPECT_EQ(0u, r.anchors->paths.value);
                EXPECT_EQ(1, std::count(r.notes.begin(), r.notes.end(),
                                        std::string(kNoteNoStopCodon)));
                check_peptide(*graph, records, mode, peptide, table,
                              make_request(Scope::ANY_OFFSET, Strands::BOTH, true));
                check_peptide(*graph, records, mode, peptide, table, make_request());
            }
        }
    }
    // in the standard code "*" (three stops, 4.4 bits) is below a floor of 24: refused
    auto graph = build(6, records, DeBruijnGraph::BASIC);
    Request request = make_request();
    request.min_information_bits = 24;
    const Result r = count_of(*graph, peptide_pattern("*", 1), request);
    ASSERT_TRUE(r.refusal);
    EXPECT_EQ("information_below_floor", r.refusal->code);
    EXPECT_EQ(0, std::count(r.notes.begin(), r.notes.end(), std::string(kNoteNoStopCodon)));
}

// a peptide read from a record with its stops: a random frame and strand, translated with the
// oracle's table, some residues replaced by an ambiguity code that admits them (never a stop
// by X), some by '*'
std::string stop_peptide_from(const std::vector<std::string> &records, const OracleCode &code,
                              size_t m, std::mt19937 &rng) {
    for (int attempt = 0; attempt < 30; ++attempt) {
        std::string s = records[rng() % records.size()];
        if (rng() % 2)
            s = rev_comp(s);
        if (s.size() < 3 * m)
            continue;
        const size_t start = rng() % (s.size() - 3 * m + 1);
        std::string peptide;
        for (size_t i = 0; i < m; ++i) {
            peptide.push_back(code.translate(std::string_view(s).substr(start + 3 * i, 3)));
        }
        if (peptide.find('?') != std::string::npos)
            continue;
        for (char &residue : peptide) {
            if (residue == '*' || rng() % 4)
                continue;
            residue = rng() % 2 ? 'X' : '*';
        }
        return peptide;
    }
    std::string peptide;
    for (size_t i = 0; i < m; ++i) {
        peptide.push_back("MKWX*"[rng() % 5]);
    }
    return peptide;
}

TEST(PatternPeptide, RandomStopPeptidesAgainstOracles) {
    // records with planted stop codons, peptides from their six-frame translations with
    // their stops, tables 1, 2, 11 and the tables without a stop (27, 28, 31), every mode,
    // short and long, every strand choice and scope
    const std::vector<std::string> stops { "TAA", "TAG", "TGA", "AGA", "AGG" };
    size_t cases = 0;
    size_t with_stop = 0;
    size_t nonempty = 0;
    for (uint32_t seed = 1; seed <= 40; ++seed) {
        std::mt19937 rng(7100 + seed);
        const size_t k = 4 + rng() % 12;
        const auto mode = static_cast<DeBruijnGraph::Mode>(rng() % 3);
        std::vector<std::string> records(1 + rng() % 3);
        for (std::string &record : records) {
            const size_t length = k + 10 + rng() % 40;
            while (record.size() < length) {
                if (rng() % 5 == 0) {
                    record += stops[rng() % stops.size()];
                } else {
                    record.push_back("ACGT"[rng() % 4]);
                }
            }
        }
        auto graph = build(k, records, mode, rng() % 2);
        const int ids[] = { 1, 2, 11, 27, 28, 31 };
        const int id = ids[rng() % 6];
        const OracleCode &code = oracle_code(id);
        for (int t = 0; t < 6; ++t) {
            const size_t m = 1 + rng() % 6;
            const std::string peptide = stop_peptide_from(records, code, m, rng);
            const bool long_peptide = 3 * m > k;
            Scope scope = Scope::ANY_OFFSET;
            if (!long_peptide && mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0)
                scope = Scope::SUFFIX;
            const auto strands = static_cast<Strands>(rng() % 3);
            SCOPED_TRACE("seed " + std::to_string(seed));
            const size_t found = check_peptide(*graph, records, mode, peptide, id,
                                               make_request(scope, strands, long_peptide));
            ++cases;
            with_stop += peptide.find('*') != std::string::npos;
            nonempty += found > 0;
        }
    }
    EXPECT_EQ(240u, cases);
    EXPECT_LT(80u, with_stop);
    EXPECT_LT(60u, nonempty);
    std::cerr << "random stop peptide cases: " << cases << ", with '*': " << with_stop
              << ", with hits: " << nonempty << std::endl;
}


// ---------------------------------------------------------------- random cases

// a peptide read from a record: a random frame and strand, translated with the oracle's
// table (no stop), some residues replaced by an ambiguity code that admits them
std::string peptide_from(const std::vector<std::string> &records, const OracleCode &code,
                         size_t m, std::mt19937 &rng) {
    for (int attempt = 0; attempt < 30; ++attempt) {
        std::string s = records[rng() % records.size()];
        if (rng() % 2)
            s = rev_comp(s);
        if (s.size() < 3 * m)
            continue;
        const size_t start = rng() % (s.size() - 3 * m + 1);
        std::string peptide;
        for (size_t i = 0; i < m; ++i) {
            peptide.push_back(code.translate(std::string_view(s).substr(start + 3 * i, 3)));
        }
        if (peptide.find_first_of("*?") != std::string::npos)
            continue;
        for (char &residue : peptide) {
            if (rng() % 5)
                continue;
            switch (residue) {
                case 'D': case 'N': residue = rng() % 2 ? 'B' : 'X'; break;
                case 'E': case 'Q': residue = rng() % 2 ? 'Z' : 'X'; break;
                case 'I': case 'L': residue = rng() % 2 ? 'J' : 'X'; break;
                default: residue = 'X'; break;
            }
        }
        return peptide;
    }
    // a random peptide
    const std::string letters = "ACDEFGHIKLMNPQRSTVWYXBZJ";
    std::string peptide;
    for (size_t i = 0; i < m; ++i) {
        peptide.push_back(letters[rng() % letters.size()]);
    }
    return peptide;
}

TEST(PatternPeptide, RandomPeptidesAgainstOracles) {
    // fixed seeds: graphs of every mode and both builders, every table, peptides drawn from
    // the records' six-frame translations, short and long, every strand choice and scope
    size_t cases = 0;
    size_t nonempty = 0;
    size_t long_cases = 0;
    size_t long_nonempty = 0;
    for (uint32_t seed = 1; seed <= 70; ++seed) {
        std::mt19937 rng(5000 + seed);
        const size_t k = 4 + rng() % 13;
        const auto mode = static_cast<DeBruijnGraph::Mode>(rng() % 3);
        const bool batch = rng() % 2;
        std::vector<std::string> records(1 + rng() % 4);
        for (size_t r = 0; r < records.size(); ++r) {
            std::string &record = records[r];
            const size_t length = k + 6 + rng() % 50;
            for (size_t i = 0; i < length; ++i) {
                // an occasional N splits a record into islands (DNA4)
                record.push_back(rng() % 40 ? "ACGT"[rng() % 4] : 'N');
            }
            // a copy of a piece of an earlier record with one change: shared k-mers, joins
            if (r && rng() % 2) {
                const std::string &earlier = records[rng() % r];
                std::string piece = earlier.substr(0, std::min<size_t>(earlier.size(), 20));
                piece[rng() % piece.size()] = "ACGT"[rng() % 4];
                record += piece;
            }
        }
        auto graph = build(k, records, mode, batch);
        const int id = rng() % 2 ? 1 : GeneticCode::ids()[rng() % GeneticCode::ids().size()];
        const OracleCode &code = oracle_code(id);

        for (int t = 0; t < 5; ++t) {
            const size_t m = 1 + rng() % 6;
            const std::string peptide = peptide_from(records, code, m, rng);
            const auto strands = static_cast<Strands>(rng() % 3);
            const bool long_peptide = 3 * m > k;
            Scope scope = Scope::ANY_OFFSET;
            if (!long_peptide && mode != DeBruijnGraph::PRIMARY && rng() % 3 == 0)
                scope = Scope::SUFFIX;
            const bool extend = long_peptide && rng() % 4;
            SCOPED_TRACE("seed " + std::to_string(seed) + " batch " + std::to_string(batch));
            const size_t found = check_peptide(*graph, records, mode, peptide, id,
                                               make_request(scope, strands, extend));
            ++cases;
            nonempty += found > 0;
            long_cases += extend;
            long_nonempty += extend && found > 0;
        }
    }
    EXPECT_EQ(350u, cases);
    // not vacuous: most peptides come from the records
    EXPECT_LT(150u, nonempty);
    EXPECT_LT(40u, long_cases);
    EXPECT_LT(20u, long_nonempty);
    std::cerr << "random peptide cases: " << cases << ", with hits: " << nonempty
              << ", long with extension: " << long_cases << ", with paths: " << long_nonempty
              << std::endl;
}


// ---------------------------------------------------------------- budgets and determinism

// |count| states a true relation to |truth| (§3)
bool holds(const Count &count, uint64_t truth) {
    switch (count.relation) {
        case Relation::EXACT: return count.value == truth;
        case Relation::AT_LEAST: return count.value <= truth;
        case Relation::BOUNDS:
            return count.lower == count.value && count.lower <= truth && truth <= count.upper;
        case Relation::UNKNOWN: return count.value == 0;
    }
    return false;
}

TEST(PatternPeptide, StepStopsStateTrueRelations) {
    // every step cap from 0 to the full run's: the counts stay true, the stop is named, and
    // a stop in the extension leaves the anchors exact and the paths at least
    const std::vector<std::string> records {
        "ATGCTGAGCAGATGGAAACATTAGGGATGTTAAGTCGG", "CCATGTTGTCCAGGTGGAAG"
    };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
        for (size_t k : { 5, 7, 13 }) {
            auto graph = build(k, records, mode);
            for (std::string peptide : { "MLSR", "MJXRW" }) {
                const bool long_peptide = 3 * peptide.size() > k;
                const Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH,
                                                     long_peptide);
                const Pattern pattern = peptide_pattern(peptide, 1);
                const Result full = count_of(*graph, pattern, request);
                ASSERT_FALSE(full.stop);
                const uint64_t steps = full.work.steps;
                const Count &full_total = long_peptide ? full.anchors->total
                                                       : full.contexts->total;
                for (uint64_t cap = 0; cap <= steps; cap += 1 + steps / 150) {
                    Result r = count_of(*graph, pattern, request, cap);
                    const Count &total = long_peptide ? r.anchors->total : r.contexts->total;
                    EXPECT_TRUE(holds(total, full_total.value)) << peptide << " cap " << cap;
                    if (cap < steps) {
                        ASSERT_TRUE(r.stop);
                        EXPECT_EQ(StopReason::MAX_STEPS, r.stop->reason);
                    }
                    if (long_peptide) {
                        EXPECT_TRUE(holds(r.anchors->paths, full.anchors->paths.value))
                            << peptide << " cap " << cap;
                        if (r.stop && r.stop->phase == StopPhase::EXTENSION) {
                            EXPECT_EQ(Relation::EXACT, r.anchors->total.relation);
                            EXPECT_EQ(Extension::STOPPED, r.anchors->extension);
                            EXPECT_EQ(Relation::AT_LEAST, r.anchors->paths.relation);
                        }
                    }
                }
            }
        }
    }
}

TEST(PatternPeptide, Deterministic) {
    const std::vector<std::string> records { "ATGCTGAGCAGATGGAAACAT", "GGATGTTAAGTCGGTGGAAG" };
    for (auto mode : { DeBruijnGraph::BASIC, DeBruijnGraph::PRIMARY }) {
        auto graph = build(7, records, mode);
        const StoredNodes stored(*graph);
        PatternSearch engine(*graph);
        for (bool extend : { false, true }) {
            Request request = make_request(Scope::ANY_OFFSET, Strands::BOTH, extend);
            Result a;
            Result b;
            if (extend) {
                auto first = run_paths(engine, stored, peptide_pattern("MLSRW", 1), request, &a);
                auto second = run_paths(engine, stored, peptide_pattern("MLSRW", 1), request, &b);
                EXPECT_EQ(first, second);
                EXPECT_FALSE(first.empty());
            } else {
                auto first = run_contexts(engine, stored, peptide_pattern("ML", 1), request, &a);
                auto second = run_contexts(engine, stored, peptide_pattern("ML", 1), request, &b);
                EXPECT_EQ(first, second);
                EXPECT_FALSE(first.empty());
            }
            EXPECT_EQ(a.work.steps, b.work.steps);
            EXPECT_EQ(a.work.ranges_visited, b.work.ranges_visited);
        }
    }
}

#endif // _DNA_GRAPH || _DNA5_GRAPH

} // namespace
