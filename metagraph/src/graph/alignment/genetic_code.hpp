#ifndef __GENETIC_CODE_HPP__
#define __GENETIC_CODE_HPP__

/**
 * The genetic codes of the pattern search's peptides (docs/DESIGN-pattern-search.md §6): every
 * NCBI translation table, and the codon sets the codon automaton of a peptide pattern is built
 * from.
 *
 * A codon is indexed in reading order over the engine's base order (pattern_search.hpp's
 * BaseSet bits: A = 0, C = 1, G = 2, T = 3): index = 16 * b0 + 4 * b1 + b2. A set of codons is
 * a 64-bit mask (CodonSet, bit c = codon index c). NCBI lists a table as 64 one-letter
 * residues in TCAG order (its "ncbieaa" string); amino_acids() returns that string as NCBI
 * writes it, and codons() and translate() convert between the two orders.
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>


namespace mtg {
namespace graph {
namespace pattern {

// a set of codons: bit c is the codon of index c (16 * b0 + 4 * b1 + b2, A C G T = 0 1 2 3)
using CodonSet = uint64_t;

constexpr CodonSet kAllCodons = ~CodonSet(0);

// the index of a base in the engine's order (A 0, C 1, G 2, T 3; case-insensitive); -1 for
// any other character
int base_index(char base);

// the index of a codon of three bases (case-insensitive); -1 when a character is not A, C, G, T
int codon_index(std::string_view codon);

// { rc(c) : c in |set| }: each codon reversed and complemented
CodonSet reverse_complement_codons(CodonSet set);

/**
 * The bases that can follow, at position |position_in_codon| (0, 1 or 2) of a codon of |set|,
 * the bases already spelled in that codon: |b0| at position 0 and |b1| at position 1, each a
 * base index, or -1 when not known (a position before a pattern window's first, §4.1). Exactly
 * the bases b for which some codon of |set| agrees with every known base and has b at
 * |position_in_codon|: no superset (Leu after C: T only; after CT: A C G T; after T: T; after
 * TT: A G). The result is a BaseSet (bit 0 A, bit 1 C, bit 2 G, bit 3 T); 0 when no codon
 * agrees.
 */
uint8_t next_bases(CodonSet set, size_t position_in_codon, int b0 = -1, int b1 = -1);

/**
 * The cylinder of |set| over the codon positions [begin, end) (0 <= begin <= end <= 3): every
 * codon that agrees with some codon of |set| at those positions, the other positions free.
 * Two codon sets admit the same strings over [begin, end) iff their cylinders are equal.
 */
CodonSet cylinder(CodonSet set, size_t begin, size_t end);

// the number of distinct strings the codons of |set| spell at the positions [begin, end)
uint32_t distinct_projections(CodonSet set, size_t begin, size_t end);


/**
 * One NCBI translation table (the NCBI Taxonomy's genetic codes, data of gc.prt version 4.6,
 * https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt): ids 1-6, 9-16 and 21-33, 1 the standard
 * code (the default). The ids NCBI retired (7 and 8, merged into 4 and 1) and the unassigned
 * ones (17-20) are not tables.
 *
 * A codon codes the residue of NCBI's ncbieaa string; '*' there is a stop. Tables 27, 28 and
 * 31 list some codons both as a residue (ncbieaa) and as a terminator in context (sncbieaa
 * '*'): Karyorelict (27) TGA W; Condylostoma (28) TAA Q, TAG Q, TGA W; Blastocrithidia (31)
 * TAA E, TAG E. They code their residue here (codons()) and are named by context_stops(): a
 * pattern can match them as that residue, which a translation that ends at them in context
 * would not show.
 */
class GeneticCode {
  public:
    static constexpr int kStandard = 1;

    // the table with NCBI id |id|; nullptr when NCBI has no table with that id
    static const GeneticCode* find(int id);
    // the table with NCBI id |id|; std::invalid_argument naming the ids otherwise
    static const GeneticCode& get(int id);
    static const GeneticCode& standard() { return get(kStandard); }
    // every id, ascending: 1-6, 9-16, 21-33
    static const std::vector<int>& ids();
    // the ids as a text for messages: "1-6, 9-16, 21-33"
    static std::string ids_text();

    int id() const { return id_; }
    // NCBI's name of the table
    const char* name() const { return name_; }
    // NCBI's ncbieaa string: 64 residues in TCAG order (TTT TTC TTA TTG TCT ... GGG), '*' a stop
    std::string_view amino_acids() const { return amino_acids_; }

    // the residue a codon (three bases, case-insensitive) codes: an upper-case letter, '*'
    // for a stop, 0 when |codon| is not three of A, C, G, T
    char translate(std::string_view codon) const;
    char translate(int codon) const { return codon >= 0 && codon < 64 ? by_index_[codon] : 0; }

    // the codons of |residue|, one of the 20 amino acids (case-insensitive); 0 for any other
    // character, the stop '*' included
    CodonSet codons(char residue) const;
    // the stop codons ('*' in ncbieaa)
    CodonSet stops() const { return stops_; }
    // the codons that code a residue and are also a terminator in context (see the class)
    CodonSet context_stops() const { return context_stops_; }

  private:
    GeneticCode(int id, const char *name, const char *amino_acids, const char *context_stops);

    int id_;
    const char *name_;
    const char *amino_acids_;
    // the residue per codon index (engine order)
    char by_index_[64];
    CodonSet stops_ = 0;
    CodonSet context_stops_ = 0;
    // per letter 'A' .. 'Z', the codons coding it
    CodonSet by_residue_[26] = {};

    static const std::vector<GeneticCode>& tables();
};

} // namespace pattern
} // namespace graph
} // namespace mtg

#endif // __GENETIC_CODE_HPP__
