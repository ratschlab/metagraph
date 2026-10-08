#include "genetic_code.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>


namespace mtg {
namespace graph {
namespace pattern {

int base_index(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

int codon_index(std::string_view codon) {
    if (codon.size() != 3)
        return -1;
    int index = 0;
    for (char c : codon) {
        int b = base_index(c);
        if (b < 0)
            return -1;
        index = 4 * index + b;
    }
    return index;
}

std::string codon_string(int index) {
    static constexpr char kBases[] = "ACGT";
    return { kBases[(index >> 4) & 3], kBases[(index >> 2) & 3], kBases[index & 3] };
}

CodonSet reverse_complement_codons(CodonSet set) {
    CodonSet result = 0;
    for (int c = 0; c < 64; ++c) {
        if (!(set >> c & 1))
            continue;
        // complement: 3 - b (A <-> T, C <-> G); reverse: b2 b1 b0
        const int b0 = 3 - (c >> 4), b1 = 3 - ((c >> 2) & 3), b2 = 3 - (c & 3);
        result |= CodonSet(1) << (16 * b2 + 4 * b1 + b0);
    }
    return result;
}

uint8_t next_bases(CodonSet set, size_t position_in_codon, int b0, int b1) {
    // the codons agreeing with the known bases, folded onto the position asked
    uint8_t result = 0;
    for (int c = 0; c < 64; ++c) {
        if (!(set >> c & 1))
            continue;
        const int bases[3] = { c >> 4, (c >> 2) & 3, c & 3 };
        if (position_in_codon > 0 && b0 >= 0 && bases[0] != b0)
            continue;
        if (position_in_codon > 1 && b1 >= 0 && bases[1] != b1)
            continue;
        result |= uint8_t(1) << bases[position_in_codon];
    }
    return result;
}

namespace {

// the bases of codon |c| at the positions [begin, end), packed two bits each
int projection(int c, size_t begin, size_t end) {
    int key = 0;
    for (size_t i = begin; i < end; ++i) {
        key = 4 * key + ((c >> (2 * (2 - i))) & 3);
    }
    return key;
}

} // namespace

CodonSet cylinder(CodonSet set, size_t begin, size_t end) {
    // the projections present (at most 4^3 keys)
    CodonSet seen = 0;
    for (int c = 0; c < 64; ++c) {
        if (set >> c & 1)
            seen |= CodonSet(1) << projection(c, begin, end);
    }
    CodonSet result = 0;
    for (int c = 0; c < 64; ++c) {
        if (seen >> projection(c, begin, end) & 1)
            result |= CodonSet(1) << c;
    }
    return result;
}

uint32_t distinct_projections(CodonSet set, size_t begin, size_t end) {
    CodonSet seen = 0;
    for (int c = 0; c < 64; ++c) {
        if (set >> c & 1)
            seen |= CodonSet(1) << projection(c, begin, end);
    }
    return __builtin_popcountll(seen);
}


GeneticCode::GeneticCode(int id, const char *name, const char *amino_acids,
                         const char *context_stops)
      : id_(id), name_(name), amino_acids_(amino_acids) {
    // NCBI's order: TCAG at each codon position
    static constexpr char kTCAG[] = "TCAG";
    for (int i = 0; i < 64; ++i) {
        const char codon[3] = { kTCAG[i >> 4], kTCAG[(i >> 2) & 3], kTCAG[i & 3] };
        const int c = codon_index(std::string_view(codon, 3));
        const char residue = amino_acids[i];
        by_index_[c] = residue;
        if (residue == '*') {
            stops_ |= CodonSet(1) << c;
        } else {
            by_residue_[residue - 'A'] |= CodonSet(1) << c;
        }
    }
    // "TAA TAG TGA": codons separated by one space
    for (std::string_view s = context_stops; s.size() >= 3;
            s.remove_prefix(std::min<size_t>(4, s.size()))) {
        context_stops_ |= CodonSet(1) << codon_index(s.substr(0, 3));
    }
}

const std::vector<GeneticCode>& GeneticCode::tables() {
    // NCBI gc.prt version 4.6 (id, name, ncbieaa; the codons sncbieaa also marks as a
    // terminator while ncbieaa codes a residue)
    static const std::vector<GeneticCode> kTables {
        { 1, "Standard",
          "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 2, "Vertebrate Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "" },
        { 3, "Yeast Mitochondrial",
          "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 4, "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; "
             "Mycoplasma; Spiroplasma",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 5, "Invertebrate Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "" },
        { 6, "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
          "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 9, "Echinoderm Mitochondrial; Flatworm Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "" },
        { 10, "Euplotid Nuclear",
          "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 11, "Bacterial, Archaeal and Plant Plastid",
          "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 12, "Alternative Yeast Nuclear",
          "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 13, "Ascidian Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "" },
        { 14, "Alternative Flatworm Mitochondrial",
          "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "" },
        { 15, "Blepharisma Macronuclear",
          "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 16, "Chlorophycean Mitochondrial",
          "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 21, "Trematode Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "" },
        { 22, "Scenedesmus obliquus Mitochondrial",
          "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 23, "Thraustochytrium Mitochondrial",
          "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 24, "Rhabdopleuridae Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "" },
        { 25, "Candidate Division SR1 and Gracilibacteria",
          "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 26, "Pachysolen tannophilus Nuclear",
          "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 27, "Karyorelict Nuclear",
          "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "TGA" },
        { 28, "Condylostoma Nuclear",
          "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
          "TAA TAG TGA" },
        { 29, "Mesodinium Nuclear",
          "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 30, "Peritrich Nuclear",
          "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 31, "Blastocrithidia Nuclear",
          "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "TAA TAG" },
        { 32, "Balanophoraceae Plastid",
          "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "" },
        { 33, "Cephalodiscidae Mitochondrial",
          "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "" },
    };
    return kTables;
}

const GeneticCode* GeneticCode::find(int id) {
    for (const GeneticCode &table : tables()) {
        if (table.id() == id)
            return &table;
    }
    return nullptr;
}

const GeneticCode& GeneticCode::get(int id) {
    if (const GeneticCode *table = find(id))
        return *table;
    throw std::invalid_argument("genetic code " + std::to_string(id)
                                + " is not an NCBI translation table (" + ids_text() + ")");
}

const std::vector<int>& GeneticCode::ids() {
    static const std::vector<int> kIds = []() {
        std::vector<int> ids;
        for (const GeneticCode &table : tables()) {
            ids.push_back(table.id());
        }
        return ids;
    }();
    return kIds;
}

std::string GeneticCode::ids_text() {
    // runs of consecutive ids
    const std::vector<int> &all = ids();
    std::string text;
    for (size_t i = 0; i < all.size();) {
        size_t j = i;
        while (j + 1 < all.size() && all[j + 1] == all[j] + 1) {
            ++j;
        }
        if (!text.empty())
            text += ", ";
        text += std::to_string(all[i]);
        if (j > i)
            text += "-" + std::to_string(all[j]);
        i = j + 1;
    }
    return text;
}

char GeneticCode::translate(std::string_view codon) const {
    return translate(codon_index(codon));
}

CodonSet GeneticCode::codons(char residue) const {
    const int c = std::toupper(static_cast<unsigned char>(residue));
    return c >= 'A' && c <= 'Z' ? by_residue_[c - 'A'] : 0;
}

} // namespace pattern
} // namespace graph
} // namespace mtg
