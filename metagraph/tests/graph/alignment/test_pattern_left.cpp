/**
 * The Model's arms inside the pattern (Model::start_at, start_left, bases_left, prev) as Pattern
 * implements them, against the oracle's reading of the pattern (pattern_segment_oracles.hpp: the
 * IUPAC codes and every NCBI genetic code decoded by the test, a string placed at an offset
 * admitted iff every residue it overlaps has a codon agreeing with it there):
 *  - for every table, every residue letter and both orientations, every codon position and every
 *    known later part of the codon: bases_left equals the bases the oracle admits there;
 *  - random peptides at every phase: from an anchor sampled from an instance, the left arm
 *    (start_left, then bases_left and prev step by step) and the right arm (start_at, then bases
 *    and next) admit exactly the oracle's bases at every step, and the right arm agrees with
 *    Pattern::allowed on the spelled bases;
 *  - DNA and IUPAC carry no state; a model without a left arm throws.
 */
#include <gtest/gtest.h>

#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "pattern_segment_oracles.hpp"

#include "graph/alignment/genetic_code.hpp"
#include "graph/alignment/pattern_search.hpp"


namespace {

using namespace mtg::graph::pattern;
using namespace mtg::test::segments;

constexpr char kResidues[] = "ACDEFGHIKLMNPQRSTVWYXBZJ*";

// the bases (BaseSet bits A C G T) b with b + x admitted at offset t - 1
BaseSet oracle_left(const OracleQ &q, const std::string &x, size_t t) {
    BaseSet set = 0;
    for (size_t b = 0; b < 4; ++b) {
        if (q.admits(std::string(1, kOracleBases[b]) + x, t - 1))
            set |= BaseSet(1) << b;
    }
    return set;
}

// the bases b with x + b admitted at offset t
BaseSet oracle_right(const OracleQ &q, const std::string &x, size_t t) {
    BaseSet set = 0;
    for (size_t b = 0; b < 4; ++b) {
        if (q.admits(x + kOracleBases[b], t))
            set |= BaseSet(1) << b;
    }
    return set;
}

char pick(BaseSet set, std::mt19937 &rng) {
    std::vector<char> bases;
    for (size_t b = 0; b < 4; ++b) {
        if (set & (1 << b))
            bases.push_back(kOracleBases[b]);
    }
    return bases[rng() % bases.size()];
}

TEST(PatternLeft, EveryResidueOfEveryTableBothOrientations) {
    for (int id : GeneticCode::ids()) {
        const GeneticCode &code = GeneticCode::get(id);
        for (const char residue : std::string(kResidues)) {
            const std::string text(1, residue);
            const Pattern p = Pattern::parse(PatternKind::PROTEIN, text, code);
            const Pattern rc = p.reverse_complement();
            const std::vector<OracleQ> q = oracle_orientations(PatternKind::PROTEIN, text,
                                                               Strands::BOTH, id);
            // a peptide of one residue is never palindromic in a table with stops; in one
            // without, X is: its oracle has one orientation, read both ways
            const OracleQ &oq_forward = q.front();
            const OracleQ &oq_reverse = q.back();
            for (int orientation = 0; orientation < 2; ++orientation) {
                const Pattern &model = orientation ? rc : p;
                const OracleQ &oracle = orientation ? oq_reverse : oq_forward;
                // every known later part y of the codon (positions 3 - |y| .. 2)
                for (size_t known = 0; known < 3; ++known) {
                    const size_t j = 2 - known;
                    for (size_t index = 0; index < (size_t(1) << (2 * known)); ++index) {
                        std::string y;
                        for (size_t i = 0; i < known; ++i) {
                            y.push_back(kOracleBases[(index >> (2 * i)) & 3]);
                        }
                        // the state reached by prepending y's bases from position 2 down,
                        // only along bases the arm admits (prev's precondition)
                        Model::State state = model.start_left(3, "");
                        bool reachable = true;
                        for (size_t i = known; i-- > 0 && reachable;) {
                            const uint32_t position = static_cast<uint32_t>(j + 1 + i);
                            const BaseSet allowed = model.bases_left(state, position);
                            const char base = y[i];
                            if (!(allowed & (1 << base_index(base)))) {
                                reachable = false;
                                break;
                            }
                            state = model.prev(state, position, base);
                        }
                        if (!reachable) {
                            // and then the oracle admits y at its place neither
                            EXPECT_FALSE(oracle.admits(y, j + 1))
                                << id << " " << residue << " " << orientation << " " << y;
                            continue;
                        }
                        EXPECT_EQ(oracle_left(oracle, y, j + 1),
                                  model.bases_left(state, static_cast<uint32_t>(j)))
                            << "table " << id << " residue " << residue << " orientation "
                            << orientation << " known " << y;
                    }
                }
            }
        }
    }
}

// a random peptide over every residue letter, and an instance of its orientation |q| (an
// empty string when a residue admits no codon)
std::string random_instance(const OracleQ &q, std::mt19937 &rng) {
    std::string instance;
    for (const auto &codons : q.codons) {
        if (codons.empty())
            return "";
        instance += codons[rng() % codons.size()];
    }
    return instance;
}

TEST(PatternLeft, ArmsAgainstTheOracleAtEveryPhase) {
    std::mt19937 rng(7);
    const std::vector<int> &ids = GeneticCode::ids();
    uint64_t steps = 0;
    for (size_t round = 0; round < 3000; ++round) {
        const int id = ids[rng() % ids.size()];
        const size_t m = 3 + rng() % 8;
        std::string text;
        for (size_t i = 0; i < m; ++i) {
            // mostly residues, sometimes the ambiguity codes, rarely the stop
            const size_t roll = rng() % 30;
            text.push_back(roll < 24 ? kResidues[roll % 20]
                                     : roll < 29 ? kResidues[20 + roll % 4] : '*');
        }
        const Pattern p = Pattern::parse(PatternKind::PROTEIN, text, GeneticCode::get(id));
        const Pattern rc = p.reverse_complement();
        const std::vector<OracleQ> qs = oracle_orientations(PatternKind::PROTEIN, text,
                                                            Strands::BOTH, id);
        for (const OracleQ &q : qs) {
            const Pattern &model = q.orientation == Orientation::REVERSE ? rc : p;
            const std::string instance = random_instance(q, rng);
            if (instance.empty())
                continue;
            const size_t L = instance.size();
            const size_t length = 2 + rng() % std::min<size_t>(11, L - 1);
            const size_t s = rng() % (L - length + 1);
            const std::string anchor = instance.substr(s, length);
            ASSERT_TRUE(q.admits(anchor, s));

            // the left arm from the anchor down to position 0
            std::string spelled = anchor;
            Model::State state = model.start_left(static_cast<uint32_t>(s), anchor);
            for (size_t t = s; t > 0; --t) {
                const BaseSet expected = oracle_left(q, spelled, t);
                const BaseSet got = model.bases_left(state, static_cast<uint32_t>(t - 1));
                ASSERT_EQ(expected, got) << text << " table " << id << " offset " << s
                                         << " at " << t - 1 << " spelled " << spelled;
                ASSERT_NE(0, got) << "an admitted string extends to the left";
                const char base = pick(got, rng);
                state = model.prev(state, static_cast<uint32_t>(t - 1), base);
                spelled.insert(spelled.begin(), base);
                ++steps;
            }

            // the right arm from the anchor up to L
            spelled = anchor;
            state = model.start_at(static_cast<uint32_t>(s), anchor);
            for (size_t position = s + length; position < L; ++position) {
                const BaseSet expected = oracle_right(q, spelled, s);
                const BaseSet got = model.bases(state, static_cast<uint32_t>(position));
                ASSERT_EQ(expected, got) << text << " table " << id << " offset " << s
                                         << " at " << position << " spelled " << spelled;
                // allowed() reads the codon prefix from the spelled bases themselves
                ASSERT_EQ(model.allowed(position, spelled), got);
                const char base = pick(got, rng);
                state = model.next(state, static_cast<uint32_t>(position), base);
                spelled.push_back(base);
                ++steps;
            }
        }
    }
    EXPECT_GT(steps, 10000u);
}

TEST(PatternLeft, StartAtOffsetZeroIsStart) {
    std::mt19937 rng(11);
    const Pattern p = Pattern::parse(PatternKind::PROTEIN, "MKVLAWQRS");
    for (size_t round = 0; round < 500; ++round) {
        std::string kmer;
        const size_t length = 1 + rng() % 12;
        for (size_t i = 0; i < length; ++i) {
            kmer.push_back("ACGTN"[rng() % 5]);
        }
        EXPECT_EQ(p.start(kmer), p.start_at(0, kmer)) << kmer;
    }
}

TEST(PatternLeft, AnchorWithoutACGTAdmitsNothing) {
    // a base other than A, C, G, T among the known bases of a codon admits nothing, both ways
    const Pattern p = Pattern::parse(PatternKind::PROTEIN, "MKVLAW");
    // position 4's codon (residue 1) holds positions 3, 4, 5: an anchor at 4 knows 4 and 5
    EXPECT_EQ(0, p.bases_left(p.start_left(4, "ANGTTC"), 3));
    EXPECT_NE(0, p.bases_left(p.start_left(4, "AAGTTC"), 3));
    // start_at(1, x) of |x| = 4 reads the codon prefix at position 5: x's last two bases
    EXPECT_EQ(0, p.bases(p.start_at(1, "TGAN"), 5));
    EXPECT_EQ(p.allowed(5, "ATGAA"), p.bases(p.start_at(1, "TGAA"), 5));
}

TEST(PatternLeft, DnaAndIupacCarryNoState) {
    std::mt19937 rng(3);
    for (size_t round = 0; round < 200; ++round) {
        std::string text;
        const size_t L = 5 + rng() % 20;
        for (size_t i = 0; i < L; ++i) {
            text.push_back("ACGTRYSWKMBDHVN"[rng() % 15]);
        }
        for (const Pattern &p : { Pattern::parse(PatternKind::IUPAC, text),
                                  Pattern::parse(PatternKind::IUPAC, text).reverse_complement() }) {
            const uint32_t s = static_cast<uint32_t>(rng() % L);
            EXPECT_EQ(0u, p.start_at(s, "ACG"));
            EXPECT_EQ(0u, p.start_left(s, "ACG"));
            for (uint32_t position = 0; position < L; ++position) {
                EXPECT_EQ(p.positions()[position], p.bases_left(0, position));
                for (size_t b = 0; b < 4; ++b) {
                    if (p.positions()[position] & (1 << b)) {
                        EXPECT_EQ(0u, p.prev(0, position, kOracleBases[b]));
                    }
                }
            }
        }
    }
}

// a model of the right arm only, as a mismatch model or a profile would be before it steps left
class RightOnly : public Model {
  public:
    size_t length() const override { return 12; }
    State start(std::string_view) const override { return 5; }
    BaseSet bases(State, uint32_t) const override { return kAllBases; }
    State next(State s, uint32_t, char) const override { return s; }
};

TEST(PatternLeft, ModelWithoutLeftArmThrows) {
    const RightOnly model;
    EXPECT_EQ(5u, model.start_at(0, "ACGT"));
    EXPECT_THROW(model.start_at(1, "ACGT"), std::logic_error);
    EXPECT_THROW(model.start_left(1, "ACGT"), std::logic_error);
    EXPECT_THROW(model.bases_left(0, 1), std::logic_error);
    EXPECT_THROW(model.prev(0, 1, 'A'), std::logic_error);
}

} // namespace
