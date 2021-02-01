#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../graph/test_aligner_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "graph/alignment/aligner_labeled.hpp"
#include "seq_io/sequence_io.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;


inline std::vector<std::string>
get_alignment_labels(const AnnotatedDBG &anno_graph,
                     const LabeledDBGAligner<>::DBGAlignment &alignment) {
    return anno_graph.get_labels(alignment.get_sequence(), 1.0);
}

template <typename GraphAnnotationPair>
class LabeledDBGAlignerTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledDBGAlignerTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledDBGAlignerTest, SimpleTangleGraph) {
    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    LabeledDBGAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::string("GAATGCAT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        for (const auto &alignment : aligner.align(query)) {
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());
                EXPECT_EQ(find->second, alignment.get_sequence());
            }
        }
    }
}

TEST(LabeledDBGAlignerTest, SimpleTangleGraphSuffixSeed) {
    size_t k = 4;
    /*  B    B                  AB   AB
       TCGA-CGAA                TGCC-GCCT
                \ BC   BC   BC /
                 GAAT-AATG-ATGC
          C    C/              \   C    C
       TGGA-GGAA                TGCA-GCAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "TCGAATGCCT",
        "TGGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto graph = std::dynamic_pointer_cast<DBGSuccinct>(build_graph<DBGSuccinct>(k, sequences));
    graph->reset_mask();
    auto anno_graph = build_anno_graph<annot::ColumnCompressed<>>(graph, sequences, labels);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
    config.min_seed_length = 2;
    LabeledDBGAligner<SuffixSeeder<ExactSeeder<>>> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("TGAAATGCAT"), {{ { std::string("C"), std::string("TGGAATGCAT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        for (const auto &alignment : aligner.align(query)) {
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());
                EXPECT_EQ(find->second, alignment.get_sequence());
            }
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(LabeledDBGAlignerTest, CanonicalTangleGraph) {
    size_t k = 5;
    /*   B     B                AB    AB
       TTAGT-TAGTC             TCGAA-CGAAA
                  \  BC   ABC /
                   AGTCG-GTCGA
          C     C /           \   C     C
       TCAGT-CAGTC             TCGAT-CGATT
        AB    AB                 B     B
       TTTCG-TTCGA             GACTA-ACTAA
                  \ ABC    BC /
                   TCGAC-CGACT
          C     C /           \   C     C
       AATCG-ATCGA             GACTG-ACTGA
    */
    const std::vector<std::string> sequences {
        "GTCGAAA", // "TTTCGAC"
        "TTAGTCGAAA", // "TTTCGACTAA"
        "TCAGTCGATT" // "AATCGACTGA"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    for (DBGMode mode : { DBGMode::CANONICAL, DBGMode::CANONICAL_WRAPPER }) {
        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k, sequences, labels, mode
        );

        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        LabeledDBGAligner<> aligner(*anno_graph, config);

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
            { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") } }} }
        }};

        for (const auto &[query, targets] : exp_alignments) {
            const auto &alignments = aligner.align(query);

            EXPECT_EQ(targets.size(), alignments.size()) << query;

            for (const auto &alignment : alignments) {
                for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                    auto find = targets.find(label);
                    ASSERT_TRUE(find != targets.end());

                    EXPECT_EQ(find->second, alignment.get_sequence())
                        << query << " " << label;
                }
            }
        }
    }
}
#endif

} // namespace
