#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../graph/test_aligner_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "kmer/alphabets.hpp"
#include "seq_io/sequence_io.hpp"


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

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
    LabeledDBGAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::string("GAATGCAT") },
                                       { std::string("B"), std::string("CGAATGCCT") },
                                       { std::string("A"), std::string("TGCCT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(targets.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end()) << label;
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment << " " << alignment.target_column;
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
        { std::string("TGAAATGCAT"), {{ { std::string("C"), std::string("TGGAATGCAT") },
                                        // { std::string("A"), std::string("TGCCT") },
                                        { std::string("B"), std::string("AATGCCT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(targets.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment << " " << alignment.target_column;
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

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
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
            auto alignments = aligner.align(query);
            EXPECT_EQ(targets.size(), alignments.size()) << query;

            for (const auto &alignment : alignments) {
                bool found = false;
                for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                    auto find = targets.find(label);
                    ASSERT_TRUE(find != targets.end());
                    if (alignment.get_sequence() == find->second) {
                        found = true;
                        break;
                    }
                }
                EXPECT_TRUE(found) << alignment << " " << alignment.target_column;
            }
        }
    }
}
#endif

} // namespace
