#include "gtest/gtest.h"

#include <unordered_map>

#include "test_annotated_dbg_helpers.hpp"

#include "graph/alignment/dbg_masked_aligner.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;


template <typename GraphAnnotationPair>
class MaskedDBGAlignerTest : public ::testing::Test {};
TYPED_TEST_SUITE(MaskedDBGAlignerTest, GraphWithNAnnotationPairTypes);


void check_alignment_labels(const AnnotatedDBG &anno_graph,
                            const QueryAlignment<DeBruijnGraph::node_index> &alignments) {
    for (const auto &alignment : alignments) {
        const auto ref_labels = anno_graph.get_labels(alignment.get_sequence(), 1.0);
        const auto &align_labels = alignment.get_labels();

        EXPECT_EQ(1u, align_labels.size());
        EXPECT_NE(ref_labels.end(),
                  std::find(ref_labels.begin(), ref_labels.end(), align_labels[0]));
    }
}


TYPED_TEST(MaskedDBGAlignerTest, SimpleTangleGraph) {
    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAC
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAC",
        "TTTTTTTTTTTTTT"
    };
    const std::vector<std::string> labels { "A", "B", "C", "D" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    MaskedDBGAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCCT"), {{ { std::string("A"), std::string("TGCCT") },
                                       { std::string("B"), std::string("CGAATGCCT") },
                                       { std::string("C"), std::string("GAATGC") } }} },
        { std::string("TTTTATTTTTATT"), {{ { std::string("D"), std::string("TTTTTTTTTTTTT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        const auto &alignments = aligner.align(query);

        EXPECT_EQ(targets.size(), alignments.size());
        check_alignment_labels(*anno_graph, alignments);

        for (const auto &alignment : alignments) {
            for (const auto &label : alignment.get_labels()) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());

                EXPECT_EQ(find->second, alignment.get_sequence())
                    << query << " " << label;
            }
        }
    }
}

TYPED_TEST(MaskedDBGAlignerTest, SimpleTangleGraphHighCoverage) {
    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAC
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAC",
        "TTTTTTTTTTTTTT"
    };
    const std::vector<std::string> labels { "A", "B", "C", "D" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.exact_kmer_match_fraction = 5.0 / 11.0;
    MaskedDBGAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCCT"), {{ { std::string("B"), std::string("CGAATGCCT") },
                                       { std::string("C"), std::string("GAATGC") } }} },
        { std::string("TTTTATTTTTATT"), {{ { std::string("D"), std::string("TTTTTTTTTTTTT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        const auto &alignments = aligner.align(query);

        EXPECT_EQ(targets.size(), alignments.size()) << query;
        check_alignment_labels(*anno_graph, alignments);

        for (const auto &alignment : alignments) {
            for (const auto &label : alignment.get_labels()) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());

                EXPECT_EQ(find->second, alignment.get_sequence())
                    << query << " " << label;
            }
        }
    }
}

TYPED_TEST(MaskedDBGAlignerTest, CanonicalTangleGraph) {
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
        "TCAGTCGATT", // "AATCGACTGA"
        "TTTTTTTTTTTTTTTTT" // "AAAAAAAAAAAAAAAAA"
    };
    const std::vector<std::string> labels { "A", "B", "C", "D" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels, true);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    MaskedDBGAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") } }} },
        { std::string("TTTTTATTTTTTTATTT"), {{ { std::string("D"), std::string("TTTTTTTTTTTTTTTTT") } }} }
    }};

    for (const auto &[query, targets] : exp_alignments) {
        const auto &alignments = aligner.align(query);

        EXPECT_EQ(targets.size(), alignments.size()) << query;
        check_alignment_labels(*anno_graph, alignments);

        for (const auto &alignment : alignments) {
            for (const auto &label : alignment.get_labels()) {
                auto find = targets.find(label);
                ASSERT_TRUE(find != targets.end());

                EXPECT_EQ(find->second, alignment.get_sequence())
                    << query << " " << label;
            }
        }
    }
}

} // namespace
