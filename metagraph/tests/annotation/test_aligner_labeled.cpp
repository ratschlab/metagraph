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
                     const LabeledAligner<>::DBGAlignment &alignment) {
    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    auto labels = anno_graph.get_labels(alignment.get_sequence(), 1.0);
    EXPECT_GE(labels.size(), alignment.target_columns.size());

    std::unordered_set<uint64_t> enc_labels;
    for (const auto &label : labels) {
        enc_labels.emplace(label_encoder.encode(label));
    }

    std::vector<std::string> dec_labels;
    for (uint64_t target : alignment.target_columns) {
        EXPECT_TRUE(enc_labels.count(target))
            << alignment << " " << label_encoder.decode(target);
        dec_labels.emplace_back(label_encoder.decode(target));
    }

    return dec_labels;
}

template <typename GraphAnnotationPair>
class LabeledAlignerTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledAlignerTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraph) {
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
    LabeledAligner<> aligner(*anno_graph, config);

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
            EXPECT_TRUE(found) << alignment << " " << alignment.target_columns.size();
        }
    }
}

TEST(LabeledAlignerTest, SimpleTangleGraphSuffixSeed) {
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
    LabeledAligner<> aligner(*anno_graph, config);

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("TGAAATGCAT"), {{ { std::string("C"), std::string("TGGAATGCAT") },
                                        { std::string("A"), std::string("TGCCT") },
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
            EXPECT_TRUE(found) << alignment << " " << alignment.target_columns.size();
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(LabeledAlignerTest, CanonicalTangleGraph) {
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
        LabeledAligner<> aligner(*anno_graph, config);

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
            // r.c. TTTGAACTAA
            { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") },
                                            { std::string("A"), std::string("GTCGAAA") },
                                            { std::string("C"), std::string("AGTCGA") } }} }
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
                EXPECT_TRUE(found) << alignment << " " << alignment.target_columns.size();
            }
        }
    }
}
#endif

} // namespace
