#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/annotation_converters.hpp"
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

inline std::vector<std::string> get_alignment_labels(const AnnotatedDBG &anno_graph,
                                                     const Alignment &alignment,
                                                     bool check_full_coverage = true) {
    const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
    auto labels = anno_graph.get_labels(alignment.get_sequence(),
                                        check_full_coverage ? 1.0 : 0.0);
    if (check_full_coverage) {
        EXPECT_GE(labels.size(), alignment.label_columns.size());
    }

    std::unordered_set<uint64_t> enc_labels;
    for (const auto &label : labels) {
        enc_labels.emplace(label_encoder.encode(label));
    }

    std::vector<std::string> dec_labels;
    for (uint64_t label : alignment.label_columns) {
        EXPECT_TRUE(enc_labels.count(label)) << alignment;
        dec_labels.emplace_back(label_encoder.decode(label));
    }

    return dec_labels;
}

template <typename GraphAnnotationPair>
class LabeledAlignerTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowDiffColumnAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledAlignerTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledAlignerTest, SimpleLinearGraph) {
    size_t k = 4;
    /*
        A    A    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCT-GCTT
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCTT"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("GCAATGCTT"), {{ { std::string("B"), std::string("AATGCTT") }, // 2S7=
                                       { std::string("A"), std::string("GCAAT") } //5=4S
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

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

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::string("GAATGCAT") }, // 1S8=
                                       { std::string("B"), std::string("CGAATGCCT") }, // 7=1X1=
                                       { std::string("A"), std::string("TGCCT") } // 4S3=1X1=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoords) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

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
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("CGAATGCAT"), {{ { std::string("C"), std::make_pair(std::string("GAATGCAT"), 1) }, // 1S8=
                                       { std::string("B"), std::make_pair(std::string("CGAATGCCT"), 0) }, // 7=1X1=
                                       { std::string("A"), std::make_pair(std::string("TGCCT"), 0) } // 4S3=1X1=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoordsMiddle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

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
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("CGAAAGCCT"), {{ { std::string("C"), std::make_pair(std::string("GAATGCAT"), 1) }, // 1S3=1X2=1X1=
                                       { std::string("B"), std::make_pair(std::string("CGAATGCCT"), 0) }, // 4=1X4=
                                       { std::string("A"), std::make_pair(std::string("GCCT"), 1) } // 5S4=
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
        }
    }
}

TYPED_TEST(LabeledAlignerTest, SimpleTangleGraphCoordsCycle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 4;
    /*
        A    A    B    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCG-GCGC-CGCA
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCGCA"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int32_t>>> exp_alignments {{
        { std::string("ATGCGCAATGCG"), {{ { std::string("B"), std::make_pair(std::string("ATGCGCA"), 1) },
                                          { std::string("A"), std::make_pair(std::string("GCAAT"), 0) },
                                     }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            ASSERT_EQ(alignment.label_columns.size(), alignment.label_coordinates.size());
            size_t label_index = 0;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                ASSERT_GT(alignment.label_coordinates[label_index].size(), 0);
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end()) << label;
                if (alignment.get_sequence() == find->second.first) {
                    found = true;
                    EXPECT_EQ(find->second.second,
                              alignment.label_coordinates[label_index][0]);
                    break;
                }
                ++label_index;
            }
            EXPECT_TRUE(found) << alignment;
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

    auto anno_graph = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.min_seed_length = 2;
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
        { std::string("TGAAATGCAT"), {{
#if ! _PROTEIN_GRAPH
            { std::string("C"), std::string("TGGAATGCAT") }, // 2=1X7=
            { std::string("B"), std::string("TCGAATGCCT") } // 1=2X5=1X1=
#else
            { std::string("C"), std::string("AATGCAT") }, // 3S7=
            { std::string("B"), std::string("AATGCCT") } // 3S5=1X1=
#endif
        }} }
    }};

    for (const auto &[query, labels] : exp_alignments) {
        auto alignments = aligner.align(query);
        EXPECT_EQ(labels.size(), alignments.size()) << query;

        for (const auto &alignment : alignments) {
            bool found = false;
            for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                auto find = labels.find(label);
                ASSERT_TRUE(find != labels.end());
                if (alignment.get_sequence() == find->second) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << alignment;
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

        DBGAlignerConfig config;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
        LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
            // r.c. TTTGAACTAA
            { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") } }} } // 5=2X3=
        }};

        for (const auto &[query, labels] : exp_alignments) {
            auto alignments = aligner.align(query);
            EXPECT_EQ(labels.size(), alignments.size()) << query;

            for (const auto &alignment : alignments) {
                bool found = false;
                for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
                    auto find = labels.find(label);
                    ASSERT_TRUE(find != labels.end());
                    if (alignment.get_sequence() == find->second) {
                        found = true;
                        break;
                    }
                }
                EXPECT_TRUE(found) << alignment;
            }
        }
    }
}
#endif

} // namespace
