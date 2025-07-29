#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_match.hpp"
#include "graph/alignment_redone/aligner_config.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/annotated_dbg.hpp"
#include "kmer/alphabets.hpp"
#include "seq_io/sequence_io.hpp"

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align_redone;
using namespace mtg::test;
using namespace mtg::kmer;

// inline std::vector<std::string> get_alignment_labels(const AnnotatedDBG &anno_graph,
//                                                      const Alignment &alignment,
//                                                      bool check_full_coverage = true) {
//     const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
//     auto labels = anno_graph.get_labels(alignment.get_path_spelling(),
//                                         check_full_coverage ? 1.0 : 0.0);
//     const auto &columns = alignment.get_label_class();
//     if (check_full_coverage) {
//         EXPECT_GE(labels.size(), columns.size());
//     }

//     std::unordered_set<uint64_t> enc_labels;
//     for (const auto &label : labels) {
//         enc_labels.emplace(label_encoder.encode(label));
//     }

//     std::vector<std::string> dec_labels;
//     for (uint64_t label : columns) {
//         EXPECT_TRUE(enc_labels.count(label)) << alignment;
//         dec_labels.emplace_back(label_encoder.decode(label));
//     }

//     return dec_labels;
// }

void run_alignment(const AnnotatedDBG &anno_graph,
                   DBGAlignerConfig config,
                   std::string_view query,
                   const std::vector<std::tuple<std::string, std::string, std::string>> &mappings,
                   size_t end_trim = 0,
                   bool needs_extension = false,
                   bool needs_extension_long_seed = false) {
    const auto &graph = anno_graph.get_graph();
    AnnotationBuffer anno_buffer(graph, anno_graph.get_annotator());

    size_t k = graph.get_k();
    if (config.min_seed_length == 0)
        config.min_seed_length = k;

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        config.max_seed_length = std::max(mx, config.min_seed_length);
        bool check_chaining = mx == std::numeric_limits<size_t>::max()
            ? !needs_extension_long_seed
            : !needs_extension;

        Query aln_query(graph, query);
        LabeledSeeder seeder(anno_buffer, aln_query, config);
        std::vector<Alignment> paths;
        LabeledExtender extender(anno_buffer, aln_query, config);
        std::vector<Alignment> paths_no_extend = seeder.get_inexact_anchors();
        for (const auto &base_path : paths_no_extend) {
            extender.extend(base_path, [&](Alignment&& path) {
                paths.emplace_back(std::move(path));
            });
        }
        std::sort(paths_no_extend.begin(), paths_no_extend.end(), [](const auto &a, const auto &b) {
            return std::make_pair(a.get_score(), b.get_orientation())
                 > std::make_pair(b.get_score(), a.get_orientation());
        });
        std::sort(paths.begin(), paths.end(), [](const auto &a, const auto &b) {
            return std::make_pair(a.get_score(), b.get_orientation())
                 > std::make_pair(b.get_score(), a.get_orientation());
        });


        ASSERT_LE(mappings.size(), paths.size()) << mx;
        paths.resize(mappings.size());

        for (size_t i = 0; i < mappings.size(); ++i) {
            const auto &[label, reference, cigar_str] = mappings[i];
            Cigar cigar(cigar_str);

            auto check_ref = [&](const Alignment &path, const std::string &type) {
                const auto &[label, reference, cigar_str] = mappings[i];
                if (reference.size()) {
                    EXPECT_EQ(reference.size() - k + 1 + end_trim, path.get_path().size()) << mx << "\t" << type;
                    EXPECT_EQ(reference, path.get_spelling()) << mx << "\t" << type;
                }

                EXPECT_EQ(end_trim, path.get_end_trim()) << mx << "\t" << type;
                EXPECT_EQ(cigar.to_string(), path.get_cigar().to_string()) << mx << "\t" << type;
                if (reference.size()) {
                    EXPECT_EQ(config.score_cigar(reference, query, cigar), path.get_score()) << mx << "\t" << type;
                }

                EXPECT_EQ(cigar.get_clipping(), path.get_clipping()) << mx << "\t" << type;
                EXPECT_EQ(cigar.get_end_clipping(), path.get_end_clipping()) << mx << "\t" << type;
            };

            check_ref(paths[i], "extend");

            if (check_chaining
                    && cigar.data()[0].first == Cigar::MATCH && cigar.data()[0].second >= config.min_seed_length
                    && cigar.data().back().first == Cigar::MATCH && cigar.data().back().second >= config.min_seed_length) {
                // this alignment should work with chaining alone
                ASSERT_LT(i, paths_no_extend.size());
                check_ref(paths_no_extend[i], "chain");
            }
        }
    }
}

template <typename GraphAnnotationPair>
class LabeledAlignerRedoneTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledAlignerRedoneTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledAlignerRedoneTest, SimpleLinearGraph) {
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
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("GCAATGCTT"), {
            std::make_tuple(std::string("B"), std::string("AATGCTT"), std::string("2S7=")),
            std::make_tuple(std::string("A"), std::string("GCAAT"), std::string("5=4S")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings);
    }
}

// TYPED_TEST(LabeledAlignerTest, SimpleTangleGraph) {
//     size_t k = 3;
//     /*  B                  AB  AB
//        CGA                 GCC-CCT
//           \ BC  BC  BC ABC/
//            GAA-AAT-ATG-TGC
//          C/               \  C   C
//        GGA                 GCA-CAT
//     */
//     const std::vector<std::string> sequences {
//         "TGCCT",
//         "CGAATGCCT",
//         "GGAATGCAT"
//     };
//     const std::vector<std::string> labels { "A", "B", "C" };

//     auto anno_graph = build_anno_graph<typename TypeParam::first_type,
//                                        typename TypeParam::second_type>(k, sequences, labels);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
//     LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

//     std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
//         { std::string("CGAATGCAT"), {{ { std::string("C"), std::string("GAATGCAT") }, // 1S8=
//                                        { std::string("B"), std::string("CGAATGCCT") }, // 7=1X1=
//                                        { std::string("A"), std::string("TGCCT") } // 4S3=1X1=
//                                      }} }
//     }};

//     for (const auto &[query, labels] : exp_alignments) {
//         auto alignments = aligner.align(query);
//         EXPECT_EQ(labels.size(), alignments.size()) << query;

//         for (const auto &alignment : alignments) {
//             bool found = false;
//             for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
//                 auto find = labels.find(label);
//                 ASSERT_TRUE(find != labels.end()) << label;
//                 if (alignment.get_sequence() == find->second) {
//                     found = true;
//                     break;
//                 }
//             }
//             EXPECT_TRUE(found) << alignment;
//         }
//     }
// }

// TEST(LabeledAlignerTest, SimpleTangleGraphSuffixSeed) {
//     size_t k = 4;
//     /*  B    B                  AB   AB
//        TCGA-CGAA                TGCC-GCCT
//                 \ BC   BC   BC /
//                  GAAT-AATG-ATGC
//           C    C/              \   C    C
//        TGGA-GGAA                TGCA-GCAT
//     */
//     const std::vector<std::string> sequences {
//         "TGCCT",
//         "TCGAATGCCT",
//         "TGGAATGCAT"
//     };
//     const std::vector<std::string> labels { "A", "B", "C" };

//     auto anno_graph = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, sequences, labels);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
//     config.min_seed_length = 2;
//     config.left_end_bonus = 5;
//     config.right_end_bonus = 5;
//     LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

//     std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
//         { std::string("TGAAATGCAT"), {{
// #if ! _PROTEIN_GRAPH
//             { std::string("C"), std::string("TGGAATGCAT") }, // 2=1X7=
//             { std::string("B"), std::string("TCGAATGCCT") } // 1=2X5=1X1=
// #else
//             { std::string("C"), std::string("AATGCAT") }, // 3S7=
//             { std::string("B"), std::string("AATGCCT") } // 3S5=1X1=
// #endif
//         }} }
//     }};

//     for (const auto &[query, labels] : exp_alignments) {
//         auto alignments = aligner.align(query);
//         EXPECT_EQ(labels.size(), alignments.size()) << query;

//         for (const auto &alignment : alignments) {
//             bool found = false;
//             for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
//                 auto find = labels.find(label);
//                 ASSERT_TRUE(find != labels.end());
//                 if (alignment.get_sequence() == find->second) {
//                     found = true;
//                     break;
//                 }
//             }
//             EXPECT_TRUE(found) << alignment;
//         }
//     }
// }

// #if ! _PROTEIN_GRAPH
// TYPED_TEST(LabeledAlignerTest, CanonicalTangleGraph) {
//     size_t k = 5;
//     /*   B     B                AB    AB
//        TTAGT-TAGTC             TCGAA-CGAAA
//                   \  BC   ABC /
//                    AGTCG-GTCGA
//           C     C /           \   C     C
//        TCAGT-CAGTC             TCGAT-CGATT
//         AB    AB                 B     B
//        TTTCG-TTCGA             GACTA-ACTAA
//                   \ ABC    BC /
//                    TCGAC-CGACT
//           C     C /           \   C     C
//        AATCG-ATCGA             GACTG-ACTGA
//     */
//     const std::vector<std::string> sequences {
//         "GTCGAAA", // "TTTCGAC"
//         "TTAGTCGAAA", // "TTTCGACTAA"
//         "TCAGTCGATT" // "AATCGACTGA"
//     };
//     const std::vector<std::string> labels { "A", "B", "C" };

//     for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
//         auto anno_graph = build_anno_graph<typename TypeParam::first_type,
//                                            typename TypeParam::second_type>(
//             k, sequences, labels, mode
//         );

//         DBGAlignerConfig config;
//         config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//         LabeledAligner<> aligner(anno_graph->get_graph(), config, anno_graph->get_annotator());

//         std::unordered_map<std::string, std::unordered_map<std::string, std::string>> exp_alignments {{
//             // r.c. TTTGAACTAA
//             { std::string("TTAGTTCAAA"), {{ { std::string("B"), std::string("TTAGTCGAAA") } }} } // 5=2X3=
//         }};

//         for (const auto &[query, labels] : exp_alignments) {
//             auto alignments = aligner.align(query);
//             EXPECT_EQ(labels.size(), alignments.size()) << query;

//             for (const auto &alignment : alignments) {
//                 bool found = false;
//                 for (const auto &label : get_alignment_labels(*anno_graph, alignment)) {
//                     auto find = labels.find(label);
//                     ASSERT_TRUE(find != labels.end());
//                     if (alignment.get_sequence() == find->second) {
//                         found = true;
//                         break;
//                     }
//                 }
//                 EXPECT_TRUE(found) << alignment;
//             }
//         }
//     }
// }
// #endif

} // namespace