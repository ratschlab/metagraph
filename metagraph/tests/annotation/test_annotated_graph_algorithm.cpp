#include "gtest/gtest.h"

#include <unordered_set>

#include "../graph/all/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/masked_graph.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/transpose.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::annot;
using namespace mtg::test;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";


template <typename GraphAnnotationPair>
class MaskedDeBruijnGraphAlgorithm : public ::testing::Test {};
// test with DBGBitmap and DBGHashFast to have k-mers in different order
typedef ::testing::Types<std::pair<DBGBitmap, ColumnCompressed<>>,
                         std::pair<DBGHashFast, ColumnCompressed<>>,
                         std::pair<DBGSuccinct, ColumnCompressed<>>,
                         std::pair<DBGSuccinct, RowDiffColumnAnnotator>
                        > GraphAnnoTypes;
TYPED_TEST_SUITE(MaskedDeBruijnGraphAlgorithm, GraphAnnoTypes);

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabel) {
    typedef typename TypeParam::first_type Graph;
    typedef typename TypeParam::second_type Annotation;
    for (size_t num_threads : { 1, 4 }) {
        size_t max_k = std::min(max_test_k<Graph>(), (size_t)15);
        for (size_t k = 3; k < max_k; ++k) {
            const std::vector<std::string> sequences {
                std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "A",
                std::string("T") + std::string(k - 1, 'A') + "G"
            };
            const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);
            auto graph_ptr = std::static_pointer_cast<const Graph>(anno_graph->get_graph_ptr());

            std::unordered_set<std::string> obs_labels, obs_kmers;
            const std::unordered_set<std::string> ref_kmers {
                std::string(k - 1, 'A') + "C"
            };
            const std::unordered_set<std::string> ref_labels {
                "B", "C"
            };

            DifferentialAssemblyConfig config {
                .test_type = "notest",
                .min_in_recurrence = 2,
                .max_out_recurrence = 0,
                .outfbase = "",
            };

            // ingroup: B, C
            // outgroup: A
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_in { "B", "C" };
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_out { "A" };

            {
                auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
                    *anno_graph,
                    labels_in, labels_out,
                    config, num_threads, test_dump_basename
                );

                masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
                    auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                    obs_labels.insert(cur_labels.begin(), cur_labels.end());
                    obs_kmers.insert(kmer);
                });

                EXPECT_EQ(ref_labels, obs_labels) << k;
                EXPECT_EQ(ref_kmers, obs_kmers) << k;
            }

            if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>) {
                // TODO
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabelTest) {
    typedef typename TypeParam::first_type Graph;
    typedef typename TypeParam::second_type Annotation;
    for (size_t num_threads : { 1, 4 }) {
        size_t max_k = std::min(max_test_k<Graph>(), (size_t)15);
        for (size_t k = 6; k < max_k; ++k) {
            const std::vector<std::string> sequences {
                std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "A",
                std::string("T") + std::string(k - 1, 'A') + "G"
            };
            const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);
            auto graph_ptr = std::static_pointer_cast<const Graph>(anno_graph->get_graph_ptr());

            std::unordered_set<std::string> obs_labels, obs_kmers;
            std::unordered_set<std::string> ref_kmers;

            DifferentialAssemblyConfig config {
                .test_type = "poisson_binom",
                .outfbase = "",
            };

            // ingroup: B, C
            // outgroup: A
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_in { "B", "C" };
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_out { "A" };

            std::unordered_set<std::string> ref_labels;

            std::for_each(labels.begin() + 2, labels.end(), [&](const auto &label) {
                if (label != "C")
                    labels_out.emplace(label);

                if (k >= 8 || label == "E" || (k >= 7 && label == "D")) {
                    ref_labels.emplace("B");
                    ref_labels.emplace("C");
                    ref_kmers.emplace(std::string(k - 1, 'A') + "C");
                }

                auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
                    *anno_graph,
                    labels_in, labels_out,
                    config, num_threads, test_dump_basename
                );

                masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
                    auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                    obs_labels.insert(cur_labels.begin(), cur_labels.end());
                    obs_kmers.insert(kmer);
                });

                EXPECT_EQ(ref_labels, obs_labels) << k << " " << label;
                EXPECT_EQ(ref_kmers, obs_kmers) << k << " " << label;
            });

            if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>) {
                // TODO
            }
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabelCounts) {
    typedef typename TypeParam::first_type Graph;
    typedef typename TypeParam::second_type Annotation;

    for (size_t num_threads : { 1, 4 }) {
        size_t max_k = std::min(max_test_k<Graph>(), (size_t)15);
        for (size_t k = 3; k < max_k; ++k) {
            std::string sig = std::string(k - 1, 'A') + "C";
            for (size_t i = 0; i < 7; ++i) {
                sig += sig;
            }

            std::vector<std::string> sequences {
                std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
                std::string("T") + sig,
                std::string("T") + sig,
                std::string("T") + std::string(k - 1, 'A') + "A",
                std::string("T") + std::string(k - 1, 'A') + "G"
            };

            std::vector<std::string> labels { "A", "B", "C", "D", "E" };
            std::unordered_set<std::string> obs_labels, obs_kmers;
            std::unordered_set<std::string> ref_kmers;
            for (size_t j = 0; j < k; ++j) {
                ref_kmers.emplace(std::string(j, 'A') + "C" + std::string(k - 1 - j, 'A'));
            }
            const std::unordered_set<std::string> ref_labels {
                "B", "C"
            };

            // ingroup: B, C
            // outgroup: A
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_in { "B", "C" };
            tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> labels_out { "A" };

            std::vector<std::string> test_types = {
                "poisson_exact",
                "nbinom_exact",
                "gnb_exact",
                "fisher",
            };

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels, DeBruijnGraph::BASIC, true);
            auto graph_ptr = std::static_pointer_cast<const Graph>(anno_graph->get_graph_ptr());

            for (const auto &test_type : test_types) {
                DifferentialAssemblyConfig config {
                    .test_type = test_type,
                    .outfbase = "",
                };

                auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
                    *anno_graph,
                    labels_in, labels_out,
                    config, num_threads, test_dump_basename
                );

                masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
                    auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                    obs_labels.insert(cur_labels.begin(), cur_labels.end());
                    obs_kmers.insert(kmer);
                });
                EXPECT_EQ(ref_labels, obs_labels) << k << " " << test_type;
                EXPECT_EQ(ref_kmers, obs_kmers) << k << " " << test_type;
            }

            if constexpr(std::is_same_v<Annotation, ColumnCompressed<>>) {
                // TODO
            }
        }
    }
}

// TODO: tests not defined when there is overlap in the groups
// TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabelOverlap) {
//     typedef typename TypeParam::first_type Graph;
//     typedef typename TypeParam::second_type Annotation;
//     for (size_t num_threads : { 1, 4 }) {
//         size_t k = 5;
//         const tsl::hopscotch_set<std::string> labels_in { "A", "B", "C" };
//         const tsl::hopscotch_set<std::string> labels_out { "A", "B", "C" };

//         const std::vector<std::string> sequences {
//             "TTTTTAAAAATTTTTTTTTT",
//             "CCCCCAAAAACCCCCCCCCC",
//             "GGGGGAAAAAGGGGGGGGGG"
//         };
//         const std::vector<std::string> labels { "A", "B", "C" };

//         auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

//         std::unordered_set<std::string> obs_kmers;
//         const std::unordered_set<std::string> ref_kmers {
//             "TTTTT", "TTTTA", "TTTAA", "TTAAA", "TAAAA", "AAAAT", "AAATT", "AATTT", "ATTTT",
//             "CCCCC", "CCCCA", "CCCAA", "CCAAA", "CAAAA", "AAAAC", "AAACC", "AACCC", "ACCCC",
//             "GGGGG", "GGGGA", "GGGAA", "GGAAA", "GAAAA", "AAAAG", "AAAGG", "AAGGG", "AGGGG",
//         };
//         DifferentialAssemblyConfig config {
//             .test_type = "poisson_binom",
//             .outfbase = "",
//         };

//         auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
//             *anno_graph,
//             labels_in, labels_out,
//             config, num_threads, test_dump_basename
//         );

//         masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
//             obs_kmers.insert(kmer);
//         });

//         EXPECT_EQ(ref_kmers, obs_kmers);
//     }
// }

template <class Graph, class Annotation = ColumnCompressed<>>
void test_mask_unitigs(const std::unordered_set<std::string> &ref_kmers) {
    for (size_t num_threads : {1, 4}) {
        const tsl::hopscotch_set<std::string> labels_in { "B", "C" };
        const tsl::hopscotch_set<std::string> labels_out { "A" };
        size_t k = 3;

        {
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

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

            std::unordered_set<std::string> obs_kmers;

            DifferentialAssemblyConfig config {
                .test_by_unitig = true,
                .test_type = "poisson_binom",
                .outfbase = "",
            };

            auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
                *anno_graph,
                labels_in, labels_out,
                config, num_threads, test_dump_basename
            );

            masked_dbg_in->call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

            EXPECT_EQ(ref_kmers, obs_kmers) << k;
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabel) {
    std::unordered_set<std::string> ref_kmers;
    test_mask_unitigs<typename TypeParam::first_type, typename TypeParam::second_type>(ref_kmers);

    // TODO: make a test where some k-mers are significant
}

#if ! _PROTEIN_GRAPH
template <class Graph, class Annotation = ColumnCompressed<>>
std::unordered_set<std::string>
test_mask_unitigs_canonical(const std::unordered_set<std::string> &ref_kmers,
                            DeBruijnGraph::Mode mode,
                            size_t k) {
    for (size_t num_threads : { 1, 4 }) {
        const tsl::hopscotch_set<std::string> labels_in { "B", "C" };
        const tsl::hopscotch_set<std::string> labels_out { "A" };

        {
            /*
                B                AB    AB
               CGAAT             ATGCC-TGCCT
                    \ BC   ABC  /
                     GAATG-AATGC
                 C  /           \  C     C
               GGAAT             ATGCA-TGCAC
                                   C  /  C                B
                                 GTGCA-TGCAT             ATTCG
                                            \ABC    BC  /
                                             GCATT-CATTC
                                 AB    AB   /           \  C
                                 AGGCA-GGCAT             ATTCC
            */
            const std::vector<std::string> sequences {
                  "AATGCCT",     // "AGGCATT"
                "CGAATGCCT",     // "AGGCATTCG"
                "GGAATGCAC",     // "GTGCATTCC"
                "TTTTTTTTTTTTTT" // "AAAAAAAAAAAAAA"
            };
            const std::vector<std::string> labels { "A", "B", "C", "D" };

            auto anno_graph = build_anno_graph<Graph, Annotation>(
                k, sequences, labels, mode
            );

            std::unordered_set<std::string> obs_kmers;

            DifferentialAssemblyConfig config {
                .test_by_unitig = true,
                .test_type = "poisson_binom",
                .outfbase = "",
            };

            auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<std::vector<uint64_t>>(
                *anno_graph,
                labels_in, labels_out,
                config, num_threads, test_dump_basename
            );

            masked_dbg_in->call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

            if (ref_kmers != obs_kmers)
                return obs_kmers;
        }
    }

    return ref_kmers;
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabelAddCanonical) {
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC,
                                      DeBruijnGraph::PRIMARY }) {
        std::string mode_str = mode == DeBruijnGraph::BASIC
            ? "BASIC"
            : (mode == DeBruijnGraph::CANONICAL ? "CANONICAL" : "PRIMARY");

        std::unordered_set<std::string> ref_kmers;

        EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
                                                typename TypeParam::second_type>(
            ref_kmers, mode, 5
        )), ref_kmers) << mode_str;
        EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
                                                typename TypeParam::second_type>(
            ref_kmers, mode, 6
        )), ref_kmers) << mode_str;
        EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
                                                typename TypeParam::second_type>(
            ref_kmers, mode, 7
        )), ref_kmers) << mode_str;

        // TODO: make a test where some k-mers are significant
    }
}

#endif

} // namespace
