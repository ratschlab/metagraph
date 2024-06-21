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
                         std::pair<DBGSuccinct, ColumnCompressed<>>
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

            const auto &annotation = static_cast<const Annotation&>(anno_graph->get_annotator());
            const auto &in_columns = annotation.get_matrix().data();
            std::vector<std::unique_ptr<const bit_vector>> columns;
            columns.reserve(in_columns.size());
            std::transform(in_columns.begin(), in_columns.end(), std::back_inserter(columns),
                           [&](const auto &a) {
                               return std::make_unique<const bit_vector_stat>(a->to_vector());
                           });

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
            std::vector<bool> groups = { true, false, false };
            columns.resize(groups.size());

            std::vector<std::unique_ptr<const sdsl::int_vector<>>> column_values;
            auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<sdsl::int_vector<>, std::vector<uint64_t>>(
                graph_ptr,
                columns,
                column_values,
                [&](uint64_t row_i, uint64_t col_j) -> uint64_t {
                    return (*columns[col_j])[row_i];
                },
                groups,
                config, num_threads
            );

            EXPECT_EQ(0u, masked_dbg_out->num_nodes());

            masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
                auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                obs_labels.insert(cur_labels.begin(), cur_labels.end());
                obs_kmers.insert(kmer);
            });

            EXPECT_EQ(ref_labels, obs_labels) << k;
            EXPECT_EQ(ref_kmers, obs_kmers) << k;
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabelCounts) {
    typedef typename TypeParam::first_type Graph;
    typedef typename TypeParam::second_type Annotation;
    for (size_t num_threads : { 1, 4 }) {
        size_t max_k = std::min(max_test_k<Graph>(), (size_t)15);
        for (size_t k = 3; k < max_k; ++k) {
            std::string sig = std::string(k - 1, 'A') + "C$";
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

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels, DeBruijnGraph::BASIC, true);
            auto graph_ptr = std::static_pointer_cast<const Graph>(anno_graph->get_graph_ptr());

            using BinMatType = matrix::ColumnMajor;
            using CoordMatType = matrix::TupleCSCMatrix<BinMatType>;
            const auto &annotation = static_cast<const ColumnCoordAnnotator&>(anno_graph->get_annotator());
            const auto &count_matrix = static_cast<const CoordMatType&>(annotation.get_matrix());
            const auto &in_binary_columns = static_cast<const BinMatType&>(count_matrix.get_binary_matrix()).data();
            std::vector<std::unique_ptr<const bit_vector>> columns;
            columns.reserve(in_binary_columns.size());
            std::transform(in_binary_columns.begin(), in_binary_columns.end(), std::back_inserter(columns),
                           [&](const auto &a) {
                               return std::make_unique<const bit_vector_stat>(a->to_vector());
                           });

            ASSERT_LT(0u, columns.size());

            std::vector<std::unique_ptr<const sdsl::int_vector<>>> column_values_all;
            std::vector<sdsl::int_vector<>::iterator> pos;
            column_values_all.reserve(columns.size());
            for (const auto &bin_col : columns) {
                column_values_all.emplace_back(std::make_unique<const sdsl::int_vector<>>(bin_col->num_set_bits()));
                pos.emplace_back(const_cast<sdsl::int_vector<>&>(*column_values_all.back()).begin());
            }

            std::vector<uint64_t> row_indices(columns[0]->size());
            std::iota(row_indices.begin(), row_indices.end(), 0);
            auto row_values = count_matrix.get_row_values(row_indices);
            for (size_t i = 0; i < row_indices.size(); ++i) {
                for (const auto &[j, c] : row_values[i]) {
                    *pos[j] = c;
                    ++pos[j];
                }
            }

            std::unordered_set<std::string> obs_labels, obs_kmers;
            const std::unordered_set<std::string> ref_kmers {
                std::string(k - 1, 'A') + "C"
            };
            const std::unordered_set<std::string> ref_labels {
                "B", "C"
            };

            // ingroup: B, C
            // outgroup: A
            std::vector<bool> groups { true, false, false };
            columns.resize(groups.size());

            ASSERT_GE(column_values_all.size(), groups.size());
            column_values_all.resize(groups.size());

            std::vector<std::string> test_types = {
                "poisson_exact",
                "nbinom_exact",
                "fisher",
            };

            for (const auto &test_type : test_types) {
                DifferentialAssemblyConfig config {
                    .test_type = test_type,
                    .outfbase = "",
                };

                auto [masked_dbg_in, masked_dbg_out, pvals, tmp_file] = mask_nodes_by_label_dual<sdsl::int_vector<>, std::vector<uint64_t>>(
                    graph_ptr,
                    columns,
                    column_values_all,
                    [&](uint64_t row_i, uint64_t col_j) -> uint64_t {
                        const auto &col = *columns[col_j];
                        const auto &col_vals = *column_values_all[col_j];
                        if (uint64_t r = col.conditional_rank1(row_i))
                            return col_vals[r - 1];

                        return 0;
                    },
                    groups,
                    config, num_threads, test_dump_basename, num_threads,
                    false
                );

                masked_dbg_in->call_kmers([&](auto, const std::string &kmer) {
                    auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                    obs_labels.insert(cur_labels.begin(), cur_labels.end());
                    obs_kmers.insert(kmer);
                });
                EXPECT_EQ(ref_labels, obs_labels) << k << " " << test_type;
                EXPECT_EQ(ref_kmers, obs_kmers) << k << " " << test_type;
            }
        }
    }
}

// TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabelOverlap) {
//     typedef typename TypeParam::first_type Graph;
//     typedef typename TypeParam::second_type Annotation;
//     for (size_t num_threads : { 1, 4 }) {
//         size_t k = 5;
//         const tsl::hopscotch_set<std::string> ingroup { "A", "B", "C" };
//         const tsl::hopscotch_set<std::string> outgroup { "A", "B", "C" };

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
//             .label_mask_in_unitig_fraction = 0.0,
//             .label_mask_in_kmer_fraction = 0.0,
//             .label_mask_out_unitig_fraction = 1.0,
//             .label_mask_out_kmer_fraction = 0.99,
//             .label_mask_other_unitig_fraction = 1.0,
//             .outfbase = "",
//         };

//         auto masked_dbg = mask_nodes_by_label(*anno_graph,
//                                               ingroup, outgroup,
//                                               {}, {},
//                                               config, num_threads);

//         masked_dbg->call_kmers([&](auto, const std::string &kmer) {
//             obs_kmers.insert(kmer);
//         });

//         EXPECT_EQ(ref_kmers, obs_kmers);
//     }
// }

// template <class Graph, class Annotation = ColumnCompressed<>>
// void test_mask_unitigs(double inlabel_fraction,
//                        double outlabel_fraction,
//                        double other_label_fraction,
//                        const std::unordered_set<std::string> &ref_kmers) {
//     for (size_t num_threads : {1, 4}) {
//         const tsl::hopscotch_set<std::string> ingroup { "B", "C" };
//         const tsl::hopscotch_set<std::string> outgroup { "A" };
//         size_t k = 3;

//         {
//             /*  B                  AB  AB
//                CGA                 GCC-CCT
//                   \ BC  BC  BC ABC/
//                    GAA-AAT-ATG-TGC
//                  C/               \  C   C
//                GGA                 GCA-CAC
//             */
//             const std::vector<std::string> sequences {
//                     "TGCCT",
//                 "CGAATGCCT",
//                 "GGAATGCAC",
//                 "TTTTTTTTTTTTTT"
//             };
//             const std::vector<std::string> labels { "A", "B", "C", "D" };

//             auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

//             std::unordered_set<std::string> obs_kmers;

//             DifferentialAssemblyConfig config {
//                 .label_mask_in_unitig_fraction = inlabel_fraction,
//                 .label_mask_in_kmer_fraction = 1.0,
//                 .label_mask_out_unitig_fraction = outlabel_fraction,
//                 .label_mask_out_kmer_fraction = 0.0,
//                 .label_mask_other_unitig_fraction = other_label_fraction,
//                 .outfbase = "",
//             };

//             auto masked_dbg = mask_nodes_by_label(*anno_graph,
//                                                   ingroup, outgroup,
//                                                   {}, {},
//                                                   config, num_threads);

//             masked_dbg->call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

//             EXPECT_EQ(ref_kmers, obs_kmers)
//                 << k << " "
//                 << inlabel_fraction << " "
//                 << outlabel_fraction << " "
//                 << other_label_fraction;
//         }
//     }
// }

// TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabel) {
//     for (double other_frac : std::vector<double>{ 0.0, 1.0 }) {
//         std::unordered_set<std::string> ref_kmers;
//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 0.0, other_frac, ref_kmers);

//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 0.24, other_frac, ref_kmers);

//         ref_kmers.insert("GAA");
//         ref_kmers.insert("AAT");
//         ref_kmers.insert("ATG");
//         ref_kmers.insert("TGC");

//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 0.25, other_frac, ref_kmers);

//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 0.50, other_frac, ref_kmers);

//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 0.75, other_frac, ref_kmers);

//         test_mask_unitigs<typename TypeParam::first_type,
//                           typename TypeParam::second_type>(1.0, 1.0, other_frac, ref_kmers);
//     }
// }

// #if ! _PROTEIN_GRAPH
// template <class Graph, class Annotation = ColumnCompressed<>>
// std::unordered_set<std::string>
// test_mask_unitigs_canonical(double inlabel_fraction,
//                             double outlabel_fraction,
//                             double other_label_fraction,
//                             const std::unordered_set<std::string> &ref_kmers,
//                             bool add_complement,
//                             DeBruijnGraph::Mode mode) {
//     for (size_t num_threads : { 1, 4 }) {
//         const tsl::hopscotch_set<std::string> ingroup { "B", "C" };
//         const tsl::hopscotch_set<std::string> outgroup { "A" };
//         size_t k = 5;

//         {
//             /*
//                 B                AB    AB
//                CGAAT             ATGCC-TGCCT
//                     \ BC   ABC  /
//                      GAATG-AATGC
//                  C  /           \  C     C
//                GGAAT             ATGCA-TGCAC
//                                    C  /  C                B
//                                  GTGCA-TGCAT             ATTCG
//                                             \ABC    BC  /
//                                              GCATT-CATTC
//                                  AB    AB   /           \  C
//                                  AGGCA-GGCAT             ATTCC
//             */
//             const std::vector<std::string> sequences {
//                   "AATGCCT",     // "AGGCATT"
//                 "CGAATGCCT",     // "AGGCATTCG"
//                 "GGAATGCAC",     // "GTGCATTCC"
//                 "TTTTTTTTTTTTTT" // "AAAAAAAAAAAAAA"
//             };
//             const std::vector<std::string> labels { "A", "B", "C", "D" };

//             // if add_complement, then the differential
//             auto anno_graph = build_anno_graph<Graph, Annotation>(
//                 k, sequences, labels, mode
//             );

//             std::unordered_set<std::string> obs_kmers;

//             DifferentialAssemblyConfig config {
//                 .label_mask_in_unitig_fraction = inlabel_fraction,
//                 .label_mask_in_kmer_fraction = 1.0,
//                 .label_mask_out_unitig_fraction = outlabel_fraction,
//                 .label_mask_out_kmer_fraction = 0.0,
//                 .label_mask_other_unitig_fraction = other_label_fraction,
//                 .add_complement = add_complement,
//                 .outfbase = "",
//             };

//             auto masked_dbg = mask_nodes_by_label(*anno_graph,
//                                                   ingroup, outgroup,
//                                                   {}, {},
//                                                   config, num_threads);

//             masked_dbg->call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

//             if (ref_kmers != obs_kmers)
//                 return obs_kmers;
//         }
//     }

//     return ref_kmers;
// }

// TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabelAddCanonical) {
//     for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC,
//                                       DeBruijnGraph::CANONICAL,
//                                       DeBruijnGraph::PRIMARY }) {
//         std::string mode_str = mode == DeBruijnGraph::BASIC
//             ? "BASIC"
//             : (mode == DeBruijnGraph::CANONICAL ? "CANONICAL" : "PRIMARY");

//         for (double other_frac : std::vector<double>{ 1.0, 0.0 }) {
//             std::unordered_set<std::string> ref_kmers;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.0, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.25, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;


//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.49, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             ref_kmers.insert("GAATG");
//             ref_kmers.insert("AATGC");

//             if (mode == DeBruijnGraph::CANONICAL || mode == DeBruijnGraph::PRIMARY) {
//                 ref_kmers.insert("CATTC");
//                 ref_kmers.insert("GCATT");
//             }

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.50, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.75, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 1.0, other_frac, ref_kmers, true, mode
//             )), ref_kmers) << other_frac << " " << mode_str;
//         }
//     }
// }

// TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabelCanonical) {
//     for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL,
//                                       DeBruijnGraph::PRIMARY }) {
//         std::string mode_str = mode == DeBruijnGraph::CANONICAL ? "CANONICAL" : "PRIMARY";
//         for (double other_frac : std::vector<double>{ 1.0, 0.0 }) {
//             std::unordered_set<std::string> ref_kmers;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.0, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.25, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;


//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.49, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             ref_kmers.insert("GAATG");
//             ref_kmers.insert("AATGC");
//             ref_kmers.insert("CATTC");
//             ref_kmers.insert("GCATT");

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.50, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 0.75, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;

//             EXPECT_EQ((test_mask_unitigs_canonical<typename TypeParam::first_type,
//                                                    typename TypeParam::second_type>(
//                 1.0, 1.0, other_frac, ref_kmers, false, mode
//             )), ref_kmers) << other_frac << " " << mode_str;
//         }
//     }
// }
// #endif

} // namespace
