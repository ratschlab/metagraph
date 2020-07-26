#include "gtest/gtest.h"

#include "../graph/all/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "common/threads/threading.hpp"
#include "graph/annotated_graph_algorithm.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::test;

template <typename GraphAnnotationPair>
class MaskedDeBruijnGraphAlgorithm : public ::testing::Test {};
TYPED_TEST_SUITE(MaskedDeBruijnGraphAlgorithm, GraphAnnotationCanonicalPairTypes);

// TYPED_TEST(MaskedDeBruijnGraphTest, CallUnitigsMaskTangle) {
//     size_t k = 4;
//     // TTGC      GCACGGGTC
//     //      TGCA
//     // ATGC      GCAGTGGTC
//     std::vector<std::string> sequences { "TTGCACGGGTC", "ATGCAGTGGTC" };
//     const std::vector<std::string> labels { "A", "B" };
//     auto anno_graph = build_anno_graph<TypeParam,
//                                        annot::ColumnCompressed<>>(
//         k, sequences, labels
//     );

//     auto masked_dbg = build_masked_graph(*anno_graph, { "A" }, {});
//     std::unordered_multiset<std::string> ref = { "TTGCACGGGTC" };
//     std::unordered_multiset<std::string> obs;

//     masked_dbg.call_unitigs([&](const auto &unitig, const auto &path) {
//         ASSERT_EQ(path, map_sequence_to_nodes(masked_dbg, unitig));
//         obs.insert(unitig);
//     });

//     EXPECT_EQ(obs, ref);
// }

template <class Graph, class Annotation = annot::ColumnCompressed<>>
void test_mask_indices(double density_cutoff) {
    for (size_t num_threads = 1; num_threads < 5; num_threads += 3) {
        set_num_threads(num_threads);
        const std::vector<std::string> ingroup { "B", "C" };
        const std::vector<std::string> outgroup { "A" };

        for (size_t k = 3; k < 15; ++k) {
            const std::vector<std::string> sequences {
                std::string("T") + std::string(k - 1, 'A') + std::string(100, 'T'),
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "C",
                std::string("T") + std::string(k - 1, 'A') + "A",
                std::string("T") + std::string(k - 1, 'A') + "G"
            };
            const std::vector<std::string> labels { "A", "B", "C", "D", "E" };

            auto anno_graph = build_anno_graph<Graph, Annotation>(k, sequences, labels);

            std::unordered_set<std::string> obs_labels, obs_kmers;
            const std::unordered_set<std::string> ref_kmers {
                std::string(k - 1, 'A') + "C"
            };
            const std::unordered_set<std::string> ref_labels {
                "B", "C"
            };

            auto masked_dbg = build_masked_graph(*anno_graph,
                                                 ingroup,
                                                 outgroup,
                                                 1.0,
                                                 0.0,
                                                 0.0,
                                                 density_cutoff);

            // FYI: num_nodes() throws exception for masked graph with lazy node mask
            // EXPECT_EQ(anno_graph->get_graph().num_nodes(), masked_dbg.num_nodes());
            ASSERT_EQ(anno_graph->get_graph().max_index(), masked_dbg.max_index());

            masked_dbg.call_kmers([&](auto i, const auto &kmer) {
                auto cur_labels = anno_graph->get_labels(i);
                obs_labels.insert(cur_labels.begin(), cur_labels.end());
                obs_kmers.insert(kmer);
            });

            EXPECT_EQ(ref_labels, obs_labels) << k << " " << density_cutoff;
            EXPECT_EQ(ref_kmers, obs_kmers) << k << " " << density_cutoff;
        }
    }
    set_num_threads(1);
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabel) {
    for (double d = 0.0; d <= 1.0; d += 0.05) {
        test_mask_indices<typename TypeParam::first_type,
                          typename TypeParam::second_type>(d);
    }
}


template <class Graph, class Annotation = annot::ColumnCompressed<>>
void
test_mask_unitigs(double inlabel_fraction,
                  double outlabel_fraction,
                  double other_label_fraction,
                  const std::unordered_set<std::string> &ref_kmers) {
    for (size_t num_threads = 1; num_threads < 5; num_threads += 3) {
        set_num_threads(num_threads);
        const std::vector<std::string> ingroup { "B", "C" };
        const std::vector<std::string> outgroup { "A" };
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

            auto masked_dbg = make_masked_graph_by_unitig_labels(*anno_graph,
                                                                 ingroup,
                                                                 outgroup,
                                                                 inlabel_fraction,
                                                                 outlabel_fraction,
                                                                 other_label_fraction);

            EXPECT_EQ(anno_graph->get_graph().max_index(), masked_dbg.max_index());

            masked_dbg.call_kmers([&](auto, const auto &kmer) { obs_kmers.insert(kmer); });

            EXPECT_EQ(ref_kmers, obs_kmers)
                << k << " "
                << inlabel_fraction << " "
                << outlabel_fraction << " "
                << other_label_fraction;
        }
    }
    set_num_threads(1);
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabel) {
    for (double other_frac : std::vector<double>{ 0.0, 1.0 }) {
        std::unordered_set<std::string> ref_kmers;
        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 0.0, other_frac, ref_kmers);

        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 0.24, other_frac, ref_kmers);

        ref_kmers.insert("GAA");
        ref_kmers.insert("AAT");
        ref_kmers.insert("ATG");
        ref_kmers.insert("TGC");

        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 0.25, other_frac, ref_kmers);


        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 0.50, other_frac, ref_kmers);

        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 0.75, other_frac, ref_kmers);

        test_mask_unitigs<typename TypeParam::first_type,
                          typename TypeParam::second_type>(1.0, 1.0, other_frac, ref_kmers);
    }
}


template <class Graph, class Annotation = annot::ColumnCompressed<>>
void
test_mask_unitigs_canonical(double inlabel_fraction,
                            double outlabel_fraction,
                            double other_label_fraction,
                            const std::unordered_set<std::string> &ref_kmers,
                            bool add_canonical = false) {
    for (size_t num_threads = 1; num_threads < 5; num_threads += 3) {
        set_num_threads(num_threads);
        const std::vector<std::string> ingroup { "B", "C" };
        const std::vector<std::string> outgroup { "A" };
        size_t k = 5;

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
                k, sequences, labels, !add_canonical
            );

            std::unordered_set<std::string> obs_kmers;

            auto masked_dbg = make_masked_graph_by_unitig_labels(
                *anno_graph,
                ingroup, outgroup,
                inlabel_fraction, outlabel_fraction,
                other_label_fraction,
                add_canonical
            );

            if (!add_canonical) {
                EXPECT_EQ(anno_graph->get_graph().max_index(), masked_dbg.max_index());
            }

            masked_dbg.call_kmers([&](auto, const auto &kmer) { obs_kmers.insert(kmer); });

            EXPECT_EQ(ref_kmers, obs_kmers)
                << k << " "
                << inlabel_fraction << " "
                << outlabel_fraction << " "
                << other_label_fraction << " "
                << add_canonical << " "
                << num_threads;
        }
    }
    set_num_threads(1);
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabelCanonical) {
    for (bool add_canonical : std::vector<bool>{ false, true }) {
        for (double other_frac : std::vector<double>{ 1.0, 0.0 }) {
            std::unordered_set<std::string> ref_kmers;
            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.0, other_frac, ref_kmers, add_canonical
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.25, other_frac, ref_kmers, add_canonical
            );


            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.49, other_frac, ref_kmers, add_canonical
            );

            ref_kmers.insert("GAATG");
            ref_kmers.insert("AATGC");
            ref_kmers.insert("GCATT");
            ref_kmers.insert("CATTC");

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.50, other_frac, ref_kmers, add_canonical
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.75, other_frac, ref_kmers, add_canonical
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 1.0, other_frac, ref_kmers, add_canonical
            );
        }
    }
}

} // namespace
