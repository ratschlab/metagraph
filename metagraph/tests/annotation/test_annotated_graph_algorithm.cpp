#include "gtest/gtest.h"

#include "../graph/all/test_dbg_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/masked_graph.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::annot;
using namespace mtg::test;


template <typename GraphAnnotationPair>
class MaskedDeBruijnGraphAlgorithm : public ::testing::Test {};
// test with DBGBitmap and DBGHashFast to have k-mers in different order
typedef ::testing::Types<std::pair<DBGBitmap, ColumnCompressed<>>,
                         std::pair<DBGHashFast, ColumnCompressed<>>,
                         std::pair<DBGBitmap, RowFlatAnnotator>,
                         std::pair<DBGHashFast, RowFlatAnnotator>
                        > GraphAnnoTypes;
TYPED_TEST_SUITE(MaskedDeBruijnGraphAlgorithm, GraphAnnoTypes);

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskIndicesByLabel) {
    typedef typename TypeParam::first_type Graph;
    typedef typename TypeParam::second_type Annotation;
    for (size_t num_threads = 1; num_threads < 5; num_threads += 3) {
        const std::vector<std::string> ingroup { "B", "C" };
        const std::vector<std::string> outgroup { "A" };

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

            std::unordered_set<std::string> obs_labels, obs_kmers;
            const std::unordered_set<std::string> ref_kmers {
                std::string(k - 1, 'A') + "C"
            };
            const std::unordered_set<std::string> ref_labels {
                "B", "C"
            };

            DifferentialAssemblyConfig config {
                .label_mask_in_unitig_fraction = 0.0,
                .label_mask_in_kmer_fraction = 1.0,
                .label_mask_out_unitig_fraction = 1.0,
                .label_mask_out_kmer_fraction = 0.0,
                .label_mask_other_unitig_fraction = 1.0
            };

            auto masked_dbg = mask_nodes_by_label(*anno_graph,
                                                  ingroup, outgroup,
                                                  {}, {},
                                                  config, num_threads);

            masked_dbg.call_kmers([&](auto, const std::string &kmer) {
                auto cur_labels = anno_graph->get_labels(kmer, 0.0);
                obs_labels.insert(cur_labels.begin(), cur_labels.end());
                obs_kmers.insert(kmer);
            });

            EXPECT_EQ(ref_labels, obs_labels) << k;
            EXPECT_EQ(ref_kmers, obs_kmers) << k;
        }
    }
}

template <class Graph, class Annotation = ColumnCompressed<>>
void test_mask_unitigs(double inlabel_fraction,
                       double outlabel_fraction,
                       double other_label_fraction,
                       const std::unordered_set<std::string> &ref_kmers) {
    for (size_t num_threads : {1, 4}) {
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

            DifferentialAssemblyConfig config {
                .label_mask_in_unitig_fraction = inlabel_fraction,
                .label_mask_in_kmer_fraction = 1.0,
                .label_mask_out_unitig_fraction = outlabel_fraction,
                .label_mask_out_kmer_fraction = 0.0,
                .label_mask_other_unitig_fraction = other_label_fraction
            };

            auto masked_dbg = mask_nodes_by_label(*anno_graph,
                                                  ingroup, outgroup,
                                                  {}, {},
                                                  config, num_threads);

            masked_dbg.call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

            EXPECT_EQ(ref_kmers, obs_kmers)
                << k << " "
                << inlabel_fraction << " "
                << outlabel_fraction << " "
                << other_label_fraction;
        }
    }
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

#if ! _PROTEIN_GRAPH
template <class Graph, class Annotation = ColumnCompressed<>>
void test_mask_unitigs_canonical(double inlabel_fraction,
                                 double outlabel_fraction,
                                 double other_label_fraction,
                                 const std::unordered_set<std::string> &ref_kmers,
                                 bool add_complement,
                                 DeBruijnGraph::Mode mode) {
    for (size_t num_threads = 1; num_threads < 5; num_threads += 3) {
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

            // if add_complement, then the differential
            auto anno_graph = build_anno_graph<Graph, Annotation>(
                k, sequences, labels, mode
            );

            std::unordered_set<std::string> obs_kmers;

            DifferentialAssemblyConfig config {
                .label_mask_in_unitig_fraction = inlabel_fraction,
                .label_mask_in_kmer_fraction = 1.0,
                .label_mask_out_unitig_fraction = outlabel_fraction,
                .label_mask_out_kmer_fraction = 0.0,
                .label_mask_other_unitig_fraction = other_label_fraction,
                .add_complement = add_complement
            };

            auto masked_dbg = mask_nodes_by_label(*anno_graph,
                                                  ingroup, outgroup,
                                                  {}, {},
                                                  config, num_threads);

            masked_dbg.call_kmers([&](auto, const std::string &kmer) { obs_kmers.insert(kmer); });

            EXPECT_EQ(ref_kmers, obs_kmers)
                << k << " "
                << inlabel_fraction << " "
                << outlabel_fraction << " "
                << other_label_fraction << " "
                << add_complement << " "
                << mode << " "
                << num_threads;
        }
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabelCanonical) {
    std::vector<std::pair<bool, DeBruijnGraph::Mode>> params {
        { true, DeBruijnGraph::BASIC },
        { false, DeBruijnGraph::CANONICAL },
        { true, DeBruijnGraph::CANONICAL },
        { false, DeBruijnGraph::PRIMARY },
        { true, DeBruijnGraph::PRIMARY },
    };
    for (const auto &[add_complement, mode] : params) {
        for (double other_frac : std::vector<double>{ 1.0, 0.0 }) {
            std::unordered_set<std::string> ref_kmers;
            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.0, other_frac, ref_kmers, add_complement, mode
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.25, other_frac, ref_kmers, add_complement, mode
            );


            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.49, other_frac, ref_kmers, add_complement, mode
            );

            ref_kmers.insert("CATTC");
            ref_kmers.insert("GCATT");
            if (mode == DeBruijnGraph::CANONICAL || add_complement) {
                ref_kmers.insert("GAATG");
                ref_kmers.insert("AATGC");
            }

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.50, other_frac, ref_kmers, add_complement, mode
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 0.75, other_frac, ref_kmers, add_complement, mode
            );

            test_mask_unitigs_canonical<typename TypeParam::first_type,
                                        typename TypeParam::second_type>(
                1.0, 1.0, other_frac, ref_kmers, add_complement, mode
            );
        }
    }
}
#endif

} // namespace
