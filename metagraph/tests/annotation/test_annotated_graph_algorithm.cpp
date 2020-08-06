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
TYPED_TEST_SUITE(MaskedDeBruijnGraphAlgorithm, GraphAnnotationPairTypes);

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
    const std::vector<std::string> ingroup { "B", "C" };
    const std::vector<std::string> outgroup { "A" };
    size_t k = 3;

    {
        /*
           CGA                 GCC-CCT
              \               /
               GAA-AAT-ATG-TGC
              /               \
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

        MaskedDeBruijnGraph masked_dbg(
            std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph->get_graph_ptr()),
            mask_nodes_by_unitig_labels(
                *anno_graph,
                ingroup,
                outgroup,
                inlabel_fraction,
                outlabel_fraction,
                other_label_fraction
            )
        );

        EXPECT_EQ(anno_graph->get_graph().max_index(), masked_dbg.max_index());

        masked_dbg.call_kmers([&](auto, const auto &kmer) { obs_kmers.insert(kmer); });

        EXPECT_EQ(ref_kmers, obs_kmers)
            << k << " "
            << inlabel_fraction << " "
            << outlabel_fraction << " "
            << other_label_fraction;
    }
}

TYPED_TEST(MaskedDeBruijnGraphAlgorithm, MaskUnitigsByLabel) {
    std::unordered_set<std::string> ref_kmers;

    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 0.0, 0.0, ref_kmers);

    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 0.24, 0.0, ref_kmers);

    ref_kmers.insert("GAA");
    ref_kmers.insert("AAT");
    ref_kmers.insert("ATG");
    ref_kmers.insert("TGC");

    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 0.25, 0.0, ref_kmers);


    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 0.50, 0.0, ref_kmers);

    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 0.75, 0.0, ref_kmers);

    test_mask_unitigs<typename TypeParam::first_type,
                      typename TypeParam::second_type>(1.0, 1.0, 0.0, ref_kmers);
}

} // namespace
