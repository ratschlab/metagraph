#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "test_aligner_helpers.hpp"

#include "graph/alignment/dbg_aligner.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;

template <typename Graph>
class DBGAlignerPostChainTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGSuccinct,
                         DBGSuccinctUnitigIndexed,
                         DBGSuccinctPathIndexed> ChainGraphTypes;
TYPED_TEST_SUITE(DBGAlignerPostChainTest, ChainGraphTypes);

inline void check_chain(const AlignmentResults &paths,
                        const DeBruijnGraph &graph,
                        const DBGAlignerConfig &config,
                        bool has_chain = true) {
    for (const auto &path : paths) {
        EXPECT_TRUE(path.is_valid(graph, &config)) << path;
        if (has_chain) {
            EXPECT_THROW(path.to_json(graph.get_k(), false, "", ""), std::runtime_error);
        } else {
            check_json_dump_load(graph, path, paths.get_query(), paths.get_query(true));
        }
    }
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_swap) {
    size_t k = 11;
    std::string reference = "ATGATATGAGGGGGGGGGGGGTTTTTTTTGACCCCGGTTTAA";
    std::string query     = "TTTTTTTTGACCCCGGTTTAAATGATATGAGGGGGGGGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = k;
    config.seed_complexity_filter = false;
    config.allow_jump = true;
    config.set_node_insertion_penalty(k);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("TTTTTTTTGACCCCGGTTTAA$ATGATATGAGGGGGGGGGGGG"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_overlap_2) {
    size_t k = 9;
    std::string reference1 = "CCCCCCTTTGAGGATCAG";
    std::string reference2 =          "CCGGATCAGCTAGCTAGCTAGC";
    std::string query      = "CCCCCCTTTGAGGATCAGCTAGCTAGCTAGC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.allow_jump = true;
    config.min_seed_length = 7;
    config.max_seed_length = 7;
    config.seed_complexity_filter = false;
    config.set_node_insertion_penalty(k);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("CCCCCCTTTGAGGATCAGCTAGCTAGCTAGC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_overlap_mismatch) {
    size_t k = 8;
    std::string reference1 = "TTTTTCCTGAGGATCCG";
    std::string reference2 =        "CCCGGATCAGCTAGCTAGCTAGC";
    std::string query      = "TTTTTCCTGAGGATCTGCTAGCTAGCTAGC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.allow_jump = true;
    config.forward_and_reverse_complement = true;
    config.min_seed_length = 5;
    config.max_seed_length = 5;
    config.set_node_insertion_penalty(k);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("TTTTTCCTGAGGATCAGCTAGCTAGCTAGC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_overlap_3_prefer_mismatch_over_gap) {
    size_t k = 11;
    std::string reference1 = "AAATTTTGAGGATCAG";
    std::string reference2 =      "CCCCGGATCAGGTTTATTTAATTAGCT";
    std::string reference3 =                      "CCCCATTAGCTTGCTAGCAAAAA";
    std::string query      = "AAATTTTGAGGATCAGCTTTATTTAATTAGCTTGCTAGCAAAAA";
    //                                        X

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2, reference3 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.min_seed_length = 7;
    config.max_seed_length = 7;
    config.seed_complexity_filter = false;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("AAATTTTGAGGATCAGGTTTATTTAATTAGCTTGCTAGCAAAAA"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_delete_no_chain_if_full_coverage) {
    size_t k = 10;
    std::string reference = "TGAGGATCAGTTCTAGCTTGCTAGC";
    std::string query     = "TGAGGATCAG""CTAGCTTGCTAGC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = k;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(reference, paths[0].get_sequence());
    check_chain(paths, *graph, config, false);
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_delete1) {
    size_t k = 10;
    std::string reference1 = "TTTTGAGGATCAGTTCTAGCTTG";
    std::string reference2 =              "CCCTAGCTTGCTAGCGCTAGCTAGATC";
    std::string query      = "TTTTGAGGATCAG""CTAGCTTGCTAGCGCTAGCTAGATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = 8;
    config.max_seed_length = 8;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("TTTTGAGGATCAGTTCTAGCTTGCTAGCGCTAGCTAGATC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_delete_mismatch) {
    size_t k = 10;
    std::string reference1 = "AAAAAGGGTTTTTGAGGATCAGTTCTGCGCTTG";
    std::string reference2 =                       "CCCTACGCTTGCTAGCGCTAGCTAGATC";
    std::string query      = "AAAAAGGGTTTTTGAGGATCAG""CTTCGCTTGCTAGCGCTAGCTAGATC";
    //                                                  X

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = 6;
    config.max_seed_length = 6;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("AAAAAGGGTTTTTGAGGATCAGTTCTGCGCTTGCTAGCGCTAGCTAGATC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_overlap_with_insert) {
    size_t k = 10;
    std::string reference1 =  "TGAGGATCAGTTCTAGCTTG";
    std::string reference2 =            "CCCTAGCTTGCTAGCGCTAGCTAGATC";
    std::string query      = "TGAGGATCAGTTCTGAGCTTGCTAGCGCTAGCTAGATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    config.min_seed_length = 6;
    config.max_seed_length = 6;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(1, -1, -1);
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("TGAGGATCAGTTCTAGCTTGCTAGCGCTAGCTAGATC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}


TYPED_TEST(DBGAlignerPostChainTest, align_chain_deletion_in_overlapping_node) {
    size_t k = 10;
    std::string reference1 = "AAATTTTTTTGAGGATCAGTTCTAAGCTTG";
    std::string reference2 =                     "CCCCAGCTTGCTAGCGCTAGCTAGATC";
    std::string query      = "AAATTTTTTTGAGGATCAG""CTAAGCTTGCTAGCGCTAGCTAGATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.allow_jump = true;
    config.min_seed_length = 5;
    config.max_seed_length = 5;
    config.set_node_insertion_penalty(k);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("AAATTTTTTTGAGGATCAGTTCTAAGCTTGCTAGCGCTAGCTAGATC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_large_overlap) {
    size_t k = 10;
    std::string reference1 = "TGAGGATCAGTTCTAGCTTG";
    std::string reference2 =      "ATCAGTTCTAGCTTGCTAGCGCTAGCTAGATC";
    std::string query      = "TGAGGATCAGTAATCTAGCTTGCTAGCGCTAGCTAGATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = k;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("TGAGGATCAGTTCTAGCTTGCTAGCGCTAGCTAGATC"), paths[0].get_sequence());
    check_chain(paths, *graph, config, false);
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_delete_in_overlap) {
    size_t k = 10;
    std::string reference1 = "AAATTTTTTTGAGGATCAGTTCTAGCTTGC";
    std::string reference2 =                     "AATAGCTTGCTAGCGCTAGCTAGATCCCCCCC";
    std::string query      = "AAATTTTTTTGAGGATCAGTTCTGCTTGCTAGCGCTAGCTAGATCCCCCCC";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = 6;
    config.max_seed_length = 6;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("AAATTTTTTTGAGGATCAGTTCTAGCTTGCTAGCGCTAGCTAGATCCCCCCC"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_disjoint) {
    size_t k = 10;
    std::string reference1 = "GGGGGGGGGGAAACCCCCCCCTGAGGATCAG";
    std::string reference2 =                                "TTCACTAGCTAGCCCCCCCCCGGGGGGGGGG";
    std::string query      = "GGGGGGGGGGAAACCCCCCCCTGAGGATCAGTTCACTAGCTAGCCCCCCCCCGGGGGGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = k;
    config.seed_complexity_filter = false;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("GGGGGGGGGGAAACCCCCCCCTGAGGATCAG$TTCACTAGCTAGCCCCCCCCCGGGGGGGGGG"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerPostChainTest, align_chain_gap) {
    size_t k = 10;
    std::string reference1 = "AAAAACCCCCTGAGGATCAG";
    std::string reference2 =                        "ACTAGCTAGCCCCCCAAAAA";
    std::string query      = "AAAAACCCCCTGAGGATCAGTTCACTAGCTAGCCCCCCAAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference1, reference2 });
    DBGAlignerConfig config;
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(1, -1, -1);
    config.min_seed_length = k;
    config.set_node_insertion_penalty(k);
    config.allow_jump = true;

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1u, paths.size());
    EXPECT_EQ(std::string("AAAAACCCCCTGAGGATCAG$ACTAGCTAGCCCCCCAAAAA"), paths[0].get_sequence());
    check_chain(paths, *graph, config);
    check_extend(graph, aligner.get_config(), paths, query);
}

} // namespace
