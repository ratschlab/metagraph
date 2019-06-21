#include "dbg_aligner.hpp"

#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "annotate_column_compressed.hpp"
#include "reverse_complement.hpp"
#include "../graph/test_dbg_helpers.hpp"

template <typename Graph>
class DBGAlignerTest : public DeBruijnGraphTest<Graph> {};

TYPED_TEST_CASE(DBGAlignerTest, GraphTypes);

TYPED_TEST(DBGAlignerTest, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto path = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(0ull, path.size());
}

TYPED_TEST(DBGAlignerTest, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", path.get_sequence());
}

TEST(DBGAlignerTest, inexact_seeding) {
    size_t k = 4;
    std::string reference = "GGCC" "TGTTTG";
    std::string query =     "ACCC" "TGTTTG";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(7ull, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(query.size() * aligner.get_match_score(), path.get_total_score());
}

TYPED_TEST(DBGAlignerTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, align_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGA" "AT" "ATTTGTT";
    std::string reference_2 = "AGCTTCGA" "CG" "ATTTGTT";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =       "AGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, variation) {
    size_t k = 4;
    std::string reference = "AGCAA" "C" "TCGAAA";
    std::string query =     "AGCAA" "T" "TCGAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAA" "CTCGAAA";
    std::string reference_2 = "TTAAGCAA" "GTCGAAA";
    std::string query =       "TTAAGCAA" "TGGGAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_TRUE(path.get_sequence().compare(reference_1) == 0 ||
                path.get_sequence().compare(reference_2) == 0)
        << "Path: " << path.get_sequence() << std::endl
        << "Ref1: " << reference_1 << std::endl
        << "Ref2: " << reference_2 << std::endl;
}

TYPED_TEST(DBGAlignerTest, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAA" "C" "TCTCTG" "A" "A" "C" "TTGT";
    std::string query =     "ACGCAA" "T" "TCTCTG" "T" "A" "T" "TTGT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAA" "CTTTTTT";
    std::string reference_2 = "AAAA" "TTGGGGG";
    std::string query =       "AAAA" "TTTTTTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAligner aligner(graph,
            /*num_top_paths =*/ 10,
            /*num_alternative_paths =*/ 2);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_TRUE(alt_paths.size() > 0);
    auto path = alt_paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, path.get_sequence());
    EXPECT_EQ((query.size() - 1) * aligner.get_match_score() - 1, path.get_total_score());
    EXPECT_EQ("4=1X6=", path.get_cigar());
}

TYPED_TEST(DBGAlignerTest, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAA" "TTTT" "TTTT",
                                           "ACAA" "TTTT" "TGTT",
                                           "ACAA" "GTTT" "TTTT",
                                           "ACAA" "GTTT" "TGTT"};
    std::string query =                    "ACAA" "CTTT" "TCTT";

    auto graph = build_graph_batch<TypeParam>(k, references);

    size_t max_num_alternative_paths = 4;
    DBGAligner aligner(graph, 10, max_num_alternative_paths);
    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));

    EXPECT_EQ(max_num_alternative_paths, alt_paths.size());

    // TODO: how to test alternative paths?
//    for (size_t i = 0; i < alt_paths.size(); ++i) {
//        auto path = alt_paths[i];
//        EXPECT_TRUE(path.get_total_score() >= min_expected_score);
//    }
}

//TYPED_TEST(DBGAlignerTest, large_search_space) {
//    size_t k = 3;
//    std::string reference = "";
//    auto alphabet = {'A', 'G', 'T'};
//    for (auto first_letter : alphabet) {
//        for (auto second_letter : alphabet) {
//            for (auto third_letter : alphabet) {
//                reference = reference + first_letter +
//                            second_letter + third_letter;
//            }
//        }
//    }
//    std::string query = "AAA";
//    uint64_t unmapped_char_length = 100;
//    for (size_t i = 0; i < unmapped_char_length; i++) {
//        query += 'C';
//    }
//    auto graph = build_graph_batch<TypeParam>(k, { reference });
//    DBGAligner aligner (graph);
//    auto alt_paths = aligner.align_by_graph_exploration(std::begin(query), std::end(query));
//
//    EXPECT_EQ(1ull, alt_paths.size());
//    auto path = alt_paths.front();
//
//    std::string aligned_query(query);
//    std::replace(aligned_query.begin(), aligned_query.end(), 'C' , 'T');
//
//    EXPECT_EQ(1ull, path.size());
//    EXPECT_EQ("AAA", path.get_sequence());
//    EXPECT_EQ("3=", path.get_cigar().substr(0, 2));
//    EXPECT_EQ(k * aligner.get_match_score(), path.get_total_score());
////    EXPECT_EQ("", path.back().get_sequence());
//}

TYPED_TEST(DBGAlignerTest, map_to_nodes_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, map_to_nodes_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCC"     "TTGTT";
    std::string query =     "TTTCC" "A" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

//TYPED_TEST(DBGAlignerTest, map_to_nodes_inverse_stitching) {
//    size_t k = 4;
//    std::string reference = "TTCGGA"     "TATGGAC";
//    std::string query =     "TTCGGA" "C" "TATGGAC";
//
//    auto graph = build_graph_batch<TypeParam>(k, { reference });
//    DBGAligner aligner(graph);
//    auto paths = aligner.map_to_nodes(query);
//
//    EXPECT_EQ(1ull, paths.size());
//    auto path = paths.front();
//
//    EXPECT_EQ(reference, path.get_sequence());
//}

TYPED_TEST(DBGAlignerTest, map_to_nodes_delete) {
    size_t k = 4;
    std::string reference = "TTCGA" "T" "TGGC";
    std::string query =     "TTCGA"     "TGGC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TYPED_TEST(DBGAlignerTest, map_to_nodes_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATA" "CCTT" "GGCGCTCTC";
    std::string query =     "TTTCTGTATA"        "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ("10=4D9=", path.get_cigar());
}

TYPED_TEST(DBGAlignerTest, map_to_nodes_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(query.size() * aligner.get_match_score(), path.get_total_score());
}

TEST(DBGAlignerTest, map_to_nodes_inexact_seed) {
    size_t k = 4;
    std::string reference = "AAA" "AGCTTCGAGGCCAA";
    std::string query =      "TT" "AGCTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference.substr(1), path.get_sequence());
    EXPECT_EQ("14=", path.get_cigar());
}

TEST(DBGAlignerTest, map_to_nodes_inexact_seed_snp) {
    size_t k = 7;
    std::string reference = "AAAAGCTT" "TCGAGGCCAA";
    std::string query =        "ACCTT" "TCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto paths = aligner.map_to_nodes(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ("13=", path.get_cigar());
}
// TODO: Back to typed test when inexact seeding with graphs other than BOSS is handled.
TEST(DBGAlignerTest, map_to_nodes_reverse_complement) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;
    std::string rev_comp_query(query);
    reverse_complement(rev_comp_query.begin(), rev_comp_query.end());

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes_forward_reverse_complement(rev_comp_query);

    EXPECT_EQ(rev_comp_query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(rev_comp_query.size() * aligner.get_match_score(), path.get_total_score());
}
