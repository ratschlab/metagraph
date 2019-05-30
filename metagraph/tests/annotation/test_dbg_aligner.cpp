#include "dbg_aligner.hpp"

#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "annotate_column_compressed.hpp"

typedef DBGSuccinct Graph;

TEST(dbg_aligner, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(0ull, path.size());
}

TEST(dbg_aligner, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(1ull, path.front().size());
    EXPECT_EQ("CAT", path.front().get_sequence());
}

TEST(dbg_aligner, inexact_seeding) {
    size_t k = 3;
    std::string reference = "CATTGTTTT";
    std::string query =     "CCCCTGTTTT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(4ull, path.front().size());
    EXPECT_EQ("TGTTTT", path.front().get_sequence());
}

TEST(dbg_aligner, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    // Query is the same as the reference.
    std::string query =     "AGCTTCGAGGCCAA";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
    EXPECT_EQ(query.size() * aligner.get_match_score(), path.front().get_total_score());
}

TEST(dbg_aligner, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    // Query is the same as the second reference.
    std::string query =       "AGCTTCGAC";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, align_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    // Query is the same as the second reference.
    std::string query =       "AGCTTCGACGATTTGTT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =       "AGGGGG";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, variation) {
    size_t k = 4;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference, path.front().get_sequence());
}

TEST(dbg_aligner, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAACTCGAAA";
    std::string reference_2 = "TTAAGCAAGTCGAAA";
    std::string query =       "TTAAGCAATGGGAAA";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_TRUE(path.front().get_sequence().compare(reference_1) == 0 ||
                path.front().get_sequence().compare(reference_2) == 0)
        << "Path: " << path.front().get_sequence() << std::endl
        << "Ref1: " << reference_1 << std::endl
        << "Ref2: " << reference_2 << std::endl;
}

TEST(dbg_aligner, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAACTCTCTGAAC";
    std::string query =     "ACGCAATTCTCTGTAT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference, path.front().get_sequence());
}

TEST(dbg_aligner, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAACTTTTTT";
    std::string reference_2 = "AAAATTGGGGG";
    std::string query =       "AAAATTTTTTT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference_1, path.front().get_sequence());
    EXPECT_EQ((query.size() - 1) * aligner.get_match_score() - 1, path.front().get_total_score());
    EXPECT_EQ("4=1X6=", path.front().get_cigar());
}

TEST(dbg_aligner, large_search_space) {
    size_t k = 3;
    std::string reference = "";
    auto alphabet = {'A', 'G', 'T'};
    for (auto first_letter : alphabet) {
        for (auto second_letter : alphabet) {
            for (auto third_letter : alphabet) {
                reference = reference + first_letter +
                            second_letter + third_letter;
            }
        }
    }
    std::string query = "AAA";
    uint64_t unmapped_char_length = 100;
    for (size_t i = 0; i < unmapped_char_length; i++) {
        query += 'C';
    }
    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph);
    auto path = aligner.align(query);

    std::string aligned_query(query);
    std::replace(aligned_query.begin(), aligned_query.end(), 'C' , 'T');

    EXPECT_EQ(2ull, path.size());
    EXPECT_EQ(1ull, path.front().size());
    EXPECT_EQ("AAA", path.front().get_sequence());
    EXPECT_EQ(k * aligner.get_match_score(), path.front().get_total_score());
    EXPECT_EQ(0ull, path.back().size());
    EXPECT_EQ("", path.back().get_sequence());
}

TEST(dbg_aligner, large_gap) {
    size_t k = 10;
    std::string reference = "AAAAAAAAAATTTTTTTTTTAAAAATTTTTATATATAATTAATTAA";
    reference +=            "TTTAAATTTTATTATTAATTAAAATATTTAGGTGTGGGGGGGGGCCCCCCCCCCCGCGCGCGCGC";
    std::string query =     "AAAAAAAAAA";
    query +=                "CCCCCCCCCC";
    query +=                "CGCGCGCGCG";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(2ull, path.size());

    auto path_seq = path.front().get_sequence();
    EXPECT_EQ(query.substr(0, k), path_seq);
    EXPECT_EQ(std::begin(query), path.front().get_query_begin_it());
    EXPECT_EQ(k * aligner.get_match_score(), path.front().get_total_score());

    path_seq = path.back().get_sequence();
    EXPECT_EQ(query.substr(query.size() - 2 * k), path_seq);
    EXPECT_EQ(std::end(query) - 2 * k, path.back().get_query_begin_it());
    EXPECT_EQ(2 * k * aligner.get_match_score(), path.back().get_total_score());
}

TEST(dbg_aligner, map_to_nodes_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAGCGGACCCTTTCCGTTAT";
    std::string query =     "AAAGGGGACCCTTTTCGTTAT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCACTTGTT";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_insert) {
    size_t k = 4;
    // NOTE: very challenging pair of reference and query!
    // Should keep it as is. Negative overlap.
    // The problem with graphs.
    std::string reference = "TTCGGATATGGAC";
    std::string query =     "TTCGGACTATGGAC";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ("TTCGGACTNNGGAC", path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_delete) {
    size_t k = 4;
    std::string reference = "TTCGATGGC";
    std::string query =     "TTCGAGGC";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

//TEST(dbg_aligner, map_to_nodes_gap) {
//    size_t k = 4;
//    std::string reference = "TTCGATTATAATTGGCGCC";
//    std::string query =     "TTCGAGGCGCC";
//
//    Graph* graph = new Graph(k);
//    graph->add_sequence(reference);
//    graph->mask_dummy_kmers(1, false);
//    DBGAligner aligner(graph);
//    auto path = aligner.map_to_nodes(query);
//
//    EXPECT_EQ(query.size() - k + 1, path.size());
//    EXPECT_EQ(query, path.get_sequence());
//}

TEST(dbg_aligner, map_to_nodes_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    // Query is the same as the reference.
    std::string query =     "AGCTTCGAGGCCAA";

    Graph* graph = new Graph(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(query.size() * aligner.get_match_score(), path.get_total_score());
}
