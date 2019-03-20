#include "dbg_aligner.hpp"

#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_aligner.hpp"
#include "dbg_succinct.hpp"
#include "annotate_column_compressed.hpp"

TEST(dbg_aligner, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query = "CAT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(0ull, path.size());
    EXPECT_EQ("", aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query = "CAT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    // Query is the same as the reference.
    std::string query = "AGCTTCGAGGCCAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    // Query is the same as the second reference.
    std::string query = "AGCTTCGAC";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, align_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    // Query is the same as the second reference.
    std::string query = "AGCTTCGACGATTTGTT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference_1 = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query = "AGGGGG";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, variation) {
    size_t k = 4;
    std::string reference_1 = "AGCAACTCGAAA";
    std::string query = "AGCAATTCGAAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AGCAACTCGAAA";
    std::string reference_2 = "AGCAAGTCGAAA";
    std::string query = "AGCAATGCGAAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_TRUE(aligner.get_path_sequence(path.get_nodes()).compare(reference_1) == 0 ||
                aligner.get_path_sequence(path.get_nodes()).compare(reference_2) == 0);
}

TEST(dbg_aligner, multiple_variations) {
    size_t k = 4;
    std::string reference_1 = "AGCAACTCGAAA";
    std::string query = "AGCAATTTGCAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, aligner.get_path_sequence(path.get_nodes()));
}

TEST(dbg_aligner, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAAACTTTTTT";
    std::string reference_2 = "AAAAATTGGGGG";
    std::string query = "AAAAATTTTTTT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, aligner.get_path_sequence(path.get_nodes()));
    EXPECT_EQ(1, path.get_total_loss());
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
    uint64_t unmapped_char_length = 50;
    for (size_t i = 0; i < unmapped_char_length; i++) {
        query += 'C';
    }
    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(query);

    std::string aligned_query(query);
    std::replace(aligned_query.begin(), aligned_query.end(), 'C' , 'T');

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(aligned_query, aligner.get_path_sequence(path.get_nodes()));
    EXPECT_EQ(unmapped_char_length, path.get_total_loss());
}
