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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
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
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, aligner.get_path_sequence(path.get_nodes()));
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
    uint32_t unmapped_char_length = 50;
    for (uint32_t i = 0; i < unmapped_char_length; i++) {
        query += 'C';
    }
    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(query);

    std::replace(query.begin(), query.end(), 'C' , 'T');

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, aligner.get_path_sequence(path.get_nodes()));
}


//TEST(dbg_aligner, variation_in_first_kmer) {
//    size_t k = 4;
//    std::string reference_1 = "AGCCTCGAAA";
//    std::string sequence_2 = "AGCTTCGAAA";
//
//    DBGSuccinct* graph = new DBGSuccinct(k);
//    graph->add_sequence(reference_1);
//    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
//    auto path = aligner.align(sequence_2);
//
//    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
//    EXPECT_EQ(sequence_2, aligner.get_path_sequence(path.get_nodes()));
//}
