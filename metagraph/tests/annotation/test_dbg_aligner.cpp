#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_aligner.hpp"
#include "dbg_succinct.hpp"
#include "annotate_column_compressed.hpp"

TEST(dbg_aligner, align_sequence_too_short) {
    size_t k = 4;
    std::string sequence = "CATTT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align("CAT");

    EXPECT_EQ(0ull, path.size());
    EXPECT_EQ("", aligner.get_path_sequence(path));
}

TEST(dbg_aligner, align_single_node) {
    size_t k = 3;
    std::string sequence = "CAT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align("CAT");

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", aligner.get_path_sequence(path));
}

TEST(dbg_aligner, align_straight) {
    size_t k = 4;
    std::string sequence = "AGCTTCGAGGCCAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence);

    EXPECT_EQ(sequence.size() - k + 1, path.size());
    EXPECT_EQ(sequence, aligner.get_path_sequence(path));
}

TEST(dbg_aligner, align_ending_branch) {
    size_t k = 4;
    std::string sequence_1 = "AGCTTCGAA";
    std::string sequence_2 = "AGCTTCGAC";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence_1);
    graph->add_sequence(sequence_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence_2);

    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
    EXPECT_EQ(sequence_2, aligner.get_path_sequence(path));
}

TEST(dbg_aligner, align_branch) {
    size_t k = 4;
    std::string sequence_1 = "AGCTTCGAATATTTGTT";
    std::string sequence_2 = "AGCTTCGACGATTTGTT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence_1);
    graph->add_sequence(sequence_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence_2);

    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
    EXPECT_EQ(sequence_2, aligner.get_path_sequence(path));
}

TEST(dbg_aligner, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string sequence_1 = "AGGGGGGGGGAAAAGGGGGGG";
    std::string sequence_2 = "AGGGGG";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence_1);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence_2);

    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
    EXPECT_EQ(sequence_2, aligner.get_path_sequence(path));
}

TEST(dbg_aligner, variation) {
    size_t k = 4;
    std::string sequence_1 = "AGCAACTCGAAA";
    std::string sequence_2 = "AGCAATTCGAAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence_1);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence_2);

    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
    EXPECT_EQ(sequence_1, aligner.get_path_sequence(path));
}

TEST(dbg_aligner, variation_in_branching_point) {
    size_t k = 4;
    std::string sequence_1 = "AGCAACTCGAAA";
    std::string sequence_2 = "AGCAATTCGAAA";
    std::string sequence_3 = "AGCAAGGCGAAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(sequence_1);
    graph->add_sequence(sequence_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
    auto path = aligner.align(sequence_3);

    EXPECT_EQ(sequence_3.size() - k + 1, path.size());
    EXPECT_TRUE(aligner.get_path_sequence(path).compare(sequence_1) == 0 ||
                aligner.get_path_sequence(path).compare(sequence_2) == 0);
}

//TEST(dbg_aligner, variation_in_first_kmer) {
//    size_t k = 4;
//    std::string sequence_1 = "AGCCTCGAAA";
//    std::string sequence_2 = "AGCTTCGAAA";
//
//    DBGSuccinct* graph = new DBGSuccinct(k);
//    graph->add_sequence(sequence_1);
//    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/1));
//    auto path = aligner.align(sequence_2);
//
//    EXPECT_EQ(sequence_2.size() - k + 1, path.size());
//    EXPECT_EQ(sequence_2, aligner.get_path_sequence(path));
//}
