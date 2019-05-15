#include "dbg_aligner.hpp"

#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "annotate_column_compressed.hpp"

TEST(dbg_aligner, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
    auto path = aligner.align(query);

    EXPECT_EQ(0ull, path.size());
}

TEST(dbg_aligner, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(1ull, path.front().size());
    EXPECT_EQ("CAT", path.front().get_sequence());
}

TEST(dbg_aligner, inexact_seeding) {
    size_t k = 3;
    std::string reference = "CATTGTTTT";
    std::string query =     "CCCCTGTTTT";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =       "AGGGGG";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
    auto path = aligner.align(query);

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, variation) {
    size_t k = 4;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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
    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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

    DBGSuccinct* graph = new DBGSuccinct(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph, new annotate::ColumnCompressed<>(/*num_rows=*/graph->num_nodes() + 1));
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
