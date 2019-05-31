#include "dbg_aligner.hpp"

#include <gtest/gtest.h>

#include "boss.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "annotate_column_compressed.hpp"
#include "reverse_complement.hpp"

typedef DBGSuccinct Graph;

TEST(dbg_aligner, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto path = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(0ull, path.size());
}

TEST(dbg_aligner, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(1ull, path.front().size());
    EXPECT_EQ("CAT", path.front().get_sequence());
}

TEST(dbg_aligner, inexact_seeding) {
    size_t k = 3;
    std::string reference = "CAT"  "TGTTTT";
    std::string query =     "CCCC" "TGTTTT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(4ull, path.front().size());
    EXPECT_EQ("TGTTTT", path.front().get_sequence());
}

TEST(dbg_aligner, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

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

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, align_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGA" "AT" "ATTTGTT";
    std::string reference_2 = "AGCTTCGA" "CG" "ATTTGTT";
    std::string query = reference_2;

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =       "AGGGGG";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(query, path.front().get_sequence());
}

TEST(dbg_aligner, variation) {
    size_t k = 4;
    std::string reference = "AGCAA" "C" "TCGAAA";
    std::string query =     "AGCAA" "T" "TCGAAA";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference, path.front().get_sequence());
}

TEST(dbg_aligner, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAA" "CTCGAAA";
    std::string reference_2 = "TTAAGCAA" "GTCGAAA";
    std::string query =       "TTAAGCAA" "TGGGAAA";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

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
    std::string reference = "ACGCAA" "C" "TCTCTG" "A" "A" "C" "TTGT";
    std::string query =     "ACGCAA" "T" "TCTCTG" "T" "A" "T" "TTGT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference, path.front().get_sequence());
}

TEST(dbg_aligner, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAA" "CTTTTTT";
    std::string reference_2 = "AAAA" "TTGGGGG";
    std::string query =       "AAAA" "TTTTTTT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ(query.size() - k + 1, path.front().size());
    EXPECT_EQ(reference_1, path.front().get_sequence());
    EXPECT_EQ((query.size() - 1) * aligner.get_match_score() - 1, path.front().get_total_score());
    EXPECT_EQ("4=1X6=", path.front().get_cigar());
}

TEST(dbg_aligner, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAA" "TTTT" "TTTT",
                                           "ACAA" "TTTT" "TGTT",
                                           "ACAA" "GTTT" "TTTT",
                                           "ACAA" "GTTT" "TGTT"};
    std::string query =                    "ACAA" "CTTT" "TCTT";

    auto graph = std::make_shared<Graph>(k);
    for (auto reference : references)
        graph->add_sequence(reference);

    size_t max_num_alternative_paths = 4;
    DBGAligner aligner(graph, 10, max_num_alternative_paths);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(max_num_alternative_paths, alt_paths.size());

    size_t min_expected_score = (query.size() - 4) * aligner.get_match_score();
    for (size_t i = 0; i < alt_paths.size(); ++i) {
        auto path = alt_paths[i];
        EXPECT_TRUE(path.front().get_total_score() >= min_expected_score);
        std::cerr << "i: " << i << " score: " << path.front().get_total_score()
            << " path sequence: " << path.front().get_sequence() << std::endl;
    }
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
    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner (graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

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
    std::string reference = "AAAAAAAAAA" "TTTTTTTTTTAAAAATTTTTATATATAATTAATTAA";
    reference +=            "TTTAAATTTTATTATTAATTAAAATATTTAGGTGTGGGGGGGGG" "CCCCCCCCCC" "CGCGCGCGCGC";
    std::string query =     "AAAAAAAAAA";
    query +=                "CCCCCCCCCC";
    query +=                "CGCGCGCGCG";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    DBGAligner aligner(graph);
    auto alt_paths = aligner.align(std::begin(query), std::end(query));

    EXPECT_EQ(1ull, alt_paths.size());
    auto path = alt_paths.front();

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
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCC"     "TTGTT";
    std::string query =     "TTTCC" "A" "TTGTT";

    auto graph = std::make_shared<Graph>(k);
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
    std::string reference = "TTCGGA"     "TATGGAC";
    std::string query =     "TTCGGA" "C" "TATGGAC";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ("TTCGGACTNNGGAC", path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_delete) {
    size_t k = 4;
    std::string reference = "TTCGA" "T" "TGGC";
    std::string query =     "TTCGA"     "TGGC";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
}

TEST(dbg_aligner, map_to_nodes_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATA" "CCTT" "GGCGCTCTC";
    std::string query =     "TTTCTGTATA"        "GGCGCTCTC";

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ("10=4D9=", path.get_cigar());
}

TEST(dbg_aligner, map_to_nodes_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(query.size() * aligner.get_match_score(), path.get_total_score());
}

TEST(dbg_aligner, map_to_nodes_inexact_seed) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    // Query starts with two unmappable kmers.
    // The rest of the query is the same as the reference.
    std::string query =     "TT";
    query += reference;

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes(query);

    EXPECT_EQ(query.size() - 2 - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(reference.size() * aligner.get_match_score(), path.get_total_score());
}

TEST(dbg_aligner, map_to_nodes_reverse_complement) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;
    std::string rev_comp_query(query);
    reverse_complement(rev_comp_query.begin(), rev_comp_query.end());

    auto graph = std::make_shared<Graph>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);
    DBGAligner aligner(graph);
    auto path = aligner.map_to_nodes_forward_reverse_complement(rev_comp_query);

    EXPECT_EQ(rev_comp_query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(rev_comp_query.size() * aligner.get_match_score(), path.get_total_score());
}
