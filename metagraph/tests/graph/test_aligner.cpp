#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "test_aligner_helpers.hpp"
#include "../test_helpers.hpp"

#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_methods.hpp"
#include "graph/representation/succinct/dbg_succinct_range.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;

const std::string test_data_dir = "../tests/data";

typedef DBGAligner<>::score_t score_t;

int8_t single_char_score(const DBGAlignerConfig &config, char a, int8_t b) {
    return config.get_row(a)[b];
}

void check_score_matrix(const DBGAlignerConfig &config,
                        const char* alphabet,
                        size_t alph_size) {
    ASSERT_EQ(strlen(alphabet), alph_size);

    for (size_t i = 0; i < alph_size; ++i) {
        if (i + 1 != alph_size) {
            EXPECT_LT(int8_t(0), single_char_score(config, alphabet[i], alphabet[i]));
        }

        for (size_t j = 0; j < alph_size; ++j) {
            //check if the match score is the greatest
            if (i + 1 != alph_size) {
                EXPECT_GE(single_char_score(config, alphabet[i], alphabet[i]),
                          single_char_score(config, alphabet[i], alphabet[j]));
            }

            // checking symmetry
            EXPECT_EQ(single_char_score(config, alphabet[i], alphabet[j]),
                      single_char_score(config, alphabet[j], alphabet[i]));
        }
    }
}

TEST(DBGAlignerTest, check_score_matrix_dna) {
    check_score_matrix(DBGAlignerConfig(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2)),
                       alphabets::kAlphabetDNA5,
                       alphabets::kSigmaDNA5);
}

TEST(DBGAlignerTest, check_score_matrix_protein) {
    check_score_matrix(DBGAlignerConfig(DBGAlignerConfig::score_matrix_blosum62),
                       alphabets::kAlphabetProtein,
                       alphabets::kSigmaProtein);
}

TEST(DBGAlignerTest, check_score_matrix_dna_unit) {
    check_score_matrix(
        DBGAlignerConfig(DBGAlignerConfig::unit_scoring_matrix(
            1, alphabets::kAlphabetDNA, alphabets::kCharToDNA
        )),
        alphabets::kAlphabetDNA5,
        alphabets::kSigmaDNA5
    );
}

TEST(DBGAlignerTest, check_score_matrix_protein_unit) {
    check_score_matrix(
        DBGAlignerConfig(DBGAlignerConfig::unit_scoring_matrix(
            1, alphabets::kAlphabetProtein, alphabets::kCharToProtein
        )),
        alphabets::kAlphabetProtein,
        alphabets::kSigmaProtein
    );
}


DBGAligner<>::DBGQueryAlignment
get_extend(std::shared_ptr<const DeBruijnGraph> graph,
           const DBGAlignerConfig &config,
           const DBGAligner<>::DBGQueryAlignment &paths,
           const std::string &query) {
    assert(graph.get());
    EXPECT_EQ(query, paths.get_query());
    auto uniconfig = config;
    uniconfig.max_seed_length = std::numeric_limits<size_t>::max();

    auto range_graph = std::dynamic_pointer_cast<const DBGSuccinctRange>(graph);
    return DBGAligner<UniMEMSeeder<>>(range_graph ? *range_graph : *graph,
                                      uniconfig).align(query);
}

void check_extend(std::shared_ptr<const DeBruijnGraph> graph,
                  const DBGAlignerConfig &config,
                  const DBGAligner<>::DBGQueryAlignment &paths,
                  const std::string &query) {
    auto unimem_paths = get_extend(graph, config, paths, query);

    ASSERT_EQ(paths.size(), unimem_paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        EXPECT_EQ(paths[i], unimem_paths[i]) << paths[i] << "\n" << unimem_paths[i];
    }
}


template <typename Graph>
class DBGAlignerTest : public DeBruijnGraphTest<Graph> {};

TYPED_TEST_SUITE(DBGAlignerTest, FewGraphTypes);

TYPED_TEST(DBGAlignerTest, bad_min_cell_score) {
    auto graph = build_graph_batch<TypeParam>(3, {});
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min();
    config.min_path_score = std::numeric_limits<score_t>::min();
    ASSERT_THROW(DBGAligner<>(*graph, config), std::runtime_error);
}

TYPED_TEST(DBGAlignerTest, align_sequence_much_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(0ull, paths.size());
}

TYPED_TEST(DBGAlignerTest, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(0ull, paths.size());
}

TYPED_TEST(DBGAlignerTest, align_big_self_loop) {
    size_t k = 3;
    std::string reference = "AAAA";
    std::string query =     "AAAAAAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(7ull, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("9=", path.get_cigar().to_string());
    EXPECT_EQ(9u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("3=", path.get_cigar().to_string());
    EXPECT_EQ(3u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight_min_path_score) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_path_score = 100;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    EXPECT_TRUE(paths.empty()) << paths.size() << "\t" << paths.front();

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight_with_N) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = "AGCTNCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(reference, query), path.get_score());
    EXPECT_EQ("4=1X9=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight_forward_and_reverse_complement) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;
    reverse_complement(query.begin(), query.end());

    auto graph = build_graph_batch<TypeParam>(k, { reference });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    auto config_fwd_and_rev = config;
    config_fwd_and_rev.forward_and_reverse_complement = true;

    DBGAligner<> aligner(*graph, config_fwd_and_rev);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);

    auto ext_paths = DBGAligner<UniMEMSeeder<>>(*graph, config_fwd_and_rev).align(query);

    EXPECT_TRUE(std::equal(paths.begin(), paths.end(),
                           ext_paths.begin(), ext_paths.end()));

    // test copy
    auto paths_copy = const_cast<const DBGAligner<>::DBGQueryAlignment&>(paths);
    for (const auto &path : paths_copy) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }

    // test move
    auto paths_move = std::move(paths);
    for (const auto &path : paths_move) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }
}


TYPED_TEST(DBGAlignerTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("9=", path.get_cigar().to_string());
    EXPECT_EQ(9u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_branch) {
    size_t k = 6;
    std::string reference_1 = "AGCTTCGA" "AT" "ATTTGTT";
    std::string reference_2 = "AGCTTCGA" "CG" "ATTTGTT";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("17=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_branch_with_cycle) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGA" "AT" "ATTTGTT";
    std::string reference_2 = "AGCTTCGA" "CG" "ATTTGTT";
    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("17=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =     "AGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("6=", path.get_cigar().to_string());
    EXPECT_EQ(6u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, variation) {
    size_t k = 4;
    std::string reference = "AGCAA" "C" "TCGAAA";
    std::string query =     "AGCAA" "T" "TCGAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_EQ("5=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(11u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAA" "CTC" "GAAA";
    std::string reference_2 = "TTAAGCAA" "GTC" "GAAA";
    std::string query =       "TTAAGCAA" "TGG" "GAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_TRUE(path.get_sequence().compare(reference_1) == 0 ||
                path.get_sequence().compare(reference_2) == 0)
        << "Path: " << path.get_sequence() << std::endl
        << "Ref1: " << reference_1 << std::endl
        << "Ref2: " << reference_2 << std::endl;
    // TODO: what about other cases?
    EXPECT_EQ("8=3X4=", path.get_cigar().to_string());
    EXPECT_EQ(12u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAA" "C" "TCTCTG" "A" "A" "C" "TTGT";
    std::string query =     "ACGCAA" "T" "TCTCTG" "T" "A" "T" "TTGT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_EQ("6=1X6=1X1=1X4=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAA" "CTTTTTT";
    std::string reference_2 = "AAAA" "TTGGGGG";
    std::string query =       "AAAA" "TTTTTTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    config.num_alternative_paths = 2;
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    ASSERT_EQ(2u, paths.size());
    EXPECT_NE(paths.front(), paths.back());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference_1), path.get_score());
    EXPECT_EQ("4=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAA" "TTTT" "TTTT",
                                           "ACAA" "TTTT" "TGTT",
                                           "ACAA" "GTTT" "TTTT",
                                           "ACAA" "GTTT" "TGTT"};
    std::string query =                    "ACAA" "CTTT" "TCTT";

    auto graph = build_graph_batch<TypeParam>(k, references);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    config.num_alternative_paths = 2;
    config.queue_size = 100;
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    EXPECT_EQ(config.num_alternative_paths, paths.size());
    for (const auto &path : paths) {
        EXPECT_EQ("4=1X4=1X2=", path.get_cigar().to_string())
            << query << "\n" << path.get_sequence();
        EXPECT_EQ(10u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());
    }

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_EQ("4=1X9=1X6=", path.get_cigar().to_string());
    EXPECT_EQ(19u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_multiple_misalignment_bandwidth) {
    size_t k = 4;
    std::string reference = "AAAG" "C" "GGACCCTTT" "C" "CGTTAT";
    std::string query =     "AAAG" "G" "GGACCCTTT" "T" "CGTTAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));

    for (uint64_t bandwidth : std::vector<uint64_t>{ 2, 5, 10, std::numeric_limits<uint64_t>::max()}) {
        auto config_bandwidth = config;
        config_bandwidth.bandwidth = bandwidth;

        DBGAligner<> aligner(*graph, config_bandwidth);
        auto paths = aligner.align(query);

        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(query.size() - k + 1, path.size());
        EXPECT_EQ(reference, path.get_sequence());
        EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
        EXPECT_EQ("4=1X9=1X6=", path.get_cigar().to_string());
        EXPECT_EQ(19u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        check_extend(graph, aligner.get_config(), paths, query);
    }
}

TYPED_TEST(DBGAlignerTest, align_delete_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCC"     "TTGTT";
    std::string query =     "TTTCC" "A" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference) + config.gap_opening_penalty, path.get_score());
    EXPECT_EQ("5=1D5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_delete_multi) {
    size_t k = 4;
    std::string reference = "TTTCC"      "TTGTT";
    std::string query =     "TTTCC" "AA" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference)
        + config.gap_opening_penalty + config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("5=2D5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_delete_long) {
    size_t k = 4;
    std::string reference = "TTTCC"             "TTGTT";
    std::string query =     "TTTCC" "AAAAAAAAA" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference) + config.gap_opening_penalty
        + score_t(8) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("5=9D5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_delete_long_offset) {
    size_t k = 4;
    std::string reference = "TTTCCG" "G"             "TTGTTA";
    std::string query =     "TTTCCG" "C" "AAAAAAAAA" "TTGTTA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(reference, "TTTCCGCTTGTTA")
        + config.gap_opening_penalty
        + score_t(8) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "6=1X9D6="
        || path.get_cigar().to_string() == "6=9D1X6=") << path.get_cigar().to_string();
    EXPECT_EQ(12u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_insert) {
    size_t k = 4;
    std::string reference = "TTCGAT" "TGGCCT";
    std::string query =     "TTCGAT"  "GGCCT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty, path.get_score());

    // TODO: the first should ideally always be true
    EXPECT_TRUE("6=1I5=" == path.get_cigar().to_string()
        || "5=1I6=" == path.get_cigar().to_string());
    // EXPECT_EQ("6=1I5=", path.get_cigar().to_string());

    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    // TODO: enable this when the above is correct
    // check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATA" "CCTT" "GGCGCTCTC";
    std::string query =     "TTTCTGTATA"        "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(3) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("10=4I9=", path.get_cigar().to_string());
    EXPECT_EQ(19u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTC" "CCTT" "GGCGCTCTC";
    std::string query =     "TTTC"        "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(3) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("4=4I9=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_loop_insertion) {
    size_t k = 4;
    std::string reference = "AAAA" "TTTT" "CGAGGCCAA";
    std::string query =     "AAAA"        "CGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::unit_scoring_matrix(
        1, alphabets::kAlphabetDNA, alphabets::kCharToDNA
    ));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(13u, path.size());
    EXPECT_EQ("AAAATTTCGAGGCCAA", path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(2) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("4=3I9=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight_long_xdrop) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAGGCCAAGCCTGACTGATCGATGCATGCTAGCTAGTCAGTCAGCGTGAGCTAGCAT";
    std::string reference_2 = "AGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    std::string query = reference_1;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.xdrop = 30;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query), path.get_score());
    EXPECT_EQ("63=", path.get_cigar().to_string());
    EXPECT_EQ(63u, path.get_num_matches());
    EXPECT_TRUE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_drop_seed) {
    size_t k = 4;
    std::string reference = "TTTCC" "CT" "GGCGCTCTC";
    std::string query =     "TTTCC" "GG" "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -6, -6));
    config.gap_opening_penalty = -10;
    config.gap_extension_penalty = -4;
    config.xdrop = 6;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(6, path.size());
    EXPECT_EQ(reference.substr(7), path.get_sequence());
    EXPECT_EQ(config.match_score(reference.substr(7)), path.get_score());
    EXPECT_EQ("7S9=", path.get_cigar().to_string());
    EXPECT_EQ(9u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(7u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    // TODO: re-enable this when gap extension is fixed
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_long_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTC" "CCTTAA" "GGCGCTCTC";
    std::string query =     "TTTC"          "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(6, path.size());
    EXPECT_EQ(reference.substr(10), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(4)), path.get_score());
    EXPECT_EQ("4S9=", path.get_cigar().to_string());
    EXPECT_EQ(9u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(4u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_repeat_sequence_no_insert_after_delete) {
    size_t k = 27;
    std::string reference = "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATG" "A"        "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATG" "ACAAATGT" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    // the alignment may be misrepresented as
    // "TTTGTGGCTAGAGCTCGAGATCGCGCG"                    "GCCACAATT" "GACAAATG" "A" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    // "TTTGTGGCTAGAGCTCGAGATCGCGCG" "GCCACAATTGACAAAT"             "GACAAATG" "T" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(67ull, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
    ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "45=7D8=1X39="
             || path.get_cigar().to_string() == "44=2D1=5D8=1X39="
             || path.get_cigar().to_string() == "44=3D1=4D8=1X39="
             || path.get_cigar().to_string() == "44=4D1=3D8=1X39=")
        << path.get_cigar().to_string();
    EXPECT_EQ(92u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    auto extends = get_extend(graph, aligner.get_config(), paths, query);
    ASSERT_EQ(1ull, extends.size());
    path = extends.front();

    EXPECT_EQ(67ull, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
    ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "45=7D8=1X39="
             || path.get_cigar().to_string() == "44=2D1=5D8=1X39="
             || path.get_cigar().to_string() == "44=3D1=4D8=1X39="
             || path.get_cigar().to_string() == "44=4D1=3D8=1X39=")
        << path.get_cigar().to_string();
    EXPECT_EQ(92u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
}

TYPED_TEST(DBGAlignerTest, align_clipping1) {
    size_t k = 4;
    std::string reference = "GGCC" "TGTTTG";
    std::string query =     "ACCC" "TGTTTG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(5ull, path.size());
    EXPECT_EQ(reference.substr(2), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
    EXPECT_EQ("2S8=", path.get_cigar().to_string())
        << reference.substr(2) << " " << path.get_sequence();
    EXPECT_EQ(8u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_clipping2) {
    size_t k = 4;
    std::string reference = "AAA" "AGCTTCGAGGCCAA";
    std::string query =      "TT" "AGCTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(11u, path.size());
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
    EXPECT_EQ("2S14=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_long_clipping) {
    size_t k = 4;
    std::string reference = "TTTTTTT" "AAAAGCTTCGAGGCCAA";
    std::string query =     "CCCCCCC" "AAAAGCTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(14u, path.size());
    EXPECT_EQ(reference.substr(7), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(7)), path.get_score());
    EXPECT_EQ("7S17=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(7u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_end_clipping) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAA" "TTTTTTT";
    std::string query =     "AAAAGCTTCGAGGCCAA" "CCCCCCC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(14u, path.size());
    EXPECT_EQ(reference.substr(0, 17), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(0, 17)), path.get_score());
    EXPECT_EQ("17=7S", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(7u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_clipping_min_cell_score) {
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(7u, path.size());
    EXPECT_EQ(reference.substr(5), path.get_sequence());
    EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
    EXPECT_EQ("2S13=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(2u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_low_similarity) {
    size_t k = 27;
    std::string reference = "CTAGAACTTAAAGTATAATAATACTAA" "TAA" "TAAAA" "TA" "A" "AATACA";
    std::string query =     "CTAGAACTTAAAGTATAATAATACTAA" "TAA"  "AAG"  "TA" "C" "AATACA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    // EXPECT_EQ(7u, path.size());
    // EXPECT_EQ(reference.substr(5), path.get_sequence());
    // EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
    // EXPECT_EQ("2S13=", path.get_cigar().to_string());
    // EXPECT_EQ(13u, path.get_num_matches());
    // EXPECT_FALSE(path.is_exact_match());
    // EXPECT_EQ(2u, path.get_clipping());
    // EXPECT_EQ(0u, path.get_end_clipping());
    // EXPECT_EQ(0u, path.get_offset());
    // EXPECT_TRUE(path.is_valid(*graph, &config));
    // check_json_dump_load(*graph,
    //                      path,
    //                      paths.get_query(),
    //                      paths.get_query_reverse_complement());

    // check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_low_similarity2) {
    size_t k = 27;
    std::string reference = "GCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
}

TYPED_TEST(DBGAlignerTest, align_low_similarity3) {
    size_t k = 27;
    std::string reference = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCTGGGATTATAGGTGTGAACCACCACACCTGGCTAATTTTTTTTGTGTGTGTGTGTGTTTTTTC";
    std::string query =     "AAAAAAAAAAAAAAAAAAAAAAAAAAACGCCAAAAAGGGGGAATAGGGGGGGGGGAACCCCAACACCGGTATGTTTTTTTGTGTGTGGGGGATTTTTTTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
}

TYPED_TEST(DBGAlignerTest, align_low_similarity4) {
    size_t k = 6;
    std::vector<std::string> seqs;
    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
                                          [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
    auto graph = build_graph_batch<TypeParam>(k, std::move(seqs));

    std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGA";

    for (size_t xdrop : { 27, 30 }) {
        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
        config.gap_opening_penalty = -5;
        config.gap_extension_penalty = -2;
        config.xdrop = xdrop;
        config.exact_kmer_match_fraction = 0.0;
        config.max_nodes_per_seq_char = 10.0;
        config.queue_size = 20;
        config.num_alternative_paths = 2;
        config.min_path_score = 0;
        config.min_cell_score = 0;

        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);

        ASSERT_EQ(2ull, paths.size());
        EXPECT_EQ(557llu, paths[0].get_score()) << paths[0];
        EXPECT_EQ(556llu, paths[1].get_score()) << paths[1];
    }
}

TEST(DBGAlignerTest, align_low_similarity5) {
    size_t k = 31;
    std::vector<std::string> seqs;
    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/refs.fa",
                                          [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
    auto base_graph = std::dynamic_pointer_cast<DBGSuccinct>(
        build_graph_batch<DBGSuccinct>(k, std::move(seqs), DBGMode::CANONICAL)
    );
    base_graph->reset_mask();
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/reads.fa", [&](auto *seq) {
        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -4, -4));
        config.gap_opening_penalty = -5;
        config.gap_extension_penalty = -5;
        config.xdrop = 27;
        config.exact_kmer_match_fraction = 0.0;
        config.max_nodes_per_seq_char = 10.0;
        config.queue_size = 20;
        config.num_alternative_paths = 1;
        config.min_seed_length = 15;
        config.min_path_score = 0;
        config.min_cell_score = 0;

        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(seq->seq.s);

        ASSERT_EQ(1ull, paths.size()) << seq->name.s;
    });
}

TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto base_graph = std::make_shared<DBGSuccinct>(k);
    base_graph->add_sequence(reference);
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    {
        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        config.min_seed_length = 2;
        config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
        config.min_path_score = std::numeric_limits<score_t>::min() + 100;
        config.max_seed_length = k;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(7u, path.size());
        EXPECT_EQ(reference.substr(5), path.get_sequence());
        EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
        EXPECT_EQ("2S13=", path.get_cigar().to_string());
        EXPECT_EQ(13u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(2u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        check_extend(graph, aligner.get_config(), paths, query);
    }
    {
        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        config.min_seed_length = 1;
        config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
        config.min_path_score = std::numeric_limits<score_t>::min() + 100;
        config.max_seed_length = k;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(9u, path.size());
        EXPECT_EQ(reference.substr(3), path.get_sequence());
        EXPECT_EQ(config.score_sequences(query, reference.substr(3)), path.get_score());
        EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
        EXPECT_EQ(14u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());

        check_extend(graph, aligner.get_config(), paths, query);
    }
}

TEST(DBGAlignerTest, align_suffix_seed_snp) {
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto base_graph = std::make_shared<DBGSuccinct>(k);
    base_graph->add_sequence(reference);
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    config.max_seed_length = k;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(9u, path.size());
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference.substr(3)), path.get_score());
    EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_suffix_seed_snp_canonical) {
    size_t k = 18;
    std::string reference = "AAAAA" "CTTTCGAGGCCAA";
    std::string query =     "GGGGG" "CTTTCGAGGCCAA";

    auto base_graph = std::make_shared<DBGSuccinct>(k, true);
    base_graph->add_sequence(reference);
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    config.max_seed_length = k;
    config.min_seed_length = 13;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    ASSERT_EQ(1u, path.size());
    EXPECT_EQ(reference.substr(5), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.substr(5), reference.substr(5)), path.get_score());
    EXPECT_EQ("5S13=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(5u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(5u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_nodummy) {
    size_t k = 7;
    std::string reference = "AAAAG" "C" "TTTCGAGGCCAA";
    std::string query =     "AAAAG" "T" "TTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.max_seed_length = k;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(6u, path.size());
    EXPECT_EQ(reference.substr(6), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.substr(6), reference.substr(6)), path.get_score());
    EXPECT_EQ("6S12=", path.get_cigar().to_string());
    EXPECT_EQ(12u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(6u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_both_directions) {
    size_t k = 7;
    std::string reference = "AAAAG" "C" "TTTCGAGGCCAA";
    std::string query =     "AAAAG" "T" "TTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DBGMode::CANONICAL);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.max_seed_length = k;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(12u, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_seed_to_end) {
    size_t k = 5;
    std::string reference = "ATCCCTTTTAAAA";
    std::string query =     "ATCCCGGGGGGGGGGGGGGGGGTTTTAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_dummy) {
    size_t k = 7;
    std::string reference = "AAAAG" "C" "TTTCGAGGCCAA";
    std::string query =     "AAAAG" "T" "TTTCGAGGCCAA";

    auto base_graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.max_seed_length = k;
    base_graph->add_sequence(reference);
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(12u, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_extended_delete_after_match) {
    size_t k = 27;
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG"    "GCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG"    "CCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
    std::string query =       "CGTGGCCCAGGCCCAGGCCCAG" "TGGGCGTTGGCCCAGGCGGCCACGGTGGCTGCGCAGGCCCGCCTGGCACAAGCCACGCTG";
    auto base_graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.max_seed_length = k;
    config.min_seed_length = 15;
    base_graph->add_sequence(reference_1);
    base_graph->add_sequence(reference_2);
    auto graph = std::make_shared<DBGSuccinctRange>(*base_graph);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();

    EXPECT_EQ(47u, path.size());
    EXPECT_EQ(reference_1, path.get_sequence());
    EXPECT_EQ("22=3D2=3X9=3I4=1D1=1X2D3=2D1=1D7=1I1=2I2=1X3=1X6=6S",
              path.get_cigar().to_string());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_suffix_seed_snp_canonical_wrapper) {
    size_t k = 18;
    std::string reference = "TTGGCCTCGAAAG" "TTTTT";
    std::string query =     "GGGGG" "CTTTCGAGGCCAA";
    std::string ref_rev =   "AAAAA" "CTTTCGAGGCCAA";

    auto dbg_succ = std::make_shared<DBGSuccinct>(k);
    dbg_succ->add_sequence(reference);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    config.min_seed_length = 13;

    for (size_t max_seed_length : { k, k + 100}) {
        config.max_seed_length = max_seed_length;
        auto base_graph = std::make_shared<const DBGSuccinctRange>(*dbg_succ);
        auto graph = std::make_shared<CanonicalDBG>(base_graph, true);
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();

        EXPECT_EQ(1u, path.size()); // includes dummy k-mers
        EXPECT_EQ(ref_rev.substr(5), path.get_sequence());
        EXPECT_EQ(config.score_sequences(query.substr(5), ref_rev.substr(5)), path.get_score());
        EXPECT_EQ("5S13=", path.get_cigar().to_string());
        EXPECT_EQ(13u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(5u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(5u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());
    }
}

TEST(DBGAlignerTest, align_range_canonical_long) {
    size_t k = 31;
    std::vector<std::string> seqs;
    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
                                          [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
    auto base_graph = build_graph_batch<DBGSuccinct>(k, std::move(seqs));
    dynamic_cast<DBGSuccinct*>(base_graph.get())->reset_mask();
    auto range_graph = std::make_shared<DBGSuccinctRange>(base_graph);
    auto graph = std::make_shared<CanonicalDBG>(*range_graph, true);

    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/long_seq.fa", [&](auto *seq) {
        std::string query = seq->seq.s;

        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
        config.gap_opening_penalty = -5;
        config.gap_extension_penalty = -2;
        config.xdrop = 27;
        config.exact_kmer_match_fraction = 0.0;
        config.max_nodes_per_seq_char = 10.0;
        config.queue_size = 20;
        config.min_path_score = 0;
        config.min_cell_score = 0;
        config.min_seed_length = 15;
        config.max_ram_per_alignment = 50000;

        DBGAligner<UniMEMSeeder<>> aligner(*graph, config);
        auto paths = aligner.align(query);

        for (auto &path : paths) {
            EXPECT_TRUE(path.is_valid(*graph, &config));
            check_json_dump_load(*graph,
                                 path,
                                 paths.get_query(),
                                 paths.get_query_reverse_complement());
        }
    });
}

TEST(DBGAlignerTest, align_canonical_nonprimary) {
    size_t k = 31;
    std::vector<std::string> seqs;
    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
                                          [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
    auto base_graph = build_graph_batch<DBGSuccinct>(k, std::move(seqs));
    auto graph = std::make_shared<CanonicalDBG>(*base_graph, true);

    std::string query = "ACACCTGTAATCCCAGCACTTTGGGAGGCCGA";

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.gap_opening_penalty = -5;
    config.gap_extension_penalty = -2;
    config.xdrop = 27;
    config.exact_kmer_match_fraction = 0.0;
    config.max_nodes_per_seq_char = 10.0;
    config.queue_size = 20;
    config.num_alternative_paths = 2;
    config.min_path_score = 0;
    config.min_cell_score = 0;

    DBGAligner<UniMEMSeeder<>> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    EXPECT_EQ(64llu, paths[0].get_score()) << paths[0];
    EXPECT_TRUE(paths[0].is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         paths[0],
                         paths.get_query(),
                         paths.get_query_reverse_complement());
}

} // namespace
