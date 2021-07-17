#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "test_aligner_helpers.hpp"
#include "../test_helpers.hpp"

#include "graph/representation/canonical_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;

const std::string test_data_dir = "../tests/data";
const bool PICK_REV_COMP = true;

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
    config.min_seed_length = k;
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
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    EXPECT_TRUE(paths.empty()) << paths.size() << "\t" << paths[0];

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_straight_with_N) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query     = "AGCTNCGAGGCCAA";
    //                           X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    EXPECT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

#if ! _PROTEIN_GRAPH

TYPED_TEST(DBGAlignerTest, align_straight_forward_and_reverse_complement) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;
    reverse_complement(query.begin(), query.end());

    auto graph = build_graph_batch<TypeParam>(k, { reference });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    auto config_fwd_and_rev = config;

    DBGAligner<> aligner(*graph, config_fwd_and_rev);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
    auto ext_paths = get_extend(graph, config_fwd_and_rev, paths, query);

    EXPECT_TRUE(std::equal(paths.begin(), paths.end(),
                           ext_paths.begin(), ext_paths.end()));

    // test copy
    auto paths_copy = paths;
    for (const auto &path : paths_copy) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }

    // test move
    auto paths_move = std::move(paths);
    for (const auto &path : paths_move) {
        EXPECT_TRUE(path.is_valid(*graph, &config));
    }
}

#endif


TYPED_TEST(DBGAlignerTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    //                                 X

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_branch) {
    size_t k = 6;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_branch_with_cycle) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, variation) {
    size_t k = 4;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAACTCGAAA";
    std::string reference_2 = "TTAAGCAAGTCGAAA";
    std::string query =       "TTAAGCAATGGGAAA";
    //                                 XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAACTCTCTGAACTTGT";
    std::string query =     "ACGCAATTCTCTGTATTTGT";
    //                             X      X X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAACTTTTTT";
    std::string reference_2 = "AAAATTGGGGG";
    std::string query =       "AAAATTTTTTT";
    //                             D

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3), -3, -1);
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    ASSERT_EQ(1u, paths.size());
    auto path = paths[0];

    EXPECT_EQ(query.size() - k + 2, path.size());
    EXPECT_EQ(reference_1 + "T", path.get_sequence());
    EXPECT_EQ("4=1D7=", path.get_cigar().to_string());
    EXPECT_EQ(11u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAATTTTTTTT",
                                           "ACAATTTTTGTT",
                                           "ACAAGTTTTTTT",
                                           "ACAAGTTTTGTT"};
    std::string query =                    "ACAACTTTTCTT";
    //                                          X    X

    auto graph = build_graph_batch<TypeParam>(k, references);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2), -3, -1);
    config.num_alternative_paths = 2;
    DBGAligner<> aligner(*graph, config);

    auto paths = aligner.align(query);

    EXPECT_EQ(config.num_alternative_paths, paths.size());
    auto path = paths[0];
    EXPECT_EQ("4=1X4=1X2=", path.get_cigar().to_string())
        << query << "\n" << path.get_sequence();
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAGCGGACCCTTTCCGTTAT";
    std::string query =     "AAAGGGGACCCTTTTCGTTAT";
    //                           X         X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCATTGTT";
    //                            I

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference) + config.gap_opening_penalty, path.get_score());
    EXPECT_EQ("5=1I5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_insert_multi) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAATTGTT";
    //                            II

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference)
        + config.gap_opening_penalty + config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("5=2I5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_insert_long) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAAAAAAAAATTGTT";
    //                            IIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference) + config.gap_opening_penalty
        + score_t(8) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("5=9I5=", path.get_cigar().to_string());
    EXPECT_EQ(10u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_insert_long_offset) {
    size_t k = 5;
    std::string reference = "TTTCCGGTTGTTA";
    std::string query =     "TTTCCGCAAAAAAAAATTGTTA";
    //                             XIIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(reference, "TTTCCGCTTGTTA")
        + config.gap_opening_penalty
        + score_t(8) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "6=1X9I6="
        || path.get_cigar().to_string() == "6=9I1X6=") << path.get_cigar().to_string();
    EXPECT_EQ(12u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_delete) {
    size_t k = 4;
    std::string reference = "TTCGATTGGCCT";
    std::string query =     "TTCGATGGCCT";
    //                             D

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty, path.get_score());

    // TODO: the first should ideally always be true
    EXPECT_TRUE("6=1D5=" == path.get_cigar().to_string()
        || "5=1D6=" == path.get_cigar().to_string());
    // EXPECT_EQ("6=1I5=", path.get_cigar().to_string());

    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    // TODO: enable this when the above is correct
    // check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATACCTTGGCGCTCTC";
    std::string query =     "TTTCTGTATAGGCGCTCTC";
    //                                 DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(3) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("10=4D9=", path.get_cigar().to_string());
    EXPECT_EQ(19u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTTGGCGCTCTC";
    std::string query =     "TTTCGGCGCTCTC";
    //                           DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(3) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("4=4D9=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_loop_deletion) {
    size_t k = 4;
    std::string reference = "AAAATTTTCGAGGCCAA";
    std::string query =     "AAAACGAGGCCAA";
    //                           DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::unit_scoring_matrix(
        1, alphabets::kAlphabetDNA, alphabets::kCharToDNA
    ));
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(13u, path.size());
    EXPECT_EQ("AAAATTTCGAGGCCAA", path.get_sequence());
    EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty
        + score_t(2) * config.gap_extension_penalty, path.get_score());
    EXPECT_EQ("4=3D9=", path.get_cigar().to_string());
    EXPECT_EQ(13u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_drop_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTGGCGCTCTC";
    std::string query =     "TTTCCGGGGCGCTCTC";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -6, -6));
    config.gap_opening_penalty = -10;
    config.gap_extension_penalty = -4;
    config.xdrop = 6;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    // TODO: re-enable this when gap extension is fixed
    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_long_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTTAAGGCGCTCTC";
    std::string query =           "TTTCGGCGCTCTC";
    //                             SSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_repeat_sequence_no_delete_after_insert) {
    size_t k = 27;
    std::string reference = "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
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
    auto path = paths[0];

    EXPECT_EQ(67ull, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
    ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "45=7I8=1X39="
             || path.get_cigar().to_string() == "45=5I1=2I7=1X39="
             || path.get_cigar().to_string() == "44=2I1=5I8=1X39="
             || path.get_cigar().to_string() == "44=3I1=4I8=1X39="
             || path.get_cigar().to_string() == "44=4I1=3I8=1X39=")
        << path.get_cigar().to_string();
    EXPECT_EQ(92u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    auto extends = get_extend(graph, aligner.get_config(), paths, query);
    ASSERT_EQ(1ull, extends.size());
    path = extends[0];

    EXPECT_EQ(67ull, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
        "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
    ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
    EXPECT_TRUE(path.get_cigar().to_string() == "45=7I8=1X39="
             || path.get_cigar().to_string() == "45=5I1=2I7=1X39="
             || path.get_cigar().to_string() == "44=2I1=5I8=1X39="
             || path.get_cigar().to_string() == "44=3I1=4I8=1X39="
             || path.get_cigar().to_string() == "44=4I1=3I8=1X39=")
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
    std::string reference = "GGCCTGTTTG";
    std::string query =     "ACCCTGTTTG";
    //                       SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_clipping2) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAA";
    std::string query =      "TTAGCTTCGAGGCCAA";
    //                        SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_long_clipping) {
    size_t k = 4;
    std::string reference = "TTTTTTTAAAAGCTTCGAGGCCAA";
    std::string query =     "CCCCCCCAAAAGCTTCGAGGCCAA";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_FALSE(paths.empty());

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_end_clipping) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAATTTTTTT";
    std::string query =     "AAAAGCTTCGAGGCCAACCCCCCC";
    //                                        SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_clipping_min_cell_score) {
    size_t k = 7;
    std::string reference = "AAAAGCTTTCGAGGCCAA";
    std::string query =        "ACCTTTCGAGGCCAA";
    //                          SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_low_similarity) {
    size_t k = 27;
    std::string reference = "CTAGAACTTAAAGTATAATAATACTAATAATAAAATAAAATACA";
    std::string query =     "CTAGAACTTAAAGTATAATAATACTAATAAAAGTACAATACA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    // check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

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
    auto path = paths[0];
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
    auto path = paths[0];
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
    std::string match = "TCGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAACGATCAAT"
                        "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
                        "TCGATCAATCGATCAACGATCAATCGATCAATCGATCAACGATCAATCGATCAAT"
                        "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
                        "TCGATCAACGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAAT"
                        "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
                        "CGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAATCGATCAAT"
                        "CGATCAATCGATCAATCGATC";

    EXPECT_TRUE(graph->find(match, 1.0));

    for (double nodes_per_seq_char : { 10.0, std::numeric_limits<double>::max() }) {
        for (size_t xdrop : { 27, 30 }) {
            for (double discovery_fraction : { 0.0, 1.0 }) {
                DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
                config.gap_opening_penalty = -5;
                config.gap_extension_penalty = -2;
                config.xdrop = xdrop;
                config.min_exact_match = discovery_fraction;
                config.max_nodes_per_seq_char = nodes_per_seq_char;
                config.num_alternative_paths = 2;
                config.min_path_score = 0;
                config.min_cell_score = 0;
                config.min_seed_length = k;

                DBGAligner<> aligner(*graph, config);
                auto paths = aligner.align(query);

                if (discovery_fraction == 0.0) {
                    ASSERT_EQ(2ull, paths.size());
                    EXPECT_NE(paths[0], paths[1]);
                    EXPECT_FALSE(paths[0].get_orientation());
                    EXPECT_FALSE(paths[1].get_orientation());
                    EXPECT_GE(paths[0].get_score(), paths[1].get_score());
                } else {
                    EXPECT_EQ(0ull, paths.size());
                }

                paths = aligner.align(match);
                ASSERT_LE(1ull, paths.size());
                EXPECT_EQ(match, paths[0].get_sequence());
                EXPECT_TRUE(paths[0].is_exact_match());
            }
        }
    }
}

TYPED_TEST(DBGAlignerTest, align_low_similarity5) {
    size_t k = 31;
    std::string reference = "GTCGTCAGATCGGAAGAGCGTCGTGTAGGGAAAGGTCTTCGCCTGTGTAGATCTCGGTGGTCG";
    std::string query =     "GTCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTCCTGGTGGTGTAGATC";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);

}

TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
    {
        size_t k = 7;
        std::string reference = "AAAAGCTTTCGAGGCCAA";
        std::string query =        "ACCTTTCGAGGCCAA";
        //                          SS

        auto graph = std::make_shared<DBGSuccinct>(k);
        graph->add_sequence(reference);

        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        config.min_seed_length = 2;
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
        config.min_path_score = std::numeric_limits<score_t>::min() + 100;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths[0];

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
        check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

        check_extend(graph, aligner.get_config(), paths, query);
    }
    {
        size_t k = 15;
        std::string reference = "TTTCGAGGCCAAAGCTTTCGAGGCCAA";
        std::string query =                 "ACCTTTCGAGGCCAA";
        //                                    X

        auto graph = std::make_shared<DBGSuccinct>(k);
        graph->add_sequence(reference);

        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -1));
        config.min_seed_length = 1;
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
        config.min_path_score = std::numeric_limits<score_t>::min() + 100;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths[0];

        EXPECT_EQ(1u, path.size()); // includes dummy k-mers
        EXPECT_EQ(reference.substr(12), path.get_sequence());
        EXPECT_EQ(config.score_sequences(query, reference.substr(12)), path.get_score());
        EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
        EXPECT_EQ(14u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

        check_extend(graph, aligner.get_config(), paths, query);
    }
}

#if ! _PROTEIN_GRAPH
TEST(DBGAlignerTest, align_suffix_seed_snp_canonical) {
    size_t k = 18;
    std::string reference = "AAAAACTTTCGAGGCCAA";
    std::string query =     "GGGGGCTTTCGAGGCCAA";
    //                       SSSSS

    std::string reference_rc = "TTGGCCTCGAAAGTTTTT";
    std::string query_rc     = "TTGGCCTCGAAAGCCCCC";

    for (auto mode : { DeBruijnGraph::PRIMARY, DeBruijnGraph::CANONICAL }) {
        auto dbg_succ = std::make_shared<DBGSuccinct>(k, mode);
        dbg_succ->add_sequence(reference_rc);

        std::shared_ptr<DeBruijnGraph> graph;
        if (mode == DeBruijnGraph::PRIMARY) {
            graph = std::make_shared<CanonicalDBG>(*dbg_succ);
        } else {
            graph = dbg_succ;
        }

        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
        config.min_path_score = std::numeric_limits<score_t>::min() + 100;
        config.min_seed_length = 13;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths[0];

        EXPECT_EQ(1u, path.size()); // includes dummy k-mers
        if (path.get_sequence() == reference.substr(5)) {
            EXPECT_EQ(config.score_sequences(query.substr(5), reference.substr(5)),
                      path.get_score());
            EXPECT_EQ("5S13=", path.get_cigar().to_string());
            EXPECT_EQ(5u, path.get_clipping());
            EXPECT_EQ(0u, path.get_end_clipping());
        } else {
            ASSERT_EQ(reference_rc.substr(0, 13), path.get_sequence());
            EXPECT_EQ(config.score_sequences(query_rc.substr(0, 13),
                                             reference_rc.substr(0, 13)),
                      path.get_score());
            EXPECT_EQ("13=5S", path.get_cigar().to_string());
            EXPECT_EQ(0u, path.get_clipping());
            EXPECT_EQ(5u, path.get_end_clipping());
        }
        EXPECT_EQ(5u, path.get_offset());
        EXPECT_EQ(13u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

        // TODO: sub-k seeds which are sink tips in the underlying graph currently can't be found
        if (mode != DeBruijnGraph::PRIMARY)
            check_extend(graph, aligner.get_config(), paths, query);
    }
}

TYPED_TEST(DBGAlignerTest, align_both_directions) {
    size_t k = 7;
    std::string reference =    "AAAAGCTTTCGAGGCCAA";
    std::string query =        "AAAAGTTTTCGAGGCCAA";
    //                               X

    std::string reference_rc = "TTGGCCTCGAAAGCTTTT";
    std::string query_rc =     "TTGGCCTCGAAAACTTTT";
    //                                      X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::CANONICAL);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(12u, path.size());

    if (path.get_sequence() == reference) {
        EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
        EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
    } else {
        ASSERT_EQ(reference_rc, path.get_sequence());
        EXPECT_EQ(config.score_sequences(query_rc, reference_rc), path.get_score());
        EXPECT_EQ("12=1X5=", path.get_cigar().to_string());
    }

    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(0u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_both_directions2) {
    size_t k = 11;
    std::string reference =    "GTAGTGCTAGCTGTAGTCGTGCTGATGC";
    std::string query =        "GTAGTGCTACCTGTAGTCGTGGTGATGC";
    //                                   X           X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::BASIC);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_EQ(18u, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TYPED_TEST(DBGAlignerTest, align_low_similarity4_rep_primary) {
    size_t k = 6;
    std::vector<std::string> seqs;
    mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
                                          [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
    auto graph = build_graph_batch<TypeParam>(k, std::move(seqs), DeBruijnGraph::PRIMARY);

    std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
                        "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
                        "CGATCGATCGATCGATCGATCGA";

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.gap_opening_penalty = -5;
    config.gap_extension_penalty = -2;
    config.xdrop = 27;
    config.min_exact_match = 0.0;
    config.max_nodes_per_seq_char = 10.0;
    config.num_alternative_paths = 3;

    DBGAligner<> aligner(*graph, config);
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(3u, aligner.align(query).size()) << i;
    }
}
#endif

TYPED_TEST(DBGAlignerTest, align_nodummy) {
    size_t k = 7;
    std::string reference = "AAAAGCTTTCGAGGCCAA";
    std::string query =     "AAAAGTTTTCGAGGCCAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));

    for (bool both_directions : { false, true }) {
#if _PROTEIN_GRAPH
        if (both_directions)
            continue;
#endif
        config.forward_and_reverse_complement = both_directions;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths[0];

        if (both_directions) {
            EXPECT_EQ(12u, path.size());
            EXPECT_EQ(reference, path.get_sequence());
            EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
            EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
            EXPECT_EQ(17u, path.get_num_matches());
        } else {
            EXPECT_EQ(6u, path.size());
            EXPECT_EQ(reference.substr(6), path.get_sequence());
            EXPECT_EQ(config.score_sequences(query.substr(6), reference.substr(6)), path.get_score());
            EXPECT_EQ("6S12=", path.get_cigar().to_string());
            EXPECT_EQ(12u, path.get_num_matches());
            EXPECT_EQ(6u, path.get_clipping());
        }
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

        check_extend(graph, aligner.get_config(), paths, query);
    }
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
    auto path = paths[0];
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_dummy) {
    size_t k = 7;
    std::string reference = "AAAAGCTTTCGAGGCCAA";
    std::string query =     "AAAAGTTTTCGAGGCCAA";
    //                            X

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_seed_length = 5;
    graph->add_sequence(reference);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

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
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_extended_insert_after_match) {
    size_t k = 27;
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
    std::string query =       "CGTGGCCCAGGCCCAGGCCCAGTGGGCGTTGGCCCAGGCGGCCACGGTGGCTGCGCAGGCCCGCCTGGCACAAGCCACGCTG";
    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -3, -3));
    config.min_seed_length = 15;
    graph->add_sequence(reference_1);
    graph->add_sequence(reference_2);

    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths[0];

    EXPECT_TRUE(path.is_valid(*graph, &config));
    EXPECT_EQ(52u, path.get_score())
        << path.get_score() << " " << path.get_cigar().to_string();
    check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

    check_extend(graph, aligner.get_config(), paths, query);
}

TEST(DBGAlignerTest, align_suffix_seed_no_full_seeds) {
    size_t k = 31;
    std::string reference = "CTGCTGCGCCATCGCAACCCACGGTTGCTTTTTGAGTCGCTGCTCACGTTAGCCATCACACTGACGTTAAGCTGGCTTTCGATGCTGTATC";
    std::string query     = "CTTACTGCTGCGCTCTTCGCAAACCCCACGGTTTCTTGTTTTGAGCTCGCCTGCTCACGATACCCATACACACTGACGTTCAAGCTGGCTTTCGATGTTGTATC";

    auto dbg_succ = std::make_shared<DBGSuccinct>(k, DeBruijnGraph::PRIMARY);
    dbg_succ->add_sequence(reference);
    auto graph = std::make_shared<CanonicalDBG>(*dbg_succ);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
    config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
    config.min_path_score = std::numeric_limits<score_t>::min() + 100;
    config.min_seed_length = 13;

    for (size_t max_seed_length : { (size_t)0, k + 100 }) {
        config.max_seed_length = max_seed_length;
        DBGAligner<> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths[0];
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));
    }
}

} // namespace
