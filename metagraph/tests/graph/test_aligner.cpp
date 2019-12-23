#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "test_aligner_helpers.hpp"
#include "../test_helpers.hpp"

#include "graph/aligner/dbg_aligner.hpp"
#include "graph/aligner/aligner_methods.hpp"

#include "annotation/column_compressed/annotate_column_compressed.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"


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


// TODO: REPLACE THIS
#if _PROTEIN_GRAPH
    const auto *alphabet = alphabets::kAlphabetProtein;
    const auto *alphabet_encoding = alphabets::kCharToProtein;
#elif _DNA_CASE_SENSITIVE_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNA;
    const auto *alphabet_encoding = alphabets::kCharToDNA;
#elif _DNA5_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNA;
    const auto *alphabet_encoding = alphabets::kCharToDNA;
#elif _DNA_GRAPH
    const auto *alphabet = alphabets::kAlphabetDNA;
    const auto *alphabet_encoding = alphabets::kCharToDNA;
#else
    static_assert(false,
        "Define an alphabet: either "
        "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
    );
#endif

void check_extend(std::shared_ptr<const DeBruijnGraph> graph,
                  const DBGAlignerConfig &config,
                  const DBGAligner<>::DBGQueryAlignment &paths,
                  const std::string &query) {
    assert(graph.get());
    EXPECT_EQ(query, paths.get_query());

    auto unimem_paths = DBGAligner<UniMEMSeeder<>>(*graph, config).align(query);
    ASSERT_EQ(paths.size(), unimem_paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        EXPECT_EQ(paths[i], unimem_paths[i]);
    }
}


template <typename Graph>
class DBGAlignerTest : public DeBruijnGraphTest<Graph> {
    void SetUp() { Cigar::initialize_opt_table(alphabet, alphabet_encoding); }
};

TYPED_TEST_CASE(DBGAlignerTest, GraphTypes);

TYPED_TEST(DBGAlignerTest, bad_min_cell_score) {
    auto graph = build_graph_batch<TypeParam>(3, {});
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.min_cell_score = std::numeric_limits<score_t>::min();
    config.min_path_score = std::numeric_limits<score_t>::min();
    ASSERT_THROW(DBGAligner<>(*graph, config), std::runtime_error);
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(1ull, path.size());
    EXPECT_EQ("CAT", path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(query, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin()),
              path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin()),
              path.get_score());
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
    std::vector<double> weights = {
        std::exp(paths[0].get_score() - config.match_score(paths[0].get_sequence().begin(),
                                                           paths[0].get_sequence().end())),
        std::exp(paths[1].get_score() - config.match_score(paths[1].get_sequence().begin(),
                                                           paths[1].get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));
    EXPECT_NE(paths.front(), paths.back());
    auto path = paths.front();

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference_1, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference_1.begin()),
              path.get_score());
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

    std::vector<double> weights(paths.size());
    std::transform(
        paths.begin(), paths.end(),
        weights.begin(),
        [&](const auto &path) {
            return std::exp(path.get_score()
                - config.match_score(path.get_sequence().begin(),
                                     path.get_sequence().end()));
        }
    );
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(query.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin()),
              path.get_score());
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
        std::vector<double> weights = {
            std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                           path.get_sequence().end()))
        };
        EXPECT_EQ(weights, paths.get_alignment_weights(config));

        EXPECT_EQ(query.size() - k + 1, path.size());
        EXPECT_EQ(reference, path.get_sequence());
        EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin()),
                  path.get_score());
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

TYPED_TEST(DBGAlignerTest, align_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCC"     "TTGTT";
    std::string query =     "TTTCC" "A" "TTGTT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(reference.begin(), reference.end())
                + config.gap_opening_penalty,
              path.get_score());
    EXPECT_EQ("5=1I5=", path.get_cigar().to_string());
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

TYPED_TEST(DBGAlignerTest, align_delete) {
    size_t k = 4;
    std::string reference = "TTCGAT" "TGGCCT";
    std::string query =     "TTCGAT"  "GGCCT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end()) + config.gap_opening_penalty,
              path.get_score());

    // TODO: the first should ideally always be true
    EXPECT_TRUE("6=1D5=" == path.get_cigar().to_string()
        || "5=1D6=" == path.get_cigar().to_string());
    // EXPECT_EQ("6=1D5=", path.get_cigar().to_string());

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
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(reference.size() - k + 1, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.end())
                + config.gap_opening_penalty
                + score_t(3) * config.gap_extension_penalty,
              path.get_score());
    EXPECT_EQ("10=4D9=", path.get_cigar().to_string());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(5ull, path.size());
    EXPECT_EQ(reference.substr(2), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(11u, path.size());
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(14u, path.size());
    EXPECT_EQ(reference.substr(7), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 7, query.end()), path.get_score());
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
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(14u, path.size());
    EXPECT_EQ(reference.substr(0, 17), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin(), query.begin() + 17), path.get_score());
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
    config.min_cell_score = std::numeric_limits<score_t>::min() + 3;
    config.min_path_score = std::numeric_limits<score_t>::min() + 3;
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);

    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(7u, path.size());
    EXPECT_EQ(reference.substr(5), path.get_sequence());
    EXPECT_EQ(config.match_score(query.begin() + 2, query.end()), path.get_score());
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

TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    {
        DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
        config.min_seed_length = 2;
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min() + 3;
        config.min_path_score = std::numeric_limits<score_t>::min() + 3;
        DBGAligner<SuffixSeeder<>> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();
        std::vector<double> weights = {
            std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                           path.get_sequence().end()))
        };
        EXPECT_EQ(weights, paths.get_alignment_weights(config));

        EXPECT_EQ(7u, path.size());
        EXPECT_EQ(reference.substr(5), path.get_sequence());
        EXPECT_EQ(config.match_score(query.begin() + 2, query.end()), path.get_score());
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
        config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
        config.min_cell_score = std::numeric_limits<score_t>::min() + 3;
        config.min_path_score = std::numeric_limits<score_t>::min() + 3;
        DBGAligner<SuffixSeeder<>> aligner(*graph, config);
        auto paths = aligner.align(query);
        ASSERT_EQ(1ull, paths.size());
        auto path = paths.front();
        std::vector<double> weights = {
            std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                           path.get_sequence().end()))
        };
        EXPECT_EQ(weights, paths.get_alignment_weights(config));

        EXPECT_EQ(15u, path.size()); // includes dummy k-mers
        EXPECT_EQ(reference.substr(3), path.get_sequence());
        EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin() + 3),
                  path.get_score());
        EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
        EXPECT_EQ(14u, path.get_num_matches());
        EXPECT_FALSE(path.is_exact_match());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(6u, path.get_offset());
        EXPECT_TRUE(path.is_valid(*graph, &config));
        check_json_dump_load(*graph,
                             path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());
    }
}

TEST(DBGAlignerTest, align_suffix_seed_snp) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "CTTTCGAGGCCAA";
    std::string query =        "AC" "CTTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    graph->add_sequence(reference);

    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
    config.min_cell_score = std::numeric_limits<score_t>::min() + 3;
    config.min_path_score = std::numeric_limits<score_t>::min() + 3;
    DBGAligner<SuffixSeeder<>> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(15u, path.size()); // includes dummy k-mers
    EXPECT_EQ(reference.substr(3), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin() + 3),
              path.get_score());
    EXPECT_EQ("1=1X13=", path.get_cigar().to_string());
    EXPECT_EQ(14u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(6u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());
}

TYPED_TEST(DBGAlignerTest, align_nodummy) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "C" "TTTCGAGGCCAA";
    std::string query =     "AAAAG" "T" "TTTCGAGGCCAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    DBGAligner<> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(6u, path.size());
    EXPECT_EQ(reference.substr(6), path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin() + 6, query.end(),
                                     reference.begin() + 6),
              path.get_score());
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

TEST(DBGAlignerTest, align_dummy) {
    Cigar::initialize_opt_table(alphabet, alphabet_encoding);
    size_t k = 7;
    std::string reference = "AAAAG" "C" "TTTCGAGGCCAA";
    std::string query =     "AAAAG" "T" "TTTCGAGGCCAA";

    auto graph = std::make_shared<DBGSuccinct>(k);
    DBGAlignerConfig config(DBGAlignerConfig::dna_scoring_matrix(2, -1, -2));
    graph->add_sequence(reference);

    DBGAligner<SuffixSeeder<>> aligner(*graph, config);
    auto paths = aligner.align(query);
    ASSERT_EQ(1ull, paths.size());
    auto path = paths.front();
    std::vector<double> weights = {
        std::exp(path.get_score() - config.match_score(path.get_sequence().begin(),
                                                       path.get_sequence().end()))
    };
    EXPECT_EQ(weights, paths.get_alignment_weights(config));

    EXPECT_EQ(14u, path.size());
    EXPECT_EQ(reference, path.get_sequence());
    EXPECT_EQ(config.score_sequences(query.begin(), query.end(), reference.begin()),
              path.get_score());
    EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
    EXPECT_EQ(17u, path.get_num_matches());
    EXPECT_FALSE(path.is_exact_match());
    EXPECT_EQ(0u, path.get_clipping());
    EXPECT_EQ(0u, path.get_end_clipping());
    EXPECT_EQ(2u, path.get_offset());
    EXPECT_TRUE(path.is_valid(*graph, &config));
    check_json_dump_load(*graph,
                         path,
                         paths.get_query(),
                         paths.get_query_reverse_complement());

    // TODO: make uni-mem seeder work with unmasked DBGSuccinct
    // check_extend(graph, aligner.get_config(), paths, query);
}
