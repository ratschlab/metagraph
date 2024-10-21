#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"

#include "graph/representation/canonical_dbg.hpp"
#include "seq_io/sequence_io.hpp"

#include "common/seq_tools/reverse_complement.hpp"
#include "kmer/alphabets.hpp"

#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_match.hpp"
#include "graph/alignment_redone/aligner_config.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align_redone;
using namespace mtg::test;
using namespace mtg::kmer;

const std::string test_data_dir = "../tests/data";
const bool PICK_REV_COMP = true;

template <typename Graph>
class DBGAlignerRedoneTest : public DeBruijnGraphTest<Graph> {};

TYPED_TEST_SUITE(DBGAlignerRedoneTest, FewGraphTypes);

void run_alignment(const DeBruijnGraph &graph,
                   DBGAlignerConfig config,
                   std::string_view query,
                   const std::vector<std::string> &reference,
                   const std::vector<std::string> &cigar_str,
                   size_t end_trim = 0) {
    size_t k = graph.get_k();
    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        config.min_seed_length = k;
        config.max_seed_length = mx;

        Query aln_query(graph, query);
        ExactSeeder seeder(aln_query, config);
        std::vector<Alignment> paths;
        Extender extender(aln_query, config);
        // for (const auto &base_path : seeder.get_inexact_anchors()) {
        for (const auto &base_path : seeder.get_alignments()) {
            extender.extend(base_path, [&](Alignment&& path) {
                paths.emplace_back(std::move(path));
            });
        }

        ASSERT_LE(reference.size(), paths.size()) << mx;
        paths.resize(reference.size());

        for (size_t i = 0; i < reference.size(); ++i) {
            auto path = paths[i];
            if (reference[i].size()) {
                EXPECT_EQ(reference[i].size() - k + 1 + end_trim, path.get_path().size());
                EXPECT_EQ(reference[i], path.get_spelling()) << mx;
            }
            EXPECT_EQ(end_trim, path.get_end_trim()) << mx;

            if (i < cigar_str.size()) {
                Cigar cigar(cigar_str[i]);
                EXPECT_EQ(cigar_str[i], path.get_cigar().to_string()) << mx;
                if (reference[i].size()) {
                    EXPECT_EQ(config.score_cigar(reference[i], query, cigar), path.get_score()) << mx;
                }

                EXPECT_EQ(cigar.get_clipping(), path.get_clipping());
                EXPECT_EQ(cigar.get_end_clipping(), path.get_end_clipping());
            }
        }
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_empty) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, align_sequence_much_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, align_sequence_too_short) {
    size_t k = 4;
    std::string reference = "CATTT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, {}, {});
}

TYPED_TEST(DBGAlignerRedoneTest, align_big_self_loop) {
    size_t k = 3;
    std::string reference = "AAAA";
    std::string query =     "AAAAAAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_single_node) {
    size_t k = 3;
    std::string reference = "CAT";
    std::string query =     "CAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "3=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = reference;

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "14=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight_with_N) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query     = "AGCTNCGAGGCCAA";
    //                           X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "4=1X9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_ending_branch) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAA";
    std::string reference_2 = "AGCTTCGAC";
    //                                 X

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_branch) {
    size_t k = 6;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "17=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_branch_with_cycle) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAATATTTGTT";
    std::string reference_2 = "AGCTTCGACGATTTGTT";
    //                                 XX

    std::string query = reference_2;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "17=" });
}

TYPED_TEST(DBGAlignerRedoneTest, repetitive_sequence_alignment) {
    size_t k = 3;
    std::string reference = "AGGGGGGGGGAAAAGGGGGGG";
    std::string query =     "AGGGGG";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { query }, { "6=" });
}

TYPED_TEST(DBGAlignerRedoneTest, variation) {
    size_t k = 4;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X6=" });
}

TYPED_TEST(DBGAlignerRedoneTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAACTCGAAA";
    std::string reference_2 = "TTAAGCAAGTCGAAA";
    std::string query =       "TTAAGCAATGGGAAA";
    //                                 XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        DBGAlignerConfig config;
        config.min_seed_length = k;
        config.max_seed_length = mx;
        config.gap_opening_penalty = -3;
        config.gap_extension_penalty = -1;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

        Query aln_query(*graph, query);
        ExactSeeder seeder(aln_query, config);
        auto paths = seeder.get_inexact_anchors();

        ASSERT_LE(1ull, paths.size());
        auto path = paths[0];

        EXPECT_EQ(query.size() - k + 1, path.size());
        EXPECT_TRUE(path.get_spelling().compare(reference_1) == 0 ||
                    path.get_spelling().compare(reference_2) == 0)
            << "Path: " << path.get_spelling() << std::endl
            << "Ref1: " << reference_1 << std::endl
            << "Ref2: " << reference_2 << std::endl;
        // TODO: what about other cases?
        EXPECT_EQ("8=3X4=", path.get_cigar().to_string());
        EXPECT_EQ(12u, path.get_cigar().get_num_matches());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_end_trim());
    }
}

TYPED_TEST(DBGAlignerRedoneTest, multiple_variations) {
    size_t k = 4;
    std::string reference = "ACGCAACTCTCTGAACTTGT";
    std::string query =     "ACGCAATTCTCTGTATTTGT";
    //                             X      X X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "6=1X6=1X1=1X4=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_noise_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "AAAACTTTTTT";
    std::string reference_2 = "AAAATTGGGGG";
    std::string query =       "AAAATTTTTTT";
    //                             D

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        DBGAlignerConfig config;
        config.gap_opening_penalty = -3;
        config.gap_extension_penalty = -1;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
        config.min_seed_length = k;
        config.max_seed_length = mx;

        Query aln_query(*graph, query);
        ExactSeeder seeder(aln_query, config);
        auto paths = seeder.get_inexact_anchors();

        ASSERT_LE(1u, paths.size());
        auto path = paths[0];

        EXPECT_EQ(query.size() - k + 2, path.size());
        if (!path.get_orientation()) {
            // the forward orientation of the read was aligned
            EXPECT_EQ("AAAACTTTTTTT", path.get_spelling());
            EXPECT_EQ("4=1D7=", path.get_cigar().to_string());
        } else {
            // the reverse complement of the read was aligned
            EXPECT_EQ("AAAAAAACTTTT", path.get_spelling());
            EXPECT_EQ("7=1D4=", path.get_cigar().to_string());
        }
        EXPECT_EQ(11u, path.get_cigar().get_num_matches());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_end_trim());
    }
}

TYPED_TEST(DBGAlignerRedoneTest, alternative_path_basic) {
    size_t k = 4;
    std::vector<std::string> references = {"ACAATTTTTTTT",
                                           "ACAATTTTTGTT",
                                           "ACAAGTTTTTTT",
                                           "ACAAGTTTTGTT"};
    std::string query =                    "ACAACTTTTCTT";
    //                                          X    X

    auto graph = build_graph_batch<TypeParam>(k, references);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { "" }, { "4=1X4=1X2=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_multiple_misalignment) {
    size_t k = 4;
    std::string reference = "AAAGCGGACCCTTTCCGTTAT";
    std::string query =     "AAAGGGGACCCTTTTCGTTAT";
    //                           X         X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "4=1X9=1X6=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_non_existent) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCATTGTT";
    //                            I

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "5=1I5=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_multi) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAATTGTT";
    //                            II

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "5=2I5=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAAAAAAAAATTGTT";
    //                            IIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference }, { "5=9I5=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long_offset) {
    size_t k = 5;
    std::string reference = "TTTCCGGTTGTTA";
    std::string query =     "TTTCCGCAAAAAAAAATTGTTA";
    //                             XIIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        DBGAlignerConfig config;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
        config.gap_opening_penalty = -1;
        config.gap_extension_penalty = -1;
        config.min_seed_length = k;
        config.max_seed_length = mx;

        Query aln_query(*graph, query);
        ExactSeeder seeder(aln_query, config);
        auto paths = seeder.get_inexact_anchors();

        ASSERT_LE(1ull, paths.size());
        auto path = paths[0];

        EXPECT_EQ(reference.size() - k + 1, path.size());
        EXPECT_EQ(reference, path.get_spelling());
        EXPECT_EQ(config.score_sequences(reference, "TTTCCGCTTGTTA")
            + config.gap_opening_penalty
            + 8 * config.gap_extension_penalty, path.get_score());
        EXPECT_TRUE(path.get_cigar().to_string() == "6=1X9I6="
            || path.get_cigar().to_string() == "6=9I1X6=") << path.get_cigar().to_string();
        EXPECT_EQ(12u, path.get_cigar().get_num_matches());
        EXPECT_EQ(0u, path.get_clipping());
        EXPECT_EQ(0u, path.get_end_clipping());
        EXPECT_EQ(0u, path.get_end_trim());
    }
}

// TYPED_TEST(DBGAlignerTest, align_delete) {
//     size_t k = 4;
//     std::string reference = "TTCGATTGGCCT";
//     std::string query =     "TTCGATGGCCT";
//     //                             D

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.gap_opening_penalty = -3;
//     config.gap_extension_penalty = -3;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(reference.size() - k + 1, path.size());
//     EXPECT_EQ(reference, path.get_sequence());
//     EXPECT_EQ(config.match_score(query) + config.gap_opening_penalty, path.get_score());

//     // TODO: the first should ideally always be true
//     EXPECT_TRUE("6=1D5=" == path.get_cigar().to_string()
//         || "5=1D6=" == path.get_cigar().to_string());
//     // EXPECT_EQ("6=1I5=", path.get_cigar().to_string());

//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(0u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     // TODO: enable this when the above is correct
//     // check_extend(graph, aligner.get_config(), paths, query);
// }

TYPED_TEST(DBGAlignerRedoneTest, align_gap) {
    size_t k = 4;
    std::string reference = "TTTCTGTATACCTTGGCGCTCTC";
    std::string query =     "TTTCTGTATAGGCGCTCTC";
    //                                 DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "10=4D9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTTGGCGCTCTC";
    std::string query =     "TT""TC""GGCGCTCTC";
    //                           II

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { "TTGGCGCTCTC" }, { "2=2I9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_loop_deletion) {
    size_t k = 4;
    std::string reference = "AAAATTTTCGAGGCCAA";
    std::string query =     "AAAACGAGGCCAA";
    //                           DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::unit_scoring_matrix(
        1, alphabets::kAlphabetDNA, alphabets::kCharToDNA
    );
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { "AAAATTTCGAGGCCAA" }, { "4=3D9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_straight_long_xdrop) {
    size_t k = 4;
    std::string reference_1 = "AGCTTCGAGGCCAAGCCTGACTGATCGATGCATGCTAGCTAGTCAGTCAGCGTGAGCTAGCAT";
    std::string reference_2 = "AGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    std::string query = reference_1;

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 30;
    config.rel_score_cutoff = 0.8;
    run_alignment(*graph, config, query, { query }, { "63=" });
}

// TYPED_TEST(DBGAlignerTest, align_drop_seed) {
//     size_t k = 4;
//     std::string reference = "TTTCCCTGGCGCTCTC";
//     std::string query =     "TTTCCGGGGCGCTCTC";
//     //                       SSSSSSS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -6, -6);
//     config.gap_opening_penalty = -10;
//     config.gap_extension_penalty = -4;
//     config.xdrop = 6;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(6, path.size());
//     EXPECT_EQ(reference.substr(7), path.get_sequence());
//     EXPECT_EQ(config.match_score(reference.substr(7)), path.get_score());
//     EXPECT_EQ("7S9=", path.get_cigar().to_string());
//     EXPECT_EQ(9u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(7u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     // TODO: re-enable this when gap extension is fixed
//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_long_gap_after_seed) {
//     size_t k = 4;
//     std::string reference = "TTTCCCTTAAGGCGCTCTC";
//     std::string query =           "TTTCGGCGCTCTC";
//     //                             SSSS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -1;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(6, path.size());
//     EXPECT_EQ(reference.substr(10), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(4)), path.get_score());
//     EXPECT_EQ("4S9=", path.get_cigar().to_string());
//     EXPECT_EQ(9u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(4u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_repeat_sequence_no_delete_after_insert) {
//     size_t k = 27;
//     std::string reference = "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
//     std::string query =     "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
//     // the alignment may be misrepresented as
//     // "TTTGTGGCTAGAGCTCGAGATCGCGCG"                    "GCCACAATT" "GACAAATG" "A" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
//     // "TTTGTGGCTAGAGCTCGAGATCGCGCG" "GCCACAATTGACAAAT"             "GACAAATG" "T" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.gap_opening_penalty = -3;
//     config.gap_extension_penalty = -3;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(67ull, path.size());
//     EXPECT_EQ(reference, path.get_sequence());
//     EXPECT_EQ(config.score_sequences(
//         "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
//         "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGAGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
//     ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
//     EXPECT_TRUE(path.get_cigar().to_string() == "45=7I8=1X39="
//              || path.get_cigar().to_string() == "45=5I1=2I7=1X39="
//              || path.get_cigar().to_string() == "44=2I1=5I8=1X39="
//              || path.get_cigar().to_string() == "44=3I1=4I8=1X39="
//              || path.get_cigar().to_string() == "44=4I1=3I8=1X39=")
//         << path.get_cigar().to_string();
//     EXPECT_EQ(92u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(0u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     auto extends = get_extend(graph, aligner.get_config(), paths, query);
//     ASSERT_EQ(1ull, extends.size());
//     path = extends[0];

//     EXPECT_EQ(67ull, path.size());
//     EXPECT_EQ(reference, path.get_sequence());
//     EXPECT_EQ(config.score_sequences(
//         "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "T" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC",
//         "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA" "GATCTAAT" "G" "AAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC"
//     ) + config.gap_opening_penalty + score_t(6) * config.gap_extension_penalty, path.get_score());
//     EXPECT_TRUE(path.get_cigar().to_string() == "45=7I8=1X39="
//              || path.get_cigar().to_string() == "45=5I1=2I7=1X39="
//              || path.get_cigar().to_string() == "44=2I1=5I8=1X39="
//              || path.get_cigar().to_string() == "44=3I1=4I8=1X39="
//              || path.get_cigar().to_string() == "44=4I1=3I8=1X39=")
//         << path.get_cigar().to_string();
//     EXPECT_EQ(92u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(0u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
// }

TYPED_TEST(DBGAlignerRedoneTest, align_clipping1) {
    size_t k = 4;
    std::string reference = "GGCCTGTTTG";
    std::string query =     "ACCCTGTTTG";
    //                       SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(2) }, { "2S8=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_clipping2) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAA";
    std::string query =      "TTAGCTTCGAGGCCAA";
    //                        SS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(3) }, { "2S14=" });
}

// TYPED_TEST(DBGAlignerTest, align_long_clipping) {
//     size_t k = 4;
//     std::string reference = "TTTTTTTAAAAGCTTCGAGGCCAA";
//     std::string query =     "CCCCCCCAAAAGCTTCGAGGCCAA";
//     //                       SSSSSSS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_FALSE(paths.empty());

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(14u, path.size());
//     EXPECT_EQ(reference.substr(7), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(7)), path.get_score());
//     EXPECT_EQ("7S17=", path.get_cigar().to_string());
//     EXPECT_EQ(17u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(7u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

TYPED_TEST(DBGAlignerRedoneTest, align_end_clipping) {
    size_t k = 4;
    std::string reference = "AAAAGCTTCGAGGCCAATTTTTTT";
    std::string query =     "AAAAGCTTCGAGGCCAACCCCCCC";
    //                                        SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(0, 17) }, { "17=7S" });
}

// TYPED_TEST(DBGAlignerTest, align_clipping_min_cell_score) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =        "ACCTTTCGAGGCCAA";
//     //                          SS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(7u, path.size());
//     EXPECT_EQ(reference.substr(5), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
//     EXPECT_EQ("2S13=", path.get_cigar().to_string());
//     EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(2u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_low_similarity) {
//     size_t k = 27;
//     std::string reference = "CTAGAACTTAAAGTATAATAATACTAATAATAAAATAAAATACA";
//     std::string query =     "CTAGAACTTAAAGTATAATAATACTAATAAAAGTACAATACA";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     // EXPECT_EQ(7u, path.size());
//     // EXPECT_EQ(reference.substr(5), path.get_sequence());
//     // EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
//     // EXPECT_EQ("2S13=", path.get_cigar().to_string());
//     // EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//     // EXPECT_FALSE(is_exact_match(path));
//     // EXPECT_EQ(2u, path.get_clipping());
//     // EXPECT_EQ(0u, path.get_end_clipping());
//     // EXPECT_EQ(0u, path.get_offset());
//     // EXPECT_TRUE(path.is_valid(*graph, &config));
//     // check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     // check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_low_similarity2) {
//     size_t k = 27;
//     std::string reference = "GCCACAATTGACAAATGAGATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
//     std::string query =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);

//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];
// }

// TYPED_TEST(DBGAlignerTest, align_low_similarity3) {
//     size_t k = 27;
//     std::string reference = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTGCTGGGATTATAGGTGTGAACCACCACACCTGGCTAATTTTTTTTGTGTGTGTGTGTGTTTTTTC";
//     std::string query =     "AAAAAAAAAAAAAAAAAAAAAAAAAAACGCCAAAAAGGGGGAATAGGGGGGGGGGAACCCCAACACCGGTATGTTTTTTTGTGTGTGGGGGATTTTTTTC";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     for (bool seed_complexity_filter : { false, true }) {
//         DBGAlignerConfig config;
//         config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//         config.seed_complexity_filter = seed_complexity_filter;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
// #if ! _PROTEIN_GRAPH
//         EXPECT_EQ(seed_complexity_filter, paths.empty());
// #else
//         EXPECT_FALSE(paths.empty());
// #endif
//     }
// }

// TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4) {
//     size_t k = 6;
//     std::vector<std::string> seqs;
//     mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
//                                           [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
//     auto graph = build_graph_batch<TypeParam>(k, std::move(seqs));

//     std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGA";
//     std::string match = "TCGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAACGATCAAT"
//                         "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
//                         "TCGATCAATCGATCAACGATCAATCGATCAATCGATCAACGATCAATCGATCAAT"
//                         "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
//                         "TCGATCAACGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAAT"
//                         "CGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAATCGATCAA"
//                         "CGATCAATCGATCAATCGATCAACGATCAATCGATCAATCGATCAATCGATCAAT"
//                         "CGATCAATCGATCAATCGATC";

//     EXPECT_TRUE(graph->find(match, 1.0));

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     run_alignment(*graph, config, query, { match }, {});
// }

// TYPED_TEST(DBGAlignerTest, align_low_similarity5) {
//     size_t k = 31;
//     std::string reference = "GTCGTCAGATCGGAAGAGCGTCGTGTAGGGAAAGGTCTTCGCCTGTGTAGATCTCGGTGGTCG";
//     std::string query =     "GTCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTTCCTGGTGGTGTAGATC";

//     auto graph = std::make_shared<DBGSuccinct>(k);
//     graph->add_sequence(reference);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);

// }

// TEST(DBGAlignerTest, align_suffix_seed_snp_min_seed_length) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =        "ACCTTTCGAGGCCAA";
//     //                          SS

//     auto graph = std::make_shared<DBGSuccinct>(k);
//     graph->add_sequence(reference);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.min_seed_length = 2;
//     config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
//     config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(7u, path.size());
//     EXPECT_EQ(reference.substr(5), path.get_sequence());
//     EXPECT_EQ(config.match_score(query.substr(2)), path.get_score());
//     EXPECT_EQ("2S13=", path.get_cigar().to_string());
//     EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(2u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// #if ! _PROTEIN_GRAPH
// TEST(DBGAlignerTest, align_suffix_seed_snp_canonical) {
//     size_t k = 18;
//     std::string reference = "AAAAACTTTCGAGGCCAA";
//     std::string query =     "GGGGGCTTTCGAGGCCAA";
//     //                       SSSSS

//     std::string reference_rc = "TTGGCCTCGAAAGTTTTT";
//     std::string query_rc     = "TTGGCCTCGAAAGCCCCC";

//     for (auto mode : { DeBruijnGraph::PRIMARY, DeBruijnGraph::CANONICAL }) {
//         auto dbg_succ = std::make_shared<DBGSuccinct>(k, mode);
//         dbg_succ->add_sequence(reference_rc);

//         std::shared_ptr<DeBruijnGraph> graph = dbg_succ;
//         if (mode == DeBruijnGraph::PRIMARY)
//             graph = std::make_shared<CanonicalDBG>(graph);

//         DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//         config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
//         config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//         config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//         config.min_seed_length = 13;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
//         ASSERT_EQ(1ull, paths.size());
//         auto path = paths[0];

//         EXPECT_EQ(1u, path.size()); // includes dummy k-mers
//         if (path.get_sequence() == reference.substr(5)) {
//             EXPECT_EQ(config.score_sequences(query.substr(5), reference.substr(5)),
//                       path.get_score());
//             EXPECT_EQ("5S13=", path.get_cigar().to_string());
//             EXPECT_EQ(5u, path.get_clipping());
//             EXPECT_EQ(0u, path.get_end_clipping());
//         } else {
//             ASSERT_EQ(reference_rc.substr(0, 13), path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query_rc.substr(0, 13),
//                                              reference_rc.substr(0, 13)),
//                       path.get_score());
//             EXPECT_EQ("13=5S", path.get_cigar().to_string());
//             EXPECT_EQ(0u, path.get_clipping());
//             EXPECT_EQ(5u, path.get_end_clipping());
//         }
//         EXPECT_EQ(5u, path.get_offset());
//         EXPECT_EQ(13u, path.get_cigar().get_num_matches());
//         EXPECT_FALSE(is_exact_match(path));
//         EXPECT_TRUE(path.is_valid(*graph, &config));
//         check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//         // TODO: sub-k seeds which are sink tips in the underlying graph currently can't be found
//         if (mode != DeBruijnGraph::PRIMARY)
//             check_extend(graph, aligner.get_config(), paths, query);
//     }
// }

// TYPED_TEST(DBGAlignerTest, align_both_directions) {
//     size_t k = 7;
//     std::string reference =    "AAAAGCTTTCGAGGCCAA";
//     std::string query =        "AAAAGTTTTCGAGGCCAA";
//     //                               X

//     std::string reference_rc = "TTGGCCTCGAAAGCTTTT";
//     std::string query_rc =     "TTGGCCTCGAAAACTTTT";
//     //                                      X

//     auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::CANONICAL);
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(12u, path.size());

//     if (path.get_sequence() == reference) {
//         EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
//         EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
//     } else {
//         ASSERT_EQ(reference_rc, path.get_sequence());
//         EXPECT_EQ(config.score_sequences(query_rc, reference_rc), path.get_score());
//         EXPECT_EQ("12=1X5=", path.get_cigar().to_string());
//     }

//     EXPECT_EQ(17u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(0u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_both_directions2) {
//     size_t k = 11;
//     std::string reference =    "GTAGTGCTAGCTGTAGTCGTGCTGATGC";
//     std::string query =        "GTAGTGCTACCTGTAGTCGTGGTGATGC";
//     //                                   X           X

//     auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::BASIC);
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(18u, path.size());
//     EXPECT_EQ(reference, path.get_sequence());
//     EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TYPED_TEST(DBGAlignerTest, align_low_similarity4_rep_primary) {
//     size_t k = 6;
//     std::vector<std::string> seqs;
//     mtg::seq_io::read_fasta_file_critical(test_data_dir + "/transcripts_100.fa",
//                                           [&](auto *seq) { seqs.emplace_back(seq->seq.s); });
//     auto graph = build_graph_batch<TypeParam>(k, std::move(seqs), DeBruijnGraph::PRIMARY);

//     std::string query = "TCGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGATCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "TCGATCGACGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
//                         "CGATCGATCGATCGATCGATCGACGATCGATCGATCGATCGATCGATCGATCGAT"
//                         "CGATCGATCGATCGATCGATCGA";

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     config.xdrop = 27;
//     config.min_exact_match = 0.0;
//     config.max_nodes_per_seq_char = 10.0;
//     config.num_alternative_paths = 3;

//     DBGAligner<> aligner(*graph, config);
//     for (size_t i = 0; i < 3; ++i) {
//         EXPECT_EQ(3u, aligner.align(query).size()) << i;
//     }
// }
// #endif

// TYPED_TEST(DBGAlignerTest, align_nodummy) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =     "AAAAGTTTTCGAGGCCAA";
//     //                            X

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

//     for (bool both_directions : { false, true }) {
// #if _PROTEIN_GRAPH
//         if (both_directions)
//             continue;
// #endif
//         config.forward_and_reverse_complement = both_directions;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
//         ASSERT_EQ(1ull, paths.size());
//         auto path = paths[0];

//         if (both_directions) {
//             EXPECT_EQ(12u, path.size());
//             EXPECT_EQ(reference, path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
//             EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
//             EXPECT_EQ(17u, path.get_cigar().get_num_matches());
//         } else {
//             EXPECT_EQ(6u, path.size());
//             EXPECT_EQ(reference.substr(6), path.get_sequence());
//             EXPECT_EQ(config.score_sequences(query.substr(6), reference.substr(6)), path.get_score());
//             EXPECT_EQ("6S12=", path.get_cigar().to_string());
//             EXPECT_EQ(12u, path.get_cigar().get_num_matches());
//             EXPECT_EQ(6u, path.get_clipping());
//         }
//         EXPECT_FALSE(is_exact_match(path));
//         EXPECT_EQ(0u, path.get_end_clipping());
//         EXPECT_EQ(0u, path.get_offset());
//         EXPECT_TRUE(path.is_valid(*graph, &config));
//         check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//         check_extend(graph, aligner.get_config(), paths, query);
//     }
// }

// TYPED_TEST(DBGAlignerTest, align_seed_to_end) {
//     size_t k = 5;
//     std::string reference = "ATCCCTTTTAAAA";
//     std::string query =     "ATCCCGGGGGGGGGGGGGGGGGTTTTAAAA";

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

TYPED_TEST(DBGAlignerRedoneTest, align_bfs_vs_dfs_xdrop) {
    size_t k = 31;
    std::string reference_1 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGGATCAGAGAAGAAAGGAATGGCAAATACTCCATAGAAGAGTACCTTCAGGACAGGGACAGATACTATGAGGAGGTGGCCAT";
    std::string reference_2 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGAATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";
    std::string query =       "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGATGAGTACGATGAGTACAAGAGCATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 27;
    run_alignment(*graph, config, query, { "" }, { "48=1X20=1X80=" });
}

// TEST(DBGAlignerTest, align_dummy) {
//     size_t k = 7;
//     std::string reference = "AAAAGCTTTCGAGGCCAA";
//     std::string query =     "AAAAGTTTTCGAGGCCAA";
//     //                            X

//     auto graph = std::make_shared<DBGSuccinct>(k);
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.min_seed_length = 5;
//     graph->add_sequence(reference);

//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_EQ(12u, path.size());
//     EXPECT_EQ(reference, path.get_sequence());
//     EXPECT_EQ(config.score_sequences(query, reference), path.get_score());
//     EXPECT_EQ("5=1X12=", path.get_cigar().to_string());
//     EXPECT_EQ(17u, path.get_cigar().get_num_matches());
//     EXPECT_FALSE(is_exact_match(path));
//     EXPECT_EQ(0u, path.get_clipping());
//     EXPECT_EQ(0u, path.get_end_clipping());
//     EXPECT_EQ(0u, path.get_offset());
//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// TEST(DBGAlignerTest, align_extended_insert_after_match) {
//     size_t k = 27;
//     std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
//     std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC";
//     std::string query =       "CGTGGCCCAGGCCCAGGCCCAGTGGGCGTTGGCCCAGGCGGCCACGGTGGCTGCGCAGGCCCGCCTGGCACAAGCCACGCTG";
//     auto graph = std::make_shared<DBGSuccinct>(k);
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.min_seed_length = 15;
//     graph->add_sequence(reference_1);
//     graph->add_sequence(reference_2);

//     DBGAligner<> aligner(*graph, config);
//     auto paths = aligner.align(query);
//     ASSERT_EQ(1ull, paths.size());
//     auto path = paths[0];

//     EXPECT_TRUE(path.is_valid(*graph, &config));
//     EXPECT_EQ(52u, path.get_score())
//         << path.get_score() << " " << path.get_cigar().to_string();
//     check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));

//     check_extend(graph, aligner.get_config(), paths, query);
// }

// #if ! _PROTEIN_GRAPH
// TEST(DBGAlignerTest, align_suffix_seed_no_full_seeds) {
//     size_t k = 31;
//     std::string reference = "CTGCTGCGCCATCGCAACCCACGGTTGCTTTTTGAGTCGCTGCTCACGTTAGCCATCACACTGACGTTAAGCTGGCTTTCGATGCTGTATC";
//     std::string query     = "CTTACTGCTGCGCTCTTCGCAAACCCCACGGTTTCTTGTTTTGAGCTCGCCTGCTCACGATACCCATACACACTGACGTTCAAGCTGGCTTTCGATGTTGTATC";

//     auto dbg_succ = std::make_shared<DBGSuccinct>(k, DeBruijnGraph::PRIMARY);
//     dbg_succ->add_sequence(reference);
//     auto graph = std::make_shared<CanonicalDBG>(dbg_succ);

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
//     config.max_num_seeds_per_locus = std::numeric_limits<size_t>::max();
//     config.min_cell_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_path_score = std::numeric_limits<score_t>::min() + 100;
//     config.min_seed_length = 13;

//     for (size_t max_seed_length : { (size_t)0, k + 100 }) {
//         config.max_seed_length = max_seed_length;
//         DBGAligner<> aligner(*graph, config);
//         auto paths = aligner.align(query);
//         ASSERT_EQ(1ull, paths.size());
//         auto path = paths[0];
//         EXPECT_TRUE(path.is_valid(*graph, &config));
//         check_json_dump_load(*graph, path, paths.get_query(), paths.get_query(PICK_REV_COMP));
//     }
// }
// #endif

} // namespace
