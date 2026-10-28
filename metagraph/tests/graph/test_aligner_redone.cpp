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

template <typename Graph>
class DBGAlignerRedoneTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGHashFast, DBGSuccinctCached, DBGSSHash> AlignGraphTypes;
TYPED_TEST_SUITE(DBGAlignerRedoneTest, AlignGraphTypes);

void run_alignment(const DeBruijnGraph &graph,
                   DBGAlignerConfig config,
                   std::string_view query,
                   const std::vector<std::string> &reference,
                   const std::vector<std::string> &cigar_str,
                   size_t end_trim = 0,
                   bool needs_extension = false,
                   bool needs_extension_long_seed = false) {
    size_t k = graph.get_k();
    if (config.min_seed_length == 0)
        config.min_seed_length = k;

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        config.max_seed_length = std::max(mx, config.min_seed_length);
        bool check_chaining = mx == std::numeric_limits<size_t>::max()
            ? !needs_extension_long_seed
            : !needs_extension;

        Query aln_query(graph, query);
        ExactSeeder seeder(aln_query, config);
        std::vector<Alignment> paths;
        Extender extender(aln_query, config);
        std::vector<Alignment> paths_no_extend = seeder.get_inexact_anchors();
        std::vector<Alignment> paths_no_align = seeder.get_inexact_anchors(false);
        for (const auto &base_path : paths_no_extend) {
            extender.extend(base_path, [&](Alignment&& path) {
                paths.emplace_back(std::move(path));
            });
        }

        auto aln_sort = [&config](const auto &a, const auto &b) {
            return std::make_pair(score_match(a, config), b.get_orientation())
                 > std::make_pair(score_match(b, config), a.get_orientation());
        };

        std::sort(paths_no_extend.begin(), paths_no_extend.end(), aln_sort);
        std::sort(paths_no_align.begin(), paths_no_align.end(), aln_sort);
        std::sort(paths.begin(), paths.end(), aln_sort);

        ASSERT_LE(reference.size(), paths.size()) << mx;
        paths.resize(reference.size());

        for (size_t i = 0; i < reference.size(); ++i) {
            auto check_ref = [&](const Alignment &path, const std::string &reference, const std::string &type) {
                if (reference.size()) {
                    EXPECT_EQ(reference.size() - k + 1 + end_trim, path.get_path().size()) << mx << "\t" << type;
                    EXPECT_EQ(reference, path.get_spelling()) << mx << "\t" << type;
                }
            };

            auto check_aln = [&](const Alignment &path, const Cigar &cigar, const std::string &reference, const std::string &type) {
                ASSERT_EQ(cigar.to_string(), path.get_cigar().to_string()) << mx << "\t" << type;
                EXPECT_EQ(end_trim, path.get_end_trim()) << mx << "\t" << type;
                if (reference.size()) {
                    ASSERT_TRUE(cigar.is_valid(reference, path.get_query()));
                    EXPECT_EQ(config.score_cigar(reference, path.get_query(), cigar),
                              score_match(path, config)) << mx << "\t" << type;
                }

                EXPECT_EQ(cigar.get_clipping(), path.get_clipping()) << mx << "\t" << type;
                EXPECT_EQ(cigar.get_end_clipping(), path.get_end_clipping()) << mx << "\t" << type;
            };

            std::string cur_reference = reference[i];
            if (paths[i].get_orientation() && !paths[i].get_end_trim())
                ::reverse_complement(cur_reference.begin(), cur_reference.end());

            check_ref(paths[i], cur_reference, "extend");

            if (i < cigar_str.size() && cigar_str[i].size()) {
                std::string cur_cigar = cigar_str[i];
                Cigar cigar(cur_cigar);
                if (paths[i].get_orientation() && !paths[i].get_end_trim()) {
                    std::reverse(cigar.data().begin(), cigar.data().end());
                    cur_cigar = cigar.to_string();
                }
                check_aln(paths[i], cigar, cur_reference, "extend");

                if (check_chaining) {
                    // this alignment should work with chaining alone
                    ASSERT_LT(i, paths_no_extend.size());
                    check_ref(paths_no_extend[i], cur_reference, "chain");
                    check_aln(paths_no_extend[i], cigar, cur_reference, "chain");

                    ASSERT_LT(i, paths_no_align.size());
                    EXPECT_EQ(cigar.to_string(), paths_no_align[i].get_cigar().to_string()) << mx << "\tnoalign";
                }
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
    size_t k = 5;
    std::string reference = "TAGCTTCGAGGCCAA";
    std::string query     = "TAGCTNCGAGGCCAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

    run_alignment(*graph, config, query, { reference }, { "5=1X9=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X9=" });
    }
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
    size_t k = 5;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";
    //                            X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X6=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X6=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, variation_in_branching_point) {
    size_t k = 4;
    std::string reference_1 = "TTAAGCAACTCGAAA";
    std::string reference_2 = "TTAAGCAAGTCGAAA";
    std::string query =       "TTAAGCAATGGGAAA";
    //                                 XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -2;

    run_alignment(*graph, config, query, { "" }, { "8=3X4=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { "" }, { "8=3X4=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, multiple_variations) {
    size_t k = 6;
    std::string reference = "TACGCAACTCTCTGAACTTGT";
    std::string query =     "TACGCAATTCTCTGTATTTGT";
    //                              X      X X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "7=1X6=1X1=1X4=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 4;
        run_alignment(*graph, config, query, { reference }, { "7=1X6=1X1=1X4=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_noise_in_branching_point) {
    size_t k = 5;
    std::string reference_1 = "AAAAACTTTTTT";
    std::string reference_2 = "AAAAATTGGGGG";
    std::string query =       "AAAAATTTTTTT";
    //                              D

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);

    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { "AAAAACTTTTTTT" }, { "5=1D7=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { "AAAAACTTTTTTT" }, { "5=1D7=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, alternative_path_basic) {
    size_t k = 5;
    std::vector<std::string> references = {
                        "AGACAATTTTTTAAA",
                        "AGACAATTTTTGAAA",
                        "AGACAAGTTTTTAAA",
                        "AGACAAGTTTTGAAA" };
    std::string query = "AGACAACTTTTCAAA";
    //                         X    X

    auto graph = build_graph_batch<TypeParam>(k, references);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);

    run_alignment(*graph, config, query, { "" }, { "6=1X4=1X3=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_multiple_misalignment) {
    size_t k = 5;
    std::string reference = "AAAAGCGGACCCTTTCCGTTAT";
    std::string query =     "AAAAGGGGACCCTTTTCGTTAT";
    //                            X         X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X9=1X6=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "5=1X9=1X6=" });
    }
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
    run_alignment(*graph, config, query, { reference }, { "5=1I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=1I5=" });
    }
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
    run_alignment(*graph, config, query, { reference }, { "5=2I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=2I5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long) {
    size_t k = 4;
    std::string reference = "TTTCCTTGTT";
    std::string query =     "TTTCCAAAAAAAAATTGTT";
    //                            IIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference }, { "5=9I5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 2;
        run_alignment(*graph, config, query, { reference }, { "5=9I5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_insert_long_offset) {
    size_t k = 5;
    std::string reference = "TTTCCGGTTGTTA";
    std::string query =     "TTTCCGCAAAAAAAAATTGTTA";
    //                             XIIIIIIIII

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.gap_opening_penalty = -1;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference }, { "6=1X9I6=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_delete) {
    size_t k = 5;
    //                               D
    std::string reference = "TTCGAT""A""GGCCT";
    std::string query =     "TTCGAT"   "GGCCT";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "6=1D5=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { reference }, { "6=1D5=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_gap) {
    size_t k = 10;
    //                                  DDDD
    std::string reference = "TTTTCTGTATACCTTGGCGCTCTCC";
    std::string query =     "TTTTCTGTATA"  "GGCGCTCTCC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "11=4D10=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference }, { "11=4D10=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_gap_after_seed) {
    size_t k = 4;
    //                           DDDD
    std::string reference = "TTTCCCAAGGCGCTCTC";
    std::string query =     "TTTC"  "GGCGCTCTC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -2;
    config.left_end_bonus = 4;
    config.right_end_bonus = 4;
    run_alignment(*graph, config, query, { reference }, { "4=4D9=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_loop_deletion) {
    size_t k = 4;
    std::string reference = "AAAATTTTCGAGGCCAA";
    std::string query =     "AAAACGAGGCCAA";
    //                           DDDD

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::unit_scoring_matrix(
        2, alphabets::kAlphabetDNA, alphabets::kCharToDNA
    );
    config.gap_opening_penalty = -2;
    config.gap_extension_penalty = -2;
    run_alignment(*graph, config, query, { "AAAATTTCGAGGCCAA" }, { "4=3D9=" }, 0, true, true);
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

TYPED_TEST(DBGAlignerRedoneTest, align_drop_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTGGCGCTCTC";
    std::string query =     "TTTCCGGGGCGCTCTC";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -6, -6);
    config.gap_opening_penalty = -10;
    config.gap_extension_penalty = -4;
    config.xdrop = 6;
    run_alignment(*graph, config, query, { reference.substr(7) }, { "7S9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_long_gap_after_seed) {
    size_t k = 4;
    std::string reference = "TTTCCCTTAAGGCGCTCTC";
    std::string query =     "TTTC"    "GGCGCTCTC";
    //                       SSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.gap_opening_penalty = -5;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference.substr(10) }, { "4S9=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_repeat_sequence_no_delete_after_insert) {
    size_t k = 27;
    std::string reference = "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "TTTGTGGCTAGAGCTCGAGATCGCGCGGCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                                                    IIIIIII        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    run_alignment(*graph, config, query, { reference }, { "45=7I8=1X39=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 25;
        run_alignment(*graph, config, query, { reference }, { "45=7I8=1X39=" }, 0, true, true);
    }
}

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

TYPED_TEST(DBGAlignerRedoneTest, align_long_clipping) {
    size_t k = 4;
    std::string reference = "TTTTTTTAAAAGCTTCGAGGCCAA";
    std::string query =     "CCCCCCCAAAAGCTTCGAGGCCAA";
    //                       SSSSSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(7) }, { "7S17=" });
}

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

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity) {
    size_t k = 27;
    std::string reference = "CTAGAACTTAAAGTATAATAATACTAATAATAAAATAAAATACA";
    std::string query =     "CTAGAACTTAAAGTATAATAATACTAATAA""AAGTACAATACA";
    //                                                     DD  X  X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, { reference }, { "30=2D2=1X2=1X6=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity2) {
    size_t k = 27;
    std::string reference = "GCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string query =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                         IIIIIII        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);

    run_alignment(*graph, config, query, { reference.substr(3) }, { "10S4=1X9=1X8=1X39=" }, 0, true, true);

    config.gap_opening_penalty = -3;
    run_alignment(*graph, config, query, { reference },           {       "18=7I8=1X39=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity2_del) {
    size_t k = 27;
    std::string query =         "GCCACAATTGACAAATGA"     "GATCTAATTAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    std::string reference =     "GCCACAATTGACAAATGACAAATGTGATCTAATGAAACTAAAGAGCTTCTGCACAGCAAAAGAAACTGTCATC";
    //                                             DDDDDDD        X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);

    config.gap_opening_penalty = -5;
    run_alignment(*graph, config, query, { reference.substr(10) }, { "3S4=1X9=1X8=1X39=" }, 0, true, true);

    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -1;
    run_alignment(*graph, config, query, { reference },           {       "18=7D8=1X39=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4) {
    GTEST_SKIP() << "Runs for minutes: the engine has no work budget";
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

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    config.gap_opening_penalty = -5;
    config.gap_extension_penalty = -2;
    run_alignment(*graph, config, query, { "" }, { "" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity5) {
    size_t k = 31;
    std::string reference = "GTCGTCAGATCGGAAGAGCGTCGTGTAGGGAAAG" "GTCTTC""GCCTGTGTAGATCTC";
    std::string query =        "GTCAGATCGGAAGAGCGTCGTGTAGGGAAAG""AGTGTTCCTGGTGGTGTAGATCTC";
    //                                                           I  X   II XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 11;
        run_alignment(*graph, config, query, { reference.substr(3, 52) }, { "31=1I2=1X3=2I1=3X11=" }, 0, true, true);
    }
}

#if ! _PROTEIN_GRAPH

TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4_rep_primary) {
    GTEST_SKIP() << "Runs for minutes: the engine has no work budget";
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

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;
    config.gap_opening_penalty = -5;
    config.gap_extension_penalty = -2;
    run_alignment(*graph, config, query, { "" }, { "" });
}
#endif

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, align_both_directions) {
    size_t k = 7;
    std::string reference =    "AAAAGCTTTCGAGGCCAA";
    std::string query =        "AAAAGTTTTCGAGGCCAA";
    //                               X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::CANONICAL);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X12=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference }, { "5=1X12=" });
    }
}
#endif

TYPED_TEST(DBGAlignerRedoneTest, align_both_directions2) {
    size_t k = 11;
    std::string reference =    "GTAGTGCTAGCTTGTAGTCGTGCTGATGCA";
    std::string query =        "GTAGTGCTACCTTGTAGTCGTGGTGATGCA";
    //                                   X            X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::BASIC);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(0, query.size()) }, { "9=1X12=1X7=" }, 0, true, true);
}

TYPED_TEST(DBGAlignerRedoneTest, align_seed_to_end) {
    size_t k = 5;
    //                       SSSSSSSSSSSSSSSSSSSSSS
    std::string reference =                  "ATCCCTTTTAAAA";
    std::string query =     "ATCCCGGGGGGGGGGGGGGGGGTTTTAAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference.substr(5) }, { "22S8=" });
}

TYPED_TEST(DBGAlignerRedoneTest, align_bfs_vs_dfs_xdrop) {
    size_t k = 31;
    std::string reference_1 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGGATCAGAGAAGAAAGGAATGGCAAATACTCCATAGAAGAGTACCTTCAGGACAGGGACAGATACTATGAGGAGGTGGCCAT";
    std::string reference_2 = "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGAAGAGTACGATGAGTACAAGAGAATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";
    std::string query =       "TCGGGGCAAGAAACACACAGCCTTCTCATCCAAGGGCCTCAGTGATGATGAGTACGATGAGTACAAGAGCATCAGAGAGGAGAGGAATGGCAAATACTCAATAGAGGAATACCTCCAAGATAGGGACAGATACTATGAAGAGCTTGCCAT";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.xdrop = 27;
    run_alignment(*graph, config, query, { "" }, { "48=1X20=1X80=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 19;
        run_alignment(*graph, config, query, { "" }, { "48=1X20=1X80=" });
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy_short) {
    size_t k = 7;
    std::string reference = "GAAAAGCT";
    std::string query =     "GAAAAGTT";
    //                             X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference }, { "6=1X1=" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy_short2) {
    size_t k = 31;
    std::string reference = "GCTAGTCTATAGCGTACGATCTAGCTAGTCT";
    std::string query =     "GCTAGTCTATAGCGTACGATCTTTTTAGTC";
    //                                             XXX

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 19;
        run_alignment(*graph, config, query, { reference.substr(0, 30) }, { "22=3X5=" }, 1, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_dummy) {
    size_t k = 7;
    std::string reference = "AAAAGCTTTCGATT";
    std::string query =     "AAAAGTTTTCGAA";
    //                            X      X

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.left_end_bonus = 2;
    config.right_end_bonus = 2;
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 5;
        run_alignment(*graph, config, query, { reference.substr(0, query.size()) }, { "5=1X6=1X" }, 0, true, true);
    }
}

TYPED_TEST(DBGAlignerRedoneTest, align_extended_insert_after_match) {
    size_t k = 27;
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG" "GCCCAGGCCCAGGC""CCA" "GGCCCAGGC""CC" "AGGCCCAGGCCCAGGCCCAAGCC";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG" "CCCCAGGCCCAGGC""CCA" "GGCCCAGGC""CC" "AGGCCCAGGCCCAGGCCCAAGCC";
    std::string query =       "CGTGGCCCAGGCCCAGGCCCAGTGGGCGTTGGCCCAGGCGGCCA""CGG" "TGGCTGCG""CAGGCCCGCCTGGCACAAGCCACGCTG";
    //                                               III  XXX         II     I  DDDX   II X  I

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, {}, {});

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 15;
        run_alignment(*graph, config, query,
                      { "CGTGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAGGCCCAAGCC" },
                      { "" }, 0, true, true);
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, align_suffix_seed_no_full_seeds) {
    size_t k = 31;
    std::string reference =     "CTGCTGCGCCATCGCAACCCACGGTTGCTTTTTGAGTCGCTGCTCACGTTAGCCATCACACTGACGTTAAGCTGGCTTTCGATGCTGTATCTTTTTTTT";
    std::string query     = "CTTACTGCTGCGCTCTTCGCAAACCCCACGGTTTCTTGTTTTGAGCTCGCCTGCTCACGATACCCATACACACTGACGTTCAAGCTGGCTTTCGATGTTGTATC";
    //                       SSSS

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::PRIMARY);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    run_alignment(*graph, config, query, {}, {});
    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
        config.min_seed_length = 9;
        run_alignment(*graph, config, query,
                      { reference.substr(0, 91) },
                      { "4S9=1I1=1X6=2I9=1X2=2I7=1I3=1I9=1X2=1X4=1I12=1I16=1X6=" }, 0, true, true);
    }
}
#endif

// Records the connect flag that align_query passes to the seeder
class RecordingSeeder : public ExactSeeder {
  public:
    using ExactSeeder::ExactSeeder;

    std::vector<Alignment> get_inexact_anchors(bool align = true) const override {
        flags.push_back(align);
        return ExactSeeder::get_inexact_anchors(align);
    }

    mutable std::vector<bool> flags;
};

TYPED_TEST(DBGAlignerRedoneTest, align_query_forwards_connect_anchors) {
    size_t k = 10;
    std::string reference = "TTTTCTGTATACCTTGGCGCTCTCC";
    std::string query =     "TTTTCTGTATA"  "GGCGCTCTCC";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.gap_opening_penalty = -3;
    config.gap_extension_penalty = -3;
    config.min_seed_length = k;
    config.max_seed_length = k;
    Query aln_query(*graph, query);
    Extender extender(aln_query, config);

    for (bool connect : { true, false }) {
        RecordingSeeder seeder(aln_query, config);
        std::vector<Alignment> paths;
        auto callback = [&](Alignment&& path) { paths.emplace_back(std::move(path)); };

        align_query(aln_query, seeder, extender, callback, connect);
        EXPECT_EQ(std::vector<bool>{ connect }, seeder.flags) << "extend";
        EXPECT_FALSE(paths.empty()) << connect;

        seeder.flags.clear();
        align_query(aln_query, seeder, callback, connect);
        EXPECT_EQ(std::vector<bool>{ connect }, seeder.flags) << "no extend";
    }
}

TYPED_TEST(DBGAlignerRedoneTest, extend_unconnected_chain) {
    size_t k = 5;
    std::string reference = "AGCAACTCGAAA";
    std::string query =     "AGCAATTCGAAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference });
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    // Anchors are chained across the mismatch with sub-k seeds (DBGSuccinct)
    // or by SSHash's own chaining
    config.min_seed_length = std::is_base_of_v<DBGSuccinct, TypeParam> ? 3 : k;
    config.max_seed_length = k;
    Query aln_query(*graph, query);
    ExactSeeder seeder(aln_query, config);
    Extender extender(aln_query, config);

    bool has_placeholder = false;
    for (const auto &chain : seeder.get_inexact_anchors(false)) {
        const auto &path = chain.get_path();
        if (std::find(path.begin(), path.end(), DeBruijnGraph::npos) == path.end())
            continue;

        has_placeholder = true;
        std::vector<Alignment> paths;
        extender.extend(chain, [&](Alignment&& path) { paths.emplace_back(std::move(path)); });
        ASSERT_EQ(1u, paths.size());
        EXPECT_EQ(chain.get_cigar().to_string(), paths[0].get_cigar().to_string());
        EXPECT_EQ(chain.get_path(), paths[0].get_path());
    }
    if constexpr(!std::is_same_v<DBGHashFast, TypeParam>) {
        EXPECT_TRUE(has_placeholder);
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(DBGAlignerRedoneTest, query_reverse_complement) {
    auto graph = build_graph_batch<TypeParam>(4, { "AGCTTCGAGGCCAA" });
    Query query(*graph, std::string("AACCG"));
    EXPECT_EQ("AACCG", query.get_query(false));
    EXPECT_EQ("CGGTT", query.get_query(true));
}

TYPED_TEST(DBGAlignerRedoneTest, align_reverse_complement_basic) {
    size_t k = 5;
    std::string reference = "AGCTTCGAGGCCAATCGT";
    std::string query = reference;
    ::reverse_complement(query.begin(), query.end());

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::BASIC);
    ASSERT_EQ(DeBruijnGraph::BASIC, graph->get_mode());
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    config.min_seed_length = k;
    config.max_seed_length = k;
    Query aln_query(*graph, query);

    ExactSeeder seeder(aln_query, config);
    Extender extender(aln_query, config);
    std::vector<Alignment> paths;
    align_query(aln_query, seeder, extender,
                [&](Alignment&& path) { paths.emplace_back(std::move(path)); });
    ASSERT_FALSE(paths.empty());
    auto best = std::max_element(paths.begin(), paths.end(),
        [&](const auto &a, const auto &b) {
            return score_match(a, config) < score_match(b, config);
        });
    EXPECT_TRUE(best->get_orientation());
    EXPECT_EQ(reference, best->get_spelling());
    EXPECT_EQ("18=", best->get_cigar().to_string());
}
#endif

} // namespace
