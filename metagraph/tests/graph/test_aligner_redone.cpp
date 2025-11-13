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
//const bool PICK_REV_COMP = true;

template <typename T>
struct template_parameter;

template <template <typename...> class C, typename T>
struct template_parameter<C<T>> {
    using type = T;
};

template <typename T>
using get_kmer_t = typename template_parameter<std::decay_t<T>>::type;

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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 3;
        run_alignment(*graph, config, query, { "AAAAACTTTTTTT" }, { "5=1D7=" });
    }

    // TODO: too many possibilities here
    // config.gap_opening_penalty = -5;
    // config.gap_extension_penalty = -2;
    // if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
    //     config.min_seed_length = 2;
    //     run_alignment(*graph, config, query, { "AAAACTTTTTT" }, { "4=1X6=" });
    // } else {
    //     run_alignment(*graph, config, query, { "AAAACTTTTTT" }, { "4=1X6=" }, 0, true, true);
    // }
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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
    // if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam>) {
    //     config.min_seed_length = 3;
    //     run_alignment(*graph, config, query, { reference }, { "6=1X9I6=" });
    // } else {
    // }
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

// TODO: double-check
// TYPED_TEST(DBGAlignerRedoneTest, align_long_gap_after_seed) {
//     size_t k = 4;
//     std::string reference = "TTTCCCTTAAGGCGCTCTC";
//     std::string query =     "TTTC"    "GGCGCTCTC";
//     //                       SSSS

//     auto graph = build_graph_batch<TypeParam>(k, { reference });
//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -1;
//     run_alignment(*graph, config, query, { reference.substr(10) }, { "4S9=" });
// }

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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

// TODO: implement min cell score
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

// TODO: implement seed complexity filter
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

//     DBGAlignerConfig config;
//     config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
//     config.left_end_bonus = 5;
//     config.right_end_bonus = 5;
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     // config.min_seed_length = 6;
//     run_alignment(*graph, config, query, { "" }, { "" });
// }

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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 11;
        run_alignment(*graph, config, query, { reference.substr(3, 52) }, { "31=1I2=1X3=2I1=3X11=" }, 0, true, true);
    }
}

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

#if ! _PROTEIN_GRAPH
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

// TYPED_TEST(DBGAlignerRedoneTest, align_low_similarity4_rep_primary) {
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
//     config.left_end_bonus = 5;
//     config.right_end_bonus = 5;
//     config.gap_opening_penalty = -5;
//     config.gap_extension_penalty = -2;
//     // config.min_seed_length = 6;
//     run_alignment(*graph, config, query, { "" }, { "" });
// }
#endif

TYPED_TEST(DBGAlignerRedoneTest, align_both_directions) {
    size_t k = 7;
    std::string reference =    "AAAAGCTTTCGAGGCCAA";
    std::string query =        "AAAAGTTTTCGAGGCCAA";
    //                               X

    auto graph = build_graph_batch<TypeParam>(k, { reference }, DeBruijnGraph::CANONICAL);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
    run_alignment(*graph, config, query, { reference }, { "5=1X12=" }, 0, true, true);

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 5;
        if constexpr(std::is_same_v<DBGSSHash, TypeParam>) {
            const auto &sshash = static_cast<const DBGSSHash&>(*graph);
            std::visit([&](const auto &dict) {
                using kmer_t = get_kmer_t<decltype(dict)>;
                auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(reference.data(), k);
                auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                // we can only find the first seed of it's also the minimizer
                // if (sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_ref, k, dict.m(), dict.seed()).second
                //         || sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_query, k, dict.m(), dict.seed()).second) {
                //     run_alignment(*graph, config, query, { reference }, { "5=1X12=" }, 0, true, true);
                // } else {
                    run_alignment(*graph, config, query, { reference }, { "5=1X12=" });
                // }
            }, sshash.data());
        } else {
            run_alignment(*graph, config, query, { reference }, { "5=1X12=" });
        }
    }
}

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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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

    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
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
    if constexpr(std::is_base_of_v<DBGSuccinct, TypeParam> || std::is_same_v<DBGSSHash, TypeParam>) {
        config.min_seed_length = 9;
        run_alignment(*graph, config, query,
                      { reference.substr(0, 91) },
                      { "4S9=1I1=1X6=2I9=1X2=2I7=1I3=1I9=1X2=1X4=1I12=1I16=1X6=" }, 0, true, true);
    }
}
#endif

} // namespace
