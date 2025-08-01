#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_match.hpp"
#include "graph/alignment_redone/aligner_config.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/annotated_dbg.hpp"
#include "kmer/alphabets.hpp"
#include "seq_io/sequence_io.hpp"

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align_redone;
using namespace mtg::test;
using namespace mtg::kmer;


void run_alignment(const AnnotatedDBG &anno_graph,
                   DBGAlignerConfig config,
                   std::string_view query,
                   const std::vector<std::tuple<std::string, std::string, std::string>> &mappings,
                   size_t end_trim = 0,
                   bool needs_extension = false,
                   bool needs_extension_long_seed = false) {
    const auto &graph = anno_graph.get_graph();
    const auto &label_encoder = anno_graph.get_annotator().get_label_encoder();
    AnnotationBuffer anno_buffer(graph, anno_graph.get_annotator());

    size_t k = graph.get_k();
    if (config.min_seed_length == 0)
        config.min_seed_length = k;

    for (auto mx : { k, std::numeric_limits<size_t>::max() }) {
        config.max_seed_length = std::max(mx, config.min_seed_length);
        bool check_chaining = mx == std::numeric_limits<size_t>::max()
            ? !needs_extension_long_seed
            : !needs_extension;

        Query aln_query(graph, query);
        LabeledSeeder seeder(anno_buffer, aln_query, config);
        std::vector<Alignment> paths;
        LabeledExtender extender(anno_buffer, aln_query, config);
        std::vector<Alignment> paths_no_extend = seeder.get_inexact_anchors();
        std::vector<Alignment> paths_no_align = seeder.get_inexact_anchors(false);
        for (const auto &base_path : paths_no_extend) {
            extender.extend(base_path, [&](Alignment&& path) {
                paths.emplace_back(std::move(path));
            });
        }

        auto aln_sort = [](const auto &a, const auto &b) {
            return std::make_pair(a.get_score(), b.get_orientation())
                 > std::make_pair(b.get_score(), a.get_orientation());
        };

        std::sort(paths_no_extend.begin(), paths_no_extend.end(), aln_sort);
        std::sort(paths_no_align.begin(), paths_no_align.end(), aln_sort);
        std::sort(paths.begin(), paths.end(), aln_sort);

        ASSERT_LE(mappings.size(), paths.size()) << mx;
        paths.resize(mappings.size());

        for (size_t i = 0; i < mappings.size(); ++i) {
            const auto &[label, reference, cigar_str] = mappings[i];
            Cigar cigar(cigar_str);

            auto check_ref = [&](const Alignment &path, const std::string &type) {
                const auto &[label, reference, cigar_str] = mappings[i];

                if (label.size()) {
                    size_t target_column = label_encoder.encode(label);
                    size_t target = anno_buffer.cache_column(target_column);
                    const auto &path_targets = path.get_label_classes();
                    std::string path_target_labels;
                    for (auto node_target : path_targets) {
                        for (auto column : anno_buffer.get_cached_column_set(node_target)) {
                            path_target_labels += label_encoder.decode(column) + ",";
                        }
                        path_target_labels += ";";
                    }
                    EXPECT_TRUE(std::all_of(path_targets.begin(), path_targets.end(),
                                            [&](auto node_target) { return node_target == target; }))
                        << path_target_labels << " vs. " << label;
                }

                ASSERT_EQ(cigar.to_string(), path.get_cigar().to_string()) << label << "\t" << mx << "\t" << type;
                if (reference.size()) {
                    ASSERT_TRUE(cigar.is_valid(reference, path.get_query()));
                    EXPECT_EQ(reference.size() - k + 1 + end_trim, path.get_path().size()) << label << "\t" << mx << "\t" << type;
                    EXPECT_EQ(reference, path.get_spelling()) << label << "\t" << mx << "\t" << type;
                }

                EXPECT_EQ(end_trim, path.get_end_trim()) << label << "\t" << mx << "\t" << type;
                if (reference.size()) {
                    EXPECT_EQ(config.score_cigar(reference, path.get_query(), cigar), path.get_score()) << label << "\t" << mx << "\t" << type;
                }

                EXPECT_EQ(cigar.get_clipping(), path.get_clipping()) << label << "\t" << mx << "\t" << type;
                EXPECT_EQ(cigar.get_end_clipping(), path.get_end_clipping()) << label << "\t" << mx << "\t" << type;
            };

            check_ref(paths[i], "extend");

            if (check_chaining) {
                // this alignment should work with chaining alone
                ASSERT_LT(i, paths_no_extend.size());
                check_ref(paths_no_extend[i], "chain");

                ASSERT_LT(i, paths_no_align.size());
                EXPECT_EQ(cigar.to_string(), paths_no_align[i].get_cigar().to_string()) << label << "\t" << mx << "\tnoalign";
            }
        }
    }
}

template <typename GraphAnnotationPair>
class LabeledAlignerRedoneTest : public ::testing::Test {};

typedef ::testing::Types<std::pair<DBGHashFast, annot::ColumnCompressed<>>,
                         std::pair<DBGSuccinct, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>> FewGraphAnnotationPairTypes;

TYPED_TEST_SUITE(LabeledAlignerRedoneTest, FewGraphAnnotationPairTypes);

TYPED_TEST(LabeledAlignerRedoneTest, SimpleLinearGraph) {
    size_t k = 4;
    /*
        A    A    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCT-GCTT
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCTT"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("GCAATGCTT"), {
            std::make_tuple(std::string("B"), std::string("AATGCTT"), std::string("2S7=")),
            std::make_tuple(std::string("A"), std::string("GCAAT"), std::string("5=4S")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings);
    }
}

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraph) {
    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC ABC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
           "ATGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -4;
    config.left_end_bonus = 6;
    config.right_end_bonus = 6;

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("CGAATGCAT"), { // RC: ATGCATTCG
            std::make_tuple(std::string("C"), std::string("GGAATGCAT"), std::string("1X8=")),
            std::make_tuple(std::string("B"), std::string("CGAATGCCT"), std::string("7=1X1=")),
            std::make_tuple(std::string("A"), std::string("ATGCCT"), std::string("4=1X1=3S")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings, 0, true, true);
    }
}

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraphSuffixSeed) {
    if constexpr(!std::is_base_of_v<DBGSuccinct, typename TypeParam::first_type>) {
        return;
    }

    size_t k = 4;
    /*   B    B    B                  AB   AB
        GTCG-TCGA-CGAA                TGCC-GCCT
                      \ BC   BC   BC /
                       GAAT-AATG-ATGC
           C    C    C/              \   C    C
        GTGG-TGGA-GGAA                TGCA-GCAT
    */
    const std::vector<std::string> sequences {
         "TGCCT",
        "GTCGAATGCCT",
        "GTGGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.min_seed_length = 2;
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;

    {
        std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
            { std::string("GTGAAATGCAT"), {
                std::make_tuple(std::string("C"), std::string("GTGGAATGCAT"), std::string("3=1X7=")),
            } }
        };

        for (const auto &[query, mappings] : exp_alignments) {
            run_alignment(*anno_graph, config, query, mappings);
        }
    }
    {
        std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
            { std::string("GTGAAATGCAT"), {
                std::make_tuple(std::string("C"), std::string("GTGGAATGCAT"), std::string("3=1X7=")),
                std::make_tuple(std::string("B"), std::string("GTCGAATGCCT"), std::string("2=2X5=1X1=")),
            } }
        };

        for (const auto &[query, mappings] : exp_alignments) {
            run_alignment(*anno_graph, config, query, mappings, 0, true, true);
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(LabeledAlignerRedoneTest, CanonicalTangleGraph) {
    if constexpr(!std::is_base_of_v<DBGSuccinct, typename TypeParam::first_type>) {
        return;
    }

    size_t k = 5;
    /*   B     B                AB    AB
       TTAGT-TAGTC             TCGAA-CGAAA
                  \  BC   ABC /
                   AGTCG-GTCGA
          C     C /           \   C     C
       TCAGT-CAGTC             TCGAT-CGATT
        AB    AB                 B     B
       TTTCG-TTCGA             GACTA-ACTAA
                  \ ABC    BC /
                   TCGAC-CGACT
          C     C /           \   C     C
       AATCG-ATCGA             GACTG-ACTGA
    */
    const std::vector<std::string> sequences {
        "GTCGAAA", // "TTTCGAC"
        "TTAGTCGAAA", // "TTTCGACTAA"
        "TCAGTCGATT" // "AATCGACTGA"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k, sequences, labels, mode
        );

        DBGAlignerConfig config;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
        config.min_seed_length = 4;
        std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
            { std::string("TTAGTTCAAA"), { // RC: TTTGAACTAA
                std::make_tuple(std::string("B"), std::string("TTAGTCGAAA"), std::string("5=2X3=")),
            } }
        };

        for (const auto &[query, mappings] : exp_alignments) {
            run_alignment(*anno_graph, config, query, mappings, 0, true, true);
        }
    }
}
#endif

} // namespace