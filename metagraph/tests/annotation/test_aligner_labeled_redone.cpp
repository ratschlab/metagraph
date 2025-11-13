#include <gtest/gtest.h>

#include "../graph/all/test_dbg_helpers.hpp"
#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include <tsl/hopscotch_set.h>

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

template <typename T>
struct template_parameter;

template <template <typename...> class C, typename T>
struct template_parameter<C<T>> {
    using type = T;
};

template <typename T>
using get_kmer_t = typename template_parameter<std::decay_t<T>>::type;


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

        auto aln_sort = [&config,&anno_buffer](const auto &a, const auto &b) {
            return std::make_tuple(score_match(a, config),
                                   anno_buffer.get_cached_column_set(b.get_label_classes()[0])[0])
                 > std::make_tuple(score_match(b, config),
                                   anno_buffer.get_cached_column_set(a.get_label_classes()[0])[0]);
        };

        std::sort(paths_no_extend.begin(), paths_no_extend.end(), aln_sort);
        std::sort(paths_no_align.begin(), paths_no_align.end(), aln_sort);
        std::sort(paths.begin(), paths.end(), aln_sort);

        ASSERT_LE(mappings.size(), paths.size()) << mx;
        paths.resize(mappings.size());

        tsl::hopscotch_set<std::string> seen_labels;

        for (size_t i = 0; i < mappings.size(); ++i) {
            const auto &[label, reference, cigar_str] = mappings[i];
            Cigar cigar(cigar_str);

            auto check_ref = [&](const Alignment &path, const std::string &type) {
                const auto &label = std::get<0>(mappings[i]);
                std::string reference = std::get<1>(mappings[i]);
                std::string cigar_str = std::get<2>(mappings[i]);

                if (path.get_orientation() && !path.get_end_trim()) {
                    ::reverse_complement(reference.begin(), reference.end());
                    cigar = Cigar(cigar_str);
                    std::reverse(cigar.data().begin(), cigar.data().end());
                    cigar_str = cigar.to_string();
                }

                if (label.size()) {
                    // make sure the labels are correct before moving on
                    size_t target_column = label_encoder.encode(label);

                    size_t old_num_column_sets = anno_buffer.num_column_sets();
                    size_t target = anno_buffer.cache_column(target_column);
                    ASSERT_EQ(old_num_column_sets, anno_buffer.num_column_sets());
                    ASSERT_LE(target, old_num_column_sets);

                    const auto &path_targets = path.get_label_classes();
                    std::string path_target_labels = anno_buffer.generate_column_set_str(path_targets[0], path.get_spelling().size());

                    if (path_targets.size() > 1) {
                        ASSERT_TRUE(std::all_of(path_targets.begin(), path_targets.end(),
                                                [&](auto node_target) { return node_target == path_targets[0]; }))
                            << path << "\t" << path_target_labels << "\t" << fmt::format("{}", fmt::join(path_targets, ","));
                    }
                    // for (auto node_target : path_targets) {
                    //     for (auto column : anno_buffer.get_cached_column_set(node_target)) {
                    //         path_target_labels += label_encoder.decode(column) + ",";
                    //     }
                    //     path_target_labels += ";";
                    // }
                    ASSERT_TRUE(std::all_of(path_targets.begin(), path_targets.end(),
                                            [&](auto node_target) { return node_target == target; }))
                        << path << "\t" << path_target_labels << " vs. " << label << "\t"
                        << fmt::format("{}", fmt::join(path_targets, ",")) << " vs. " << target;
                }

                ASSERT_EQ(cigar.to_string(), path.get_cigar().to_string()) << label << "\t" << mx << "\t" << type;
                if (reference.size()) {
                    ASSERT_TRUE(cigar.is_valid(reference, path.get_query()));
                    EXPECT_EQ(reference.size() - k + 1 + end_trim, path.get_path().size()) << label << "\t" << mx << "\t" << type;
                    EXPECT_EQ(reference, path.get_spelling()) << label << "\t" << mx << "\t" << type;
                }

                EXPECT_EQ(end_trim, path.get_end_trim()) << label << "\t" << mx << "\t" << type;
                if (reference.size()) {
                    EXPECT_EQ(config.score_cigar(reference, path.get_query(), cigar),
                              score_match(path, config)) << label << "\t" << mx << "\t" << type;
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
                         std::pair<DBGSSHash, annot::ColumnCompressed<>>,
                         std::pair<DBGHashFast, annot::RowFlatAnnotator>,
                         std::pair<DBGSuccinct, annot::RowFlatAnnotator>,
                         std::pair<DBGSSHash, annot::RowFlatAnnotator>> FewGraphAnnotationPairTypes;

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

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));

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

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));
    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("B"),
              anno_graph->get_annotator().get_label_encoder().encode("C"));

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -2, -2);
    config.gap_opening_penalty = -4;
    config.gap_extension_penalty = -4;
    config.left_end_bonus = 6;
    config.right_end_bonus = 6;

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("CGAATGCAT"), { // RC: ATGCATTCG
            std::make_tuple(std::string("B"), std::string("CGAATGCCT"), std::string("7=1X1=")),
            std::make_tuple(std::string("C"), std::string("GGAATGCAT"), std::string("1X8=")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings, 0, true, true);
    }
}

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraphSuffixSeed) {
    if constexpr(!std::is_base_of_v<DBGSuccinct, typename TypeParam::first_type>
                    && !std::is_base_of_v<DBGSSHash, typename TypeParam::first_type>) {
        return;
    }

    size_t k = 7;
    /*
         B       B       B       B       B       B               B       B      AB      AB
        TATTGCC-ATTGCCG-TTGCCGA-TGCCGAT-GCCGATT-CCGATTG         GATTGCC-ATTGCCT-TTGCCTT-TGCCTTG
                                                       \ BC    /
                                                        CGATTGC
           C      C       C       C       C       C    /       \  C       C       C       C
        TATTGGC-ATTGGCG-TTGGCGA-TGGCGAT-GGCGATT-GCGATTG         GATTGCA-ATTGCAT-TTGCATT-TGCATTG
    */
    const std::vector<std::string> sequences {
                 "TTGCCTTG",
        "TATTGCCGATTGCCTTG",
        "TATTGGCGATTGCATTG"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(k, sequences, labels);

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));
    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("B"),
              anno_graph->get_annotator().get_label_encoder().encode("C"));

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.min_seed_length = 5;
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;

    {
        std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
            { std::string("TATTGACGATTGCATTG"), {
                std::make_tuple(std::string("C"), std::string("TATTGGCGATTGCATTG"), std::string("5=1X11=")),
            } }
        };

        for (size_t i = 0; i < exp_alignments.size(); ++i) {
            bool first_seed_matches = true;
            if constexpr(std::is_same_v<DBGSSHash, typename TypeParam::first_type>) {
                const auto &sshash = static_cast<const DBGSSHash&>(anno_graph->get_graph());
                std::visit([&](const auto &dict) {
                    ASSERT_GE(config.min_seed_length, dict.m());
                    const auto &[query, mappings] = exp_alignments[i];
                    ASSERT_LT(0u, mappings.size());
                    using kmer_t = get_kmer_t<decltype(dict)>;
                    auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(std::get<1>(mappings[0]).data(), k);
                    auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                    // we can only find the first seed of it's also the minimizer
                    // if (sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_ref, k, dict.m(), dict.seed()).second
                            // || sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_query, k, dict.m(), dict.seed()).second) {
                        first_seed_matches = false;
                    // }
                }, sshash.data());
            }

            const auto &[query, mappings] = exp_alignments[i];
            if (first_seed_matches) {
                run_alignment(*anno_graph, config, query, mappings);
            } else {
                run_alignment(*anno_graph, config, query, mappings, 0, true, true);
            }
        }
    }
    {
        std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
            { std::string("TATTGCCGATTGCATTG"), {
                std::make_tuple(std::string("B"), std::string("TATTGCCGATTGCCTTG"), std::string("13=1X3=")),
                std::make_tuple(std::string("C"), std::string("TATTGGCGATTGCATTG"), std::string("5=1X11=")),
            } }
        };

        for (const auto &[query, mappings] : exp_alignments) {
            run_alignment(*anno_graph, config, query, mappings, 0, true, true);
        }
    }
}

#if ! _PROTEIN_GRAPH
TYPED_TEST(LabeledAlignerRedoneTest, CanonicalTangleGraph) {
    if constexpr(!std::is_base_of_v<DBGSuccinct, typename TypeParam::first_type>
                    && !std::is_base_of_v<DBGSSHash, typename TypeParam::first_type>) {
        return;
    }

    size_t k = 7;
    const std::vector<std::string> sequences {
                 "CAAGGCAA",
        "CAAGGCAATCGGCAATA",
        "CAATGCAATCGCCAATA"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k, sequences, labels, mode
        );

        ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
                  anno_graph->get_annotator().get_label_encoder().encode("B"));
        ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("B"),
                  anno_graph->get_annotator().get_label_encoder().encode("C"));

        DBGAlignerConfig config;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
        config.min_seed_length = 5;
        config.left_end_bonus = 5;
        config.right_end_bonus = 5;

        {
            std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
                { std::string("TATTGACGATTGCATTG"), {
                    std::make_tuple(std::string("C"), std::string("TATTGGCGATTGCATTG"), std::string("5=1X11=")),
                } }
            };

            for (size_t i = 0; i < exp_alignments.size(); ++i) {
                bool first_seed_matches = true;
                if constexpr(std::is_same_v<DBGSSHash, typename TypeParam::first_type>) {
                    const auto &sshash = static_cast<const DBGSSHash&>(anno_graph->get_graph());
                    std::visit([&](const auto &dict) {
                        ASSERT_GE(config.min_seed_length, dict.m());
                        const auto &[query, mappings] = exp_alignments[i];
                        ASSERT_LT(0u, mappings.size());
                        using kmer_t = get_kmer_t<decltype(dict)>;
                        auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(std::get<1>(mappings[0]).data(), k);
                        auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                        // we can only find the first seed of it's also the minimizer
                        // if (sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_ref, k, dict.m(), dict.seed()).second
                                // || sshash::util::compute_minimizer_pos<kmer_t>(first_kmer_query, k, dict.m(), dict.seed()).second) {
                            first_seed_matches = false;
                        // }
                    }, sshash.data());
                }

                const auto &[query, mappings] = exp_alignments[i];
                if (first_seed_matches) {
                    run_alignment(*anno_graph, config, query, mappings);
                } else {
                    run_alignment(*anno_graph, config, query, mappings, 0, true, true);
                }
            }
        }
        {
            std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
                { std::string("TATTGCCGATTGCATTG"), {
                    std::make_tuple(std::string("B"), std::string("TATTGCCGATTGCCTTG"), std::string("13=1X3=")),
                    std::make_tuple(std::string("C"), std::string("TATTGGCGATTGCATTG"), std::string("5=1X11=")),
                } }
            };

            for (const auto &[query, mappings] : exp_alignments) {
                run_alignment(*anno_graph, config, query, mappings, 0, true, true);
            }
        }
    }
}
#endif

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraphCoords) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C
       GGA                 GCA-CAT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCAT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));
    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("B"),
              anno_graph->get_annotator().get_label_encoder().encode("C"));

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.left_end_bonus = 5;
    config.right_end_bonus = 5;

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("CGAATGCAT"), { // RC: ATGCATTCG
            std::make_tuple(std::string("B"), std::string("CGAATGCCT"), std::string("7=1X1=")),
            std::make_tuple(std::string("C"), std::string("GGAATGCAT"), std::string("1X8=")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings, 0, true, true);
    }
}

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraphCoordsMiddle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 3;
    /*  B                  AB  AB
       CGA                 GCC-CCT
          \ BC  BC  BC ABC/
           GAA-AAT-ATG-TGC
         C/               \  C   C   C   C
       GGA                 GCA-CAT-ATT-TTT
    */
    const std::vector<std::string> sequences {
        "TGCCT",
        "CGAATGCCT",
        "GGAATGCATTTT"
    };
    const std::vector<std::string> labels { "A", "B", "C" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));
    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("B"),
              anno_graph->get_annotator().get_label_encoder().encode("C"));

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("CGAAAGCCTTTT"), { // RC: AAAGGCTTTCG
            std::make_tuple(std::string("C"), std::string("GAATGCATTTT"), std::string("1S3=1X2=1X4=")),
            std::make_tuple(std::string("B"), std::string("CGAATGCCT"), std::string("4=1X4=3S")),
            std::make_tuple(std::string("A"), std::string("GCCT"), std::string("5S4=3S")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings, 0, true, true);
    }
}

TYPED_TEST(LabeledAlignerRedoneTest, SimpleTangleGraphCoordsCycle) {
    // TODO: for now, not implemented for other annotators
    if constexpr(!std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>
                    && !std::is_same_v<typename TypeParam::second_type, annot::RowDiffColumnAnnotator>) {
        return;
    }

    size_t k = 4;
    /*
        A    A    B    B    B    B    B
        GCAA-CAAT-AATG-ATGC-TGCG-GCGC-CGCA
    */
    const std::vector<std::string> sequences {
        "GCAAT",
        "AATGCGCA"
    };
    const std::vector<std::string> labels { "A", "B" };

    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        k, sequences, labels, DeBruijnGraph::BASIC, true
    );

    ASSERT_LT(anno_graph->get_annotator().get_label_encoder().encode("A"),
              anno_graph->get_annotator().get_label_encoder().encode("B"));

    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);

    std::vector<std::pair<std::string, std::vector<std::tuple<std::string, std::string, std::string>>>> exp_alignments {
        { std::string("ATGCGCAATGCG"), { // RC: CGCATTGCGCAT
            std::make_tuple(std::string("B"), std::string("ATGCGCA"), std::string("7=5S")),
        } }
    };

    for (const auto &[query, mappings] : exp_alignments) {
        run_alignment(*anno_graph, config, query, mappings, 0, true, true);
    }
}

} // namespace