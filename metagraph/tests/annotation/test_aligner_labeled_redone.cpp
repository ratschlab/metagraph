#include <gtest/gtest.h>
#include <map>
#include <random>
#include <set>

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
#include "graph/alignment/aligner_labeled.hpp"
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

TEST(LabeledFastAlignment, RelativeCutoffUsesEachLabelsBestScore) {
    const std::string sequence = "AACGTTGCAATGCCGATTCAGGCTAACGTCTAGCGTACGA";
    std::vector<std::string> references{sequence, sequence};
    references[0][20] = references[1][20] = 'A';
    references[1][31] = 'C';
    for (bool coordinates : {false, true}) {
        auto annotated = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            7, references, {"first", "second"}, DeBruijnGraph::BASIC, coordinates);
        const auto &graph = annotated->get_graph();
        AnnotationBuffer buffer(graph, annotated->get_annotator());
        DBGAlignerConfig config;
        config.min_seed_length = 7;
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.forward_and_reverse_complement = false;
        config.seed_complexity_filter = false;
        config.rel_score_cutoff = 0.95;
        config.xdrop = 10000;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
        Query query(graph, sequence);
        LabeledSeeder seeder(buffer, query, config);
        LabeledExtender extender(buffer, query, config);
        std::map<std::string, int> scores;
        align_query(query, seeder, extender, [&](Alignment&& alignment) {
            EXPECT_TRUE(alignment.is_spelling_valid(graph));
            for (auto column : buffer.get_cached_column_set(alignment.get_label_classes()[0])) {
                const auto &label = annotated->get_annotator().get_label_encoder().decode(column);
                scores[label] = score_match(alignment, config);
                EXPECT_EQ(references[label == "first" ? 0 : 1], alignment.get_spelling());
            }
        });
        EXPECT_EQ((std::map<std::string, int>{{"first", 75}, {"second", 70}}), scores)
            << coordinates;
    }
}

TEST(LabeledFastAlignment, ReusesAnnotationsWithIndependentSearchConfiguration) {
    const std::vector<std::string> references{
        "AACGTTGCAATGCCGATTCAGGCTAAC", "TTGATAGCCTTGACCGATAGCTTCCA"
    };
    for (bool coordinates : {false, true}) {
        auto annotated = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            5, references, {"first", "second"}, DeBruijnGraph::BASIC, coordinates);
        const auto &graph = annotated->get_graph();
        AnnotationBuffer shared(graph, annotated->get_annotator());
        for (size_t round = 0; round < 4; ++round) {
            DBGAlignerConfig config;
            config.min_seed_length = 5;
            config.max_seed_length = std::numeric_limits<size_t>::max();
            config.forward_and_reverse_complement = false;
            config.seed_complexity_filter = false;
            config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2 + round, -3, -3);
            std::string sequence = references[round % 2];
            sequence.at(18) = sequence.at(18) == 'A' ? 'C' : 'A';
            Query query(graph, sequence);
            using Result = std::tuple<int, std::string, std::string,
                                      AnnotationBuffer::Columns, AnnotationBuffer::CoordinateSet>;
            auto run = [&](AnnotationBuffer &buffer) {
                std::vector<Result> results;
                LabeledSeeder seeder(buffer, query, config);
                LabeledExtender extender(buffer, query, config);
                align_query(query, seeder, extender, [&](Alignment&& alignment) {
                    results.emplace_back(score_match(alignment, config), alignment.get_spelling(),
                        alignment.get_cigar().to_string(),
                        buffer.get_cached_column_set(alignment.get_label_classes()[0]),
                        alignment.get_coords());
                });
                return results;
            };
            AnnotationBuffer fresh(graph, annotated->get_annotator());
            auto expected = run(fresh);
            ASSERT_FALSE(expected.empty());
            EXPECT_EQ(expected, run(shared)) << coordinates << " " << round;
        }
        mtg::graph::align::DBGAlignerConfig legacy{};
        legacy.min_seed_length = legacy.max_seed_length = 5;
        legacy.score_matrix = mtg::graph::align::DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
        for (bool chain : {false, true}) {
            for (bool two_sided : {false, true}) {
                legacy.chain_alignments = chain;
                legacy.two_sided = two_sided;
                mtg::graph::align::LabeledAligner<> aligner(graph, legacy, annotated->get_annotator());
                EXPECT_EQ(chain, aligner.get_config().chain_alignments);
                EXPECT_EQ(two_sided, aligner.get_config().two_sided);
            }
        }
    }
}

#if ! _PROTEIN_GRAPH
TEST(LabeledFastAlignment, XdropDoesNotPruneOtherLabels) {
    const std::vector<std::string> references{
        "AGCCCCTACACTACACTACAATTATGTGTGGCAACGCATGTGTCGGCATTATGGCTGACACATGGGGAAT"
        "TGGCTTTCTCCATTCGGTGAGCCAATTGGCGTTTGCCGTGCACTTACTCTTCTGT",
        "AGCCCCTACACTACACTACAATTATGTGTGGCAACGCATGTGTCGGCATTATGGCTGTCACATGGGGAAT"
        "AGGCTTTCTCCATTCGGTGAGCCAGTTGGCGTTTGCCATGCACTTACTCTTCTGT"
    };
    const std::string sequence =
        "AGCCCCTACACTACACTACAATTATGTGTGGTAACGCATGAGTCGGCATTATGGCTGTCACATGGGGAAT"
        "AGGGCTTTCTCCATTCGGTGAGCCAGTTGGCGTTTGCCATGCACTTACTCTTCTGT";
    for (bool coordinates : {false, true}) {
        auto annotated = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            15, references, {"first", "second"}, DeBruijnGraph::BASIC, coordinates);
        const auto &graph = annotated->get_graph();
        AnnotationBuffer buffer(graph, annotated->get_annotator());
        DBGAlignerConfig config;
        config.min_seed_length = 15;
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.forward_and_reverse_complement = false;
        config.rel_score_cutoff = 0;
        config.xdrop = 27;
        config.gap_opening_penalty = -6;
        config.left_end_bonus = config.right_end_bonus = 5;
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
        for (auto [nodes, ram] : {std::pair{5.0, 200.0}, {0.1, 200.0}, {5.0, 0.001}}) {
            SCOPED_TRACE(::testing::Message() << coordinates << " " << nodes << " " << ram);
            const bool full_search = nodes == 5.0 && ram == 200.0;
            config.max_nodes_per_seq_char = nodes;
            config.max_ram_per_alignment = ram;
            // Let capped searches report the exact seeds they retain.
            config.min_exact_match = full_search ? 0.7 : 0;
            Query query(graph, sequence);
            LabeledSeeder seeder(buffer, query, config);
            LabeledExtender extender(buffer, query, config);
            std::map<std::string, int> scores;
            align_query(query, seeder, extender, [&](Alignment&& alignment) {
                EXPECT_TRUE(alignment.is_spelling_valid(graph));
                EXPECT_TRUE(alignment.get_cigar().is_valid(alignment.get_spelling(), sequence));
                for (auto column : buffer.get_cached_column_set(alignment.get_label_classes()[0])) {
                    const auto &label = annotated->get_annotator().get_label_encoder().decode(column);
                    EXPECT_TRUE(scores.emplace(label, score_match(alignment, config)).second);
                }
            });
            if (full_search) {
                EXPECT_EQ((std::map<std::string, int>{{"first", 224}, {"second", 244}}), scores);
            } else {
                EXPECT_EQ((std::map<std::string, int>{{"first", 67}, {"second", 113}}), scores);
            }
        }
    }
}

TEST(LabeledAlignment, PartialNodeSuffixKeepsCoordinates) {
    const std::string reference = "CGCCACCGGCCGACCGACCGCGGCTGCCCCAGTTGATG";
    auto annotated = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        31, {reference}, {"reference"}, DeBruijnGraph::BASIC, true);
    const auto &graph = annotated->get_graph();
    AnnotationBuffer buffer(graph, annotated->get_annotator());
    const auto path = map_to_nodes_sequentially(graph, reference);
    buffer.queue_path(path);
    buffer.fetch_queued_annotations();
    DBGAlignerConfig config;
    config.min_seed_length = 19;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.min_exact_match = 0.7;
    config.xdrop = 27;
    config.gap_opening_penalty = -6;
    config.left_end_bonus = config.right_end_bonus = 5;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    auto check = [&](const Alignment &alignment, const char *cigar, int score) {
        EXPECT_EQ(reference, alignment.get_spelling());
        EXPECT_EQ(cigar, alignment.get_cigar().to_string());
        EXPECT_EQ(score, score_match(alignment, config));
        EXPECT_TRUE(alignment.is_spelling_valid(graph));
        const auto &columns = buffer.get_cached_column_set(alignment.get_label_classes()[0]);
        ASSERT_EQ(1u, columns.size());
        EXPECT_EQ("reference", annotated->get_annotator().get_label_encoder().decode(columns[0]));
        ASSERT_EQ(1u, alignment.get_coords().size());
        EXPECT_EQ(AnnotationBuffer::Tuple({0}), alignment.get_coords()[0]);
    };
    Query exact(graph, reference);
    Anchor seed(exact.get_query(), 0, 19, false, {path.front()}, reference.substr(19, 12));
    seed.set_label_class(buffer.get_labels_id(path.front()));
    seed.set_coord(0);
    size_t count = 0;
    // Consume the remaining first-node suffix before advancing reference coordinates.
    LabeledExtender(buffer, exact, config).extend(Alignment(seed), [&](Alignment&& alignment) {
        ++count;
        check(alignment, "38=", 86);
    }, true);
    EXPECT_EQ(1u, count);

    Query query(graph, "CGCCACCGGCCGACCGACCGCGCGCTGCCACCAGTTGATG");
    LabeledSeeder seeder(buffer, query, config);
    LabeledExtender extender(buffer, query, config);
    count = 0;
    // The default route must retry partial-node seeds when the column search misses.
    align_query(query, seeder, extender, [&](Alignment&& alignment) {
        ++count;
        check(alignment, "22=1I6=1I10=", 74);
    });
    EXPECT_EQ(1u, count);
}
#endif

TEST(LabeledChaining, SharedSeedsKeepDistinctConnectingPaths) {
    const std::vector<std::string> references{
        "TCCCATTACATCCCCCCACTCATCCTCCATGTTCCCCCCACTCATAGACCCCAAC",
        "TCCCATTACATCCCCCCTCTCATCCATCTACATGTTCCTCCCACTCATAGACCCCAAC",
        "TCCCATTACATCCCCCCAGTCATCCTTGTTCCCACCACTCATAGACCCCAAC"
    };
    const std::string sequence = "TCCCATTACATCCCCCCACTAACCTCCTCCATTTTCCCCTCTCATAGAACCCAAC";
    auto annotated = build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
        3, references, {"0", "1", "2"});
    const auto &graph = annotated->get_graph();
    const auto &encoder = annotated->get_annotator().get_label_encoder();
    AnnotationBuffer buffer(graph, annotated->get_annotator());
    DBGAlignerConfig config;
    config.min_seed_length = 3;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.chain_alignments = true;
    config.two_sided = true;
    config.forward_and_reverse_complement = false;
    config.seed_complexity_filter = false;
    config.min_exact_match = 0;
    config.rel_score_cutoff = 0;
    config.num_alternative_paths = 3;
    config.xdrop = 27;
    config.gap_opening_penalty = config.gap_extension_penalty = -2;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    Query query(graph, sequence);
    LabeledSeeder seeder(buffer, query, config);
    LabeledExtender extender(buffer, query, config);
    std::set<std::string> found_labels;
    align_query(query, seeder, extender, [&](Alignment&& alignment) {
        EXPECT_TRUE(alignment.is_spelling_valid(graph));
        EXPECT_TRUE(alignment.get_cigar().is_valid(alignment.get_spelling(), alignment.get_query()));
        for (auto column : buffer.get_cached_column_set(alignment.get_label_classes()[0])) {
            const auto &label = encoder.decode(column);
            found_labels.insert(label);
            const auto &reference = references.at(std::stoul(label));
            // Column annotations constrain membership in a label's k-mer
            // subgraph. Reference-order continuity needs coordinate annotations.
            for (auto node : alignment.get_path()) {
                ASSERT_NE(DeBruijnGraph::npos, node); // anchor connection is enabled
                EXPECT_NE(std::string::npos, reference.find(graph.get_node_sequence(node)))
                    << "label " << label << " in " << alignment;
            }
        }
    });
    EXPECT_EQ((std::set<std::string>{"0", "1", "2"}), found_labels);
}

template <class Graph>
class MatchFractionAlignment : public ::testing::Test {};
using MatchFractionGraphs = ::testing::Types<DBGSuccinct, DBGSSHash>;
TYPED_TEST_SUITE(MatchFractionAlignment, MatchFractionGraphs);

std::pair<std::string, std::string> sparse_seed_sequences() {
    std::mt19937 random(73);
    std::string reference(73, 'A');
    for (char &c : reference)
        c = "ACGT"[random() % 4];
    std::string sequence = reference;
    for (size_t i : { 25, 35, 45, 55, 65 })
        sequence[i] = sequence[i] == 'A' ? 'C' : 'A';
    return {reference, sequence};
}

TYPED_TEST(MatchFractionAlignment, UsesCompletedAlignment) {
    const auto [reference, sequence] = sparse_seed_sequences();
    auto anno_graph = build_anno_graph<TypeParam, annot::ColumnCompressed<>>(
        11, { reference }, { "reference" });
    const auto &graph = anno_graph->get_graph();
    Query query(graph, sequence);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.min_seed_length = 10;
    config.max_seed_length = 11;
    config.two_sided = true;
    config.forward_and_reverse_complement = false;

    std::vector<bool> covered(sequence.size());
    for (const auto &anchor : ExactSeeder(query, config).get_anchors())
        std::fill(covered.begin() + anchor.get_clipping(),
                  covered.begin() + anchor.get_clipping() + anchor.get_seed().size(), true);
    ASSERT_LT(std::count(covered.begin(), covered.end(), true), 0.7 * sequence.size());

    for (bool labeled : { false, true }) {
        for (double fraction : { 0.7, 1.0 }) {
            SCOPED_TRACE(labeled);
            SCOPED_TRACE(fraction);
            config.min_exact_match = fraction;
            std::vector<Alignment> results;
            auto collect = [&](Alignment&& aln) { results.push_back(std::move(aln)); };
            AnnotationBuffer buffer(graph, anno_graph->get_annotator());
            if (labeled) {
                LabeledSeeder seeder(buffer, query, config);
                LabeledExtender extender(buffer, query, config);
                align_query(query, seeder, extender, collect);
            } else {
                ExactSeeder seeder(query, config);
                Extender extender(query, config);
                align_query(query, seeder, extender, collect);
            }
            if (fraction == 1.0) {
                EXPECT_TRUE(results.empty());
            } else {
                ASSERT_FALSE(results.empty());
                EXPECT_EQ(reference, results.front().get_spelling());
                EXPECT_EQ(reference.size() - 5, results.front().get_cigar().get_num_matches());
            }
        }
    }
}

TYPED_TEST(MatchFractionAlignment, LowCoverageLabelSurvivesAnExactOtherLabel) {
    const auto [reference, sequence] = sparse_seed_sequences();
    auto anno_graph = build_anno_graph<TypeParam, annot::ColumnCompressed<>>(
        11, {reference, sequence}, {"reference", "exact"});
    const auto &graph = anno_graph->get_graph();
    const auto &encoder = anno_graph->get_annotator().get_label_encoder();
    Query query(graph, sequence);
    DBGAlignerConfig config;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);
    config.min_seed_length = 10;
    config.max_seed_length = 11;
    config.min_exact_match = 0.7;
    config.forward_and_reverse_complement = false;
    config.two_sided = true;
    AnnotationBuffer buffer(graph, anno_graph->get_annotator());
    LabeledSeeder seeder(buffer, query, config);
    LabeledExtender extender(buffer, query, config);
    std::set<std::string> found;
    align_query(query, seeder, extender, [&](Alignment&& alignment) {
        for (auto column : buffer.get_cached_column_set(alignment.get_label_classes()[0])) {
            auto label = encoder.decode(column);
            found.insert(label);
            EXPECT_EQ(label == "reference" ? reference : sequence, alignment.get_spelling());
        }
    });
    EXPECT_EQ((std::set<std::string>{"reference", "exact"}), found);
}




#if ! _PROTEIN_GRAPH
TEST(AnnotationBufferRedone, PendingReverseComplementAliases) {
    auto anno_graph = build_anno_graph<DBGSSHash, annot::ColumnCompressed<>>(
        3, { "AACGT" }, { "reference" }, DeBruijnGraph::CANONICAL);
    const auto &graph = dynamic_cast<const DBGSSHash&>(anno_graph->get_graph());
    // Queue every base node before its reverse complement, in both insertion
    // orders and with a fetch either between orientations or after both.
    for (bool reverse_first : { false, true }) {
        for (bool fetch_between : { false, true }) {
            AnnotationBuffer buffer(graph, anno_graph->get_annotator());
            std::vector<DeBruijnGraph::node_index> nodes;
            graph.call_nodes([&](auto node) {
                if (graph.get_base_node(node) == node)
                    nodes.push_back(node);
            });
            ASSERT_FALSE(nodes.empty());
            for (auto node : nodes)
                buffer.queue_path({ reverse_first ? graph.reverse_complement(node) : node });
            if (fetch_between)
                buffer.fetch_queued_annotations();
            for (auto node : nodes)
                buffer.queue_path({ reverse_first ? node : graph.reverse_complement(node) });
            buffer.fetch_queued_annotations();
            for (auto node : nodes) {
                const auto *forward = buffer.get_labels(node);
                const auto *reverse = buffer.get_labels(graph.reverse_complement(node));
                ASSERT_NE(nullptr, forward);
                ASSERT_NE(nullptr, reverse);
                ASSERT_FALSE(forward->empty());
                EXPECT_EQ(*forward, *reverse);
            }
        }
    }
}
#endif

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

    config.chain_alignments = true;
    config.two_sided = true;
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

        // These fixtures specify the best alignment for each label. Extending
        // isolated seeds can also produce lower-scoring paths for that label.
        for (auto *candidates : { &paths, &paths_no_extend, &paths_no_align }) {
            // Check every result, including alternatives discarded below:
            // a valid CIGAR alone cannot detect a traceback into another label.
            for (const auto &path : *candidates) {
                for (auto node : path.get_path()) {
                    if (node == DeBruijnGraph::npos) continue; // unconnected anchor gap
                    // Sequence lookup maps canonical aliases to annotation rows.
                    const auto labels = anno_graph.get_labels(graph.get_node_sequence(node), 1.0);
                    for (auto column : anno_buffer.get_cached_column_set(path.get_label_classes()[0])) {
                        const auto &label = label_encoder.decode(column);
                        EXPECT_NE(labels.end(), std::find(labels.begin(), labels.end(), label))
                            << label << " at node " << node << " in " << path;
                    }
                }
            }
            tsl::hopscotch_set<Anchor::label_class_t> labels;
            std::erase_if(*candidates, [&](const auto &path) {
                return !labels.insert(path.get_label_classes()[0]).second;
            });
        }

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
                    std::vector<std::string> decoded_labels;
                    for (auto column : anno_buffer.get_cached_column_set(path_targets[0]))
                        decoded_labels.push_back(label_encoder.decode(column));
                    const std::string path_target_labels = fmt::format("{}", fmt::join(decoded_labels, ";"));

                    if (path_targets.size() > 1) {
                        ASSERT_TRUE(std::all_of(path_targets.begin(), path_targets.end(),
                                                [&](auto node_target) { return node_target == path_targets[0]; }))
                            << path << "\t" << path_target_labels << "\t" << fmt::format("{}", fmt::join(path_targets, ","));
                    }
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

TYPED_TEST(LabeledAlignerRedoneTest, SharedSeedKeepsAllLabels) {
    const std::string sequence = "AACGT";
    const std::vector<std::string> labels { "A", "B" };

    for (bool coordinates : { false, true }) {
        if (coordinates
                && !std::is_same_v<typename TypeParam::second_type, annot::ColumnCompressed<>>)
            continue;
        SCOPED_TRACE(coordinates);

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            sequence.size(), { sequence, sequence }, labels, DeBruijnGraph::BASIC, coordinates);
        AnnotationBuffer buffer(anno_graph->get_graph(), anno_graph->get_annotator());
        Query query(anno_graph->get_graph(), sequence);
        DBGAlignerConfig config{};
        config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -2);
        config.min_seed_length = sequence.size();
        config.max_seed_length = sequence.size();
        LabeledSeeder seeder(buffer, query, config);

        // The shared seed introduces {A, B}; splitting it into singleton label
        // sets grows the cache while get_anchors() is iterating its labels.
        auto anchors = seeder.get_anchors();
        ASSERT_EQ(labels.size(), anchors.size());
        std::vector<std::string> actual_labels;
        for (const auto &anchor : anchors) {
            const auto &columns = buffer.get_cached_column_set(anchor.get_label_class());
            ASSERT_EQ(1u, columns.size());
            actual_labels.push_back(anno_graph->get_annotator().get_label_encoder().decode(columns[0]));
            EXPECT_EQ(sequence, anchor.get_seed());
            EXPECT_EQ(0u, anchor.get_clipping());
            EXPECT_FALSE(anchor.get_orientation());
            EXPECT_EQ(coordinates ? 0 : Anchor::ncoord, anchor.get_coord());
        }
        std::sort(actual_labels.begin(), actual_labels.end());
        EXPECT_EQ(labels, actual_labels);

        // Public emission merges the two identical paths, retaining both
        // labels and their coordinates despite the one-alternative limit.
        LabeledExtender extender(buffer, query, config);
        std::vector<Alignment> results;
        align_query(query, seeder, extender, [&](Alignment&& aln) {
            results.push_back(std::move(aln));
        });
        ASSERT_EQ(1u, results.size());
        const auto &columns = buffer.get_cached_column_set(results[0].get_label_classes()[0]);
        EXPECT_EQ(2u, columns.size());
        if (coordinates) {
            ASSERT_EQ(2u, results[0].get_coords().size());
            EXPECT_EQ(AnnotationBuffer::Tuple({ 0 }), results[0].get_coords()[0]);
            EXPECT_EQ(AnnotationBuffer::Tuple({ 0 }), results[0].get_coords()[1]);
        }
    }
}

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

TYPED_TEST(LabeledAlignerRedoneTest, RelativeCutoffKeepsLabelsSeparate) {
    auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                       typename TypeParam::second_type>(
        4, { "GCAAT", "AATGCTT" }, { "A", "B" });
    AnnotationBuffer buffer(anno_graph->get_graph(), anno_graph->get_annotator());
    Query query(anno_graph->get_graph(), "GCAATGCTT");
    DBGAlignerConfig config{};
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -1, -1);
    config.min_seed_length = config.max_seed_length = 4;
    config.forward_and_reverse_complement = false;
    config.rel_score_cutoff = 0.95;
    LabeledSeeder seeder(buffer, query, config);
    LabeledExtender extender(buffer, query, config);
    std::map<std::string, std::string> actual;
    align_query(query, seeder, extender, [&](Alignment&& aln) {
        for (auto column : buffer.get_cached_column_set(aln.get_label_classes()[0])) {
            auto label = anno_graph->get_annotator().get_label_encoder().decode(column);
            actual.emplace(label, aln.get_cigar().to_string());
        }
    });
    EXPECT_EQ((std::map<std::string, std::string>{{ "A", "5=4S" }, { "B", "2S7=" }}), actual);
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
                    using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
                    auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(std::get<1>(mappings[0]).data(), k);
                    auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                    // we can only find the first seed of it's also the minimizer
                    if (sshash::util::compute_minimizer<kmer_t>(first_kmer_ref, k, dict.m(), dict.hasher()).pos_in_seq
                            || sshash::util::compute_minimizer<kmer_t>(first_kmer_query, k, dict.m(), dict.hasher()).pos_in_seq) {
                        first_seed_matches = false;
                    }
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
                        using kmer_t = typename std::decay_t<decltype(dict)>::kmer_type;
                        auto first_kmer_ref = sshash::util::string_to_uint_kmer<kmer_t>(std::get<1>(mappings[0]).data(), k);
                        auto first_kmer_query = sshash::util::string_to_uint_kmer<kmer_t>(query.data(), k);
                        // we can only find the first seed of it's also the minimizer
                        if (sshash::util::compute_minimizer<kmer_t>(first_kmer_ref, k, dict.m(), dict.hasher()).pos_in_seq
                                || sshash::util::compute_minimizer<kmer_t>(first_kmer_query, k, dict.m(), dict.hasher()).pos_in_seq) {
                            first_seed_matches = false;
                        }
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

TEST(LabeledChaining, ConnectionUsesBestScoreAcrossReferenceLengths) {
    const std::string sequence = "ACGTTAGGACTTTCAGAC";
    const std::string longer = "ACGTTAGTATTTTCAGAC";
    const std::string shorter = "ACGTTAGGTTTCAGAC";
    auto annotated = build_anno_graph<DBGHashFast, annot::ColumnCompressed<>>(
        7, {longer, shorter}, {"sample", "sample"}, DeBruijnGraph::BASIC, true, {0, 100});
    const auto &graph = annotated->get_graph();
    AnnotationBuffer buffer(graph, annotated->get_annotator());
    Query query(graph, sequence);
    DBGAlignerConfig config{};
    config.min_seed_length = 7;
    config.max_seed_length = std::numeric_limits<size_t>::max();
    config.forward_and_reverse_complement = false;
    config.seed_complexity_filter = false;
    config.chain_alignments = true;
    config.score_matrix = DBGAlignerConfig::dna_scoring_matrix(2, -3, -3);

    const auto results = LabeledSeeder(buffer, query, config).get_inexact_anchors();
    ASSERT_FALSE(results.empty());
    const auto best = std::max_element(results.begin(), results.end(), [&](const auto &a, const auto &b) {
        return score_match(a, config) < score_match(b, config);
    });
    // The shorter path reaches the target first (cost 18, score 25), but two
    // mismatches on the longer path score 26 despite their greater cost 20.
    EXPECT_EQ(26, score_match(*best, config));
    EXPECT_EQ("7=1X1=1X8=", best->get_cigar().to_string());
    EXPECT_EQ(longer, best->get_spelling());
    EXPECT_EQ(AnnotationBuffer::CoordinateSet({{0}}), best->get_coords());
    for (const auto &alignment : results) {
        EXPECT_TRUE(alignment.is_spelling_valid(graph));
        EXPECT_TRUE(alignment.get_cigar().is_valid(alignment.get_spelling(), alignment.get_query()));
    }
}

} // namespace
