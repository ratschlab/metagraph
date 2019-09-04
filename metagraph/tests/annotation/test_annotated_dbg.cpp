#include <stdio.h>

#include <set>
#include "gtest/gtest.h"

#include "test_annotated_dbg_helpers.hpp"

#include "annotation_converters.hpp"


void check_labels(const AnnotatedDBG &anno_graph,
                  const std::string &sequence,
                  const std::vector<std::string> labels_present,
                  const std::vector<std::string> labels_not_present) {
    std::set<SequenceGraph::node_index> indices;
    anno_graph.get_graph().map_to_nodes(
        sequence,
        [&](const auto &index) {
            ASSERT_NE(SequenceGraph::npos, index);
            indices.insert(index);
            EXPECT_EQ(labels_present.size(), anno_graph.get_labels(index).size());
            EXPECT_EQ(convert_to_set(labels_present),
                      convert_to_set(anno_graph.get_labels(index)));

            for (const auto &label : labels_present) {
                EXPECT_TRUE(anno_graph.has_label(index, label));
            }

            for (const auto &label : labels_not_present) {
                EXPECT_FALSE(anno_graph.has_label(index, label));
            }
        }
    );

    for (const auto &label : labels_present) {
        std::set<SequenceGraph::node_index> cur_indices;
        anno_graph.call_annotated_nodes(
            label,
            [&](const auto &index) {
                ASSERT_NE(SequenceGraph::npos, index);
                cur_indices.insert(index);
                EXPECT_TRUE(anno_graph.has_label(index, label));
            }
        );
        std::vector<SequenceGraph::node_index> diff;
        std::set_difference(indices.begin(), indices.end(),
                            cur_indices.begin(), cur_indices.end(),
                            diff.begin());
        EXPECT_EQ(0u, diff.size())
            << diff.front()
            << anno_graph.get_graph().get_node_sequence(diff.front());
    }
}


TEST(AnnotatedDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        AnnotatedDBG anno_graph(std::make_shared<DBGSuccinct>(k + 1),
                                std::make_unique<annotate::ColumnCompressed<>>(1));

        ASSERT_EQ(anno_graph.get_graph().num_nodes(),
                  anno_graph.get_annotation().num_objects());

        std::string sequence(100, 'A');

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);
        anno_graph.graph_->add_sequence(sequence, &inserted_edges);

        ASSERT_EQ(k + 2, anno_graph.get_graph().num_nodes());
        EXPECT_EQ(1u, anno_graph.get_annotation().num_objects());

        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_edges.size());

        EXPECT_FALSE(anno_graph.label_exists("Label"));
        EXPECT_FALSE(anno_graph.label_exists("NotLabel"));

        anno_graph.annotate_sequence(sequence, { "Label" });

        EXPECT_EQ(std::vector<std::string> { "Label" },
                  anno_graph.get_labels(sequence, 1));

        EXPECT_TRUE(anno_graph.label_exists("Label"));
        EXPECT_FALSE(anno_graph.label_exists("NotLabel"));

        check_labels(anno_graph, sequence, { "Label" }, { "NotLabel" });
    }
}

TEST(AnnotatedDBG, ExtendGraphAddPath) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(2, 'N');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));

        check_labels(anno_graph, seq_first, { "First" }, { "Second", "Third" });

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        check_labels(anno_graph, seq_first, { "First" }, { "Second", "Third" });

        anno_graph.annotate_sequence(seq_second, { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));

#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 2)
                                                            / (seq_second.size() - (k + 1) + 1)));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));

        check_labels(anno_graph, seq_second, { "Second" }, { "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, Transform) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(2, 'N');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        uint64_t num_nodes = graph->num_nodes();
        auto anno_graph = std::make_unique<AnnotatedDBG>(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );
        EXPECT_EQ(num_nodes, anno_graph->get_graph().num_nodes());

        EXPECT_FALSE(anno_graph->label_exists("First"));
        EXPECT_FALSE(anno_graph->label_exists("Second"));
        EXPECT_FALSE(anno_graph->label_exists("Third"));

        anno_graph->annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_FALSE(anno_graph->label_exists("Second"));
        EXPECT_FALSE(anno_graph->label_exists("Third"));

        check_labels(*anno_graph, seq_first, { "First" }, { "Second", "Third" });

        bit_vector_dyn inserted_edges(anno_graph->get_graph().num_nodes() + 1, 0);

        anno_graph->graph_->add_sequence(seq_second, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph->annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));

        check_labels(*anno_graph, seq_first, { "First" }, { "Second", "Third" });

        anno_graph->annotate_sequence(seq_second, { "Second" });

        anno_graph = std::make_unique<AnnotatedDBG>(
            std::move(anno_graph->graph_),
            std::unique_ptr<AnnotatedDBG::Annotator>(
                annotate::convert<annotate::RowFlatAnnotator>(
                    std::move(dynamic_cast<annotate::ColumnCompressed<>&>(
                        *anno_graph->annotator_
                    )
                )).release()
            )
        );

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_FALSE(anno_graph->label_exists("Third"));

#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 2)
                                                            / (seq_second.size() - (k + 1) + 1)));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 1));

        check_labels(*anno_graph, seq_second, { "Second" }, { "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_third, { "Third" });
        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPathsParallel) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPathsWithoutDummy) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(0, false);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_third, { "Third" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPathsWithoutDummyParallel) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(10, false);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes),
            10
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPathsPruneDummy) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(k, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(0, true);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes)
        );

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_third, { "Third" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TEST(AnnotatedDBG, ExtendGraphAddTwoPathsPruneDummyParallel) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(k, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(100, 'N');

        std::string seq_third = std::string(100, 'G')
                                    + std::string(2, 'N')
                                    + std::string(100, 'A');

        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(10, true);

        uint64_t num_nodes = graph->num_nodes();
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<annotate::ColumnCompressed<>>(num_nodes),
            10
        );

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        anno_graph.graph_->add_sequence(seq_third, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_TRUE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
#endif

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

#ifndef _DNA_GRAPH
        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}


template <typename GraphAnnotationPair>
class AnnotatedDBGTest : public ::testing::Test {};

template <typename GraphAnnotationPair>
class AnnotatedDBGWithNTest : public ::testing::Test {};

template <typename GraphAnnotationPair>
class AnnotatedDBGNoNTest : public ::testing::Test {};

TYPED_TEST_CASE(AnnotatedDBGTest, GraphAnnotationPairTypes);
TYPED_TEST_CASE(AnnotatedDBGWithNTest, GraphWithNAnnotationPairTypes);
TYPED_TEST_CASE(AnnotatedDBGNoNTest, GraphNoNAnnotationPairTypes);

TYPED_TEST(AnnotatedDBGWithNTest, check_labels) {
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'),
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')
        };

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        check_labels(*anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(*anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(*anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(*anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });
#ifndef _DNA_GRAPH
        check_labels(*anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
#endif
    }
}

TYPED_TEST(AnnotatedDBGNoNTest, check_labels) {
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'),
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')
        };

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        check_labels(*anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(*anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(*anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(*anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });
    }
}

TYPED_TEST(AnnotatedDBGWithNTest, get_labels) {
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'),
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')
        };
        const auto& seq_first = sequences[0];
        const auto& seq_second = sequences[1];
        const auto& seq_third = sequences[2];

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));
#if _DNA_GRAPH
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph->get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
#else
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph->get_labels(seq_third, 1));
#endif
    }
}

TYPED_TEST(AnnotatedDBGNoNTest, get_labels) {
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'),
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')
        };
        const auto& seq_first = sequences[0];
        const auto& seq_second = sequences[1];
        const auto& seq_third = sequences[2];

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));

        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 1. * (seq_second.size() - (k + 1) + 1 - 100)
                                                            / (seq_second.size() - (k + 1) + 1) - 1e-9));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_third, 1));
        EXPECT_EQ(std::vector<std::string> {},
                  anno_graph->get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) + 1e-9));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph->get_labels(seq_third, 1. * (seq_third.size() - (k + 1) + 1 - (k + 2))
                                                            / (seq_third.size() - (k + 1) + 1) - 1e-9));
    }
}

TYPED_TEST(AnnotatedDBGWithNTest, get_top_labels) {
    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'), //TTTTTGGNNNNN
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')  //GGGGGNNAAAAA
        };

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        std::vector<VectorCounts> results {
            {
                std::make_pair("First", (100 + k) - (k + 1) + 1),
                std::make_pair("Third", 100 - (k + 1) + 1)
            },
            {
#if _DNA_GRAPH
                std::make_pair("Second", 102 - (k + 1) + 1)
#else
                std::make_pair("Second", 202 - (k + 1) + 1),
#endif
            },
            {
#if _DNA_GRAPH
                std::make_pair("Third", 2 * (100 - (k + 1) + 1)),
                std::make_pair("First", 100 - (k + 1) + 1)
#else
                std::make_pair("Third", 202 - (k + 1) + 1),
                std::make_pair("First", 100 - (k + 1) + 1),
#endif
            }
        };
        if (k == 1) {
#if _DNA_GRAPH
            results[1].emplace_back("Third", 1);
            results[2].emplace_back("Second", 100 - (k + 1) + 1);
#else
            results[1].emplace_back("Third", 2 + 100 - (k + 1) + 1);
            results[2].emplace_back("Second", 2 + 100 - (k + 1) + 1);
            std::swap(results[2][1], results[2][2]);
#endif
        }

#ifndef _DNA_GRAPH
        switch (k) {
            case 2:
                results[1].emplace_back("Third", 2);
                results[2].emplace_back("Second", 2);
                break;
            case 3:
                results[1].emplace_back("Third", 1);
                results[2].emplace_back("Second", 1);
                break;
        }
#endif

        std::vector<double> percentages;
        for (size_t i = 0; i < results.size(); ++i) {
            percentages.clear();
            std::transform(results[i].begin(), results[i].end(),
                           std::back_inserter(percentages),
                           [&](const auto &pair) {
                               return 1. * pair.second / (sequences[i].size() - (k + 1) + 1);
                           });
            percentages.emplace_back(0.0);

            for (size_t j = 1; j <= results[i].size(); ++j) {
#if _DNA_GRAPH
                if (k == 1 && i == 2 && j == 2) {
                    auto label_counts = anno_graph->get_top_labels(sequences[j], j);
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second);
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                              anno_graph->get_top_labels(sequences[i], j))
                        << k << " " << i << " " << j;
                }
#else
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                          anno_graph->get_top_labels(sequences[i], j))
                    << k << " " << i << " " << j;
#endif

                for (size_t m = 1; m <= j; ++m) {
#ifdef _DNA_GRAPH
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;
#endif

                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                              anno_graph->get_top_labels(sequences[i],
                                                         j,
                                                         percentages[m] + 1e-9))
                        << k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
#ifdef _DNA_GRAPH
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;
#endif

                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                          anno_graph->get_top_labels(sequences[i],
                                                     results[i].size() + 1,
                                                     percentages[m] + 1e-9))
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

#ifdef _DNA_GRAPH
        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                          anno_graph->get_top_labels(sequences[i],
                                                     j,
                                                     percentages[m] + 1e-9))
                    << k << " " << i << " " << j << " " << m;
            }

            EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                      anno_graph->get_top_labels(sequences[i],
                                                 results[i].size() + 1,
                                                 percentages[m] + 1e-9))
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
#endif
    }
}

TYPED_TEST(AnnotatedDBGNoNTest, get_top_labels) {
    typedef std::vector<std::pair<std::string, size_t>> VectorCounts;
    for (size_t k = 1; k < 10; ++k) {
        const std::vector<std::string> sequences {
            std::string(100, 'A') + std::string(k, 'C'),
            std::string(100, 'T') + std::string(2, 'G') + std::string(100, 'N'), //TTTTTGGNNNNN
            std::string(100, 'G') + std::string(2, 'N') + std::string(100, 'A')  //GGGGGNNAAAAA
        };

        auto anno_graph = build_anno_graph<typename TypeParam::first_type,
                                           typename TypeParam::second_type>(
            k + 1, sequences, { "First", "Second" , "Third" }
        );

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        std::vector<VectorCounts> results {
            {
                std::make_pair("First", (100 + k) - (k + 1) + 1),
                std::make_pair("Third", 100 - (k + 1) + 1)
            },
            {
                std::make_pair("Second", 102 - (k + 1) + 1)
            },
            {
                std::make_pair("Third", 2 * (100 - (k + 1) + 1)),
                std::make_pair("First", 100 - (k + 1) + 1)
            }
        };
        if (k == 1) {
            results[1].emplace_back("Third", 1);
            results[2].emplace_back("Second", 100 - (k + 1) + 1);
        }

        std::vector<double> percentages;
        for (size_t i = 0; i < results.size(); ++i) {
            percentages.clear();
            std::transform(results[i].begin(), results[i].end(),
                           std::back_inserter(percentages),
                           [&](const auto &pair) {
                               return 1. * pair.second / (sequences[i].size() - (k + 1) + 1);
                           });
            percentages.emplace_back(0.0);

            for (size_t j = 1; j <= results[i].size(); ++j) {
                if (k == 1 && i == 2 && j == 2) {
                    auto label_counts = anno_graph->get_top_labels(sequences[j], j);
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second);
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                              anno_graph->get_top_labels(sequences[i], j))
                        << k << " " << i << " " << j;
                }

                for (size_t m = 1; m <= j; ++m) {
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;

                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                              anno_graph->get_top_labels(sequences[i],
                                                         j,
                                                         percentages[m] + 1e-9))
                        << k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;

                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                          anno_graph->get_top_labels(sequences[i],
                                                     results[i].size() + 1,
                                                     percentages[m] + 1e-9))
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                          anno_graph->get_top_labels(sequences[i],
                                                     j,
                                                     percentages[m] + 1e-9))
                    << k << " " << i << " " << j << " " << m;
            }

            EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                      anno_graph->get_top_labels(sequences[i],
                                                 results[i].size() + 1,
                                                 percentages[m] + 1e-9))
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
    }
}
