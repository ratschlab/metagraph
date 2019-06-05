#include <stdio.h>

#include "gtest/gtest.h"

#define protected public
#define private public

#include <set>

#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "annotated_dbg.hpp"
#include "annotate_column_compressed.hpp"


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

        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));

        check_labels(anno_graph, seq_second, { "Second" }, { "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'C'), { "First" }, { "Second", "Third" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(100, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
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
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));

        check_labels(anno_graph, std::string(100, 'A'), { "First", "Third" }, { "Second" });
        check_labels(anno_graph, std::string(100, 'T'), { "Second" }, { "First", "Third" });
        check_labels(anno_graph, std::string(k, 'A') + std::string(k, 'C'), { "First" }, { "Second", "Third" });

        check_labels(anno_graph, std::string(100, 'G'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Third" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Second" });

        check_labels(anno_graph, std::string(100, 'N'),
                     k == 1 ? std::vector<std::string>{ "Second", "Third" } : std::vector<std::string>{ "Second" },
                     k == 1 ? std::vector<std::string>{ "First" } : std::vector<std::string>{ "First", "Third" });
    }
}
