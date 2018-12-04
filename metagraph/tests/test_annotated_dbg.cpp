#include <stdio.h>

#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "annotated_dbg.hpp"
#include "annotate_column_compressed.hpp"


TEST(AnnotatedDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        AnnotatedDBG anno_graph(new DBGSuccinct(k + 1),
                                new annotate::ColumnCompressed<>(1));

        ASSERT_EQ(anno_graph.get_graph().num_nodes(),
                  anno_graph.get_annotation().num_objects());

        std::string sequence(100, 'A');

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);
        anno_graph.graph_->add_sequence(sequence, &inserted_edges);

        ASSERT_EQ(k + 2, anno_graph.get_graph().num_nodes());

        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_edges.size());

        anno_graph.annotate_sequence(sequence, { "Label" });

        EXPECT_EQ(std::vector<std::string> { "Label" },
                  anno_graph.get_labels(sequence, 1));
    }
}

TEST(AnnotatedDBG, ExtendGraphAddPath) {
    for (size_t k = 1; k < 10; ++k) {

        std::string seq_first = std::string(100, 'A')
                                    + std::string(100, 'C');

        std::string seq_second = std::string(100, 'T')
                                    + std::string(2, 'G')
                                    + std::string(2, 'N');

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_nodes() + 1, 0);

        anno_graph.graph_->add_sequence(seq_second, &inserted_edges);
        AnnotatedDBG::insert_zero_rows(anno_graph.annotator_.get(), inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });

        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

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

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

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

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(0, false);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes())
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        anno_graph.annotate_sequence(seq_first, { "First" });

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

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(10, false);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes()),
            10
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

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

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(0, true);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

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

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
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

        auto graph = std::make_unique<DBGSuccinct>(k + 1);
        graph->add_sequence(seq_first);
        graph->mask_dummy_kmers(10, true);

        AnnotatedDBG anno_graph(
            graph.release(),
            new annotate::ColumnCompressed<>(graph->num_nodes()),
            10
        );

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

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

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
    }
}
