#include <stdio.h>

#include "gtest/gtest.h"

#include "annotated_dbg.hpp"
#include "annotate_column_compressed.hpp"


TEST(AnnotatedDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        AnnotatedDBG anno_graph(new DBG_succ(k),
                                new annotate::ColumnCompressed<>(1));

        ASSERT_EQ(anno_graph.get_graph().num_edges(),
                  anno_graph.get_annotation().num_objects());

        std::string sequence(100, 'A');

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);
        anno_graph.get_graph().add_sequence(sequence, true, &inserted_edges);

        ASSERT_EQ(k + 1, anno_graph.get_graph().num_nodes());
        ASSERT_EQ(k + 2, anno_graph.get_graph().num_edges());

        anno_graph.adjust_annotation(inserted_edges);
        EXPECT_EQ(anno_graph.get_graph().num_edges() + 1, inserted_edges.size());

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

        AnnotatedDBG anno_graph(new DBG_succ(k));
        anno_graph.get_graph().add_sequence(seq_first, false);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.get_graph().num_edges())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

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

        AnnotatedDBG anno_graph(new DBG_succ(k));
        anno_graph.get_graph().add_sequence(seq_first, false);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.get_graph().num_edges())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

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

        AnnotatedDBG anno_graph(new DBG_succ(k), 10);
        anno_graph.get_graph().add_sequence(seq_first, false);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.get_graph().num_edges())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

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

        AnnotatedDBG anno_graph(new DBG_succ(k));
        anno_graph.get_graph().add_sequence(seq_first, false);
        anno_graph.initialize_annotation_mask(0, false);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.num_anno_rows())
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

        anno_graph.annotate_sequence(seq_first, { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

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

        AnnotatedDBG anno_graph(new DBG_succ(k), 10);
        anno_graph.get_graph().add_sequence(seq_first, false);
        anno_graph.initialize_annotation_mask(10, false);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.num_anno_rows())
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

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

        AnnotatedDBG anno_graph(new DBG_succ(k));
        anno_graph.get_graph().add_sequence(seq_first, false);
        anno_graph.initialize_annotation_mask(0, true);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.num_anno_rows())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

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

        AnnotatedDBG anno_graph(new DBG_succ(k), 10);
        anno_graph.get_graph().add_sequence(seq_first, false);
        anno_graph.initialize_annotation_mask(10, true);

        anno_graph.set_annotation(
            new annotate::ColumnCompressed<>(anno_graph.num_anno_rows())
        );

        anno_graph.annotate_sequence(seq_first, { "First" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        anno_graph.get_graph().add_sequence(seq_second, false, &inserted_edges);
        anno_graph.get_graph().add_sequence(seq_third, false, &inserted_edges);
        anno_graph.adjust_annotation(inserted_edges);

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(seq_second, { "Second" });
        anno_graph.annotate_sequence(seq_third, { "Third" });
        anno_graph.join();

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < anno_graph.get_graph().num_edges()) << anno_graph.get_graph();

        EXPECT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));
        EXPECT_EQ(std::vector<std::string> { "Second" },
                  anno_graph.get_labels(seq_second, 1));
        EXPECT_EQ(std::vector<std::string> { "Third" },
                  anno_graph.get_labels(seq_third, 1));
    }
}
