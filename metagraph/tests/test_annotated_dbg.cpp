#include <stdio.h>

#include "gtest/gtest.h"

#include "annotated_dbg.hpp"
#include "annotate_column_compressed.hpp"


TEST(AnnotatedDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        auto *graph = new DBG_succ(k);
        AnnotatedDBG anno_graph(
            graph,
            new annotate::ColumnCompressed<>(graph->num_edges() + 1)
        );

        ASSERT_EQ(anno_graph.get_graph().num_edges() + 1,
                  anno_graph.get_annotation().num_objects());
        bit_vector_dyn inserted_edges(anno_graph.get_graph().num_edges() + 1, 0);

        std::string sequence(100, 'A');
        anno_graph.get_graph().add_sequence(sequence, true, &inserted_edges);

        ASSERT_EQ(k + 1, anno_graph.get_graph().num_nodes());
        ASSERT_EQ(k + 2, anno_graph.get_graph().num_edges());

        std::vector<uint64_t> inserted_edge_idx;
        for (uint64_t j = 1; j <= inserted_edges.num_set_bits(); ++j) {
            inserted_edge_idx.push_back(inserted_edges.select1(j));
        }

        anno_graph.get_annotation().insert_rows(inserted_edge_idx);
        EXPECT_EQ(anno_graph.get_graph().num_edges() + 1, inserted_edges.size());

        anno_graph.annotate_sequence(sequence, { "Label" });

        EXPECT_EQ(std::vector<std::string> { "Label" },
                  anno_graph.get_labels(sequence, 1));
    }
}
