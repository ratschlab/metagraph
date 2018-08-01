#include <stdio.h>

#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "annotate_color_compressed.hpp"
#include "config.hpp"


TEST(ColoredDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k);
        annotate::ColorCompressed<> annotation(graph.num_edges() + 1);
        bit_vector_dyn inserted_edges(graph.num_edges() + 1, 0);

        std::string sequence(100, 'A');
        graph.add_sequence(sequence, true, &inserted_edges);

        ASSERT_EQ(k + 1, graph.num_nodes());
        ASSERT_EQ(k + 2, graph.num_edges());

        std::vector<uint64_t> inserted_edge_idx;
        for (uint64_t j = 1; j <= inserted_edges.num_set_bits(); ++j) {
            inserted_edge_idx.push_back(inserted_edges.select1(j));
        }

        annotation.insert_rows(inserted_edge_idx);

        graph.map_to_nodes(sequence, [&](uint64_t i) {
            if (i > 0)
                annotation.add_color(i, "Label");
        });

        EXPECT_EQ(1u, annotation.num_colors());
    }
}
