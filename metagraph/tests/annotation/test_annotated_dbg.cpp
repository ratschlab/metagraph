#include <stdio.h>

#include <set>
#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "test_annotated_dbg_helpers.hpp"

#include "common/threads/threading.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "common/vectors/vector_algorithm.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

// this next #include includes AnnotatedDBG. we need access to its protected
// members to modify the underlying annotator
#define protected public
#include "annotation/annotation_converters.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;
using namespace mtg::graph;
using namespace mtg::annot;

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

std::vector<uint64_t> edge_to_row_idx(const bitmap &edge_mask) {
    // transform indexes of k-mers to the annotation format
    std::vector<uint64_t> rows;
    edge_mask.call_ones([&](auto i) {
        rows.push_back(AnnotatedDBG::graph_to_anno_index(i));
    });
    return rows;
}

TEST(AnnotatedDBG, ExtendGraphWithSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        auto graph = std::make_shared<DBGSuccinct>(k + 1);
        AnnotatedDBG anno_graph(graph, std::make_unique<ColumnCompressed<>>(1));

        ASSERT_EQ(anno_graph.get_graph().num_nodes(),
                  anno_graph.get_annotation().num_objects());

        std::string sequence(100, 'A');

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(sequence,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        ASSERT_EQ(k + 2, anno_graph.get_graph().num_nodes());
        EXPECT_EQ(1u, anno_graph.get_annotation().num_objects());

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        EXPECT_FALSE(anno_graph.label_exists("Label"));
        EXPECT_FALSE(anno_graph.label_exists("NotLabel"));

        anno_graph.annotate_sequence(std::string(sequence), { "Label" });

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
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));

        anno_graph.annotate_sequence(std::string(seq_first), { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));

        check_labels(anno_graph, seq_first, { "First" }, { "Second", "Third" });

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        check_labels(anno_graph, seq_first, { "First" }, { "Second", "Third" });

        anno_graph.annotate_sequence(std::string(seq_second), { "Second" });

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
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );
        EXPECT_EQ(num_nodes, anno_graph->get_graph().num_nodes());

        EXPECT_FALSE(anno_graph->label_exists("First"));
        EXPECT_FALSE(anno_graph->label_exists("Second"));
        EXPECT_FALSE(anno_graph->label_exists("Third"));

        anno_graph->annotate_sequence(std::string(seq_first), { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_FALSE(anno_graph->label_exists("Second"));
        EXPECT_FALSE(anno_graph->label_exists("Third"));

        check_labels(*anno_graph, seq_first, { "First" }, { "Second", "Third" });

        bit_vector_dyn inserted_nodes(anno_graph->get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph->annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph->get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph->get_labels(seq_first, 1));

        check_labels(*anno_graph, seq_first, { "First" }, { "Second", "Third" });

        anno_graph->annotate_sequence(std::string(seq_second), { "Second" });

        anno_graph = std::make_unique<AnnotatedDBG>(
            graph,
            std::unique_ptr<AnnotatedDBG::Annotator>(
                convert<RowFlatAnnotator>(
                    dynamic_cast<ColumnCompressed<>&&>(*anno_graph->annotator_)
                ).release()
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
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_first), { "First" });

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(std::string(seq_second), { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_third), { "Third" });
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

        ThreadPool thread_pool(10);
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_first), std::vector<std::string> { "First" }
        );
        thread_pool.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_second), std::vector<std::string> { "Second" }
        );
        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_third), std::vector<std::string> { "Third" }
        );
        thread_pool.join();

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
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );
        EXPECT_EQ(num_nodes, anno_graph.get_graph().num_nodes());

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_first), { "First" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(std::string(seq_second), { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_third), { "Third" });

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

        ThreadPool thread_pool(10);
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + k
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_first), std::vector<std::string> { "First" }
        );
        thread_pool.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_second), std::vector<std::string> { "Second" }
        );
        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_third), std::vector<std::string> { "Third" }
        );
        thread_pool.join();

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

        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_first), { "First" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        anno_graph.annotate_sequence(std::string(seq_second), { "Second" });

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_TRUE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        anno_graph.annotate_sequence(std::string(seq_third), { "Third" });

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

        ThreadPool thread_pool(10);
        AnnotatedDBG anno_graph(
            graph,
            std::make_unique<ColumnCompressed<>>(graph->max_index())
        );

        EXPECT_FALSE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_first), std::vector<std::string> { "First" }
        );
        thread_pool.join();

        EXPECT_TRUE(anno_graph.label_exists("First"));
        EXPECT_FALSE(anno_graph.label_exists("Second"));
        EXPECT_FALSE(anno_graph.label_exists("Third"));
        EXPECT_FALSE(anno_graph.label_exists("Fourth"));

        EXPECT_TRUE(anno_graph.get_annotation().num_objects() + 1
                        < dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss().num_edges())
            << dynamic_cast<const DBGSuccinct&>(anno_graph.get_graph()).get_boss();

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        bit_vector_dyn inserted_nodes(anno_graph.get_graph().max_index() + 1, 0);
        graph->add_sequence(seq_second,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );
        graph->add_sequence(seq_third,
            [&](auto new_node) { inserted_nodes.insert_bit(new_node, true); }
        );

        anno_graph.annotator_->insert_rows(edge_to_row_idx(inserted_nodes));
        EXPECT_EQ(anno_graph.get_graph().num_nodes() + 1, inserted_nodes.size());

        ASSERT_EQ(std::vector<std::string> { "First" },
                  anno_graph.get_labels(seq_first, 1));

        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_second), std::vector<std::string> { "Second" }
        );
        thread_pool.enqueue(
            [&anno_graph](const std::string &sequence, const auto &labels) {
                anno_graph.annotate_sequence(sequence, labels);
            },
            std::string(seq_third), std::vector<std::string> { "Third" }
        );
        thread_pool.join();

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
typedef ::testing::Types<std::pair<DBGBitmap, ColumnCompressed<>>,
                         std::pair<DBGHashString, ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, ColumnCompressed<>>,
                         std::pair<DBGHashFast, ColumnCompressed<>>,
                         std::pair<DBGSuccinct, ColumnCompressed<>>,
                         std::pair<DBGBitmap, RowFlatAnnotator>,
                         std::pair<DBGHashString, RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, RowFlatAnnotator>,
                         std::pair<DBGHashFast, RowFlatAnnotator>,
                         std::pair<DBGSuccinct, RowFlatAnnotator>
                        > GraphAnnotationPairTypes;
TYPED_TEST_SUITE(AnnotatedDBGTest, GraphAnnotationPairTypes);


template <typename GraphAnnotationPair>
class AnnotatedDBGWithNTest : public ::testing::Test {};
typedef ::testing::Types<std::pair<DBGHashString, ColumnCompressed<>>,
                         std::pair<DBGSuccinct, ColumnCompressed<>>,
                         std::pair<DBGHashString, RowFlatAnnotator>,
                         std::pair<DBGSuccinct, RowFlatAnnotator>
                        > GraphWithNAnnotationPairTypes;
TYPED_TEST_SUITE(AnnotatedDBGWithNTest, GraphWithNAnnotationPairTypes);


#if ! _PROTEIN_GRAPH
template <typename GraphAnnotationPair>
class AnnotatedDBGNoNTest : public ::testing::Test {};
typedef ::testing::Types<std::pair<DBGBitmap, ColumnCompressed<>>,
                         std::pair<DBGHashOrdered, ColumnCompressed<>>,
                         std::pair<DBGHashFast, ColumnCompressed<>>,
                         std::pair<DBGBitmap, RowFlatAnnotator>,
                         std::pair<DBGHashOrdered, RowFlatAnnotator>,
                         std::pair<DBGHashFast, RowFlatAnnotator>
                        > GraphNoNAnnotationPairTypes;
TYPED_TEST_SUITE(AnnotatedDBGNoNTest, GraphNoNAnnotationPairTypes);
#endif



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

#if ! _PROTEIN_GRAPH
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
#endif

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

        EXPECT_EQ(std::vector<std::string>({ "First", "Third" }),
                  anno_graph->get_labels(seq_first, 0));
#if _DNA_GRAPH
        EXPECT_EQ(k == 1 ? std::vector<std::string>({ "Second", "Third" })
                         : std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 0));
        EXPECT_EQ(k == 1 ? std::vector<std::string>({ "First", "Second", "Third" })
                         : std::vector<std::string>({ "First", "Third" }),
                  anno_graph->get_labels(seq_third, 0));
#else
        EXPECT_EQ(k <= 3 ? std::vector<std::string>({ "Second", "Third" })
                         : std::vector<std::string> { "Second" },
                  anno_graph->get_labels(seq_second, 0));
        EXPECT_EQ(k <= 3 ? std::vector<std::string>({ "First", "Second", "Third" })
                         : std::vector<std::string>({ "First", "Third" }),
                  anno_graph->get_labels(seq_third, 0));
#endif
    }
}

TYPED_TEST(AnnotatedDBGWithNTest, get_top_label_signatures) {
    typedef std::vector<std::pair<std::string, sdsl::bit_vector>> VectorSignature;

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

        const auto &label_encoder = anno_graph->get_annotation().get_label_encoder();
        auto comp = [&](const std::pair<std::string, sdsl::bit_vector> &a,
                        const std::pair<std::string, sdsl::bit_vector> &b) {
            size_t a_cnt = sdsl::util::cnt_one_bits(a.second);
            size_t b_cnt = sdsl::util::cnt_one_bits(b.second);
            return a_cnt > b_cnt
                || (a_cnt == b_cnt
                        && label_encoder.encode(a.first) < label_encoder.encode(b.first));
        };

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        std::vector<std::vector<bool>> temps(6);
        temps[0].assign(100 - (k + 1) + 1, true);
        temps[0].insert(temps[0].end(), 100 + k - (k + 1) + 1 - (100 - (k + 1) + 1), false);

        temps[1].assign(102 - (k + 1) + 1, true);
        temps[1].insert(temps[1].end(), 202 - (k + 1) + 1 - (102 - (k + 1) + 1), false);

        temps[2].assign(100 - (k + 1) + 1, true);
        temps[2].insert(temps[2].end(), 202 - (k + 1) + 1 - 2 * (100 - (k + 1) + 1), false);
        temps[2].insert(temps[2].end(), 100 - (k + 1) + 1, true);

        temps[3].assign(202 - (k + 1) + 1 - (100 - (k + 1) + 1), false);
        temps[3].insert(temps[3].end(), 100 - (k + 1) + 1, true);

        std::vector<VectorSignature> results {
            {
                std::make_pair("First", sdsl::bit_vector((100 + k) - (k + 1) + 1, true)),
                std::make_pair("Third", to_sdsl(temps[0]))
            },
            {
#if _DNA_GRAPH
                std::make_pair("Second", to_sdsl(temps[1]))
#else
                std::make_pair("Second", sdsl::bit_vector(202 - (k + 1) + 1, true))
#endif
            },
            {
#if _DNA_GRAPH
                std::make_pair("Third", to_sdsl(temps[2])),
                std::make_pair("First", to_sdsl(temps[3]))
#else
                std::make_pair("Third", sdsl::bit_vector(202 - (k + 1) + 1, true)),
                std::make_pair("First", to_sdsl(temps[3]))
#endif
            }
        };
        if (k == 1) {
#if _DNA_GRAPH
            temps[4].assign(100 - (k + 1) + 1 + 1, false);
            temps[4].insert(temps[4].end(), 1, true);
            temps[4].insert(temps[4].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1 + 1) - 1, false);

            temps[5].assign(100 - (k + 1) + 1, true);
            temps[5].insert(temps[5].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1), false);

            results[1].emplace_back("Third", to_sdsl(temps[4]));
            results[2].emplace_back("Second", to_sdsl(temps[5]));
#else
            temps[4].assign(100 - (k + 1) + 1 + 1, false);
            temps[4].insert(temps[4].end(), 102 - (k + 1) + 1, true);
            temps[4].insert(temps[4].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1 + 1) - (102 - (k + 1) + 1), false);

            temps[5].assign(102 - (k + 1) + 1, true);
            temps[5].insert(temps[5].end(), 202 - (k + 1) + 1 - (102 - (k + 1) + 1), false);

            results[1].emplace_back("Third", to_sdsl(temps[4]));
            results[2].emplace_back("Second", to_sdsl(temps[5]));
            std::swap(results[2][1], results[2][2]);
#endif
        }

#ifndef _DNA_GRAPH
        switch (k) {
            case 2:
                temps[4].assign(100 - (k + 1) + 1 + 2, false);
                temps[4].insert(temps[4].end(), 2, true);
                temps[4].insert(temps[4].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1 + 2) - 2, false);

                temps[5].assign(100 - (k + 1) + 1, false);
                temps[5].insert(temps[5].end(), 2, true);
                temps[5].insert(temps[5].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1) - 2, false);

                results[1].emplace_back("Third", to_sdsl(temps[4]));
                results[2].emplace_back("Second", to_sdsl(temps[5]));
                break;
            case 3:
                temps[4].assign(100 - (k + 1) + 1 + 3, false);
                temps[4].insert(temps[4].end(), 1, true);
                temps[4].insert(temps[4].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1 + 3) - 1, false);

                temps[5].assign(100 - (k + 1) + 2, false);
                temps[5].insert(temps[5].end(), 1, true);
                temps[5].insert(temps[5].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 2) - 1, false);

                results[1].emplace_back("Third", to_sdsl(temps[4]));
                results[2].emplace_back("Second", to_sdsl(temps[5]));
                break;
        }
#endif

        std::vector<double> percentages;
        for (size_t i = 0; i < results.size(); ++i) {
            percentages.clear();
            std::transform(results[i].begin(), results[i].end(),
                           std::back_inserter(percentages),
                           [&](const auto &pair) {
                               return 1. * sdsl::util::cnt_one_bits(pair.second)
                                   / (sequences[i].size() - (k + 1) + 1);
                           });
            percentages.emplace_back(0.0);

            for (size_t j = 1; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_label_signatures(sequences[i], j);
                std::sort(label_counts.begin(), label_counts.end(), comp);
                ASSERT_GE(j, label_counts.size());
#if _DNA_GRAPH
                if (k == 1 && i == 2 && j == 2) {
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second);
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + j),
                              label_counts) << k << " " << i << " " << j;
                }
#else
                EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + j),
                          label_counts) << k << " " << i << " " << j;
#endif

                for (size_t m = 1; m <= j; ++m) {
#ifdef _DNA_GRAPH
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;
#endif

                    auto label_counts = anno_graph->get_top_label_signatures(
                        sequences[i],
                        j,
                        percentages[m] + 1e-9
                    );
                    ASSERT_GE(j, label_counts.size());
                    std::sort(label_counts.begin(), label_counts.end(), comp);
                    EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m),
                              label_counts) << k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
#ifdef _DNA_GRAPH
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;
#endif

                auto label_counts = anno_graph->get_top_label_signatures(
                    sequences[i],
                    results[i].size() + 1,
                    percentages[m] + 1e-9
                );
                ASSERT_GE(results[i].size() + 1, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m),
                          label_counts)
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

#ifdef _DNA_GRAPH
        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_label_signatures(
                    sequences[i],
                    j,
                    percentages[m] + 1e-9
                );
                ASSERT_GE(j, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m - 1),
                          label_counts) << k << " " << i << " " << j << " " << m;
            }

            auto label_counts = anno_graph->get_top_label_signatures(
                sequences[i],
                results[i].size() + 1,
                percentages[m] + 1e-9
            );
            ASSERT_GE(results[i].size() + 1, label_counts.size());
            std::sort(label_counts.begin(), label_counts.end(), comp);
            EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m - 1),
                      label_counts)
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
#endif
    }
}

#if ! _PROTEIN_GRAPH
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

TYPED_TEST(AnnotatedDBGNoNTest, get_top_label_signatures) {
    typedef std::vector<std::pair<std::string, sdsl::bit_vector>> VectorSignature;

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

        const auto &label_encoder = anno_graph->get_annotation().get_label_encoder();
        auto comp = [&](const std::pair<std::string, sdsl::bit_vector> &a,
                        const std::pair<std::string, sdsl::bit_vector> &b) {
            size_t a_cnt = sdsl::util::cnt_one_bits(a.second);
            size_t b_cnt = sdsl::util::cnt_one_bits(b.second);
            return a_cnt > b_cnt
                || (a_cnt == b_cnt
                        && label_encoder.encode(a.first) < label_encoder.encode(b.first));
        };

        EXPECT_TRUE(anno_graph->label_exists("First"));
        EXPECT_TRUE(anno_graph->label_exists("Second"));
        EXPECT_TRUE(anno_graph->label_exists("Third"));
        EXPECT_FALSE(anno_graph->label_exists("Fourth"));

        std::vector<std::vector<bool>> temps(6);
        temps[0].assign(100 - (k + 1) + 1, true);
        temps[0].insert(temps[0].end(), 100 + k - (k + 1) + 1 - (100 - (k + 1) + 1), false);

        temps[1].assign(102 - (k + 1) + 1, true);
        temps[1].insert(temps[1].end(), 202 - (k + 1) + 1 - (102 - (k + 1) + 1), false);

        temps[2].assign(100 - (k + 1) + 1, true);
        temps[2].insert(temps[2].end(), 202 - (k + 1) + 1 - 2 * (100 - (k + 1) + 1), false);
        temps[2].insert(temps[2].end(), 100 - (k + 1) + 1, true);

        temps[3].assign(202 - (k + 1) + 1 - (100 - (k + 1) + 1), false);
        temps[3].insert(temps[3].end(), 100 - (k + 1) + 1, true);

        std::vector<VectorSignature> results {
            {
                std::make_pair("First", sdsl::bit_vector((100 + k) - (k + 1) + 1, true)),
                std::make_pair("Third", to_sdsl(temps[0]))
            },
            {
                std::make_pair("Second", to_sdsl(temps[1]))
            },
            {
                std::make_pair("Third", to_sdsl(temps[2])),
                std::make_pair("First", to_sdsl(temps[3]))
            }
        };

        if (k == 1) {
            temps[4].assign(100 - (k + 1) + 1 + 1, false);
            temps[4].insert(temps[4].end(), 1, true);
            temps[4].insert(temps[4].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1 + 1) - 1, false);

            temps[5].assign(100 - (k + 1) + 1, true);
            temps[5].insert(temps[5].end(), 202 - (k + 1) + 1 - (100 - (k + 1) + 1), false);

            results[1].emplace_back("Third", to_sdsl(temps[4]));
            results[2].emplace_back("Second", to_sdsl(temps[5]));
        }

        std::vector<double> percentages;
        for (size_t i = 0; i < results.size(); ++i) {
            percentages.clear();
            for (size_t j = 0; j < results[i].size(); ++j) {
                ASSERT_EQ(sequences[i].size() - (k + 1) + 1, results[i][j].second.size())
                    << k << " " << i << " " << j;
            }

            std::transform(results[i].begin(), results[i].end(),
                           std::back_inserter(percentages),
                           [&](const auto &pair) {
                               return 1. * sdsl::util::cnt_one_bits(pair.second)
                                   / (sequences[i].size() - (k + 1) + 1);
                           });
            percentages.emplace_back(0.0);

            for (size_t j = 1; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_label_signatures(sequences[i], j);
                std::sort(label_counts.begin(), label_counts.end(), comp);
                ASSERT_GE(j, label_counts.size());
                if (k == 1 && i == 2 && j == 2) {
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second)
                        << results[i][1].second << " " << label_counts[1].second;
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + j),
                              label_counts) << k << " " << i << " " << j;
                }

                for (size_t m = 1; m <= j; ++m) {
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;

                    auto label_counts = anno_graph->get_top_label_signatures(
                        sequences[i],
                        j,
                        percentages[m] + 1e-9
                    );
                    ASSERT_GE(j, label_counts.size());
                    std::sort(label_counts.begin(), label_counts.end(), comp);
                    EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m),
                              label_counts)<< k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;

                auto label_counts = anno_graph->get_top_label_signatures(
                    sequences[i],
                    results[i].size() + 1,
                    percentages[m] + 1e-9
                );
                ASSERT_GE(results[i].size() + 1, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m),
                          label_counts)
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_label_signatures(
                    sequences[i],
                    j,
                    percentages[m] + 1e-9
                );
                ASSERT_GE(j, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m - 1),
                          label_counts) << k << " " << i << " " << j << " " << m;
            }

            auto label_counts = anno_graph->get_top_label_signatures(
                sequences[i],
                results[i].size() + 1,
                percentages[m] + 1e-9
            );
            ASSERT_GE(results[i].size() + 1, label_counts.size());
            std::sort(label_counts.begin(), label_counts.end(), comp);
            EXPECT_EQ(VectorSignature(results[i].begin(), results[i].begin() + m - 1),
                      label_counts)
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
    }
}
#endif

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

        const auto &label_encoder = anno_graph->get_annotation().get_label_encoder();
        auto comp = [&](const std::pair<std::string, size_t> &a,
                        const std::pair<std::string, size_t> &b) {
            return a.second > b.second
                || (a.second == b.second
                    && label_encoder.encode(a.first) < label_encoder.encode(b.first));
        };

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
                auto label_counts = anno_graph->get_top_labels(sequences[i], j);
                std::sort(label_counts.begin(), label_counts.end(), comp);
#if _DNA_GRAPH
                if (k == 1 && i == 2 && j == 2) {
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second);
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                              label_counts)
                        << k << " " << i << " " << j;
                }
#else
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                          label_counts)
                    << k << " " << i << " " << j;
#endif

                for (size_t m = 1; m <= j; ++m) {
#ifdef _DNA_GRAPH
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;
#endif

                    auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                                   j,
                                                                   percentages[m] + 1e-9);
                    std::sort(label_counts.begin(), label_counts.end(), comp);
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                              label_counts)
                        << k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
#ifdef _DNA_GRAPH
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;
#endif

                auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                               results[i].size() + 1,
                                                               percentages[m] + 1e-9);
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                          label_counts)
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

#ifdef _DNA_GRAPH
        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                               j,
                                                               percentages[m] + 1e-9);
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                          label_counts)
                    << k << " " << i << " " << j << " " << m;
            }

            auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                           results[i].size() + 1,
                                                           percentages[m] + 1e-9);
            std::sort(label_counts.begin(), label_counts.end(), comp);
            EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                      label_counts)
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
#endif
    }
}

#if ! _PROTEIN_GRAPH
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

        const auto &label_encoder = anno_graph->get_annotation().get_label_encoder();
        auto comp = [&](const std::pair<std::string, size_t> &a,
                        const std::pair<std::string, size_t> &b) {
            return a.second > b.second
                || (a.second == b.second
                    && label_encoder.encode(a.first) < label_encoder.encode(b.first));
        };

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
                auto label_counts = anno_graph->get_top_labels(sequences[i], j);
                std::sort(label_counts.begin(), label_counts.end(), comp);
                ASSERT_GE(j, label_counts.size());
                if (k == 1 && i == 2 && j == 2) {
                    ASSERT_EQ(2u, label_counts.size());
                    EXPECT_EQ(results[i][0], label_counts[0]);
                    EXPECT_EQ(results[i][1].second, label_counts[1].second);
                    EXPECT_TRUE(results[i][1].first == "First"
                        || results[i][1].first == "Second");
                } else {
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + j),
                              label_counts)
                        << k << " " << i << " " << j;
                }

                for (size_t m = 1; m <= j; ++m) {
                    // Special case to handle later
                    if (k == 1 && i == 2 && m == 2)
                        continue;

                    auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                                   j,
                                                                   percentages[m] + 1e-9);
                    ASSERT_GE(j, label_counts.size());
                    std::sort(label_counts.begin(), label_counts.end(), comp);
                    EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                              label_counts)
                        << k << " " << i << " " << j << " " << m;
                }
            }

            for (size_t m = 1; m <= results[i].size(); ++m) {
                // Special case to handle later
                if (k == 1 && i == 2 && m == 2)
                    continue;

                auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                               results[i].size() + 1,
                                                               percentages[m] + 1e-9);
                ASSERT_GE(results[i].size() + 1, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m),
                          label_counts)
                    << k << " " << i << " " << results[i].size() + 1 << " " << m;
            }
        }

        if (k == 1) {
            // special case for third sequence (First and Second are equally good matches)
            size_t i = 2;
            size_t m = 2;
            for (size_t j = 2; j <= results[i].size(); ++j) {
                auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                               j,
                                                               percentages[m] + 1e-9);
                ASSERT_GE(j, label_counts.size());
                std::sort(label_counts.begin(), label_counts.end(), comp);
                EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                          label_counts)
                    << k << " " << i << " " << j << " " << m;
            }

            auto label_counts = anno_graph->get_top_labels(sequences[i],
                                                           results[i].size() + 1,
                                                           percentages[m] + 1e-9);
            ASSERT_GE(results[i].size() + 1, label_counts.size());
            std::sort(label_counts.begin(), label_counts.end(), comp);
            EXPECT_EQ(VectorCounts(results[i].begin(), results[i].begin() + m - 1),
                      label_counts)
                << k << " " << i << " " << results[i].size() + 1 << " " << m;
        }
    }
}
#endif

TEST(AnnotatedDBG, score_kmer_presence_mask) {
    auto anno_graph = build_anno_graph<DBGSuccinct, ColumnCompressed<>>(31);
    std::vector<std::pair<sdsl::bit_vector, int32_t>> results {
       { sdsl::bit_vector(), 0},
       { sdsl::bit_vector({
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }), 0 },
       { sdsl::bit_vector({
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }), 97 },
       { sdsl::bit_vector({
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1 }), 924 },
       { sdsl::bit_vector({
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }), 0 },
       { sdsl::bit_vector({
          0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0,
          1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
          1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1,
          0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1,
          0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1,
          0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0,
          1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0,
          1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0,
          0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0,
          1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
          0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0,
          1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1,
          1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0,
          0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1,
          1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
          1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
          0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0,
          0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1,
          0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1,
          1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0,
          0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0,
          0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0,
          0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
          1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1,
          0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1,
          0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1,
          1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
          1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
          1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1,
          0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1,
          0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1,
          0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0,
          1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0,
          0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1,
          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1 }), 869 },
       { sdsl::bit_vector({
          1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1,
          0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
          1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
          0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,
          0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0,
          0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0,
          0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
          0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
          0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1,
          0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0,
          1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0,
          1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
          1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
          1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
          0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1,
          1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1,
          1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0,
          1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
          1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1,
          1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
          0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0,
          1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
          1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0,
          0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0,
          0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0,
          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0,
          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1,
          1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0,
          1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0,
          0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
          1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
          1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1,
          0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
          1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1,
          0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
          1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1 }), 822 }
    };

    size_t i = 0;
    for (const auto &[vector, score] : results) {
        EXPECT_EQ(score, anno_graph->score_kmer_presence_mask(vector)) << i;
        ++i;
    }
}

} // namespace
