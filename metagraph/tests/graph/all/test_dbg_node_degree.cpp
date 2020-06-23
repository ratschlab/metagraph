#include "gtest/gtest.h"

#define private public
#define protected public

#include "../../test_helpers.hpp"
#include "test_dbg_helpers.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

TYPED_TEST_SUITE(DeBruijnGraphTest, GraphTypes);


TYPED_TEST(DeBruijnGraphTest, get_outdegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A') + 'C' });
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->outdegree(graph->kmer_to_node(std::string(k - 1, 'A') + 'C')));
    }
}

TYPED_TEST(DeBruijnGraphTest, get_maximum_outdegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k - 1, 'A') + 'A',
            std::string(k - 1, 'A') + 'C',
            std::string(k - 1, 'A') + 'G',
            std::string(k - 1, 'A') + 'T'
        });

        auto max_outdegree_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_EQ(4ull, graph->num_nodes());
        graph->call_nodes([&](auto i) {
            if (i == max_outdegree_node_index) {
                EXPECT_EQ(4ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(0ull, graph->outdegree(i));
            }
        });
    }
}

void check_degree_functions(const DeBruijnGraph &graph) {
    graph.call_nodes([&](auto i) {

        size_t outdegree = graph.outdegree(i);
        if (outdegree == 0) {
            ASSERT_FALSE(graph.has_single_outgoing(i));
            ASSERT_FALSE(graph.has_multiple_outgoing(i));
        } else if (outdegree == 1) {
            ASSERT_TRUE(graph.has_single_outgoing(i));
            ASSERT_FALSE(graph.has_multiple_outgoing(i));
        } else {
            ASSERT_FALSE(graph.has_single_outgoing(i));
            ASSERT_TRUE(graph.has_multiple_outgoing(i));
        }
    });

    graph.call_nodes([&](auto i) {

        size_t indegree = graph.indegree(i);

        if (indegree == 0) {
            ASSERT_TRUE(graph.has_no_incoming(i));
            ASSERT_FALSE(graph.has_single_incoming(i));
        } else if (indegree == 1) {
            ASSERT_FALSE(graph.has_no_incoming(i));
            ASSERT_TRUE(graph.has_single_incoming(i));
        } else {
            ASSERT_FALSE(graph.has_no_incoming(i));
            ASSERT_FALSE(graph.has_single_incoming(i));
        }

    });
}

TYPED_TEST(DeBruijnGraphTest, get_outdegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k - 1, 'A') + std::string(k - 1, 'C') +
                std::string(k - 1, 'G') + std::string(k, 'T'),
            std::string(k, 'A')
        });

        auto loop_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_TRUE(graph->num_nodes() > 1);
        graph->call_nodes([&](auto i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(1ull, graph->outdegree(i));
            }
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, get_indegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A') + 'C' });

        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->indegree(graph->kmer_to_node(std::string(k - 1, 'A') + 'C')));

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, get_maximum_indegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            'A' + std::string(k - 1, 'A'),
            'C' + std::string(k - 1, 'A'),
            'G' + std::string(k - 1, 'A'),
            'T' + std::string(k - 1, 'A')
        });

        auto max_indegree_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_EQ(4ull, graph->num_nodes());
        graph->call_nodes([&](auto i) {
            if (i == max_indegree_node_index) {
                EXPECT_EQ(4ull, graph->indegree(i));
            } else {
                EXPECT_EQ(0ull, graph->indegree(i));
            }
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, get_indegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k, 'T')
        });

        auto loop_node_index = graph->kmer_to_node(std::string(k, 'T'));

        ASSERT_TRUE(graph->num_nodes() > 1);
        graph->call_nodes([&](auto i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->indegree(i));
            } else {
                EXPECT_EQ(1ull, graph->indegree(i));
            }
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree1) {
    auto graph = build_graph<TypeParam>(3, { "ACTAAGCCC",
                                             "AAAGC",
                                             "TAAGCA" });
    ASSERT_EQ(9u, graph->num_nodes());

    EXPECT_EQ(2u, graph->indegree(graph->kmer_to_node("CCC")));
    EXPECT_EQ(1u, graph->outdegree(graph->kmer_to_node("CCC")));

    EXPECT_EQ(2u, graph->indegree(graph->kmer_to_node("AAA")));
    EXPECT_EQ(2u, graph->outdegree(graph->kmer_to_node("AAA")));

    check_degree_functions(*graph);
}

TYPED_TEST(DeBruijnGraphTest, get_degree1) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k, 'T')
        });

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node(std::string(k, 'T'));
        ASSERT_NE(DeBruijnGraph::npos, node_T);
        EXPECT_EQ(1ull, graph->outdegree(node_T));
        EXPECT_EQ(2ull, graph->indegree(node_T));

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, get_degree2) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k - 1, 'T')
        });

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node('G' + std::string(k - 1, 'T'));
        ASSERT_NE(DeBruijnGraph::npos, node_T);
        EXPECT_EQ(0ull, graph->outdegree(node_T));
        EXPECT_EQ(1ull, graph->indegree(node_T));

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_incoming_indegree) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        graph->call_nodes([&](auto node) {
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, [&](auto i) { incoming_nodes.push_back(i); });
            EXPECT_EQ(graph->indegree(node), incoming_nodes.size())
                << "adjacent_incoming_nodes and indegree are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_indegree_traverse_back) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        graph->call_nodes([&](auto node) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            EXPECT_EQ(graph->indegree(node), num_incoming_edges)
                << "traverse_back and indegree are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_traverse_back_incoming) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        graph->call_nodes([&](auto node) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, [&](auto i) { incoming_nodes.push_back(i); });
            EXPECT_EQ(num_incoming_edges, incoming_nodes.size())
                << "adjacent_incoming_nodes and traverse_back are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, is_single_outgoing_simple) {
    size_t k = 4;
    std::string reference = "CATC";

    auto graph = build_graph<TypeParam>(k, { reference });

    uint64_t single_outgoing_counter = 0;
    graph->call_nodes([&](auto i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    });

    EXPECT_EQ(0u, single_outgoing_counter);
}

TYPED_TEST(DeBruijnGraphTest, is_single_outgoing_for_multiple_valid_edges) {
    size_t k = 4;
    std::string reference = "AGGGGTC";

    auto graph = build_graph<TypeParam>(k, { reference });

    uint64_t single_outgoing_counter = 0;
    graph->call_nodes([&](auto i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    });

    EXPECT_EQ(1u, single_outgoing_counter);
}

} // namespace
