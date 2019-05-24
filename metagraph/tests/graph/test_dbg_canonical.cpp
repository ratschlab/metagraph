#include "gtest/gtest.h"

#include "test_dbg_helpers.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


template <typename Graph>
class DeBruijnGraphCanonicalTest : public DeBruijnGraphTest<Graph> { };
typedef ::testing::Types<DBGBitmap, DBGHashOrdered, DBGSuccinct> CanonicalGraphTypes;
TYPED_TEST_CASE(DeBruijnGraphCanonicalTest, CanonicalGraphTypes);

template <typename Graph>
class DeBruijnGraphWithNTest : public DeBruijnGraphTest<Graph> { };
typedef ::testing::Types<DBGBitmap, DBGHashOrdered> NoNGraphTypes;
TYPED_TEST_CASE(DeBruijnGraphWithNTest, NoNGraphTypes);


TYPED_TEST(DeBruijnGraphCanonicalTest, CheckGraph) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGT", true));
}

TYPED_TEST(DeBruijnGraphWithNTest, CheckGraph) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGTN", false));
}

TYPED_TEST(DeBruijnGraphWithNTest, CheckGraphCanonical) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGTN", true));
}


TYPED_TEST(DeBruijnGraphCanonicalTest, InitializeEmpty) {
    auto graph = build_graph<TypeParam>(2, {}, true);

    EXPECT_EQ(0u, graph->num_nodes());
    EXPECT_FALSE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, SerializeEmpty) {
    {
        auto graph = build_graph<TypeParam>(20, {}, true);
        ASSERT_EQ(0u, graph->num_nodes());
        graph->serialize(test_dump_basename);
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(20u, graph.get_k());
    EXPECT_TRUE(graph.is_canonical_mode());

    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, Serialize) {
    {
        auto graph = build_graph<TypeParam>(20, {
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "CATGTACTAGCTGATCGTAGCTAGCTAGC"
        }, true);

        EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph->serialize(test_dump_basename);
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, InsertSequence) {
    auto graph = build_graph<TypeParam>(20, {
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "CATGTACTAGCTGATCGTAGCTAGCTAGC"
    }, true);

    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_TRUE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, ReverseComplement) {
    auto graph1 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" }, true);
    auto graph2 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                               "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT" }, true);

    EXPECT_EQ(graph1->num_nodes(), graph2->num_nodes());

    auto graph = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                              "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                                              "CATGTACTAGCTGATCGTAGCTAGCTAGC" }, true);
    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_TRUE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, Traversals1) {
    for (size_t k = 2; k <= 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C')
        }, true);

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        auto it = DeBruijnGraph::npos;
        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        map_to_nodes_sequentially(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );

        auto it2 = DeBruijnGraph::npos;
        map_to_nodes_sequentially(
            std::string(k - 1, 'A') + "C",
            [&](auto i) { it2 = i; }
        );
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it2, 'G'));

        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );

        map_to_nodes_sequentially(
            std::string(k - 1, 'G') + "T",
            [&](auto i) { it2 = i; }
        );
        ASSERT_NE(DeBruijnGraph::npos, it2);
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'A'));
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        EXPECT_EQ(it2, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'G'));
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, Traversals2) {
    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C'),
            std::string(100, 'G') + std::string(100, 'T')
        }, true);

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        auto it = DeBruijnGraph::npos;

        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));

        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));

        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));

        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it + 1, 'G'));

        // reverse complement
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'G'));
        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'T'), 'G'));
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, Traversals3) {
    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C'),
            std::string(100, 'G') + std::string(100, 'T')
        }, true);

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        auto it = DeBruijnGraph::npos;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        map_to_nodes_sequentially(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );

        DeBruijnGraph::node_index it2;
        graph->map_to_nodes(
            std::string(k - 1, 'A') + "C",
            [&](auto i) { it2 = i; }
        );
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it2, 'G'));

        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );

        graph->map_to_nodes(
            std::string(k - 1, 'G') + "T",
            [&](auto i) { it2 = i; }
        );
        ASSERT_NE(DeBruijnGraph::npos, it2);
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it2, 'T'));
        EXPECT_NE(DeBruijnGraph::npos, graph->traverse(it2, 'C'));

        map_to_nodes_sequentially(
            std::string(k - 1, 'G') + "T",
            [&](auto i) { it2 = i; }
        );
        ASSERT_NE(DeBruijnGraph::npos, it2);
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it2, 'A'));
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        EXPECT_EQ(it2, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'G'));
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, Traversals4) {
    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C')
        }, true);

        auto it = DeBruijnGraph::npos;

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it + 1, 'G'));


        // reverse complement
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'T'), 'G'));


        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'T'));
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, map_to_nodes_canonical) {
    for (size_t k = 2; k <= 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C'),
            std::string(100, 'G') + std::string(100, 'T')
        });
        auto graph_can = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C'),
            std::string(100, 'G') + std::string(100, 'T')
        }, true);

        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 2, 'A')
                                        + std::string(2 * (k - 1), 'C');
        auto rev_seq = sequence_to_map;
        reverse_complement(rev_seq.begin(), rev_seq.end());

        std::vector<DeBruijnGraph::node_index> indices_forward;
        graph->map_to_nodes(
            sequence_to_map,
            [&](auto i) { indices_forward.push_back(i); }
        );
        std::vector<DeBruijnGraph::node_index> indices_reverse;
        graph->map_to_nodes(
            rev_seq,
            [&](auto i) { indices_reverse.push_back(i); }
        );
        std::reverse(indices_reverse.begin(), indices_reverse.end());

        std::vector<DeBruijnGraph::node_index> indices_canonical;
        graph_can->map_to_nodes(
            sequence_to_map,
            [&](auto i) { indices_canonical.push_back(i); }
        );

        ASSERT_EQ(indices_forward.size(), indices_reverse.size());
        ASSERT_EQ(indices_forward.size(), indices_canonical.size());

        for (size_t i = 0; i < indices_forward.size(); ++i) {
            if (!indices_forward[i]) {
                EXPECT_EQ(indices_canonical[i], indices_reverse[i]);
            } else if (!indices_reverse[i]) {
                EXPECT_EQ(indices_canonical[i], indices_forward[i]);
            } else {
                EXPECT_EQ(indices_canonical[i], std::min(indices_forward[i], indices_reverse[i]));
            }
        }
    }
}
