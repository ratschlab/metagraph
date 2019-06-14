#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "dbg_hash_string.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


TEST(DBGHashString, InitializeEmpty) {
    DBGHashString graph(20);

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
}

TEST(DBGHashString, SerializeEmpty) {
    {
        DBGHashString graph(20);

        ASSERT_EQ(0u, graph.num_nodes());

        graph.serialize(test_dump_basename);
    }

    DBGHashString graph(0);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(20u, graph.get_k());

    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
}

TEST(DBGHashString, InsertSequence) {
    DBGHashString graph(20);

    graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TEST(DBGHashString, CheckGraph) {
    DBGHashString graph(20);

    const std::string alphabet = "ACGTN";

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        graph.add_sequence(seq);
    }

    for (uint64_t i = 1; i <= graph.num_nodes(); ++i) {
        ASSERT_EQ(i, graph.kmer_to_node(graph.node_to_kmer(i)));
    }
}

TEST(DBGHashString, Traversals) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashString graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        graph.map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });

        uint64_t it2 = graph.kmer_to_node(std::string(k - 1, 'A') + "C");
        ASSERT_EQ(graph.kmer_to_node(std::string(k, 'A')), it);
        EXPECT_EQ(it, graph.traverse(it, 'A'));
        EXPECT_EQ(it2, graph.traverse(it, 'C'));
        EXPECT_EQ(it, graph.traverse_back(it2, 'A'));
        EXPECT_EQ(DBGHashString::npos, graph.traverse(it, 'G'));
        EXPECT_EQ(DBGHashString::npos, graph.traverse_back(it2, 'G'));
    }
}

TEST(DBGHashString, OutgoingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashString graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        std::vector<DBGHashString::node_index> adjacent_nodes;

        // AA, AAAAA
        graph.map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph.adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ it, graph.traverse(it, 'C') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph.traverse(it, 'C');
        graph.adjacent_outgoing_nodes(it, &adjacent_nodes);
        auto outset = convert_to_set(std::vector<uint64_t>{ graph.traverse(it, 'C') });
        if (k == 2) {
            outset.insert(graph.traverse(it, 'G'));
            ASSERT_EQ(2u, adjacent_nodes.size());
        } else {
            ASSERT_EQ(1u, adjacent_nodes.size());
        }

        EXPECT_EQ(
            outset,
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CC, CCCCC
        graph.map_to_nodes(std::string(k, 'C'), [&](auto i) { it = i; });
        graph.adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph.traverse(it, 'G')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CG, CCCCG
        it = graph.traverse(it, 'G');
        graph.adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph.traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GG, GGGGG
        graph.map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        graph.adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph.traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TEST(DBGHashString, IncomingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashString graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        std::vector<DBGHashString::node_index> adjacent_nodes;

        // AA, AAAAA
        graph.map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph.adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ it }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph.traverse(it, 'C');
        graph.adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph.traverse_back(it, 'A') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CC, CCCCC
        graph.map_to_nodes(std::string(k, 'C'), [&](auto i) { it = i; });
        graph.adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph.traverse_back(it, 'A')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CG, CCCCG
        it = graph.traverse(it, 'C');
        graph.adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                graph.traverse_back(it, 'A'),
                graph.traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GG, GGGGG
        graph.map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        graph.adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph.traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TEST(DBGHashString, Serialize) {
    {
        DBGHashString graph(20);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));

        graph.serialize(test_dump_basename);
    }

    DBGHashString graph(0);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TEST(DBGHashString, get_outdegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);
        graph->add_sequence(std::string(k - 1, 'A') + 'C');
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->outdegree(1));
    }
}

TEST(DBGHashString, get_maximum_outdegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);
        graph->add_sequence(std::string(k - 1, 'A') + 'A');
        graph->add_sequence(std::string(k - 1, 'A') + 'C');
        graph->add_sequence(std::string(k - 1, 'A') + 'G');
        graph->add_sequence(std::string(k - 1, 'A') + 'T');

        auto max_outdegree_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_EQ(4ull, graph->num_nodes());
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == max_outdegree_node_index) {
                EXPECT_EQ(4ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(0ull, graph->outdegree(i));
            }
        }
    }
}

TEST(DBGHashString, get_outdegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);
        graph->add_sequence(std::string(k - 1, 'A') + std::string(k - 1, 'C') +
                            std::string(k - 1, 'G') + std::string(k, 'T'));
        graph->add_sequence(std::string(k, 'A'));

        auto loop_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_TRUE(graph->num_nodes() > 1);
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(1ull, graph->outdegree(i));
            }
        }
    }
}

TEST(DBGHashString, get_indegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);
        graph->add_sequence(std::string(k - 1, 'A') + 'C');
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->indegree(1));
    }
}

TEST(DBGHashString, get_maximum_indegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);
        graph->add_sequence('A' + std::string(k - 1, 'A'));
        graph->add_sequence('C' + std::string(k - 1, 'A'));
        graph->add_sequence('G' + std::string(k - 1, 'A'));
        graph->add_sequence('T' + std::string(k - 1, 'A'));

        auto max_indegree_node_index = graph->kmer_to_node(std::string(k, 'A'));

        ASSERT_EQ(4ull, graph->num_nodes());
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == max_indegree_node_index) {
                EXPECT_EQ(4ull, graph->indegree(i));
            } else {
                EXPECT_EQ(0ull, graph->indegree(i));
            }
        }
    }
}

TEST(DBGHashString, get_indegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGHashString>(k);

        graph->add_sequence(std::string(k, 'A')
                                + std::string(k - 1, 'C')
                                + std::string(k - 1, 'G')
                                + std::string(k, 'T'));

        auto loop_node_index = graph->kmer_to_node(std::string(k, 'T'));

        ASSERT_TRUE(graph->num_nodes() > 1);
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->indegree(i));
            } else {
                EXPECT_EQ(1ull, graph->indegree(i));
            }
        }
    }
}
