#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "dbg_hash_ordered.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


TEST(DBGHashOrdered, InitializeEmpty) {
    {
        DBGHashOrdered graph(20, false);

        EXPECT_EQ(0u, graph.num_nodes());
        EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGHashOrdered graph(20, true);

        EXPECT_EQ(0u, graph.num_nodes());
        EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGHashOrdered, SerializeEmpty) {
    {
        DBGHashOrdered graph(20, false);

        ASSERT_EQ(0u, graph.num_nodes());

        graph.serialize(test_dump_basename);
    }

    {
        DBGHashOrdered graph(0, false);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(0u, graph.num_nodes());
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGHashOrdered graph(20, true);

        ASSERT_EQ(0u, graph.num_nodes());

        graph.serialize(test_dump_basename);
    }

    {
        DBGHashOrdered graph(0, true);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(0u, graph.num_nodes());
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGHashOrdered, InsertSequence) {
    {
        DBGHashOrdered graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGHashOrdered graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGHashOrdered, ReverseComplement) {
    {
        DBGHashOrdered graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        uint64_t num_nodes = graph.num_nodes();
        graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        EXPECT_EQ(num_nodes * 2, graph.num_nodes());

        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");


        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGHashOrdered graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        uint64_t num_nodes = graph.num_nodes();
        graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        EXPECT_EQ(num_nodes, graph.num_nodes());

        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGHashOrdered, CheckGraph) {
    {
        DBGHashOrdered graph(20, false);

        const std::string alphabet = "ACGT";

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

    {
        DBGHashOrdered graph(20, true);

        const std::string alphabet = "ACGT";

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
}

TEST(DBGHashOrdered, Traversals) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashOrdered graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        graph.map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });

        uint64_t it2 = graph.kmer_to_node(std::string(k - 1, 'A') + "C");
        ASSERT_EQ(graph.kmer_to_node(std::string(k, 'A')), it);
        EXPECT_EQ(it, graph.traverse(it, 'A'));
        EXPECT_EQ(it2, graph.traverse(it, 'C'));
        EXPECT_EQ(it, graph.traverse_back(it2, 'A'));
        EXPECT_EQ(DBGHashOrdered::npos, graph.traverse(it, 'G'));
        EXPECT_EQ(DBGHashOrdered::npos, graph.traverse_back(it2, 'G'));
    }
}

TEST(DBGHashOrdered, OutgoingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashOrdered graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        std::vector<DBGHashOrdered::node_index> adjacent_nodes;

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

TEST(DBGHashOrdered, IncomingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGHashOrdered graph(k);

        graph.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                 + std::string(100, 'G'));

        uint64_t it = 0;
        std::vector<DBGHashOrdered::node_index> adjacent_nodes;

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

TEST(DBGHashOrdered, Serialize) {
    {
        DBGHashOrdered graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGHashOrdered graph(0, false);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGHashOrdered graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGHashOrdered graph(0, true);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGHashOrdered, get_outdegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        std::unique_ptr<DBGHashOrdered> graph { new DBGHashOrdered(k) };
        graph->add_sequence(std::string(k - 1, 'A') + 'C');
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->outdegree(1));
    }
}

TEST(DBGHashOrdered, get_maximum_outdegree) {
    for (size_t k = 2; k < 10; ++k) {
        std::unique_ptr<DBGHashOrdered> graph { new DBGHashOrdered(k) };
        graph->add_sequence(std::string(k - 1, 'A') + 'A');
        graph->add_sequence(std::string(k - 1, 'A') + 'C');
        graph->add_sequence(std::string(k - 1, 'A') + 'G');
        graph->add_sequence(std::string(k - 1, 'A') + 'T');

        DBGHashOrdered::node_index max_outdegree_node_index;
        graph->map_to_nodes(std::string(k, 'A'), [&](DBGHashOrdered::node_index node) {
                                                    max_outdegree_node_index = node; });

        EXPECT_EQ(4ull, graph->num_nodes());
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == max_outdegree_node_index)
                EXPECT_EQ(4ull, graph->outdegree(i));
            else
                EXPECT_EQ(0ull, graph->outdegree(i));
        }
    }
}

TEST(DBGHashOrdered, get_outdegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        std::unique_ptr<DBGHashOrdered> graph { new DBGHashOrdered(k) };
        graph->add_sequence(std::string(k - 1, 'A') + std::string(k - 1, 'C') +
                            std::string(k - 1, 'G') + std::string(k, 'T'));
        graph->add_sequence(std::string(k, 'A'));

        DBGHashOrdered::node_index loop_node_index;
        graph->map_to_nodes(std::string(k, 'A'), [&](DBGHashOrdered::node_index node) { loop_node_index = node; });
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == loop_node_index)
                EXPECT_EQ(2ull, graph->outdegree(i));
            else
                EXPECT_EQ(1ull, graph->outdegree(i));
        }
    }
}
