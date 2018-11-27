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

    {
        DBGHashOrdered graph(20, true);

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
