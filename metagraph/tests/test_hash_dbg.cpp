#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "dbg_hash.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


TEST(DBGHash, InitializeEmpty) {
    DBGHash graph(20);

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
}

TEST(DBGHash, SerializeEmpty) {
    {
        DBGHash graph(20);
    
        ASSERT_EQ(0u, graph.num_nodes());
    
        graph.serialize(test_dump_basename);
    }

    DBGHash graph(0);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(20u, graph.get_k());

    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
}

TEST(DBGHash, InsertSequence) {
    DBGHash graph(20);

    graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TEST(DBGHash, CheckGraph) {
    DBGHash graph(20);

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

TEST(DBGHash, Serialize) {
    {
        DBGHash graph(20);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));

        graph.serialize(test_dump_basename);
    }

    DBGHash graph(0);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}
