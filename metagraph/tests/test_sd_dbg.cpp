#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "dbg_sd.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


TEST(DBGSD, InitializeEmpty) {
    {
        DBGSD graph(20, false);
        ASSERT_EQ(1u, graph.kmer_to_node("AAAAAAAAAAAAAAAAAAAA"));
        ASSERT_EQ(graph.num_nodes(), graph.kmer_to_node("TTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGSD graph(20, true);
        ASSERT_EQ(1u, graph.kmer_to_node("AAAAAAAAAAAAAAAAAAAA"));
        ASSERT_EQ(1u, graph.kmer_to_node("TTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGSD, SerializeEmpty) {
    {
        DBGSD graph(20, false);

        graph.serialize(test_dump_basename);
    }

    {
        DBGSD graph(0, false);

        ASSERT_TRUE(graph.load(test_dump_basename));
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGSD graph(20, true);

        graph.serialize(test_dump_basename);
    }

    {
        DBGSD graph(0, true);

        ASSERT_TRUE(graph.load(test_dump_basename));
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGSD, InsertSequence) {
    {
        DBGSD graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGSD graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGSD, ReverseComplement) {
    {
        DBGSD graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        //uint64_t num_nodes = graph.num_nodes();
        graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");


        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGSD graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        //uint64_t num_nodes = graph.num_nodes();
        graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGSD, CheckGraph) {
    {
        DBGSD graph(20, false);

        const std::string alphabet = "ACGT";

        for (size_t i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            graph.add_sequence(seq);
        }

        for (uint64_t i = 1; i <= 1000; ++i) {
            auto kmer = graph.node_to_kmer(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = graph.kmer_to_node(kmer);
            ASSERT_EQ(i, node);
        }
    }

    {
        DBGSD graph(20, true);

        const std::string alphabet = "ACGT";

        for (size_t i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            graph.add_sequence(seq);
        }

        for (uint64_t i = 1; i <= 1000; ++i) {
            auto kmer = graph.node_to_kmer(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = graph.kmer_to_node(kmer);
            ASSERT_EQ(i, node);
        }
    }
}

TEST(DBGSD, Serialize) {
    {
        DBGSD graph(20, false);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGSD graph(0, false);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGSD graph(20, true);

        graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGSD graph(0, true);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}
