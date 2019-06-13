#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"
#include "../test_helpers.hpp"


KSEQ_INIT(gzFile, gzread);

#define protected public
#define private public

#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";

const KmerExtractor2Bit kmer_extractor;

uint64_t kmer_string_to_index(const DBGBitmap &graph, const std::string &kmer) {
    DBGBitmap::node_index node;
    graph.map_to_nodes_sequentially(
        kmer.begin(), kmer.end(),
        [&](const auto i) { node = i; }
    );
    return graph.node_to_index(node);
}


TEST(DBGBitmapComplete, InitializeComplete) {
    {
        DBGBitmap graph(20, false);
        ASSERT_TRUE(graph.is_complete());
        EXPECT_EQ(std::string("AAAAAAAAAAAAAAAAAAAA"), graph.get_node_sequence(1));
        EXPECT_EQ(1u, kmer_string_to_index(graph, "AAAAAAAAAAAAAAAAAAAA"));

    // #if _DNA4_GRAPH
        EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "TTTTTTTTTTTTTTTTTTTT"));
    // #elif _DNA_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "NNNNNNNNNNNNNNNNNNNN"));
    // #elif _DNA_CASE_SENSITIVE_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "tttttttttttttttttttt"));
    // #elif _PROTEIN_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "XXXXXXXXXXXXXXXXXXXX"));
    // #endif

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGBitmap graph(20, true);
        ASSERT_TRUE(graph.is_complete());
        EXPECT_EQ(std::string("AAAAAAAAAAAAAAAAAAAA"), graph.get_node_sequence(1));
        EXPECT_EQ(1u, kmer_string_to_index(graph, "AAAAAAAAAAAAAAAAAAAA"));

    // #if _DNA4_GRAPH
        EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "TTTTTTTTTTTTTTTTTTTT"));
    // #elif _DNA_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "NNNNNNNNNNNNNNNNNNNN"));
    // #elif _DNA_CASE_SENSITIVE_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "tttttttttttttttttttt"));
    // #elif _PROTEIN_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), kmer_string_to_index(graph, "XXXXXXXXXXXXXXXXXXXX"));
    // #endif

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGBitmapComplete, SerializeComplete) {
    {
        DBGBitmap graph(20, false);

        graph.serialize(test_dump_basename);
    }

    {
        DBGBitmap graph(2, false);

        ASSERT_TRUE(graph.load(test_dump_basename));
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGBitmap graph(20, true);

        graph.serialize(test_dump_basename);
    }

    {
        DBGBitmap graph(2, true);

        ASSERT_TRUE(graph.load(test_dump_basename));
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGBitmapComplete, InsertSequence) {
    {
        DBGBitmap graph(20, false);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);
        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGBitmap graph(20, true);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);
        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGBitmapComplete, ReverseComplement) {
    {
        DBGBitmap graph(20, false);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);

        //uint64_t num_nodes = graph.num_nodes();
        ASSERT_THROW(graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), std::runtime_error);

        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);


        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGBitmap graph(20, true);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);

        //uint64_t num_nodes = graph.num_nodes();
        ASSERT_THROW(graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), std::runtime_error);

        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGBitmapComplete, CheckGraph) {
    {
        DBGBitmap graph(20, false);

        const std::string alphabet = "ACGT";

        for (DBGBitmap::node_index i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            ASSERT_THROW(graph.add_sequence(seq), std::runtime_error);
        }

        for (DBGBitmap::node_index i = 1; i <= 1000; ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = kmer_string_to_index(graph, kmer);
            ASSERT_EQ(i, node);
        }
    }

    {
        DBGBitmap graph(20, true);

        const std::string alphabet = "ACGT";

        for (DBGBitmap::node_index i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            ASSERT_THROW(graph.add_sequence(seq), std::runtime_error);
        }

        for (DBGBitmap::node_index i = 1; i <= 1000; ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = kmer_string_to_index(graph, kmer);
            ASSERT_EQ(i, node) << kmer;
        }
    }
}

TEST(DBGBitmapComplete, Serialize) {
    {
        DBGBitmap graph(20, false);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);
        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGBitmap graph(2, false);

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
        DBGBitmap graph(20, true);

        ASSERT_THROW(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), std::runtime_error);
        ASSERT_THROW(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), std::runtime_error);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph.serialize(test_dump_basename);
    }

    {
        DBGBitmap graph(2, true);

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

/*
// TODO
TEST(DBGBitmap, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    DBGBitmapConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_sequences({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBGBitmap *graph = new DBGBitmap(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_nodes() + 1);

    EXPECT_EQ(1u, graph->outgoing(1, DBGBitmap::kSentinelCode));

    for (size_t i = 1; i <= graph->num_edges(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBGBitmap::kSentinelCode) {
            uint64_t node_idx = graph->rank_last(i - 1) + 1;

            EXPECT_EQ(outgoing_edges[i],
                graph->select_last(graph->outgoing(node_idx, graph->get_W(i)))
            ) << "Edge index: " << i << "\n"
              << "Outgoing: "->outgoing(node_idx, graph->get_W(i)) << "\n"
              << *graph;

            EXPECT_EQ(node_idx,
                graph->incoming(graph->outgoing(node_idx, graph->get_W(i)),
                                graph->get_minus_k_value(i, graph->get_k() - 1).first)
            );
            for (int c = 0; c < graph->alph_size; ++c) {
                uint64_t prev_node = graph->incoming(node_idx, c);
                if (prev_node) {
                    EXPECT_EQ(node_idx,
                        graph->outgoing(prev_node, graph->get_node_last_value(i))
                    );
                }
            }
        }
    }

    delete graph;
}
*/
