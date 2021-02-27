#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include <gtest/gtest.h>

#include "../test_helpers.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";

const kmer::KmerExtractor2Bit kmer_extractor;


TEST(DBGBitmapComplete, InitializeComplete) {
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        DBGBitmap graph(12, mode);
        ASSERT_TRUE(graph.is_complete());
        EXPECT_EQ(std::string(12, 'A'), graph.get_node_sequence(1));
        EXPECT_EQ(1u, graph.kmer_to_node(std::string(12, 'A')));

        #if _DNA_GRAPH
            EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node(std::string(12, 'T')));
        // #elif _DNA5_GRAPH
        //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node(std::string(12, 'N')));
        // #elif _DNA_CASE_SENSITIVE_GRAPH
        //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node(std::string(12, 't')));
        // #elif _PROTEIN_GRAPH
        //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node(std::string(12, 'X')));
        #endif

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGBitmapComplete, SerializeComplete) {
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        {
            DBGBitmap graph(12, mode);

            graph.serialize(test_dump_basename);
        }

        {
            DBGBitmap graph(2);

            ASSERT_TRUE(graph.load(test_dump_basename));
            EXPECT_EQ(12u, graph.get_k());

            EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
            EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
            EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
            EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        }
    }
}

TEST(DBGBitmapComplete, InsertSequence) {
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        DBGBitmap graph(12, mode);

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
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        DBGBitmap graph(12, mode);

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
    for (DeBruijnGraph::Mode mode : { DeBruijnGraph::BASIC, DeBruijnGraph::CANONICAL }) {
        DBGBitmap graph(12, mode);

        const std::string alphabet = "ACGT";

        for (DBGBitmap::node_index i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % alphabet.size()];
            }
            ASSERT_THROW(graph.add_sequence(seq), std::runtime_error);
        }
        for (DBGBitmap::node_index i = 1; i <= graph.alphabet().size(); ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(12u, kmer.size());
            EXPECT_EQ(i, graph.kmer_to_node(kmer)) << kmer;
        }
// complete BitmapDBG uses contiguous indexing only for DNA4
#if _DNA_GRAPH
        for (DBGBitmap::node_index i = 1; i <= 1000; ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(12u, kmer.size());
            EXPECT_EQ(i, graph.kmer_to_node(kmer)) << kmer;
        }
#endif
    }
}

TEST(DBGBitmapComplete, Serialize) {
    {
        DBGBitmap graph(12, DeBruijnGraph::BASIC);

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
        DBGBitmap graph(2, DeBruijnGraph::BASIC);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(12u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGBitmap graph(12, DeBruijnGraph::CANONICAL);

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
        DBGBitmap graph(2, DeBruijnGraph::CANONICAL);

        ASSERT_TRUE(graph.load(test_dump_basename));

        EXPECT_EQ(12u, graph.get_k());

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

} // namespace
