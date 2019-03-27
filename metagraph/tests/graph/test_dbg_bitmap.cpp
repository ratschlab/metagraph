#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef ASSERT_DEATH
#undef ASSERT_DEATH
#define ASSERT_DEATH(a, b) (void)0
#endif
#endif

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


TEST(DBGBitmapConstructedFull, InitializeComplete) {
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

TEST(DBGBitmapConstructedFull, SerializeComplete) {
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

TEST(DBGBitmapConstructedFull, InsertSequence) {
    {
        DBGBitmap graph(20, false);

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");
        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGBitmap graph(20, true);

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");
        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGBitmapConstructedFull, ReverseComplement) {
    {
        DBGBitmap graph(20, false);

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");

        //uint64_t num_nodes = graph.num_nodes();
        ASSERT_DEATH(graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), "");

        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");


        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }

    {
        DBGBitmap graph(20, true);

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");

        //uint64_t num_nodes = graph.num_nodes();
        ASSERT_DEATH(graph.add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), "");

        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_TRUE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_TRUE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    }
}

TEST(DBGBitmapConstructedFull, CheckGraph) {
    {
        DBGBitmap graph(20, false);

        const std::string alphabet = "ACGT";

        for (size_t i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            ASSERT_DEATH(graph.add_sequence(seq), "");
        }

        for (uint64_t i = 1; i <= 1000; ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = kmer_string_to_index(graph, kmer);
            ASSERT_EQ(i, node);
        }
    }

    {
        DBGBitmap graph(20, true);

        const std::string alphabet = "ACGT";

        for (size_t i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            ASSERT_DEATH(graph.add_sequence(seq), "");
        }

        for (uint64_t i = 1; i <= 1000; ++i) {
            auto kmer = graph.get_node_sequence(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = kmer_string_to_index(graph, kmer);
            ASSERT_EQ(i, node) << kmer;
        }
    }
}

TEST(DBGBitmapConstructedFull, Serialize) {
    {
        DBGBitmap graph(20, false);

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");
        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");

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

        ASSERT_DEATH(graph.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), "");
        ASSERT_DEATH(graph.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC"), "");

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

TEST(DBGBitmapConstructed, InsertSequence) {
    DBGBitmapConstructor constructor(20);
    constructor.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    constructor.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");
    DBGBitmap graph(20);
    constructor.build_graph(&graph);

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TEST(DBGBitmapConstructed, CheckGraph) {
    DBGBitmapConstructor constructor(20);

    const std::string alphabet = "ACGTN";

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        constructor.add_sequence(seq);
    }
    DBGBitmap graph(20);
    constructor.build_graph(&graph);

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        EXPECT_TRUE(graph.find(seq));
    }
}

TEST(DBGBitmapConstructed, Serialize) {
    {
        DBGBitmapConstructor constructor(20);
        constructor.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        constructor.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");
        DBGBitmap graph(20);
        constructor.build_graph(&graph);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));

        graph.serialize(test_dump_basename);
    }

    DBGBitmap graph(20);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

constexpr double kEps = std::numeric_limits<double>::epsilon();

//TODO
/*
void test_graph(DBGBitmap *graph, const bit_vector_sd &ref, bool canonical_only) {
    DBGBitmap ngraph(2);
    ngraph.k_ = ngraph.infer_k();
    EXPECT_EQ(ref, graph->kmers_);
    EXPECT_EQ(canonical_only, graph->canonical_only_);
}
*/


TEST(DBGBitmap, GraphDefaultConstructor) {
    DBGBitmap *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new DBGBitmap(2);
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        DBGBitmap graph(2);
    });
}

/*

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

TEST(DBGBitmap, AddSequenceSimplePath) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        DBGBitmap graph(&constructor);

        EXPECT_EQ(1u, graph.num_nodes());
    }
}

TEST(DBGBitmap, AddSequenceSimplePaths) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'C'));
        DBGBitmap graph(&constructor);

        EXPECT_EQ(2u, graph.num_nodes());
    }
}

TEST(DBGBitmap, NonASCIIStrings) {
    DBGBitmapConstructor constructor_first(6);
    constructor_first.add_sequences({
        // cyrillic A and C
        "АСАСАСАСАСАСА",
        "плохая строка",
        "АСАСАСАСАСАСА"
    });
    DBGBitmap graph(&constructor_first);
    ASSERT_EQ(1u, graph.num_nodes());
}

TEST(DBGBitmap, AddSequence) {
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("CAAC");
        DBGBitmap graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("CAAC");
        constructor.add_sequence("GAAC");
        DBGBitmap graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("AACG");
        DBGBitmap graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(5);
        constructor.add_sequence("AGACT");
        constructor.add_sequence("GACTT");
        constructor.add_sequence("ACTAT");
        DBGBitmap graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
}

TEST(DBGBitmap, AddSequences) {
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequences({
            "AAAC",
            "CAAC"
        });
        DBGBitmap graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequences({
           "AAAC",
           "CAAC",
           "GAAC"
        });
        DBGBitmap graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(4);
        constructor.add_sequences({
            "AAAC",
            "AACG"
        });
        DBGBitmap graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGBitmapConstructor constructor(5);
        constructor.add_sequences({
           "AGACT",
           "GACTT",
           "ACTAT"
        });
        DBGBitmap graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
}

TEST(DBGBitmap, CallPathsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor empty_const(k);
        DBGBitmap empty(&empty_const);
        DBGBitmapConstructor empty_reconst(k);

        uint64_t nseq = 0;
        empty.call_sequences([&](const auto &sequence) {
            empty_reconst.add_sequence(sequence);
            ++nseq;
        }, false);
        ASSERT_EQ(0u, nseq);
        DBGBitmap reconstructed(&empty_reconst);

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(DBGBitmap, CallContigsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor empty_const(k);
        DBGBitmap empty(&empty_const);
        DBGBitmapConstructor empty_reconst(k);

        uint64_t nseq = 0;
        empty.call_sequences([&](const auto &sequence) {
            empty_reconst.add_sequence(sequence);
        }, true);
        ASSERT_EQ(0u, nseq);
        DBGBitmap reconstructed(&empty_reconst);

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(DBGBitmap, CallPathsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; }, false);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGBitmap, CallContigsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_sequences([&](const auto &) { num_sequences++; }, true);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGBitmap, CallPathsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(3u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; }, false);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGBitmap, CallContigsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(3u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_sequences([&](const auto &) { num_sequences++; }, true);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGBitmap, CallPaths) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACACTAG");
            constructor.add_sequence("AACGACATG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);

            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AGACACTGA");
            constructor.add_sequence("GACTACGTA");
            constructor.add_sequence("ACTAACGTA");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AGACACAGT");
            constructor.add_sequence("GACTTGCAG");
            constructor.add_sequence("ACTAGTCAG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACTCGTAGC");
            constructor.add_sequence("AAATGCGTAGC");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACT");
            constructor.add_sequence("AAATG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
    }
}

TEST(DBGBitmap, CallContigs) {
    for (size_t k = 2; k <= 10; ++k) {
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACACTAG");
            constructor.add_sequence("AACGACATG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);

            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AGACACTGA");
            constructor.add_sequence("GACTACGTA");
            constructor.add_sequence("ACTAACGTA");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AGACACAGT");
            constructor.add_sequence("GACTTGCAG");
            constructor.add_sequence("ACTAGTCAG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACTCGTAGC");
            constructor.add_sequence("AAATGCGTAGC");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGBitmapConstructor constructor(k);
            constructor.add_sequence("AAACT");
            constructor.add_sequence("AAATG");
            DBGBitmap graph(&constructor);

            DBGBitmapConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGBitmap reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
    }
}

TEST(DBGBitmap, CallKmersEmptyGraph) {
    for (size_t k = 2; k <= 30; ++k) {
        DBGBitmapConstructor constructor(k);
        DBGBitmap empty(&constructor);

        size_t num_kmers = 0;
        empty.call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGBitmap graph(&constructor);

        ASSERT_EQ(3u, graph.num_nodes());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersFourLoopsDynamic) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'G'));
        constructor.add_sequence(std::string(100, 'C'));
        DBGBitmap graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersTestPath) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(k, 'C'));
        DBGBitmap graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k + 1, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersTestPathACA) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A')
                            + std::string(k, 'C')
                            + std::string(100, 'A'));
        DBGBitmap graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersTestPathDisconnected) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'T'));
        DBGBitmap graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TEST(DBGBitmap, CallKmersTestPathDisconnected2) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'G'));
        constructor.add_sequence(std::string(k, 'A') + "T");
        DBGBitmap graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(3u, num_kmers);
    }
}

/*
void test_pred_kmer(const DBGBitmap &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::vector<DBGBitmap::TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(),
                   [&graph](char c) {
                       return c == DBGBitmap::kSentinel
                                   ? DBGBitmap::kSentinelCode
                                   : graph.encode(c);
                   });
    EXPECT_EQ(expected_idx, graph.select_last(graph.pred_kmer(kmer)))
        << kmer_s << std::endl
       ;
}

TEST(DBGBitmap, PredKmer) {
    {
        DBGBitmap graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        DBGBitmap graph(5);
        graph.add_sequence("AAAAAA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGBitmap graph(5);
        graph.add_sequence("ACACAA");

        test_pred_kmer(graph, "ACGCG", 8);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 8);
        test_pred_kmer(graph, "NNNNN", 8);
        test_pred_kmer(graph, "$$$$$", 2);
    }
#ifndef _PROTEIN_GRAPH
    {
        DBGBitmap graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGBitmap graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGBitmap graph(2);
        graph.add_sequence("ATAATATCC");
        graph.add_sequence("ATACGC");
        graph.add_sequence("ATACTC");
        graph.add_sequence("ATACTA");
        graph.add_sequence("CATT");
        graph.add_sequence("CCC");
        graph.add_sequence("GGGC");
        graph.add_sequence("GGTGTGAC");
        graph.add_sequence("GGTCT");
        graph.add_sequence("GGTA");

        test_pred_kmer(graph, "$$", 4);
        test_pred_kmer(graph, "$A", 5);
        test_pred_kmer(graph, "$T", 26);
        test_pred_kmer(graph, "AT", 29);
        test_pred_kmer(graph, "TT", 35);
        test_pred_kmer(graph, "NT", 35);
        test_pred_kmer(graph, "TN", 35);
    }
#endif
}

TEST(DBGBitmap, PredKmerRandomTest) {
    srand(1);

    for (size_t k = 1; k < 8; ++k) {
        DBGBitmap graph(k);

        for (size_t p = 0; p < 10; ++p) {
            size_t length = rand() % 400;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = graph.alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = graph.alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, true);
        }

        auto all_kmer_str = utils::generate_strings("ACGT", k);
        for (size_t i = 1; i < k; ++i) {
            auto kmer_str_suffices = utils::generate_strings("ACGT", i);
            for (size_t j = 0; j < kmer_str_suffices.size(); ++j) {
                all_kmer_str.push_back(std::string(k - i, '$')
                                        + kmer_str_suffices[j]);
            }
        }

        for (const auto &kmer_str : all_kmer_str) {
            std::vector<DBGBitmap::TAlphabet> kmer = graph.encode(kmer_str);

            uint64_t lower_bound = graph.select_last(graph.pred_kmer(kmer));

            EXPECT_FALSE(
                utils::colexicographically_greater(
                    graph.get_node_seq(lower_bound), kmer
                )
            )
              << "kmer: " << kmer_str << std::endl
              << "lower bound: " << lower_bound << std::endl
              << "which is: ".get_node_str(lower_bound) << std::endl;

            if (lower_bound < graph.get_W().size() - 1) {
                EXPECT_TRUE(
                    utils::colexicographically_greater(
                        graph.get_node_seq(lower_bound + 1), kmer
                    )
                );
            }
        }
    }
}
*/

TEST(DBGBitmap, FindSequenceDBG) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        EXPECT_FALSE(graph->find(std::string(k - 1, 'A')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0));
        EXPECT_TRUE(graph->find(std::string(k, 'A')));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 1));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 0.75));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 0.5));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k, 'A'), 0));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A')));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 1));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 0.75));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 0.5));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'A'), 0));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A')));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A'), 1));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A'), 0.75));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A'), 0.5));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'A'), 0));

        EXPECT_FALSE(graph->find(std::string(k - 1, 'C')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'C'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'C'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'C'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'C'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'C'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'C')));
        EXPECT_FALSE(graph->find(std::string(k, 'C'), 1));
        EXPECT_FALSE(graph->find(std::string(k, 'C'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k, 'C'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k, 'C'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k, 'C'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'C')));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'C'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'C'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'C'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'C'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'C'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'C')));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'C'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'C'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'C'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'C'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'C'), 0));

        std::string pattern = std::string(k - 1, 'A') + std::string(k - 1, 'C');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 2)));
        EXPECT_FALSE(graph->find(pattern, 0.0 / (k - 1) + kEps));
        EXPECT_TRUE(graph->find(pattern, 0.0 / (k - 1) - kEps));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k, 'A') + std::string(k, 'C');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 1)));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 1) + kEps));
        // we map (k+1)-mers to the graph edges
        // just 1 out of k+2 (k+1)-mers can be mapped here
        EXPECT_TRUE(graph->find(pattern, 1.0 / (k + 1)));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k + 1, 'A') + std::string(k + 1, 'C');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 3.0 / (k + 3)));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 3) + kEps));
        EXPECT_TRUE(graph->find(pattern, 2.0 / (k + 3)));
        EXPECT_TRUE(graph->find(pattern, 0));
    }
}

// TODO
/*
TEST(DBGBitmap, KmerMappingMode) {
    for (size_t k = 1; k < 10; ++k) {
        DBGBitmap graph(k);

        const std::string check(k + 1, 'A');
        graph.add_sequence(std::string(100, 'A'));

        // kmer mapping modes
        for (size_t h = 0; h < 3; ++h) {
            EXPECT_TRUE(graph.find(check, 1, h));

            std::string query(100, 'A');

            for (double a : { 1.0, 0.75, 0.5, 0.25, 0.0 }) {
                EXPECT_TRUE(graph.find(query, a, h)) << "Mode: " << h;
            }

            // number of mutations
            for (size_t n = 1; n <= 100; ++n) {
                std::string query(100, 'A');
                for (size_t i = 0; i < query.size(); i += query.size() / n) {
                    query[i] = 'B';
                }

                uint64_t num_good_kmers = 0;
                for (size_t i = 0; i + k + 1 <= query.size(); ++i) {
                    if (query.substr(i, k + 1) == check)
                        num_good_kmers++;
                }
                double good_fraction = static_cast<double>(num_good_kmers)
                                            / (query.size() - k);

                for (double a : { 1.0, 0.75, 0.5, 0.25, 0.0 }) {
                    EXPECT_EQ(a <= good_fraction, graph.find(query, a, h))
                        << "Mode: " << h;
                }
            }
        }
    }
}
*/

TEST(DBGBitmap, Traversals) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k);

        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });

        uint64_t it2;
        graph->map_to_nodes(std::string(k - 1, 'A') + "C", [&](auto i) { it2 = i; });
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse_back(it2, 'G'));
    }
}

TEST(DBGBitmap, TraversalsCanonical) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k, true);

        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::unique_ptr<DBGBitmap> graph { new DBGBitmap(&constructor) };

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph->map_to_nodes(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        map_to_nodes_sequentially(
            std::string(k, 'A'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        map_to_nodes_sequentially(
            std::string(k, 'T'),
            [&](auto i) {
                EXPECT_NE(i, it);
            }
        );

        uint64_t it2 = 0;
        graph->map_to_nodes(
            std::string(k - 1, 'A') + "C",
            [&](auto i) { it2 = i; }
        );
        map_to_nodes_sequentially(
            std::string(k - 1, 'A') + "C",
            [&](auto i) { EXPECT_EQ(it2, i); }
        );
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse_back(it2, 'G'));

        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DBGBitmap::npos, it);
        uint64_t it3 = 0;
        graph->map_to_nodes(
            std::string(k, 'C'),
            [&](auto i) {
                it3 = i;
                EXPECT_EQ(i, it);
            }
        );
        map_to_nodes_sequentially(
            std::string(k, 'C'),
            [&](auto i) {
                EXPECT_EQ(i, it);
            }
        );
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'A'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'T'));

        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DBGBitmap::npos, it);
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
        ASSERT_NE(DBGBitmap::npos, it2);
        EXPECT_EQ(DBGBitmap::npos, graph->traverse(it, 'A'));
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        EXPECT_EQ(it2, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'G'));

        graph->map_to_nodes(
            std::string(k - 1, 'G') + "T",
            [&](auto i) { it2 = i; }
        );
        ASSERT_NE(DBGBitmap::npos, it2);
        EXPECT_EQ(DBGBitmap::npos, graph->traverse_back(it2, 'G'));
        EXPECT_EQ(it3, graph->traverse(it2, 'C'));
        EXPECT_NE(DBGBitmap::npos, graph->traverse_back(it2, 'A'));
    }
}

TEST(DBGBitmap, TraversalsDBG) {
    const auto npos = DeBruijnGraph::npos;

    for (size_t k = 2; k < 11; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        uint64_t it = 0;

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_EQ(npos, it);

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));

        ASSERT_NE(npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));

        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));

        EXPECT_EQ(npos, graph->traverse(it, 'G'));
        EXPECT_EQ(npos, graph->traverse_back(it + 1, 'G'));
    }
}

TEST(DBGBitmap, OutgoingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                       + std::string(100, 'G'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        uint64_t it = 0;
        std::vector<DBGBitmap::node_index> adjacent_nodes;

        // AA, AAAAA
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ it, graph->traverse(it, 'C') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        auto outset = convert_to_set(std::vector<uint64_t>{ graph->traverse(it, 'C') });
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        if (k == 2) {
            outset.insert(graph->traverse(it, 'G'));
            ASSERT_EQ(2u, adjacent_nodes.size());
        } else {
            ASSERT_EQ(1u, adjacent_nodes.size());
        }

        EXPECT_EQ(outset, convert_to_set(adjacent_nodes));
        adjacent_nodes.clear();

        // CC, CCCCC
        graph->map_to_nodes(std::string(k, 'C'), [&](auto i) { it = i; });
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph->traverse(it, 'G')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CG, CCCCG
        it = graph->traverse(it, 'G');
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph->traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GGGGG
        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph->traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TEST(DBGBitmap, IncomingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C')
                                                       + std::string(100, 'G'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        uint64_t it = 0;
        std::vector<DBGBitmap::node_index> adjacent_nodes;

        // AA, AAAAA
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ it }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{ graph->traverse_back(it, 'A') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CC, CCCCC
        graph->map_to_nodes(std::string(k, 'C'), [&](auto i) { it = i; });
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph->traverse_back(it, 'A')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CG, CCCCG
        it = graph->traverse(it, 'G');
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                graph->traverse_back(it, 'A'),
                graph->traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GG, GGGGG
        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<uint64_t>{
                it,
                graph->traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TEST(DBGBitmap, TraversalsDBGCanonical) {
    const auto npos = DeBruijnGraph::npos;

    for (size_t k = 2; k < 11; ++k) {
        DBGBitmapConstructor constructor(k, true);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGBitmap(&constructor) };

        uint64_t it = 0;

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);
        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));
        ASSERT_NE(npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));
        EXPECT_EQ(npos, graph->traverse(it, 'G'));
        EXPECT_EQ(npos, graph->traverse_back(it + 1, 'G'));


        // reverse complement
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        ASSERT_NE(npos, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'T'), 'G'));


        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        EXPECT_EQ(npos, graph->traverse(it, 'G'));
        EXPECT_EQ(npos, graph->traverse(it, 'T'));
    }
}

TEST(DBGBitmap, map_to_nodes) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DBGBitmap> graph { new DBGBitmap(&constructor) };


        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 3, 'A')
                                        + std::string(2 * k, 'C');
        std::vector<uint64_t> expected_result {
            SequenceGraph::npos,
            SequenceGraph::npos
        };

        for (size_t i = 2; i + k <= sequence_to_map.size(); ++i) {
            graph->map_to_nodes(
                sequence_to_map.substr(i, k),
                [&](auto i) { expected_result.push_back(i);}
            );
        }

        std::vector<uint64_t> observed_result;
        graph->map_to_nodes(sequence_to_map,
            [&](const auto &index) {
                observed_result.emplace_back(index);
            }
        );
        EXPECT_EQ(expected_result, observed_result);

        size_t pos = 0;
        graph->map_to_nodes(sequence_to_map,
                            [&](auto i) { EXPECT_EQ(expected_result[pos++], i); });
    }
}

TEST(DBGBitmap, get_outdegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(k - 1, 'A') + 'C');
        std::unique_ptr<DBGBitmap> graph { new DBGBitmap(&constructor) };
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->outdegree(1));
    }
}

TEST(DBGBitmap, get_maximum_outdegree) {
    for (size_t k = 2; k < 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(k - 1, 'A') + 'A');
        constructor.add_sequence(std::string(k - 1, 'A') + 'C');
        constructor.add_sequence(std::string(k - 1, 'A') + 'G');
        constructor.add_sequence(std::string(k - 1, 'A') + 'T');
        std::unique_ptr<DBGBitmap> graph { new DBGBitmap(&constructor) };

        DBGBitmap::node_index max_outdegree_node_index;
        graph->map_to_nodes(std::string(k, 'A'), [&](DBGBitmap::node_index node) {
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

TEST(DBGBitmap, get_outdegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        DBGBitmapConstructor constructor(k);
        constructor.add_sequence(std::string(k - 1, 'A') + std::string(k - 1, 'C') +
                                 std::string(k - 1, 'G') + std::string(k, 'T'));
        constructor.add_sequence(std::string(k, 'A'));
        std::unique_ptr<DBGBitmap> graph { new DBGBitmap(&constructor) };

        DBGBitmap::node_index loop_node_index;
        graph->map_to_nodes(std::string(k, 'A'),
                            [&](DBGBitmap::node_index node) { loop_node_index = node; });
        for (size_t i = 1; i <= graph->num_nodes(); ++i) {
            if (i == loop_node_index)
                EXPECT_EQ(2ull, graph->outdegree(i));
            else
                EXPECT_EQ(1ull, graph->outdegree(i));
        }
    }
}
