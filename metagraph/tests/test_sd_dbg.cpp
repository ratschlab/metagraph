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

#include "dbg_sd.hpp"
#include "dbg_construct.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";

const KmerExtractor2Bit kmer_extractor;


TEST(DBGSDConstructedFull, InitializeComplete) {
    {
        DBGSD graph(20, false);
        EXPECT_EQ(std::string("AAAAAAAAAAAAAAAAAAAA"), graph.node_to_kmer(1));
        EXPECT_EQ(1u, graph.kmer_to_node("AAAAAAAAAAAAAAAAAAAA"));

    // #if _DNA4_GRAPH
        EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("TTTTTTTTTTTTTTTTTTTT"));
    // #elif _DNA_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("NNNNNNNNNNNNNNNNNNNN"));
    // #elif _DNA_CASE_SENSITIVE_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("tttttttttttttttttttt"));
    // #elif _PROTEIN_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("XXXXXXXXXXXXXXXXXXXX"));
    // #endif

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }

    {
        DBGSD graph(20, true);
        EXPECT_EQ(std::string("AAAAAAAAAAAAAAAAAAAA"), graph.node_to_kmer(1));
        EXPECT_EQ(1u, graph.kmer_to_node("AAAAAAAAAAAAAAAAAAAA"));

    // #if _DNA4_GRAPH
        EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("TTTTTTTTTTTTTTTTTTTT"));
    // #elif _DNA_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("NNNNNNNNNNNNNNNNNNNN"));
    // #elif _DNA_CASE_SENSITIVE_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("tttttttttttttttttttt"));
    // #elif _PROTEIN_GRAPH
    //     EXPECT_EQ(graph.num_nodes(), graph.kmer_to_node("XXXXXXXXXXXXXXXXXXXX"));
    // #endif

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGSDConstructedFull, SerializeComplete) {
    {
        DBGSD graph(20, false);

        graph.serialize(test_dump_basename);
    }

    {
        DBGSD graph(2, false);

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
        DBGSD graph(2, true);

        ASSERT_TRUE(graph.load(test_dump_basename));
        EXPECT_EQ(20u, graph.get_k());

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_TRUE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    }
}

TEST(DBGSDConstructedFull, InsertSequence) {
    {
        DBGSD graph(20, false);

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
        DBGSD graph(20, true);

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

TEST(DBGSDConstructedFull, ReverseComplement) {
    {
        DBGSD graph(20, false);

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
        DBGSD graph(20, true);

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

TEST(DBGSDConstructedFull, CheckGraph) {
    {
        DBGSD graph(20, false);

        const std::string alphabet = "ACGT";

        for (size_t i = 0; i < 100; ++i) {
            std::string seq(1'000, 'A');
            for (size_t j = 0; j < seq.size(); ++j) {
                seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
            }
            ASSERT_DEATH(graph.add_sequence(seq), "");
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
            ASSERT_DEATH(graph.add_sequence(seq), "");
        }

        for (uint64_t i = 1; i <= 1000; ++i) {
            auto kmer = graph.node_to_kmer(i);
            ASSERT_EQ(20u, kmer.size());
            auto node = graph.kmer_to_node(kmer);
            ASSERT_EQ(i, node) << kmer;
        }
    }
}

TEST(DBGSDConstructedFull, Serialize) {
    {
        DBGSD graph(20, false);

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
        DBGSD graph(2, false);

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
        DBGSD graph(2, true);

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

TEST(DBGSDConstructed, InsertSequence) {
    DBGSDConstructor constructor(20);
    constructor.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    constructor.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");
    DBGSD graph(20);
    constructor.build_graph(&graph);

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TEST(DBGSDConstructed, CheckGraph) {
    DBGSDConstructor constructor(20);

    const std::string alphabet = "ACGTN";

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        constructor.add_sequence(seq);
    }
    DBGSD graph(20);
    constructor.build_graph(&graph);

    for (size_t i = 0; i < 100; ++i) {
        std::string seq(1'000, 'A');
        for (size_t j = 0; j < seq.size(); ++j) {
            seq[j] = alphabet[(i * i + j + 17 * j * j) % 5];
        }
        EXPECT_TRUE(graph.find(seq));
    }
}

TEST(DBGSDConstructed, Serialize) {
    {
        DBGSDConstructor constructor(20);
        constructor.add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        constructor.add_sequence("CATGTACTAGCTGATCGTAGCTAGCTAGC");
        DBGSD graph(20);
        constructor.build_graph(&graph);

        EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));

        graph.serialize(test_dump_basename);
    }

    DBGSD graph(20);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

constexpr double kEps = std::numeric_limits<double>::epsilon();

//TODO
/*
void test_graph(DBGSD *graph, const bit_vector_sd &ref, bool canonical_only) {
    DBGSD ngraph(2);
    ngraph.k_ = ngraph.infer_k();
    EXPECT_EQ(ref, graph->kmers_);
    EXPECT_EQ(canonical_only, graph->canonical_only_);
}
*/


TEST(DBGSD, GraphDefaultConstructor) {
    DBGSD *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new DBGSD(2);
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        DBGSD graph(2);
    });
}

/*

TEST(DBGSD, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    DBGSDConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_sequences({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBGSD *graph = new DBGSD(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_nodes() + 1);

    EXPECT_EQ(1u, graph->outgoing(1, DBGSD::kSentinelCode));

    for (size_t i = 1; i <= graph->num_edges(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBGSD::kSentinelCode) {
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

TEST(DBGSD, AddSequenceSimplePath) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        DBGSD graph(&constructor);

        EXPECT_EQ(1u, graph.num_nodes());
    }
}

TEST(DBGSD, AddSequenceSimplePaths) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'C'));
        DBGSD graph(&constructor);

        EXPECT_EQ(2u, graph.num_nodes());
    }
}

TEST(DBGSD, NonASCIIStrings) {
    DBGSDConstructor constructor_first(6);
    constructor_first.add_sequences({
        // cyrillic A and C
        "АСАСАСАСАСАСА",
        "плохая строка",
        "АСАСАСАСАСАСА"
    });
    DBGSD graph(&constructor_first);
    ASSERT_EQ(1u, graph.num_nodes());
}

TEST(DBGSD, AddSequence) {
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("CAAC");
        DBGSD graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("CAAC");
        constructor.add_sequence("GAAC");
        DBGSD graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequence("AAAC");
        constructor.add_sequence("AACG");
        DBGSD graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(5);
        constructor.add_sequence("AGACT");
        constructor.add_sequence("GACTT");
        constructor.add_sequence("ACTAT");
        DBGSD graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
}

TEST(DBGSD, AddSequences) {
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequences({
            "AAAC",
            "CAAC"
        });
        DBGSD graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequences({
           "AAAC",
           "CAAC",
           "GAAC"
        });
        DBGSD graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(4);
        constructor.add_sequences({
            "AAAC",
            "AACG"
        });
        DBGSD graph(&constructor);
        EXPECT_EQ(2u, graph.num_nodes());
    }
    {
        DBGSDConstructor constructor(5);
        constructor.add_sequences({
           "AGACT",
           "GACTT",
           "ACTAT"
        });
        DBGSD graph(&constructor);
        EXPECT_EQ(3u, graph.num_nodes());
    }
}

TEST(DBGSD, CallPathsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor empty_const(k);
        DBGSD empty(&empty_const);
        DBGSDConstructor empty_reconst(k);

        uint64_t nseq = 0;
        empty.call_sequences([&](const auto &sequence) {
            empty_reconst.add_sequence(sequence);
            ++nseq;
        }, false);
        ASSERT_EQ(0u, nseq);
        DBGSD reconstructed(&empty_reconst);

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(DBGSD, CallContigsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor empty_const(k);
        DBGSD empty(&empty_const);
        DBGSDConstructor empty_reconst(k);

        uint64_t nseq = 0;
        empty.call_sequences([&](const auto &sequence) {
            empty_reconst.add_sequence(sequence);
        }, true);
        ASSERT_EQ(0u, nseq);
        DBGSD reconstructed(&empty_reconst);

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(DBGSD, CallPathsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGSD graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; }, false);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGSD, CallContigsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGSD graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_sequences([&](const auto &) { num_sequences++; }, true);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGSD, CallPathsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGSD graph(&constructor);

        ASSERT_EQ(3u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; }, false);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGSD, CallContigsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGSD graph(&constructor);

        ASSERT_EQ(3u, graph.num_nodes());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_sequences([&](const auto &) { num_sequences++; }, true);

        EXPECT_EQ(graph.num_nodes(), num_paths);
        EXPECT_EQ(graph.num_nodes(), num_sequences);
    }
}

TEST(DBGSD, CallPaths) {
    for (size_t k = 3; k <= 10; ++k) {
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACACTAG");
            constructor.add_sequence("AACGACATG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);

            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AGACACTGA");
            constructor.add_sequence("GACTACGTA");
            constructor.add_sequence("ACTAACGTA");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AGACACAGT");
            constructor.add_sequence("GACTTGCAG");
            constructor.add_sequence("ACTAGTCAG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACTCGTAGC");
            constructor.add_sequence("AAATGCGTAGC");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACT");
            constructor.add_sequence("AAATG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, false);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
    }
}

TEST(DBGSD, CallContigs) {
    for (size_t k = 2; k <= 10; ++k) {
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACACTAG");
            constructor.add_sequence("AACGACATG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);

            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AGACACTGA");
            constructor.add_sequence("GACTACGTA");
            constructor.add_sequence("ACTAACGTA");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AGACACAGT");
            constructor.add_sequence("GACTTGCAG");
            constructor.add_sequence("ACTAGTCAG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACTCGTAGC");
            constructor.add_sequence("AAATGCGTAGC");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
        {
            DBGSDConstructor constructor(k);
            constructor.add_sequence("AAACT");
            constructor.add_sequence("AAATG");
            DBGSD graph(&constructor);

            DBGSDConstructor reconst(k);
            graph.call_sequences([&](const auto &sequence) {
                reconst.add_sequence(sequence);
            }, true);
            DBGSD reconstructed(&reconst);

            ASSERT_EQ(graph, reconstructed);
        }
    }
}

TEST(DBGSD, CallKmersEmptyGraph) {
    for (size_t k = 2; k <= 30; ++k) {
        DBGSDConstructor constructor(k);
        DBGSD empty(&constructor);

        size_t num_kmers = 0;
        empty.call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TEST(DBGSD, CallKmersTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        DBGSD graph(&constructor);

        ASSERT_EQ(1u, graph.num_nodes());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TEST(DBGSD, CallKmersFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        DBGSD graph(&constructor);

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

TEST(DBGSD, CallKmersFourLoopsDynamic) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'G'));
        constructor.add_sequence(std::string(100, 'C'));
        DBGSD graph(&constructor);

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

TEST(DBGSD, CallKmersTestPath) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(k, 'C'));
        DBGSD graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k + 1, num_kmers);
    }
}

TEST(DBGSD, CallKmersTestPathACA) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A')
                            + std::string(k, 'C')
                            + std::string(100, 'A'));
        DBGSD graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k, num_kmers);
    }
}

TEST(DBGSD, CallKmersTestPathDisconnected) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        constructor.add_sequence(std::string(100, 'T'));
        DBGSD graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TEST(DBGSD, CallKmersTestPathDisconnected2) {
    for (size_t k = 2; k <= 20; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'G'));
        constructor.add_sequence(std::string(k, 'A') + "T");
        DBGSD graph(&constructor);

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(3u, num_kmers);
    }
}

/*
void test_pred_kmer(const DBGSD &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::vector<DBGSD::TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(),
                   [&graph](char c) {
                       return c == DBGSD::kSentinel
                                   ? DBGSD::kSentinelCode
                                   : graph.encode(c);
                   });
    EXPECT_EQ(expected_idx, graph.select_last(graph.pred_kmer(kmer)))
        << kmer_s << std::endl
       ;
}

TEST(DBGSD, PredKmer) {
    {
        DBGSD graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        DBGSD graph(5);
        graph.add_sequence("AAAAAA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGSD graph(5);
        graph.add_sequence("ACACAA");

        test_pred_kmer(graph, "ACGCG", 8);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 8);
        test_pred_kmer(graph, "NNNNN", 8);
        test_pred_kmer(graph, "$$$$$", 2);
    }
#ifndef _PROTEIN_GRAPH
    {
        DBGSD graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGSD graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBGSD graph(2);
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

TEST(DBGSD, PredKmerRandomTest) {
    srand(1);

    for (size_t k = 1; k < 8; ++k) {
        DBGSD graph(k);

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
            std::vector<DBGSD::TAlphabet> kmer = graph.encode(kmer_str);

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

TEST(DBGSD, FindSequenceDBG) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGSD(&constructor) };

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
TEST(DBGSD, KmerMappingMode) {
    for (size_t k = 1; k < 10; ++k) {
        DBGSD graph(k);

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

TEST(DBGSD, Traversals) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k);

        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        auto graph = std::make_unique<DBGSD>(&constructor);

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });

        uint64_t it2 = graph->kmer_to_node(std::string(k - 1, 'A') + "C");
        ASSERT_EQ(graph->kmer_to_node(std::string(k, 'A')), it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DBGSD::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBGSD::npos, graph->traverse_back(it2, 'G'));
    }
}

TEST(DBGSD, TraversalsCanonical) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k, true);

        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        auto graph = std::make_unique<DBGSD>(&constructor);

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(graph->kmer_to_node(std::string(k, 'T')), it);

        uint64_t it2 = graph->kmer_to_node(std::string(k - 1, 'A') + "C");
        ASSERT_EQ(graph->kmer_to_node(std::string(k, 'A')), it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DBGSD::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBGSD::npos, graph->traverse_back(it2, 'G'));

        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_EQ(graph->kmer_to_node(std::string(k, 'G')), it);
        ASSERT_NE(graph->kmer_to_node(std::string(k, 'C')), it);
        ASSERT_NE(DBGSD::npos, it);
        it2 = graph->kmer_to_node(std::string(k - 1, 'G') + "T");
        ASSERT_NE(DBGSD::npos, it2);
        EXPECT_EQ(DBGSD::npos, graph->traverse(it, 'A'));
        EXPECT_EQ(it, graph->traverse(it, 'G'));
        EXPECT_EQ(it2, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'G'));
    }
}

TEST(DBGSD, TraversalsDBG) {
    const auto npos = DeBruijnGraph::npos;

    for (size_t k = 2; k < 11; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGSD(&constructor) };

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

TEST(DBGSD, TraversalsDBGCanonical) {
    const auto npos = DeBruijnGraph::npos;

    for (size_t k = 2; k < 11; ++k) {
        DBGSDConstructor constructor(k, true);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DeBruijnGraph> graph { new DBGSD(&constructor) };

        uint64_t it = 0;

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));

        ASSERT_NE(npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));

        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));

        EXPECT_EQ(npos, graph->traverse(it, 'G'));
        EXPECT_EQ(npos, graph->traverse_back(it + 1, 'G'));

        // reverse complement
        graph->map_to_nodes(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'G'))
            << *dynamic_cast<DBGSD*>(graph.get());
        EXPECT_NE(npos, graph->traverse(it, 'T'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'T'), 'G'));
    }
}

TEST(DBGSD, map_to_nodes) {
    for (size_t k = 2; k <= 10; ++k) {
        DBGSDConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        std::unique_ptr<DBGSD> graph { new DBGSD(&constructor) };


        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 3, 'A')
                                        + std::string(2 * k, 'C');
        std::vector<uint64_t> expected_result {
            SequenceGraph::npos,
            SequenceGraph::npos
        };

        for (size_t i = 2; i + k <= sequence_to_map.size(); ++i) {
            expected_result.push_back(
                graph->kmer_to_node(sequence_to_map.substr(i, k))
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
