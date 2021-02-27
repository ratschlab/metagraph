#include "gtest/gtest.h"

#include "../test_helpers.hpp"
#include "all/test_dbg_helpers.hpp"

#include "common/seq_tools/reverse_complement.hpp"


namespace {

#if ! _PROTEIN_GRAPH

using namespace mtg;
using namespace mtg::test;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


template <typename Graph>
class DeBruijnGraphCanonicalTest : public DeBruijnGraphTest<Graph> { };
// TODO: add support for canonical mode in DBGHashString
typedef ::testing::Types<DBGBitmap,
                         DBGHashOrdered,
                         DBGHashFast,
                         DBGSuccinct,
                         DBGSuccinctBloom<4, 1>> CanonicalGraphTypes;
TYPED_TEST_SUITE(DeBruijnGraphCanonicalTest, CanonicalGraphTypes);


TYPED_TEST(DeBruijnGraphCanonicalTest, CheckGraph) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGT", DeBruijnGraph::CANONICAL, true));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CheckGraphInputWithN) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGTN", DeBruijnGraph::CANONICAL, false));
    EXPECT_EQ(TypeParam(3).alphabet().find('N') != std::string::npos,
              check_graph<TypeParam>("ACGTN", DeBruijnGraph::CANONICAL, true));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, InitializeEmpty) {
    auto graph = build_graph<TypeParam>(2, {}, DeBruijnGraph::CANONICAL);

    EXPECT_EQ(0u, graph->num_nodes());
    EXPECT_FALSE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, SerializeEmpty) {
    {
        auto graph = build_graph<TypeParam>(20, {}, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(0u, graph->num_nodes());
        graph->serialize(test_dump_basename);
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(20u, graph.get_k());
    EXPECT_EQ(DeBruijnGraph::CANONICAL, graph.get_mode());

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
        }, DeBruijnGraph::CANONICAL);

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
    }, DeBruijnGraph::CANONICAL);

    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_TRUE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphCanonicalTest, ReverseComplement) {
    auto graph1 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" }, DeBruijnGraph::CANONICAL);
    auto graph2 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                               "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT" }, DeBruijnGraph::CANONICAL);

    EXPECT_EQ(graph1->num_nodes(), graph2->num_nodes());

    auto graph = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                              "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                                              "CATGTACTAGCTGATCGTAGCTAGCTAGC" }, DeBruijnGraph::CANONICAL);
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
        }, DeBruijnGraph::CANONICAL);

        auto it = DeBruijnGraph::npos;
        graph->map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->map_to_nodes_sequentially(
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

        graph->map_to_nodes_sequentially(
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
        }, DeBruijnGraph::CANONICAL);

        auto it = DeBruijnGraph::npos;

        graph->map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        graph->map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));

        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));

        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));

        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it, 'G'));

        // reverse complement
        graph->map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
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
        }, DeBruijnGraph::CANONICAL);

        auto it = DeBruijnGraph::npos;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(
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
        graph->map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->map_to_nodes_sequentially(
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

        graph->map_to_nodes_sequentially(
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
        }, DeBruijnGraph::CANONICAL);

        auto it = DeBruijnGraph::npos;

        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_NE(DeBruijnGraph::npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it + 1, 'G'));


        // reverse complement
        graph->map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { it = i; });
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
        }, DeBruijnGraph::CANONICAL);

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

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsEmptyGraphCanonical) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            auto empty = build_graph<TypeParam>(k, {}, DeBruijnGraph::CANONICAL);
            std::vector<std::string> sequences;
            std::mutex seq_mutex;
            empty->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*empty, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                sequences.push_back(sequence);
            }, num_threads);
            ASSERT_EQ(0u, sequences.size());

            EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL));
            EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL));
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsEmptyGraph) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            auto empty = build_graph<TypeParam>(k, {}, DeBruijnGraph::CANONICAL);
            std::vector<std::string> sequences;
            std::mutex seq_mutex;
            empty->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*empty, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                sequences.push_back(sequence);
            }, num_threads);
            ASSERT_EQ(0u, sequences.size());

            EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL));
            EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL));
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsOneSelfLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 20; ++k) {
            std::vector<std::string> sequences { std::string(100, 'A') };
            auto graph = build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            ASSERT_EQ(2u, graph->num_nodes());
            ASSERT_EQ(2u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsOneSelfLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 20; ++k) {
            std::vector<std::string> sequences { std::string(100, 'A') };
            auto graph = build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            ASSERT_EQ(2u, graph->num_nodes());
            ASSERT_EQ(2u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsThreeSelfLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 20; ++k) {
            std::vector<std::string> sequences { std::string(100, 'A'),
                                                 std::string(100, 'G'),
                                                 std::string(100, 'C') };
            auto graph = build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            ASSERT_EQ(4u, graph->num_nodes());
            ASSERT_EQ(4u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsExtractsLongestOneLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 6; k < 14; ++k) {
            std::vector<std::string> sequences { "ATGCCGTACTCAG",
                                                 "GGGGGGGGGGGGG" };
            auto graph = build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

            std::vector<std::string> contigs;
            std::mutex seq_mutex;
            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                contigs.push_back(sequence);
            }, num_threads);

            EXPECT_EQ(4u, contigs.size());
            EXPECT_EQ(convert_to_set({ "ATGCCGTACTCAG", std::string(k, 'G'),
                                       "CTGAGTACGGCAT", std::string(k, 'C') }),
                      convert_to_set(contigs)) << k;
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallContigsUniqueKmers) {
    for (size_t num_threads : { 1, 4 }) {
        std::string sequence = "GCAAATAAC";
        auto graph = build_graph<TypeParam>(3, { sequence }, DeBruijnGraph::CANONICAL);

        std::atomic<size_t> num_kmers = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_kmers += sequence.size() - 2;
        }, num_threads);

        EXPECT_EQ((sequence.size() - 2) * 2, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsUniqueKmersCycle) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 4;
        std::string sequence = "AAACCCGGGTTTAAA";
        auto graph = build_graph<TypeParam>(k, { sequence }, DeBruijnGraph::CANONICAL);

        std::atomic<size_t> num_unitigs = 0;
        std::atomic<size_t> num_kmers = 0;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_unitigs++;
            num_kmers += sequence.size() - k + 1;
        }, num_threads);

        EXPECT_EQ(1u, num_unitigs);
        EXPECT_EQ(sequence.size() - k + 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallContigsUniqueKmersCycle) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 4;
        std::string sequence = "AAACCCGGGTTTAAA";
        auto graph = build_graph<TypeParam>(k, { sequence }, DeBruijnGraph::CANONICAL);

        std::atomic<size_t> num_contigs = 0;
        std::atomic<size_t> num_kmers = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_contigs++;
            num_kmers += sequence.size() - k + 1;
        }, num_threads);

        EXPECT_EQ(1u, num_contigs);
        EXPECT_EQ(sequence.size() - k + 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsFourLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 20; ++k) {
            std::vector<std::string> sequences { std::string(100, 'A'),
                                                 std::string(100, 'G'),
                                                 std::string(100, 'C') };
            auto graph = build_graph<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            auto graph_batch = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);
            ASSERT_EQ(4u, graph->num_nodes());
            ASSERT_EQ(4u, graph_batch->num_nodes());

            std::atomic<size_t> num_sequences = 0;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                num_sequences++;
            }, num_threads);
            std::atomic<size_t> num_sequences_batch = 0;
            graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
                num_sequences_batch++;
            }, num_threads);

            EXPECT_EQ(graph->num_nodes(), num_sequences);
            EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
            EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPaths) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::CANONICAL);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_sequences([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads);
                    },
                    DeBruijnGraph::CANONICAL
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::CANONICAL);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_sequences([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads, true);
                    },
                    DeBruijnGraph::CANONICAL
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallPathsCheckHalfSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 15; k += 2) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::atomic<size_t> num_kmers_both = 0;
                graph->call_sequences([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    num_kmers_both += path.size();
                }, num_threads);

                std::atomic<size_t> num_kmers = 0;
                graph->call_sequences([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    num_kmers += path.size();
                }, num_threads, true);

                EXPECT_EQ(num_kmers_both, num_kmers * 2);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigs) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::CANONICAL);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_unitigs([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads);
                    },
                    DeBruijnGraph::CANONICAL
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 2; k <= 10; ++k) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                // in stable graphs the order of input sequences
                // does not change the order of k-mers and their indexes
                auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences, DeBruijnGraph::CANONICAL);

                std::mutex seq_mutex;
                auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                    k,
                    [&](const auto &callback) {
                        graph->call_unitigs([&](const auto &sequence, const auto &path) {
                            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                            std::unique_lock<std::mutex> lock(seq_mutex);
                            callback(sequence);
                        }, num_threads, 1, true);
                    },
                    DeBruijnGraph::CANONICAL
                );

                EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsCheckHalfSingleKmerForm) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 3; k <= 15; k += 2) {
            for (const std::vector<std::string> &sequences
                    : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                        std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                        std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                        std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                        std::vector<std::string>({ "AAACT", "AAATG" }),
                        std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

                auto graph = build_graph_batch<TypeParam>(k, sequences, DeBruijnGraph::CANONICAL);

                std::atomic<size_t> num_kmers_both = 0;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    num_kmers_both += path.size();
                }, num_threads);

                std::atomic<size_t> num_kmers = 0;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    num_kmers += path.size();
                }, num_threads, 1, true);

                EXPECT_EQ(num_kmers_both, num_kmers * 2);
            }
        }
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsWithoutTips) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 3;
        std::mutex seq_mutex;
        std::set<std::string> unitigs;

        auto graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                                 "TCTAAGC" }, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(12u, graph->num_nodes());

        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGA", "AGCT", "AGT", "CTA", "CTTA",
                                          "TAAG", "TAG", "TCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGA", "AGCT", "AGT", "CTA", "CTTA",
                                          "TAAG", "TAG", "TCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGA", "AGCT", "AGT", "CTA", "CTTA",
                                          "TAAG", "TAG", "TCT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGA", "AGCT", "AGT", "CTA", "CTTA",
                                          "TAAG", "TAG", "TCT" }), unitigs);

        graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                            "ACTAAGT" }, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(10u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGCT", "AGT", "CTA", "CTTA", "TAAG", "TAG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGCT", "AGT", "CTA", "CTTA", "TAAG", "TAG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGCT", "AGT", "CTA", "CTTA", "TAAG", "TAG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "ACT", "AGCT", "AGT", "CTA", "CTTA", "TAAG", "TAG" }), unitigs);

        graph = build_graph<TypeParam>(k, { "ACTAAGCCC",
                                            "AAAGC",
                                            "TAAGCA" }, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(18u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAA", "AAG", "ACT", "AGC", "AGT", "CCC",
                                          "CTA", "CTT", "GCA", "GCC", "GCT", "GGC",
                                          "GGG", "TAA", "TAG", "TGC", "TTA", "TTT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAA", "AAG", "ACT", "AGC", "AGT", "CCC",
                                          "CTA", "CTT", "GCA", "GCC", "GCT", "GGC",
                                          "GGG", "TAA", "TAG", "TGC", "TTA", "TTT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC", "GCC", "CCC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAA", "AAG", "ACT", "AGC", "AGT", "CCC",
                                          "CTA", "CTT", "GCA", "GCC", "GCT", "GGC",
                                          "GGG", "TAA", "TAG", "TGC", "TTA", "TTT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAA", "AAG", "ACT", "AGC", "AGT", "CCC",
                                          "CTA", "CTT", "GCA", "GCC", "GCT", "GGC",
                                          "GGG", "TAA", "TAG", "TGC", "TTA", "TTT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        // EXPECT_EQ(std::set<std::string>({ "AAGC" }), unitigs);
        EXPECT_EQ(std::set<std::string>({ "AAA", "AAG", "ACT", "AGC", "AGT", "CCC",
                                          "CTA", "CTT", "GCA", "GCC", "GCT", "GGC",
                                          "GGG", "TAA", "TAG", "TGC", "TTA", "TTT" }), unitigs);

        graph = build_graph<TypeParam>(k, { "ACGAAGCCT",
                                            "AAGC",
                                            "TAAGCA" }, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(18u, graph->num_nodes());

        // TODO: make DBGSuccinct work properly even if it has redundant source dummy edges
        if (dynamic_cast<DBGSuccinct*>(graph.get()))
            dynamic_cast<DBGSuccinct&>(*graph).mask_dummy_kmers(1, true);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "AAG", "ACG", "AGC", "AGGC", "CGAA",
                                          "CGT", "CTT", "GCA", "GCCT", "GCT",
                                          "TGC", "TTAA", "TTCG" }), unitigs) << *graph;

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "AAG", "ACG", "AGC", "AGGC", "CGAA",
                                          "CGT", "CTT", "GCA", "GCCT", "GCT",
                                          "TGC", "TTAA", "TTCG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "AAG", "ACG", "AGC", "AGGC", "CGAA",
                                          "CGT", "CTT", "GCA", "GCCT", "GCT",
                                          "TGC", "TTAA", "TTCG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        EXPECT_EQ(std::set<std::string>({ "AAG", "ACG", "AGC", "AGGC", "CGAA",
                                          "CGT", "CTT", "GCA", "GCCT", "GCT",
                                          "TGC", "TTAA", "TTCG" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "AAG", "ACG", "AGC", "AGGC", "CGAA",
                                          "CGT", "CTT", "GCA", "GCCT", "GCT",
                                          "TGC", "TTAA", "TTCG" }), unitigs);

        graph = build_graph<TypeParam>(k, { "TCTAAGCCG",
                                            "CATAAGCCG",
                                            "CATAACCGA" }, DeBruijnGraph::CANONICAL);
        ASSERT_EQ(24u, graph->num_nodes());

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "AACC", "AAG", "AGA", "AGC", "ATA", "ATG",
                                          "CAT", "CCG", "CGA", "CGG", "CTA", "CTT",
                                          "GCC", "GCT", "GGC", "GGTT", "TAA", "TAG",
                                          "TAT", "TCG", "TCT", "TTA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "AACC", "AAG", "AGA", "AGC", "ATA", "ATG",
                                          "CAT", "CCG", "CGA", "CGG", "CTA", "CTT",
                                          "GCC", "GCT", "GGC", "GGTT", "TAA", "TAG",
                                          "TAT", "TCG", "TCT", "TTA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "AACC", "AAG", "AGA", "AGC", "ATA", "ATG",
                                          "CAT", "CCG", "CGA", "CGG", "CTA", "CTT",
                                          "GCC", "GCT", "GGC", "GGTT", "TAA", "TAG",
                                          "TAT", "TCG", "TCT", "TTA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 3);
        EXPECT_EQ(std::set<std::string>({ "AACC", "AAG", "AGA", "AGC", "ATA", "ATG",
                                          "CAT", "CCG", "CGA", "CGG", "CTA", "CTT",
                                          "GCC", "GCT", "GGC", "GGTT", "TAA", "TAG",
                                          "TAT", "TCG", "TCT", "TTA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "AACC", "AAG", "AGA", "AGC", "ATA", "ATG",
                                          "CAT", "CCG", "CGA", "CGG", "CTA", "CTT",
                                          "GCC", "GCT", "GGC", "GGTT", "TAA", "TAG",
                                          "TAT", "TCG", "TCT", "TTA" }), unitigs);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsWithoutTips2) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 5;
        auto graph = build_graph<TypeParam>(k, { "ACTATAGCTAGTCTATGCGA",
                                                 "ACTATAGCTAGTCTAA",
                                                 "ACTATAGCTA",
                                                 "ACTATAGCTT",
                                                 "ACTATC", }, DeBruijnGraph::CANONICAL);
        std::mutex seq_mutex;
        ASSERT_EQ(34u, graph->num_nodes());
        std::set<std::string> unitigs;
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 0);
        EXPECT_EQ(std::set<std::string>({ "AAGCT", "ACTAG", "ACTAT", "AGCTA", "AGCTT",
                                          "ATAGA", "ATAGC", "ATAGT", "CTAGC", "CTAGT",
                                          "CTATAG", "CTATC", "CTATGCGA", "GATAG",
                                          "GCTAG", "GCTAT", "TAGACTA", "TAGCT",
                                          "TAGTCTA", "TCGCATAG", "TCTAA", "TCTAT",
                                          "TTAGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 1);
        EXPECT_EQ(std::set<std::string>({ "AAGCT", "ACTAG", "ACTAT", "AGCTA", "AGCTT",
                                          "ATAGA", "ATAGC", "ATAGT", "CTAGC", "CTAGT",
                                          "CTATAG", "CTATC", "CTATGCGA", "GATAG",
                                          "GCTAG", "GCTAT", "TAGACTA", "TAGCT",
                                          "TAGTCTA", "TCGCATAG", "TCTAA", "TCTAT",
                                          "TTAGA" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 2);
        EXPECT_EQ(std::set<std::string>({ "AAGCT", "ACTAG", "ACTAT", "AGCTA", "AGCTT",
                                          "ATAGA", "ATAGC", "ATAGT", "CTAGC", "CTAGT",
                                          "CTATAG", "CTATC", "CTATGCGA", "GATAG",
                                          "GCTAG", "GCTAT", "TAGACTA", "TAGCT",
                                          "TAGTCTA", "TCGCATAG", "TCTAT" }), unitigs);

        unitigs.clear();
        graph->call_unitigs([&](const auto &unitig, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
            std::unique_lock<std::mutex> lock(seq_mutex);
            unitigs.insert(unitig);
        }, num_threads, 10);
        EXPECT_EQ(std::set<std::string>({ "AAGCT", "ACTAG", "ACTAT", "AGCTA", "AGCTT",
                                          "ATAGA", "ATAGC", "ATAGT", "CTAGC", "CTAGT",
                                          "CTATAG", "CTATC", "CTATGCGA", "GATAG",
                                          "GCTAG", "GCTAT", "TAGACTA", "TAGCT",
                                          "TAGTCTA", "TCGCATAG", "TCTAT" }), unitigs);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallKmersEmptyGraph) {
    for (size_t k = 2; k <= 30; ++k) {
        auto empty = build_graph<TypeParam>(k, {}, DeBruijnGraph::CANONICAL);
        size_t num_kmers = 0;
        empty->call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallKmersTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A') }, DeBruijnGraph::CANONICAL);

        ASSERT_EQ(2u, graph->num_nodes());

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'T') == sequence);
            num_kmers++;
        });
        EXPECT_EQ(2u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsCheckDegree) {
    for (size_t num_threads : { 1, 4 }) {
        std::vector<std::string> sequences {
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGG",
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGC",
            "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATAACTGGAAAA",
            "CCGAACTAGTGAAACTGCAACAGACATACGCTGCTCTGAACTCTAAGGC",
            "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
            "ATCGGAAGAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGC",
            "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAATTGGCTGAAAAGC",
            "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAATCGCTCCACCAA"
        };

        auto graph = build_graph_batch<TypeParam>(9, sequences, DeBruijnGraph::CANONICAL);

        std::multiset<std::string> unitigs {
            "AAATATTTACACAGAGTAGGAGACAAAT",
            "AAATATTTTTCCCTCCG",
            "AGACAAATCGCTCCACCAA",
            "AGACAAATTGGCTGAAAAGC",
            "AGTTCAGACGTGTGCTCTTCCGAT",
            "AGTTCAGAGCAGCGTATGTCTG",
            "ATCGGAAGAGCACACGTCTGAACT",
            "ATTTGTCTCCTACTCTGTGTAAATATTT",
            "ATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGG",
            "CAGACATAACTGGAAAA",
            "CAGACATACGCTGCTCTGAACT",
            "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATA",
            "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAAT",
            "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCG",
            "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
            "CCGAACTAGTGAAACTGCAACAGACATA",
            "CGCACCCTGTCTGACTACAACATCCAGAAAGAGTCCACCCTGCACCTGG",
            "CGCCCGAATCTGGCTTGGCGGAATATCTCTTTGACAAGCACACCCTGG",
            "CGGAGGGAAAAATATTT",
            "CTGAACTCCAGACACTAAGGCATCTCGTATGC",
            "CTGAACTCTAAGGC",
            "GAGTTCAGA",
            "GCATACGAGATGCCTTAGTGTCTGGAGTTCAG",
            "GCCTTAGAGTTCAG",
            "GCTTTTCAGCCAATTTGTCT",
            "TATGTCTGATTAAGAGTTAAAACTGAAGGTTTCATTTTGG",
            "TATGTCTGTTGCAGTTTCACTAGTTCGG",
            "TCTGAACTC",
            "TTGGTGGAGCGATTTGTCT",
            "TTTTCCAGTTATGTCTG"
        };

        std::mutex seq_mutex;
        std::multiset<std::string> obs_unitigs;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            std::unique_lock<std::mutex> lock(seq_mutex);
            obs_unitigs.insert(sequence);
        }, num_threads, 2);

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsIndegreeFirstNodeIsZero) {
    for (size_t num_threads : { 1, 4 }) {
        std::vector<std::string> sequences {
            "AGAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
            "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
            "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAATTAG"
        };

        auto graph = build_graph_batch<TypeParam>(31, sequences, DeBruijnGraph::CANONICAL);

        std::multiset<std::string> unitigs {
            "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
            "ATTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGCTGGTCAGGC",
            "CGCCACCACTCCCGGCTAATTTTGTATTTTTAGTAGAGACGGGGTTTC",
            "GAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
            "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAAT",
            "GTCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCT"
        };

        std::multiset<std::string> obs_unitigs;
        std::mutex seq_mutex;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            std::unique_lock<std::mutex> lock(seq_mutex);
            obs_unitigs.insert(sequence);
        }, num_threads, 2);

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

TYPED_TEST(DeBruijnGraphCanonicalTest, CallUnitigsCross) {
    for (size_t num_threads : { 1, 4 }) {
        // AATTT - ATTTT           TTTAA - TTAAA
        //               > TTTTA <
        // GGTTT - GTTTT           TTTAG - TTAGG

        // build graph from k-mers added in different order
        std::mutex seq_mutex;
        for (const auto &sequences : {
            std::vector<std::string>({ "AATTTTAAA",
                                       "GGTTTTAGG", }),
            std::vector<std::string>({ "GGTTTTAGG",
                                       "AATTTTAAA", }),
            std::vector<std::string>({ "TTTTAAA",
                                       "TTTTAGG",
                                       "AATTTTA",
                                       "GGTTTTA", }),
            std::vector<std::string>({ "AATTTTA",
                                       "GGTTTTA",
                                       "TTTTAAA",
                                       "TTTTAGG", }) }) {
            auto graph = build_graph_batch<TypeParam>(5, sequences, DeBruijnGraph::CANONICAL);

            std::multiset<std::string> unitigs {
                "AAAACC",
                "AAAATTTT",
                "CCTAAA",
                "GGTTTT",
                "TAAAA",
                "TTTAAA",
                "TTTAGG",
                "TTTTA",
            };

            for (size_t t = 0; t <= 2; ++t) {
                std::multiset<std::string> obs_unitigs;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    obs_unitigs.insert(sequence);
                }, num_threads, t);
                EXPECT_EQ(unitigs, obs_unitigs) << t;
            }

            std::multiset<std::string> long_unitigs {
                "AAAATTTT",
                "TAAAA",
                "TTTAAA",
                "TTTTA",
            };

            for (size_t t = 3; t <= 10; ++t) {
                std::multiset<std::string> obs_long_unitigs;
                graph->call_unitigs([&](const auto &sequence, const auto &path) {
                    ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    obs_long_unitigs.insert(sequence);
                }, num_threads, 3);
                EXPECT_EQ(long_unitigs, obs_long_unitigs);
            }
        }
    }
}

#endif // ! _PROTEIN_GRAPH

} // namespace
