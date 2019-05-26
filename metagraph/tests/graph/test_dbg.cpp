#include "gtest/gtest.h"

#include "test_dbg_helpers.hpp"
#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";


template <typename Graph>
class DeBruijnGraphTest : public ::testing::Test { };
typedef ::testing::Types<DBGBitmap,
                         DBGHashString,
                         DBGHashOrdered,
                         DBGSuccinct> GraphTypes;
TYPED_TEST_CASE(DeBruijnGraphTest, GraphTypes);

TYPED_TEST(DeBruijnGraphTest, GraphDefaultConstructor) {
    TypeParam *graph = nullptr;

    ASSERT_NO_THROW({ graph = new TypeParam(2); });
    ASSERT_NE(nullptr, graph);
    delete graph;
    ASSERT_NO_THROW(TypeParam(2));
}

TYPED_TEST(DeBruijnGraphTest, InitializeEmpty) {
    auto graph = build_graph<TypeParam>(2);

    EXPECT_EQ(0u, graph->num_nodes());
    EXPECT_FALSE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphTest, SerializeEmpty) {
    {
        auto graph = build_graph<TypeParam>(20);
        ASSERT_EQ(0u, graph->num_nodes());
        graph->serialize(test_dump_basename);
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(20u, graph.get_k());

    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphTest, Serialize) {
    {
        auto graph = build_graph<TypeParam>(20, {
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "CATGTACTAGCTGATCGTAGCTAGCTAGC"
        });

        EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        EXPECT_FALSE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
        EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
        EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
        EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
        EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));

        graph->serialize(test_dump_basename);
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(20u, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphTest, InsertSequence) {
    auto graph = build_graph<TypeParam>(20, {
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "CATGTACTAGCTGATCGTAGCTAGCTAGC"
    });

    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TYPED_TEST(DeBruijnGraphTest, ReverseComplement) {
    auto graph1 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" });
    auto graph2 = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                               "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT" });

    EXPECT_EQ(graph1->num_nodes() * 2, graph2->num_nodes());

    auto graph = build_graph<TypeParam>(20, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                              "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                                              "CATGTACTAGCTGATCGTAGCTAGCTAGC" });
    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphTest, CheckGraph) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGT", false));
}

TYPED_TEST(DeBruijnGraphTest, Alphabet) {
    for (size_t k = 2; k <= 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {});
        std::set<char> alphabet(graph->alphabet().begin(), graph->alphabet().end());
        EXPECT_TRUE(alphabet.count('A'));
        EXPECT_TRUE(alphabet.count('C'));
        EXPECT_TRUE(alphabet.count('G'));
        EXPECT_TRUE(alphabet.count('T'));
    }
}

TYPED_TEST(DeBruijnGraphTest, AddSequenceSimplePath) {
    for (size_t k = 2; k <= 10; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A') };
        EXPECT_EQ(1u, build_graph<TypeParam>(k, sequences)->num_nodes());
        EXPECT_EQ(1u, build_graph_batch<TypeParam>(k, sequences)->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, AddSequenceSimplePaths) {
    for (size_t k = 2; k <= 10; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A'),
                                             std::string(100, 'C') };
        EXPECT_EQ(2u, build_graph<TypeParam>(k, sequences)->num_nodes());
        EXPECT_EQ(2u, build_graph_batch<TypeParam>(k, sequences)->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, TestNonASCIIStrings) {
    std::vector<std::string> sequences { // cyrillic A and C
                                         "АСАСАСАСАСАСА",
                                         "плохая строка",
                                         "АСАСАСАСАСАСА" };
    EXPECT_EQ(1u, build_graph<TypeParam>(6, sequences)->num_nodes());
    EXPECT_EQ(1u, build_graph_batch<TypeParam>(6, sequences)->num_nodes());
}

TYPED_TEST(DeBruijnGraphTest, AddSequences) {
    {
        std::vector<std::string> sequences { "AAAC", "CAAC" };
        EXPECT_EQ(2u, build_graph<TypeParam>(4, sequences)->num_nodes());
        EXPECT_EQ(2u, build_graph_batch<TypeParam>(4, sequences)->num_nodes());
    }
    {
        std::vector<std::string> sequences { "AAAC", "CAAC", "GAAC" };
        EXPECT_EQ(3u, build_graph<TypeParam>(4, sequences)->num_nodes());
        EXPECT_EQ(3u, build_graph_batch<TypeParam>(4, sequences)->num_nodes());
    }
    {
        std::vector<std::string> sequences { "AAAC", "AACG" };
        EXPECT_EQ(2u, build_graph<TypeParam>(4, sequences)->num_nodes());
        EXPECT_EQ(2u, build_graph_batch<TypeParam>(4, sequences)->num_nodes());
    }
    {
        // TODO: add version with N at the ends of these
        std::vector<std::string> sequences { "AGACT", "GACTT", "ACTAT" };
        EXPECT_EQ(3u, build_graph<TypeParam>(5, sequences)->num_nodes());
        EXPECT_EQ(3u, build_graph_batch<TypeParam>(5, sequences)->num_nodes());
    }
    {
        std::vector<std::string> sequences { "AGAC", "GACT", "ACTA" };
        EXPECT_EQ(3u, build_graph<TypeParam>(4, sequences)->num_nodes());
        EXPECT_EQ(3u, build_graph_batch<TypeParam>(4, sequences)->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        auto empty = build_graph<TypeParam>(k);
        std::vector<std::string> sequences;
        empty->call_sequences([&](const auto &sequence) {
            sequences.push_back(sequence);
        });
        ASSERT_EQ(0u, sequences.size());

        EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences));
        EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences));
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        auto empty = build_graph<TypeParam>(k);
        std::vector<std::string> sequences;
        empty->call_unitigs([&](const auto &sequence) {
            sequences.push_back(sequence);
        });
        ASSERT_EQ(0u, sequences.size());

        EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences));
        EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences));
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(1u, graph->num_nodes());
        ASSERT_EQ(1u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_sequences([&](auto&) { num_sequences++; });
        size_t num_sequences_batch = 0;
        graph_batch->call_sequences([&](auto&) { num_sequences_batch++; });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(1u, graph->num_nodes());
        ASSERT_EQ(1u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_unitigs([&](auto&) { num_sequences++; });
        size_t num_sequences_batch = 0;
        graph_batch->call_unitigs([&](auto&) { num_sequences_batch++; });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A'),
                                             std::string(100, 'G'),
                                             std::string(100, 'C') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(3u, graph->num_nodes());
        ASSERT_EQ(3u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_sequences([&](auto&) { num_sequences++; });
        size_t num_sequences_batch = 0;
        graph_batch->call_sequences([&](auto&) { num_sequences_batch++; });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A'),
                                             std::string(100, 'G'),
                                             std::string(100, 'C') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(3u, graph->num_nodes());
        ASSERT_EQ(3u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_unitigs([&](auto&) { num_sequences++; });
        size_t num_sequences_batch = 0;
        graph_batch->call_unitigs([&](auto&) { num_sequences_batch++; });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPaths) {
    for (size_t k = 2; k <= 10; ++k) {
        {
            auto graph = build_graph<TypeParam>(k, { "AAACACTAG",
                                                     "AACGACATG" });

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AGACACTGA",
                                                     "GACTACGTA",
                                                     "ACTAACGTA" });

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AGACACAGT",
                                                     "GACTTGCAG",
                                                     "ACTAGTCAG" });

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AAACTCGTAGC",
                                                     "AAATGCGTAGC" });

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AAACT",
                                                     "AAATG" });

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigs) {
    for (size_t k = 2; k <= 10; ++k) {
        {
            auto graph = build_graph<TypeParam>(k, { "AAACACTAG",
                                                     "AACGACATG" });

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AGACACTGA",
                                                     "GACTACGTA",
                                                     "ACTAACGTA" });

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AGACACAGT",
                                                     "GACTTGCAG",
                                                     "ACTAGTCAG" });

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AAACTCGTAGC",
                                                     "AAATGCGTAGC" });

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
        {
            auto graph = build_graph<TypeParam>(k, { "AAACT",
                                                     "AAATG" });

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence) {
                reconst.push_back(sequence);
            });
            auto reconstructed = build_graph<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersEmptyGraph) {
    for (size_t k = 2; k <= 30; ++k) {
        auto empty = build_graph<TypeParam>(k);
        size_t num_kmers = 0;
        empty->call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTwoLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A') });

        ASSERT_EQ(1u, graph->num_nodes());

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence);
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersFourLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A'),
                                                 std::string(100, 'G'),
                                                 std::string(100, 'C') });
        ASSERT_EQ(3u, graph->num_nodes());

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersFourLoopsDynamic) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph_batch<TypeParam>(k, { std::string(100, 'A'),
                                                       std::string(100, 'G'),
                                                       std::string(100, 'C') });

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPath) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(k - 1, 'C')
        });

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathACA) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(k - 1, 'C') + std::string(100, 'A')
        });

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k - 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathDisconnected) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A'),
                                                 std::string(100, 'T') });

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathDisconnected2) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'G'),
                                                 std::string(k - 1, 'A') + "T" });

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, FindSequence1) {
    for (size_t k = 2; k <= 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A') });

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

        constexpr double kEps = std::numeric_limits<double>::epsilon();
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

TYPED_TEST(DeBruijnGraphTest, Traversals) {
    const auto npos = DeBruijnGraph::npos;

    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C')
        });

        auto it = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(npos, it);

        it = graph->kmer_to_node(std::string(k, 'G'));
        ASSERT_EQ(npos, it);

        it = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(npos, it);

        EXPECT_EQ(it, graph->traverse(it, 'A'));

        ASSERT_NE(npos, graph->traverse(it, 'C'));
        EXPECT_NE(it, graph->traverse(it, 'C'));

        EXPECT_EQ(it, graph->traverse_back(graph->traverse(it, 'C'), 'A'));

        EXPECT_EQ(npos, graph->traverse(it, 'G'));
        EXPECT_EQ(npos, graph->traverse_back(it + 1, 'G'));
    }
}

TYPED_TEST(DeBruijnGraphTest, Traversals2) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(100, 'C') + std::string(100, 'G')
        });

        auto it = graph->kmer_to_node(std::string(k, 'A'));
        auto it2 = graph->kmer_to_node(std::string(k - 1, 'A') + "C");

        ASSERT_EQ(graph->kmer_to_node(std::string(k, 'A')), it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it2, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it2, 'A'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DeBruijnGraph::npos, graph->traverse_back(it2, 'G'));
    }
}

TYPED_TEST(DeBruijnGraphTest, CallOutgoingEdges) {
    for (size_t k = 3; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C')
                                               + std::string(k - 1, 'G') });
        // AAA -> AAA
        // AAA -> AAC
        auto it = graph->kmer_to_node(std::string(k, 'A'));
        std::set<char> set = { 'A', 'C' };
        graph->call_outgoing_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        // AAC -> ACC
        it = graph->traverse(it, 'C');
        set = { 'C', };
        graph->call_outgoing_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        // CCC -> CCC
        // CCC -> CCG
        it = graph->kmer_to_node(std::string(k, 'C'));
        set = { 'C', 'G' };
        graph->call_outgoing_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        // CCG -> CGG
        it = graph->traverse(it, 'G');
        set = { 'G', };
        graph->call_outgoing_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        // CGG -> {}
        it = graph->kmer_to_node("C" + std::string(k - 1, 'G'));
        set = {};
        graph->call_outgoing_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse(it, c))
                    << graph->get_node_sequence(i) << '\n' << c;
            }
        );
        ASSERT_TRUE(set.empty());

        // GGG does not exist
        it = graph->kmer_to_node(std::string(k, 'G'));
        ASSERT_EQ(DeBruijnGraph::npos, it);
    }
}

TYPED_TEST(DeBruijnGraphTest, OutgoingAdjacent) {
    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C')
                                               + std::string(100, 'G') });
        std::vector<DeBruijnGraph::node_index> adjacent_nodes;

        auto map_to_nodes_sequentially = [&](const auto &seq, auto callback) {
            graph->map_to_nodes_sequentially(seq.begin(), seq.end(), callback);
        };

        // AA, AAAAA
        auto it = graph->kmer_to_node(std::string(k, 'A'));
        map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { EXPECT_EQ(it, i); });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<SequenceGraph::node_index>{ it, graph->traverse(it, 'C') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        auto outset = convert_to_set(std::vector<SequenceGraph::node_index>{ graph->traverse(it, 'C') });
        if (k == 2) {
            outset.insert(graph->traverse(it, 'G'));
            ASSERT_EQ(2u, adjacent_nodes.size());
        } else {
            ASSERT_EQ(1u, adjacent_nodes.size());
        }

        EXPECT_EQ(outset, convert_to_set(adjacent_nodes));
        adjacent_nodes.clear();

        // CC, CCCCC
        it = graph->kmer_to_node(std::string(k, 'C'));
        map_to_nodes_sequentially(std::string(k, 'C'), [&](auto i) { EXPECT_EQ(it, i); });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<SequenceGraph::node_index>{
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
            convert_to_set(std::vector<SequenceGraph::node_index>{ graph->traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GGGGG
        it = graph->kmer_to_node(std::string(k, 'G'));
        map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { EXPECT_EQ(it, i); });
        ASSERT_NE(DeBruijnGraph::npos, it);
        graph->adjacent_outgoing_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<SequenceGraph::node_index>{ graph->traverse(it, 'G') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TYPED_TEST(DeBruijnGraphTest, IncomingAdjacent) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C')
                                               + std::string(100, 'G') });

        std::vector<DeBruijnGraph::node_index> adjacent_nodes;

        // AA, AAAAA
        auto it = graph->kmer_to_node(std::string(k, 'A'));
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{ it }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{ graph->traverse_back(it, 'A') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CC, CCCCC
        it = graph->kmer_to_node(std::string(k, 'C'));
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{
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
            convert_to_set(std::vector<DeBruijnGraph::node_index>{
                graph->traverse_back(it, 'A'),
                graph->traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // GG, GGGGG
        it = graph->kmer_to_node(std::string(k, 'G'));
        graph->adjacent_incoming_nodes(it, &adjacent_nodes);
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{
                it,
                graph->traverse_back(it, 'C')
            }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();
    }
}

TYPED_TEST(DeBruijnGraphTest, map_to_nodes) {
    for (size_t k = 2; k <= 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C') });

        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 2, 'A')
                                        + std::string(2 * (k - 1), 'C');
        std::vector<DeBruijnGraph::node_index> expected_result {
            DeBruijnGraph::npos,
            DeBruijnGraph::npos
        };

        for (size_t i = 2; i + k <= sequence_to_map.size(); ++i) {
            graph->map_to_nodes(
                sequence_to_map.substr(i, k),
                [&](auto i) { expected_result.push_back(i);}
            );
        }

        std::vector<DeBruijnGraph::node_index> observed_result;
        graph->map_to_nodes(
            sequence_to_map,
            [&](const auto &index) { observed_result.emplace_back(index); }
        );
        EXPECT_EQ(expected_result, observed_result);

        size_t pos = 0;
        graph->map_to_nodes(sequence_to_map,
                            [&](auto i) { EXPECT_EQ(expected_result[pos++], i); });
    }
}

TYPED_TEST(DeBruijnGraphTest, get_outdegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A') + 'C' });
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->outdegree(1));
    }
}

TYPED_TEST(DeBruijnGraphTest, get_maximum_outdegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k - 1, 'A') + 'A',
            std::string(k - 1, 'A') + 'C',
            std::string(k - 1, 'A') + 'G',
            std::string(k - 1, 'A') + 'T'
        });

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

TYPED_TEST(DeBruijnGraphTest, get_outdegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k - 1, 'A') + std::string(k - 1, 'C') +
                std::string(k - 1, 'G') + std::string(k, 'T'),
            std::string(k, 'A')
        });

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

TYPED_TEST(DeBruijnGraphTest, get_indegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A') + 'C' });
        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->indegree(1));
    }
}

TYPED_TEST(DeBruijnGraphTest, get_maximum_indegree) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            'A' + std::string(k - 1, 'A'),
            'C' + std::string(k - 1, 'A'),
            'G' + std::string(k - 1, 'A'),
            'T' + std::string(k - 1, 'A')
        });

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

TYPED_TEST(DeBruijnGraphTest, get_indegree_loop) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k, 'T')
        });

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

TYPED_TEST(DeBruijnGraphTest, get_degree1) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k, 'T')
        });

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node(std::string(k, 'T'));
        ASSERT_NE(DeBruijnGraph::npos, node_T);
        EXPECT_EQ(1ull, graph->outdegree(node_T));
        EXPECT_EQ(2ull, graph->indegree(node_T));
    }
}

TYPED_TEST(DeBruijnGraphTest, get_degree2) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k - 1, 'T')
        });

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DeBruijnGraph::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node('G' + std::string(k - 1, 'T'));
        ASSERT_NE(DeBruijnGraph::npos, node_T);
        EXPECT_EQ(0ull, graph->outdegree(node_T));
        EXPECT_EQ(1ull, graph->indegree(node_T));
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_incoming_indegree) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        for (uint64_t node = 1; node <= graph->num_nodes(); node++) {
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, &incoming_nodes);
            EXPECT_EQ(graph->indegree(node), incoming_nodes.size())
                << "adjacent_incoming_nodes and indegree are inconsistent for node: " << node;
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_indegree_traverse_back) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        for (uint64_t node = 1; node <= graph->num_nodes(); node++) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            EXPECT_EQ(graph->indegree(node), num_incoming_edges)
                << "traverse_back and indegree are inconsistent for node: " << node;
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree_identity_traverse_back_incoming) {
    for (int k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            "ATGCGATCGATATGCGAGA",
            "ATGCGATCGAGACTACGAG",
            "GTACGATAGACATGACGAG",
            "ACTGACGAGACACAGATGC"
        });

        for (uint64_t node = 1; node <= graph->num_nodes(); node++) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, &incoming_nodes);
            EXPECT_EQ(num_incoming_edges, incoming_nodes.size())
                << "adjacent_incoming_nodes and traverse_back are inconsistent for node: " << node;
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, get_node_sequence) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = "AGCT";

    auto graph = build_graph<TypeParam>(k, { reference });

    std::string mapped_query = "";
    graph->map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
        mapped_query += graph->get_node_sequence(node);
    });

    EXPECT_EQ(query, mapped_query);
}

TYPED_TEST(DeBruijnGraphTest, is_single_outgoing_simple) {
    size_t k = 4;
    std::string reference = "CATC";

    auto graph = build_graph<TypeParam>(k, { reference });

    uint64_t single_outgoing_counter = 0;
    for (DeBruijnGraph::node_index i = 1; i <= graph->num_nodes(); ++i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    }

    EXPECT_EQ(0u, single_outgoing_counter);
}

TYPED_TEST(DeBruijnGraphTest, is_single_outgoing_for_multiple_valid_edges) {
    size_t k = 4;
    std::string reference = "AGGGGTC";

    auto graph = build_graph<TypeParam>(k, { reference });

    uint64_t single_outgoing_counter = 0;
    for (DeBruijnGraph::node_index i = 1; i <= graph->num_nodes(); ++i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    }

    EXPECT_EQ(1u, single_outgoing_counter);
}

TYPED_TEST(DeBruijnGraphTest, CallNodes) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(k, 'A')
                + std::string(k - 1, 'C')
                + std::string(k - 1, 'G')
                + std::string(k, 'T')
        });

        std::vector<DeBruijnGraph::node_index> nodes;
        graph->call_nodes([&](const auto &index) { nodes.push_back(index);});
        EXPECT_EQ(graph->num_nodes(), nodes.size());

        std::sort(nodes.begin(), nodes.end());

        size_t i = 1;
        for (const auto &index : nodes) {
            EXPECT_EQ(i++, index);
        }
    }
}
