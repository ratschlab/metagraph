#include "gtest/gtest.h"

#define private public
#define protected public

#include <set>

#include "../../test_helpers.hpp"
#include "test_dbg_helpers.hpp"
#include "graph/graph_extensions/node_weights.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";

const size_t kBitsPerCount = 8;

#if _PROTEIN_GRAPH
const size_t maxK = 12;
#else
const size_t maxK = 20;
#endif

TYPED_TEST_SUITE(DeBruijnGraphTest, GraphTypes);


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
    EXPECT_TRUE(check_graph_nodes(*graph));
    EXPECT_FALSE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
}

TYPED_TEST(DeBruijnGraphTest, SerializeEmpty) {
    {
        auto graph = build_graph<TypeParam>(maxK);
        ASSERT_EQ(0u, graph->num_nodes());
        graph->serialize(test_dump_basename);
        EXPECT_TRUE(check_graph_nodes(*graph));
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(0u, graph.num_nodes());
    EXPECT_EQ(maxK, graph.get_k());

    EXPECT_FALSE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_FALSE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_TRUE(check_graph_nodes(graph));
}

TYPED_TEST(DeBruijnGraphTest, Serialize) {
    {
        auto graph = build_graph<TypeParam>(maxK, {
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
        EXPECT_TRUE(check_graph_nodes(*graph));
    }

    TypeParam graph(2);

    ASSERT_TRUE(graph.load(test_dump_basename));

    EXPECT_EQ(maxK, graph.get_k());

    EXPECT_TRUE(graph.find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_FALSE(graph.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph.find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph.find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph.find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph.find("GCTAAAAATATATATATTAAAAAAACATG"));
    EXPECT_TRUE(check_graph_nodes(graph));
}

template <class Graph>
void test_graph_serialization(size_t k_max) {
    for (size_t k = 2; k <= k_max; ++k) {
        std::vector<std::string> data = { std::string(k, 'A'),
                                          std::string(k, 'C') + 'G', };
        {
            auto graph = build_graph<Graph>(k, data);

            EXPECT_TRUE(graph->find(std::string(k, 'A')));
            EXPECT_TRUE(graph->find(std::string(k, 'C')));
            EXPECT_TRUE(graph->find(std::string(k - 1, 'C') + 'G'));
            EXPECT_FALSE(graph->find(std::string(k, 'G')));

            graph->serialize(test_dump_basename);
        }
        {
            Graph graph(2);

            ASSERT_TRUE(graph.load(test_dump_basename)) << "k: " << k;

            EXPECT_EQ(k, graph.get_k());

            EXPECT_TRUE(graph.find(std::string(k, 'A')));
            EXPECT_TRUE(graph.find(std::string(k, 'C')));
            EXPECT_TRUE(graph.find(std::string(k - 1, 'C') + 'G'));
            EXPECT_FALSE(graph.find(std::string(k, 'G')));
        }
    }
}

#if _PROTEIN_GRAPH
TYPED_TEST(DeBruijnGraphTest, SerializeAnyK) {
    test_graph_serialization<TypeParam>(maxK);
}
#else
TYPED_TEST(DeBruijnGraphTest, SerializeAnyK) {
    test_graph_serialization<TypeParam>(31);
}
#endif

TEST(DBGHashString, SerializeAnyK) {
    test_graph_serialization<DBGHashString>(200);
}

TEST(DBGHashOrdered, SerializeAnyK) {
    test_graph_serialization<DBGHashOrdered>(256 / kmer::KmerExtractor2Bit::bits_per_char);
}

TEST(DBGHashFast, SerializeAnyK) {
    test_graph_serialization<DBGHashFast>(256 / kmer::KmerExtractor2Bit::bits_per_char);
}

TEST(DBGSuccinct, SerializeAnyK) {
    test_graph_serialization<DBGSuccinct>(256 / (kmer::KmerExtractor2Bit::bits_per_char + 1));
}

TYPED_TEST(DeBruijnGraphTest, InsertSequence) {
    auto graph = build_graph<TypeParam>(maxK, {
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "CATGTACTAGCTGATCGTAGCTAGCTAGC"
    });
    EXPECT_TRUE(check_graph_nodes(*graph));

    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
}

TYPED_TEST(DeBruijnGraphTest, Weighted) {
    for (size_t k = 2; k < 10; ++k) {
        auto sequences = {
            std::string(100, 'A'),
            std::string(k, 'G'),
            std::string(50, 'C')
        };
        auto graph = build_graph<TypeParam>(k, sequences);

        graph->add_extension(std::make_shared<NodeWeights>(graph->max_index() + 1, kBitsPerCount));
        auto &weights = *graph->template get_extension<NodeWeights>();

        for (const auto &sequence : sequences) {
            graph->map_to_nodes(sequence, [&](auto node) { weights.add_weight(node, 1); });
        }

        auto node_idx = graph->kmer_to_node(std::string(k, 'A'));
        EXPECT_EQ(100u - k + 1, weights[node_idx]);

        node_idx = graph->kmer_to_node(std::string(k, 'C'));
        EXPECT_EQ(50u - k + 1, weights[node_idx]);

        node_idx = graph->kmer_to_node(std::string(k, 'G'));
        EXPECT_EQ(1u, weights[node_idx]);

        //TODO should throw if weights extension is present
        //bit_vector_dyn nodes_inserted(graph->max_index() + 1, 0);
        //graph->add_sequence(std::string(25, 'T'), &nodes_inserted);

        EXPECT_TRUE(check_graph_nodes(*graph));
    }
}

TYPED_TEST(DeBruijnGraphTest, ReverseComplement) {
    auto graph1 = build_graph<TypeParam>(maxK, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA" });
    auto graph2 = build_graph<TypeParam>(maxK, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                                 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT" });

    EXPECT_EQ(graph1->num_nodes() * 2, graph2->num_nodes());

    auto graph = build_graph<TypeParam>(maxK, { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                                                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                                                "CATGTACTAGCTGATCGTAGCTAGCTAGC" });
    EXPECT_TRUE(check_graph_nodes(*graph));
    EXPECT_TRUE(graph->find("AAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    EXPECT_TRUE(graph->find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    EXPECT_TRUE(graph->find("CATGTACTAGCTGATCGTAGCTAGCTAGC"));
    EXPECT_FALSE(graph->find("GCTAGCTAGCTACGATCAGCTAGTACATG"));
    EXPECT_FALSE(graph->find("CATGTTTTTTTAATATATATATTTTTAGC"));
    EXPECT_FALSE(graph->find("GCTAAAAATATATATATTAAAAAAACATG"));
}

TYPED_TEST(DeBruijnGraphTest, CheckGraph) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGT", DBGMode::NORMAL, true));
}

TYPED_TEST(DeBruijnGraphTest, CheckGraphInputWithN) {
    EXPECT_TRUE(check_graph<TypeParam>("ACGTN", DBGMode::NORMAL, false));
    EXPECT_EQ(TypeParam(3).alphabet().find('N') != std::string::npos,
              check_graph<TypeParam>("ACGTN", DBGMode::NORMAL, true));
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
    if (TypeParam(2).alphabet().find('N') == std::string::npos) {
        EXPECT_EQ(0u, build_graph<TypeParam>(6, sequences)->num_nodes());
        EXPECT_EQ(0u, build_graph_batch<TypeParam>(6, sequences)->num_nodes());
    } else {
        EXPECT_EQ(1u, build_graph<TypeParam>(6, sequences)->num_nodes());
        EXPECT_EQ(1u, build_graph_batch<TypeParam>(6, sequences)->num_nodes());
    }
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

TYPED_TEST(DeBruijnGraphTest, CallKmersEmptyGraph) {
#if _PROTEIN_GRAPH
    for (size_t k = 2; k <= maxK; ++k) {
#else
    for (size_t k = 2; k <= 30; ++k) {
#endif
        auto empty = build_graph<TypeParam>(k);

        EXPECT_TRUE(check_graph_nodes(*empty));

        size_t num_kmers = 0;
        empty->call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTwoLoops) {
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A') });

        ASSERT_EQ(1u, graph->num_nodes());
        EXPECT_TRUE(check_graph_nodes(*graph));

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence);
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersFourLoops) {
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A'),
                                                 std::string(100, 'G'),
                                                 std::string(100, 'C') });
        ASSERT_EQ(3u, graph->num_nodes());
        EXPECT_TRUE(check_graph_nodes(*graph));

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
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph_batch<TypeParam>(k, { std::string(100, 'A'),
                                                       std::string(100, 'G'),
                                                       std::string(100, 'C') });
        EXPECT_TRUE(check_graph_nodes(*graph));

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
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(k - 1, 'C')
        });
        EXPECT_TRUE(check_graph_nodes(*graph));

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathACA) {
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, {
            std::string(100, 'A') + std::string(k - 1, 'C') + std::string(100, 'A')
        });
        EXPECT_TRUE(check_graph_nodes(*graph));

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k - 1, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathDisconnected) {
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A'),
                                                 std::string(100, 'T') });
        EXPECT_TRUE(check_graph_nodes(*graph));

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, CallKmersTestPathDisconnected2) {
    for (size_t k = 2; k <= maxK; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'G'),
                                                 std::string(k - 1, 'A') + "T" });
        EXPECT_TRUE(check_graph_nodes(*graph));

        size_t num_kmers = 0;
        graph->call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers);
    }
}

TYPED_TEST(DeBruijnGraphTest, call_source_nodes) {
    const std::vector<std::string> sequences {
        "ATGCGATCGATATGCGAGA",
        "ATGCGATCGAGACTACGAG",
        "GTACGATAGACATGACGAG",
        "ACTGACGAGACACAGATGC",
        "AAAAAAAAAAAAAAAAAAA"
    };
    for (int k = 8; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, sequences);

        EXPECT_TRUE(check_graph_nodes(*graph));

        std::multiset<std::string> start_nodes {
            sequences[0].substr(0, k),
            sequences[2].substr(0, k),
            sequences[3].substr(0, k)
        };

        std::multiset<std::string> obs_start_nodes;
        graph->call_source_nodes([&](auto start) {
            obs_start_nodes.emplace(graph->get_node_sequence(start));
        });

        EXPECT_EQ(start_nodes, obs_start_nodes) << *graph;
    }
}

TYPED_TEST(DeBruijnGraphTest, get_node_sequence) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = "AGCT";

    auto graph = build_graph<TypeParam>(k, { reference });
    EXPECT_TRUE(check_graph_nodes(*graph));

    std::string mapped_query = "";
    graph->map_to_nodes(query, [&](DeBruijnGraph::node_index node) {
        mapped_query += graph->get_node_sequence(node);
    });

    EXPECT_EQ(query, mapped_query);
}

TYPED_TEST(DeBruijnGraphTest, CallStartNodes) {
    {
        std::vector<std::string> sequences = { "AAACACTAG",
                                               "AACGACATG" };

        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(2, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAACA", "AACGA" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAACAC", "AACGAC" }), nodes);
        }
    }
    {
        std::vector<std::string> sequences = { "AGACACTGA",
                                               "GACTACGTA",
                                               "ACTAACGTA" };

        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(2, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGA" }), nodes) << *graph;
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACA", "GACTA" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACAC", "GACTAC", "ACTAAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(7, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACACT", "GACTACG", "ACTAACG" }), nodes);
        }
    }
    {
        std::vector<std::string> sequences = { "AGACACAGT",
                                               "GACTTGCAG",
                                               "ACTAGTCAG" };

        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(2, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACA", "GACTT", "ACTAG" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACAC", "GACTTG", "ACTAGT" }), nodes);
        }

    }
    {
        std::vector<std::string> sequences = { "AAACTCGTAGC",
                                               "AAATGCGTAGC" };

        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(2, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAAC", "AAAT" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            EXPECT_TRUE(check_graph_nodes(*graph));
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAACT", "AAATG" }), nodes);
        }

    }
}

} // namespace
