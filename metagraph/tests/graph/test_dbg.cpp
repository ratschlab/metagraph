#include "gtest/gtest.h"

#define private public
#define protected public

#include "../test_helpers.hpp"
#include "test_dbg_helpers.hpp"
#include "node_weights.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test_graph";

const size_t kBitsPerCount = 8;

TYPED_TEST_CASE(DeBruijnGraphTest, GraphTypes);
TYPED_TEST_CASE(StableDeBruijnGraphTest, StableGraphTypes);

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

TYPED_TEST(DeBruijnGraphTest, SerializeAnyK) {
    test_graph_serialization<TypeParam>(31);
}

TEST(DBGHashString, SerializeAnyK) {
    test_graph_serialization<DBGHashString>(200);
}

TEST(DBGHashOrdered, SerializeAnyK) {
    test_graph_serialization<DBGHashOrdered>(256 / KmerExtractor2Bit::bits_per_char);
}

TEST(DBGHashFast, SerializeAnyK) {
    test_graph_serialization<DBGHashFast>(256 / KmerExtractor2Bit::bits_per_char);
}

TEST(DBGSuccinct, SerializeAnyK) {
    test_graph_serialization<DBGSuccinct>(256 / (KmerExtractor2Bit::bits_per_char + 1));
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

TYPED_TEST(DeBruijnGraphTest, Weighted) {
    for (size_t k = 2; k < 10; ++k) {
        auto sequences = {
            std::string(100, 'A'),
            std::string(k, 'G'),
            std::string(50, 'C')
        };
        auto graph = build_graph<TypeParam>(k, sequences, false);

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
    }
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

TYPED_TEST(DeBruijnGraphTest, traverse_string) {
    for (size_t k = 2; k < 11; ++k) {
        std::string sequence = "AGCTTCGAAGGCCTT";
        auto graph = build_graph_batch<TypeParam>(k, { sequence });

        for (size_t i = 0; i + k <= sequence.size(); ++i) {
            auto cur_node = graph->kmer_to_node(sequence.substr(i, k));
            ASSERT_NE(DeBruijnGraph::npos, cur_node);

            std::string path;
            graph->traverse(
                cur_node,
                sequence.data() + i + k, sequence.data() + sequence.size(),
                [&](auto node) {
                    path += graph->get_node_sequence(node).back();
                }
            );
            EXPECT_EQ(std::string(sequence.begin() + i + k, sequence.end()), path);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, traverse_string_stop_when_no_edge) {
    size_t k = 4;
    std::string sequence = "AGGCCTGTTTG";
    auto graph = build_graph_batch<TypeParam>(k, { sequence });

    std::string query = "CCCTGTTTG";
    graph->traverse(
        graph->kmer_to_node("AGGC"),
        query.data() + 4,
        query.data() + query.size(),
        [&](auto node) {
            EXPECT_FALSE(true) << node << " " << graph->get_node_sequence(node);
        },
        []() { return false; }
    );
}

TYPED_TEST(DeBruijnGraphTest, CallPathsEmptyGraph) {
    for (size_t k = 2; k <= 10; ++k) {
        auto empty = build_graph<TypeParam>(k);
        std::vector<std::string> sequences;
        empty->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*empty, sequence));
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
        empty->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*empty, sequence));
            sequences.push_back(sequence);
        });
        ASSERT_EQ(0u, sequences.size());

        EXPECT_EQ(*empty, *build_graph<TypeParam>(k, sequences));
        EXPECT_EQ(*empty, *build_graph_batch<TypeParam>(k, sequences));
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsOneSelfLoop) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(1u, graph->num_nodes());
        ASSERT_EQ(1u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_sequences++;
        });
        size_t num_sequences_batch = 0;
        graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
            num_sequences_batch++;
        });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsOneSelfLoop) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(1u, graph->num_nodes());
        ASSERT_EQ(1u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_sequences++;
        });
        size_t num_sequences_batch = 0;
        graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
            num_sequences_batch++;
        });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsThreeSelfLoops) {
    for (size_t k = 2; k <= 20; ++k) {
        std::vector<std::string> sequences { std::string(100, 'A'),
                                             std::string(100, 'G'),
                                             std::string(100, 'C') };
        auto graph = build_graph<TypeParam>(k, sequences);
        auto graph_batch = build_graph_batch<TypeParam>(k, sequences);
        ASSERT_EQ(3u, graph->num_nodes());
        ASSERT_EQ(3u, graph_batch->num_nodes());

        size_t num_sequences = 0;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_sequences++;
        });
        size_t num_sequences_batch = 0;
        graph_batch->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
            num_sequences_batch++;
        });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsExtractsLongestOneLoop) {
    for (size_t k = 4; k < 14; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "GGGGGGGGGGGGG" };
        auto graph = build_graph<TypeParam>(k, sequences);

        std::vector<std::string> contigs;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            contigs.push_back(sequence);
        });

        EXPECT_EQ(2u, contigs.size());
        EXPECT_EQ(convert_to_set({ "ATGCAGTACTCAG", std::string(k, 'G') }),
                  convert_to_set(contigs)) << k;
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPathsExtractsLongestTwoLoops) {
    for (size_t k = 4; k < 14; ++k) {
        std::vector<std::string> sequences { "ATGCAGTACTCAG",
                                             "ATGCAGTACTGAG",
                                             "GGGGGGGGGGGGG" };
        auto graph = build_graph<TypeParam>(k, sequences);

        std::vector<std::string> contigs;
        graph->call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            contigs.push_back(sequence);
        });

        EXPECT_EQ(3u, contigs.size());
    }
}

TYPED_TEST(DeBruijnGraphTest, CallContigsUniqueKmers) {
    std::string sequence = "GCAAATAAC";
    auto graph = build_graph<TypeParam>(3, { sequence });

    size_t num_kmers = 0;
    graph->call_sequences([&](const auto &sequence, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
        num_kmers += sequence.size() - 2;
    });

    EXPECT_EQ(sequence.size() - 2, num_kmers);
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsUniqueKmersCycle) {
    size_t k = 4;
    std::string sequence = "AAACCCGGGTTTAA";
    auto graph = build_graph<TypeParam>(k, { sequence });

    size_t num_unitigs = 0;
    size_t num_kmers = 0;
    graph->call_unitigs([&](const auto &sequence, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
        num_unitigs++;
        num_kmers += sequence.size() - k + 1;
    });

    EXPECT_EQ(1u, num_unitigs);
    EXPECT_EQ(sequence.size() - k + 1, num_kmers);
}

TYPED_TEST(DeBruijnGraphTest, CallContigsUniqueKmersCycle) {
    size_t k = 4;
    std::string sequence = "AAACCCGGGTTTAAA";
    auto graph = build_graph<TypeParam>(k, { sequence });

    size_t num_contigs = 0;
    size_t num_kmers = 0;
    graph->call_sequences([&](const auto &sequence, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
        num_contigs++;
        num_kmers += sequence.size() - k + 1;
    });

    EXPECT_EQ(1u, num_contigs);
    EXPECT_EQ(sequence.size() - k + 1, num_kmers);
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
        graph->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
            num_sequences++;
        });
        size_t num_sequences_batch = 0;
        graph_batch->call_unitigs([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, map_sequence_to_nodes(*graph_batch, sequence));
            num_sequences_batch++;
        });

        EXPECT_EQ(graph->num_nodes(), num_sequences);
        EXPECT_EQ(graph_batch->num_nodes(), num_sequences_batch);
        EXPECT_EQ(graph->num_nodes(), graph_batch->num_nodes());
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallPaths) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences);

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            });
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigs) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences);

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            });
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallPathsFromCanonical) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences, true);

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            });
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, true);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallPathsFromCanonicalSingleKmerForm) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences, true);

            std::vector<std::string> reconst;

            graph->call_sequences([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            }, true);
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, true);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigsFromCanonical) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences, true);

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            });
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, true);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(StableDeBruijnGraphTest, CallUnitigsFromCanonicalSingleKmerForm) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences, true);

            std::vector<std::string> reconst;

            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                reconst.push_back(sequence);
            }, 0, true);
            auto reconstructed_graph = build_graph_batch<TypeParam>(k, reconst, true);

            EXPECT_EQ(*graph, *reconstructed_graph);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallPaths) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences);

            // in stable graphs the order of input sequences
            // does not change the order of k-mers and their indexes
            auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

            auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                k,
                [&](const auto &callback) {
                    graph->call_sequences([&](const auto &sequence, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                        callback(sequence);
                    });
                }
            );

            EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigs) {
    for (size_t k = 2; k <= 10; ++k) {
        for (const std::vector<std::string> &sequences
                : { std::vector<std::string>({ "AAACACTAG", "AACGACATG" }),
                    std::vector<std::string>({ "AGACACTGA", "GACTACGTA", "ACTAACGTA" }),
                    std::vector<std::string>({ "AGACACAGT", "GACTTGCAG", "ACTAGTCAG" }),
                    std::vector<std::string>({ "AAACTCGTAGC", "AAATGCGTAGC" }),
                    std::vector<std::string>({ "AAACT", "AAATG" }),
                    std::vector<std::string>({ "ATGCAGTACTCAG", "ATGCAGTAGTCAG", "GGGGGGGGGGGGG" }) }) {

            auto graph = build_graph_batch<TypeParam>(k, sequences);

            // in stable graphs the order of input sequences
            // does not change the order of k-mers and their indexes
            auto stable_graph = build_graph_batch<DBGSuccinct>(k, sequences);

            auto reconstructed_stable_graph = build_graph_iterative<DBGSuccinct>(
                k,
                [&](const auto &callback) {
                    graph->call_unitigs([&](const auto &sequence, const auto &path) {
                        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                        callback(sequence);
                    });
                }
            );

            EXPECT_EQ(*stable_graph, *reconstructed_stable_graph);
        }
    }
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsWithoutTips) {
    size_t k = 3;
    auto graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                             "TCTAAGC" });
    ASSERT_EQ(6u, graph->num_nodes());

    std::set<std::string> unitigs;
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    EXPECT_EQ(std::set<std::string>({ "ACT", "TCT", "CTAAGC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    EXPECT_EQ(std::set<std::string>({ "ACT", "TCT", "CTAAGC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    EXPECT_EQ(std::set<std::string>({ "CTAAGC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    EXPECT_EQ(std::set<std::string>({ "CTAAGC" }), unitigs);


    graph = build_graph<TypeParam>(k, { "ACTAAGC",
                                        "ACTAAGT" });
    ASSERT_EQ(6u, graph->num_nodes());

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    EXPECT_EQ(std::set<std::string>({ "ACTAAG", "AGC", "AGT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    EXPECT_EQ(std::set<std::string>({ "ACTAAG", "AGC", "AGT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    EXPECT_EQ(std::set<std::string>({ "ACTAAG" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    EXPECT_EQ(std::set<std::string>({ "ACTAAG" }), unitigs);


    graph = build_graph<TypeParam>(k, { "ACTAAGCCC",
                                        "AAAGC",
                                        "TAAGCA" });
    ASSERT_EQ(9u, graph->num_nodes());

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
    EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCC", "CCC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCCC" }), unitigs);
    EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCA", "GCC", "CCC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC", "GCC", "CCC" }), unitigs);
    EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCC", "CCC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 3);
    // EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAGC" }), unitigs);
    EXPECT_EQ(std::set<std::string>({ "ACTAA", "AAA", "AAGC", "GCC", "CCC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    // EXPECT_EQ(std::set<std::string>({ "AAGC" }), unitigs);
    EXPECT_EQ(std::set<std::string>({ "AAGC", "AAA", "ACTAA", "GCC", "CCC" }), unitigs);


    graph = build_graph<TypeParam>(k, { "ACGAAGCCT",
                                        "AAGC",
                                        "TAAGCA" });
    ASSERT_EQ(9u, graph->num_nodes());

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    EXPECT_EQ(std::set<std::string>({ "ACGAA", "TAA", "AAGC", "GCA", "GCCT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    EXPECT_EQ(std::set<std::string>({ "ACGAA", "TAA", "AAGC", "GCA", "GCCT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    EXPECT_EQ(std::set<std::string>({ "ACGAA", "AAGC", "GCCT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 3);
    EXPECT_EQ(std::set<std::string>({ "ACGAA", "AAGC" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    EXPECT_EQ(std::set<std::string>({ "AAGC" }), unitigs);


    graph = build_graph<TypeParam>(k, { "TCTAAGCCG",
                                        "CATAAGCCG",
                                        "CATAACCGA" });
    ASSERT_EQ(12u, graph->num_nodes());

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    EXPECT_EQ(std::set<std::string>({ "CATA", "TCTA", "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 3);
    EXPECT_EQ(std::set<std::string>({ "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    EXPECT_EQ(std::set<std::string>({ "TAA", "AAGCC", "AACC", "CCGA" }), unitigs);
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsWithoutTips2) {
    size_t k = 5;
    auto graph = build_graph<TypeParam>(k, { "ACTATAGCTAGTCTATGCGA",
                                             "ACTATAGCTAGTCTAA",
                                             "ACTATAGCTA",
                                             "ACTATAGCTT",
                                             "ACTATC", });
    ASSERT_EQ(19u, graph->num_nodes());
    std::set<std::string> unitigs;
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 0);
    EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTT", "AGCTAGTCTA", "TCTAA", "TCTAT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 1);
    EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTT", "AGCTAGTCTA", "TCTAA", "TCTAT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 2);
    EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT" }), unitigs);

    unitigs.clear();
    graph->call_unitigs([&](const auto &unitig, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, unitig));
        unitigs.insert(unitig);
    }, 10);
    EXPECT_EQ(std::set<std::string>({ "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT" }), unitigs);
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

TYPED_TEST(DeBruijnGraphTest, CallUnitigsCheckDegree) {
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

    auto graph = build_graph_batch<TypeParam>(9, sequences);

    std::multiset<std::string> unitigs {
        "AGACAAATCGCTCCACCAA",
        "AGACAAATTGGCTGAAAAGC",
        "ATCGGAAGAGCACACGTCTGAACT",
        "CAGACATAACTGGAAAA",
        "CAGACATACGCTGCTCTGAACT",
        "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATA",
        "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAAT",
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCG",
        "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
        "CCGAACTAGTGAAACTGCAACAGACATA",
        "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAAT",
        "CTGAACTCCAGACACTAAGGCATCTCGTATGC",
        "CTGAACTCTAAGGC",
        "TCTGAACTC"
    };

    std::multiset<std::string> obs_unitigs;
    graph->call_unitigs([&](const auto &sequence, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
        obs_unitigs.insert(sequence);
    }, 2);

    EXPECT_EQ(unitigs, obs_unitigs);
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsIndegreeFirstNodeIsZero) {
    std::vector<std::string> sequences {
        "AGAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAATTAG"
    };

    auto graph = build_graph_batch<TypeParam>(31, sequences);

    std::multiset<std::string> unitigs {
        "GAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAAT"
    };

    std::multiset<std::string> obs_unitigs;
    graph->call_unitigs([&](const auto &sequence, const auto &path) {
        ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
        obs_unitigs.insert(sequence);
    }, 2);

    EXPECT_EQ(unitigs, obs_unitigs);
}

TYPED_TEST(DeBruijnGraphTest, CallUnitigsCross) {
    // AATTT - ATTTT           TTTAA - TTAAA
    //               > TTTTA <
    // GGTTT - GTTTT           TTTAG - TTAGG

    // build graph from k-mers added in different order
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
        auto graph = build_graph_batch<TypeParam>(5, sequences);

        std::multiset<std::string> unitigs {
            "AATTTT",
            "GGTTTT",
            "TTTTA",
            "TTTAAA",
            "TTTAGG",
        };

        for (size_t t = 0; t <= 2; ++t) {
            std::multiset<std::string> obs_unitigs;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                obs_unitigs.insert(sequence);
            }, t);
            EXPECT_EQ(unitigs, obs_unitigs) << t;
        }

        std::multiset<std::string> long_unitigs {
            "TTTTA",
        };

        for (size_t t = 3; t <= 10; ++t) {
            std::multiset<std::string> obs_long_unitigs;
            graph->call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, map_sequence_to_nodes(*graph, sequence));
                obs_long_unitigs.insert(sequence);
            }, 3);
            EXPECT_EQ(long_unitigs, obs_long_unitigs);
        }
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
        graph->adjacent_outgoing_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
        ASSERT_EQ(2u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<SequenceGraph::node_index>{ it, graph->traverse(it, 'C') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        graph->adjacent_outgoing_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_outgoing_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_outgoing_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_outgoing_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_incoming_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{ it }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // AC, AAAAC
        it = graph->traverse(it, 'C');
        graph->adjacent_incoming_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
        ASSERT_EQ(1u, adjacent_nodes.size());
        EXPECT_EQ(
            convert_to_set(std::vector<DeBruijnGraph::node_index>{ graph->traverse_back(it, 'A') }),
            convert_to_set(adjacent_nodes)
        );
        adjacent_nodes.clear();

        // CC, CCCCC
        it = graph->kmer_to_node(std::string(k, 'C'));
        graph->adjacent_incoming_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_incoming_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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
        graph->adjacent_incoming_nodes(it, [&](auto i) { adjacent_nodes.push_back(i); });
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

TYPED_TEST(DeBruijnGraphTest, RankIncomingEdge) {
    for (size_t k = 2; k <= 20; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C')
                                               + std::string(100, 'G') });

        EXPECT_EQ(0u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node(std::string(k, 'A')),
                                         graph->kmer_to_node(std::string(k, 'A'))));
        EXPECT_EQ(0u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node(std::string(k, 'A')),
                                         graph->kmer_to_node(std::string(k - 1, 'A') + 'C')));
        EXPECT_EQ(0u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node('A' + std::string(k - 1, 'C')),
                                         graph->kmer_to_node(std::string(k, 'C'))));
        EXPECT_EQ(1u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node(std::string(k, 'C')),
                                         graph->kmer_to_node(std::string(k, 'C'))));
        EXPECT_EQ(0u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node('C' + std::string(k - 1, 'G')),
                                         graph->kmer_to_node(std::string(k, 'G'))));
        EXPECT_EQ(1u, incoming_edge_rank(*graph,
                                         graph->kmer_to_node(std::string(k, 'G')),
                                         graph->kmer_to_node(std::string(k, 'G'))));
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
        EXPECT_EQ(0ull, graph->outdegree(graph->kmer_to_node(std::string(k - 1, 'A') + 'C')));
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
        graph->call_nodes([&](auto i) {
            if (i == max_outdegree_node_index) {
                EXPECT_EQ(4ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(0ull, graph->outdegree(i));
            }
        });
    }
}

void check_degree_functions(const DeBruijnGraph &graph) {
    graph.call_nodes([&](auto i) {

        size_t outdegree = graph.outdegree(i);
        if (outdegree == 0) {
            ASSERT_FALSE(graph.has_single_outgoing(i));
            ASSERT_FALSE(graph.has_multiple_outgoing(i));
        } else if (outdegree == 1) {
            ASSERT_TRUE(graph.has_single_outgoing(i));
            ASSERT_FALSE(graph.has_multiple_outgoing(i));
        } else {
            ASSERT_FALSE(graph.has_single_outgoing(i));
            ASSERT_TRUE(graph.has_multiple_outgoing(i));
        }
    });

    graph.call_nodes([&](auto i) {

        size_t indegree = graph.indegree(i);

        if (indegree == 0) {
            ASSERT_TRUE(graph.has_no_incoming(i));
            ASSERT_FALSE(graph.has_single_incoming(i));
        } else if (indegree == 1) {
            ASSERT_FALSE(graph.has_no_incoming(i));
            ASSERT_TRUE(graph.has_single_incoming(i));
        } else {
            ASSERT_FALSE(graph.has_no_incoming(i));
            ASSERT_FALSE(graph.has_single_incoming(i));
        }

    });
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
        graph->call_nodes([&](auto i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->outdegree(i));
            } else {
                EXPECT_EQ(1ull, graph->outdegree(i));
            }
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, get_indegree_single_node) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A') + 'C' });

        EXPECT_EQ(1ull, graph->num_nodes());
        EXPECT_EQ(0ull, graph->indegree(graph->kmer_to_node(std::string(k - 1, 'A') + 'C')));

        check_degree_functions(*graph);
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
        graph->call_nodes([&](auto i) {
            if (i == max_indegree_node_index) {
                EXPECT_EQ(4ull, graph->indegree(i));
            } else {
                EXPECT_EQ(0ull, graph->indegree(i));
            }
        });

        check_degree_functions(*graph);
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
        graph->call_nodes([&](auto i) {
            if (i == loop_node_index) {
                EXPECT_EQ(2ull, graph->indegree(i));
            } else {
                EXPECT_EQ(1ull, graph->indegree(i));
            }
        });

        check_degree_functions(*graph);
    }
}

TYPED_TEST(DeBruijnGraphTest, indegree1) {
    auto graph = build_graph<TypeParam>(3, { "ACTAAGCCC",
                                             "AAAGC",
                                             "TAAGCA" });
    ASSERT_EQ(9u, graph->num_nodes());

    EXPECT_EQ(2u, graph->indegree(graph->kmer_to_node("CCC")));
    EXPECT_EQ(1u, graph->outdegree(graph->kmer_to_node("CCC")));

    EXPECT_EQ(2u, graph->indegree(graph->kmer_to_node("AAA")));
    EXPECT_EQ(2u, graph->outdegree(graph->kmer_to_node("AAA")));

    check_degree_functions(*graph);
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

        check_degree_functions(*graph);
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

        check_degree_functions(*graph);
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

        graph->call_nodes([&](auto node) {
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, [&](auto i) { incoming_nodes.push_back(i); });
            EXPECT_EQ(graph->indegree(node), incoming_nodes.size())
                << "adjacent_incoming_nodes and indegree are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
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

        graph->call_nodes([&](auto node) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            EXPECT_EQ(graph->indegree(node), num_incoming_edges)
                << "traverse_back and indegree are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
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

        graph->call_nodes([&](auto node) {
            size_t num_incoming_edges = 0;
            for (auto c : graph->alphabet()) {
                if (graph->traverse_back(node, c))
                    num_incoming_edges++;
            }
            std::vector<DeBruijnGraph::node_index> incoming_nodes;
            graph->adjacent_incoming_nodes(node, [&](auto i) { incoming_nodes.push_back(i); });
            EXPECT_EQ(num_incoming_edges, incoming_nodes.size())
                << "adjacent_incoming_nodes and traverse_back are inconsistent for node: " << node;
        });

        check_degree_functions(*graph);
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
    graph->call_nodes([&](auto i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    });

    EXPECT_EQ(0u, single_outgoing_counter);
}

TYPED_TEST(DeBruijnGraphTest, is_single_outgoing_for_multiple_valid_edges) {
    size_t k = 4;
    std::string reference = "AGGGGTC";

    auto graph = build_graph<TypeParam>(k, { reference });

    uint64_t single_outgoing_counter = 0;
    graph->call_nodes([&](auto i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    });

    EXPECT_EQ(1u, single_outgoing_counter);
}

TYPED_TEST(DeBruijnGraphTest, CallStartNodes) {
    {
        std::vector<std::string> sequences = { "AAACACTAG",
                                               "AACGACATG" };

        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(2, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAACA", "AACGA" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
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
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGA" }), nodes) << *graph;
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACA", "GACTA" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACAC", "GACTAC", "ACTAAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(7, sequences);
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
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGAC" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AGACA", "GACTT", "ACTAG" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(6, sequences);
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
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(3, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>(), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(4, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAAC", "AAAT" }), nodes);
        }
        {
            std::multiset<std::string> nodes;
            auto graph = build_graph_batch<TypeParam>(5, sequences);
            graph->call_source_nodes([&](const auto &node) {
                nodes.insert(graph->get_node_sequence(node));
            });
            EXPECT_EQ(std::multiset<std::string>({ "AAACT", "AAATG" }), nodes);
        }

    }
}
