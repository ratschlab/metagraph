#include "gtest/gtest.h"

#define private public
#define protected public

#include <set>

#include "../../test_helpers.hpp"
#include "test_dbg_helpers.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

TYPED_TEST_SUITE(DeBruijnGraphTest, GraphTypes);


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
        EXPECT_EQ(npos, graph->traverse_back(it, 'G'));
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
        auto it = graph->kmer_to_node(std::string(k, 'A'));
        // AAA -> AAA
        // AAA -> AAC
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

        it = graph->traverse(it, 'C');
        // AAC -> ACC
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

        it = graph->kmer_to_node(std::string(k, 'C'));
        // CCC -> CCC
        // CCC -> CCG
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

        it = graph->traverse(it, 'G');
        // CCG -> CGG
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

        it = graph->kmer_to_node("C" + std::string(k - 1, 'G'));
        // CGG -> {}
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
        it = graph->traverse(it, 'G');
        ASSERT_EQ(DeBruijnGraph::npos, it);
        ASSERT_EQ(DeBruijnGraph::npos, graph->kmer_to_node(std::string(k, 'G')));
    }
}

TYPED_TEST(DeBruijnGraphTest, CallIncomingEdges) {
    for (size_t k = 3; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(k - 1, 'A')
                                               + std::string(100, 'C')
                                               + std::string(k - 1, 'G') });
        // GGG does not exist
        auto it = graph->kmer_to_node(std::string(k, 'G'));
        ASSERT_EQ(DeBruijnGraph::npos, it);

        it = graph->kmer_to_node('C' + std::string(k - 1, 'G'));
        // CCG <- CGG
        std::set<char> set = { 'C' };
        graph->call_incoming_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse_back(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        it = graph->kmer_to_node(std::string(k - 1, 'C') + 'G');
        // ACC <- CCG
        // CCC <- CCG
        set = { 'A', 'C' };
        graph->call_incoming_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse_back(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        it = graph->traverse_back(it, 'C');
        // CCC <- CCC
        // ACC <- CCC
        set = { 'C', 'A' };
        graph->call_incoming_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse_back(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        it = graph->traverse_back(it, 'A');
        // AAC <- ACC
        set = { 'A' };
        graph->call_incoming_kmers(it,
            [&](auto i, char c) {
                ASSERT_TRUE(set.count(c))
                    << k
                    << "\n" << c
                    << "\n" << graph->get_node_sequence(it)
                    << "\n" << graph->get_node_sequence(i);
                set.erase(c);
                EXPECT_EQ(i, graph->traverse_back(it, c));
            }
        );
        ASSERT_TRUE(set.empty());

        it = graph->kmer_to_node(std::string(k - 1, 'A') + 'C');
        // {} <- AAC
        set = {};
        graph->call_incoming_kmers(it,
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

        // AAA does not exist
        it = graph->traverse_back(it, 'A');
        ASSERT_EQ(DeBruijnGraph::npos, it);
        ASSERT_EQ(DeBruijnGraph::npos, graph->kmer_to_node(std::string(k, 'A')));
    }
}

TYPED_TEST(DeBruijnGraphTest, OutgoingAdjacent) {
    for (size_t k = 2; k < 11; ++k) {
        auto graph = build_graph<TypeParam>(k, { std::string(100, 'A')
                                               + std::string(100, 'C')
                                               + std::string(100, 'G') });
        std::vector<DeBruijnGraph::node_index> adjacent_nodes;

        // AA, AAAAA
        auto it = graph->kmer_to_node(std::string(k, 'A'));
        graph->map_to_nodes_sequentially(std::string(k, 'A'), [&](auto i) { EXPECT_EQ(it, i); });
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
        graph->map_to_nodes_sequentially(std::string(k, 'C'), [&](auto i) { EXPECT_EQ(it, i); });
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
        graph->map_to_nodes_sequentially(std::string(k, 'G'), [&](auto i) { EXPECT_EQ(it, i); });
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

} // namespace
