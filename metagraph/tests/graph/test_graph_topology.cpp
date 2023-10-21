#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"

#include "graph/graph_extensions/graph_topology.hpp"

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::test;

template <typename Graph>
class DBGTopologyTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGSuccinctTopology> TopologyGraphTypes;

TYPED_TEST_SUITE(DBGTopologyTest, TopologyGraphTypes);

void test_get_coords(const DeBruijnGraph &graph,
                     const std::vector<std::tuple<std::string, size_t, size_t>> &expected_results) {
    const auto *topology = graph.get_extension_threadsafe<GraphTopology>();
    ASSERT_TRUE(topology);

    for (const auto &[kmer, exp_unitig_id, exp_cluster_id] : expected_results) {
        auto nodes = map_to_nodes(graph, kmer);
        ASSERT_EQ(1u, nodes.size());
        ASSERT_NE(DeBruijnGraph::npos, nodes[0]);

        auto coords = topology->get_coords(nodes);
        ASSERT_EQ(1u, coords.size());
        ASSERT_EQ(1u, coords[0].size());
        auto &[coord, topo] = coords[0][0];
        EXPECT_EQ(exp_unitig_id, topo.first) << kmer;
        EXPECT_EQ(exp_cluster_id, topo.second) << kmer;
    }
}

void test_unclustered(const DeBruijnGraph &graph,
                      const std::vector<std::string> &test_kmers) {
    const auto *topology = graph.get_extension_threadsafe<GraphTopology>();
    ASSERT_TRUE(topology);

    for (const auto &kmer : test_kmers) {
        auto nodes = map_to_nodes(graph, kmer);
        ASSERT_EQ(1u, nodes.size());
        ASSERT_NE(DeBruijnGraph::npos, nodes[0]);

        auto coords = topology->get_coords(nodes);
        ASSERT_EQ(1u, coords.size());
        ASSERT_EQ(1u, coords[0].size());
        auto &[coord, topo] = coords[0][0];
        EXPECT_EQ(topo.first, topo.second) << kmer;
    }
}

TYPED_TEST(DBGTopologyTest, superbubble) {
    size_t k = 7;

    std::string reference1 = "ATGTGATCATGTAGC";
    std::string reference2 = "ATGTGATGATGTAGC";

    auto graph = build_graph<TypeParam>(k, { reference1, reference2 });

    std::vector<std::tuple<std::string, size_t, size_t>> expected_results {
        { "ATGTGAT", 0, 0 }, { "ATGTAGC", 3, 0 }
    };

    test_get_coords(*graph, expected_results);
}

TYPED_TEST(DBGTopologyTest, superbubble_chain) {
    size_t k = 7;

    std::string reference1 = "ATGTGATCATGTAGCAGAGCTAG";
    std::string reference2 = "ATGTGATGATGTAGCTGAGCTAG";

    auto graph = build_graph<TypeParam>(k, { reference1, reference2 });

    std::vector<std::tuple<std::string, size_t, size_t>> expected_results {
        { "ATGTGAT", 0, 0 }, { "ATGTAGC", 3, 0 }, { "GAGCTAG", 6, 0 }
    };

    test_get_coords(*graph, expected_results);
}

TYPED_TEST(DBGTopologyTest, superbubble_chain_cycle_nochain) {
    size_t k = 7;

    std::string reference1 = "ATGTGATCATGTAGCGATGTGAT";
    std::string reference2 = "ATGTGATGATGTAGCCATGTGAT";

    auto graph = build_graph<TypeParam>(k, { reference1, reference2 });

    std::vector<std::string> test_kmers;
    for (size_t i = 0; i < reference1.size() - k + 1; ++i) {
        test_kmers.emplace_back(reference1.begin() + i, reference1.begin() + i + k);
    }
    for (size_t i = 0; i < reference2.size() - k + 1; ++i) {
        test_kmers.emplace_back(reference2.begin() + i, reference2.begin() + i + k);
    }
    test_unclustered(*graph, test_kmers);
}

TYPED_TEST(DBGTopologyTest, superbubble_chain_cycle_maybechain) {
    size_t k = 7;

    std::string reference1 = "TCCCCCCAATGTGATCATGTAGCGATGTGAT";
    std::string reference2 = "TCCCCCCTATGTGATGATGTAGCCATGTGAT";

    auto graph = build_graph<TypeParam>(k, { reference1, reference2 });

    std::vector<std::tuple<std::string, size_t, size_t>> expected_results {
        { "TCCCCCC", 0, 0 }, { "ATGTGAT", 3, 3 }, { "ATGTAGC", 6, 6 }
    };

    test_get_coords(*graph, expected_results);

    std::vector<std::string> test_kmers;
    for (size_t i = 0; i < reference1.size() - k + 1; ++i) {
        test_kmers.emplace_back(reference1.begin() + i, reference1.begin() + i + k);
    }
    for (size_t i = 0; i < reference2.size() - k + 1; ++i) {
        test_kmers.emplace_back(reference2.begin() + i, reference2.begin() + i + k);
    }
    test_unclustered(*graph, test_kmers);
}

} // namespace