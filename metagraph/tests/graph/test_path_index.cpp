#include "all/test_dbg_helpers.hpp"

#include "graph/graph_extensions/path_index.hpp"

namespace mtg::test {

using namespace mtg::graph;
using node_index = DeBruijnGraph::node_index;
using TypeParam = DBGSuccinctUnitigIndexed;

TEST(PathIndex, single_unitig) {
    size_t k = 11;
    std::string reference = "ACGATGCGATG";
    auto graph = build_graph_batch<TypeParam>(k, { reference });
    auto path_index = graph->get_extension_threadsafe<IPathIndex>();
    ASSERT_NE(nullptr, path_index);

    auto nodes = map_to_nodes(*graph, reference);
    ASSERT_EQ(1u, nodes.size());

    node_index node = nodes[0];
    auto coords = path_index->get_coords({ node });
    ASSERT_EQ(1u, coords.size());
    ASSERT_EQ(1u, coords[0].size());

    EXPECT_EQ(1u, coords[0][0].first);
    EXPECT_EQ(0u, path_index->get_dist(coords[0][0].first, coords[0][0].first));
}

TEST(PathIndex, simple_superbubble_ends) {
    size_t k = 22;
    /*
     *    o
     *   / \
     *  o   o
     *   \ /
     *    o
     */
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCGCATCATCTAGCTACGATCTA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCGCATCATCTAGCTACGATCTA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    auto path_index = graph->get_extension_threadsafe<IPathIndex>();
    ASSERT_NE(nullptr, path_index);

    auto nodes_1 = map_to_nodes(*graph, reference_1.substr(0, k));
    ASSERT_EQ(1u, nodes_1.size());

    auto nodes_2 = map_to_nodes(*graph, reference_1.substr(reference_1.size() - k, k));
    ASSERT_EQ(1u, nodes_2.size());

    EXPECT_NE(nodes_1[0], nodes_2[0]);

    auto coords = path_index->get_coords({ nodes_1[0], nodes_2[0] });
    ASSERT_EQ(2u, coords.size());
    ASSERT_EQ(1u, coords[0].size());
    ASSERT_EQ(1u, coords[1].size());

    size_t source = coords[0][0].first;
    size_t terminus = coords[1][0].first;
    EXPECT_NE(source, terminus);
    EXPECT_NE(coords[0][0].second, coords[1][0].second);

    EXPECT_EQ(k + 1, path_index->get_dist(source, terminus));
}

void test_traversal_distances(const DeBruijnGraph &graph,
                              const std::vector<std::tuple<size_t, size_t, size_t>> &tests,
                              const std::string &reference_1,
                              const std::string &reference_2) {
    size_t k = graph.get_k();
    auto path_index = graph.get_extension_threadsafe<IPathIndex>();
    ASSERT_NE(nullptr, path_index);
    for (auto [i1, i2, exp_dist] : tests) {
        auto nodes_1 = map_to_nodes(graph, reference_1.substr(i1, k));
        ASSERT_EQ(1u, nodes_1.size());

        auto nodes_2 = map_to_nodes(graph, reference_2.substr(i2, k));
        ASSERT_EQ(1u, nodes_2.size());

        auto coords = path_index->get_coords({ nodes_1[0], nodes_2[0] });
        ASSERT_EQ(2u, coords.size());
        ASSERT_EQ(1u, coords[0].size());
        ASSERT_EQ(1u, coords[1].size());

        size_t internal1 = coords[0][0].first;
        size_t internal2 = coords[1][0].first;
        ASSERT_EQ(nodes_1[0] == nodes_2[0],
                  coords[0][0].second == coords[1][0].second);

        ASSERT_EQ(exp_dist, path_index->get_dist(internal1, internal2)) << i1 << "\t" << i2;
    }
}

TEST(PathIndex, simple_superbubble) {
    size_t k = 22;
    /*
     *    o
     *   / \
     *  o   o
     *   \ /
     *    o
     */
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCGCATCATCTAGCTACGATCTA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCGCATCATCTAGCTACGATCTA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });

    std::vector<std::tuple<size_t, size_t, size_t>> tests;
    tests.emplace_back(0, 0, 0);
    tests.emplace_back(reference_1.size() - k, reference_1.size() - k, 0);
    tests.emplace_back(0, reference_1.size() - k, k + 1);
    for (size_t i = 1; i <= k; ++i) {
        tests.emplace_back(i, i, std::numeric_limits<size_t>::max());
        tests.emplace_back(0, i, 1);
        tests.emplace_back(i, 0, std::numeric_limits<size_t>::max());
        tests.emplace_back(i, reference_2.size() - k, k);
        tests.emplace_back(reference_1.size() - k, i, std::numeric_limits<size_t>::max());
    }

    test_traversal_distances(*graph, tests, reference_1, reference_2);
}

TEST(PathIndex, double_superbubble) {
    size_t k = 22;
    /*
     *    o   o
     *   / \ / \
     *  o   o   o
     *   \ / \ /
     *    o   o
     */
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCGCATCATCTAGCTACGATCTAGTGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCGCATCATCTAGCTACGATCTACTGATCGATCGGAGTGACGTGAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 });
    std::vector<std::tuple<size_t, size_t, size_t>> tests;
    tests.emplace_back(0, reference_2.size() - k, (k + 1) * 2);
    for (size_t i = 1; i <= k; ++i) {
        tests.emplace_back(i, reference_2.size() - k, (k + 1) * 2 - 1);
    }

    test_traversal_distances(*graph, tests, reference_1, reference_2);
}

TEST(PathIndex, nested_superbubbles) {
    size_t k = 22;
    /*
     *    o   o
     *   / \ / \
     *  o   o   o
     *  |\ / \ /|
     *  | o   o |
     *  |       |
     *  ----o----
     */
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCGCATCATCTAGCTACGATCTAGTGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGCCGCATCATCTAGCTACGATCTACTGATCGATCGGAGTGACGTGAA";
    std::string reference_3 = "CGTGGCCCAGGCCCAGGCCCAGAGCTAGTCGATGCCTAGCTGGCGATGATCGATCGGAGTGACGTGAA";

    auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2, reference_3 });
    std::vector<std::tuple<size_t, size_t, size_t>> tests;
    for (size_t i = 1; i <= k; ++i) {
        tests.emplace_back(i, reference_2.size() - k, (k + 1) * 2 - 1);
        tests.emplace_back(i, reference_2.size() - k - 1, std::numeric_limits<size_t>::max());
    }

    test_traversal_distances(*graph, tests, reference_3, reference_2);

    test_traversal_distances(*graph,
        std::vector<std::tuple<size_t, size_t, size_t>>{
            { 1, reference_2.size() - k - 1, k + 1 }
        },
        reference_1, reference_2
    );
}

} // namespace mtg::test
