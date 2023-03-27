#include "all/test_dbg_helpers.hpp"

#include "graph/graph_extensions/path_index.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::test {

using namespace mtg::graph;
using node_index = DeBruijnGraph::node_index;

template <typename Graph>
class PathIndexTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGSuccinctUnitigIndexed,
                         DBGSuccinctPathIndexed> ChainGraphTypes;

TYPED_TEST_SUITE(PathIndexTest, ChainGraphTypes);

static const std::vector<DeBruijnGraph::Mode> MODES_TO_TEST {
    DeBruijnGraph::BASIC
};

TYPED_TEST(PathIndexTest, selfloop) {
    size_t k = 4;
    std::string reference = "AAAA";
    auto graph = build_graph_batch<TypeParam>(k, { reference });
    EXPECT_TRUE(graph);
}

TYPED_TEST(PathIndexTest, padded_selfloop) {
    size_t k = 4;
    std::string reference = "TTTAAAACCC";
    auto graph = build_graph_batch<TypeParam>(k, { reference });
    EXPECT_TRUE(graph);
}

TYPED_TEST(PathIndexTest, single_unitig) {
    size_t k = 11;
    std::string reference = "ACGATGCGATG";
    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference }, mode);
        auto path_index = graph->template get_extension_threadsafe<IPathIndex>();
        ASSERT_NE(nullptr, path_index);

        auto nodes = map_to_nodes(*graph, reference);
        ASSERT_EQ(1u, nodes.size());

        node_index node = nodes[0];
        auto coords = path_index->get_coords({ node });
        ASSERT_EQ(1u, coords.size());
        ASSERT_EQ(1u, coords[0].size());

        bool found = false;
        path_index->call_dists(coords[0][0].first, coords[0][0].first, [&](size_t dist) {
            found = true;
            EXPECT_EQ(0u, dist);
        });
        EXPECT_TRUE(found);
    }
}

template <class TypeParam>
void test_traversal_distances(const DeBruijnGraph &graph,
                              const std::vector<std::tuple<std::string::const_iterator,
                                                           std::string::const_iterator, size_t>> &tests) {
    size_t k = graph.get_k();
    auto path_index = graph.get_extension_threadsafe<IPathIndex>();
    ASSERT_NE(nullptr, path_index);
    for (auto [i1, i2, exp_dist] : tests) {
        auto nodes_1 = map_to_nodes_sequentially(graph, std::string_view(&*i1, k));
        ASSERT_EQ(1u, nodes_1.size());

        auto nodes_2 = map_to_nodes_sequentially(graph, std::string_view(&*i2, k));
        ASSERT_EQ(1u, nodes_2.size());

        const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
        if (canonical) {
            ASSERT_EQ(canonical->get_base_node(nodes_1[0]) == nodes_1[0],
                      canonical->get_base_node(nodes_2[0]) == nodes_2[0]);
            if (canonical->get_base_node(nodes_1[0]) != nodes_1[0]) {
                std::swap(nodes_1[0], nodes_2[0]);
                nodes_1[0] = canonical->get_base_node(nodes_1[0]);
                nodes_2[0] = canonical->get_base_node(nodes_2[0]);
            }
        }

        auto coords = path_index->get_coords({ nodes_1[0], nodes_2[0] });
        ASSERT_EQ(2u, coords.size());

        if constexpr(std::is_same_v<TypeParam, DBGSuccinctUnitigIndexed>) {
            ASSERT_EQ(1u, coords[0].size());
            ASSERT_EQ(1u, coords[1].size());
        } else {
            ASSERT_LE(1u, coords[0].size());
            ASSERT_LE(1u, coords[1].size());
        }

        size_t internal1 = coords[0][0].first;
        size_t internal2 = coords[1][0].first;
        ASSERT_EQ(nodes_1[0] == nodes_2[0],
                  coords[0][0].second == coords[1][0].second);

        bool found = false;
        path_index->call_dists(internal1, internal2, [&](size_t dist) {
            found = true;
            ASSERT_EQ(exp_dist, dist)
                << (canonical ? "PRIMARY\t" : "BASIC\t")
                << std::string_view(&*i1, k) << "\t" << std::string_view(&*i2, k);
        }, 100);
        ASSERT_EQ(found, exp_dist != std::numeric_limits<size_t>::max())
                << (canonical ? "PRIMARY\t" : "BASIC\t")
                << std::string_view(&*i1, k) << "\t" << std::string_view(&*i2, k);
    }
}

TYPED_TEST(PathIndexTest, simple_superbubble) {
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

    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 }, mode);

        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        tests.emplace_back(reference_1.begin(), reference_1.begin(), 0);
        tests.emplace_back(reference_1.end() - k, reference_1.end() - k, 0);
        tests.emplace_back(reference_1.begin(), reference_1.end() - k, k + 1);
        for (size_t i = 1; i <= k; ++i) {
            tests.emplace_back(reference_1.begin() + i, reference_2.begin() + i, std::numeric_limits<size_t>::max());
            tests.emplace_back(reference_1.begin(), reference_2.begin() + i, 1);
            tests.emplace_back(reference_1.begin() + i, reference_2.begin(), std::numeric_limits<size_t>::max());
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k, k);
            tests.emplace_back(reference_1.end() - k, reference_2.begin() + i, std::numeric_limits<size_t>::max());
        }

        test_traversal_distances<TypeParam>(*graph, tests);
    }
}

TYPED_TEST(PathIndexTest, double_superbubble_chain) {
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

    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        tests.emplace_back(reference_1.begin(), reference_2.end() - k, (k + 1) * 2);
        for (size_t i = 1; i <= k; ++i) {
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k, (k + 1) * 2 - 1);
        }

        test_traversal_distances<TypeParam>(*graph, tests);
    }
}

TYPED_TEST(PathIndexTest, nested_superbubbles) {
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

    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2, reference_3 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        for (size_t i = 1; i <= k; ++i) {
            tests.emplace_back(reference_3.begin() + i, reference_2.end() - k, (k + 1) * 2 - 1);
            tests.emplace_back(reference_3.begin() + i, reference_2.end() - k - 1, std::numeric_limits<size_t>::max());
        }

        // the lower long unitig forces this entire graph to be a superbubble
        tests.emplace_back(reference_1.begin() + 1, reference_2.end() - k - 1, std::numeric_limits<size_t>::max());

        test_traversal_distances<TypeParam>(*graph, tests);
    }
}

TYPED_TEST(PathIndexTest, disconnect) {
    size_t k = 22;
    /*
     *      o
     *     / \
     *  o-o   o
     *   \ \ /
     *    o o
     */
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAGGCGCATCATCTAGCTACGATCTAGTGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAGC";
    std::string reference_3 =                        "CGCATCATCTAGCTACGATCTACTGATCGATCGGAGTGACGTGAA";
    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2, reference_3 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        tests.emplace_back(reference_1.begin(), reference_2.begin() + 1, 1);
        tests.emplace_back(reference_2.begin() + 1, reference_1.begin() + 1, std::numeric_limits<size_t>::max());
        tests.emplace_back(reference_1.begin(), reference_3.end() - k, (k + 1) * 2);
        tests.emplace_back(reference_1.begin() + 1, reference_3.end() - k, (k + 1) * 2 - 1);
        tests.emplace_back(reference_1.begin() + 1, reference_3.end() - k - 1, (k + 1) * 2 - 1 - k);
        test_traversal_distances<TypeParam>(*graph, tests);
    }
}

TYPED_TEST(PathIndexTest, simple_cycle) {
    size_t k = 11;
    std::string reference = "GACTGTAGCTAGACTGTAGCTA";
    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
    }
}

} // namespace mtg::test
