#include "all/test_dbg_helpers.hpp"

#include "graph/graph_extensions/path_index.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::test {

using namespace mtg::graph;
using node_index = DeBruijnGraph::node_index;

template <typename Graph>
class PathIndexTest : public DeBruijnGraphTest<Graph> {};

typedef ::testing::Types<DBGSuccinctUnitigIndexed> ChainGraphTypes;

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
        auto path_index = graph->template get_extension_threadsafe<ColumnPathIndex>();
        ASSERT_NE(nullptr, path_index);

        auto nodes = map_to_nodes(*graph, reference);
        ASSERT_EQ(1u, nodes.size());

        node_index node = nodes[0];
        auto chain_info_all = path_index->get_chain_info({ node });
        ASSERT_EQ(1u, chain_info_all.size());
        auto &chain_info = chain_info_all[0].second;
        ASSERT_EQ(1u, chain_info.size());

        bool found = false;
        path_index->call_distances("", chain_info[0], chain_info[0], [&](int64_t dist) {
            found = true;
            EXPECT_EQ(0, dist);
        });
        EXPECT_TRUE(found);
    }
}

void test_traversal_distances(const DeBruijnGraph &graph,
                              std::vector<std::tuple<std::string::const_iterator,
                                                     std::string::const_iterator,
                                                     tsl::hopscotch_set<size_t>>>&& tests,
                              int64_t max_distance = std::numeric_limits<int64_t>::max(),
                              size_t max_steps = 100) {
    size_t k = graph.get_k();
    auto path_index = graph.get_extension_threadsafe<ColumnPathIndex>();
    ASSERT_NE(nullptr, path_index);
    for (auto &[i1, i2, exp_dists] : tests) {
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

        auto chain_info_all = path_index->get_chain_info({ nodes_1[0], nodes_2[0] });
        ASSERT_EQ(1u, chain_info_all.size());
        auto &chain_info = chain_info_all[0].second;
        ASSERT_EQ(2u, chain_info.size());

        tsl::hopscotch_set<size_t> extra_dists;
        tsl::hopscotch_set<size_t> found_dists;

        path_index->call_distances("", chain_info[0], chain_info[1], [&](int64_t dist) {
            ASSERT_FALSE(found_dists.count(dist))
                << (canonical ? "PRIMARY\t" : "BASIC\t")
                << std::string_view(&*i1, k) << "\t" << std::string_view(&*i2, k) << "\n"
                << fmt::format("{}\t{}", dist, fmt::join(found_dists.begin(), found_dists.end(), ","));

            found_dists.emplace(dist);

            if (exp_dists.count(dist)) {
                exp_dists.erase(dist);
            } else {
                extra_dists.emplace(dist);
            }
        }, max_distance, max_steps);

        ASSERT_TRUE(exp_dists.empty())
            << (canonical ? "PRIMARY\t" : "BASIC\t")
            << std::string_view(&*i1, k) << "\t" << std::string_view(&*i2, k) << "\n"
            << fmt::format("{}", fmt::join(exp_dists.begin(), exp_dists.end(), ","));
        ASSERT_TRUE(extra_dists.empty())
            << (canonical ? "PRIMARY\t" : "BASIC\t")
            << std::string_view(&*i1, k) << "\t" << std::string_view(&*i2, k) << "\n"
            << fmt::format("{}", fmt::join(extra_dists.begin(), extra_dists.end(), ","));
    }
}

void test_traversal_distances(const DeBruijnGraph &graph,
                              std::vector<std::tuple<std::string::const_iterator,
                                                     std::string::const_iterator,
                                                     size_t>>&& tests,
                              int64_t max_distance = std::numeric_limits<int64_t>::max(),
                              size_t max_steps = 100) {
    std::vector<std::tuple<std::string::const_iterator,
                           std::string::const_iterator,
                           tsl::hopscotch_set<size_t>>> long_tests;
    long_tests.reserve(tests.size());
    for (const auto &[begin, end, d] : tests) {
        auto &[i, j, d_s] = long_tests.emplace_back(begin, end, tsl::hopscotch_set<size_t>{});
        if (d != std::numeric_limits<size_t>::max())
            d_s.emplace(d);
    }

    test_traversal_distances(graph, std::move(long_tests), max_distance, max_steps);
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
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG""G""CGCATCATCTAGCTACGATCTA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG""C""CGCATCATCTAGCTACGATCTA";

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
            tests.emplace_back(reference_1.begin() + i, reference_2.begin(), std::numeric_limits<size_t>::max());
            tests.emplace_back(reference_1.end() - k, reference_2.begin() + i, std::numeric_limits<size_t>::max());

            tests.emplace_back(reference_1.begin(), reference_2.begin() + i, i);
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k, k + 1 - i);
        }

        test_traversal_distances(*graph, std::move(tests));
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
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG""G""CGCATCATCTAGCTACGATCTA""G""TGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG""C""CGCATCATCTAGCTACGATCTA""C""TGATCGATCGGAGTGACGTGAA";

    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        tests.emplace_back(reference_1.begin(), reference_2.end() - k, (k + 1) * 2);
        for (size_t i = 1; i <= k; ++i) {
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k, (k + 1) * 2 - i);
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k - 1, (k + 1) * 2 - i - 1);
        }

        test_traversal_distances(*graph, std::move(tests));
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
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG""G""CGCATCATCTAGCTACGATCTA""G""TGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG""C""CGCATCATCTAGCTACGATCTA""C""TGATCGATCGGAGTGACGTGAA";

    std::string reference_3 = "CGTGGCCCAGGCCCAGGCCCAG""""AGCTAGTCGATGCCTAGCTGGCGA""""TGATCGATCGGAGTGACGTGAA";

    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2, reference_3 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        for (size_t i = 1; i <= k; ++i) {
            tests.emplace_back(reference_1.begin() + i, reference_2.end() - k, (k + 1) * 2 - i);
            tests.emplace_back(reference_3.begin() + i, reference_2.end() - k, (k + 1) * 2 - i);

            tests.emplace_back(reference_3.begin() + i, reference_2.end() - k - 1, std::numeric_limits<size_t>::max());
        }

        test_traversal_distances(*graph, std::move(tests));
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
    std::string reference_1 = "CGTGGCCCAGGCCCAGGCCCAG""G""CGCATCATCTAGCTACGATCTA""G""TGATCGATCGGAGTGACGTGAA";
    std::string reference_2 = "CGTGGCCCAGGCCCAGGCCCAG""C";
    std::string reference_3 =                            "CGCATCATCTAGCTACGATCTA""C""TGATCGATCGGAGTGACGTGAA";
    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference_1, reference_2, reference_3 }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               size_t>> tests;
        tests.emplace_back(reference_1.begin(), reference_2.begin() + 1, 1);
        tests.emplace_back(reference_2.begin() + 1, reference_1.begin() + 1, std::numeric_limits<size_t>::max());
        tests.emplace_back(reference_1.begin(), reference_3.end() - k, (k + 1) * 2);
        tests.emplace_back(reference_1.begin() + 1, reference_3.end() - k, (k + 1) * 2 - 1);
        tests.emplace_back(reference_1.begin() + 1, reference_3.end() - k - 1, (k + 1) * 2 - 2);
        test_traversal_distances(*graph, std::move(tests));
    }
}

TYPED_TEST(PathIndexTest, simple_cycle) {
    size_t k = 11;
    std::string reference = "GACTGTAGCTAGACTGTAGCTA";
    for (DeBruijnGraph::Mode mode : MODES_TO_TEST) {
        auto graph = build_graph_batch<TypeParam>(k, { reference }, mode);
        std::vector<std::tuple<std::string::const_iterator,
                               std::string::const_iterator,
                               tsl::hopscotch_set<size_t>>> tests;
        tests.emplace_back(reference.begin(), reference.begin(), tsl::hopscotch_set<size_t>{
            0, 11, 22, 33
        });
        test_traversal_distances(*graph, std::move(tests), k * 3, 3);
    }
}

} // namespace mtg::test
