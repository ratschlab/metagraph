#include "gtest/gtest.h"

#define private public
#define protected public

#include "../../test_helpers.hpp"
#include "test_dbg_helpers.hpp"


namespace {

using namespace mtg;
using namespace mtg::test;

TYPED_TEST_SUITE(DeBruijnGraphTest, GraphTypes);


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

} // namespace
