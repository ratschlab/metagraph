#include "graph/representation/succinct/dbg_succinct.hpp"

#include "graph/representation/base/sequence_graph.hpp"

#include <gtest/gtest.h>


namespace {

using namespace mtg;
using namespace mtg::graph;

TEST(DBGSuccinct, get_degree_with_source_dummy) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGSuccinct>(k);

        graph->add_sequence(std::string(k, 'A')
                                + std::string(k - 1, 'C')
                                + std::string(k - 1, 'G')
                                + std::string(k, 'T'));

        // dummy source k-mer: '$$$$$'
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(1));
        EXPECT_EQ(1ull, graph->outdegree(1));
        EXPECT_EQ(1ull, graph->indegree(1));

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DBGSuccinct::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(2ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node(std::string(k, 'T'));
        ASSERT_NE(DBGSuccinct::npos, node_T);
        EXPECT_EQ(1ull, graph->outdegree(node_T));
        EXPECT_EQ(2ull, graph->indegree(node_T));


        // Now mask out all dummy k-mers

        graph->mask_dummy_kmers(1, false);
        // dummy source k-mer: '$$$$$'
        EXPECT_NE(std::string(k, '$'), graph->get_node_sequence(1));

        // 'AAAAA'
        node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DBGSuccinct::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        node_T = graph->kmer_to_node(std::string(k, 'T'));
        ASSERT_NE(DBGSuccinct::npos, node_T);
        EXPECT_EQ(1ull, graph->outdegree(node_T));
        EXPECT_EQ(2ull, graph->indegree(node_T));
    }
}

TEST(DBGSuccinct, get_degree_with_source_and_sink_dummy) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGSuccinct>(k);

        graph->add_sequence(std::string(k, 'A')
                                + std::string(k - 1, 'C')
                                + std::string(k - 1, 'G')
                                + std::string(k - 1, 'T'));

        // dummy source k-mer: '$$$$$'
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(1));
        EXPECT_EQ(1ull, graph->outdegree(1));
        EXPECT_EQ(1ull, graph->indegree(1));

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DBGSuccinct::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(2ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node('G' + std::string(k - 1, 'T'));
        ASSERT_NE(DBGSuccinct::npos, node_T);
        EXPECT_EQ(0ull, graph->outdegree(node_T));
        EXPECT_EQ(1ull, graph->indegree(node_T));


        // Now mask out all dummy k-mers

        graph->mask_dummy_kmers(1, false);
        // dummy source k-mer: '$$$$$'
        EXPECT_NE(std::string(k, '$'), graph->get_node_sequence(1));

        // 'AAAAA'
        node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DBGSuccinct::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(1ull, graph->indegree(node_A));

        node_T = graph->kmer_to_node('G' + std::string(k - 1, 'T'));
        ASSERT_NE(DBGSuccinct::npos, node_T);
        EXPECT_EQ(0ull, graph->outdegree(node_T));
        EXPECT_EQ(1ull, graph->indegree(node_T));
    }
}

TEST(DBGSuccinct, is_single_outgoing_simple) {
    size_t k = 4;
    std::string reference = "CATC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);

    uint64_t single_outgoing_counter = 0;
    for (DBGSuccinct::node_index i = 1; i <= graph->num_nodes(); ++i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    }

    // All nodes except the last dummy node should be single outgoing.
    EXPECT_EQ(reference.size(), single_outgoing_counter);
}

TEST(DBGSuccinct, is_single_outgoing_for_multiple_valid_edges) {
    size_t k = 4;
    std::string reference = "AGGGGTC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);

    uint64_t single_outgoing_counter = 0;
    for (DBGSuccinct::node_index i = 1; i <= graph->num_nodes(); ++i) {
        if (graph->outdegree(i) == 1)
            single_outgoing_counter++;
    }

    // All nodes except the source dummy, the sink dummy, 'AGGG',
    // and 'GGGG' have a single outgoing edge.
    EXPECT_EQ(reference.size() - 2, single_outgoing_counter);
}

TEST(DBGSuccinct, call_outgoing_kmers_source) {
    size_t k = 4;
    std::vector<std::string> sequences {
        "AATGG", "CCGAA"
    };

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);

    std::multiset<std::pair<DBGSuccinct::node_index, char>> ref {
        { 2, 'A' }, { 3, 'C' }
    };

    std::multiset<std::pair<DBGSuccinct::node_index, char>> obs;

    graph->call_outgoing_kmers(
        1,
        [&](auto node, char c) { obs.emplace(node, c); }
    );

    EXPECT_EQ(ref, obs);
}

} // namespace
