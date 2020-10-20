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

TEST(DBGSuccinct, CallNodesWithSuffix) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("AGGG"),
        graph->kmer_to_node("CAGG"),
        graph->kmer_to_node("GGGG")
    };

    std::multiset<std::string> ref_node_str {
        "AGGG",
        "CAGG",
        "GGGG"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithSuffixMinLength) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "CAGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(4)) },
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        },
        4
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithSuffixK) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithSuffixK_v2) {
    size_t k = 4;


    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence("AGCCC");
    graph->add_sequence("CGCC");
    graph->add_sequence("GGCC");
    graph->add_sequence("TGCC");
    // graph->mask_dummy_kmers(1, false);

    std::string query = "AGCC";

    size_t nodes_called = 0;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](DBGSuccinct::node_index /* i */, size_t length) {
            EXPECT_EQ(query.size(), length);
            nodes_called++;
        },
        k
    );

    EXPECT_EQ(1u, nodes_called) << *graph;
}

TEST(DBGSuccinct, CallNodesWithSuffixKEarlyCutoff) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(2)) },
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 2), ins->substr(ins->size() - 2));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("AGGG"),
        graph->kmer_to_node("CAGG"),
        graph->kmer_to_node("GGGG")
    };

    std::multiset<std::string> ref_node_str {
        "AGGG",
        "CAGG",
        "GGGG"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithSuffixKEarlyCutoffMaxNum) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(2)) },
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 2), ins->substr(ins->size() - 2));
        },
        4,
        [&]() {
            if (nodes.size() > 2) {
                nodes.clear();
                node_str.clear();
                return true;
            }

            return false;
        }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithSuffixEarlyCutoffKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGGG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(3)) },
        [&](auto node, auto length) {
            EXPECT_EQ(3u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 3), ins->substr(ins->size() - 3));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("AGGG"),
        graph->kmer_to_node("GGGG")
    };

    std::multiset<std::string> ref_node_str {
        "AGGG",
        "GGGG"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithSuffixKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithSuffixKMinusOneBeginning) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithSuffixMinusTwoBeginning) {
    size_t k = 4;
    std::string reference = "TGCCCAGGGGTC";

    std::string query = "TG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithSuffixMultipleOut) {
    size_t k = 3;
    std::vector<std::string> sequences {
        "GGGGGGATGTAG",
        "GGGGGGATGCCTAATTAA"
    };

    std::string query = "TGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> ref_nodes { graph->kmer_to_node("TGC") };
    std::multiset<std::string> ref_node_str { "TGC" };

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

TEST(DBGSuccinct, CallNodesWithSuffixMultipleInOut) {
    size_t k = 4;
    std::vector<std::string> sequences {
        "AAAAAAAAATGC",
        "GGGGGGGGATGG",
        "GGGGGGGGTTGC",
        "AAAAAAAATTGG"
    };

    std::string query = "TGA";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);
    graph->add_sequence(sequences[2]);
    graph->add_sequence(sequences[3]);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("AATG"),
        graph->kmer_to_node("GATG"),
        graph->kmer_to_node("GTTG"),
        graph->kmer_to_node("ATTG")
    };
    std::multiset<std::string> ref_node_str {
        "AATG",
        "GATG",
        "GTTG",
        "ATTG"
    };

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    graph->call_nodes_with_suffix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, length), ins->substr(ins->size() - length, length));
        }
    );

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

TEST(DBGSuccinct, CallNodesWithPrefix) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC"),
        graph->kmer_to_node("GGGG"),
        graph->kmer_to_node("GGGT"),
        graph->kmer_to_node("GGTC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC", "GGGG", "GGGT", "GGTC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixMinLength) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "CAGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(4)) },
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        },
        4
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithPrefixK) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixK_v2) {
    size_t k = 4;

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence("CAGCC");
    graph->add_sequence("AGCCC");
    graph->add_sequence("AGCCG");
    graph->add_sequence("AGCCT");
    // graph->mask_dummy_kmers(1, false);

    std::string query = "AGCC";

    size_t nodes_called = 0;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](DBGSuccinct::node_index /* i */, size_t length) {
            EXPECT_EQ(query.size(), length);
            nodes_called++;
        },
        k
    );

    EXPECT_EQ(1u, nodes_called) << *graph;
}

TEST(DBGSuccinct, CallNodesWithPrefixKEarlyCutoff) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(2)) },
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 2), ins->substr(0, 2));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC"),
        graph->kmer_to_node("GGGG"),
        graph->kmer_to_node("GGGT"),
        graph->kmer_to_node("GGTC"),
    };

    std::multiset<std::string> ref_node_str {
        "GGCC",
        "GGGG",
        "GGGT",
        "GGTC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixKEarlyCutoffMaxNum) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(2)) },
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 2), ins->substr(0, 2));
        },
        4,
        [&]() { return nodes.size() >= 2; }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithPrefixEarlyCutoffKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGGG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        { query.data(), std::min(size_t(query.size()), size_t(3)) },
        [&](auto node, auto length) {
            EXPECT_EQ(3u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, 3), ins->substr(0, 3));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGGT"),
        graph->kmer_to_node("GGGG")
    };

    std::multiset<std::string> ref_node_str {
        "GGGT",
        "GGGG"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GCCC")
    };

    std::multiset<std::string> ref_node_str {
        "GCCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixKMinusOneBeginning) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixMinusTwoBeginning) {
    size_t k = 4;
    std::string reference = "TGCCCAGGGGTC";

    std::string query = "TG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("TGCC")
    };

    std::multiset<std::string> ref_node_str {
        "TGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinct, CallNodesWithPrefixKMinusOneEnd) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GTC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithPrefixMinusTwoEnd) {
    size_t k = 4;
    std::string reference = "TGCCCAGGGGTC";

    std::string query = "TC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinct, CallNodesWithPrefixMultipleOut) {
    size_t k = 3;
    std::vector<std::string> sequences {
        "GGGGGGATGTAG",
        "GGGGGGATGCCTAATTAA"
    };

    std::string query = "TGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> ref_nodes { graph->kmer_to_node("TGC") };
    std::multiset<std::string> ref_node_str { "TGC" };

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(query.size(), length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(0, length));
        }
    );

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

TEST(DBGSuccinct, CallNodesWithPrefixMultipleInOutEnd) {
    size_t k = 4;
    std::vector<std::string> sequences {
        "AAAAAAAAATGC",
        "GGGGGGGGATGG",
        "GGGGGGGGTTGC",
        "AAAAAAAATTGG"
    };

    std::string query = "TGA";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);
    graph->add_sequence(sequences[2]);
    graph->add_sequence(sequences[3]);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, length), ins->substr(0, length));
        }
    );

    EXPECT_TRUE(nodes.empty()) << *graph;
    EXPECT_TRUE(node_str.empty()) << *graph;
}

TEST(DBGSuccinct, CallNodesWithPrefixMultipleInOut) {
    size_t k = 4;
    std::vector<std::string> sequences {
        "ATGCAAAAAAAA",
        "ATGGGGGGGGGG",
        "TTGCGGGGGGGG",
        "TTGGAAAAAAAA"
    };

    std::string query = "TGA";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(sequences[0]);
    graph->add_sequence(sequences[1]);
    graph->add_sequence(sequences[2]);
    graph->add_sequence(sequences[3]);
    graph->mask_dummy_kmers(1, false);

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("TGCA"),
        graph->kmer_to_node("TGGG"),
        graph->kmer_to_node("TGCG"),
        graph->kmer_to_node("TGGA")
    };
    std::multiset<std::string> ref_node_str {
        "TGCA",
        "TGGG",
        "TGCG",
        "TGGA"
    };

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    graph->call_nodes_with_prefix_matching_longest_prefix(
        query,
        [&](auto node, auto length) {
            EXPECT_EQ(2u, length);
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query.substr(0, length), ins->substr(0, length));
        }
    );

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

} // namespace
