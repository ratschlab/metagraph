#include "graph/representation/succinct/dbg_succinct_range.hpp"

#include "graph/representation/base/sequence_graph.hpp"

#include <gtest/gtest.h>

namespace {

using namespace mtg;
using namespace mtg::graph;


TEST(DBGSuccinctRange, CallNodesInRange) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);
    ASSERT_EQ(2u, suffix_matches.size());
    ASSERT_NE(DeBruijnGraph::npos, suffix_matches[0]);
    ASSERT_NE(DeBruijnGraph::npos, suffix_matches[1]);
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ("$$GG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$$G", range_graph.get_node_sequence(suffix_matches[1]));

    {
        std::multiset<DBGSuccinct::node_index> nodes;
        std::multiset<std::string> node_str;
        range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
            nodes.insert(node);
            auto ins = node_str.insert(graph->get_node_sequence(node));
            EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
        });

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
    {
        std::multiset<DBGSuccinct::node_index> nodes;
        std::multiset<std::string> node_str;
        range_graph.call_nodes_in_range(suffix_matches[1], [&](auto node) {
            nodes.insert(node);
            node_str.insert(graph->get_node_sequence(node));
        });

        std::multiset<DBGSuccinct::node_index> ref_nodes {
            graph->kmer_to_node("AGGG"),
            graph->kmer_to_node("CAGG"),
            graph->kmer_to_node("GGGG"),
            graph->kmer_to_node("CCAG")
        };

        std::multiset<std::string> ref_node_str {
            "AGGG",
            "CAGG",
            "GGGG",
            "CCAG"
        };

        EXPECT_EQ(ref_nodes, nodes) << *graph;
        EXPECT_EQ(ref_node_str, node_str);
    }
}

TEST(DBGSuccinctRange, CheckOffset) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "CAGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);
    auto suffix_matches = map_sequence_to_nodes(
        range_graph,
        { query.data(), std::min(size_t(query.size()), size_t(4)) }
    );

    ASSERT_EQ(4u, suffix_matches.size());
    EXPECT_EQ(1u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[2]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[3]));

    EXPECT_EQ("$CAG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$AG", range_graph.get_node_sequence(suffix_matches[1]));
    EXPECT_EQ("$$GC", range_graph.get_node_sequence(suffix_matches[2]));
    EXPECT_EQ("$$$C", range_graph.get_node_sequence(suffix_matches[3]));
}

TEST(DBGSuccinctRange, CallNodesWithSuffixK) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);
    ASSERT_EQ(1u, suffix_matches.size());
    EXPECT_EQ(0u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(map_sequence_to_nodes(*graph, query)[0], suffix_matches[0]);

    std::multiset<DBGSuccinct::node_index> nodes(suffix_matches.begin(),
                                                 suffix_matches.end());
    std::multiset<std::string> node_str {
        range_graph.get_node_sequence(suffix_matches[0])
    };

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinctRange, CallNodesWithSuffixK_v2) {
    size_t k = 4;


    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence("AGCCC");
    graph->add_sequence("CGCC");
    graph->add_sequence("GGCC");
    graph->add_sequence("TGCC");
    // graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    std::string query = "AGCC";

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);
    ASSERT_EQ(1u, suffix_matches.size());
    EXPECT_EQ(0u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(map_sequence_to_nodes(*graph, query)[0], suffix_matches[0]);
    EXPECT_EQ(query, range_graph.get_node_sequence(suffix_matches[0]));

    bool found = false;
    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        EXPECT_FALSE(found);
        EXPECT_EQ(suffix_matches[0], node);
        found = true;
    });
}

TEST(DBGSuccinctRange, CallNodesWithSuffixKEarlyCutoff) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(
        range_graph,
        { query.data(), std::min(size_t(query.size()), size_t(2)) }
    );

    ASSERT_EQ(2u, suffix_matches.size());
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[1]));

    EXPECT_EQ("$$GG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$$G", range_graph.get_node_sequence(suffix_matches[1]));

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query.substr(0, 2), ins->substr(ins->size() - 2));
    });

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

TEST(DBGSuccinctRange, CallNodesWithSuffixEarlyCutoffKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGGG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(
        range_graph,
        { query.data(), std::min(size_t(query.size()), size_t(3)) }
    );

    ASSERT_EQ(3u, suffix_matches.size());
    EXPECT_EQ(1u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[2]));

    EXPECT_EQ("$GGG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$GG", range_graph.get_node_sequence(suffix_matches[1]));
    EXPECT_EQ("$$$G", range_graph.get_node_sequence(suffix_matches[2]));

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;
    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query.substr(0, 3), ins->substr(ins->size() - 3));
    });

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

TEST(DBGSuccinctRange, CallNodesWithSuffixKMinusOne) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GCC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);

    ASSERT_EQ(3u, suffix_matches.size());
    EXPECT_EQ(1u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[2]));

    EXPECT_EQ("$GCC", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$CC", range_graph.get_node_sequence(suffix_matches[1]));
    EXPECT_EQ("$$$C", range_graph.get_node_sequence(suffix_matches[2]));

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
    });

    std::multiset<DBGSuccinct::node_index> ref_nodes {
        graph->kmer_to_node("GGCC")
    };

    std::multiset<std::string> ref_node_str {
        "GGCC"
    };

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str);
}

TEST(DBGSuccinctRange, CallNodesWithSuffixKMinusOneBeginning) {
    size_t k = 4;
    std::string reference = "GGCCCAGGGGTC";

    std::string query = "GGC";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);

    ASSERT_EQ(3u, suffix_matches.size());
    EXPECT_EQ(1u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[2]));

    EXPECT_EQ("$GGC", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$GC", range_graph.get_node_sequence(suffix_matches[1]));
    EXPECT_EQ("$$$C", range_graph.get_node_sequence(suffix_matches[2]));

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
    });

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinctRange, CallNodesWithSuffixMinusTwoBeginning) {
    size_t k = 4;
    std::string reference = "TGCCCAGGGGTC";

    std::string query = "TG";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);
    graph->mask_dummy_kmers(1, false);

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);

    ASSERT_EQ(2u, suffix_matches.size());
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[1]));

    EXPECT_EQ("$$TG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$$G", range_graph.get_node_sequence(suffix_matches[1]));

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query, ins->substr(ins->size() - query.size()));
    });

    EXPECT_TRUE(nodes.empty());
    EXPECT_TRUE(node_str.empty());
}

TEST(DBGSuccinctRange, CallNodesWithSuffixMultipleOut) {
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

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);

    ASSERT_EQ(1u, suffix_matches.size());
    EXPECT_EQ(map_sequence_to_nodes(*graph, query)[0], suffix_matches[0]);
    EXPECT_EQ(0u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ("TGC", range_graph.get_node_sequence(suffix_matches[0]));

    std::multiset<DBGSuccinct::node_index> ref_nodes { graph->kmer_to_node("TGC") };
    std::multiset<std::string> ref_node_str { "TGC" };

    std::multiset<DBGSuccinct::node_index> nodes;
    std::multiset<std::string> node_str;

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        size_t length = graph->get_k() - range_graph.get_offset(suffix_matches[0]);
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query, ins->substr(0, length));
    });

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

TEST(DBGSuccinctRange, CallNodesWithSuffixMultipleInOut) {
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

    DBGSuccinctRange range_graph(*graph);

    auto suffix_matches = map_sequence_to_nodes(range_graph, query);

    ASSERT_EQ(3u, suffix_matches.size());
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[0]));
    EXPECT_EQ(2u, range_graph.get_offset(suffix_matches[1]));
    EXPECT_EQ(3u, range_graph.get_offset(suffix_matches[2]));

    EXPECT_EQ("$$TG", range_graph.get_node_sequence(suffix_matches[0]));
    EXPECT_EQ("$$GA", range_graph.get_node_sequence(suffix_matches[1]));
    EXPECT_EQ("$$$A", range_graph.get_node_sequence(suffix_matches[2]));

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

    range_graph.call_nodes_in_range(suffix_matches[0], [&](auto node) {
        size_t length = graph->get_k() - range_graph.get_offset(suffix_matches[0]);
        nodes.insert(node);
        auto ins = node_str.insert(graph->get_node_sequence(node));
        EXPECT_EQ(query.substr(0, length), ins->substr(ins->size() - length, length));
    });

    EXPECT_EQ(ref_nodes, nodes) << *graph;
    EXPECT_EQ(ref_node_str, node_str) << *graph;
}

} // namespace
