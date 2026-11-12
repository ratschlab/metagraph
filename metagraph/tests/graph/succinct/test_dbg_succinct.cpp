#include "graph/representation/succinct/dbg_succinct.hpp"

#include "graph/representation/base/sequence_graph.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <set>


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
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(graph->select_node(1)));
        EXPECT_EQ(1ull, graph->outdegree(graph->select_node(1)));
        EXPECT_EQ(1ull, graph->indegree(graph->select_node(1)));

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
        EXPECT_NE(std::string(k, '$'), graph->get_node_sequence(graph->select_node(1)));

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
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(graph->select_node(1)));
        EXPECT_EQ(1ull, graph->outdegree(graph->select_node(1)));
        EXPECT_EQ(1ull, graph->indegree(graph->select_node(1)));

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
        EXPECT_NE(std::string(k, '$'), graph->get_node_sequence(graph->select_node(1)));

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
        if (graph->outdegree(graph->select_node(i)) == 1)
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
        if (graph->outdegree(graph->select_node(i)) == 1)
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

    for (auto [node, c] : graph->outgoing_kmers(1)) {
        obs.emplace(node, c);
    }

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

TEST(DBGSuccinct, CallNodesSmallGraphManyThreads) {
    // Regression: call_nodes with max_index() < num_threads caused
    // block_size = 0, triggering SIGFPE in the OpenMP trip count computation.
    for (size_t k = 2; k < 5; ++k) {
        auto graph = std::make_unique<DBGSuccinct>(k);
        graph->add_sequence(std::string(k, 'A'));
        graph->mask_dummy_kmers(1, false);

        ASSERT_GT(graph->max_index(), 0u);

        size_t num_nodes = 0;
        graph->call_nodes([&](auto) { num_nodes++; },
                          []() { return false; },
                          graph->max_index() + 10);
        EXPECT_EQ(num_nodes, graph->num_nodes());
    }
}

TEST(DBGSuccinct, CallNodesEmptyGraphManyThreads) {
    auto graph = std::make_unique<DBGSuccinct>(3);
    graph->mask_dummy_kmers(1, false);

    size_t num_nodes = 0;
    graph->call_nodes([&](auto) { num_nodes++; },
                      []() { return false; },
                      8);
    EXPECT_EQ(0u, num_nodes);
}

#if ! _PROTEIN_GRAPH

// (k-mer, offset)
typedef std::set<std::pair<DBGSuccinct::node_index, size_t>> Contexts;

const std::chrono::milliseconds kNoTimeLimit(1'000'000);

std::string random_sequence(std::mt19937 &gen, size_t length, const std::string &alphabet) {
    std::string s(length, 'A');
    for (char &c : s) {
        c = alphabet[gen() % alphabet.size()];
    }
    return s;
}

// The graph without a mask. Sequences shorter than 2k put dummy source k-mers in range.
std::unique_ptr<DBGSuccinct> random_graph(std::mt19937 &gen, size_t k, DeBruijnGraph::Mode mode,
                                          size_t num_sequences = 6) {
    std::vector<std::string> sequences;
    for (size_t i = 0; i < num_sequences; ++i) {
        sequences.push_back(random_sequence(gen, k + gen() % (2 * k), "ACGT"));
    }
    if (mode == DeBruijnGraph::PRIMARY) {
        DBGSuccinct canonical(k, DeBruijnGraph::CANONICAL);
        for (const std::string &sequence : sequences) {
            canonical.add_sequence(sequence);
        }
        sequences.clear();
        canonical.call_sequences([&](const std::string &contig, const auto &) {
            sequences.push_back(contig);
        }, 1, true);
    }
    auto graph = std::make_unique<DBGSuccinct>(k, mode);
    for (const std::string &sequence : sequences) {
        graph->add_sequence(sequence);
    }
    return graph;
}

// scan the spelling of every edge, and skip the dummy k-mers
Contexts brute_force(const DBGSuccinct &graph, const boss::Pattern &pattern) {
    const boss::BOSS &boss = graph.get_boss();
    Contexts contexts;
    for (uint64_t i = 1; i <= boss.num_edges(); ++i) {
        std::string kmer = boss.get_node_str(i) + boss.decode(boss.get_W(i) % boss.alph_size);
        if (kmer.find('$') != std::string::npos)
            continue;

        for (size_t offset = 0; offset + pattern.size() <= kmer.size(); ++offset) {
            bool match = true;
            for (size_t j = 0; j < pattern.size(); ++j) {
                match &= (pattern[j] >> boss.encode(kmer[offset + j])) & 1;
            }
            if (match)
                contexts.emplace(i, offset);
        }
    }
    return contexts;
}

// check the list and the count of a pattern against the scan
void test_pattern(const DBGSuccinct &graph, const std::string &sequence) {
    SCOPED_TRACE("k: " + std::to_string(graph.get_k()) + ", pattern: " + sequence
                    + (graph.get_mask() ? ", mask" : ", no mask"));
    boss::Pattern pattern = boss::parse_pattern(graph.get_boss(), sequence);
    Contexts expected = brute_force(graph, pattern);

    Contexts contexts;
    boss::SearchBudget list_budget(1'000'000, kNoTimeLimit);
    EXPECT_TRUE(graph.call_pattern_contexts(pattern, &list_budget,
        [&](DBGSuccinct::node_index node, const std::string &kmer, size_t offset) {
            EXPECT_EQ(graph.get_node_sequence(node), kmer);
            EXPECT_TRUE(contexts.emplace(node, offset).second) << "a context called twice";
        }
    ));
    EXPECT_EQ(expected, contexts);

    uint64_t min, max;
    boss::SearchBudget count_budget(1'000'000, kNoTimeLimit);
    EXPECT_TRUE(graph.count_pattern(pattern, &count_budget, &min, &max));
    if (graph.get_mask()) {
        EXPECT_EQ(expected.size(), min);
        EXPECT_EQ(expected.size(), max);
    } else {
        // without a mask the count is exact at offset 0
        EXPECT_EQ(min, static_cast<uint64_t>(std::count_if(expected.begin(), expected.end(),
            [](const auto &context) { return !context.second; }
        )));
        EXPECT_LE(expected.size(), max);
    }
}

void test_patterns(DeBruijnGraph::Mode mode) {
    std::mt19937 gen(1);
    for (size_t k : { 2, 3, 4, 6, 9 }) {
        for (size_t i = 0; i < 4; ++i) {
            auto graph = random_graph(gen, k, mode);
            std::vector<std::string> patterns;
            for (size_t length = 1; length <= k; ++length) {
                patterns.push_back(std::string(length, 'N'));
                for (size_t j = 0; j < 3; ++j) {
                    patterns.push_back(random_sequence(gen, length, "ACGT"));
                    patterns.push_back(random_sequence(gen, length, "ACGTRYSWKMBDHVN"));
                }
            }
            for (const std::string &pattern : patterns) {
                test_pattern(*graph, pattern);
            }
            graph->mask_dummy_kmers(1, false);
            for (const std::string &pattern : patterns) {
                test_pattern(*graph, pattern);
            }
        }
    }
}

TEST(DBGSuccinct, PatternBasic) {
    test_patterns(DeBruijnGraph::BASIC);
}

TEST(DBGSuccinct, PatternCanonical) {
    test_patterns(DeBruijnGraph::CANONICAL);
}

TEST(DBGSuccinct, PatternPrimary) {
    test_patterns(DeBruijnGraph::PRIMARY);
}

TEST(DBGSuccinct, PatternDummySourceKmers) {
    // k-mers of the graph: ACGT, CGTA. Dummy source k-mers: $$$A, $$AC, $ACG.
    DBGSuccinct graph(4);
    graph.add_sequence("ACGTA");
    for (bool mask : { false, true }) {
        if (mask)
            graph.mask_dummy_kmers(1, false);

        boss::Pattern pattern = boss::parse_pattern(graph.get_boss(), "AC");
        std::vector<std::pair<std::string, size_t>> contexts;
        boss::SearchBudget budget(1'000, kNoTimeLimit);
        EXPECT_TRUE(graph.call_pattern_contexts(pattern, &budget,
            [&](auto, const std::string &kmer, size_t offset) {
                contexts.emplace_back(kmer, offset);
            }
        ));
        EXPECT_EQ((std::vector<std::pair<std::string, size_t>>{ { "ACGT", 0 } }), contexts);

        uint64_t min, max;
        EXPECT_TRUE(graph.count_pattern(pattern, &budget, &min, &max));
        EXPECT_EQ(1u, min);
        // $$AC and $ACG are candidates without a mask
        EXPECT_EQ(mask ? 1u : 3u, max);
    }
}

TEST(DBGSuccinct, PatternManyMaskedEdges) {
    // more than 10,000 masked edges in the range of all edges, and fewer in the others
    std::mt19937 gen(2);
    auto graph = random_graph(gen, 12, DeBruijnGraph::BASIC, 2'000);
    graph->mask_dummy_kmers(1, false);
    ASSERT_LT(10'000u, graph->max_index() - graph->num_nodes());

    for (const std::string sequence : { "A", "R", "N", "AC", "ACN" }) {
        boss::Pattern pattern = boss::parse_pattern(graph->get_boss(), sequence);
        uint64_t expected = brute_force(*graph, pattern).size();
        uint64_t min, max;
        boss::SearchBudget budget(10'000'000, kNoTimeLimit);
        EXPECT_TRUE(graph->count_pattern(pattern, &budget, &min, &max));
        EXPECT_LE(min, expected) << sequence;
        EXPECT_LE(expected, max) << sequence;
        // one symbol that is not N: its masked edges are too many to visit
        if (sequence == "A" || sequence == "R") {
            EXPECT_LT(min, max) << sequence;
        } else {
            EXPECT_EQ(min, max) << sequence;
        }
    }
}

TEST(DBGSuccinct, PatternStepBound) {
    std::mt19937 gen(3);
    for (bool mask : { false, true }) {
        auto graph = random_graph(gen, 9, DeBruijnGraph::BASIC);
        if (mask)
            graph->mask_dummy_kmers(1, false);

        boss::Pattern pattern = boss::parse_pattern(graph->get_boss(), "NN");
        uint64_t min, max;
        boss::SearchBudget count_budget(5, kNoTimeLimit);
        EXPECT_FALSE(graph->count_pattern(pattern, &count_budget, &min, &max));
        EXPECT_STREQ("max_steps", count_budget.stop());

        size_t called = 0;
        boss::SearchBudget list_budget(5, kNoTimeLimit);
        EXPECT_FALSE(graph->call_pattern_contexts(pattern, &list_budget,
            [&](auto, const std::string &, size_t) { ++called; }
        ));
        EXPECT_GE(5u, called);
    }
}

#endif // ! _PROTEIN_GRAPH

} // namespace
