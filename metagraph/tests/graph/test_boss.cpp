#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "dbg_succinct.hpp"
#include "boss.hpp"
#include "boss_construct.hpp"
#include "utils.hpp"
#include "reverse_complement.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

constexpr double kEps = std::numeric_limits<double>::epsilon();


void test_graph(BOSS *graph, const std::string &last,
                             const std::vector<uint64_t> &W,
                             const std::string &F,
                             Config::StateType state) {
    Config::StateType old_state = graph->get_state();
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << graph->get_last();
    EXPECT_EQ(last, ostr.str()) << "state: " << state
                                << ", old state: " << old_state;

    for (size_t i = 1; i < graph->get_last().size(); ++i) {
        EXPECT_EQ((graph->get_last())[i], graph->get_last(i));
        EXPECT_EQ((graph->get_W())[i], graph->get_W(i));

        auto last_outgoing = graph->succ_last(i);
        graph->call_adjacent_incoming_edges(i, [&](auto incoming) {
            EXPECT_EQ(last_outgoing, graph->fwd(incoming));
        });
    }

    ostr.clear();
    ostr.str("");

    auto W_vector = graph->get_W().to_vector();
    EXPECT_EQ(W, std::vector<uint64_t>(W_vector.begin(), W_vector.end()))
        << "state: " << state
        << ", old state: " << old_state;

    ostr.clear();
    ostr.str("");

    for (size_t i = 0; i < graph->get_F().size(); ++i) {
        ostr << graph->get_F()[i] << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state
                             << ", old state: " << old_state;
}


void test_graph(BOSS *graph, const std::string &last,
                             const std::vector<uint64_t> &W,
                             const std::string &F) {
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::STAT);
    test_graph(graph, last, W, F, Config::STAT);
    test_graph(graph, last, W, F, Config::SMALL);
    test_graph(graph, last, W, F, Config::SMALL);
    test_graph(graph, last, W, F, Config::STAT);
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::SMALL);
    test_graph(graph, last, W, F, Config::DYN);
}


TEST(BOSS, GraphDefaultConstructor) {
    BOSS *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new BOSS;
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        graph_ = new BOSS();
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        BOSS graph;
    });
}

#if _DNA5_GRAPH
TEST(BOSS, EmptyGraph) {
    BOSS *graph = new BOSS(3);
    test_graph(graph, "01", { 0, 0 }, "0 1 1 1 1 1 ");
    delete graph;
}

TEST(BOSS, SwitchState) {
    BOSS *graph = new BOSS(3);
    test_graph(graph, "01", { 0, 0 }, "0 1 1 1 1 1 ");
    delete graph;
}

TEST(BOSS, AddSequenceFast) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    BOSSConstructor constructor(3);

    std::vector<std::string> names;

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_sequences({ read_stream->seq.s });
        names.emplace_back(read_stream->name.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    EXPECT_EQ(4llu, names.size());
    EXPECT_EQ("1", names[0]);
    EXPECT_EQ("2", names[1]);
    EXPECT_EQ("3", names[2]);
    EXPECT_EQ("4", names[3]);

    BOSS *graph = new BOSS(&constructor);

    //test graph construction
    test_graph(graph, "00011101101111111111111",
                      { 0, 0, 1, 3, 1, 1, 2, 4,
                        4, 3, 4, 0, 1, 0, 1, 4,
                        1, 7, 2, 0, 4, 3, 3 },
                      "0 3 11 13 17 22 ");
    delete graph;
}
#endif

TEST(BOSS, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    BOSSConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_sequences({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    BOSS *graph = new BOSS(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    EXPECT_EQ(1u, graph->outgoing(1, BOSS::kSentinelCode));

    for (size_t i = 1; i <= graph->num_edges(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != BOSS::kSentinelCode) {
            uint64_t node_idx = graph->rank_last(i - 1) + 1;

            EXPECT_EQ(outgoing_edges[i],
                graph->select_last(graph->outgoing(node_idx, graph->get_W(i) % graph->alph_size))
            ) << "Edge index: " << i << "\n"
              << "Outgoing: " << graph->outgoing(node_idx, graph->get_W(i) % graph->alph_size) << "\n"
              << *graph;

            EXPECT_EQ(node_idx,
                graph->incoming(graph->outgoing(node_idx, graph->get_W(i) % graph->alph_size),
                                graph->get_minus_k_value(i, graph->get_k() - 1).first)
            );
            for (int c = 0; c < graph->alph_size; ++c) {
                uint64_t prev_node = graph->incoming(node_idx, c);
                if (prev_node) {
                    EXPECT_EQ(node_idx,
                        graph->outgoing(prev_node, graph->get_node_last_value(i))
                    );
                }
            }

            //test FM index property
            EXPECT_TRUE(graph->get_last(graph->fwd(i)));
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_value(graph->fwd(i)));
        }
    }

    delete graph;
}

TEST(BOSS, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    BOSSConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_sequences({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    BOSS *graph = new BOSS(&constructor);

    graph->serialize(test_dump_basename);

    BOSS loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";
    EXPECT_FALSE(BOSS() == loaded_graph);

    delete graph;
}

TEST(BOSS, AddSequenceSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        EXPECT_EQ(k + 1, graph.num_nodes());
        EXPECT_EQ(k + 2, graph.num_edges());
    }
}

TEST(BOSS, CountDummyEdgesSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'G');
        EXPECT_EQ(2u, graph.num_edges()
                        - graph.mark_sink_dummy_edges()
                        - graph.mark_source_dummy_edges()) << graph;
    }
}

TEST(BOSS, CountDummyEdgesSimplePathParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'G');
        EXPECT_EQ(2u, graph.num_edges()
                        - graph.mark_sink_dummy_edges()
                        - graph.mark_source_dummy_edges(NULL, 10)) << graph;
    }
}

TEST(BOSS, CountDummyEdgesTwoPaths) {
    for (size_t k = 1; k < 40; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'G');
        graph.add_sequence(std::string(100, 'C') + 'T');
        EXPECT_EQ(4u, graph.num_edges()
                        - graph.mark_sink_dummy_edges()
                        - graph.mark_source_dummy_edges()) << graph;
    }
}

TEST(BOSS, CountDummyEdgesTwoPathsParallel) {
    for (size_t k = 1; k < 40; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'G');
        graph.add_sequence(std::string(100, 'C') + 'T');
        EXPECT_EQ(4u, graph.num_edges()
                        - graph.mark_sink_dummy_edges()
                        - graph.mark_source_dummy_edges(NULL, 10)) << graph;
    }
}

TEST(BOSS, MarkDummySinkEdgesSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'G');
        std::vector<bool> sink_nodes(graph.num_edges() + 1);
        sink_nodes.back() = true;
        std::vector<bool> sink_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(1u, graph.mark_sink_dummy_edges(&sink_nodes_result));
        EXPECT_EQ(sink_nodes, sink_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySinkEdgesTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'T');
        graph.add_sequence(std::string(100, 'A') + 'G');
        std::vector<bool> sink_nodes(graph.num_edges() + 1);
        sink_nodes[sink_nodes.size() - 2] = true;
        sink_nodes[sink_nodes.size() - 1] = true;
        std::vector<bool> sink_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(2u, graph.mark_sink_dummy_edges(&sink_nodes_result));
        EXPECT_EQ(sink_nodes, sink_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySourceEdgesSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        std::vector<bool> source_nodes(graph.num_edges() + 1, true);
        source_nodes.front() = false;
        source_nodes.back() = false;
        std::vector<bool> source_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(int(k + 1), int(std::count(source_nodes.begin(),
                                             source_nodes.end(), true)));
        ASSERT_EQ(k + 1, graph.mark_source_dummy_edges(&source_nodes_result));
        EXPECT_EQ(source_nodes, source_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySourceEdgesTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        graph.add_sequence(std::string(100, 'C'));
        std::vector<bool> source_nodes(graph.num_edges() + 1, true);
        source_nodes.front() = false;
        source_nodes[1 + 2 + k] = false;
        source_nodes[1 + 2 + 2 * k] = false;
        std::vector<bool> source_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(int(2 * k + 1), int(std::count(source_nodes.begin(),
                                                 source_nodes.end(), true)));
        ASSERT_EQ(2 * k + 1, graph.mark_source_dummy_edges(&source_nodes_result));
        EXPECT_EQ(source_nodes, source_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySourceEdgesSimplePathParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        std::vector<bool> source_nodes(graph.num_edges() + 1, true);
        source_nodes.front() = false;
        source_nodes.back() = false;
        std::vector<bool> source_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(int(k + 1), int(std::count(source_nodes.begin(),
                                             source_nodes.end(), true)));
        ASSERT_EQ(k + 1, graph.mark_source_dummy_edges(&source_nodes_result, 10));
        EXPECT_EQ(source_nodes, source_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySourceEdgesTwoPathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        graph.add_sequence(std::string(100, 'C'));
        std::vector<bool> source_nodes(graph.num_edges() + 1, true);
        source_nodes.front() = false;
        source_nodes[1 + 2 + k] = false;
        source_nodes[1 + 2 + 2 * k] = false;
        std::vector<bool> source_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(int(2 * k + 1), int(std::count(source_nodes.begin(),
                                                 source_nodes.end(), true)));
        ASSERT_EQ(2 * k + 1, graph.mark_source_dummy_edges(&source_nodes_result, 10));
        EXPECT_EQ(source_nodes, source_nodes_result) << graph;
    }
}

TEST(BOSS, RemoveDummyEdgesForClearGraph) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<BOSS> first_ptr;
        std::unique_ptr<BOSS> second_ptr;

        {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'C'),
                                        std::string(100, 'G'),
                                        "ACATCTAGTAGTCGATCGTACG",
                                        "ATTAGTAGTAGTAGTGATGTAG", });
            first_ptr.reset(new BOSS(&constructor));
        }

        {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'C'),
                                        std::string(100, 'G'),
                                        "ACATCTAGTAGTCGATCGTACG",
                                        "ATTAGTAGTAGTAGTGATGTAG", });
            second_ptr.reset(new BOSS(&constructor));
        }

        auto &first = *first_ptr;
        auto &second = *second_ptr;

        ASSERT_TRUE(first.equals_internally(second)) << first;

        std::vector<bool> source_dummy_edges(second.num_edges() + 1, false);
        auto to_remove = second.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(second.num_edges() + 1, false);
        second.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_EQ(0u, std::count(to_remove.begin(), to_remove.end(), true));
        EXPECT_TRUE(first.equals_internally(second)) << first;
    }
}

TEST(BOSS, RemoveDummyEdgesLinear) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(20, 'A'));
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesThreePaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'), });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesFourPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    std::string(20, 'T') + 'A' });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesFivePaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    std::string(20, 'T'),
                                    std::string(20, 'N') + 'A', });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T'));
        dynamic_graph.add_sequence(std::string(20, 'N') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdges) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence("ACATCTAGTAGTCGATCGTACG");
        dynamic_graph.add_sequence("ATTAGTAGTAGTAGTGATGTAG");

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesForClearGraphParallel) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<BOSS> first_ptr;
        std::unique_ptr<BOSS> second_ptr;

        {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'C'),
                                        std::string(100, 'G'),
                                        "ACATCTAGTAGTCGATCGTACG",
                                        "ATTAGTAGTAGTAGTGATGTAG", });
            first_ptr.reset(new BOSS(&constructor));
        }

        {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'C'),
                                        std::string(100, 'G'),
                                        "ACATCTAGTAGTCGATCGTACG",
                                        "ATTAGTAGTAGTAGTGATGTAG", });
            second_ptr.reset(new BOSS(&constructor));
        }

        auto &first = *first_ptr;
        auto &second = *second_ptr;

        ASSERT_TRUE(first.equals_internally(second)) << first;

        std::vector<bool> source_dummy_edges(second.num_edges() + 1, false);
        auto to_remove = second.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(second.num_edges() + 1, false);
        second.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_EQ(0u, std::count(to_remove.begin(), to_remove.end(), true));
        EXPECT_TRUE(first.equals_internally(second)) << first;
    }
}

TEST(BOSS, RemoveDummyEdgesLinearParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(20, 'A'));
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesThreePathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'), });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesFourPathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    std::string(20, 'T') + 'A' });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesFivePathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    std::string(20, 'T'),
                                    std::string(20, 'N') + 'A', });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T'));
        dynamic_graph.add_sequence(std::string(20, 'N') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, RemoveDummyEdgesParallel) {
    for (size_t k = 1; k < 10; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(20, 'A'),
                                    std::string(20, 'C'),
                                    std::string(20, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
        BOSS graph(&constructor);

        BOSS dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence("ACATCTAGTAGTCGATCGTACG");
        dynamic_graph.add_sequence("ATTAGTAGTAGTAGTGATGTAG");

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        std::vector<bool> source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        std::vector<bool> source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << to_sdsl(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(BOSS, AddSequenceBugRevealingTestcase) {
    BOSS graph(1);
    graph.add_sequence("CTGAG", false);
}

// TODO: these are duplicates of tests in test_dbg
TEST(BOSS, NonASCIIStrings) {
    BOSSConstructor constructor_first(5);
    constructor_first.add_sequences({
        // cyrillic A and C
        "АСАСАСАСАСАСА",
        "плохая строка",
        "АСАСАСАСАСАСА"
    });
    BOSS graph(&constructor_first);
#if _DNA_GRAPH
    ASSERT_EQ(1u, graph.num_edges()) << graph;
#else
    ASSERT_EQ(2u, graph.num_edges()) << graph;
#endif
}

TEST(BOSS, AddSequence) {
    {
        BOSS graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
    {
        BOSS graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        graph.add_sequence("GAAC");
        EXPECT_EQ(11u, graph.num_nodes());
        EXPECT_EQ(14u, graph.num_edges());
    }
    {
        BOSS graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("AACG");
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
#ifndef _DNA_GRAPH
    {
        BOSS graph(4);
        graph.add_sequence("AGACN");
        graph.add_sequence("GACTN");
        graph.add_sequence("ACTAN");
        EXPECT_EQ(15u, graph.num_nodes());
        EXPECT_EQ(18u, graph.num_edges());
    }
#endif
}

TEST(BOSS, AppendSequence) {
    {
        BOSS graph(3);
        graph.add_sequence("AAAC", true);
        graph.add_sequence("AACG", true);
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(7u, graph.num_edges());
    }
    {
        BOSS graph(3);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(7u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
#ifndef _DNA_GRAPH
    {
        BOSS graph(4);
        graph.add_sequence("AGACN", true);
        graph.add_sequence("GACTN", true);
        graph.add_sequence("ACTAN", true);
        EXPECT_EQ(15u, graph.num_nodes());
        EXPECT_EQ(18u, graph.num_edges());
    }
#endif
    {
        BOSS graph(3);
        graph.add_sequence("AAACT", true);
        graph.add_sequence("AAATG", true);
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
}

TEST(BOSS, AppendSequenceAnyKmerSize) {
    for (size_t k = 1; k < 10; ++k) {
        {
            BOSS graph(k);
            graph.add_sequence("AAAC", true);
            graph.add_sequence("AACG", true);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACT", true);
            graph.add_sequence("AAATG", true);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
        }
    }
}

TEST(BOSS, CallPathsEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        BOSS empty(k);
        BOSS reconstructed(k);

        empty.call_sequences([&](const auto &sequence) {
            reconstructed.add_sequence(sequence);
        });

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(BOSS, CallUnitigsEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        BOSS empty(k);
        BOSS reconstructed(k);

        empty.call_unitigs([&](const auto &sequence) {
            reconstructed.add_sequence(sequence);
        });

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(BOSS, CallPathsOneLoop) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);

        ASSERT_EQ(1u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallUnitigsOneLoop) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);

        ASSERT_EQ(1u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_unitigs([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallPathsTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallUnitigsTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_unitigs([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallPathsFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        BOSS graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, false);
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallUnitigsFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        BOSS graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; }, true);
        graph.call_unitigs([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(BOSS, CallPaths) {
    for (size_t k = 1; k < 10; ++k) {
        {
            BOSS graph(k);
            graph.add_sequence("AAACACTAG", true);
            graph.add_sequence("AACGACATG", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGACACTGA", true);
            graph.add_sequence("GACTACGTA", true);
            graph.add_sequence("ACTAACGTA", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGACACAGT", true);
            graph.add_sequence("GACTTGCAG", true);
            graph.add_sequence("ACTAGTCAG", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACTCGTAGC", true);
            graph.add_sequence("AAATGCGTAGC", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
    }
}

TEST(BOSS, CallUnitigs) {
    for (size_t k = 1; k < 10; ++k) {
        {
            BOSS graph(k);
            graph.add_sequence("AAACACTAG", true);
            graph.add_sequence("AACGACATG", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_unitigs([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGACACTGA", true);
            graph.add_sequence("GACTACGTA", true);
            graph.add_sequence("ACTAACGTA", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_unitigs([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AGACACAGT", true);
            graph.add_sequence("GACTTGCAG", true);
            graph.add_sequence("ACTAGTCAG", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_unitigs([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACTCGTAGC", true);
            graph.add_sequence("AAATGCGTAGC", true);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_unitigs([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            BOSS graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
            graph.switch_state(Config::STAT);

            BOSS reconstructed(k);

            graph.call_unitigs([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
    }
}

TEST(BOSS, CallUnitigs1) {
    BOSSConstructor constructor(3);
    constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                "ACTCT" });
    BOSS graph(&constructor);

    size_t num_contigs = 0;

    std::set<std::string> contigs {
        "$$$$",
        "$$$ACT",
        "ACTCT$",
        "ACTA",
        "CTAGCTA",
    };

    graph.call_paths(
        [&](const auto &, const auto &seq) {
            auto str = graph.decode(seq);
            EXPECT_TRUE(contigs.count(str)) << str;
            num_contigs++;
        },
        true
    );

    EXPECT_EQ(contigs.size(), num_contigs);
}

TEST(BOSS, CallUnitigsDisconnected1) {
    BOSSConstructor constructor(3);
    constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                "ACTCT",
                                "ATCATCATCATCATCATCAT" });
    BOSS graph(&constructor);

    size_t num_contigs = 0;

    std::set<std::string> contigs {
        "$$$$",
        "$$$ACT",
        "ACTCT$",
        "ACTA",
        "CTAGCTA",
        "TCATCA",
    };

    graph.call_paths(
        [&](const auto &, const auto &seq) {
            auto str = graph.decode(seq);
            EXPECT_TRUE(contigs.count(str)) << str;
            num_contigs++;
        },
        true
    );

    EXPECT_EQ(contigs.size(), num_contigs);
}

#ifndef _DNA_GRAPH
TEST(BOSS, CallUnitigsDisconnected2) {
    BOSSConstructor constructor(3);
    constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                "ACTCT",
                                "ATCATCATCATCATCATCAT",
                                "ATNATNATNATNATNATNAT" });
    BOSS graph(&constructor);

    size_t num_contigs = 0;

    std::set<std::string> contigs {
        "$$$$",
        "$$$ACT",
        "ACTCT$",
        "ACTA",
        "CTAGCTA",
        "TCATCA",
        "TNATNA",
    };

    graph.call_paths(
        [&](const auto &, const auto &seq) {
            auto str = graph.decode(seq);
            EXPECT_TRUE(contigs.count(str)) << str;
            num_contigs++;
        },
        true
    );

    EXPECT_EQ(contigs.size(), num_contigs);
}

TEST(BOSS, CallUnitigsTwoComponents) {
    BOSSConstructor constructor(3);
    constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                "ACTCT",
                                "ATCATCATCATCATCATCAT",
                                "ATNATNATNATNATNATNAT",
                                "ATCNATCNATCNATCNATCNATCNAT" });
    BOSS graph(&constructor);

    size_t num_contigs = 0;

    std::set<std::string> contigs {
        "$$$$",
        "$$$ACT",
        "ACTCT$",
        "ACTA",
        "CTAGCTA",
        "NATC",
        "ATCNAT",
        "NATNAT",
        "ATCATC",
    };

    graph.call_paths(
        [&](const auto &, const auto &seq) {
            auto str = graph.decode(seq);
            EXPECT_TRUE(contigs.count(str)) << str;
            num_contigs++;
        },
        true
    );

    EXPECT_EQ(contigs.size(), num_contigs);
}
#endif

TEST(BOSS, CallUnitigsWithPruning) {
    BOSSConstructor constructor(4);
    constructor.add_sequences({
        "ACTATAGCTAGTCTATGCGA",
        "ACTATAGCTAGTCTAG",
        "ACTATAGCTAN",
        "ACTATAGCTT",
        "ACTATT",
    });
    BOSS graph(&constructor);

#ifndef _DNA_GRAPH
    {
        std::set<std::string> contigs {
            "ACTAT",
            "CTATT",
            "CTATGCGA",
            "CTATAGCT",
            "AGCTT",
            "AGCTA",
            "GCTAN",
            "GCTAG",
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        );
        EXPECT_EQ(contigs.size(), num_contigs);
    }
#endif
    {
        std::set<std::string> contigs {
            "CTATGCGA",
            "CTATAGCT",
#if _DNA_GRAPH
            "AGCTAG",
#else
            "AGCTA",
            "GCTAG",
#endif
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        , 1);
        EXPECT_EQ(contigs.size(), num_contigs);
    }
    {
        std::set<std::string> contigs {
            "CTATGCGA",
            "CTATAGCT",
#if _DNA_GRAPH
            "AGCTAG",
#else
            "AGCTA",
            "GCTAG",
#endif
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        , 2);
        EXPECT_EQ(contigs.size(), num_contigs);
    }
    {
        std::set<std::string> contigs {
            "CTATGCGA",
            "CTATAGCT",
#if _DNA_GRAPH
            "AGCTAG",
#else
            "AGCTA",
            "GCTAG",
#endif
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        , 3);
        EXPECT_EQ(contigs.size(), num_contigs);
    }
    {
        std::set<std::string> contigs {
            "CTATAGCT",
#if _DNA_GRAPH
            "AGCTAG",
#else
            "AGCTA",
            "GCTAG",
#endif
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        , 4);
        EXPECT_EQ(contigs.size(), num_contigs);
    }
    {
        std::set<std::string> contigs {
            "CTATAGCT",
#if _DNA_GRAPH
            "AGCTAG",
#else
            "AGCTA",
            "GCTAG",
#endif
            "CTAGTCTA",
            "TCTAT",
            "TCTAG",
        };
        size_t num_contigs = 0;
        graph.call_unitigs(
            [&](const auto &str) {
                EXPECT_TRUE(contigs.count(str)) << str;
                num_contigs++;
            }
        , 5);
        EXPECT_EQ(contigs.size(), num_contigs);
    }
}

TEST(BOSS, CallEdgesEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        BOSS empty(k);

        size_t num_edges = 0;
        empty.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(empty.get_k() + 1, edge.size()) << empty;
            num_edges++;
        });

        EXPECT_EQ(empty.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        BOSS graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesFourLoopsDynamic) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        graph.add_sequence(std::string(100, 'G'));
        graph.add_sequence(std::string(100, 'C'));

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesTestPath) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + std::string(k, 'C'));

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesTestPathACA) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A')
                            + std::string(k, 'C')
                            + std::string(100, 'A'));

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(BOSS, CallEdgesTestPathDisconnected) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        BOSS graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(100, 'T'));

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges) << graph;
    }
}

TEST(BOSS, CallEdgesTestPathDisconnected2) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'G'));
        BOSS graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_edges = 0;
        graph.call_edges([&](auto edge_idx, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            EXPECT_EQ(graph.get_node_seq(edge_idx),
                      std::vector<BOSS::TAlphabet>(edge.begin(), edge.end() - 1))
                << edge_idx << "\n" << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges) << graph;
    }
}

TEST(BOSS, CallKmersEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        BOSS empty(k);

        size_t num_kmers = 0;
        empty.call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TEST(BOSS, CallKmersTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A') });
        BOSS graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TEST(BOSS, CallKmersFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequences({ std::string(100, 'A'),
                                    std::string(100, 'G'),
                                    std::string(100, 'C') });
        BOSS graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TEST(BOSS, CallKmersFourLoopsDynamic) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        graph.add_sequence(std::string(100, 'G'));
        graph.add_sequence(std::string(100, 'C'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_TRUE(std::string(k, 'A') == sequence
                        || std::string(k, 'G') == sequence
                        || std::string(k, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers);
    }
}

TEST(BOSS, CallKmersTestPath) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + std::string(k, 'C'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k + 1, num_kmers);
    }
}

TEST(BOSS, CallKmersTestPathACA) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A')
                            + std::string(k, 'C')
                            + std::string(100, 'A'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k, num_kmers);
    }
}

TEST(BOSS, CallKmersTestPathDisconnected) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        BOSS graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(100, 'T'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers) << graph;
    }
}

TEST(BOSS, CallKmersTestPathDisconnected2) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'G'));
        BOSS graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(3u, num_kmers) << graph;
    }
}

void test_pred_kmer(const BOSS &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::vector<BOSS::TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(),
                   [&graph](char c) {
                       return c == BOSS::kSentinel
                                   ? BOSS::kSentinelCode
                                   : graph.encode(c);
                   });
    EXPECT_EQ(expected_idx, graph.select_last(graph.pred_kmer(kmer)))
        << kmer_s << std::endl
        << graph;
}

TEST(BOSS, PredKmer) {
    {
        BOSS graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        BOSS graph(5);
        graph.add_sequence("AAAAAA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        BOSS graph(5);
        graph.add_sequence("ACACAA");

        test_pred_kmer(graph, "ACGCG", 8);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 8);
        test_pred_kmer(graph, "NNNNN", 8);
        test_pred_kmer(graph, "$$$$$", 2);
    }
#ifndef _PROTEIN_GRAPH
    {
        BOSS graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        BOSS graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        BOSS graph(2);
        graph.add_sequence("ATAATATCC");
        graph.add_sequence("ATACGC");
        graph.add_sequence("ATACTC");
        graph.add_sequence("ATACTA");
        graph.add_sequence("CATT");
        graph.add_sequence("CCC");
        graph.add_sequence("GGGC");
        graph.add_sequence("GGTGTGAC");
        graph.add_sequence("GGTCT");
        graph.add_sequence("GGTA");

        test_pred_kmer(graph, "$$", 4);
        test_pred_kmer(graph, "$A", 5);
        test_pred_kmer(graph, "$T", 26);
        test_pred_kmer(graph, "AT", 29);
        test_pred_kmer(graph, "TT", 35);
        test_pred_kmer(graph, "NT", 35);
        test_pred_kmer(graph, "TN", 35);
    }
#endif
}

TEST(BOSS, PredKmerRandomTest) {
    srand(1);

    for (size_t k = 1; k < 8; ++k) {
        BOSS graph(k);

        for (size_t p = 0; p < 10; ++p) {
            size_t length = rand() % 400;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = graph.alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = graph.alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, true);
        }

        auto all_kmer_str = utils::generate_strings("ACGT", k);
        for (size_t i = 1; i < k; ++i) {
            auto kmer_str_suffices = utils::generate_strings("ACGT", i);
            for (size_t j = 0; j < kmer_str_suffices.size(); ++j) {
                all_kmer_str.push_back(std::string(k - i, '$')
                                        + kmer_str_suffices[j]);
            }
        }

        for (const auto &kmer_str : all_kmer_str) {
            std::vector<BOSS::TAlphabet> kmer = graph.encode(kmer_str);

            uint64_t lower_bound = graph.select_last(graph.pred_kmer(kmer));

            EXPECT_FALSE(
                utils::colexicographically_greater(
                    graph.get_node_seq(lower_bound), kmer
                )
            ) << graph
              << "kmer: " << kmer_str << std::endl
              << "lower bound: " << lower_bound << std::endl
              << "which is: " << graph.get_node_str(lower_bound) << std::endl;

            if (lower_bound < graph.get_W().size() - 1) {
                EXPECT_TRUE(
                    utils::colexicographically_greater(
                        graph.get_node_seq(lower_bound + 1), kmer
                    )
                );
            }
        }
    }
}

TEST(BOSS, FindSequence) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS *graph = new BOSS(k);

        graph->add_sequence(std::string(100, 'A'));

        uint64_t index = 777;
        graph->map_to_nodes(std::string(k - 1, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(777u, index) << *graph;
        graph->map_to_nodes(std::string(k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 1, index) << *graph;

        index = 777;
        graph->map_to_nodes(std::string(2 * k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 1, index) << *graph;

        index = 777;
        graph->map_to_edges(std::string(k - 1, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(777u, index);
        graph->map_to_edges(std::string(k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(777u, index);
        graph->map_to_edges(std::string(k + 1, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 2, index) << *graph;

        index = 777;
        graph->map_to_edges(std::string(2 * k, 'A'), [&](uint64_t i) { index = i; });
        EXPECT_EQ(k + 2, index) << *graph;

        EXPECT_FALSE(graph->find(std::string(k - 1, 'A')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'A')));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0));
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

        EXPECT_FALSE(graph->find(std::string(k - 1, 'B')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'B')));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B')));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B')));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'B'), 0));

        std::string pattern = std::string(k, 'A') + std::string(k, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 3)));
        EXPECT_FALSE(graph->find(pattern, 0.0 / k + kEps));
        EXPECT_TRUE(graph->find(pattern, 0.0 / k - kEps));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k + 1, 'A') + std::string(k + 1, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 2)));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 2) + kEps));
        // we map (k+1)-mers to the graph edges
        // just 1 out of k+2 (k+1)-mers can be mapped here
        EXPECT_TRUE(graph->find(pattern, 1.0 / (k + 2)));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k + 2, 'A') + std::string(k + 2, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 3.0 / (k + 4)));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 4) + kEps));
        EXPECT_TRUE(graph->find(pattern, 2.0 / (k + 4)));
        EXPECT_TRUE(graph->find(pattern, 0));

        delete graph;
    }
}

TEST(BOSS, FindSequenceDBG) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(k)) };

        graph->add_sequence(std::string(100, 'A'));

        EXPECT_FALSE(graph->find(std::string(k - 1, 'A')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'A'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'A')));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 1));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k, 'A'), 0));
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

        EXPECT_FALSE(graph->find(std::string(k - 1, 'B')));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k - 1, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k, 'B')));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0.25));
        EXPECT_FALSE(graph->find(std::string(k, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B')));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 1, 'B'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 1, 'B'), 0));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B')));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 1));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.75));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.5));
        EXPECT_FALSE(graph->find(std::string(k + 2, 'B'), 0.25));
        EXPECT_TRUE(graph->find(std::string(k + 2, 'B'), 0));

        std::string pattern = std::string(k, 'A') + std::string(k, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 3)));
        EXPECT_FALSE(graph->find(pattern, 0.0 / k + kEps));
        EXPECT_TRUE(graph->find(pattern, 0.0 / k - kEps));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k + 1, 'A') + std::string(k + 1, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 2)));
        EXPECT_FALSE(graph->find(pattern, 1.0 / (k + 2) + kEps));
        // we map (k+1)-mers to the graph edges
        // just 1 out of k+2 (k+1)-mers can be mapped here
        EXPECT_TRUE(graph->find(pattern, 1.0 / (k + 2)));
        EXPECT_TRUE(graph->find(pattern, 0));

        pattern = std::string(k + 2, 'A') + std::string(k + 2, 'B');
        EXPECT_FALSE(graph->find(pattern));
        EXPECT_FALSE(graph->find(pattern, 1));
        EXPECT_FALSE(graph->find(pattern, 3.0 / (k + 4)));
        EXPECT_FALSE(graph->find(pattern, 2.0 / (k + 4) + kEps));
        EXPECT_TRUE(graph->find(pattern, 2.0 / (k + 4)));
        EXPECT_TRUE(graph->find(pattern, 0));
    }
}

TEST(BOSS, KmerMappingMode) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);

        const std::string check(k + 1, 'A');
        graph.add_sequence(std::string(100, 'A'));

        // kmer mapping modes
        for (size_t h = 0; h < 3; ++h) {
            EXPECT_TRUE(graph.find(check, 1, h));

            std::string query(100, 'A');

            for (double a : { 1.0, 0.75, 0.5, 0.25, 0.0 }) {
                EXPECT_TRUE(graph.find(query, a, h)) << "Mode: " << h;
            }

            // number of mutations
            for (size_t n = 1; n <= 100; ++n) {
                std::string query(100, 'A');
                for (size_t i = 0; i < query.size(); i += query.size() / n) {
                    query[i] = 'B';
                }

                uint64_t num_good_kmers = 0;
                for (size_t i = 0; i + k + 1 <= query.size(); ++i) {
                    if (query.substr(i, k + 1) == check)
                        num_good_kmers++;
                }
                double good_fraction = static_cast<double>(num_good_kmers)
                                            / (query.size() - k);

                for (double a : { 1.0, 0.75, 0.5, 0.25, 0.0 }) {
                    EXPECT_EQ(a <= good_fraction, graph.find(query, a, h))
                        << "Mode: " << h;
                }
            }
        }
    }
}

TEST(BOSS, Traversals) {
    for (size_t k = 1; k < 10; ++k) {
        auto graph = std::make_unique<BOSS>(k);

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_EQ(k + 1, it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it + 1, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it + 1, 'A'));
        EXPECT_EQ(BOSS::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(BOSS::npos, graph->traverse_back(it + 1, 'G'));
    }
}

TEST(BOSS, map_to_nodes) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<BOSS> graph { new BOSS(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::vector<BOSS::node_index> expected_result {
            SequenceGraph::npos, SequenceGraph::npos,
            k + 1, k + 1, k + 1, k + 1
        };
        for (size_t i = 1; i <= k; ++i) {
            expected_result.push_back(k + 1 + i);
        }
        for (size_t i = 0; i < k; ++i) {
            expected_result.push_back(k + 1 + k);
        }

        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 3, 'A')
                                        + std::string(2 * k, 'C');

        EXPECT_EQ(expected_result, graph->map_to_nodes(sequence_to_map));

        size_t pos = 0;
        graph->map_to_nodes(sequence_to_map,
                            [&](auto i) { EXPECT_EQ(expected_result[pos++], i); });
    }
}

TEST(BOSS, map_to_edges) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<BOSS> graph { new BOSS(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::vector<BOSS::node_index> expected_result {
            SequenceGraph::npos, SequenceGraph::npos,
            k + 2, k + 2, k + 2
        };
        for (size_t i = 1; i <= k; ++i) {
            expected_result.push_back(k + 2 + i);
        }
        for (size_t i = 0; i < k; ++i) {
            expected_result.push_back(k + 2 + k + 1);
        }

        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 3, 'A')
                                        + std::string(2 * k, 'C');

        EXPECT_EQ(expected_result, graph->map_to_edges(sequence_to_map));

        size_t pos = 0;
        graph->map_to_edges(sequence_to_map,
                            [&](auto i) { EXPECT_EQ(expected_result[pos++], i); });
    }
}

TEST(BOSS, map_to_nodes_BOSS_vs_DBG) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(k)) };
        std::unique_ptr<BOSS> boss { new BOSS(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));
        boss->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::string sequence_to_map = std::string(2, 'T')
                                        + std::string(k + 3, 'A')
                                        + std::string(2 * k, 'C');

        std::vector<BOSS::edge_index> boss_indexes;

        boss->map_to_edges(
            sequence_to_map,
            [&boss_indexes](auto i) { boss_indexes.push_back(i); }
        );
        std::vector<BOSS::edge_index> dbg_indexes;

        graph->map_to_nodes(
            sequence_to_map,
            [&dbg_indexes](auto i) { dbg_indexes.push_back(i); }
        );

        EXPECT_EQ(boss_indexes, dbg_indexes);
    }
}

TEST(BOSS, get_degree_with_source_dummy) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGSuccinct>(k);

        graph->add_sequence(std::string(k, 'A')
                                + std::string(k - 1, 'C')
                                + std::string(k - 1, 'G')
                                + std::string(k, 'T'));

        // dummy source k-mer: '$$$$$'
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(1));
        EXPECT_EQ(2ull, graph->outdegree(1));
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

TEST(BOSS, get_degree_with_source_and_sink_dummy) {
    for (size_t k = 2; k < 10; ++k) {
        auto graph = std::make_unique<DBGSuccinct>(k);

        graph->add_sequence(std::string(k, 'A')
                                + std::string(k - 1, 'C')
                                + std::string(k - 1, 'G')
                                + std::string(k - 1, 'T'));

        // dummy source k-mer: '$$$$$'
        EXPECT_EQ(std::string(k, '$'), graph->get_node_sequence(1));
        EXPECT_EQ(2ull, graph->outdegree(1));
        EXPECT_EQ(1ull, graph->indegree(1));

        // 'AAAAA'
        auto node_A = graph->kmer_to_node(std::string(k, 'A'));
        ASSERT_NE(DBGSuccinct::npos, node_A);
        EXPECT_EQ(2ull, graph->outdegree(node_A));
        EXPECT_EQ(2ull, graph->indegree(node_A));

        auto node_T = graph->kmer_to_node('G' + std::string(k - 1, 'T'));
        ASSERT_NE(DBGSuccinct::npos, node_T);
        EXPECT_EQ(1ull, graph->outdegree(node_T));
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

TEST(BOSS, get_node_sequence) {
    size_t k = 4;
    std::string reference = "AGCTTCGAGGCCAA";
    std::string query = "AGCT";

    auto graph = std::make_unique<DBGSuccinct>(k);
    graph->add_sequence(reference);

    std::string mapped_query = "";
    graph->map_to_nodes(query, [&](DBGSuccinct::node_index node) {
        mapped_query += graph->get_node_sequence(node);
    });

    EXPECT_EQ(query, mapped_query);
}

TEST(BOSS, is_single_outgoing_simple) {
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

TEST(BOSS, is_single_outgoing_for_multiple_valid_edges) {
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

TEST(BOSS, IndegreeIncomingIdentity) {
    int k_kmer = 5;

    auto graph_constructor = BOSSConstructor(k_kmer - 1);
    for (const auto &read : { "ATGCGATCGATATGCGAGA",
                              "ATGCGATCGAGACTACGAG",
                              "GTACGATAGACATGACGAG",
                              "ACTGACGAGACACAGATGC" }) {
        graph_constructor.add_sequence(read);
    }
    DBGSuccinct graph(new BOSS(&graph_constructor));
    graph.mask_dummy_kmers(1, false);

    for (uint64_t node = 1; node <= graph.num_nodes(); node++) {
        size_t computed_indegree = 0;
        for (auto c : {'A', 'C', 'G', 'T', 'N'}) {
            if (graph.traverse_back(node, c))
                computed_indegree++;
        }
        ASSERT_EQ(graph.indegree(node), computed_indegree);

        std::vector<DBGSuccinct::node_index> incoming_nodes;
        graph.adjacent_incoming_nodes(node, &incoming_nodes);
        ASSERT_EQ(graph.indegree(node), incoming_nodes.size());
    }
}

TEST(BOSS, EraseEdgesDynSingle) {
    for (size_t k = 1; k < 40; ++k) {
        BOSS graph2(k);
        graph2.add_sequence("AGACACAGT", true);
        graph2.add_sequence("GACTTGCAG", true);
        graph2.add_sequence("ACTAGTCAG", true);
        if (true) {
            for(uint64_t m = 1; m <= graph2.num_edges(); ++m) {
                BOSS graph(k);
                graph.add_sequence("AGACACAGT", true);
                graph.add_sequence("GACTTGCAG", true);
                graph.add_sequence("ACTAGTCAG", true);
                //graph.print_internal_representation();
                //graph.print();
                EXPECT_TRUE(graph.is_valid());
                graph.erase_edges_dyn({m});
                //graph.print_internal_representation();
                //graph.print();
                //std::cout << "m: " << m << std::endl;
                EXPECT_TRUE(graph.is_valid());
                //dump unitigs and visualize
                //std::cout << "viz" << std::endl;
                //std::cout << k+1 << std::endl;
                //graph.call_paths([&](const auto &sequence) {
                //    std::cout << sequence << std::endl;
                //});
                graph.call_paths(
                    [&](const auto &, const auto &seq) {
                        auto str = graph.decode(seq);
                        //std::cout << str << std::endl;
                    },
                    false
                );
            }
        }
    }
}
