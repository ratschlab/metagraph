#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <zlib.h>
#include <htslib/kseq.h>
#include <unordered_set>

#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "common/algorithms.hpp"


namespace {

KSEQ_INIT(gzFile, gzread);

using namespace mtg;
using namespace mtg::graph::boss;

const std::string test_data_dir = TEST_DATA_DIR;
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";


uint64_t fwd(const BOSS &boss, uint64_t i) {
    return boss.fwd(i, boss.get_W(i) % boss.alph_size);
}

void test_graph(BOSS *graph, const std::string &last,
                             const std::vector<uint64_t> &W,
                             const std::string &F,
                             BOSS::State state) {
    BOSS::State old_state = graph->get_state();
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << graph->get_last();
    EXPECT_EQ(last, ostr.str()) << "state: " << state
                                << ", old state: " << old_state;

    for (size_t i = 1; i < graph->get_last().size(); ++i) {
        EXPECT_EQ((graph->get_last())[i], graph->get_last(i));
        EXPECT_EQ((graph->get_W())[i], graph->get_W(i));

        auto last_outgoing = graph->succ_last(i);
        graph->call_incoming_to_target(graph->bwd(i), graph->get_node_last_value(i),
            [&](auto incoming) {
                EXPECT_EQ(last_outgoing, fwd(*graph, incoming));
            }
        );
    }

    ostr.clear();
    ostr.str("");

    auto W_vector = graph->get_W().to_vector();
    EXPECT_EQ(W, std::vector<uint64_t>(W_vector.begin(), W_vector.end()))
        << "state: " << state
        << ", old state: " << old_state;

    ostr.clear();
    ostr.str("");

    for (size_t c = 0; c < graph->alph_size; ++c) {
        ostr << graph->get_F(c) << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state
                             << ", old state: " << old_state;
}


void test_graph(BOSS *graph, const std::string &last,
                             const std::vector<uint64_t> &W,
                             const std::string &F) {
    test_graph(graph, last, W, F, BOSS::State::DYN);
    test_graph(graph, last, W, F, BOSS::State::DYN);
    test_graph(graph, last, W, F, BOSS::State::STAT);
    test_graph(graph, last, W, F, BOSS::State::STAT);
    test_graph(graph, last, W, F, BOSS::State::SMALL);
    test_graph(graph, last, W, F, BOSS::State::SMALL);
    test_graph(graph, last, W, F, BOSS::State::STAT);
    test_graph(graph, last, W, F, BOSS::State::DYN);
    test_graph(graph, last, W, F, BOSS::State::SMALL);
    test_graph(graph, last, W, F, BOSS::State::DYN);
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

TEST(BOSS, EmptyGraph) {
    BOSS *graph = new BOSS(3);
    std::string F_str = "0 ";
    for (size_t i = 1; i < graph->alphabet.size(); ++i) {
        F_str += "1 ";
    }
    test_graph(graph, "01", { 0, 0 }, F_str, graph->get_state());
    delete graph;
}

TEST(BOSS, SwitchState) {
    BOSS *graph = new BOSS(3);
    std::string F_str = "0 ";
    for (size_t i = 1; i < graph->alphabet.size(); ++i) {
        F_str += "1 ";
    }
    test_graph(graph, "01", { 0, 0 }, F_str);
    delete graph;
}

#if _DNA5_GRAPH
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

    BOSSConstructor constructor(3 /* k-mer length */);

    while (kseq_read(read_stream) >= 0) {
        constructor.add_sequences({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    BOSS *graph = new BOSS(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    uint64_t dummy_edge = graph->select_last(1);
    EXPECT_EQ(1u, graph->pick_edge(dummy_edge, BOSS::kSentinelCode));
    EXPECT_EQ(dummy_edge, fwd(*graph, 1));

    for (size_t i = 1; i <= graph->num_edges(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != BOSS::kSentinelCode) {
            EXPECT_EQ(outgoing_edges[i],
                fwd(*graph, graph->pick_edge(graph->succ_last(i),
                                             graph->get_W(i) % graph->alph_size))
            ) << "Edge index: " << i << "\n"
              << *graph;

            EXPECT_EQ(i,
                graph->pick_incoming_edge(
                    graph->bwd(fwd(*graph, i)),
                    graph->get_minus_k_value(i, graph->get_k() - 1).first
                )
            );
            for (int c = 0; c < graph->alph_size; ++c) {
                uint64_t prev_edge = graph->pick_incoming_edge(graph->bwd(i), c);
                if (prev_edge) {
                    EXPECT_EQ(i, graph->pick_edge(fwd(*graph, prev_edge),
                                                  graph->get_W(i) % graph->alph_size));
                }
            }

            //test FM index property
            EXPECT_TRUE(graph->get_last(fwd(*graph, i)));
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_value(fwd(*graph, i)));
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
        sdsl::bit_vector sink_nodes(graph.num_edges() + 1);
        sink_nodes[sink_nodes.size() - 1] = true;
        sdsl::bit_vector sink_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(1u, graph.mark_sink_dummy_edges(&sink_nodes_result));
        EXPECT_EQ(sink_nodes, sink_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySinkEdgesTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + 'T');
        graph.add_sequence(std::string(100, 'A') + 'G');
        sdsl::bit_vector sink_nodes(graph.num_edges() + 1);
        sink_nodes[sink_nodes.size() - 2] = true;
        sink_nodes[sink_nodes.size() - 1] = true;
        sdsl::bit_vector sink_nodes_result(graph.num_edges() + 1, false);
        ASSERT_EQ(2u, graph.mark_sink_dummy_edges(&sink_nodes_result));
        EXPECT_EQ(sink_nodes, sink_nodes_result) << graph;
    }
}

TEST(BOSS, MarkDummySourceEdgesSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A'));
        sdsl::bit_vector source_nodes(graph.num_edges() + 1, true);
        source_nodes[0] = false;
        source_nodes[source_nodes.size() - 1] = false;
        sdsl::bit_vector source_nodes_result(graph.num_edges() + 1, false);
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
        sdsl::bit_vector source_nodes(graph.num_edges() + 1, true);
        source_nodes[0] = false;
        source_nodes[1 + 2 + k] = false;
        source_nodes[1 + 2 + 2 * k] = false;
        sdsl::bit_vector source_nodes_result(graph.num_edges() + 1, false);
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
        sdsl::bit_vector source_nodes(graph.num_edges() + 1, true);
        source_nodes[0] = false;
        source_nodes[source_nodes.size() - 1] = false;
        sdsl::bit_vector source_nodes_result(graph.num_edges() + 1, false);
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
        sdsl::bit_vector source_nodes(graph.num_edges() + 1, true);
        source_nodes[0] = false;
        source_nodes[1 + 2 + k] = false;
        source_nodes[1 + 2 + 2 * k] = false;
        sdsl::bit_vector source_nodes_result(graph.num_edges() + 1, false);
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

        BOSS &first = *first_ptr;
        BOSS &second = *second_ptr;

        ASSERT_TRUE(first.equals_internally(second)) << first << second;

        sdsl::bit_vector source_dummy_edges(second.num_edges() + 1, false);
        auto to_remove = second.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(second.num_edges() + 1, false);
        second.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_EQ(0u, std::count(to_remove.begin(), to_remove.end(), true));
        EXPECT_TRUE(first.equals_internally(second)) << first << second;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(second.num_edges() + 1, false);
        auto to_remove = second.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(second.num_edges() + 1, false);
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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

        sdsl::bit_vector source_dummy_edges(dynamic_graph.num_edges() + 1, false);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(&source_dummy_edges, 10);

        sdsl::bit_vector source_dummy_edges_result(dynamic_graph.num_edges() + 1, false);
        dynamic_graph.mark_source_dummy_edges(&source_dummy_edges_result, 1);
        EXPECT_EQ(source_dummy_edges_result, source_dummy_edges);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << redundant_edges;
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
        empty.call_sequences([&](const auto &sequence, const auto &path) {
            ASSERT_EQ(path, empty.map_to_edges(sequence));
            reconstructed.add_sequence(sequence);
        });

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(BOSS, CallUnitigsEmptyGraph) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 30; ++k) {
            BOSS empty(k);

            BOSS reconstructed(k);
            std::mutex seq_mutex;
            empty.call_unitigs([&](const auto &sequence, const auto &path) {
                ASSERT_EQ(path, empty.map_to_edges(sequence));
                std::unique_lock<std::mutex> lock(seq_mutex);
                reconstructed.add_sequence(sequence);
            }, num_threads);

            EXPECT_EQ(empty, reconstructed);
        }
    }
}

TEST(BOSS, CallSequencesRowDiff_EmptyGraph) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 30; ++k) {
            BOSS empty(k);

            sdsl::bit_vector terminal(empty.get_last().size(), false);
            empty.row_diff_traverse(num_threads, 1, empty.get_last(), &terminal);
            ASSERT_EQ(sdsl::bit_vector(2, false), terminal)
                << "Empty graph must have no anchors";
        }
    }
}

TEST(BOSS, CallPathsOneLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSS graph(k);

            ASSERT_EQ(1u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads);
            graph.call_sequences([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallUnitigsOneLoop) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSS graph(k);

            ASSERT_EQ(1u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads,
                             true /* unitigs */);
            graph.call_unitigs([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallPathsTwoLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A') });
            BOSS graph(&constructor);

            ASSERT_EQ(2u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads);
            graph.call_sequences([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallUnitigsTwoLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A') });
            BOSS graph(&constructor);

            ASSERT_EQ(2u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads,
                             true /* unitigs */);
            graph.call_unitigs([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallSequenceRowDiff_TwoLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A') });
            BOSS graph(&constructor);

            ASSERT_EQ(2u, graph.num_edges());

            sdsl::bit_vector terminal(graph.get_last().size(), false);
            graph.row_diff_traverse(num_threads, 1, graph.get_last(), &terminal);
            ASSERT_EQ(graph.num_edges() + 1, terminal.size());
            ASSERT_EQ(sdsl::bit_vector({ 0, 0, 1 }), terminal);
        }
    }
}

TEST(BOSS, CallUnitigsTwoBigLoops) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 10;
        std::vector<std::string> sequences {
            "ATCGGAAGAGCACACGTCTG" "AACTCCAGACA" "CTAAGGCATCTCGTATGCATCGGAAGAGC",
            "GTGAGGCGTCATGCATGCAT" "TGTCTGGAGTT" "TCGTAGCGGCGGCTAGTGCGCGTAGTGAGGCGTCA"
        };
        BOSSConstructor constructor(k);
        constructor.add_sequences(std::vector<std::string>(sequences));
        BOSS graph(&constructor);

        std::atomic<size_t> num_sequences = 0;
        std::atomic<size_t> num_kmers = 0;

        graph.call_unitigs([&](const auto &seq, const auto &path) {
            ASSERT_EQ(path, graph.map_to_edges(seq));
            num_sequences++;
            num_kmers += path.size();
        }, num_threads);

        EXPECT_EQ(2, num_sequences);
        EXPECT_EQ(sequences[0].size() - k - 1 + sequences[1].size() - k - 1,
                  num_kmers);

#if ! _PROTEIN_GRAPH
        num_sequences = 0;
        num_kmers = 0;

        graph.call_unitigs([&](const auto &seq, const auto &path) {
            ASSERT_EQ(path, graph.map_to_edges(seq)) << seq;
            num_sequences++;
            num_kmers += path.size();
        }, num_threads, 0, true);

        // There is no guarantee on the number of contigs in the primary mode
        // EXPECT_EQ(2, num_sequences);
        EXPECT_EQ(sequences[0].size() - k - 1 + sequences[1].size() - k - 1 - 1,
                  num_kmers);
#endif
    }
}

TEST(BOSS, CallSequenceRowDiff_TwoBigLoops) {
    for (size_t num_threads : { 1, 4 }) {
        size_t k = 10;
        std::vector<std::string> sequences {
                "ATCGGAAGAGCACACGTCTG" "AACTCCAGACA" "CTAAGGCATCTCGTATGCATCGGAAGAGC",
                "GTGAGGCGTCATGCATGCAT" "TGTCTGGAGTT" "TCGTAGCGGCGGCTAGTGCGCGTAGTGAGGCGTCA"
        };
        BOSSConstructor constructor(k);
        constructor.add_sequences(std::vector<std::string>(sequences));
        BOSS graph(&constructor);

        sdsl::bit_vector terminal(graph.get_last().size(), false);
        graph.row_diff_traverse(num_threads, 100, graph.get_last(), &terminal);
        ASSERT_EQ(graph.num_edges() + 1, terminal.size());
        ASSERT_EQ(2, std::accumulate(terminal.begin() + 1, terminal.end(), 0U));

        graph.row_diff_traverse(num_threads, 1, graph.get_last(), &terminal);
        ASSERT_EQ(graph.num_edges() + 1, terminal.size());
        ASSERT_EQ(104, std::accumulate(terminal.begin() + 1, terminal.end(), 0U));
    }
}

TEST(BOSS, CallPathsFourLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'G'),
                                        std::string(100, 'C') });
            BOSS graph(&constructor);

            ASSERT_EQ(4u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads);
            graph.call_sequences([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallUnitigsFourLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'G'),
                                        std::string(100, 'C') });
            BOSS graph(&constructor);

            ASSERT_EQ(4u, graph.num_edges());

            std::atomic<size_t> num_paths = 0;
            std::atomic<size_t> num_sequences = 0;

            graph.call_paths([&](const auto &, const auto &) { num_paths++; },
                             num_threads,
                             true /* unitigs */);
            graph.call_unitigs([&](const auto &seq, const auto &path) {
                ASSERT_EQ(path, graph.map_to_edges(seq));
                num_sequences++;
            }, num_threads);

            EXPECT_EQ(graph.num_edges(), num_paths);
            EXPECT_EQ(graph.num_edges() - 1, num_sequences);
        }
    }
}

TEST(BOSS, CallSequenceRowDiff_FourLoops) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 20; ++k) {
            BOSSConstructor constructor(k);
            constructor.add_sequences({ std::string(100, 'A'),
                                        std::string(100, 'G'),
                                        std::string(100, 'C'),
                                        std::string(100, 'T')});
            BOSS graph(&constructor);

            sdsl::bit_vector terminal(graph.get_last().size(), false);
            graph.row_diff_traverse(num_threads, 1, graph.get_last(), &terminal);
            ASSERT_EQ(graph.num_edges() + 1, terminal.size());
            ASSERT_EQ(4, std::accumulate(terminal.begin() + 1, terminal.end(), 0U));
        }
    }
}

TEST(BOSS, CallSequenceRowDiff_FourPaths) {
    constexpr size_t k = 5;
    BOSSConstructor constructor(k);
    std::vector<std::string> sequences
            = { "ATCGGAAGA", "TTTAAACCCGGG", "ATACAGCTCGCT", "AAAAAA" };
    constructor.add_sequences(std::vector(sequences.begin(), sequences.end()));
    BOSS graph(&constructor);
    for (size_t num_threads : { 1, 4 }) {
        sdsl::bit_vector terminal(graph.get_last().size(), false);
        graph.row_diff_traverse(num_threads, 20, graph.get_last(), &terminal);
        ASSERT_EQ(graph.num_edges() + 1, terminal.size());
        ASSERT_EQ(4, std::accumulate(terminal.begin() + 1, terminal.end(), 0U));

        graph.row_diff_traverse(num_threads, 1, graph.get_last(), &terminal);
        ASSERT_EQ(graph.num_edges() + 1, terminal.size());
        ASSERT_EQ(19, std::accumulate(terminal.begin() + 1, terminal.end(), 0U));
    }
}

TEST(BOSS, CallPaths) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 10; ++k) {
            std::mutex seq_mutex;
            {
                BOSS graph(k);
                graph.add_sequence("AAACACTAG", true);
                graph.add_sequence("AACGACATG", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_sequences([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AGACACTGA", true);
                graph.add_sequence("GACTACGTA", true);
                graph.add_sequence("ACTAACGTA", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_sequences([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AGACACAGT", true);
                graph.add_sequence("GACTTGCAG", true);
                graph.add_sequence("ACTAGTCAG", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_sequences([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AAACTCGTAGC", true);
                graph.add_sequence("AAATGCGTAGC", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_sequences([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AAACT", false);
                graph.add_sequence("AAATG", false);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_sequences([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
        }
    }
}

TEST(BOSS, CallUnitigs) {
    for (size_t num_threads : { 1, 4 }) {
        for (size_t k = 1; k < 10; ++k) {
            std::mutex seq_mutex;
            {
                BOSS graph(k);
                graph.add_sequence("AAACACTAG", true);
                graph.add_sequence("AACGACATG", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_unitigs([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AGACACTGA", true);
                graph.add_sequence("GACTACGTA", true);
                graph.add_sequence("ACTAACGTA", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_unitigs([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AGACACAGT", true);
                graph.add_sequence("GACTTGCAG", true);
                graph.add_sequence("ACTAGTCAG", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_unitigs([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AAACTCGTAGC", true);
                graph.add_sequence("AAATGCGTAGC", true);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_unitigs([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
            {
                BOSS graph(k);
                graph.add_sequence("AAACT", false);
                graph.add_sequence("AAATG", false);
                graph.switch_state(BOSS::State::STAT);

                BOSS reconstructed(k);

                graph.call_unitigs([&](const auto &sequence, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(sequence));
                    std::unique_lock<std::mutex> lock(seq_mutex);
                    reconstructed.add_sequence(sequence);
                }, num_threads);

                EXPECT_EQ(graph, reconstructed);
            }
        }
    }
}

TEST(BOSS, CallUnitigs1) {
    for (size_t num_threads : { 1, 4 }) {
        BOSSConstructor constructor(3);
        constructor.add_sequences(std::vector<std::string> {
            "ACTAGCTAGCTAGCTAGCTAGC",
            "ACTCT"
        });
        BOSS graph(&constructor);

        std::multiset<std::string> contigs {
            "$$$$",
            "$$$ACT",
            "ACTCT$",
            "ACTA",
            "CTAGCTA",
        };

        std::multiset<std::string> obs_contigs;
        std::mutex seq_mutex;
        graph.call_paths(
            [&](const auto &, const auto &seq) {
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_contigs.insert(graph.decode(seq));
            },
            num_threads,
            true // unitigs
        );

        EXPECT_EQ(contigs, obs_contigs) << graph;
    }
}

TEST(BOSS, CallUnitigsDisconnected1) {
    for (size_t num_threads : { 1, 4 }) {
        BOSSConstructor constructor(3);
        constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                    "ACTCT",
                                    "ATCATCATCATCATCATCAT" });
        BOSS graph(&constructor);

        std::multiset<std::string> contigs {
            "$$$$",
            "$$$ACT",
            "ACTCT$",
            "ACTA",
            "CTAGCTA",
            "TCATCA",
        };

        std::multiset<std::string> obs_contigs;
        std::mutex seq_mutex;

        graph.call_paths(
            [&](const auto &, const auto &seq) {
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_contigs.insert(graph.decode(seq));
            },
            num_threads,
            true // unitigs
        );

        EXPECT_EQ(contigs, obs_contigs);
    }
}

#ifndef _DNA_GRAPH
TEST(BOSS, CallUnitigsDisconnected2) {
    for (size_t num_threads : { 1, 4 }) {
        BOSSConstructor constructor(3);
        constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                    "ACTCT",
                                    "ATCATCATCATCATCATCAT",
                                    "ATNATNATNATNATNATNAT" });
        BOSS graph(&constructor);

        std::multiset<std::string> contigs {
            "$$$$",
            "$$$ACT",
            "ACTCT$",
            "ACTA",
            "CTAGCTA",
            "TCATCA",
            "TNATNA",
        };

        std::multiset<std::string> obs_contigs;
        std::mutex seq_mutex;

        graph.call_paths(
            [&](const auto &, const auto &seq) {
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_contigs.insert(graph.decode(seq));
            },
            num_threads,
            true // unitigs
        );

        EXPECT_EQ(contigs, obs_contigs);
    }
}

TEST(BOSS, CallUnitigsTwoComponents) {
    for (size_t num_threads : { 1, 4 }) {
        BOSSConstructor constructor(3);
        constructor.add_sequences({ "ACTAGCTAGCTAGCTAGCTAGC",
                                    "ACTCT",
                                    "ATCATCATCATCATCATCAT",
                                    "ATNATNATNATNATNATNAT",
                                    "ATCNATCNATCNATCNATCNATCNAT" });
        BOSS graph(&constructor);

        std::multiset<std::string> contigs {
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

        std::multiset<std::string> obs_contigs;
        std::mutex seq_mutex;

        graph.call_paths(
            [&](const auto &, const auto &seq) {
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_contigs.insert(graph.decode(seq));
            },
            num_threads,
            true // unitigs
        );

        EXPECT_EQ(contigs, obs_contigs);
    }
}
#endif

TEST(BOSS, CallUnitigsWithPruning) {
    for (size_t num_threads : { 1, 4 }) {
        BOSSConstructor constructor(4);
        constructor.add_sequences({
            "ACTATAGCTAGTCTATGCGA",
            "ACTATAGCTAGTCTAA",
            "ACTATAGCTA",
            "ACTATAGCTT",
            "ACTATC",
        });
        BOSS graph(&constructor);
        ASSERT_EQ(4u, graph.get_k());

        // BOSS constructs unitigs from its edges as k-mers
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTT", "AGCTAGTCTA", "TCTAA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads);
            EXPECT_EQ(contigs.size(), num_contigs);
            num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 1 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 2 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 3 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 4 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 5 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
        {
            std::set<std::string> contigs {
                "ACTAT", "CTATC", "CTATGCGA", "CTATAGCT", "AGCTAGTCTA", "TCTAT"
            };
            std::atomic<size_t> num_contigs = 0;
            graph.call_unitigs(
                [&](const auto &str, const auto &path) {
                    EXPECT_EQ(path, graph.map_to_edges(str));
                    EXPECT_TRUE(contigs.count(str)) << str;
                    num_contigs++;
                }
            , num_threads, 6 /* max_pruned_dead_end_size */);
            EXPECT_EQ(contigs.size(), num_contigs);
        }
    }
}

TEST(BOSS, CallUnitigsCheckDegree) {
    BOSSConstructor constructor(8);
    constructor.add_sequences({
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGG",
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGC",
        "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATAACTGGAAAA",
        "CCGAACTAGTGAAACTGCAACAGACATACGCTGCTCTGAACTCTAAGGC",
        "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
        "ATCGGAAGAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGC",
        "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAATTGGCTGAAAAGC",
        "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAATCGCTCCACCAA"
    });
    BOSS graph(&constructor);
    ASSERT_EQ(8u, graph.get_k());
    std::mutex seq_mutex;

    std::multiset<std::string> unitigs {
        "AGACAAATCGCTCCACCAA",
        "AGACAAATTGGCTGAAAAGC",
        "ATCGGAAGAGCACACGTCTGAACT",
        "CAGACATAACTGGAAAA",
        "CAGACATACGCTGCTCTGAACT",
        "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATA",
        "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAAT",
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCG",
        "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
        "CCGAACTAGTGAAACTGCAACAGACATA",
        "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAAT",
        "CTGAACTCCAGACACTAAGGCATCTCGTATGC",
        "CTGAACTCTAAGGC",
        "TCTGAACTC"
    };

    for (size_t num_threads : { 1, 4 }) {
        std::multiset<std::string> obs_unitigs;
        graph.call_unitigs(
            [&](const auto &unitig, const auto &path) {
                EXPECT_EQ(path, graph.map_to_edges(unitig));
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_unitigs.insert(unitig);
            },
            num_threads,
            2 // max_pruned_dead_end_size
        );

        EXPECT_EQ(unitigs, obs_unitigs) << num_threads;
    }
}

TEST(BOSS, CallUnitigsWithEdgesCheckDegree) {
    BOSSConstructor constructor(8);
    constructor.add_sequences({
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGG",
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCGC",
        "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATAACTGGAAAA",
        "CCGAACTAGTGAAACTGCAACAGACATACGCTGCTCTGAACTCTAAGGC",
        "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
        "ATCGGAAGAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGC",
        "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAATTGGCTGAAAAGC",
        "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAATCGCTCCACCAA"
    });
    BOSS graph(&constructor);
    ASSERT_EQ(8u, graph.get_k());
    std::mutex seq_mutex;

    std::multiset<std::string> unitigs {
        "AGACAAATCGCTCCACCAA",
        "AGACAAATTGGCTGAAAAGC",
        "ATCGGAAGAGCACACGTCTGAACT",
        "CAGACATAACTGGAAAA",
        "CAGACATACGCTGCTCTGAACT",
        "CCAAAATGAAACCTTCAGTTTTAACTCTTAATCAGACATA",
        "CCAGAGTCTCGTTCGTTATCGGAATTAACCAGACAAAT",
        "CCAGGGTGTGCTTGTCAAAGAGATATTCCGCCAAGCCAGATTCGGGCG",
        "CCAGGTGCAGGGTGGACTCTTTCTGGATGTTGTAGTCAGACAGGGTGCG",
        "CCGAACTAGTGAAACTGCAACAGACATA",
        "CGGAGGGAAAAATATTTACACAGAGTAGGAGACAAAT",
        "CTGAACTCCAGACACTAAGGCATCTCGTATGC",
        "CTGAACTCTAAGGC",
        "TCTGAACTC"
    };

    for (size_t num_threads : { 1, 4 }) {
        std::multiset<std::string> obs_unitigs;
        graph.call_unitigs(
            [&](const auto &unitig, const auto &path) {
                EXPECT_EQ(path, graph.map_to_edges(unitig));
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_unitigs.insert(unitig);
            },
            num_threads,
            2 // max_pruned_dead_end_size
        );

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

TEST(BOSS, CallUnitigsIndegreeFirstNodeIsZero) {
    BOSSConstructor constructor(30);
    constructor.add_sequences({
        "AGAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAATTAG"
    });
    BOSS graph(&constructor);
    ASSERT_EQ(30u, graph.get_k());
    std::mutex seq_mutex;

    std::multiset<std::string> unitigs {
        "GAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAAT"
    };

    for (size_t num_threads : { 1, 4 }) {
        std::multiset<std::string> obs_unitigs;
        graph.call_unitigs(
            [&](const auto &unitig, const auto &path) {
                EXPECT_EQ(path, graph.map_to_edges(unitig));
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_unitigs.insert(unitig);
            },
            num_threads,
            2 // max_pruned_dead_end_size
        );

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

TEST(BOSS, CallUnitigsWithEdgesIndegreeFirstNodeIsZero) {
    BOSSConstructor constructor(30);
    constructor.add_sequences({
        "AGAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAATTAG"
    });
    BOSS graph(&constructor);
    ASSERT_EQ(30u, graph.get_k());
    std::mutex seq_mutex;

    std::multiset<std::string> unitigs {
        "GAAACCCCGTCTCTACTAAAAATACAAAATTAGCCGGGAGTGGTGGCG",
        "AGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCAGGTGTGGTGAC",
        "GCCTGACCAGCATGGTGAAACCCCGTCTCTACTAAAAATACAAAAT"
    };

    for (size_t num_threads : { 1, 4 }) {
        std::multiset<std::string> obs_unitigs;
        graph.call_unitigs(
            [&](const auto &unitig, const auto &path) {
                EXPECT_EQ(path, graph.map_to_edges(unitig));
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs_unitigs.insert(unitig);
            },
            num_threads,
            2 // max_pruned_dead_end_size
        );

        EXPECT_EQ(unitigs, obs_unitigs);
    }
}

TEST(BOSS, CallUnitigsMasked) {
    size_t k = 3;
    // TTGC      GCACGGGTC
    //      TGCA
    // ATGC      GCAGTGGTC
    std::vector<std::string> sequences { "TTGCACGGGTC", "ATGCAGTGGTC" };
    BOSSConstructor constructor(k);
    constructor.add_sequences(std::vector<std::string>(sequences));
    BOSS graph(&constructor);
    std::mutex seq_mutex;

    sdsl::bit_vector mask_bv(graph.num_edges() + 1, false);
    graph.map_to_edges(
        sequences[0],
        [&](auto edge) { mask_bv[edge] = true; }
    );
    bit_vector_stat mask(std::move(mask_bv));

    std::unordered_multiset<std::string> ref = { "TTGCACGGGTC" };

    for (size_t num_threads : { 1, 4 }) {
        std::unordered_multiset<std::string> obs;

        graph.call_unitigs(
            [&](const auto &unitig, const auto &path) {
                EXPECT_EQ(path, graph.map_to_edges(unitig));
                std::unique_lock<std::mutex> lock(seq_mutex);
                obs.insert(unitig);
            },
            num_threads,
            0, // max_pruned_dead_end_size
            false, // kmers_in_single_form
            &mask // subgraph_mask
        );

        EXPECT_EQ(obs, ref);
    }
}

#if ! _PROTEIN_GRAPH
TEST(BOSS, CallUnitigsSingleKmer) {
    size_t k = 6;
    std::vector<std::string> sequences {
        "ATCGGAAGAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGCATCGGAA",
        "GGGGGGGTGCTCTTTTTTT"
    };
    BOSSConstructor constructor(k);
    constructor.add_sequences(std::move(sequences));
    BOSS graph(&constructor);
    graph.prune_and_mark_all_dummy_edges(1);
    std::multiset<std::string> unitigs {
        "GGGGGGTGCTCTTTTTT",
        "GGGGGGG",
        "GAGCACACGTCTGAACTCCAGACACTAAGGCATCTCGTATGCATCGGAAGAGC",
        "TTTTTTT"
    };
    std::mutex seq_mutex;
    for (size_t num_threads : { 1, 4 }) {
        std::multiset<std::string> obs_unitigs;
        graph.call_unitigs([&](auto&& seq, auto&&) {
            std::lock_guard<std::mutex> lock(seq_mutex);
            obs_unitigs.emplace(seq);
        }, num_threads, 0, true);
        EXPECT_EQ(unitigs, obs_unitigs);
    }
}
#endif

template <class Callback>
void call_edges(const BOSS &boss, Callback callback) {
    boss.call_paths([&](auto&& edges, auto&& path) {
        assert(path.size() == edges.size() + boss.get_k());

        for (size_t i = 0; i < edges.size(); ++i) {
            callback(edges[i],
                     std::vector<BOSS::TAlphabet>(path.begin() + i,
                                                  path.begin() + i + boss.get_k() + 1));
        }
    });
}

TEST(BOSS, CallEdgesEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        BOSS empty(k);

        size_t num_edges = 0;
        call_edges(empty, [&](auto, const auto &edge) {
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
        call_edges(graph, [&](auto, const auto &edge) {
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
        call_edges(graph, [&](auto, const auto &edge) {
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
        call_edges(graph, [&](auto, const auto &edge) {
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
        call_edges(graph, [&](auto, const auto &edge) {
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
        call_edges(graph, [&](auto, const auto &edge) {
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
        graph.switch_state(BOSS::State::DYN);

        graph.add_sequence(std::string(100, 'T'));

        size_t num_edges = 0;
        call_edges(graph, [&](auto, const auto &edge) {
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
        graph.switch_state(BOSS::State::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_edges = 0;
        call_edges(graph, [&](auto edge_idx, const auto &edge) {
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
            EXPECT_EQ(std::string(k + 1, 'A'), sequence) << sequence;
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
            EXPECT_TRUE(std::string(k + 1, 'A') == sequence
                        || std::string(k + 1, 'G') == sequence
                        || std::string(k + 1, 'C') == sequence) << sequence;
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
            EXPECT_TRUE(std::string(k + 1, 'A') == sequence
                        || std::string(k + 1, 'G') == sequence
                        || std::string(k + 1, 'C') == sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(3u, num_kmers) << graph;
    }
}

TEST(BOSS, CallKmersTestPath) {
    for (size_t k = 1; k < 20; ++k) {
        BOSS graph(k);
        graph.add_sequence(std::string(100, 'A') + std::string(k, 'C'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k + 1, num_kmers) << graph;
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
        EXPECT_EQ(k + 1 + k, num_kmers) << graph;
    }
}

TEST(BOSS, CallKmersTestPathDisconnected) {
    for (size_t k = 1; k < 20; ++k) {
        BOSSConstructor constructor(k);
        constructor.add_sequence(std::string(100, 'A'));
        BOSS graph(&constructor);
        graph.switch_state(BOSS::State::DYN);

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
        graph.switch_state(BOSS::State::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers) << graph;
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
#if ! _PROTEIN_GRAPH
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
            auto kmer_str_suffixes = utils::generate_strings("ACGT", i);
            for (size_t j = 0; j < kmer_str_suffixes.size(); ++j) {
                all_kmer_str.push_back(std::string(k - i, '$')
                                        + kmer_str_suffixes[j]);
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

TEST(BOSS, map_to_edges) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<BOSS> graph { new BOSS(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::vector<BOSS::node_index> expected_result
            = { graph::SequenceGraph::npos, graph::SequenceGraph::npos, k + 2, k + 2, k + 2 };

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

} // namespace
