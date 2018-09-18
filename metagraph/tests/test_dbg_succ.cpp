#include <stdio.h>
#include <string>
#include <sstream>
#include <mutex>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define protected public
#define private public

#include "dbg_succinct.hpp"
#include "dbg_succinct_construct.hpp"
#include "utils.hpp"

KSEQ_INIT(gzFile, gzread);

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";

double kEps = std::numeric_limits<double>::epsilon();


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::vector<uint64_t> &W,
                                 const std::string &F,
                                 Config::StateType state) {
    Config::StateType old_state = graph->state;
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << *graph->last_;
    EXPECT_EQ(last, ostr.str()) << "state: " << state
                                << ", old state: " << old_state;

    ostr.clear();
    ostr.str("");

    auto W_vector = graph->W_->to_vector();
    EXPECT_EQ(W, std::vector<uint64_t>(W_vector.begin(), W_vector.end()))
        << "state: " << state
        << ", old state: " << old_state;

    ostr.clear();
    ostr.str("");

    for (size_t i = 0; i < graph->F_.size(); ++i) {
        ostr << graph->F_[i] << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state
                             << ", old state: " << old_state;
}


void test_graph(DBG_succ *graph, const std::string &last,
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


TEST(DBGSuccinct, GraphDefaultConstructor) {
    DBG_succ *graph_ = NULL;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ;
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        graph_ = new DBG_succ();
    });
    ASSERT_TRUE(graph_ != NULL);
    delete graph_;

    ASSERT_NO_THROW({
        DBG_succ graph;
    });

    ASSERT_NO_THROW({
        DBG_succ graph();
    });
}

#if _DNA_GRAPH
TEST(DBGSuccinct, EmptyGraph) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", { 0, 0 }, "0 1 1 1 1 1 ");
    delete graph;
}

TEST(DBGSuccinct, SwitchState) {
    DBG_succ *graph = new DBG_succ(3);
    test_graph(graph, "01", { 0, 0 }, "0 1 1 1 1 1 ");
    delete graph;
}

TEST(DBGSuccinct, AddSequenceFast) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    std::vector<std::string> names;

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
        names.emplace_back(read_stream->name.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    EXPECT_EQ(4llu, names.size());
    EXPECT_EQ("1", names[0]);
    EXPECT_EQ("2", names[1]);
    EXPECT_EQ("3", names[2]);
    EXPECT_EQ("4", names[3]);

    DBG_succ *graph = new DBG_succ(&constructor);

    //test graph construction
    test_graph(graph, "00011101101111111111111",
                      { 0, 0, 1, 3, 1, 1, 2, 4,
                        4, 3, 4, 0, 1, 0, 1, 4,
                        1, 7, 2, 0, 4, 3, 3 },
                      "0 3 11 13 17 22 ");
    delete graph;
}
#endif

TEST(DBGSuccinct, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(&constructor);

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    EXPECT_EQ(1u, graph->outgoing(1, DBG_succ::kSentinelCode));

    for (size_t i = 1; i <= graph->num_edges(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBG_succ::kSentinelCode) {
            uint64_t node_idx = graph->rank_last(i - 1) + 1;

            EXPECT_EQ(outgoing_edges[i],
                graph->select_last(graph->outgoing(node_idx, graph->get_W(i)))
            ) << "Edge index: " << i << "\n"
              << "Outgoing: " << graph->outgoing(node_idx, graph->get_W(i)) << "\n"
              << *graph;

            EXPECT_EQ(node_idx,
                graph->incoming(graph->outgoing(node_idx, graph->get_W(i)),
                                graph->get_minus_k_value(i, graph->get_k() - 1).first)
            );
            for (TAlphabet c = 0; c < DBG_succ::alph_size; ++c) {
                uint64_t prev_node = graph->incoming(node_idx, c);
                if (prev_node) {
                    EXPECT_EQ(node_idx,
                        graph->outgoing(prev_node, graph->get_node_last_value(i))
                    );
                }
            }
        }

        //test FM index property
        EXPECT_TRUE(graph->get_last(graph->fwd(i)));
        if (graph->get_W(i)) {
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_value(graph->fwd(i)));
        }
    }

    delete graph;
}

TEST(DBGSuccinct, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);

    KMerDBGSuccConstructor constructor(3);

    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        constructor.add_reads({ read_stream->seq.s });
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    DBG_succ *graph = new DBG_succ(&constructor);

    graph->serialize(test_dump_basename);

    DBG_succ loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";
    EXPECT_FALSE(DBG_succ() == loaded_graph);

    delete graph;
}

TEST(DBGSuccinct, AddSequenceSimplePath) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A'));
        EXPECT_EQ(k + 1, graph.num_nodes());
        EXPECT_EQ(k + 2, graph.num_edges());
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesForClearGraph) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DBG_succ> first_ptr;
        std::unique_ptr<DBG_succ> second_ptr;

        {
            KMerDBGSuccConstructor constructor(k);
            constructor.add_reads({ std::string(100, 'A'),
                                    std::string(100, 'C'),
                                    std::string(100, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
            first_ptr.reset(new DBG_succ(&constructor));
        }

        {
            KMerDBGSuccConstructor constructor(k);
            constructor.add_reads({ std::string(100, 'A'),
                                    std::string(100, 'C'),
                                    std::string(100, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
            second_ptr.reset(new DBG_succ(&constructor));
        }

        auto &first = *first_ptr;
        auto &second = *second_ptr;

        ASSERT_TRUE(first.equals_internally(second)) << first;
        auto to_remove = second.erase_redundant_dummy_edges();
        EXPECT_EQ(0u, std::count(to_remove.begin(), to_remove.end(), true));
        EXPECT_TRUE(first.equals_internally(second)) << first;
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesLinear) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(20, 'A'));
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges();
        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesThreePaths) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'), });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges();

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesFourPaths) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                std::string(20, 'T') + 'A' });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges();

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesFivePaths) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                std::string(20, 'T'),
                                std::string(20, 'N') + 'A', });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T'));
        dynamic_graph.add_sequence(std::string(20, 'N') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges();

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdges) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                "ACATCTAGTAGTCGATCGTACG",
                                "ATTAGTAGTAGTAGTGATGTAG", });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence("ACATCTAGTAGTCGATCGTACG");
        dynamic_graph.add_sequence("ATTAGTAGTAGTAGTGATGTAG");

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges();
        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesForClearGraphParallel) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DBG_succ> first_ptr;
        std::unique_ptr<DBG_succ> second_ptr;

        {
            KMerDBGSuccConstructor constructor(k);
            constructor.add_reads({ std::string(100, 'A'),
                                    std::string(100, 'C'),
                                    std::string(100, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
            first_ptr.reset(new DBG_succ(&constructor));
        }

        {
            KMerDBGSuccConstructor constructor(k);
            constructor.add_reads({ std::string(100, 'A'),
                                    std::string(100, 'C'),
                                    std::string(100, 'G'),
                                    "ACATCTAGTAGTCGATCGTACG",
                                    "ATTAGTAGTAGTAGTGATGTAG", });
            second_ptr.reset(new DBG_succ(&constructor));
        }

        auto &first = *first_ptr;
        auto &second = *second_ptr;

        ASSERT_TRUE(first.equals_internally(second)) << first;
        auto to_remove = second.erase_redundant_dummy_edges(10);
        EXPECT_EQ(0u, std::count(to_remove.begin(), to_remove.end(), true));
        EXPECT_TRUE(first.equals_internally(second)) << first;
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesLinearParallel) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(20, 'A'));
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(10);
        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesThreePathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'), });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(10);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesFourPathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                std::string(20, 'T') + 'A' });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(10);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesFivePathsParallel) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                std::string(20, 'T'),
                                std::string(20, 'N') + 'A', });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence(std::string(20, 'T'));
        dynamic_graph.add_sequence(std::string(20, 'N') + 'A');

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);

        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(10);

        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, RemoveDummyEdgesParallel) {
    for (size_t k = 1; k < 10; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(20, 'A'),
                                std::string(20, 'C'),
                                std::string(20, 'G'),
                                "ACATCTAGTAGTCGATCGTACG",
                                "ATTAGTAGTAGTAGTGATGTAG", });
        DBG_succ graph(&constructor);

        DBG_succ dynamic_graph(k);
        dynamic_graph.add_sequence(std::string(20, 'A'));
        dynamic_graph.add_sequence(std::string(20, 'C'));
        dynamic_graph.add_sequence(std::string(20, 'G'));
        dynamic_graph.add_sequence("ACATCTAGTAGTCGATCGTACG");
        dynamic_graph.add_sequence("ATTAGTAGTAGTAGTGATGTAG");

        ASSERT_FALSE(graph.equals_internally(dynamic_graph));
        ASSERT_EQ(graph, dynamic_graph);
        auto redundant_edges = dynamic_graph.erase_redundant_dummy_edges(10);
        EXPECT_TRUE(graph.equals_internally(dynamic_graph))
            << "Clear graph\n" << graph
            << "Cleaned up graph\n" << dynamic_graph
            << "Removed edges\n" << bit_vector_stat(redundant_edges);
        EXPECT_EQ(graph, dynamic_graph);
    }
}

TEST(DBGSuccinct, AddSequenceBugRevealingTestcase) {
    DBG_succ graph(1);
    graph.add_sequence("CTGAG", false);
}

TEST(DBGSuccinct, NonASCIIStrings) {
    KMerDBGSuccConstructor constructor_first(5);
    constructor_first.add_reads({ "АСАСАСАСАСАСА",
                                  "плохая строка",
                                  "АСАСАСАСАСАСА" });
    DBG_succ graph(&constructor_first);
    ASSERT_EQ(2u, graph.num_edges()) << graph;
}

TEST(DBGSuccinct, AddSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("CAAC");
        graph.add_sequence("GAAC");
        EXPECT_EQ(11u, graph.num_nodes());
        EXPECT_EQ(14u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC");
        graph.add_sequence("AACG");
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC");
        graph.add_sequence("GACT");
        graph.add_sequence("ACTA");
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequence) {
    {
        DBG_succ graph(3);
        graph.add_sequence("AAAC", true);
        graph.add_sequence("AACG", true);
        EXPECT_EQ(6u, graph.num_nodes());
        EXPECT_EQ(7u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(7u, graph.num_nodes());
        EXPECT_EQ(8u, graph.num_edges());
    }
    {
        DBG_succ graph(4);
        graph.add_sequence("AGAC", true);
        graph.add_sequence("GACT", true);
        graph.add_sequence("ACTA", true);
        EXPECT_EQ(12u, graph.num_nodes());
        EXPECT_EQ(15u, graph.num_edges());
    }
    {
        DBG_succ graph(3);
        graph.add_sequence("AAACT", true);
        graph.add_sequence("AAATG", true);
        EXPECT_EQ(8u, graph.num_nodes());
        EXPECT_EQ(10u, graph.num_edges());
    }
}

TEST(DBGSuccinct, AppendSequenceAnyKmerSize) {
    for (size_t k = 1; k < 10; ++k) {
        {
            DBG_succ graph(k);
            graph.add_sequence("AAAC", true);
            graph.add_sequence("AACG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGAC", true);
            graph.add_sequence("GACT", true);
            graph.add_sequence("ACTA", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", true);
            graph.add_sequence("AAATG", true);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
        }
    }
}

TEST(DBGSuccinct, CallSimplePathsEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        DBG_succ empty(k);
        DBG_succ reconstructed(k);

        empty.call_sequences([&](const auto &sequence) {
            reconstructed.add_sequence(sequence);
        });

        EXPECT_EQ(empty, reconstructed);
    }
}

TEST(DBGSuccinct, CallSimplePathsOneLoop) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);

        ASSERT_EQ(1u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; });
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(DBGSuccinct, CallSimplePathsTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; });
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(DBGSuccinct, CallSimplePathsFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A'),
                                std::string(100, 'G'),
                                std::string(100, 'C') });
        DBG_succ graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_paths = 0;
        size_t num_sequences = 0;

        graph.call_paths([&](const auto &, const auto &) { num_paths++; });
        graph.call_sequences([&](const auto &) { num_sequences++; });

        EXPECT_EQ(graph.num_edges(), num_paths);
        EXPECT_EQ(graph.num_edges() - 1, num_sequences);
    }
}

TEST(DBGSuccinct, CallSimplePaths) {
    for (size_t k = 1; k < 10; ++k) {
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACACTAG", true);
            graph.add_sequence("AACGACATG", true);
            graph.switch_state(Config::STAT);

            DBG_succ reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGACACTGA", true);
            graph.add_sequence("GACTACGTA", true);
            graph.add_sequence("ACTAACGTA", true);
            graph.switch_state(Config::STAT);

            DBG_succ reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AGACACAGT", true);
            graph.add_sequence("GACTTGCAG", true);
            graph.add_sequence("ACTAGTCAG", true);
            graph.switch_state(Config::STAT);

            DBG_succ reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACTCGTAGC", true);
            graph.add_sequence("AAATGCGTAGC", true);
            graph.switch_state(Config::STAT);

            DBG_succ reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
        {
            DBG_succ graph(k);
            graph.add_sequence("AAACT", false);
            graph.add_sequence("AAATG", false);
            graph.switch_state(Config::STAT);

            DBG_succ reconstructed(k);

            graph.call_sequences([&](const auto &sequence) {
                reconstructed.add_sequence(sequence);
            });

            EXPECT_EQ(graph, reconstructed);
        }
    }
}

TEST(DBGSuccinct, CallEdgesEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        DBG_succ empty(k);

        size_t num_edges = 0;
        empty.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(empty.get_k() + 1, edge.size()) << empty;
            num_edges++;
        });

        EXPECT_EQ(empty.num_edges(), num_edges);
    }
}

TEST(DBGSuccinct, CallEdgesTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(DBGSuccinct, CallEdgesFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A'),
                                std::string(100, 'G'),
                                std::string(100, 'C') });
        DBG_succ graph(&constructor);

        ASSERT_EQ(4u, graph.num_edges());

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(DBGSuccinct, CallEdgesFourLoopsDynamic) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
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

TEST(DBGSuccinct, CallEdgesTestPath) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A') + std::string(k, 'C'));

        size_t num_edges = 0;
        graph.call_edges([&](auto, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges);
    }
}

TEST(DBGSuccinct, CallEdgesTestPathACA) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
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

TEST(DBGSuccinct, CallEdgesTestPathDisconnected) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(100, 'A'));
        DBG_succ graph(&constructor);
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

TEST(DBGSuccinct, CallEdgesTestPathDisconnected2) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(100, 'G'));
        DBG_succ graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_edges = 0;
        graph.call_edges([&](auto edge_idx, const auto &edge) {
            EXPECT_EQ(graph.get_k() + 1, edge.size()) << graph;
            EXPECT_EQ(graph.get_node_seq(edge_idx),
                      std::vector<TAlphabet>(edge.begin(), edge.end() - 1))
                << edge_idx << "\n" << graph;
            num_edges++;
        });
        EXPECT_EQ(graph.num_edges(), num_edges) << graph;
    }
}

TEST(DBGSuccinct, CallKmersEmptyGraph) {
    for (size_t k = 1; k < 30; ++k) {
        DBG_succ empty(k);

        size_t num_kmers = 0;
        empty.call_kmers([&](auto, const auto &sequence) {
            EXPECT_FALSE(true) << sequence;
            num_kmers++;
        });

        EXPECT_EQ(0u, num_kmers);
    }
}

TEST(DBGSuccinct, CallKmersTwoLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A') });
        DBG_succ graph(&constructor);

        ASSERT_EQ(2u, graph.num_edges());

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto &sequence) {
            EXPECT_EQ(std::string(k, 'A'), sequence) << sequence;
            num_kmers++;
        });
        EXPECT_EQ(1u, num_kmers);
    }
}

TEST(DBGSuccinct, CallKmersFourLoops) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_reads({ std::string(100, 'A'),
                                std::string(100, 'G'),
                                std::string(100, 'C') });
        DBG_succ graph(&constructor);

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

TEST(DBGSuccinct, CallKmersFourLoopsDynamic) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
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

TEST(DBGSuccinct, CallKmersTestPath) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A') + std::string(k, 'C'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(k + 1, num_kmers);
    }
}

TEST(DBGSuccinct, CallKmersTestPathACA) {
    for (size_t k = 1; k < 20; ++k) {
        DBG_succ graph(k);
        graph.add_sequence(std::string(100, 'A')
                            + std::string(k, 'C')
                            + std::string(100, 'A'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2 * k, num_kmers);
    }
}

TEST(DBGSuccinct, CallKmersTestPathDisconnected) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(100, 'A'));
        DBG_succ graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(100, 'T'));

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(2u, num_kmers) << graph;
    }
}

TEST(DBGSuccinct, CallKmersTestPathDisconnected2) {
    for (size_t k = 1; k < 20; ++k) {
        KMerDBGSuccConstructor constructor(k);
        constructor.add_read(std::string(100, 'G'));
        DBG_succ graph(&constructor);
        graph.switch_state(Config::DYN);

        graph.add_sequence(std::string(k, 'A') + "T");

        size_t num_kmers = 0;
        graph.call_kmers([&](auto, const auto&) { num_kmers++; });
        EXPECT_EQ(3u, num_kmers) << graph;
    }
}

void test_pred_kmer(const DBG_succ &graph,
                    const std::string &kmer_s,
                    uint64_t expected_idx) {
    std::vector<TAlphabet> kmer(kmer_s.size());
    std::transform(kmer_s.begin(), kmer_s.end(), kmer.begin(),
                   [](char c) {
                       return c == DBG_succ::kSentinel
                                   ? DBG_succ::kSentinelCode
                                   : DBG_succ::encode(c);
                   });
    EXPECT_EQ(expected_idx, graph.select_last(graph.pred_kmer(kmer))) << kmer_s << std::endl << graph;
}

TEST(DBGSuccinct, PredKmer) {
    {
        DBG_succ graph(5);

        test_pred_kmer(graph, "ACGCG", 1);
        test_pred_kmer(graph, "$$$$A", 1);
        test_pred_kmer(graph, "TTTTT", 1);
        test_pred_kmer(graph, "NNNNN", 1);
        test_pred_kmer(graph, "$$$$$", 1);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAAAA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("ACACA");

        test_pred_kmer(graph, "ACGCG", 7);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 7);
        test_pred_kmer(graph, "NNNNN", 7);
        test_pred_kmer(graph, "$$$$$", 2);
    }
#ifndef _PROTEIN_GRAPH
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGTAGTATGTAGC");

        test_pred_kmer(graph, "ACGCG", 13);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 18);
        test_pred_kmer(graph, "NNNNN", 18);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(5);
        graph.add_sequence("AAACGAAGGAAGTACGC");

        test_pred_kmer(graph, "ACGCG", 17);
        test_pred_kmer(graph, "$$$$A", 3);
        test_pred_kmer(graph, "TTTTT", 19);
        test_pred_kmer(graph, "NNNNN", 19);
        test_pred_kmer(graph, "$$$$$", 2);
    }
    {
        DBG_succ graph(2);
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

TEST(DBGSuccinct, PredKmerRandomTest) {
    for (size_t k = 1; k < 8; ++k) {
        DBG_succ graph(k);

        for (size_t p = 0; p < 10; ++p) {
            size_t length = rand() % 400;
            std::string sequence(length, 'A');

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
            }
            graph.add_sequence(sequence, false);

            for (size_t s = 0; s < sequence.size(); ++s) {
                sequence[s] = DBG_succ::alphabet[1 + rand() % 4];
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
            std::vector<TAlphabet> kmer = DBG_succ::encode(kmer_str);

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

TEST(DBGSuccinct, FindSequence) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ *graph = new DBG_succ(k);

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

TEST(DBGSuccinct, KmerMappingMode) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k);

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

TEST(DBGSuccinct, Traversals) {
    for (size_t k = 1; k < 10; ++k) {
        SequenceGraph *graph = new DBG_succ(k);

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        uint64_t it = 0;
        graph->map_to_nodes(std::string(k, 'A'), [&](auto i) { it = i; });
        ASSERT_EQ(k + 1, it);
        EXPECT_EQ(it, graph->traverse(it, 'A'));
        EXPECT_EQ(it + 1, graph->traverse(it, 'C'));
        EXPECT_EQ(it, graph->traverse_back(it + 1, 'A'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse(it, 'G'));
        EXPECT_EQ(DBG_succ::npos, graph->traverse_back(it + 1, 'G'));

        delete graph;
    }
}

TEST(DBGSuccinct, map_to_nodes) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DBG_succ> graph { new DBG_succ(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::vector<uint64_t> expected_result {
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

TEST(DBGSuccinct, map_to_edges) {
    for (size_t k = 1; k < 10; ++k) {
        std::unique_ptr<DBG_succ> graph { new DBG_succ(k) };

        graph->add_sequence(std::string(100, 'A') + std::string(100, 'C'));

        std::vector<uint64_t> expected_result {
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
