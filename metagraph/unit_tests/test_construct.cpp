#include <stdio.h>
#include <string>
#include <sstream>

#include <zlib.h>
#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define private public

#include "dbg_succinct.hpp"

KSEQ_INIT(gzFile, gzread)

const std::string test_data_dir = "../unit_tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F,
                                 Config::StateType state) {
    graph->switch_state(state);

    std::ostringstream ostr;

    ostr << *graph->last;
    EXPECT_EQ(last, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    ostr << *(graph->W);
    EXPECT_EQ(W, ostr.str()) << "state: " << state;

    ostr.clear();
    ostr.str("");

    for (size_t i = 0; i < graph->F.size(); ++i) {
        ostr << graph->F[i] << " ";
    }
    EXPECT_EQ(F, ostr.str()) << "state: " << state;
}


void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F) {
    test_graph(graph, last, W, F, Config::DYN);
    test_graph(graph, last, W, F, Config::STAT);
}


TEST(Construct, GraphDefaultConstructor) {
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

TEST(DBGSuccinct, EmptyGraph) {
    DBG_succ *graph = new DBG_succ(3);
    graph->construct_succ();
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ");
    delete graph;
}

TEST(DBGSuccinct, SwitchState) {
    DBG_succ *graph = new DBG_succ(3);
    graph->construct_succ();
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::STAT);
    test_graph(graph, "01", "00", "0 1 1 1 1 1 ", Config::DYN);
    delete graph;
}

TEST(DBGSuccinct, AddSequenceFast) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = new DBG_succ(3);
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        graph->add_sequence_fast(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    graph->construct_succ();

    //test graph construction
    test_graph(graph, "00011101101111111111111",
                      "00131124434010141720433",
                      "0 3 11 13 17 22 ");
}

TEST(Construct, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = new DBG_succ(3);
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        graph->add_sequence_fast(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    graph->construct_succ();

    //traversal
    std::vector<size_t> outgoing_edges = { 0, 3, 4, 14, 5, 7, 12, 18, 19, 15, 20, 0,
                                           8, 0, 10, 21, 11, 11, 13, 0, 22, 16, 17 };
    ASSERT_EQ(outgoing_edges.size(), graph->num_edges() + 1);

    EXPECT_EQ(outgoing_edges[1], graph->outgoing(1, DBG_succ::encode('$')));

    for (size_t i = 1; i < graph->W->size(); ++i) {
        //test forward traversal given an output edge label
        if (graph->get_W(i) != DBG_succ::encode('$')) {
            EXPECT_EQ(outgoing_edges[i], graph->outgoing(i, graph->get_W(i)))
                << "Edge index: " << i;
        }

        //test FM index property
        EXPECT_TRUE(graph->get_last(graph->fwd(i)));
        if (graph->get_W(i)) {
            EXPECT_EQ(graph->get_W(i) % graph->alph_size,
                      graph->get_node_last_char(graph->fwd(i)));
        }
    }

    delete graph;
}

TEST(DBGSuccinct, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = new DBG_succ(3);
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        graph->add_sequence_fast(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    graph->construct_succ();

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

TEST(DBGSuccinct, TraversalMergeWithEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinct, TraversalMergeEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        second.add_sequence(std::string(100, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinct, TraversalMergeEqualPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'A'));
        first.merge(second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinct, TraversalMergeTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        first.merge(second);
        EXPECT_EQ(2 * k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(2 * k + 3, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinct, TraversalMergeSinglePathWithTwo) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        first.merge(second);
        EXPECT_EQ(3 * k + 1, first.num_nodes());
        EXPECT_EQ(2 * k + 1, second.num_nodes());
        EXPECT_EQ(3 * k + 4, first.num_edges());
        EXPECT_EQ(2 * k + 3, second.num_edges());
    }
}

TEST(DBGSuccinct, TraversalMergeTwoGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        first.add_sequence(std::string(50, 'C'));
        first.add_sequence(std::string(60, 'G'));
        first.add_sequence("AAAGT");
        second.add_sequence("AAACT", true);
        second.add_sequence("AAATG", true);
        second.add_sequence("ACTGA", true);
        DBG_succ merged(k);
        merged.merge(second);
        merged.merge(first);
        first.merge(second);
        EXPECT_EQ(first, merged);
    }
}

TEST(DBGSuccinct, ParallelMergeTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));

        std::vector<const DBG_succ*> graphs = { &first, &second };
        std::vector<uint64_t> kv;
        std::vector<uint64_t> nv;
        for (size_t i = 0; i < graphs.size(); ++i) {
            kv.push_back(1);
            nv.push_back(graphs[i]->get_W().size());
        }
        DBG_succ *merged = merge::merge(graphs, kv, nv);

        first.merge(second);

        EXPECT_EQ(first, *merged);
    }
}

TEST(DBGSuccinct, ParallelMergeSinglePathWithTwo) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));

        std::vector<const DBG_succ*> graphs = { &first, &second };
        std::vector<uint64_t> kv;
        std::vector<uint64_t> nv;
        for (size_t i = 0; i < graphs.size(); ++i) {
            kv.push_back(1);
            nv.push_back(graphs[i]->get_W().size());
        }
        DBG_succ *merged = merge::merge(graphs, kv, nv);

        first.merge(second);

        EXPECT_EQ(first, *merged);
    }
}

TEST(DBGSuccinct, ParallelMergeThreeGraphs) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        DBG_succ third(k, true);
        first.add_sequence("AAACT", true);
        first.add_sequence("ACTATG", true);
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        third.add_sequence(std::string(60, 'A'));
        third.add_sequence(std::string(60, 'T'));

        std::vector<const DBG_succ*> graphs = { &first, &second, &third };
        std::vector<uint64_t> kv;
        std::vector<uint64_t> nv;
        for (size_t i = 0; i < graphs.size(); ++i) {
            kv.push_back(1);
            nv.push_back(graphs[i]->get_W().size());
        }
        DBG_succ *merged = merge::merge(graphs, kv, nv);

        first.merge(second);
        first.merge(third);

        EXPECT_EQ(first, *merged);
    }
}
