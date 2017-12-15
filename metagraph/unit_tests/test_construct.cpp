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



DBG_succ* init_graph(size_t k = 3) {
    DBG_succ *graph = new DBG_succ(k);
    assert(graph);
    graph->add_sink();
    graph->switch_state(Config::CSTR);
    return graph;
}

void construct_succ(DBG_succ *graph) {
    graph->construct_succ();
    graph->switch_state(Config::DYN);
}

void test_graph(DBG_succ *graph, const std::string &last,
                                 const std::string &W,
                                 const std::string &F) {
    std::ostringstream ostr;
    ostr << *(graph->last);
    EXPECT_EQ(ostr.str(), last);
    ostr.clear();
    ostr.str("");
    ostr << *(graph->W);
    EXPECT_EQ(ostr.str(), W);
    ostr.clear();
    ostr.str("");
    for (size_t i = 0; i < graph->F.size(); ++i) {
        ostr << graph->F[i] << " ";
    }
    EXPECT_EQ(ostr.str(), F);

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
    DBG_succ *graph = init_graph();
    construct_succ(graph);
    graph->switch_state(Config::DYN);
    test_graph(graph, "01111", "06660", "0 1 1 1 1 1 1 ");
    delete graph;
}

TEST(Construct, SmallGraphTraversal) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = init_graph();
    graph->add_sink();
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        graph->add_sequence_fast(read_stream->seq.s);
    }
    graph->construct_succ();
    graph->switch_state(Config::DYN);

    //test graph construction
    test_graph(graph, "01011110111111111111111111111001",
                      "06241463411641812643366666131313013",
                      "0 1 9 11 15 20 20 ");

    //traversal
    std::vector<size_t> outgoing_edges = { 21, 10, 16, 3, 17, 22, 12, 18, 4, 5,
                                           23, 19, 6, 6, 8, 11, 24, 20, 13, 14,
                                           25, 26, 27, 28, 31, 31, 31, 31, 1, 9, 15 };
    assert(outgoing_edges.size() == graph->num_edges());
    std::vector<size_t> incoming_edges = {};
    for (size_t i = 0; i < outgoing_edges.size(); ++i) {
        //test forward traversal given an output edge label
        EXPECT_EQ(outgoing_edges[i], graph->outgoing(i + 1, graph->get_W(i + 1)))
            << "Edge index: " << i + 1;

        //test FM index property
        if (graph->get_W(i + 1) < graph->alph_size) {
            EXPECT_TRUE(graph->get_last(graph->fwd(i + 1)));
            EXPECT_EQ(graph->get_W(i + 1),
                      graph->get_node_last_char(graph->fwd(i + 1)));
        }
    }

    delete graph;
    kseq_destroy(read_stream);
    gzclose(input_p);
}

TEST(DBGSuccinct, Serialization) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = init_graph();
    graph->add_sink();
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        graph->add_sequence_fast(read_stream->seq.s);
    }
    kseq_destroy(read_stream);
    gzclose(input_p);

    graph->construct_succ();
    graph->switch_state(Config::DYN);

    //test graph construction
    test_graph(graph, "01011110111111111111111111111001",
                      "06241463411641812643366666131313013",
                      "0 1 9 11 15 20 20 ");

    graph->serialize(test_dump_basename);

    DBG_succ loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";

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

TEST(DBGSuccinct, MergeWithEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        merge::merge(&first, second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinct, MergeEmpty) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        second.add_sequence(std::string(100, 'A'));
        merge::merge(&first, second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
    }
}

TEST(DBGSuccinct, MergeEqualPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'A'));
        merge::merge(&first, second);
        EXPECT_EQ(k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(k + 2, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinct, MergeTwoPaths) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        merge::merge(&first, second);
        EXPECT_EQ(2 * k + 1, first.num_nodes());
        EXPECT_EQ(k + 1, second.num_nodes());
        EXPECT_EQ(2 * k + 3, first.num_edges());
        EXPECT_EQ(k + 2, second.num_edges());
    }
}

TEST(DBGSuccinct, MergeSinglePathWithTwo) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ first(k, true);
        DBG_succ second(k, true);
        first.add_sequence(std::string(100, 'A'));
        second.add_sequence(std::string(50, 'C'));
        second.add_sequence(std::string(60, 'G'));
        merge::merge(&first, second);
        EXPECT_EQ(3 * k + 1, first.num_nodes());
        EXPECT_EQ(2 * k + 1, second.num_nodes());
        EXPECT_EQ(3 * k + 4, first.num_edges());
        EXPECT_EQ(2 * k + 3, second.num_edges());
    }
}
