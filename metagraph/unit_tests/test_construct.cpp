#include <stdio.h>
#include <string>
#include <sstream>

#include <htslib/kseq.h>
#include "gtest/gtest.h"

#define private public

#include "dbg_succinct.hpp"

KSEQ_INIT(gzFile, gzread)

const std::string test_data_dir = "../unit_tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";
const std::string test_dump_basename = test_data_dir + "/graph_dump_test";



DBG_succ* init_graph(size_t k = 3) {
    DBG_succ *graph = new DBG_succ(k, NULL);
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
                                 const std::string &F, const size_t &p) {
    std::ostringstream ostr;
    ostr << *(graph->last);
    EXPECT_EQ(ostr.str(), last);
    ostr.clear();
    ostr.str("");
    ostr << *(graph->W);
    EXPECT_EQ(ostr.str(), W);
    ostr.clear();
    ostr.str("");
    EXPECT_EQ(graph->p_, p);
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
    test_graph(graph, "01111", "06660", "0 1 1 1 1 1 1 ", 4u);
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
                      "0 1 9 11 15 20 20 ",
                      29u);

    //traversal
    std::vector<size_t> outgoing_edges = { 21, 10, 16, 3, 17, 22, 12, 18, 4, 5,
                                           23, 19, 6, 6, 8, 11, 24, 20, 13, 14,
                                           25, 26, 27, 28, 31, 31, 31, 31, 1, 9, 15 };
    assert(outgoing_edges.size() == graph->num_edges());
    std::vector<size_t> incoming_edges = {};
    for (size_t i = 0; i < outgoing_edges.size(); ++i) {
        //test forward traversal given an output edge label
        EXPECT_EQ(graph->outgoing(i + 1, graph->get_W(i + 1)), outgoing_edges[i]);
        //test that there is only one terminus
        auto sink_rank = graph->get_equal_node_range(graph->p_);
        if (i + 1 < sink_rank.first || i + 1 > sink_rank.second) {
            EXPECT_EQ(graph->outgoing(i + 1, 0), 0u);
        } else {
            EXPECT_EQ(graph->outgoing(i + 1, 0), 1u);
        }
        //test FM index property
        if (graph->get_W(i + 1) < graph->alph_size) {
            EXPECT_TRUE(graph->get_last(graph->fwd(i + 1)));
            EXPECT_EQ(graph->get_W(i + 1),
                      graph->get_node_end_value(graph->fwd(i + 1)));
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
                      "0 1 9 11 15 20 20 ",
                      29u);

    graph->serialize(test_dump_basename);

    DBG_succ loaded_graph;
    ASSERT_TRUE(loaded_graph.load(test_dump_basename)) << "Can't load the graph";
    EXPECT_EQ(*graph, loaded_graph) << "Loaded graph differs";

    delete graph;
}

TEST(DBGSuccinct, AddSequence) {
    for (size_t k = 1; k < 10; ++k) {
        DBG_succ graph(k, true);
        graph.add_sequence(std::string(100, 'A'));
        EXPECT_EQ(3 * k, graph.num_nodes());
        EXPECT_EQ(3 * k + 2, graph.num_edges());
    }
}
