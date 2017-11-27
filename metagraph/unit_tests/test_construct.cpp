#include <stdio.h>
#include <string>
#include <sstream>

#include <htslib/kseq.h>
#include "gtest/gtest.h"

#include "construct.hpp"
#include "dbg_succinct.hpp"

KSEQ_INIT(gzFile, gzread)

const std::string test_data_dir = "../unit_tests/data";
const std::string test_fasta = test_data_dir + "/test_construct.fa";


DBG_succ* init_graph() {
    Config *config = new Config();
    DBG_succ *graph = new DBG_succ(config->k, config);
    delete config;
    assert(graph);
    construct::add_sink(graph);
    graph->switch_state(Config::CSTR);
    return graph;
}

void construct_succ(DBG_succ *graph) {
    construct::construct_succ(graph);
    graph->switch_state(Config::DYN);
}

void test_graph(DBG_succ *graph, const std::string& last, const std::string& W, const std::string& F, const size_t& p) {
    std::ostringstream ostr;
    ostr << *(graph->last);
    EXPECT_EQ(ostr.str(), last);
    ostr.clear();
    ostr.str("");
    ostr << *(graph->W);
    EXPECT_EQ(ostr.str(), W);
    ostr.clear();
    ostr.str("");
    EXPECT_EQ(graph->p, p);
    for (size_t i = 0; i < graph->F.size(); ++i) {
        ostr << graph->F[i] << " ";
    }
    EXPECT_EQ(ostr.str(), F);

}

TEST(Construct, EmptyGraph) {
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
    construct::add_sink(graph);
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        construct::add_seq_fast(graph, read_stream->seq.s, std::string(read_stream->name.s, read_stream->name.l), true);
    }
    construct::construct_succ(graph);
    graph->switch_state(Config::DYN);

    //test graph construction
    test_graph(graph, "01011110111111111111111111111001", "06241463411641812643366666131313013", "0 1 9 11 15 20 20 ", 29u);

    //traversal
    std::vector<size_t> outgoing_edges = {21, 10, 16, 3, 17, 22, 12, 18, 4, 5, 23, 19, 6, 6, 8, 11, 24, 20, 13, 14, 25, 26, 27, 28, 31, 31, 31, 31, 1, 9, 15};
    assert(outgoing_edges.size() == graph->get_edge_count());
    std::vector<size_t> incoming_edges = {};
    for (size_t i = 0; i < outgoing_edges.size(); ++i) {
        //test forward traversal given an output edge label
        EXPECT_EQ(graph->outgoing(i + 1, graph->get_W(i + 1)), outgoing_edges[i]);
        //test that there is only one terminus
        auto sink_rank = graph->get_equal_node_range(graph->p);
        if (i + 1 < sink_rank.first || i + 1 > sink_rank.second) {
            EXPECT_EQ(graph->outgoing(i + 1, 0), 0u);
        } else {
            EXPECT_EQ(graph->outgoing(i + 1, 0), 1u);
        }
        //test FM index property
        if (graph->get_W(i + 1) < graph->alph_size) {
            EXPECT_TRUE(graph->get_last(graph->fwd(i + 1)));
            EXPECT_EQ(graph->get_W(i + 1), graph->get_node_end_value(graph->fwd(i + 1)));
        }
    }

    delete graph;
    kseq_destroy(read_stream);
    gzclose(input_p);
}

