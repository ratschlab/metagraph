#include <stdio.h>
#include <string>
#include <sstream>

#include "gtest/gtest.h"

#include "kseq.h"
#include "construct.hpp"
#include "dbg_succinct_libmaus.hpp"
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

TEST(Construct, EmptyGraph) {
    DBG_succ *graph = init_graph();
    construct_succ(graph);
    graph->switch_state(Config::DYN);
    std::ostringstream ostr;
    ostr << *(graph->last);
    EXPECT_EQ(ostr.str(), "01111");
    EXPECT_EQ(graph->p, 4u);
    EXPECT_EQ(graph->W->size(), 5u);
    EXPECT_EQ(graph->last->size(), 5u);
    delete graph;
}

TEST(Construct, ShortGraphWithLoop) {
    gzFile input_p = gzopen(test_fasta.c_str(), "r");
    kseq_t *read_stream = kseq_init(input_p);
    ASSERT_TRUE(read_stream);
    DBG_succ *graph = init_graph();
    construct::add_sink(graph);
    for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
        construct::add_seq_fast(graph, read_stream->seq, std::string(read_stream->name.s, read_stream->name.l), true);
    }
    construct::construct_succ(graph);
    graph->switch_state(Config::DYN);
    std::ostringstream ostr;
    ostr << *(graph->last);
    EXPECT_EQ(ostr.str(), "01011110111111111111111111111001");
    ostr.clear();
    ostr.str("");
    ostr << *(graph->W);
    EXPECT_EQ(ostr.str(), "06241463411641812643366666131313013");
    delete graph;
    kseq_destroy(read_stream);
    gzclose(input_p);
}
