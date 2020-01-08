//
// Created by Jan Studený on 2019-04-26.
//

#ifndef __GRAPH_PATCH_HPP__
#define __GRAPH_PATCH_HPP__

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <prettyprint.hpp>

#include "dbg_succinct.hpp"
using node_index = SequenceGraph::node_index;

template <typename Graph=DBGSuccinct>
char get_outgoing_base(const Graph& g,node_index node) {
    char base;
    assert(g.outdegree(node) == 1);
    g.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) {
        PRINT_VAR(node,edge_label);
        base = edge_label;});
    return base;
}

template <typename Graph=DBGSuccinct>
bool is_join_node(const Graph& g, node_index node)  {
    return g.indegree(node) > 1;
}
template <typename Graph=DBGSuccinct>
bool is_split_node(const Graph& g,node_index node)  {
    return g.outdegree(node) > 1;
}

#endif // __GRAPH_PATCH_HPP__
