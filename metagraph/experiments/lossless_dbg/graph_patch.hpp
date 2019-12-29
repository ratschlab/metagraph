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
#include "cxx-prettyprint.hpp"

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


//void _call_incoming_kmers_mine(node_index node,const std::function<void(node_index, char)>& callback) const {
//    // doesn't work with call_incoming_kmers because of unknown reasons
//    // speculation: maybe the order of bases changes
//    // addition of $ as a first letter didn't have any effect
//    // TODO: try: vector & sort call_incoming_kmers
//    int count = 0;
//    for(auto c : "ACGTN"s) {
//        auto new_node = traverse_back(node,c);
//        if (new_node) {
//            callback(new_node,c);
//            count++;
//        }
//    }
//#ifdef MASK_DUMMY_KMERS
//    alt_assert((count == indegree(node) || [&]() {
//            PRINT_VAR(node,count,this->indegree(node));
//            std::vector<node_index> adjacent_nodes;
//            adjacent_nodes.reserve(10);
//
//            this->adjacent_incoming_nodes(node, &adjacent_nodes);
//            for(auto& e : adjacent_nodes) {
//                cout << this->get_node_sequence(e) << endl;
//            }
//            return false; }()));
//#endif
//
//    if (count != indegree(node)) {
//        PRINT_VAR("has dummy input kmer",node,this->get_node_sequence(node),traverse_back(node,'$'));
//        callback(0,'$'); // TODO: 0 is false node index (find the right one)
//    }
//}




class BetterDBGSuccinct : public DBGSuccinct {
    // make all member methods functions
public:
    using DBGSuccinct::DBGSuccinct;







    bool is_join(node_index node) const {
        return indegree(node) > 1;
    }
    bool is_split(node_index node) const {
        return outdegree(node) > 1;
    }

    int8_t encode(char c) const {
        if (c == '#') return get_boss().alph_size;
        if (c == '$') return 0;
        return get_boss().encode(c);
    }
    char decode(int8_t c) const {
        if (c == get_boss().alph_size) return '#';
        if (c == 0) return '$';
        return get_boss().decode(c);
    }

};

#define DBGSuccinct BetterDBGSuccinct

#endif // __GRAPH_PATCH_HPP__
