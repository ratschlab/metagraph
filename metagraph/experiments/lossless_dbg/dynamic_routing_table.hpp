//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
#define METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP

#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include "path_database.hpp"
#include "path_database_baseline.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "dbg_succinct.hpp"
using node_index = SequenceGraph::node_index;

class DynamicRoutingTable {
public:
    DynamicRoutingTable() = default;
    template<class Container>
    explicit DynamicRoutingTable(const Container& routing_table_array) {
    }

//    int select(node_index node, int occurrence, char symbol) const {
//    }

    int rank(node_index node, char symbol, int position) const {
        assert(position < routing_table.size());
        int result = 0;
        int i = 0;
        auto& node_entry = routing_table.at(node);
        for(auto it=begin(node_entry);i<position and it != end(node_entry);it++,i++) {
            result += *it == symbol;
        }
        return result;
    }

    char get(node_index node, int position) const {
        return routing_table.at(node).at(position);
    }

    void insert(node_index node, int position, char symbol) {
        auto rt_index = routing_table[node].begin();
        assert(position <= size(node));
        advance(rt_index,position);
        routing_table[node].insert(rt_index,symbol);
    }

    int size(node_index node) const {
        return routing_table.at(node).size();
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for (int i=0;i<table_size;i++) {
            out << get(node, i);
        }
        out << endl;
        cout << out.str();
        return out.str();
    }

    map<node_index,vector<char>> routing_table;

//    int serialize(std::ostream& out)const {
//    }
//    int load(std::istream& in) {
//    }
};





#endif //METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
