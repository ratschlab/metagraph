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
#include <tsl/hopscotch_map.h>

#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "dbg_succinct.hpp"

using node_index = SequenceGraph::node_index;

class DynamicRoutingTable {
public:
    DynamicRoutingTable() = default;
    template<class Container>
    explicit DynamicRoutingTable(const Container& routing_table_array) {}

//    int select(node_index node, int occurrence, char symbol) const {
//    }

    int rank(node_index node, char symbol, int position) const {
        assert(position < routing_table.size());
        int result = 0;
        const auto &node_entry = routing_table.at(node);
        assert(position <= node_entry.size());
        for (size_t i = 0; i < position; ++i) {
            result += node_entry[i] == symbol;
        }
        return result;
    }

    char get(node_index node, int position) const {
        return routing_table.at(node).at(position);
    }

    void insert(node_index node, int position, char symbol) {
        routing_table[node].insert(routing_table[node].begin() + position, symbol);
    }

    int size(node_index node) const {
        if (routing_table.count(node))
         return routing_table.at(node).size();
        else
            return 0;
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


    tsl::hopscotch_map<node_index, vector<char>> routing_table;

    int new_relative_position(node_index node, int position) {
        return rank(node,get(node,position),position);
    }


    //map<node_index,vector<char>> routing_table;

//    int serialize(std::ostream& out)const {
//    }
//    int load(std::istream& in) {
//    }
};

#endif //METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
