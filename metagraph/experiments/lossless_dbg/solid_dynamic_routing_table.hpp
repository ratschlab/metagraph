//
// Created by Jan Studený on 2019-06-02.
//

#ifndef METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP
#define METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP

#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>

#include "graph_patch.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"

template <typename EntryT=wavelet_tree_dyn>
class SolidDynamicRoutingTable {
public:
    using BaseT = EntryT;
    SolidDynamicRoutingTable(ll size) : routing_table(7) {
        vector<int8_t> initial_content(size+1,encode('#'));
        for(ll i=0;i<=size;i++) {
            routing_table.insert(0,encode('#'));
        }

        //routing_table = decltype(routing_table)(7,initial_content);
    }

//    int64_t select(ll block, int64_t occurrence, char symbol) const {
//    }


    int64_t offset(node_index node) const { return routing_table.select( delimiter_encoded,node+1) + 1; }

    int64_t select_unchecked(node_index node, int64_t occurrence, int64_t encoded_symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(encoded_symbol,routing_table_block-inclusive);
        return routing_table.select(encoded_symbol,occurrences_of_symbol_before_block+occurrence) - routing_table_block;
    }

    int64_t select(node_index node, int64_t occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(encode(symbol),routing_table_block-inclusive);
        if (occurrence > rank(node,size(node)+1,symbol)) {
            PRINT_VAR(node,occurrence,symbol);
            print_content(node);
        }
        assert(occurrence <= rank(node,size(node)+1,symbol));
        return routing_table.select(encode(symbol),occurrences_of_symbol_before_block+occurrence) - routing_table_block;
    }

    int64_t rank(node_index node, int64_t position, char symbol) const {
        auto routing_table_block = offset(node);
        auto absolute_position = routing_table_block+position;
        auto occurrences_of_base_before_block = routing_table.rank(encode(symbol),routing_table_block-inclusive);
        return routing_table.rank(encode(symbol),absolute_position-1) - occurrences_of_base_before_block;
    }

    char get(node_index node, int64_t position) const {
        auto routing_table_block = offset(node);
        return decode(routing_table[routing_table_block+position]);
    }

    int64_t size(node_index node) const {
        return select_unchecked(node, 1, delimiter_encoded);
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for (int64_t i=0;i<table_size;i++) {
            out << get(node, i);
        }
        out << endl;
        cerr << out.str();
        return out.str();
    }

    char traversed_base(node_index node, int64_t position) const {
        return get(node,position);
    }

    int64_t new_relative_position(node_index node, int64_t position) const {
        auto base = get(node,position);
        auto base_rank = rank(node,position,base);
        return base_rank;
    }

    void insert(ll block, int64_t position, char symbol) {
        int64_t encoded = encode(symbol);
        assert(position <= size(block));
        routing_table.insert(offset(block)+position,encoded);
    }

    // Todo change from magic number to some define elsewhere
    const static int64_t delimiter_encoded = 6;
    const static bool inclusive = true;
    wavelet_tree_dyn routing_table;
};

#endif //METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP
