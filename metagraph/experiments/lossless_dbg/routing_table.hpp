//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_ROUTING_TABLE_HPP
#define METAGRAPH_ROUTING_TABLE_HPP

#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"

const char RoutingTableAlphabet[] = {'$','A','C','G','T','N','#','?'};
// improvement (constexpr use https://github.com/serge-sans-paille/frozen)
// TODO: get rid of std::map
const map<char, int> RoutingTableInverseAlphabet = {{'$',0},{'A',1},{'C',2},{'G',3},{'T',4},{'N',5},{'#',6},{'?',7}};
const auto& rte2int = RoutingTableInverseAlphabet;

using routing_character_t = int;

int rc(char c) {
    auto it = RoutingTableInverseAlphabet.find(c);
    if (it != RoutingTableInverseAlphabet.end())
        return it->second;

    return RoutingTableInverseAlphabet.at('N');
}

int operator""_rc(char c) {
    return rc(c);
}

char tochar(routing_character_t rc) {
    return RoutingTableAlphabet[rc];
}

template<class Wavelet = sdsl::wt_rlmn<>>
class RoutingTable {
public:
    RoutingTable() = default;
    template<class Container>
    explicit RoutingTable(const Container& routing_table_array) {
        sdsl::int_vector<0> routing_table_array_encoded(routing_table_array.size());
        for(int i=0;i<routing_table_array.size();i++) {
            routing_table_array_encoded[i] = RoutingTableInverseAlphabet.at(routing_table_array[i]);
        }
        construct_im(routing_table,routing_table_array_encoded,0);
    }

    int offset(node_index node) const { return routing_table.select(node, '#'_rc) + 1; }

    int select(node_index node, int occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,rc(symbol));
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,rc(symbol)) - routing_table_block;
    }

    int rank(node_index node, int position, char symbol) const {
        auto routing_table_block = offset(node);
        auto absolute_position = routing_table_block+position;
        auto occurrences_of_base_before_block = routing_table.rank(routing_table_block,rc(symbol));
        return routing_table.rank(absolute_position,rc(symbol)) - occurrences_of_base_before_block;
    }

    char get(node_index node, int position) const {
        auto routing_table_block = offset(node);
        return tochar(routing_table[routing_table_block+position]);
    }

    int size(node_index node) const {
        return select(node, 1, '#');
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

    Wavelet routing_table;
    int serialize(std::ostream& out)const {
        return routing_table.serialize(out);
    }
    int load(std::istream& in) {
        return routing_table.load(in);
    }
};

#endif //METAGRAPH_ROUTING_TABLE_HPP
