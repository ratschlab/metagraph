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

#include "utils.hpp"
#include "alphabets.hpp"

#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "utilities.hpp"
#include "routing_table_transformation.hpp"

const char RoutingTableCoreAlphabet[] = {'$','A','C','G','T','N','#','?'};
// improvement (constexpr use https://github.com/serge-sans-paille/frozen)
// TODO: get rid of std::map
const map<char, int> RoutingTableCoreInverseAlphabet = {{'$',0},{'A',1},{'C',2},{'G',3},{'T',4},{'N',5},{'#',6},{'?',7}};
const auto& rte2int = RoutingTableCoreInverseAlphabet;

using routing_character_t = int;

int rc(char c) {
    auto it = RoutingTableCoreInverseAlphabet.find(c);
    if (it != RoutingTableCoreInverseAlphabet.end())
        return it->second;
    return RoutingTableCoreInverseAlphabet.at('N');
}

int operator""_rc(char c) {
    return rc(c);
}

char tochar(routing_character_t rc) {
    return RoutingTableCoreAlphabet[rc];
}

template<class Wavelet = sdsl::wt_rlmn<>>
class RoutingTableCore {
public:
    RoutingTableCore(const DBGSuccinct& graph) :
        graph(graph), delimiter_encoded(graph.get_boss().alph_size)
    {};

    template<class Container>
    void initialize_content(const Container& routing_table_array) {
        sdsl::int_vector<0> routing_table_array_encoded(routing_table_array.size());
        for(int i=0;i<routing_table_array.size();i++) {
            routing_table_array_encoded[i] = routing_table_array[i] == '#' ? 
                                             delimiter_encoded : 
                                             graph.encode(routing_table_array[i]);
        }
        construct_im(routing_table,routing_table_array_encoded,0);
    }

    template<class Container>
    explicit RoutingTableCore(const DBGSuccinct& graph,const Container& routing_table_array) :
        RoutingTableCore(graph)
     {
        initialize_content(routing_table_array);
    }

    int offset(node_index node) const { return routing_table.select(node, delimiter_encoded) + 1; }

    int select_unchecked(node_index node, int occurrence, int encoded_symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,encoded_symbol);
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,encoded_symbol) - routing_table_block;
    }

    int select(node_index node, int occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,graph.encode(symbol));
        if (occurrence > rank(node,size(node)+1,symbol)) {
            PRINT_VAR(node,occurrence,symbol);
            print_content(node);
        }
        assert(occurrence <= rank(node,size(node)+1,symbol));
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,graph.encode(symbol)) - routing_table_block;
    }

    int rank(node_index node, int position, char symbol) const {
        auto routing_table_block = offset(node);
        auto absolute_position = routing_table_block+position;
        auto occurrences_of_base_before_block = routing_table.rank(routing_table_block,graph.encode(symbol));
        return routing_table.rank(absolute_position,graph.encode(symbol)) - occurrences_of_base_before_block;
    }

    char get(node_index node, int position) const {
        auto routing_table_block = offset(node);
        return tochar(routing_table[routing_table_block+position]);
    }

    int size(node_index node) const {
        return select_unchecked(node, 1, delimiter_encoded);
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for (int i=0;i<table_size;i++) {
            out << get(node, i);
        }
        out << endl;
        cerr << out.str();
        return out.str();
    }

    char traversed_base(node_index node, int position) const {
        return get(node,position);
    }

    int new_relative_position(node_index node, int position) const {
        auto base = get(node,position);
        auto base_rank = rank(node,position,base);
        return base_rank;
    }

    Wavelet routing_table;
    int serialize(std::ostream& out)const {
        return routing_table.serialize(out);
    }

    void load(std::istream& in) {
        return routing_table.load(in);
    }
    const int delimiter_encoded;
    const DBGSuccinct& graph;
};

template <typename Wavelet>
class RoutingTable : public TransformationsEnabler<RoutingTableCore<Wavelet>> {
    using TransformationsEnabler<RoutingTableCore<Wavelet>>::TransformationsEnabler;
};

#endif //METAGRAPH_ROUTING_TABLE_HPP
