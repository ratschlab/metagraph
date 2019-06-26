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

#include "utilities.hpp"
#include "routing_table_transformation.hpp"


template<class Wavelet = sdsl::wt_rlmn<>>
class RoutingTableCore {
public:
    RoutingTableCore() = default;
    RoutingTableCore(const shared_ptr<const DBGSuccinct> graph) :
        graph(graph), delimiter_encoded(graph->get_boss().alph_size)
    {};

    template<class Container>
    void initialize_content(const Container& routing_table_array) {
        sdsl::int_vector<3> routing_table_array_encoded(routing_table_array.size());
        for(int64_t i=0;i<routing_table_array.size();i++) {
            routing_table_array_encoded[i] = routing_table_array[i] == '#' ? 
                                             delimiter_encoded : 
                                             encode(routing_table_array[i]);
        }
        construct_im(routing_table,routing_table_array_encoded,0);
    }

    template<class Container>
    explicit RoutingTableCore(const shared_ptr<const DBGSuccinct> graph,const Container& routing_table_array) :
        RoutingTableCore(graph)
     {
        initialize_content(routing_table_array);
    }

    int64_t offset(node_index node) const { return routing_table.select(node, delimiter_encoded) + 1; }

    int64_t select_unchecked(node_index node, int64_t occurrence, int64_t encoded_symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,encoded_symbol);
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,encoded_symbol) - routing_table_block;
    }

    int64_t select(node_index node, int64_t occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,encode(symbol));
        if (occurrence > rank(node,size(node)+1,symbol)) {
            PRINT_VAR(node,occurrence,symbol);
            print_content(node);
        }
        assert(occurrence <= rank(node,size(node)+1,symbol));
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,encode(symbol)) - routing_table_block;
    }

    int64_t rank(node_index node, int64_t position, char symbol) const {
        auto routing_table_block = offset(node);
        auto absolute_position = routing_table_block+position;
        auto occurrences_of_base_before_block = routing_table.rank(routing_table_block,encode(symbol));
        return routing_table.rank(absolute_position,encode(symbol)) - occurrences_of_base_before_block;
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

    Wavelet routing_table;
    int64_t serialize(std::ostream& out)const {
        return routing_table.serialize(out);
    }

    void load(std::istream& in) {
        return routing_table.load(in);
    }

    json get_statistics(int64_t verbosity=0) const {
        return {{"dummy", "true"}}; // throws exception if updating with empty value
    }

    const int64_t delimiter_encoded = 6;
    const shared_ptr<const DBGSuccinct> graph;
};

template <typename Wavelet=sdsl::wt_rlmn<>>
using RoutingTable = TransformationsEnabler<RoutingTableCore<Wavelet>>;

#endif //METAGRAPH_ROUTING_TABLE_HPP
