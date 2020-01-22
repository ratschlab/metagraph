//
// Created by Jan Studený on 2019-04-26.
//

#ifndef __ROUTING_TABLE_HPP__
#define __ROUTING_TABLE_HPP__

#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <sdsl/construct.hpp>

//#include "utils.hpp"
#include "configuration.hpp"
#include "alphabets.hpp"

#include "utilities.hpp"
#include "routing_table_transformation.hpp"

template <class Wavelet = WaveletTreeRLMN>
class RoutingTableCore {
  public:
    RoutingTableCore() = default;
    RoutingTableCore(const shared_ptr<const DBGSuccinct> graph) : graph(graph) {};

    template <class Container>
    void initialize_content(const Container &routing_table_array) {
        // can't use int_vector<3> as the construction of routing table fails (routing table is malformed)
        // now even (routing_table_array.size(),0,3) is malformed
        // int8_t* routing_table_array_encoded = new int8_t[routing_table_array.size()]; // also fails
        // DEP_TODO: change to use only 3 bits (dep: https://github.com/simongog/sdsl-lite/issues/422
        // needs to be solved OR https://github.com/xxsds/sdsl-lite/issues/71 solved & change to xxsds/sdsl-lite)
        sdsl::int_vector<8> routing_table_array_encoded(
                routing_table_array.size()); //,0,8); // 3->8
        for (uint64_t i = 0; i < routing_table_array.size(); i++) {
            routing_table_array_encoded[i] = encode(routing_table_array[i]);
        }

        construct_im(routing_table, routing_table_array_encoded, 0);
        assert(routing_table.size() == routing_table_array_encoded.size());
        for (uint64_t i = 0; i < routing_table_array.size(); i++) {
            assert(routing_table[i] == routing_table_array_encoded[i]);
        }
    }

    template <class Container>
    explicit RoutingTableCore(const shared_ptr<const DBGSuccinct> graph,
                              const Container &routing_table_array)
        : RoutingTableCore(graph) {
        initialize_content(routing_table_array);
    }

    int64_t offset(node_index node) const {
        return routing_table.select(node, delimiter_encoded) + 1;
    }

    int64_t select_unchecked(node_index node, int64_t occurrence, int64_t encoded_symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block
                = routing_table.rank(routing_table_block, encoded_symbol);
        return routing_table.select(occurrences_of_symbol_before_block + occurrence,
                                    encoded_symbol)
                - routing_table_block;
    }

    int64_t select(node_index node, int64_t occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block
                = routing_table.rank(routing_table_block, encode(symbol));
        if (occurrence > rank(node, size(node) + 1, symbol)) {
            PRINT_VAR(graph->get_node_sequence(node))(std::string());
            PRINT_VAR(node, occurrence, symbol)(std::string());
            print_content(node);
        }
        assert(occurrence <= rank(node, size(node) + 1, symbol));
        return routing_table.select(occurrences_of_symbol_before_block + occurrence,
                                    encode(symbol))
                - routing_table_block;
    }

    int64_t rank(node_index node, int64_t position, char symbol) const {
        auto routing_table_block = offset(node);
        auto absolute_position = routing_table_block + position;
        auto occurrences_of_base_before_block
                = routing_table.rank(routing_table_block, encode(symbol));
        return routing_table.rank(absolute_position, encode(symbol))
                - occurrences_of_base_before_block;
    }

    char get(node_index node, int64_t position) const {
        auto routing_table_block = offset(node);
        return decode(routing_table[routing_table_block + position]);
    }

    int64_t size(node_index node) const {
        return select_unchecked(node, 1, delimiter_encoded);
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for (int64_t i = 0; i < table_size; i++) {
            out << get(node, i);
        }
        mg::common::logger->debug(out.str());
        return out.str();
    }

    string print_content() const {
        stringstream out;
        for (uint64_t i = 0; i < routing_table.size(); i++) {
            out << decode(routing_table[i]);
        }

        mg::common::logger->debug(out.str());
        return out.str();
    }


    char traversed_base(node_index node, int64_t position) const {
        return get(node, position);
    }

    int64_t new_relative_position(node_index node, int64_t position) const {
        auto base = get(node, position);
        auto base_rank = rank(node, position, base);
        return base_rank;
    }

    Wavelet routing_table;
    int64_t serialize(std::ostream &out) const { return routing_table.serialize(out); }

    void load(std::istream &in) { return routing_table.load(in); }

    json get_statistics(int64_t) const {
        return { { "dummy", "true" } }; // throws exception if updating with empty value
    }

    const int64_t delimiter_encoded = 6;
    const shared_ptr<const DBGSuccinct> graph;
};

template <typename Wavelet = WaveletTreeRLMN>
using RoutingTable = TransformationsEnabler<RoutingTableCore<Wavelet>>;

#endif // __ROUTING_TABLE_HPP__
