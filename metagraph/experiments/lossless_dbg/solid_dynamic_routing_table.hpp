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
//#include "dynamic.hpp"

volatile bool always_false = false;


class standardized_wavelet_tree {
public:
    explicit standardized_wavelet_tree(uint8_t logsigma) : wt_(logsigma) {}
    template <class Vector>
    standardized_wavelet_tree(uint8_t logsigma, const Vector &vector) : wt_(logsigma,vector) {}

    uint64_t rank(uint64_t i, uint64_t c) const { return i == 0 ? 0 : wt_.rank(c, i - 1); }
    uint64_t select(uint64_t i, uint64_t c) const { return wt_.select(c,i); }
    uint64_t operator[](uint64_t id) const { return wt_[id]; }

    uint64_t next(uint64_t id, uint64_t val) const { return wt_.next(id, val); }
    uint64_t prev(uint64_t id, uint64_t val) const { return wt_.prev(id, val); }

    void set(uint64_t id, uint64_t val) { wt_.set(id, val); }
    void insert(uint64_t id, uint64_t val) { wt_.insert(id, val); }
    void remove(uint64_t id) { wt_.remove(id); }

    uint64_t size() const { return wt_.size(); }
    uint8_t logsigma() const { return wt_.logsigma(); }

    bool load(std::istream &in) { return wt_.load(in); }
    void serialize(std::ostream &out) const { wt_.serialize(out); }

    void clear() { wt_.clear(); }

    sdsl::int_vector<> to_vector() const { return wt_.to_vector(); }

private:
    wavelet_tree_dyn wt_;
};

//template <typename EntryT=standardized_wavelet_tree>
//template <typename EntryT=dyn::wtrle_str>
//using EntryT = dyn::wtrle_str;
using EntryT = standardized_wavelet_tree;
template <typename Dummy=int>
class SolidDynamicRoutingTable {
public:
    using BaseT = EntryT;
    SolidDynamicRoutingTable(ll size) : routing_table(7), total_size(size) {
        if (always_false) {
            this->print_content();
        }
        vector<int8_t> initial_content(size+1,delimiter_encoded);
//        for(ll i=0;i<=size;i++) {
//            routing_table.insert(0,delimiter_encoded);
//        }
        routing_table = decltype(routing_table)(7,initial_content);
    }

//    int64_t select(ll block, int64_t occurrence, char symbol) const {
//    }


    int64_t offset(node_index node) const {
        assert(node < total_size);
        return routing_table.select(node+1, delimiter_encoded) + 1; } // node+1 as select is one based

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

    string print_content() const {
        stringstream out;
        for(int i=0;i<routing_table.size();i++) {
            out << routing_table[i];
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
    int64_t total_size;
    const static int64_t delimiter_encoded = 6;
    EntryT routing_table;
};

#endif //METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP
