//
// Created by Jan Studený on 2019-06-02.
//

#ifndef METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP
#define METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP

#include <iostream>
#include <set>
#include <map>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>

#include "graph_patch.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "wavelet_tree.hpp"
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
        vector<uint8_t> initial_content(size+1, delimiter_encoded);
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
        auto result = routing_table.select(occurrences_of_symbol_before_block+occurrence,encoded_symbol) - routing_table_block;
        assert(result>=0);
        return result;
    }

    int64_t select(node_index node, int64_t occurrence, char symbol) const {
        auto routing_table_block = offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,encode(symbol));
        assert((occurrence <= rank(node,size(node)+1,symbol) || [&](){
            PRINT_VAR(node,occurrence,symbol);
            print_content(node);
            return false; }()) );
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,encode(symbol)) - routing_table_block;
    }

    int64_t rank(node_index node, int64_t position, char symbol, ll hint_block_offset) const {
        assert(hint_block_offset == offset(node));
        auto absolute_position = hint_block_offset+position;
        auto occurrences_of_base_before_block = routing_table.rank(hint_block_offset,encode(symbol));
        return routing_table.rank(absolute_position,encode(symbol)) - occurrences_of_base_before_block;
    }

    int64_t rank(node_index node, int64_t position, char symbol) const {
        return rank(node,position,symbol,offset(node));
    }

    char get(node_index node, int64_t position, ll hint_block_offset) const {
        assert(offset(node) == hint_block_offset);
        return decode(routing_table[hint_block_offset+position]);
    }

    char get(node_index node, int64_t position) const {
        return get(node,position,offset(node));
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

    int64_t new_relative_position(node_index node, int64_t position, int hint_block_offset) const {
        auto base = get(node,position,hint_block_offset);
        auto base_rank = rank(node,position,base,hint_block_offset);
        return base_rank;
    }

    int64_t new_relative_position(node_index node, int64_t position) const {
        return new_relative_position(node,position,offset(node));
    }

    ll insert(ll block, int64_t position, char symbol) {
        int64_t encoded = encode(symbol);
        assert((position <= size(block) || [&](){ PRINT_VAR(block,position,symbol,size(block),encoded,symbol); return false; }()));
        auto block_offset = offset(block);
        routing_table.insert(block_offset+position,encoded);
        return block_offset;
    }

    // Todo change from magic number to some define elsewhere
    int64_t total_size;
    const static int64_t delimiter_encoded = 6;
    EntryT routing_table;
};

#endif //METAGRAPH_SOLID_DYNAMIC_ROUTING_TABLE_HPP
