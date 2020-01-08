//
// Created by Jan Studený on 2019-06-04.
//

#ifndef __REFERENCE_DYNAMIC_ROUTING_TABLE_HPP__
#define __REFERENCE_DYNAMIC_ROUTING_TABLE_HPP__

#include "utilities.hpp"
#include <variant>
#include <filesystem>

// Verified implementation of routing table (to test against)
template<typename DummyT=int>
class ReferenceDynamicRoutingTable {
public:
    ReferenceDynamicRoutingTable() = default;
    ReferenceDynamicRoutingTable(shared_ptr<const DBGSuccinct> ) {}

    template<typename MockA,typename MockB> // TODO rename MockA to decl...bit_vector
    ReferenceDynamicRoutingTable(shared_ptr<const DBGSuccinct> /* graph */, MockA* is_element,MockB* rank_element, ll /* chunks = DefaultChunks */) : routing_table(is_element,rank_element) {

    }

//    int select(node_index node, int occurrence, char symbol) const {
//    }


protected:
    // rank [0..position)
    int rank(node_index node, int position, char symbol) const {
        int result = 0;
        const auto &node_entry = routing_table.at(node);
        assert(position <= node_entry.size());
        for (size_t i = 0; i < position; ++i) {
            result += node_entry[i] == symbol;
        }
        return result;
    }
public:
    int select(node_index node,int rank,char symbol) const {
        int crank = 1;
        int i = 0;
        for(auto& c : routing_table.at(node)) {
            if (c == symbol) {
                if (crank == rank) {
                    return i;
                }
                else {
                    crank++;
                }
            }
            i++;
        }
        return INT_MAX;
    }
    char get(node_index node, int position) const {
        return routing_table.at(node).at(position);
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
        assert(position >= 0 && "traversing o");
        return get(node,position);
    }

    int new_relative_position(node_index node, int position) const {
        auto base = get(node,position);
        auto base_rank = rank(node,position,base);
        return base_rank;
    }

    int new_relative_position(node_index node, int position, ll /* hint_block_offset */) const {
        return new_relative_position(node,position);
    }

    int size(node_index node) const {
        if (routing_table.count(node))
            return routing_table.at(node).size();
        else
            return 0;
    }

    ll insert(node_index node, int position, char symbol) {
        assert(position <= size(node));
        routing_table[node].insert(routing_table[node].begin() + position, symbol);
        return -1;
    }
    DenseHashMap<vector<char>> routing_table;
};

#endif // __REFERENCE_DYNAMIC_ROUTING_TABLE_HPP__
