#include <utility>

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
#include "graph_patch.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "graph_preprocessor.hpp"
#include "routing_table_transformation.hpp"

using node_index = SequenceGraph::node_index;

template<typename DummyT=int>
class DynamicRoutingTableCore {
public:
    DynamicRoutingTableCore() = default;
    DynamicRoutingTableCore(const DBGSuccinct& graph) {}

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

    int size(node_index node) const {
        if (routing_table.count(node))
            return routing_table.at(node).size();
        else
            return 0;
    }

    void insert(node_index node, int position, char symbol) {
        assert(position <= size(node));
        routing_table[node].insert(routing_table[node].begin() + position, symbol);
    }
    tsl::hopscotch_map<node_index, vector<char>> routing_table;

};

template <typename DummyT=int>
class DynamicRoutingTable : public TransformationsEnabler<DynamicRoutingTableCore<DummyT>> {
    using TransformationsEnabler<DynamicRoutingTableCore<DummyT>>::TransformationsEnabler;
};

#endif //METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
