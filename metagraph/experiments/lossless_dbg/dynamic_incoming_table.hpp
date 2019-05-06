//
// Created by Jan Studený on 2019-05-06.
//

#ifndef METAGRAPH_INCOMING_TABLE_HPP
#define METAGRAPH_INCOMING_TABLE_HPP

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include "graph_patch.hpp"
#include "path_database.hpp"
#include "path_database_baseline.hpp"
#include "cxx-prettyprint.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "routing_table.hpp"

using default_bit_vector = bit_vector_small;

template <typename GraphT=DBGSuccinct,typename edge_identifier_t=char>
class DynamicIncomingTable {
public:
    explicit DynamicIncomingTable(const GraphT & graph) : graph(graph) {}
    int branch_offset(node_index node,edge_identifier_t incoming) const {
        int result = 0;
        for(auto& base : "$ACGTN") {
            if (base < incoming) {
                result += branch_size(node,base);
            }
        }
        return result;
    }

    bool is_join(node_index node) const {
        return size(node); // as 1
    }


    int branch_size(node_index node,edge_identifier_t incoming) const {
        return incoming_table.at(node).at(incoming);
    }

    int size(node_index node) const {
        int result = 0;
        for(auto& base : "$ACGTN") {
            result += branch_size(node,base);
        }
        return result;
    }

    void increment(node_index node,edge_identifier_t incoming)  {
        incoming_table[node][incoming]++;
    }

    bool has_new_reads(node_index node) const {
        return branch_size(node,'$');
    }

    std::map<node_index,map<char,int>> incoming_table;
    const GraphT & graph;
};



#endif //METAGRAPH_INCOMING_TABLE_HPP
