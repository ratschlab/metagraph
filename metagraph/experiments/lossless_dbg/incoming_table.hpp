//
// Created by Jan Studený on 2019-05-06.
//

#ifndef METAGRAPH_INCOMING_TABLE_HPP
#define METAGRAPH_INCOMING_TABLE_HPP


using default_bit_vector = bit_vector_small;

template <typename BitVector=default_bit_vector,typename GraphT=DBGSuccinct>
class DynamicIncomingTable {
public:
    explicit DynamicIncomingTable(const GraphT & graph) : graph(graph) {}
    using bit_vector_t = BitVector;
    BitVector joins;
    sdsl::enc_vector<> edge_multiplicity_table;
    int branch_offset(node_index node,node_index prev_node) const {
        int branch_offset = relative_offset(node,prev_node);
        int result = 0;
        for(int i=0;i<branch_offset;i++) {
            result += branch_size_rank(node,i);
        }
        return result;
    }

    bool is_join(node_index node) const {
        return size(node); // as 1
    }

    int branch_size_rank(node_index node,int offset) const {
        if (offset < 0) { return 0; }
        int joins_position = joins.select1(node);
        int table_offset = joins.rank0(joins_position);
        return edge_multiplicity_table[table_offset+offset];
    }

    int branch_size(node_index node,node_index prev_node) const {
        return branch_size_rank(node,relative_offset(node,prev_node));
    }

    int size(node_index node) const {
        return joins.select1(node+1)-joins.select1(node)-1;
    }

    bool has_new_reads(node_index node) const {
        // todo : remove or after the indegree bug is fixed and use assertion
        // assert(!(graph.indegree(node) < 2 and size(node)) or size(node) > graph.indegree(node));
        // indegree smaller than two with nonempty incoming table implies that node has new reads
        return (size(node) > graph.indegree(node)) or (graph.indegree(node) < 2 and size(node));
    }

    int relative_offset(node_index node,node_index prev_node) const {
        bool increment = has_new_reads(node);
        int result;
        if (prev_node) {
            result = graph.branch_id(node,prev_node);
        }
        else {
            result = -1; // starting branch is before all other branches
        }
        return result + increment;
    }

    const GraphT & graph;
};

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

#endif //METAGRAPH_INCOMING_TABLE_HPP
