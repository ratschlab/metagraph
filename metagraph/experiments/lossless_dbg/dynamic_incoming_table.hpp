//
// Created by Jan Studený on 2019-05-06.
//

#ifndef METAGRAPH_DYNAMIC_INCOMING_TABLE_HPP
#define METAGRAPH_DYNAMIC_INCOMING_TABLE_HPP

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <tsl/hopscotch_map.h>

#include "graph_patch.hpp"
#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "cxx-prettyprint.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "routing_table.hpp"

using default_bit_vector = bit_vector_small;

template <typename GraphT=DBGSuccinct,typename _edge_identifier_t=char>
class DynamicIncomingTable {
public:
    using edge_identifier_t = _edge_identifier_t;
    explicit DynamicIncomingTable(const GraphT &graph) : graph(graph) {}

    int branch_offset(node_index node, edge_identifier_t incoming) const {
        int result = 0;
        for (char base : "$ACGTN") {
            if (base < incoming) {
                result += branch_size(node, base);
            }
        }
        return result;
    }

    bool is_join(node_index node) const {
        return size(node); // as 1
    }

    int branch_size(node_index node, edge_identifier_t incoming) const {
        auto it = incoming_table.find(node);
        if (it == incoming_table.end() || !it->second.count(incoming))
            return 0;

        return it->second.at(incoming);
    }

    int size(node_index node) const {
        int result = 0;
        for(char base : "$ACGTN") {
            result += branch_size(node,base);
        }
        return result;
    }

    int branch_offset_and_increment(node_index node,
                                    edge_identifier_t incoming) {
        int result = 0;

        auto it = incoming_table.find(node);

        if (it != incoming_table.end()) {
            const auto &offsets = it->second;

            for (char base : "$ACGTN") {
                if (base < incoming) {
                    if (offsets.count(base))
                        result += offsets.at(base);
                }
            }
            it.value()[incoming]++;
        } else {
            incoming_table[node][incoming]++;
        }

        return result;
    }

    bool has_new_reads(node_index node) const {
        return branch_size(node, '$');
    }

    // TODO: replace map with array
    tsl::hopscotch_map<node_index, map<char, int>> incoming_table;
    const GraphT &graph;
};

#endif //METAGRAPH_DYNAMIC_INCOMING_TABLE_HPP
