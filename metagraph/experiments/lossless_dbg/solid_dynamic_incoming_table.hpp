//
// Created by Jan Studený on 2019-06-02.
//

#ifndef METAGRAPH_SOLID_DYNAMIC_INCOMING_TABLE_HPP
#define METAGRAPH_SOLID_DYNAMIC_INCOMING_TABLE_HPP

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <tsl/hopscotch_map.h>

#include "graph_patch.hpp"
#include "cxx-prettyprint.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"

using default_bit_vector = bit_vector_small;

template <typename GraphT=DBGSuccinct,typename _edge_identifier_t=char>
class SolidDynamicIncomingTable {
public:
    using edge_identifier_t = _edge_identifier_t;
    explicit SolidDynamicIncomingTable(ll size) {}

    int64_t branch_offset(node_index node, edge_identifier_t incoming) const {
        int64_t result = 0;
        for (char base : "$ACGTN") {
            if (base < incoming) {
                result += branch_size(node, base);
            }
        }
        return result;
    }

    bool is_join(node_index node) const {
        return incoming_table.count(node);
    }

    int64_t branch_size(node_index node, edge_identifier_t incoming) const {
        int64_t result = 0;
        int64_t encoded = encode(incoming);
        auto it = incoming_table.find(node);
        if (it == incoming_table.end())
            return 0;

        return it->second.at(encoded);
    }

    ll size(node_index node) const {
        ll result = 0;
        for(char base : "$ACGTN") {
            result += branch_size(node,base);
        }
        return result;
    }

    int64_t branch_offset_and_increment(node_index node,
                                    edge_identifier_t incoming) {
        assert(incoming == '$' or
               incoming == 'A' or
               incoming == 'C' or
               incoming == 'G' or
               incoming == 'T' or
               incoming == 'N'
        );
        int result = 0;

        auto it = incoming_table.find(node);
        int64_t encoded = encode(incoming);

        if (it != incoming_table.end()) {
            const auto& array = it->second;

            for(auto i=0;i<6;i++) {
                if (i >= encoded) break;
                result += array[i];
            }
            it.value()[encoded]++;
        } else {
            incoming_table[node][encoded]++;
        }

        return result;
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for(char c : "$ACGTN") {
            out << c << ":" << branch_size(node,c) << endl;
        }
        cerr << out.str();
        return out.str();
    }

    bool has_new_reads(node_index node) const {
        return branch_size(node, '$');
    }

    tsl::hopscotch_map<node_index, array<int32_t,6>> incoming_table;
};

#endif //METAGRAPH_SOLID_DYNAMIC_INCOMING_TABLE_HPP
