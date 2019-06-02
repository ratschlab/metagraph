//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_GRAPH_PATCH_HPP
#define METAGRAPH_GRAPH_PATCH_HPP

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>

#include "dbg_succinct.hpp"

class BetterDBGSuccinct : public DBGSuccinct {
public:
    using DBGSuccinct::DBGSuccinct;

    void call_incoming_kmers_mine(node_index node,const std::function<void(node_index, char)>& callback) const {
        for(auto c : "ACGTN"s) {
            auto new_node = traverse_back(node,c);
            if (new_node) {
                callback(new_node,c);
            }
        }
    }

    char get_outgoing_base(node_index node) const {
        char base;
        assert(outdegree(node) == 1);
        call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
        return base;
    }

    int64_t branch_id(node_index node,node_index prev_node) const {
        int64_t result = INT_MIN;
        int64_t i = 0;
        call_incoming_kmers_mine(node,[&i,&result,&prev_node](node_index node,char base) {
            if (node==prev_node) {
                result = i;
            }
            i++;
        });
        return result;
    }

    bool is_join(node_index node) const {
        return indegree(node) > 1;
    }
    bool is_split(node_index node) const {
        return outdegree(node) > 1;
    }

    int encode(char c) const {
        if (c == '#') return get_boss().alph_size;
        if (c == '$') return 0;
        return get_boss().encode(c);
    }
    char decode(int c) const {
        if (c == get_boss().alph_size) return '#';
        if (c == 0) return '$';
        return get_boss().decode(c);
    }

};
#define DBGSuccinct BetterDBGSuccinct

#endif //METAGRAPH_GRAPH_PATCH_HPP
