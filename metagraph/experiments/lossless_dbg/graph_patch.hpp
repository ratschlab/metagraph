//
// Created by Jan Studený on 2019-04-26.
//

#ifndef __GRAPH_PATCH_HPP__
#define __GRAPH_PATCH_HPP__

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
        int count = 0;
        for(auto c : "ACGTN"s) {
            auto new_node = traverse_back(node,c);
            if (new_node) {
                callback(new_node,c);
                count++;
            }
        }
#ifdef MASK_DUMMY_KMERS
        alt_assert((count == indegree(node) || [&]() {
            PRINT_VAR(node,count,this->indegree(node));
            std::vector<node_index> adjacent_nodes;
            adjacent_nodes.reserve(10);

            this->adjacent_incoming_nodes(node, &adjacent_nodes);
            for(auto& e : adjacent_nodes) {
                cout << this->get_node_sequence(e) << endl;
            }
            return false; }()));
#endif
        if (count != indegree(node)) {
            //PRINT_VAR("has dummy input kmer",node,this->get_node_sequence(node),traverse_back(node,'$'));
            callback(0,'$'); // TODO: 0 is false node index (find the right one)
        }
    }

    char get_outgoing_base(node_index node) const {
        char base;
        assert(outdegree(node) == 1);
        call_outgoing_kmers(node,[&base](node_index node,char edge_label ) {
            PRINT_VAR(node,edge_label);
            base = edge_label;});
        return base;
    }

    size_t incoming_edge_rank(DeBruijnGraph::node_index source,
                             DeBruijnGraph::node_index target) const {
        assert(source && source <= num_nodes());
        assert(target && target <= num_nodes());

        assert(get_node_sequence(source).substr(1)
               == get_node_sequence(target).substr(0, get_k() - 1));

        std::vector<node_index> adjacent_nodes;
        adjacent_nodes.reserve(10);

        adjacent_incoming_nodes(target, &adjacent_nodes);

        uint64_t edge_rank = 0;

        for (node_index node : adjacent_nodes) {
            if (node == source)
                return edge_rank;

            edge_rank++;
        }

        throw std::runtime_error("the edge does not exist in graph");
    }

    int64_t branch_id(node_index node,node_index prev_node) const {
        return incoming_edge_rank(prev_node, node);
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

    int8_t encode(char c) const {
        if (c == '#') return get_boss().alph_size;
        if (c == '$') return 0;
        return get_boss().encode(c);
    }
    char decode(int8_t c) const {
        if (c == get_boss().alph_size) return '#';
        if (c == 0) return '$';
        return get_boss().decode(c);
    }

};
#define DBGSuccinct BetterDBGSuccinct

#endif // __GRAPH_PATCH_HPP__
