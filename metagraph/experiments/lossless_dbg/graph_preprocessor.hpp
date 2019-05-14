//
// Created by Jan Studený on 2019-05-13.
//

#ifndef METAGRAPH_GRAPH_PREPROCESSOR_HPP
#define METAGRAPH_GRAPH_PREPROCESSOR_HPP

#include <tsl/hopscotch_map.h>
#include "sequence_graph.hpp"


using transformations_t = tsl::hopscotch_map<node_index,pair<char,char>>;

template<typename Graph>
class GraphPreprocessor {
public:
    const int delay_gap = 0;
    GraphPreprocessor(const Graph &graph) : graph(graph) {}
    transformations_t find_weak_splits() {
        // find nodes that will
        transformations_t result;
        for(auto node=1;node<=graph.num_nodes();node++) {
            if (graph.is_split(node)) {
                auto possible_weak_split = is_weak_split(node);
                if (possible_weak_split) {
                    result[node] = *possible_weak_split;
                }
            }
        }
        return result;
    }

    optional<pair<char,char>> is_weak_split(node_index node) {
        // == is substitution split for now
        tsl::hopscotch_map<node_index,vector<char>> joins;
        graph.call_outgoing_kmers(node,[&](node_index next_node, char base) {
            auto possible_join = first_join(next_node,graph.get_k());
            if (possible_join) {
                joins[node].push_back(base);
            }
        });
        for(auto& [join_node,bases] : joins) {
            if (bases.size() > 1) {
                // currently just randomly pick one;
                return {{bases[0],bases[1]}};
            }
        }
        return {};
    };

    node_index first_join(node_index node, int delay_steps) {
        if (delay_steps < -delay_gap) {
            return 0;
        }
        if (delay_steps <= 0 and graph.is_join(node)) {
            return node;
        }
        if (graph.is_split(node) or graph.outdegree(node) == 0) {
            return 0;
        }
        return first_join(graph.traverse(node,graph.get_outgoing_base(node)),delay_steps-1);
    }
    const Graph& graph;

};

#endif //METAGRAPH_GRAPH_PREPROCESSOR_HPP
