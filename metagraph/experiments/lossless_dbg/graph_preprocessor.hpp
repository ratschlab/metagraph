//
// Created by Jan Studený on 2019-05-13.
//

#ifndef __GRAPH_PREPROCESSOR_HPP__
#define __GRAPH_PREPROCESSOR_HPP__

#include <tsl/hopscotch_map.h>
#include "sequence_graph.hpp"

using node_index = SequenceGraph::node_index;

using transformations_t = tsl::hopscotch_map<node_index,pair<char,char>>;

template<typename Graph>
class GraphPreprocessor {
public:
    const int64_t delay_gap = 0;
    GraphPreprocessor(const Graph &graph) : graph(graph) {}
    transformations_t find_weak_splits() {
        Timer timer;
        cerr << "Started finding weak splits" << endl;
        // find nodes that will
        vector<optional<pair<char,char>>> result_intermediate(graph.num_nodes()+1);
        #pragma omp parallel for
        for(uint64_t node=1;node<=graph.num_nodes();node++) {
            if (is_split_node(graph,node)) {
                auto possible_weak_split = is_weak_split(node);
                if (possible_weak_split) {
                    result_intermediate[node] = *possible_weak_split;
                }
            }
        }
        transformations_t result;
        int64_t node = 0;
        for(auto& pos_transformation : result_intermediate) {
            if (pos_transformation) {
                result[node] = *pos_transformation;
            }
            node++;
        }
        cerr << "Finished finding weak splits in " << timer.elapsed() << " sec." << endl;

        return result;
    }

    optional<pair<char,char>> is_weak_split(node_index node) {
        // == is substitution split for now
        tsl::hopscotch_map<node_index,vector<char>> joins;
        graph.call_outgoing_kmers(node,[&](node_index next_node, char base) {
            auto possible_join = first_join(next_node,graph.get_k());
            if (possible_join) {
                joins[possible_join].push_back(base);
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

    node_index first_join(node_index node, int64_t delay_steps) {
        if (delay_steps < -delay_gap) {
            return 0;
        }
        if (delay_steps <= 0 and is_join_node(graph,node)) {
            return node;
        }
        if (is_split_node(graph,node) or graph.outdegree(node) == 0) {
            return 0;
        }
        return first_join(graph.traverse(node,get_outgoing_base(graph,node)),delay_steps-1);
    }
    const Graph& graph;

};

#endif // __GRAPH_PREPROCESSOR_HPP__
