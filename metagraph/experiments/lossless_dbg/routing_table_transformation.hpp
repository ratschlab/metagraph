//
// Created by Jan Studený on 2019-05-21.
//

#ifndef __ROUTING_TABLE_TRANSFORMATION_HPP__
#define __ROUTING_TABLE_TRANSFORMATION_HPP__

#include <utility>


#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>
#include <nlohmann/json.hpp>

#include "utils.hpp"
#include "alphabets.hpp"

#include "utilities.hpp"
#include "graph_preprocessor.hpp"
#include "graph_patch.hpp"

bool transform_done = 0;

using node_index = SequenceGraph::node_index;
using json = nlohmann::json;

template<typename RoutingTable>
class TransformationsEnabler : public RoutingTable {
public:

//    explicit TransformationsEnabler(const shared_ptr<const DBGSuccinct> graph) : RoutingTable(graph) {
//#ifndef DISABLE_TRANSFORMATIONS
//        transformations = GraphPreprocessor(*graph).find_weak_splits();
//#endif
//    }

    template<typename ...Args>
    TransformationsEnabler(shared_ptr<const DBGSuccinct> graph, Args... args) : RoutingTable(graph,args...){
        #ifndef DISABLE_TRANSFORMATIONS
        transformations = GraphPreprocessor(*graph).find_weak_splits();
        #endif
    }

    //using RoutingTable::RoutingTable;

    //    template <typename Container>
//    TransformationsEnabler(const shared_ptr<const DBGSuccinct> graph,const Container& routing_table_array) :
//                                            RoutingTable(graph,routing_table_array) {
//        transformations = GraphPreprocessor(graph).find_weak_splits();
//                                            }
    explicit TransformationsEnabler(transformations_t transformations) : transformations(std::move(transformations)) {

    }

    char traversed_base(node_index node, int64_t position) const {
        char base = this->get(node,position);
        return transform(node,base);
    }

//TODO: Merge with hinted version
    int64_t new_relative_position(node_index node, int64_t position) const {
        auto base = this->get(node,position);
        auto base_rank = this->rank(node,position,base);
        if (transformations.count(node)) {
            auto transformation = transformations.at(node);
            if (is_affected(transformation,base)) {
                char opposite = opposite_element(transformation,base);
#ifdef ASSUME_NO_TRANSFORMATIONS
                assert(!this->rank(node,position,opposite));
#endif
                base_rank += this->rank(node,position,opposite);
            }
        }
        return base_rank;
    }

    int64_t new_relative_position(node_index node, int64_t position, ll hint_block_offset) const {
        auto base = this->get(node,position,hint_block_offset);
        auto base_rank = this->rank(node,position,base,hint_block_offset);
        if (transformations.count(node)) {
            auto transformation = transformations.at(node);
            if (is_affected(transformation,base)) {
                char opposite = opposite_element(transformation,base);
#ifdef ASSUME_NO_TRANSFORMATIONS
                assert(!this->rank(node,position,opposite));
#endif
                base_rank += this->rank(node,position,opposite,hint_block_offset);
            }
        }
        return base_rank;
    }

    char transform(node_index node, char base) const {
        if (transformations.count(node)) {
            auto transformation = transformations.at(node);
            return transform(transformation,base);
        }
        return base;
    }

    static bool is_affected(pair<char,char> transformation, char base) {
        return transformation.first == base or transformation.second == base;
    }
    static char opposite_element(pair<char,char> transformation, char base) {
        return transformation.first == base ? transformation.second : transformation.first;
    }
    static char transform(pair<char,char> transformation, char base) {
        return transformation.first == base ? transformation.second : base;
    }
    json get_statistics(int64_t verbosity=0) const {
        vector<string> transformed_nodes;
        for (auto&[node,transformation] : transformations) {
            if (this->graph) {
                transformed_nodes.push_back(this->graph->get_node_sequence(node));
            }
        }
        return {
            {"transformations", transformations.size()},
            {"transformed_nodes", transformed_nodes}
        };
    }

     transformations_t transformations;
};


#endif // __ROUTING_TABLE_TRANSFORMATION_HPP__
