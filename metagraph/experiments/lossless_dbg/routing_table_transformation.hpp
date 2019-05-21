//
// Created by Jan Studený on 2019-05-21.
//

#ifndef METAGRAPH_ROUTING_TABLE_TRANSFORMATION_HPP
#define METAGRAPH_ROUTING_TABLE_TRANSFORMATION_HPP

#include <utility>


#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>

#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "dbg_succinct.hpp"
#include "graph_preprocessor.hpp"

using node_index = SequenceGraph::node_index;

template<typename RoutingTable>
class TransformationsEnabler : public RoutingTable {
public:
    using RoutingTable::RoutingTable;

    explicit TransformationsEnabler(transformations_t transformations) : transformations(std::move(transformations)) {

    }

    char traversed_base(node_index node, int position) const {
        char base = this->get(node,position);
        return transform(node,base);
    }

    int new_relative_position(node_index node, int position) const {
        auto base = this->get(node,position);
        auto base_rank = this->rank(node,position,base);
        if (transformations.count(node)) {
            auto transformation = transformations.at(node);
            if (is_affected(transformation,base)) {
                char opposite = opposite_element(transformation,base);
                base_rank += this->rank(node,position,opposite);
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

     transformations_t transformations;
};


#endif //METAGRAPH_ROUTING_TABLE_TRANSFORMATION_HPP
