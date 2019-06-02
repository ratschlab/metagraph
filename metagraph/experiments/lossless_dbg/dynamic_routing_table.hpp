#include <utility>

//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
#define METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP

#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>

#include "path_database.hpp"
#include "graph_patch.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "graph_preprocessor.hpp"
#include "routing_table_transformation.hpp"
#include "dense_hashmap.hpp"
#include "solid_dynamic_routing_table.hpp"

using node_index = SequenceGraph::node_index;



//template<typename EntryT=wavelet_tree_dyn>
template <typename EntryT=SolidDynamicRoutingTable<>,typename BitVector=sdsl::bit_vector,typename RankSupport=typename BitVector::rank_1_type>
class DynamicRoutingTableCore {
public:
    using BaseT = typename EntryT::BaseT;
    DynamicRoutingTableCore() = default;
    DynamicRoutingTableCore(shared_ptr<const DBGSuccinct> graph) : graph(graph) {}
    DynamicRoutingTableCore(shared_ptr<const DBGSuccinct> graph, BitVector* is_element,RankSupport* rank_element, ll chunks = 1000)
            : graph(graph), chunks(is_element,rank_element,chunks) {
    }

    ll select(int64_t location, ll occurrence, char symbol) const {
        auto &chunk = chunks.at(location);
        return chunk.select(chunks.position_in_chunk(location), occurrence, symbol);
    }

    ll rank(int64_t location, ll position, char symbol) const {
        auto &chunk = chunks.at(location);
        return chunk.rank(chunks.position_in_chunk(location), position, symbol);
    }

    char get(int64_t location, ll position) const {
        auto &chunk = chunks.at(location);
        return chunk.get(chunks.position_in_chunk(location), position);
    }

    ll size(int64_t location) const {
        auto &chunk = chunks.at(location);
        return chunk.size(chunks.position_in_chunk(location));
    }

    string print_content(int64_t location) {
        auto &chunk = chunks.at(location);
        return chunk.print_content(chunks.position_in_chunk(location));
    }

    char traversed_base(int64_t location, ll position) const {
        auto &chunk = chunks.at(location);
        return chunk.traversed_base(chunks.position_in_chunk(location), position);
    }

    ll new_relative_position(int64_t location, ll position) {
        auto &chunk = chunks.at(location);
        return chunk.new_relative_position(chunks.position_in_chunk(location), position);
    }

    void insert(int64_t location, ll position, char symbol) {
        auto &chunk = chunks[location];
        chunk.insert(chunks.position_in_chunk(location), position, symbol);
    }


    ChunkedDenseHashMap<EntryT,BitVector,RankSupport> chunks;
    shared_ptr<const BetterDBGSuccinct> graph;
};



//template <typename EntryT=wavelet_tree_dyn>
template <typename EntryT=SolidDynamicRoutingTable<>,typename BitVector=sdsl::bit_vector,typename RankSupport=typename BitVector::rank_1_type>
class DynamicRoutingTable : public TransformationsEnabler<DynamicRoutingTableCore<EntryT,BitVector,RankSupport>> {
    using TransformationsEnabler<DynamicRoutingTableCore<EntryT,BitVector,RankSupport>>::TransformationsEnabler;
};

#endif //METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
