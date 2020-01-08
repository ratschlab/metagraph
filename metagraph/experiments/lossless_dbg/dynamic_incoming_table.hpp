//
// Created by Jan Studený on 2019-05-06.
//

#ifndef __DYNAMIC_INCOMING_TABLE_HPP__
#define __DYNAMIC_INCOMING_TABLE_HPP__

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <tsl/hopscotch_map.h>

#include "graph_patch.hpp"
#include <prettyprint.hpp>
//#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "solid_dynamic_incoming_table.hpp"

using default_bit_vector = bit_vector_small;


template <typename EntryT=SolidDynamicIncomingTable<>,typename BitVector=sdsl::bit_vector,typename RankSupport=typename BitVector::rank_1_type>
class DynamicIncomingTable {
public:
    using edge_identifier_t = char;
    DynamicIncomingTable() = default;
    DynamicIncomingTable(shared_ptr<const DBGSuccinct> graph) : graph(graph) {}
    DynamicIncomingTable(shared_ptr<const DBGSuccinct> graph, BitVector* is_element,RankSupport* rank_element, ll chunks = DefaultChunks)
            : graph(graph), chunks(is_element,rank_element,chunks) {
    }

    bool is_join(int64_t location) const {
        auto &chunk = chunks.at(location);
        return chunk.is_join(chunks.position_in_chunk(location));
    }

    ll branch_offset(int64_t location, edge_identifier_t incoming) const {
        auto &chunk = chunks.at(location);
        return chunk.branch_offset(chunks.position_in_chunk(location), incoming);
    }
    
    ll branch_size(int64_t location, edge_identifier_t incoming) const {
        auto &chunk = chunks.at(location);
        return chunk.branch_size(chunks.position_in_chunk(location), incoming);
    }

    ll size(int64_t location) const {
        auto &chunk = chunks.at(location);
        return chunk.size(chunks.position_in_chunk(location));
    }

    ll branch_offset_and_increment(int64_t location, edge_identifier_t incoming) {
        auto &chunk = chunks[location];
        return chunk.branch_offset_and_increment(chunks.position_in_chunk(location), incoming);
    }

    string print_content(int64_t location) const {
        auto &chunk = chunks.at(location);
        return chunk.print_content(chunks.position_in_chunk(location));
    }

    shared_ptr<const DBGSuccinct> graph;
    ChunkedDenseHashMap<EntryT,BitVector,RankSupport> chunks;
};


#endif // __DYNAMIC_INCOMING_TABLE_HPP__
