//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_GRAPH_PATCH_HPP
#define METAGRAPH_GRAPH_PATCH_HPP

#include "path_database.hpp"
#include "path_database_baseline.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include <iostream>
#include <set>
#include <functional>
#include <map>
#include "alphabets.hpp"
#include "routing_table.hpp"
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
};
#define DBGSuccinct BetterDBGSuccinct

#endif //METAGRAPH_GRAPH_PATCH_HPP
