//
//  path_encoder_baseline.cpp
//  lossless_dbg
//
//  Created by Jan Studený on 20/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#include "path_encoder.hpp"
#include "utils.hpp"
#include <iostream>
#include <set>
#include <map>

using namespace std;
// todo find tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this


class PathDatabaseBaseline : public PathDatabase {
public:
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaseline(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabase(graph), graph(*graph_) {
        
    }
    
    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering
        sort(all(sequences));
        for(auto& sequence : sequences) {
            increment_joins(sequence);
        }
        return {};
    }
    
    void increment_joins(string sequence) {
        auto kmer = sequence.substr(0,graph.get_k());
        joins[kmer]['$']++;
        
    }
    
    void compress() {
        
    }
    
    size_t num_paths() const override;
    
    std::string decode(PathDatabase::path_id path) const override;
    
    std::vector<path_id> get_paths_going_through(const std::string &str) const override;
    
    std::vector<path_id> get_paths_going_through(PathDatabase::node_index node) const override;
    
    PathDatabase::node_index get_next_node(PathDatabase::node_index node, PathDatabase::path_id path) const override;
    
    PathDatabase::node_index get_next_consistent_node(const std::string &history) const override;
    
private:
    // denote how many reads are joining from every branch (ATCGN$) ($ denotes start of a new read)
    std::map<string,map<char,int>> joins;
    // denote where the reads should go (ATCGN$) ($ denodes the end of particular read)
    std::map<string,map<int,char>> splits;
    
    std::set<string> additional_bifurcations;
    
    const DeBruijnGraph & graph;
    
};

//
//int main() {
//
//}
