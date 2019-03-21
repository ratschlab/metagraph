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
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaseline(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabase(graph), graph(*graph_) {
        
    }
    
    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering
        //sort(all(sequences));
        for(auto& sequence : sequences) {
            // add additional bifurcation
            additional_joins.insert(graph.kmer_to_node(sequence.substr(0,graph.get_k())));
            additional_splits.insert(graph.kmer_to_node(sequence.substr(sequence.length()-graph.get_k())));
        }
        
        for(auto& sequence : sequences) {
            route_sequence(sequence);
        }
//        for(auto& sequence : sequences) {
//            populate_splits(sequence);
//        }
        return {};
    }
    
    // when parallelizing beware that we need to account also artificial joins that will be populatated gradually
    // optimize : don't always copy the string when creating a substring (simple replacement with string_view doesn't work as maps don't work with them)
    // use node instead of kmer
//    void increment_joins(const string& sequence) {
//        auto kmer = sequence.substr(0,graph.get_k());
//        additional_bifurcations.insert(kmer);
//        joins[kmer]['~']++;
//        auto node = graph.kmer_to_node(kmer);
//        int kmer_position = 0;
//        for(auto& base : sequence.substr(graph.get_k())) {
//            node = graph.traverse(node,base);
//            kmer_position++;
//            kmer = sequence.substr(kmer_position,graph.get_k());
//            if (node_is_join(kmer, node)) {
//                joins[kmer][sequence[kmer_position-1]]++;
//            }
//        }
//    }
    
    void route_sequence(const string& sequence) {
        auto kmer = sequence.substr(0,graph.get_k());
        auto node = graph.kmer_to_node(kmer);
        // always putting new read above all other reads
        int relative_position = offset_for_symbol(joins[node],'~');
        joins[node]['~']++;

        int kmer_position = 0;
        for(auto& base : sequence.substr(graph.get_k())) {
            if (node_is_split(node)) {
                auto& routing_table = splits[node];
                auto rt_index = routing_table.begin();
                advance(rt_index,relative_position);
                routing_table.insert(rt_index,base);
                relative_position = rank(routing_table,base,relative_position)-1;
            }
            node = graph.traverse(node,base);
            kmer_position++;
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += offset_for_symbol(joins[node],join_symbol);
                joins[node][join_symbol]++;
            }
        }
        
        auto& routing_table = splits[node];
        auto rt_index = routing_table.begin();
        advance(rt_index,relative_position);
        routing_table.insert(rt_index,'~');
        
    }
    
//    void populate_splits(const string& sequence) {
//        auto kmer = sequence.substr(0,graph.get_k());
//        auto node = graph.kmer_to_node(kmer);
//        // always putting new read above all other reads
//        int relative_position = number_of_reads_at_join(kmer)
//                                - joins[kmer]['~'] // excluding the newly created kmers
//                                + number_of_reads_starting_at_kmer[kmer]; // + relative offset
//        // todo: currently the relative offset is based on ordering of the reads in file
//        // better to sort them
//        number_of_reads_starting_at_kmer[kmer]++;
//        for(auto& base : sequence.substr(graph.get_k())) {
//
//        }
//    }
    
//    int number_of_reads_at_join(const string& kmer) {
//        int result = 0;
//        for(auto& [branch_base_character,count] : joins[kmer]) {
//            result += count;
//        }
//        return result;
//    }
    
    
    // returns the number of reads that were already stored and should have lower index
    int offset_for_symbol(const map<char,int>& join,char symbol) {
        int result = 0;
        for(auto&[base,count] : join) {
            result += count;
            if (base == symbol) break;
        }
        return result;
    }
    
    bool node_is_join(node_index node) const {
        return !graph.is_single_incoming(node) or additional_joins.count(node);
    }
    bool node_is_split(node_index node) const {
        return !graph.is_single_outgoing(node) or additional_splits.count(node);
    }
    
    int rank(const routing_table_t& routing_table, char symbol, int position) {
        assert(position < routing_table.size());
        int result = 0;
        int i = 0;
        for(auto it=begin(routing_table);i<=position;it++,i++) {
            result += *it == symbol;
        }
        return result;
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
    // denote how many reads are joining from every branch (ATCGN~) (~ denotes start of a new read)
    std::map<node_index,map<char,int>> joins;
    // denote where the reads should go (ATCGN~) (~ denodes the end of particular read)
    
    std::map<node_index,routing_table_t> splits;
    
    std::set<node_index> additional_joins;
    std::set<node_index> additional_splits;
    
    const DeBruijnGraph & graph;
    
};


int main() {

}
