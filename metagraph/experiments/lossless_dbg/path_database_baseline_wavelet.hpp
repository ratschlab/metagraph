//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_wavelet_hpp
#define path_database_baseline_wavelet_hpp

#include "path_database.hpp"
#include "utils.hpp"
#include <iostream>
#include <set>
#include <map>
#include "alphabets.hpp"

#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"

using namespace std;
using alphabets::log2;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this


using routing_character_t = int;
const char RoutingTableAlphabet[] = {'$','A','C','G','T','N','#','?'};
// improvement (constexpr use https://github.com/serge-sans-paille/frozen)
const map<char,int> RoutingTableInverseAlphabet = {{'$',0},{'A',1},{'C',2},{'G',3},{'T',4},{'N',5},{'#',6},{'?',7}};
const auto& rte2int = RoutingTableInverseAlphabet;

int operator""_rc(char c) {
    return RoutingTableInverseAlphabet.at(c);
}
int rc(char c) {
    return RoutingTableInverseAlphabet.at(c);
}

char tochar(routing_character_t rc) {
    return RoutingTableAlphabet[rc];
}



class PathDatabaseBaselineWavelet : public PathDatabaseBaseline {
public:
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaselineWavelet(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabaseBaseline(graph),
                                                                              routing_table(log2(sizeof(RoutingTableAlphabet)))
                                                                              {}

    PathDatabaseBaselineWavelet(const vector<string> &raw_reads, size_t k_kmer) : PathDatabaseBaseline(raw_reads,k_kmer),
                                                                                  routing_table(log2(sizeof(RoutingTableAlphabet)))
                                                                                  {}

    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {

        vector<path_id> encoded = PathDatabaseBaseline::encode(sequences);
            // add additional bifurcation
        construct_routing_table();
        return encoded;
    }

//    void routing_table_add_delimiters(const std::vector<std::string> &sequences) {
//        split_size.resize(graph.num_nodes());
//        for(auto& sequence : sequences) {
//            increment_splits(sequence);
//        }
//
//        int routing_table_entries = 1;
//        for(int node=1;node<=graph.num_nodes();node++) {
//            routing_table_entries++; // for # symbol
//            routing_table_entries += split_size[node];
//        }
//
//        vector<int> table(routing_table_entries,'?'_rc);
//    }
//
//    void increment_splits(const string& sequence) {
//        auto kmer = sequence.substr(0,graph.get_k());
//        auto node = graph.kmer_to_node(kmer);
//        int kmer_position = 0;
//        for(auto& base : sequence.substr(graph.get_k())) {
//            if (node_is_split(node)) {
//                split_size[node]++;
//            }
//            node = graph.traverse(node,base);
//            kmer_position++;
//        }
//        split_size[node]++; // additional split (end)
//    }

    void construct_routing_table() {
        vector<char> routing_table_array;
        for(int node=1;node<=graph.num_nodes();node++) {
            routing_table_array.push_back('#');
            if (splits.count(node)) {
                for(auto& choice : splits[node]) {
                    routing_table_array.push_back(choice);
                }
            }
        }
        vector<int> routing_table_array_encoded;
        transform(all(routing_table_array),back_inserter(routing_table_array_encoded),[](char c) {
            return RoutingTableInverseAlphabet.at(c);}
            );
        routing_table = wavelet_tree_stat(sizeof(RoutingTableAlphabet),routing_table_array_encoded);
    }


    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;

        int relative_starting_position = path.second;
        int relative_position = branch_starting_offset(joins.at(node),'$') + relative_starting_position;

        int kmer_position = 0;
        int base;
        while(true) {
            if (node_is_split(node)) {
                auto routing_table_block = routing_table.select('#'_rc, node)+1;
                auto absolute_position = routing_table_block+relative_position;
                base = routing_table[absolute_position];
                auto occurrences_of_base_before_block = routing_table.rank(base,routing_table_block-1);

                // checkers
//                auto& routing_table_naive = splits.at(node);
//                auto rt_index = routing_table_naive.begin();
//                advance(rt_index,relative_position);
//                char base_check = *rt_index;
//                auto new_relative_position = rank(routing_table_naive,tochar(base),relative_position)-1;

                // -1 as rank is inclusive of the absolute position
                auto rank_of_base = routing_table.rank(base,absolute_position) - occurrences_of_base_before_block - 1;

//                assert(base_check == tochar(base));
//                assert(rank_of_base == new_relative_position);

                relative_position = rank_of_base;
            }
            else {
                assert(graph.is_single_outgoing(node));
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = rc(edge_label);});
            }
            if (base == '$'_rc) break;
            node = graph.traverse(node,tochar(base));
            assert(node);
            kmer_position++;
            sequence.append(1,tochar(base)); // 1 times base
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += branch_starting_offset(joins.at(node),join_symbol);
            }
        }

        return sequence;
    }


private:
    wavelet_tree_stat routing_table;
};

#endif /* path_database_baseline_hpp */
