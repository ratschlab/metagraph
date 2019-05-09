//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_wavelet_deprecated_hpp
#define path_database_baseline_wavelet_deprecated_hpp

#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "path_database_wavelet.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include <iostream>
#include <set>
#include <map>
#include "alphabets.hpp"
#include "routing_table.hpp"


#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"

using namespace std;
using alphabets::log2;

using GraphT = DBGSuccinct;
class PathDatabaseBaselineWaveletDeprecated : public PathDatabaseDynamic<GraphT> {
public:
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaselineWaveletDeprecated(std::shared_ptr<const GraphT> graph) : PathDatabaseDynamic(graph),
                                                                              routing_table(sizeof(RoutingTableAlphabet))
                                                                              {}

    PathDatabaseBaselineWaveletDeprecated(const vector<string> &raw_reads,
                                size_t k_kmer = 21 /* default */) : PathDatabaseDynamic(raw_reads,k_kmer),
                                                                    routing_table(sizeof(RoutingTableAlphabet))
                                                                    {}


    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {

        vector<path_id> encoded = PathDatabaseDynamic::encode(sequences);
            // add additional bifurcation
        construct_routing_table();
        construct_edge_multiplicity_table();
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
            if (PathDatabaseDynamic::node_is_split(node)) {
                auto dynamic_table = PathDatabaseDynamic::routing_table;
                for(int i=0;i<dynamic_table.size(node);i++) {
                    routing_table_array.push_back(dynamic_table.get(node,i));
                }
            }
        }
        vector<int> routing_table_array_encoded;
        transform(all(routing_table_array),back_inserter(routing_table_array_encoded),[](char c) {
            return RoutingTableInverseAlphabet.at(c);}
            );
        routing_table = wavelet_tree_stat(sizeof(RoutingTableAlphabet),routing_table_array_encoded);
    }

    void construct_edge_multiplicity_table() {
        sdsl::bit_vector is_join_node(graph.num_nodes());
        for(int node=1;node<=graph.num_nodes();node++) {
            is_join_node[node-1] = PathDatabaseDynamic::node_is_join(node);
            if (PathDatabaseDynamic::node_is_join(node)) {
                for(int rc=0;rc<'N'_rc;rc++) { // don't need to store last branch as we only compute prefix sum excluding
                                               // the branch which we came from (N in this case)
                    edge_multiplicity_table.push_back(PathDatabaseDynamic::incoming_table.branch_size(node,tochar(rc)));
                }
            }
        }
        joins = bit_vector_small(is_join_node);

    }


    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;

        int relative_starting_position = path.second;
        int relative_position = branch_starting_offset(node,'$') + relative_starting_position;

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
                assert(graph.outdegree(node) == 1);
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
//                auto tv = PathDatabaseDynamic::branch_starting_offset(node,join_symbol);
//                auto cv = branch_starting_offset(node,join_symbol);
//                assert(tv==cv);
                relative_position += branch_starting_offset(node,join_symbol);
            }
        }

        return sequence;
    }

    vector<string> decode_all_reads() const {
        auto reads = vector<string>();
        for(node_index node=1;node<=graph.num_nodes();node++) {
            auto count = number_of_reads_starting_at_node(node);
            for(auto relative_index=0;relative_index<count;relative_index++) {
                reads.push_back(decode({node,relative_index}));
            }
        }
        return reads;
    }

    int number_of_reads_starting_at_node(node_index node) const {
        int result = 0;
        if (node_is_join(node)) {
            int starting_offset = (joins.rank1(node-1)-1)*'N'_rc;
            result = edge_multiplicity_table[starting_offset + '$'_rc];
        }
        return result;
    }

    bool node_is_join(node_index node) const {
        return joins[node-1]; // 0 based
    }

    bool node_is_split(node_index node) const {
        auto offset = routing_table.select('#'_rc,node);
        // is not the last element of routing table and next character is not starting of new node
        return routing_table.size() != (offset + 1) and routing_table[offset+1] != '#'_rc;
    }

    int branch_starting_offset(node_index node,char branch_label) const {
        //node-1 as we are indexing from 0
        //rank1 - 1 because the rank is inclusive
        // * 'N' as we write multiple values
        int starting_offset = (joins.rank1(node-1)-1)*'N'_rc;
        int result = 0;
        for(int previous_branch=0;previous_branch<rc(branch_label);previous_branch++) {
            result += edge_multiplicity_table[starting_offset+previous_branch];
        }
        return result;
    }


    void serialize(const fs::path& folder) const {
        fs::create_directories(folder / "xxx.bin");
        ofstream edge_multiplicity_file(folder / "edge_multiplicity.bin", ios_base::trunc | ios_base::out);
        ofstream routing_table_file(folder / "routing_table.bin", ios_base::trunc | ios_base::out);
        ofstream joins_file(folder / "joins.bin", ios_base::trunc | ios_base::out);
        string graph_filename = folder / "graph.bin";

        ::serialize(edge_multiplicity_file,edge_multiplicity_table);
        routing_table.serialize(routing_table_file);
        joins.serialize(joins_file);
        graph.serialize(graph_filename);
    }

    static PathDatabaseBaselineWaveletDeprecated deserialize(const fs::path& folder) {
        ifstream edge_multiplicity_file(folder / "edge_multiplicity.bin");
        ifstream routing_table_file(folder / "routing_table.bin");
        ifstream joins_file(folder / "joins.bin");
        string graph_filename = folder / "graph.bin";

        auto graph = std::shared_ptr<GraphT>{
                new DBGSuccinct(21)
                };
        graph->load(graph_filename);
        auto db = PathDatabaseBaselineWaveletDeprecated(graph);
        ::deserialize(edge_multiplicity_file,db.edge_multiplicity_table);
        db.routing_table.load(routing_table_file);
        db.joins.load(joins_file);
        return db;
    }

private:
    wavelet_tree_stat routing_table;
    bit_vector_small joins;
    vector<int> edge_multiplicity_table;
};

#endif /* path_database_baseline_hpp */
