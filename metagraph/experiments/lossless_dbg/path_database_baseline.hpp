//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_hpp
#define path_database_baseline_hpp

#include "path_database.hpp"
#include "utils.hpp"
#include <iostream>
#include <set>
#include <map>

#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"

using namespace std;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this
using node_index = DeBruijnGraph::node_index;
//using path_id = pair<node_index,int>;

class PathDatabaseBaseline : public PathDatabase<pair<node_index,int>> {
public:
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseBaseline(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabase(graph), graph(*graph_) {}

    explicit PathDatabaseBaseline(const vector<string> &raw_reads,
                 size_t k_kmer = 21 /* default kmer */) : PathDatabase(raw_reads,k_kmer), graph(*graph_) {}

    node_index starting_node(const string& sequence) {
        return graph.kmer_to_node(sequence.substr(0,graph.get_k()));
    }

    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering
        vector<path_id> encoded;
        for(auto& sequence : sequences) {
            // add additional bifurcation
            additional_joins.insert(starting_node(sequence));
            additional_splits.insert(graph.kmer_to_node(sequence.substr(sequence.length()-graph.get_k())));
        }

        for(auto& sequence : sequences) {
            int relative_order = route_sequence(sequence);
            encoded.push_back({starting_node(sequence),relative_order});
            encoded_paths++;
        }

        return encoded;
    }


    int route_sequence(const string& sequence) {
        auto kmer = sequence.substr(0,graph.get_k());
        auto node = graph.kmer_to_node(kmer);
        // always putting new read above all other reads
        int relative_starting_position = joins[node]['$'];
        int relative_position = branch_starting_offset(joins[node],'$') + relative_starting_position;
        joins[node]['$']++;

        int kmer_position = 0;
        for(auto& base : sequence.substr(graph.get_k())) {
            if (node_is_split(node)) {
                auto& routing_table = splits[node];
                auto rt_index = routing_table.begin();
                assert(relative_position <= routing_table.size());
                advance(rt_index,relative_position);
                routing_table.insert(rt_index,base);
                relative_position = rank(routing_table,base,relative_position)-1;
            }
            node = graph.traverse(node,base);
            kmer_position++;
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += branch_starting_offset(joins[node],join_symbol);
                joins[node][join_symbol]++;
            }
        }

        auto& routing_table = splits[node];
        auto rt_index = routing_table.begin();
        advance(rt_index,relative_position);
        routing_table.insert(rt_index,'$');

        return relative_starting_position;
    }

    // returns the number of reads that were already stored and should have lower index
    int offset_for_symbol(const map<char,int>& join,char symbol) const {
        int result = 0;
        for(auto&[base,count] : join) {
            result += count;
            if (base >= symbol) break;
        }
        return result;
    }

    int branch_starting_offset(const map<char,int>& join,char symbol) const {
        int result = 0;
        for(auto&[base,count] : join) {
            if (base >= symbol) break;
            result += count;
        }
        return result;
    }

    bool node_is_join(node_index node) const {
        return graph.indegree(node) > 1 or additional_joins.count(node);
    }
    bool node_is_split(node_index node) const {
        return graph.outdegree(node) > 1 or additional_splits.count(node);
    }

    int rank(const routing_table_t& routing_table, char symbol, int position) const {
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

    size_t num_paths() const override {
        return encoded_paths;
    };

    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;

        int relative_starting_position = path.second;
        int relative_position = branch_starting_offset(joins.at(node),'$') + relative_starting_position;

        int kmer_position = 0;
        char base;
        while(true) {
            if (node_is_split(node)) {
                auto& routing_table = splits.at(node);
                auto rt_index = routing_table.begin();
                advance(rt_index,relative_position);
                base = *rt_index;
                relative_position = rank(routing_table,base,relative_position)-1;
            }
            else {
                assert(graph.outdegree(node) == 1);
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
            }
            if (base == '$') break;
            node = graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,base); // 1 times base

            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
                relative_position += branch_starting_offset(joins.at(node),join_symbol);
            }
        }

        return sequence;
    }

    std::vector<path_id> get_paths_going_through(const std::string &str) const override {};

    std::vector<path_id> get_paths_going_through(PathDatabase::node_index node) const override {};

    PathDatabase::node_index get_next_node(PathDatabase::node_index node, path_id path) const override {};

    PathDatabase::node_index get_next_consistent_node(const std::string &history) const override {};

protected:
    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int encoded_paths = 0;
    std::map<node_index,map<char,int>> joins;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)

    std::map<node_index,routing_table_t> splits;

    std::set<node_index> additional_joins;
    std::set<node_index> additional_splits;

    const DeBruijnGraph & graph;

};

#endif /* path_database_baseline_hpp */
