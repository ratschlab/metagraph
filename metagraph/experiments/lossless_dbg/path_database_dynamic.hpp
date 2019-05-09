//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_hpp
#define path_database_baseline_hpp

#include <iostream>
#include <set>
#include <map>

#include "path_database.hpp"
#include "dynamic_routing_table.hpp"
#include "dynamic_incoming_table.hpp"
#include "utils.hpp"


using namespace std;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this
using node_index = DeBruijnGraph::node_index;

using path_id = pair<node_index,int>;

template<typename GraphT=DBGSuccinct>
class PathDatabaseDynamic : public PathDatabase<pair<node_index,int>,GraphT> {
public:
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseDynamic(std::shared_ptr<const GraphT> graph) :
                PathDatabase<pair<node_index,int>,GraphT>(graph),
                graph(*(this->graph_)),
                incoming_table(*(this->graph_)) {}

    explicit PathDatabaseDynamic(const vector<string> &raw_reads,
                 size_t k_kmer = 21 /* default kmer */) :
                         PathDatabase<pair<node_index,int>,GraphT>(raw_reads,k_kmer),
                         graph(*(this->graph_)),
                         incoming_table(*(this->graph_)) {}

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
#ifdef MEMOIZE
        is_split = vector<bool>(graph.num_nodes()+1);
        is_join = vector<bool>(graph.num_nodes()+1);
        for (int node=1;node<=graph.num_nodes();node++) {
            is_split[node] = node_is_split_raw(node);
            is_join[node] = node_is_join_raw(node);
        }
#endif

        for(auto& sequence : sequences) {
            int relative_order = route_sequence(sequence);
            encoded.push_back({starting_node(sequence),relative_order});
            encoded_paths++;
        }

        return encoded;
    }


    int route_sequence(const string& sequence) {
        // always putting new read above all other reads
        int relative_position = 0;
        int relative_starting_position = -1;
        int kmer_left_border = 0;
        int kmer_right_border = graph.get_k();
        bool first_node = 1;
        graph.map_to_nodes(sequence,[&](node_index node){
            if (this->node_is_join(node)) {
                auto join_symbol = kmer_left_border ? sequence[kmer_left_border-1] : '$';
                relative_position += incoming_table.branch_offset(node,join_symbol);
                if (join_symbol == '$') {
                    relative_starting_position = incoming_table.branch_size(node,'$');
                    relative_position += relative_starting_position;
                }
                incoming_table.increment(node,join_symbol);
            }
            if (this->node_is_split(node)) {
                auto split_symbol = kmer_right_border < sequence.size() ? sequence[kmer_right_border] : '$';
                routing_table.insert(node,relative_position,split_symbol);
                relative_position = routing_table.rank(node,split_symbol,relative_position);
            }
            kmer_left_border++;
            kmer_right_border++;
        });

        return relative_starting_position;
    }

#ifdef MEMOIZE

    bool node_is_join_raw(node_index node) const {
        return graph.indegree(node) > 1 or additional_joins.count(node);
    }

    bool node_is_split_raw(node_index node) const {
        return graph.outdegree(node) > 1 or additional_splits.count(node);
    }

    bool node_is_join(node_index node) const {
        return is_join[node];
    }

    bool node_is_split(node_index node) const {
        return is_split[node];
    }

#else
    bool node_is_join(node_index node) const {
        return graph.indegree(node) > 1 or additional_joins.count(node);
    }
    bool node_is_split(node_index node) const {
        return graph.outdegree(node) > 1 or additional_splits.count(node);
    }
#endif

    void compress() {

    }

    size_t num_paths() const override {
        return encoded_paths;
    };

    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;
        int relative_position = path.second;

        int kmer_position = 0;
        char base = '\0';
        while(true) {
            if (node_is_split(node)) {
                base = routing_table.get(node,relative_position);
                relative_position = routing_table.rank(node,base,relative_position);
            }
            else {
                assert(graph.outdegree(node) == 1);
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
            }
            assert(base);
            if (base == '$') break;
            node = graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,base); // 1 times base

            if (node_is_join(node)) {
                auto join_symbol = sequence[kmer_position-1];
                relative_position += incoming_table.branch_offset(node,join_symbol);
            }
        }

        return sequence;
    }

    std::vector<path_id> get_paths_going_through(const std::string &str) const override {};

    std::vector<path_id> get_paths_going_through(node_index node) const override {};

    node_index get_next_node(node_index node, path_id path) const override {};

    node_index get_next_consistent_node(const std::string &history) const override {};

    void serialize(const fs::path& folder) const {};

protected:
    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)

    DynamicRoutingTable routing_table;
    DynamicIncomingTable<> incoming_table;

    std::set<node_index> additional_joins;
    std::set<node_index> additional_splits;

    const GraphT & graph;

#ifdef MEMOIZE
    vector<bool> is_join;
    vector<bool> is_split;
#endif
};

#endif /* path_database_baseline_hpp */
