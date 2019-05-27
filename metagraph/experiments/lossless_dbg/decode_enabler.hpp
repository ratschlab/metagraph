//
// Created by Jan Studený on 2019-05-20.
//

#ifndef METAGRAPH_DECODE_ENABLER_HPP
#define METAGRAPH_DECODE_ENABLER_HPP

#include <iostream>
#include <set>
#include <map>
#include <tsl/hopscotch_set.h>
#include <optional>

#include "path_database.hpp"
#include "dynamic_routing_table.hpp"
#include "dynamic_incoming_table.hpp"
#include "utils.hpp"
#include "unix_tools.hpp"
#include "threading.hpp"

#define USE_LOCKS

// TODO: Never use 'using namespace std;' in .hpp files
// todo find a tool that removes this relative namespacing issue
using namespace std;

// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this
using node_index = DeBruijnGraph::node_index;

using path_id = pair<node_index,int>;

template<typename Database>
class DecodeEnabler : public Database {
public:
    using Database::Database;
    using edge_id_t = typename decltype(Database::incoming_table)::edge_identifier_t;
    static constexpr bool use_char = is_same<edge_id_t,char>::value;
    static constexpr edge_id_t origin_node_symbol = use_char ? '$' : 0;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k

    std::string decode(path_id path) const {
        auto prev_node = 0;
        auto node = path.first;
        auto kmer = this->graph.get_node_sequence(node);
        string sequence = kmer;
        string sequence_path = kmer;
        int relative_position = path.second;

        int kmer_position = 0;
        char base = '\0';
        char encoded_base = '\0';
        while (true) {
            if (this->node_is_split(node)) {
                encoded_base = this->routing_table.get(node,relative_position); // maybe different
                base = this->routing_table.traversed_base(node,relative_position);
                relative_position = this->routing_table.new_relative_position(node,relative_position);
                if (base != encoded_base) {
                    cout << encoded_base << base << endl;
                    this->routing_table.print_content(node);
                    this->incoming_table.print_content(node);
                }
            }
            else {
                assert(this->graph.outdegree(node) == 1);
                this->graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
                encoded_base = base; // same
            }
            assert(base);
            if (base == '$') break;
            prev_node = node;
            node = this->graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,encoded_base); // 1 times base
            sequence_path.append(1,base); // 1 times base

            if (this->node_is_join(node)) {
                if constexpr (use_char) {
                    auto join_symbol = sequence_path[kmer_position-1];
                    relative_position += this->incoming_table.branch_offset(node,join_symbol);
                }
                else {
                    relative_position += this->incoming_table.branch_offset(node,prev_node);
                }
            }
        }
        return sequence;
    }

    vector<string> decode_all_reads() const {
        auto reads = vector<string>();
        for(node_index node=1;node<=this->graph.num_nodes();node++) {
            auto count = number_of_reads_starting_at_node(node);
            for(auto relative_index=0;relative_index<count;relative_index++) {
                reads.push_back(decode({node,relative_index}));
            }
        }
        return reads;
    }

    vector<string> decode_all_reads_inverse() const {
        auto reads = vector<string>();
        for(node_index node=1;node<=this->graph.num_nodes();node++) {
            auto count = number_of_reads_ending_at_node(node);
            for(auto read_index=0;read_index<count;read_index++) {
                auto relative_index = this->routing_table.select(node,read_index+1,'$');
                auto id = get_global_path_id(node,relative_index);
                reads.push_back(decode(id));
            }
        }
        return reads;
    }

    int number_of_reads_starting_at_node(node_index node) const {
        int result = 0;
        if (this->node_is_join(node)) {
            result = this->incoming_table.branch_size(node,origin_node_symbol);
        }
        return result;
    }

    int number_of_reads_ending_at_node(node_index node) const {
        int result = 0;
        if (this->node_is_split(node)) {
            result = this->routing_table.rank(node,this->routing_table.size(node),'$');
        }
        return result;
    }

    path_id get_global_path_id(node_index node, int relative_position) const {
        node_index prev_node = 0;
        int prev_offset = 0;
        if (this->node_is_join(node)) {
            this->graph.call_incoming_kmers_mine(node,[&](node_index possible_node,char c) {
                auto offset = use_char ?
                        this->incoming_table.branch_offset(node,c) :
                              this->incoming_table.branch_offset(node,possible_node);
                if (offset <= relative_position && offset >= prev_offset) {
                    prev_node = possible_node;
                    prev_offset = offset;
                }
            });
            relative_position -= prev_offset;
        }
        else {
            assert(this->graph.indegree(node) == 1);
            this->graph.call_incoming_kmers_mine(node,[&](node_index possible_node,char c) {
                prev_node = possible_node;
            });
        }
        if (!prev_node) {
            assert(this->is_valid_path_id({node,relative_position}));
            return {node,relative_position};
        }
        assert(prev_node);
        if (this->node_is_split(prev_node)) {
            relative_position = this->routing_table.select(prev_node,relative_position+1,this->node_get_last_char(node));// +1 as relative_position is 0-based
        }
        return get_global_path_id(prev_node,relative_position);
    }

    bool is_valid_path_id(path_id path_id) const {
        return this->node_is_join(path_id.first) && path_id.second < this->incoming_table.branch_size(path_id.first,origin_node_symbol);
    }

    char node_get_last_char(node_index node) const {
        auto kmer = this->graph.get_node_sequence(node);
        return kmer[kmer.size()-1];
    }

    char node_get_first_char(node_index node) const {
        auto kmer = this->graph.get_node_sequence(node);
        return kmer.front();
    }

};



#endif //METAGRAPH_DECODE_ENABLER_HPP
