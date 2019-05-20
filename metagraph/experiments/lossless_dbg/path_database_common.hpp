//
// Created by Jan Studený on 2019-05-20.
//

#ifndef METAGRAPH_PATH_DATABASE_COMMON_HPP
#define METAGRAPH_PATH_DATABASE_COMMON_HPP

#include <iostream>
#include <set>
#include <map>
#include <tsl/hopscotch_set.h>
#include <progress_bar.hpp>
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

template<typename GraphT=DBGSuccinct,
        typename RoutingTable=DynamicRoutingTable,
        typename IncomingTable=DynamicIncomingTable<>>
class PathDatabaseCommon : public PathDatabase<pair<node_index,int>,GraphT> {
public:
    static constexpr bool use_char = is_same<typename IncomingTable::edge_identifier_t,char>::value;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseCommon(std::shared_ptr<const GraphT> graph) :
            PathDatabase<pair<node_index,int>,GraphT>(graph),
            graph(*(this->graph_)),
            incoming_table(*(this->graph_)),
            routing_table(nullptr) {}


    explicit PathDatabaseCommon(const vector<string> &filenames,
                                 size_t k_kmer = 21 /* default kmer */) :
            PathDatabase<pair<node_index,int>,GraphT>(filenames, k_kmer),
            graph(*(this->graph_)),
            incoming_table(*(this->graph_)),
            routing_table()
    {}

    virtual ~PathDatabaseCommon() {}


    virtual bool node_is_split(node_index node) const = 0;
    virtual bool node_is_join(node_index node) const = 0;

    std::string decode(path_id path) const override {
        auto prev_node = 0;
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;
        string sequence_path = kmer;
        int relative_position = path.second;

        int kmer_position = 0;
        char base = '\0';
        char encoded_base = '\0';
        while (true) {
            if (node_is_split(node)) {
                encoded_base = routing_table.get(node,relative_position); // maybe different
                base = routing_table.traversed_base(node,relative_position);
                relative_position = routing_table.new_relative_position(node,relative_position);
                if (base != encoded_base) {
                    cout << encoded_base << base << endl;
                }
            }
            else {
                assert(graph.outdegree(node) == 1);
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
                encoded_base = base; // same
            }
            assert(base);
            if (base == '$') break;
            prev_node = node;
            node = graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,encoded_base); // 1 times base
            sequence_path.append(1,base); // 1 times base

            if (node_is_join(node)) {
                if constexpr (use_char) {
                    auto join_symbol = sequence_path[kmer_position-1];
                    relative_position += incoming_table.branch_offset(node,join_symbol);
                }
                else {
                    relative_position += incoming_table.branch_offset(node,prev_node);
                }
            }
        }
        return sequence;
    }

    std::vector<path_id> get_paths_going_through(const std::string &str) const override { throw std::runtime_error("not implemented"); };

    std::vector<path_id> get_paths_going_through(node_index node) const override { throw std::runtime_error("not implemented"); };

    node_index get_next_node(node_index node, path_id path) const override { throw std::runtime_error("not implemented"); };

    node_index get_next_consistent_node(const std::string &history) const override { throw std::runtime_error("not implemented"); };

    void serialize(const fs::path& folder) const {};

protected:
    json statistics;

    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)

    RoutingTable routing_table;
    IncomingTable incoming_table;

    const GraphT & graph;

};



#endif //METAGRAPH_PATH_DATABASE_COMMON_HPP
