//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_wavelet_hpp
#define path_database_baseline_wavelet_hpp

#include <iostream>
#include <set>
#include <functional>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>

#include "utils.hpp"
#include "alphabets.hpp"
#include "cxx-prettyprint.hpp"

#include "path_database.hpp"
#include "path_database_dynamic.hpp"
#include "routing_table.hpp"
#include "incoming_table.hpp"
#include "utilities.hpp"
#include "query_enabler.hpp"

#include "graph_patch.hpp"
//#define CHECK_CORECTNESS 1

#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"

const unsigned int STATS_JOINS_HISTOGRAM (1u << 0u);
const unsigned int STATS_SPLITS_HISTOGRAM (1u << 1u);

using namespace std;
using alphabets::log2;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this



//template <class Wavelet = sdsl::wt_rlmn<sdsl::sd_vector<>>
template<class Wavelet = sdsl::wt_rlmn<>,class BitVector=default_bit_vector>
class PathDatabaseWaveletCore : public PathDatabaseDynamicCore<> {
public:
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseWaveletCore(std::shared_ptr<const DBGSuccinct> graph) : PathDatabaseDynamicCore(graph),
                                                                              incoming_table(*graph)
                                                                              {}

    PathDatabaseWaveletCore(const vector<string> &filenames,
                        size_t kmer_length = 21 /* default */) : PathDatabaseDynamicCore(filenames,kmer_length),
                                                            incoming_table(graph) {}



    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        auto encoded = PathDatabaseDynamicCore::encode(sequences);

        // convert dynamic_(routing_table/incoming_table) to routing_table/incoming_table
        construct_routing_table();
        construct_incoming_table();
        return encoded;
    }

    void construct_routing_table() {
        Timer timer;
        cerr << "Started transforming routing_table." << endl;
        vector<char> routing_table_array;
        for(int node=1;node<=graph.num_nodes();node++) {
            routing_table_array.push_back('#');// to always start a block with #
            if (PathDatabaseDynamicCore::node_is_split(node)) {
                auto& dynamic_table = PathDatabaseDynamicCore::routing_table;
                for(int i=0;i<dynamic_table.size(node);i++) {
                    routing_table_array.push_back(dynamic_table.get(node,i));
                }
            }
        }
        routing_table_array.push_back('#'); // to also always end a block with #
        routing_table = decltype(routing_table)(graph,routing_table_array);
        routing_table.transformations = PathDatabaseDynamicCore::routing_table.transformations;
        statistics["transformation_routing_table_time"] = timer.elapsed();
        cerr << "Transformation finished in " << statistics["transformation_routing_table_time"] << endl;
    }

    void construct_incoming_table() {
        Timer timer;
        cerr << "Started transforming incoming_table." << endl;
        vector<int> incoming_table_builder;
        vector<bool> is_join_node;
        for(int node=1;node<=graph.num_nodes();node++) {
            is_join_node.push_back(1);
            if (PathDatabaseDynamicCore::node_is_join(node)) {
                auto new_reads = PathDatabaseDynamicCore::incoming_table.branch_size(node,'$');
                if (new_reads) {
                    is_join_node.push_back(0);
                    incoming_table_builder.push_back(new_reads);
                }
#ifdef ALL_EDGES_COVERED
                for(auto& base : "ACGTN") {
                    auto branch_size = PathDatabaseDynamicCore::incoming_table.branch_size(node,base);
                    if (branch_size) {
                        // so it is an actual edge in a graph (because all edges are covered)
                        is_join_node.push_back(0);
                        incoming_table_builder.push_back(branch_size);
                    }
                }
#else
                graph.call_incoming_kmers_mine(node,[&node,&incoming_table_builder,
                        &is_join_node,this](node_index prev_node,char c) {
                    auto branch_size = PathDatabaseDynamicCore::incoming_table.branch_size(node,c);
                    is_join_node.push_back(0);
                    incoming_table_builder.push_back(branch_size);
                });
#endif
        #ifndef FULL_INCOMING_TABLE
                assert(is_join_node.back() == 0);
                if (*(is_join_node.end()-2) == 0) {
                    is_join_node.pop_back();
                    incoming_table_builder.pop_back();
                }
        #endif
            }

        }
        is_join_node.push_back(1); // to also always end a block with 1
        incoming_table.edge_multiplicity_table = sdsl::enc_vector<>(incoming_table_builder);
        sdsl::bit_vector temporary_representation(is_join_node.size());
        for(int i=0;i<is_join_node.size();i++) {
            temporary_representation[i] = is_join_node[i];
        }
        incoming_table.joins = BitVector(temporary_representation);
        statistics["transformation_incoming_table_time"] = timer.elapsed();
        cerr << "Transformation finished in " << statistics["transformation_incoming_table_time"] << endl;
    }










    // is join or start of the read
    bool node_is_join(node_index node) const {
        return incoming_table.is_join(node);
    }

    // is split or end of the read
    bool node_is_split(node_index node) const {
        // is not the last element of routing table and next character is not starting of new node
        return routing_table.size(node);
    }



    void serialize(const fs::path& folder) const {
        fs::create_directories(folder / "xxx.bin");
        ofstream edge_multiplicity_file(folder / "edge_multiplicity.bin", ios_base::trunc | ios_base::out);
        ofstream routing_table_file(folder / "routing_table.bin", ios_base::trunc | ios_base::out);
        ofstream joins_file(folder / "joins.bin", ios_base::trunc | ios_base::out);
        string graph_filename = folder / "graph.bin";

        incoming_table.edge_multiplicity_table.serialize(edge_multiplicity_file);
        routing_table.serialize(routing_table_file);
        incoming_table.joins.serialize(joins_file);
        graph.serialize(graph_filename);
    }

    static PathDatabaseWaveletCore deserialize(const fs::path& folder) {
        ifstream edge_multiplicity_file(folder / "edge_multiplicity.bin");
        ifstream routing_table_file(folder / "routing_table.bin");
        ifstream joins_file(folder / "joins.bin");
        string graph_filename = folder / "graph.bin";

        auto graph = std::shared_ptr<DBGSuccinct>{
                new DBGSuccinct(21)
                };
        graph->load(graph_filename);
        auto db = PathDatabaseWaveletCore(graph);
        db.incoming_table.edge_multiplicity_table.load(edge_multiplicity_file);
        db.routing_table.load(routing_table_file);
        db.incoming_table.joins.load(joins_file);
        return db;
    }

    json get_statistics(unsigned int verbosity = ~0u) const {
        json result = PathDatabaseDynamicCore::get_statistics(verbosity);
        json routing_table_stats = routing_table.get_statistics(verbosity);
        result.update(statistics);
        result.update(routing_table_stats);
        int true_joins = 0;
        int added_joins = 0;
        int true_splits = 0;
        int added_splits = 0;
        std::map<int, int> joins_diff_symbols_histogram;
        std::map<int, int> splits_size_histogram;
        std::map<int, int> splits_diff_symbols_histogram;

        for (int node = 1; node <= graph.num_nodes();node++) {
            if (node_is_join(node)) {
                if (graph.indegree(node) > 1) {
                    true_joins++;
                }
                else {
                    added_joins++;
                }
                if (verbosity & STATS_JOINS_HISTOGRAM) {
                    int prev = 0;
                    int cardinality = incoming_table.size(node);
                    joins_diff_symbols_histogram[cardinality]++;
                }
            }
            if (node_is_split(node)) {
                if (graph.outdegree(node) > 1) {
                    true_splits++;
                }
                else {
                    added_splits++;
                }
                if (verbosity & STATS_SPLITS_HISTOGRAM) {
                    set<int> diff_symbols;
                    for (int i=0; i < routing_table.size(node); i++) {
                        diff_symbols.insert(routing_table.get(node,i));
                    }
                    splits_diff_symbols_histogram[diff_symbols.size()]++;
                    splits_size_histogram[routing_table.size(node)]++;
                }
            }
        }
        json addition = {{"true_joins", true_joins},
                       {"added_joins", added_joins},
                       {"true_splits", true_splits},
                       {"added_splits", added_splits},
                       {"num_of_nodes", graph.num_nodes()}
                      };
        result.update(addition);
        if (verbosity & STATS_SPLITS_HISTOGRAM) {
            result["splits_diff_symbols_histogram"] = splits_diff_symbols_histogram;
            result["splits_size_histogram"] = splits_size_histogram;
        }
        if (verbosity & STATS_JOINS_HISTOGRAM) {
            result["joins_diff_symbols_histogram"] = joins_diff_symbols_histogram;
        }
        return result;
    }



//protected:
    json statistics;
    RoutingTable<Wavelet> routing_table;
    IncomingTable<BitVector> incoming_table;


};
template<class Wavelet = sdsl::wt_rlmn<>,class BitVector=default_bit_vector>
class PathDatabaseWavelet : public QueryEnabler<DecodeEnabler<PathDatabaseWaveletCore<>>> {
    using QueryEnabler<DecodeEnabler<PathDatabaseWaveletCore<>>>::QueryEnabler;
};

#endif /* path_database_baseline_hpp */
