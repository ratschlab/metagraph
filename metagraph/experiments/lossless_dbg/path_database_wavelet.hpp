#pragma ide diagnostic ignored "openmp-use-default-none"
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


using namespace std;

const uint64_t STATS_JOINS_HISTOGRAM (1u << 0u);
const uint64_t STATS_SPLITS_HISTOGRAM (1u << 1u);

using alphabets::log2;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this

template <typename StaticTable>
class DynamicVersion {

};
template <>
class DynamicVersion<RoutingTable<>> {
public:
    using Type = DynamicRoutingTable<>;
};
template <>
class DynamicVersion<RoutingTableCore<>> {
public:
    using Type = DynamicRoutingTableCore<>;
};
template <>
class DynamicVersion<IncomingTable<>> {
public:
    using Type = DynamicIncomingTable<>;
};

using DefaultRoutingTable = RoutingTable<>;
using DefaultIncomingTable = IncomingTable<>;
//template <class Wavelet = sdsl::wt_rlmn<sdsl::sd_vector<>>
template<class RoutingTableT = DefaultRoutingTable,class IncomingTableT=DefaultIncomingTable>
class PathDatabaseWaveletCore : public PathDatabaseDynamicCore<typename DynamicVersion<RoutingTableT>::Type,typename DynamicVersion<IncomingTableT>::Type> {
public:
    using DRT = typename DynamicVersion<RoutingTableT>::Type;
    using DIT = typename DynamicVersion<IncomingTableT>::Type;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseWaveletCore(std::shared_ptr<const DBGSuccinct> graph,int chunks=DefaultChunks) : PathDatabaseDynamicCore<DRT,DIT>(graph,chunks)
                                                                              {}

    PathDatabaseWaveletCore(const vector<string> &filenames,
                        size_t kmer_length = 21 /* default */) : PathDatabaseDynamicCore<DRT,DIT>(filenames,kmer_length)
                                                         {}



    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        auto encoded = PathDatabaseDynamicCore<DRT,DIT>::encode(sequences);

        // convert dynamic_(routing_table/incoming_table) to routing_table/incoming_table
        construct_routing_table();
        construct_incoming_table();
        return encoded;
    }

    void construct_routing_table() {
        Timer timer;
        cerr << "Started transforming routing_table." << endl;
        vector<char> routing_table_array;
        #pragma omp parallel for reduction(append: routing_table_array)
        for(int64_t node=1;node<=this->graph.num_nodes();node++) {
            routing_table_array.push_back('#');// to always start a block with #
            if (PathDatabaseDynamicCore<DRT,DIT>::node_is_split(node)) {
                auto& dynamic_table = PathDatabaseDynamicCore<DRT,DIT>::routing_table;
                for(int64_t i=0;i<dynamic_table.size(node);i++) {
                    routing_table_array.push_back(dynamic_table.get(node,i));
                }
            }
        }
        routing_table_array.push_back('#'); // to also always end a block with #
        routing_table.initialize_content(routing_table_array);
        if constexpr (std::is_base_of<TransformationsEnabler<RoutingTableCore<>>,RoutingTableT>::value) {
            routing_table.transformations = PathDatabaseDynamicCore< DRT, DIT>::routing_table.transformations;
        }
        statistics["transformation_routing_table_time"] = timer.elapsed();
        cerr << "Transformation finished in " << statistics["transformation_routing_table_time"] << endl;
    }

    void construct_incoming_table() {
        Timer timer;
        cerr << "Started transforming incoming_table." << endl;

        vector<int64_t> incoming_table_builder;
        vector<bool> delimiter_vector;


        #pragma omp parallel for reduction(append: incoming_table_builder, delimiter_vector)
        for(int64_t node=1;node<=this->graph.num_nodes();node++) {
            delimiter_vector.push_back(true);
            if (PathDatabaseDynamicCore<DRT,DIT>::node_is_join(node)) {
                auto new_reads = PathDatabaseDynamicCore<DRT,DIT>::incoming_table.branch_size(node,'$');
                int current_table_size = 0;
                if (new_reads) {
                    incoming_table_builder.push_back(new_reads);
                    current_table_size++;
                    delimiter_vector.push_back(false);
                }
#ifdef ALL_EDGES_COVERED
                for(auto& base : "ACGTN") {
                    auto branch_size = PathDatabaseDynamicCore<DRT,DIT>::incoming_table.branch_size(node,base);
                    if (branch_size) {
                        // so it is an actual edge in a this->graph (because all edges are covered)
                        incoming_table_builder.push_back(branch_size);
                        current_table_size++;
                        delimiter_vector.push_back(false);
                    }
                }
#else
                this->graph.call_incoming_kmers_mine(node,[&,this](node_index prev_node,char c) {
                    auto branch_size = PathDatabaseDynamicCore<DRT,DIT>::incoming_table.branch_size(node,c);
                    incoming_table_builder.push_back(branch_size);
                    current_table_size++;
                    delimiter_vector.push_back(false);
                });
#endif
                assert((current_table_size>0 || [&](){
                    PathDatabaseDynamicCore<DRT,DIT>::incoming_table.print_content(node);
                    PRINT_VAR(PathDatabaseDynamicCore<DRT,DIT>::node_is_join(node));
#ifdef MEMOIZE
                    PRINT_VAR(PathDatabaseDynamicCore<DRT,DIT>::node_is_join_raw(node));
#endif
                    return false;
                }()));


#ifndef FULL_INCOMING_TABLE
                if (current_table_size > 1) {
                    incoming_table_builder.pop_back();
                }
        #endif
            }
        }
        delimiter_vector.push_back(true);
        sdsl::bit_vector temporary_representation(delimiter_vector.size());
        for(int64_t i=0; i <delimiter_vector.size();i++) {
            temporary_representation[i] = delimiter_vector[i];
        }
        using joins_type = typename decltype(incoming_table)::BitVector;
        incoming_table = decltype(incoming_table)(this->graph_,joins_type(temporary_representation),incoming_table_builder);
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
        Timer timer;
        cerr << "Started serializing the path encoder." << endl;
        fs::create_directories(folder / "path_encoder.flag");
        ofstream edge_multiplicity_file(folder / "edge_multiplicity.bin", ios_base::trunc | ios_base::out);
        ofstream routing_table_file(folder / "routing_table.bin", ios_base::trunc | ios_base::out);
        ofstream joins_file(folder / "joins.bin", ios_base::trunc | ios_base::out);
        string graph_filename = folder / "graph.bin";

        incoming_table.edge_multiplicity_table.serialize(edge_multiplicity_file);
        routing_table.serialize(routing_table_file);
        incoming_table.joins.serialize(joins_file);
        this->graph.serialize(graph_filename);
        cerr << "Finished serializing the path encoder in " << timer.elapsed() << " sec." << endl;
    }

    static auto deserialize(const fs::path& folder) {
        ifstream edge_multiplicity_file(folder / "edge_multiplicity.bin");
        ifstream routing_table_file(folder / "routing_table.bin");
        ifstream joins_file(folder / "joins.bin");
        string graph_filename = folder / "graph.bin";

        auto graph = std::shared_ptr<DBGSuccinct>{
                new DBGSuccinct(21)
                };
        graph->load(graph_filename);
        auto db = QueryEnabler<DecodeEnabler<PathDatabaseWaveletCore<RoutingTableT,IncomingTableT>>>(graph);
        db.incoming_table.edge_multiplicity_table.load(edge_multiplicity_file);
        db.routing_table.load(routing_table_file);
        db.incoming_table.joins.load(joins_file);
        return db;
    }

    json get_statistics(uint64_t verbosity = ~0u) const {
        json result = PathDatabaseDynamicCore<DRT,DIT>::get_statistics(verbosity);
        json routing_table_stats = routing_table.get_statistics(verbosity);
        result.update(statistics);
        result.update(routing_table_stats);
        int64_t added_joins = 0;
        int64_t added_splits = 0;
        std::map<int64_t, int64_t> joins_diff_symbols_histogram;
        std::map<int64_t, int64_t> splits_size_histogram;
        std::map<int64_t, int64_t> splits_diff_symbols_histogram;

        for (int64_t node = 1; node <= this->graph.num_nodes();node++) {
            if (node_is_join(node)) {
                if (this->graph.indegree(node) <= 1) {
                    added_joins++;
                }
                if (verbosity & STATS_JOINS_HISTOGRAM) {
                    int64_t cardinality = incoming_table.size(node);
                    joins_diff_symbols_histogram[cardinality]++;
                }
            }
            if (node_is_split(node)) {
                if (this->graph.outdegree(node) <= 1) {
                    added_splits++;
                }
                if (verbosity & STATS_SPLITS_HISTOGRAM) {
                    //set<int64_t> diff_symbols;
                    //for (int64_t i=0; i < routing_table.size(node); i++) {
                    //    diff_symbols.insert(routing_table.get(node,i));
                    //}
                    //splits_diff_symbols_histogram[diff_symbols.size()]++;
                    splits_size_histogram[routing_table.size(node)]++;
                }
            }
        }
        json addition = {
                       {"added_joins", added_joins},
                       {"added_splits", added_splits},
                       {"num_of_nodes", this->graph.num_nodes()}
                      };
        result.update(addition);
        if (verbosity & STATS_SPLITS_HISTOGRAM) {
            //result["splits_diff_symbols_histogram"] = splits_diff_symbols_histogram;
            result["splits_size_histogram"] = splits_size_histogram;
        }
        if (verbosity & STATS_JOINS_HISTOGRAM) {
            result["joins_diff_symbols_histogram"] = joins_diff_symbols_histogram;
        }
        return result;
    }



//protected:
    mutable json statistics;
    RoutingTableT routing_table;
    IncomingTableT incoming_table;


};

template<class RoutingTableT = DefaultRoutingTable,class IncomingTableT=DefaultIncomingTable>
using PathDatabaseWavelet = QueryEnabler<DecodeEnabler<PathDatabaseWaveletCore<RoutingTableT,IncomingTableT>>>;

template<typename DummyT=int64_t>
using PathDatabaseWaveletWithtoutTransformation = QueryEnabler<DecodeEnabler<PathDatabaseWaveletCore<RoutingTableCore<>,DefaultIncomingTable>>>;

#endif /* path_database_baseline_hpp */

