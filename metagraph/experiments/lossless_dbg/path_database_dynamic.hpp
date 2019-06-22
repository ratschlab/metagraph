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
#include <tsl/hopscotch_set.h>
#include <optional>
#include <omp.h>

#include "decode_enabler.hpp"
#include "exit_barrier.hpp"
#include "dynamic_routing_table.hpp"
#include "dynamic_incoming_table.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "unix_tools.hpp"
#include "threading.hpp"

// TODO: Never use 'using namespace std;' in .hpp files
// todo find a tool that removes this relative namespacing issue
using namespace std;

// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this
using node_index = DeBruijnGraph::node_index;

class VerboseTimer {
public:
    VerboseTimer(string procedure_name) : procedure_name(procedure_name) {
        cerr << "Started " << procedure_name << endl;
    }
    double finished() {
        double elapsed = timer.elapsed();
        cerr << "Finished " << procedure_name << " in " << elapsed << " sec." << endl;
        return elapsed;
    }
    string procedure_name;
    Timer timer;
};


using DefaultDynamicRoutingTable = DynamicRoutingTable<>;
using DefaultDynamicIncomingTable = DynamicIncomingTable<>;
CREATE_MEMBER_CHECK(transformations);


template<typename RoutingTableT=DefaultDynamicRoutingTable,typename IncomingTableT=DefaultDynamicIncomingTable>
class PathDatabaseDynamicCore {
public:
    using GraphT = DBGSuccinct;
    using path_id = pair<node_index,int64_t>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseDynamicCore(std::shared_ptr<const GraphT> graph, int64_t chunk_size=DefaultChunks) :
            graph_(graph),
            graph(*graph_),
            incoming_table(graph_),
            routing_table(graph_),// weak
            chunk_size(chunk_size)
            {}


    explicit PathDatabaseDynamicCore(const vector<string> &reads,
                 size_t kmer_length = 21 /* default kmer */) :
            graph_(shared_ptr<const GraphT>(buildGraph(this,reads,kmer_length))),
            graph(*graph_),
            incoming_table(graph_),
            routing_table(graph_),
            chunk_size(DefaultChunks)
            {}

    virtual ~PathDatabaseDynamicCore() {}



    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        VerboseTimer encoding_timer("encoding reads");
        VerboseTimer preprocessing_timer("preprocessing step");
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering

        // use std::conditional<has transformations,vector<string>&,vector<string>> transformed_sequences = transformed_sequences;
        // problem is that transform_sequence will additionally copy & replace it so we are getting additional overhead
//        vector<string> transformed_sequences = sequences;

        //decide if the transformed sequences will be different or not and either copy or point to original sequences.
        typename std::conditional<has_member_transformations<RoutingTableT>::value,
                         vector<string>,
                         const vector<string>&>::type transformed_sequences = sequences;
        if constexpr (has_member_transformations<RoutingTableT>::value) {
            #pragma omp parallel for
            for (ll i = 0; i < sequences.size(); i++) {
                transform_sequence_inplace(transformed_sequences[i]);
            }
        }

        populate_additional_joins(transformed_sequences);

// Mark split and join-nodes in graph for faster queries and construction
        populate_bifurcation_bitvectors();

        auto alloc_routing_table = VerboseTimer("allocation of routing & incoming table");
        routing_table = decltype(routing_table)(graph_,&is_bifurcation,&rank_is_bifurcation,chunk_size);
        incoming_table = decltype(incoming_table)(graph_,&is_bifurcation,&rank_is_bifurcation,chunk_size);
        alloc_routing_table.finished();

        #pragma omp parallel for num_threads(get_num_threads())
        for (node_index node = 1; node <= graph.num_nodes(); node++) {
            assert(is_split[node] == node_is_split_raw(node));
            assert(is_join[node] == node_is_join_raw(node));
        }

        vector<path_id> encoded(sequences.size());

#ifndef DISABLE_PARALELIZATION

        auto alloc_lock_t = VerboseTimer("memory allocation of locks");
        ChunkedDenseHashMap<omp_lock_t,decltype(is_bifurcation), decltype(rank_is_bifurcation),false> node_locks(&is_bifurcation,&rank_is_bifurcation,chunk_size);
        alloc_lock_t.finished();

        //using Barrier = ReferenceExitBarrier<>;
        //using Barrier = IdentityComparator<ExitBarrier<>,ReferenceExitBarrier<>>;




        auto lock_init_timer = VerboseTimer("initializing locks & barriers");
        for(int64_t i=0;i<node_locks.elements.size();i++) {
            omp_init_lock(&node_locks.elements[i]);
        }
        lock_init_timer.finished();

#endif
        auto threads_map_t = VerboseTimer("memory allocation of exit barriers");
        using Barrier = ExitBarrier<>;
        auto exit_barrier = Barrier(&is_bifurcation,&rank_is_bifurcation,chunk_size);
        threads_map_t.finished();


        statistics["preprocessing_time"] = preprocessing_timer.finished();

        cerr << "Cumulative sizes in bytes" << endl;
        PRINT_VAR(is_join.size()/8);
        PRINT_VAR(is_split.size()/8);
        VerboseTimer routing_timer("routing step");
#ifndef DISABLE_PARALELIZATION
        #pragma omp parallel for num_threads(get_num_threads()) //default(none) shared(encoded,node_locks,routing_table,incoming_table)
#endif
        for (size_t i = 0; i < sequences.size(); i++) {
            int16_t tid = omp_get_thread_num();
            std::vector<std::tuple<node_index, char, char>> bifurcations;

            const auto &sequence = sequences[i];
            const auto& path_for_sequence = transformed_sequences[i];

            size_t kmer_begin = 0;
            auto kmer_end = graph.get_k();
#ifdef DEBUG_ADDITIONAL_INFORMATION
            vector<node_index> debug_list_of_nodes;
            vector<int64_t> debug_relative_position_history;
            bool debug_transformation_used = sequence != path_for_sequence;
#endif
            // TODO: use map_to_nodes_sequentially when implemented
            graph.map_to_nodes(path_for_sequence, [&](node_index node) {
                assert(node || [&](){
                    auto pn = [&](int64_t offset) {
                        PRINT_VAR(sequence.substr(offset,graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(sequence.substr(offset,graph.get_k())));
                        PRINT_VAR(path_for_sequence.substr(offset,graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(path_for_sequence.substr(offset,graph.get_k())));
                    };
                    pn(kmer_begin);
                    return false; }());

#ifdef DEBUG_ADDITIONAL_INFORMATION
                debug_list_of_nodes.push_back(node);
#endif
                char join_char = node_is_join(node)
                                 ? (kmer_begin ? path_for_sequence[kmer_begin-1] : '$')
                                 : '\0';
                char split_char = node_is_split(node)
                                  ? (kmer_end < sequence.size() ?
                                     sequence[kmer_end] :
                                     '$')
                                  : '\0';

                if (join_char || split_char) {
                    bifurcations.emplace_back(node, join_char, split_char);
                }
                kmer_begin++;
                kmer_end++;
            });

            auto &[first_node,first_join,first_split] = bifurcations.front();
            auto &[last_node,last_join,last_split] = bifurcations.back();
            assert(first_join == '$');
            assert(last_split == '$');

            int64_t relative_position = INT_MIN;
#ifdef DEBUG_ADDITIONAL_INFORMATION
            int64_t debug_bifurcation_idx = 0;
#endif
            int64_t previous_node = 0;
            char traversed_edge = '\0';
            for (const auto &[node, join_symbol, split_symbol] : bifurcations) {
#ifndef DISABLE_PARALELIZATION
                omp_lock_t *node_lock;
                node_lock = node_locks.ptr_to(node);
                omp_set_lock(node_lock);
#endif
                if (previous_node) {
                    exit_barrier.exit(previous_node,tid);
                }

                if (join_symbol) {
                    if (join_symbol == '$') {
                        // always putting new read above all other reads
                        relative_position = incoming_table.branch_size(node, '$');
                        encoded[i] = {node, relative_position};
                    }
                    assert(relative_position>=0);
                    assert((relative_position <= incoming_table.branch_size(node,join_symbol) || [&,node=node,join_symbol=join_symbol]{
                        using TT = DecodeEnabler<PathDatabaseDynamicCore<>>;
                        auto self = reinterpret_cast<TT*>(this);
                        cerr << "current" << endl;
                        PRINT_VAR(graph.get_node_sequence(node));
                        PRINT_VAR(node,tid,relative_position,previous_node,traversed_edge);
                        return false;
                    }()));

#ifdef DEBUG_ADDITIONAL_INFORMATION
                    debug_relative_position_history.push_back(relative_position);
#endif
                    auto lower_routes = incoming_table.branch_offset_and_increment(node, join_symbol);
                    assert(lower_routes>=0);
                    relative_position += lower_routes;
                }

                if (split_symbol) {
                    assert(relative_position>=0);
                    routing_table.insert(node, relative_position, split_symbol);
#ifdef DEBUG_ADDITIONAL_INFORMATION
                    debug_relative_position_history.push_back(relative_position);
#endif
                    relative_position = routing_table.new_relative_position(node, relative_position);
                }
                traversed_edge = split_symbol ? split_symbol : 'X'; // X as it doesn't matter
                if (traversed_edge != '$') { // doesn't need any update
                    exit_barrier.enter(node,traversed_edge,relative_position,tid);
                }
#ifndef DISABLE_PARALELIZATION
                omp_unset_lock(node_lock);
#endif
#ifdef DEBUG_ADDITIONAL_INFORMATION
                debug_bifurcation_idx++;
#endif
                previous_node = node;
            }
        }

#ifndef DISABLE_PARALELIZATION
        for(int64_t i=0;i<node_locks.elements.size();i++) {
            omp_destroy_lock(&node_locks.elements[i]);

        }
#endif
        encoded_paths += encoded.size();
        statistics["routing_time"] = routing_timer.finished();
        encoding_timer.finished();

        return encoded;
    }

#ifndef DISABLE_PARALELIZATION


#endif

    void populate_additional_joins(const vector<std::string> &transformed_sequences) {
        auto additional_splits_t = VerboseTimer("computing additional splits and joins");
        vector<node_index> additional_joins_vec(transformed_sequences.size());
        vector<node_index> additional_splits_vec(transformed_sequences.size());
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < transformed_sequences.size(); ++i) {
            auto& transformed_sequence = transformed_sequences[i];
            // add additional bifurcation
            additional_joins_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(0, graph.get_k())
            );
            assert(additional_joins_vec[i]); // node has to be in graph
            additional_splits_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(transformed_sequence.length() - graph.get_k())
            );
            assert(additional_splits_vec[i]); // node has to be in graph
        }
        additional_joins = decltype(additional_joins)(additional_joins_vec.begin(),
                                                      additional_joins_vec.end());
        additional_joins_vec.clear();
        additional_splits = decltype(additional_splits)(additional_splits_vec.begin(),
                                                        additional_splits_vec.end());
        additional_splits_vec.clear();
        additional_splits_t.finished();
    }

    void populate_bifurcation_bitvectors() {
        auto bifurcation_timer = VerboseTimer("construction of bifurcation bit_vectors");
        is_split = decltype(is_split)(graph.num_nodes() + 1); // bit
        is_join = decltype(is_join)(graph.num_nodes() + 1);
        is_bifurcation = decltype(is_join)(graph.num_nodes() + 1);

        #pragma omp parallel for num_threads(get_num_threads())
        for (uint64_t node = 0; node <= graph.num_nodes(); node += 8) {
            for (int64_t i = node; i < node + 8 && i <= graph.num_nodes(); ++i) {
                if (!i)
                    continue;
                is_split[i] = node_is_split_raw(i);
                is_join[i] = node_is_join_raw(i);
                is_bifurcation[i] = is_split[i] || is_join[i];
            }
        }

        rank_is_split = decltype(rank_is_split)(&is_split);
        rank_is_join = decltype(rank_is_join)(&is_join);
        rank_is_bifurcation = decltype(rank_is_bifurcation)(&is_bifurcation);

        bifurcation_timer.finished();
    }


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

    void compress() {

    }

    size_t num_paths() const {
        return encoded_paths;
    };

    json get_statistics(uint64_t verbosity = ~0u) const {
        //auto result = PathDatabase<pair<node_index,int>,GraphT>::get_statistics(verbosity);
        //result.update(statistics);
        //return result;
        return statistics;
    }


    //todo put to proper place
    string transform_sequence(const string& sequence) const {
        string result = sequence;
        size_t kmer_end = graph.get_k();
        static_assert(has_member_transformations<RoutingTableT>::value == std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,RoutingTableT>::value);
        if constexpr (std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,RoutingTableT>::value) {
            graph.map_to_nodes(sequence, [&](node_index node) {
                if (kmer_end < sequence.size()) {
                    result[kmer_end] = routing_table.transform(node, result[kmer_end]);
                }
                kmer_end++;
            });
        }
        return result;
    }

    void transform_sequence_inplace(string& sequence) const {
        size_t kmer_end = graph.get_k();
        static_assert(has_member_transformations<RoutingTableT>::value == std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,RoutingTableT>::value);
        if constexpr (std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,RoutingTableT>::value) {
            graph.map_to_nodes(sequence, [&](node_index node) {
                if (kmer_end < sequence.size()) {
                    sequence[kmer_end] = routing_table.transform(node, sequence[kmer_end]);
                }
                kmer_end++;
            });
        }
    }

    void serialize(const fs::path& folder) const {};

//protected:
    json statistics;

    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int64_t encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)


    tsl::hopscotch_set<node_index> additional_joins;
    tsl::hopscotch_set<node_index> additional_splits;

    sdsl::bit_vector is_join;
    sdsl::bit_vector is_split;
    sdsl::bit_vector is_bifurcation;
    typename decltype(is_join)::rank_1_type rank_is_join;
    typename decltype(is_split)::rank_1_type rank_is_split;
    typename decltype(is_bifurcation)::rank_1_type rank_is_bifurcation;

    std::shared_ptr<const GraphT> graph_;
    const GraphT& graph;
    RoutingTableT routing_table;
    IncomingTableT incoming_table;
    int64_t chunk_size;

    static DBGSuccinct* buildGraph(PathDatabaseDynamicCore* self,vector<string> reads,int64_t kmer_length) {
        Timer timer;
        cerr << "Started building the graph" << endl;
        auto graph = new DBGSuccinct(dbg_succ_graph_constructor(reads, kmer_length));
        graph->mask_dummy_kmers(1, false);
        auto elapsed = timer.elapsed();
        cerr << "Building finished in " << elapsed << " sec." << endl;
        self->statistics["graph_build_time"] = elapsed;
        return graph;
    }

};

template<typename RoutingTableT=DefaultDynamicRoutingTable,typename IncomingTableT=DefaultDynamicIncomingTable>
class PathDatabaseDynamic : public DecodeEnabler<PathDatabaseDynamicCore<RoutingTableT,IncomingTableT>> {
    using DecodeEnabler<PathDatabaseDynamicCore<RoutingTableT,IncomingTableT>>::DecodeEnabler;
};




#endif /* path_database_baseline_hpp */
