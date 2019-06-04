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
#include "dynamic_routing_table.hpp"
#include "dynamic_incoming_table.hpp"
#include "utils.hpp"
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
    void finished() {
        cerr << "Finished " << procedure_name << " in " << timer.elapsed() << " sec." << endl;
    }
    string procedure_name;
    Timer timer;
};

struct waiting_thread_info_t {
    int16_t thread_id;
    int64_t relative_position;
    uint64_t node;
    char traversed_edge;
    friend ostream& operator<<(ostream& os, const waiting_thread_info_t& dt);
};
ostream& operator<<(ostream& os, const waiting_thread_info_t& wt)
{
    os << "[tid=" << wt.thread_id << ",rp=" << wt.relative_position << ",node=" << wt.node << ",tedg=" << wt.traversed_edge << "]";
    return os;
}

template<typename GraphT=DBGSuccinct,typename RoutingTableT=DynamicRoutingTable<>,typename IncomingTableT=DynamicIncomingTable<>>
class PathDatabaseDynamicCore {
public:

    using path_id = pair<node_index,int64_t>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseDynamicCore(std::shared_ptr<const GraphT> graph, int64_t chunk_size=1000) :
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
            chunk_size(1000)
            {}

    virtual ~PathDatabaseDynamicCore() {}



    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        // memory expected size
        PRINT_VAR(sizeof(vector<char>));
        PRINT_VAR(sizeof(omp_lock_t));
        PRINT_VAR(sizeof(deque<tuple<int64_t,int64_t,int64_t>>));

        Timer timer;
        cerr << "Started encoding reads" << endl;
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering

        populate_additional_joins(sequences);

// Mark split and join-nodes in graph for faster queries and construction
        populate_bifurcation_bitvectors();

        auto alloc_routing_table = VerboseTimer("allocation of routing & incoming table");
        routing_table = decltype(routing_table)(graph_,&is_split,&rank_is_split,chunk_size);
        incoming_table = decltype(incoming_table)(graph_,&is_join,&rank_is_join,chunk_size);
        alloc_routing_table.finished();

        #pragma omp parallel for num_threads(get_num_threads())
        for (node_index node = 1; node <= graph.num_nodes(); node++) {

            assert(is_split[node] == node_is_split_raw(node));
            assert(is_join[node] == node_is_join_raw(node));
        }

        vector<path_id> encoded(sequences.size());

#ifdef _OPENMP
        auto alloc_lock_t = VerboseTimer("memory allocation of locks");

        ChunkedDenseHashMap<omp_lock_t,decltype(is_bifurcation), decltype(rank_is_bifurcation),false> node_locks(&is_bifurcation,&rank_is_bifurcation,chunk_size);
        ChunkedDenseHashMap<omp_lock_t,decltype(is_bifurcation), decltype(rank_is_bifurcation),false> outgoing_locks(&is_bifurcation,&rank_is_bifurcation,chunk_size);
        alloc_lock_t.finished();
        auto threads_map_t = VerboseTimer("memory allocation of thread deques");
        ChunkedDenseHashMap<deque<waiting_thread_info_t>,decltype(is_bifurcation), decltype(rank_is_bifurcation),false> waiting_threads(&is_bifurcation,&rank_is_bifurcation,chunk_size);
        threads_map_t.finished();
        auto lock_init_timer = VerboseTimer("initializing locks");
        for(int64_t i=0;i<node_locks.elements.size();i++) {
            omp_init_lock(&node_locks.elements[i]);
            omp_init_lock(&outgoing_locks.elements[i]);
        }
        lock_init_timer.finished();

#endif

        statistics["preprocessing_time"] = timer.elapsed();
        std::cerr << "Finished preprocessing in " << statistics["preprocessing_time"] << " sec" << std::endl;

        timer.reset();
        cerr << "Cumulative sizes in bytes" << endl;
        PRINT_VAR(is_join.size()/8);
        PRINT_VAR(is_split.size()/8);

        #pragma omp parallel for num_threads(get_num_threads()) //default(none) shared(encoded,node_locks,routing_table,incoming_table)
        for (size_t i = 0; i < sequences.size(); i++) {
            int16_t tid = omp_get_thread_num() + 1;
            std::vector<std::tuple<node_index, char, char>> bifurcations;

            const auto &sequence = sequences[i];
            auto path_for_sequence = transform_sequence(sequence);

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
                omp_lock_t *node_lock;
                node_lock = node_locks.ptr_to(node);
                omp_set_lock(node_lock);
                if (previous_node) {
                    auto outgoing_lock = outgoing_locks.ptr_to(previous_node);
                    omp_set_lock(outgoing_lock);
                    auto& waiting_queue = waiting_threads[previous_node];
                    relative_position = update_position_after_wakeup(tid, previous_node, traversed_edge,
                                                                     waiting_queue, relative_position);
                    omp_unset_lock(outgoing_lock);
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
#ifdef DEBUG_ADDITIONAL_INFORMATION
                        PRINT_VAR(debug_bifurcation_idx);
#endif
                        incoming_table.print_content(node);
                        PRINT_VAR(join_symbol);
                        cerr << "previous" << endl;
#ifdef DEBUG_ADDITIONAL_INFORMATION
                        auto& prev_bifurcation = bifurcations[debug_bifurcation_idx-1];
                        const auto &[previous_node_2, prev_join_symbol, prev_split_symbol] = prev_bifurcation;
                        PRINT_VAR(graph.get_node_sequence(previous_node));
                        if (prev_split_symbol) {
                            routing_table.print_content(previous_node);
                            PRINT_VAR(debug_relative_position_history.back());
                            PRINT_VAR(prev_split_symbol);
                        }
                        if (prev_join_symbol) {
                            incoming_table.print_content(previous_node);
                            int64_t prev_position = !prev_split_symbol ?
                                                    debug_relative_position_history.back()
                                                                       :
                                                    *(debug_relative_position_history.end()-2)
                            ;
                            auto path_id = self->get_global_path_id(previous_node,prev_position-1);
                            PRINT_VAR(transform_sequence(decode_from_input(sequences,
                                                                           {graph.get_node_sequence(path_id.first),path_id.second}
                                    ,graph.get_k())
                            ));
                            PRINT_VAR(prev_position);
                            PRINT_VAR(prev_join_symbol);
                        }
#endif
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
#ifdef _OPENMP
                traversed_edge = split_symbol ? split_symbol : 'X'; // X as it doesn't matter
                if (traversed_edge != '$') { // doesn't need any update
                    auto outgoing_lock = outgoing_locks.ptr_to(node);
                    omp_set_lock(outgoing_lock);
                    waiting_threads[node].push_back({tid, relative_position, node, traversed_edge});
                    omp_unset_lock(outgoing_lock);
                }
                omp_unset_lock(node_lock);
#endif
#ifdef DEBUG_ADDITIONAL_INFORMATION
                debug_bifurcation_idx++;
#endif
                previous_node = node;
            }
        }

#ifdef _OPENMP
        for(int64_t i=0;i<node_locks.elements.size();i++) {
            omp_destroy_lock(&node_locks.elements[i]);
            omp_destroy_lock(&outgoing_locks.elements[i]);
        }
#endif
        encoded_paths += encoded.size();
        statistics["routing_time"] = timer.elapsed();
        std::cerr << "Finished routing in " << statistics["routing_time"] << " sec" << std::endl;

        return encoded;
    }

#ifdef _OPENMP

    int64_t update_position_after_wakeup(int tid,
                                 int64_t previous_node, char traversed_edge,
                                 deque<waiting_thread_info_t> &waiting_queue,
                                 int64_t relative_position) const {
            assert(previous_node);
            bool first_it = 1;
            bool me_first = 0;
            bool was_me = 0;
            int64_t past_offset = 0;
            // todo_wrap in define or better in macro
            int64_t debug_my_id = 0;
            int64_t debug_idx = 0;
            for(auto&other : waiting_queue) {
                if (was_me) {
                    // future
                    if (other.node == previous_node and
                            other.traversed_edge == traversed_edge and
                            other.thread_id < 0 and // job has already finished
                            other.relative_position <= relative_position) {
#ifdef DEBUG_ORDER_CORRECTION
                        cerr << waiting_queue << endl;
                        PRINT_VAR(tid, node, other.relative_position);
#endif
                        relative_position++;
                    }
                }
                else if (tid == other.thread_id) {
                    debug_my_id = debug_idx;
                    assert(!was_me); // can appear only once
                    // present
                    was_me = 1;
                    if (first_it) {
                        me_first = 1;
                    }
                    other.thread_id = -tid;
                }  else {
                    // past
                    if (other.node == previous_node and
                            other.traversed_edge == traversed_edge and
                            other.thread_id > 0 and// job is waiting
                            other.relative_position < relative_position) {
                            past_offset--;
                    }
                }
                first_it = 0;
                debug_idx++;
            }
            assert(was_me);
            relative_position += past_offset;
            assert(relative_position >= 0 || [&](){
                cerr << waiting_queue << endl;
                PRINT_VAR(tid,traversed_edge,debug_my_id,past_offset,relative_position);
                return false;
            }());

            if (me_first) {
                while (!waiting_queue.empty() and waiting_queue.front().thread_id < 0) {
                    waiting_queue.pop_front();
                }
            }
            return relative_position;
        }
#endif

    void populate_additional_joins(const vector<std::string> &sequences) {
        auto additional_splits_t = VerboseTimer("computing additional splits and joins");
        vector<node_index> additional_joins_vec(sequences.size());
        vector<node_index> additional_splits_vec(sequences.size());
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < sequences.size(); ++i) {
            const auto &sequence = sequences[i];
            auto transformed_sequence = transform_sequence(sequence);
            // add additional bifurcation
            additional_joins_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(0, graph.get_k())
            );
            assert(additional_joins_vec[i]); // node has to be in graph
            additional_splits_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(sequence.length() - graph.get_k())
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

    static BOSS* dbg_succ_graph_constructor(const vector<string> &reads,
                                            size_t kmer_length) {

        auto graph_constructor = BOSSConstructor(kmer_length - 1);// because BOSS has smaller kmers

        for (const auto &read : reads) {
            graph_constructor.add_sequence(read);
        }

        return new BOSS(&graph_constructor);
    }

};

template<typename GraphT=DBGSuccinct>
class PathDatabaseDynamic : public DecodeEnabler<PathDatabaseDynamicCore<GraphT>> {
    using DecodeEnabler<PathDatabaseDynamicCore<GraphT>>::DecodeEnabler;
};




#endif /* path_database_baseline_hpp */
