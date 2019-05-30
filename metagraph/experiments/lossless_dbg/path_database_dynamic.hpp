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



template<typename GraphT=DBGSuccinct>
class PathDatabaseDynamicCore {
public:
    using path_id = pair<node_index,int>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseDynamicCore(std::shared_ptr<const GraphT> graph) :
            graph(*graph),
            incoming_table(*graph)
            {}


    explicit PathDatabaseDynamicCore(const vector<string> &reads,
                 size_t kmer_length = 21 /* default kmer */) :
            graph(*buildGraph(this,reads,kmer_length)),
            incoming_table(graph),
            routing_table(graph)
            {}

    virtual ~PathDatabaseDynamicCore() {}

    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        Timer timer;
        cerr << "Started encoding reads" << endl;
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering
        // #pragma omp parallel for num_threads(num_threads)

        std::vector<node_index> additional_joins_vec(sequences.size());
        std::vector<node_index> additional_splits_vec(sequences.size());

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < sequences.size(); ++i) {
            const auto &sequence = sequences[i];
            //todo transform sequence here
            auto transformed_sequence = transform_sequence(sequence);
            // add additional bifurcation
            additional_joins_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(0, graph.get_k())
            );
            additional_splits_vec[i] = graph.kmer_to_node(
                    transformed_sequence.substr(sequence.length() - graph.get_k())
            );
        }

        additional_joins = decltype(additional_joins)(additional_joins_vec.begin(),
                                                      additional_joins_vec.end());
        additional_joins_vec = std::vector<node_index>();

        additional_splits = decltype(additional_splits)(additional_splits_vec.begin(),
                                                        additional_splits_vec.end());
        additional_splits_vec = std::vector<node_index>();

// Mark split and join-nodes in graph for faster queries and construction
#ifdef MEMOIZE
        is_split = decltype(is_split)(graph.num_nodes() + 1); // bit
        is_join = decltype(is_join)(graph.num_nodes() + 1);

        #pragma omp parallel for num_threads(get_num_threads())
        for (uint64_t node = 0; node <= graph.num_nodes(); node += 8) {
            for (int i = node; i < node + 8 && i <= graph.num_nodes(); ++i) {
                if (!i)
                    continue;
                is_split[i] = node_is_split_raw(i);
                is_join[i] = node_is_join_raw(i);
            }
        }
#endif
        //prepare hash maps
        for (uint64_t node = 0; node <= graph.num_nodes(); node += 8) {
            for (int i = node; i < node + 8 && i <= graph.num_nodes(); ++i) {
                if (!i)
                    continue;

                if (node_is_split(i)) {
                    routing_table.routing_table[i];
                }
                if (node_is_join(i)) {
                    incoming_table.incoming_table[i];
                }
            }
        }

#ifdef MEMOIZE
        #pragma omp parallel for num_threads(get_num_threads())
        for (node_index node = 1; node <= graph.num_nodes(); node++) {
            assert(is_split[node] == node_is_split_raw(node));
            assert(is_join[node] == node_is_join_raw(node));
        }
#endif

        vector<path_id> encoded(sequences.size());

        statistics["preprocessing_time"] = timer.elapsed();
        std::cerr << "Finished preprocessing in " << statistics["preprocessing_time"] << " sec" << std::endl;

        //ProgressBar progress_bar(sequences.size(), "Building dRT and dEM");

        const size_t batch_size = std::max(static_cast<size_t>(1),
                                           static_cast<size_t>(sequences.size() / 200));

        double join_time = 0;
        double split_time = 0;

#ifdef _OPENMP
        vector<omp_lock_t> node_locks(graph.num_nodes()+1);
        vector<omp_lock_t> outgoing_locks(graph.num_nodes()+1);
        vector<deque<tuple<int,int,int>>> waiting_threads(graph.num_nodes()+1);
        for(int i=0;i<node_locks.size();i++) {
            omp_init_lock(&node_locks[i]);
            omp_init_lock(&outgoing_locks[i]);
        }
#endif

        timer.reset();

        #pragma omp parallel for num_threads(get_num_threads()) //default(none) shared(encoded,node_locks,routing_table,incoming_table)
        for (size_t i = 0; i < sequences.size(); i++) {
            int tid = omp_get_thread_num() + 1;
            std::vector<std::tuple<node_index, char, char>> bifurcations;

            const auto &sequence = sequences[i];
            //auto &path_for_sequence = path_for_sequences[i];
            auto path_for_sequence = transform_sequence(sequence);


            size_t kmer_end = graph.get_k();
            graph.map_to_nodes(sequence, [&](node_index node) {
                if (kmer_end < sequence.size()) {
                    path_for_sequence[kmer_end] = routing_table.transform(node, path_for_sequence[kmer_end]);
                }
                kmer_end++;
            });

            size_t kmer_begin = 0;
            kmer_end = graph.get_k();

            vector<node_index> debug_list_of_nodes;
            vector<int> debug_relative_position_history;
            bool debug_transformation_used = sequence != path_for_sequence;
            // TODO: use map_to_nodes_sequentially when implemented
            graph.map_to_nodes(path_for_sequence, [&](node_index node) {
                if (not node) {
                    auto pn = [&](int offset) {
                        PRINT_VAR(sequence.substr(offset,graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(sequence.substr(offset,graph.get_k())));
                        PRINT_VAR(path_for_sequence.substr(offset,graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(path_for_sequence.substr(offset,graph.get_k())));
                    };
                    pn(kmer_begin);
                    pn(kmer_begin-6);
                    auto gp = GraphPreprocessor(graph);
                    PRINT_VAR(bool(gp.is_weak_split(271541)));
                }
                assert(node);
                debug_list_of_nodes.push_back(node);
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

            int relative_position = INT_MIN;
#ifdef _OPENMP
            optional<omp_lock_t*> prev_outgoing_lock;
#endif
            int debug_bifurcation_idx = 0;
            int prev_node = 0;
            char traversed_edge = '\0';
            for (const auto &[node, join_symbol, split_symbol] : bifurcations) {
#ifdef _OPENMP
                omp_set_lock(&node_locks[node]);
                if (prev_node) {
                    bool first_it = 1;
                    bool me_first = 0;
                    bool was_me = 0;
                    omp_set_lock(&outgoing_locks[prev_node]);
                    auto& target_queue = waiting_threads[prev_node];
                    for(auto&[their_thread_id,their_relative_position,their_traversed_edge] : target_queue) {
                        if (tid == their_thread_id) {
                            was_me = 1;
                            if (first_it) {
                                me_first = 1;
                            }
                            their_thread_id = -tid;
                        }
                        else if (was_me and
                                their_traversed_edge == traversed_edge and
                                their_thread_id < 0 and
                                their_relative_position <= relative_position) {
                            cerr << target_queue << endl;
                            PRINT_VAR(tid,node,their_relative_position);
                            relative_position++;
                        }
                        first_it = 0;
                    }
                    if (me_first) {
                        while (!target_queue.empty() and get<0>(target_queue.front()) < 0) {
                            target_queue.pop_front();
                        }
                    }
                    omp_unset_lock(&outgoing_locks[prev_node]);
                }
#endif
                if (join_symbol) {
                    if (join_symbol == '$') {
                        // always putting new read above all other reads
                        relative_position = incoming_table.branch_size(node, '$');
                        encoded[i] = {node, relative_position};
                    }
                    assert(relative_position>=0);
                    if (relative_position > incoming_table.branch_size(node,join_symbol)) {
                        using TT = DecodeEnabler<PathDatabaseDynamicCore<>>;
                        auto self = reinterpret_cast<TT*>(this);

                        cerr << "current" << endl;
                        PRINT_VAR(graph.get_node_sequence(node));
                        incoming_table.print_content(node);
                        PRINT_VAR(join_symbol);
                        cerr << "previous" << endl;
                        auto& prev_bifurcation = bifurcations[debug_bifurcation_idx-1];
                        const auto &[prev_node_2, prev_join_symbol, prev_split_symbol] = prev_bifurcation;
                        PRINT_VAR(graph.get_node_sequence(prev_node));
                        if (prev_split_symbol) {
                            routing_table.print_content(prev_node);
                            PRINT_VAR(debug_relative_position_history.back());
                            PRINT_VAR(prev_split_symbol);
                        }
                        if (prev_join_symbol) {
                            incoming_table.print_content(prev_node);
                            int prev_position = !prev_split_symbol ?
                                    debug_relative_position_history.back()
                                    :
                                    *(debug_relative_position_history.end()-2)
                                    ;
                            auto path_id = self->get_global_path_id(prev_node,prev_position-1);
                            PRINT_VAR(transform_sequence(decode_from_input(sequences,
                                    {graph.get_node_sequence(path_id.first),path_id.second}
                                    ,graph.get_k())
                                    ));
                            PRINT_VAR(prev_position);
                            PRINT_VAR(prev_join_symbol);
                        }
                    }
                    assert(relative_position <= incoming_table.branch_size(node,join_symbol));
                    debug_relative_position_history.push_back(relative_position);
                    relative_position += incoming_table.branch_offset_and_increment(node, join_symbol);
                }


                if (split_symbol) {
                    assert(relative_position>=0);
                    routing_table.insert(node, relative_position, split_symbol);
                    debug_relative_position_history.push_back(relative_position);
                    relative_position = routing_table.new_relative_position(node, relative_position);
                    if (split_symbol == '$') {
                        //++progress_bar;
                    }
                }
#ifdef _OPENMP
                traversed_edge = split_symbol ? split_symbol : 'X'; // X as it doesn't matter
                if (traversed_edge != '$') { // doesn't need any update
                    omp_set_lock(&outgoing_locks[node]);
                    waiting_threads[node].push_back({tid, relative_position, traversed_edge});
                    omp_unset_lock(&outgoing_locks[node]);
                }
                omp_unset_lock(&node_locks[node]);
#endif
                debug_bifurcation_idx++;
                prev_node = node;
            }
        }

#ifdef _OPENMP
        for(int i=0;i<node_locks.size();i++) {
            omp_destroy_lock(&node_locks[i]);
        }
#endif
        encoded_paths += encoded.size();
        statistics["routing_time"] = timer.elapsed();
        std::cerr << "Finished routing in " << statistics["routing_time"] << " sec" << std::endl;

        return encoded;
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

    size_t num_paths() const {
        return encoded_paths;
    };

    json get_statistics(unsigned int verbosity = ~0u) const {
        //auto result = PathDatabase<pair<node_index,int>,GraphT>::get_statistics(verbosity);
        //result.update(statistics);
        //return result;
        return statistics;
    }


    //todo put to proper place
    string transform_sequence(const string& sequence) const {
        transform_done = 0;
        string result = sequence;
        size_t kmer_end = graph.get_k();
        graph.map_to_nodes(sequence, [&](node_index node) {
            if (kmer_end < sequence.size()) {
                result[kmer_end] = routing_table.transform(node, result[kmer_end]);
            }
            kmer_end++;
        });
        transform_done = 1;
        return result;
    }

    void serialize(const fs::path& folder) const {};

//protected:
    json statistics;

    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)


    tsl::hopscotch_set<node_index> additional_joins;
    tsl::hopscotch_set<node_index> additional_splits;

#ifdef MEMOIZE
    vector<char> is_join;
    vector<char> is_split;
#endif
    const GraphT& graph;
    DynamicRoutingTable<> routing_table;
    DynamicIncomingTable<> incoming_table;

    static DBGSuccinct* buildGraph(PathDatabaseDynamicCore* self,vector<string> reads,int kmer_length) {
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
