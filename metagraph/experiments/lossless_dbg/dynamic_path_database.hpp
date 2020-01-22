//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef __PATH_DATABASE_DYNAMIC_HPP__
#define __PATH_DATABASE_DYNAMIC_HPP__

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
#include "utilities.hpp"
#include "unix_tools.hpp"
#include "threading.hpp"

using node_index = DeBruijnGraph::node_index;


using DefaultDynamicRoutingTable = DynamicRoutingTable<>;
using DefaultDynamicIncomingTable = DynamicIncomingTable<>;
CREATE_MEMBER_CHECK(transformations);


template <typename RoutingTableT = DefaultDynamicRoutingTable,
          typename IncomingTableT = DefaultDynamicIncomingTable,
          typename ExitBarrierT = ExitBarrier<>>
class DynamicPathDatabaseCore {
  public:
    using GraphT = DBGSuccinct;
    using path_id = pair<node_index, int64_t>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit DynamicPathDatabaseCore(std::shared_ptr<const GraphT> graph,
                                     int64_t chunks = DefaultChunks)
        : graph_(graph),
          graph(*graph_),
          incoming_table(graph_),
          routing_table(graph_), // weak
          chunks(chunks) {
        assert(chunks > 0);
    }


    explicit DynamicPathDatabaseCore(const vector<string> &reads,
                                     size_t kmer_length = 21 /* default kmer */)
        : graph_(shared_ptr<const GraphT>(buildGraph(this, reads, kmer_length))),
          graph(*graph_),
          incoming_table(graph_),
          routing_table(graph_),
          chunks(DefaultChunks) {}

    virtual ~DynamicPathDatabaseCore() {}


    std::vector<path_id> encode(const std::vector<std::string> &sequences) {
        statistics["encode_initial_ram"] = get_used_memory();
        VerboseTimer encoding_timer("encoding reads");
        VerboseTimer preprocessing_timer("preprocessing step");
        // improvement
        // - when multiple reads start at a same symbol, sort them so they share longest prefix
        // sort them so we have fixed relative ordering
        // probably not need to sort globally as we want only relative ordering

        // decide if the transformed sequences will be different or not and either copy or point to original sequences.
        typename std::conditional<has_member_transformations<RoutingTableT>::value,
                                  vector<string>, const vector<string> &>::type transformed_sequences
                = sequences;
        if constexpr (has_member_transformations<RoutingTableT>::value) {
#pragma omp parallel for
            for (uint64_t i = 0; i < sequences.size(); i++) {
                transform_sequence_inplace(transformed_sequences[i]);
            }
        }
        statistics["transform_sequences_ram"] = get_used_memory();
        populate_additional_joins(transformed_sequences);

        auto alloc_routing_table = VerboseTimer("allocation of routing & incoming table");

        routing_table = decltype(routing_table)(graph_, &is_bifurcation,
                                                &rank_is_bifurcation, chunks);
        incoming_table = decltype(incoming_table)(graph_, &is_bifurcation,
                                                  &rank_is_bifurcation, chunks);
        alloc_routing_table.finished();


        vector<path_id> encoded(sequences.size());

#ifndef DISABLE_PARALELIZATION

        auto alloc_lock_t = VerboseTimer("memory allocation of locks");
        ChunkedDenseHashMap<omp_lock_t, decltype(is_bifurcation), decltype(rank_is_bifurcation), false>
                node_locks(&is_bifurcation, &rank_is_bifurcation, chunks);
        alloc_lock_t.finished();


        auto lock_init_timer = VerboseTimer("initializing locks & barriers");
        for (uint64_t i = 0; i < node_locks.elements.size(); i++) {
            omp_init_lock(&node_locks.elements[i]);
        }
        lock_init_timer.finished();

#endif
        auto threads_map_t = VerboseTimer("memory allocation of exit barriers");

        auto exit_barrier = ExitBarrierT(&is_bifurcation, &rank_is_bifurcation, chunks);
        threads_map_t.finished();


        statistics["preprocessing_time"] = preprocessing_timer.finished();
        statistics["preprocessing_ram"] = get_used_memory();

        PRINT_VAR("Cumulative size in bytes", is_join.size() / 8)
        PRINT_VAR("Cumulative size in bytes", is_split.size() / 8)
        VerboseTimer routing_timer("routing step");
#ifndef DISABLE_PARALELIZATION
#pragma omp parallel for num_threads( \
        get_num_threads()) // default(none) shared(encoded,node_locks,routing_table,incoming_table)
#endif
        for (size_t i = 0; i < sequences.size(); i++) {
            int16_t tid = omp_get_thread_num();
            std::vector<std::tuple<node_index, char, char>> bifurcations;

            const auto &sequence = sequences[i];
            const auto &path_for_sequence = transformed_sequences[i];

            size_t kmer_begin = 0;
            auto kmer_end = graph.get_k();
#ifdef DEBUG_ADDITIONAL_INFORMATION
            vector<node_index> debug_list_of_nodes;
            vector<int64_t> debug_relative_position_history;
            bool debug_transformation_used = sequence != path_for_sequence;
#endif

            graph.map_to_nodes_sequentially(all(path_for_sequence), [&](node_index node) {
                assert(node || [&]() {
                    auto pn = [&](int64_t offset) {
                        PRINT_VAR(sequence.substr(offset, graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(sequence.substr(offset, graph.get_k())));
                        PRINT_VAR(path_for_sequence.substr(offset, graph.get_k()));
                        PRINT_VAR(graph.kmer_to_node(
                                path_for_sequence.substr(offset, graph.get_k())));
                    };
                    pn(kmer_begin);
                    return false;
                }() && "sequence doesn't represent a walk");

#ifdef DEBUG_ADDITIONAL_INFORMATION
                debug_list_of_nodes.push_back(node);
#endif
                char join_char = node_is_join(node)
                        ? (kmer_begin ? path_for_sequence[kmer_begin - 1] : '$')
                        : '\0';
                char split_char = node_is_split(node)
                        ? (kmer_end < sequence.size() ? sequence[kmer_end] : '$')
                        : '\0';

                if (join_char || split_char) {
                    bifurcations.emplace_back(node, join_char, split_char);
                }
                kmer_begin++;
                kmer_end++;
            });

            auto &[first_node, first_join, first_split] = bifurcations.front();
            auto &[last_node, last_join, last_split] = bifurcations.back();
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
                    relative_position = exit_barrier.exit(previous_node, tid);
                }

                if (join_symbol) {
                    if (join_symbol == '$') {
                        // always putting new read above all other reads
                        relative_position = incoming_table.branch_size(node, '$');
                        encoded[i] = { node, relative_position };
                    }
                    assert(relative_position >= 0);
                    assert((relative_position <= incoming_table.branch_size(node, join_symbol) ||
                            [&, node = node, join_symbol = join_symbol] {
                                using TT = DecodeEnabler<DynamicPathDatabaseCore<>>;
                                auto self = reinterpret_cast<TT *>(this);
                                (void)self; // so it is always used (because of conditional compilation)
                                mg::common::logger->info("current");
                                PRINT_VAR(relative_position,
                                          incoming_table.branch_size(node, join_symbol));
                                PRINT_VAR(graph.get_node_sequence(node));
                                PRINT_VAR(node, tid, relative_position, previous_node,
                                          traversed_edge);
#ifdef DEBUG_ADDITIONAL_INFORMATION
                                PRINT_VAR(debug_bifurcation_idx);
#endif
                                PRINT_VAR("current incoming table");
                                incoming_table.print_content(node);
                                PRINT_VAR(join_symbol);
                                mg::common::logger->info("previous");
                                PRINT_VAR("exit barrier");
                                exit_barrier.print_content(previous_node);
#ifdef DEBUG_ADDITIONAL_INFORMATION
                                auto &prev_bifurcation
                                        = bifurcations[debug_bifurcation_idx - 1];
                                const auto &[previous_node_2, prev_join_symbol,
                                             prev_split_symbol]
                                        = prev_bifurcation;
                                PRINT_VAR(graph.get_node_sequence(previous_node));
                                if (debug_bifurcation_idx > 0) {
                                    PRINT_VAR(debug_relative_position_history.back());
                                }
                                if (prev_split_symbol) {
                                    PRINT_VAR("previous routing table");
                                    routing_table.print_content(previous_node);
                                    PRINT_VAR(debug_relative_position_history.back());
                                    PRINT_VAR(prev_split_symbol);
                                }
                                if (prev_join_symbol) {
                                    PRINT_VAR("previous incoming table");
                                    incoming_table.print_content(previous_node);
                                    int64_t prev_position = !prev_split_symbol
                                            ? debug_relative_position_history.back()
                                            : *(debug_relative_position_history.end() - 2);
                                    PRINT_VAR(prev_position);
                                    PRINT_VAR(prev_join_symbol);
                                    auto path_id
                                            = self->get_global_path_id(previous_node,
                                                                       prev_position - 1);
                                    PRINT_VAR(transform_sequence(decode_from_input(
                                            sequences,
                                            { graph.get_node_sequence(path_id.first),
                                              path_id.second },
                                            graph.get_k())));
                                }
#endif
                                return false;
                            }()));

#ifdef DEBUG_ADDITIONAL_INFORMATION
                    debug_relative_position_history.push_back(relative_position);
#endif
                    auto lower_routes
                            = incoming_table.branch_offset_and_increment(node, join_symbol);
                    assert(lower_routes >= 0);
                    relative_position += lower_routes;
                }

                if (split_symbol) {
                    assert(relative_position >= 0);
                    auto block_offset
                            = routing_table.insert(node, relative_position, split_symbol);
#ifdef DEBUG_ADDITIONAL_INFORMATION
                    debug_relative_position_history.push_back(relative_position);
#endif
                    relative_position
                            = routing_table.new_relative_position(node, relative_position,
                                                                  block_offset);
                }
                traversed_edge = split_symbol ? split_symbol : 'X'; // X as it doesn't matter
                if (traversed_edge != '$') { // doesn't need any update
                    exit_barrier.enter(node, traversed_edge, relative_position, tid);
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
        for (uint64_t i = 0; i < node_locks.elements.size(); i++) {
            omp_destroy_lock(&node_locks.elements[i]);
        }
#endif
        encoded_paths += encoded.size();
        statistics["routing_time"] = routing_timer.finished();
        statistics["routing_ram"] = get_used_memory();
        encoding_timer.finished();

        return encoded;
    }

#ifndef DISABLE_PARALELIZATION


#endif

    void populate_additional_joins(const vector<std::string> &transformed_sequences) {
        // change
        ll bits_to_set = transformed_sequences.size();

        // vector<node_index> debug_join_ids(bits_to_set);
        // vector<node_index> debug_split_ids(bits_to_set);
        auto bifurcation_timer = VerboseTimer("construction of bifurcation bit_vectors");

        auto additional_bifurcations_timer
                = VerboseTimer("computing additional splits and joins");
        is_split = decltype(is_split)(graph.num_nodes() + 1); // bit
        is_join = decltype(is_join)(graph.num_nodes() + 1);
        is_bifurcation = decltype(is_join)(graph.num_nodes() + 1);
        uint64_t chunk_size = (graph.num_nodes() + 1) / chunks + 64ull;
        chunk_size &= ~63ull;
        vector<omp_lock_t> node_locks(chunks);
        for (uint64_t i = 0; i < node_locks.size(); i++) {
            omp_init_lock(&node_locks[i]);
        }
#pragma omp parallel for
        for (int64_t i = 0; i < bits_to_set; ++i) {
            auto &transformed_sequence = transformed_sequences[i];
            omp_lock_t *lock_ptr;
            ll start_node
                    = graph.kmer_to_node(transformed_sequence.substr(0, graph.get_k()));
            assert(start_node);
            lock_ptr = &node_locks[start_node / chunk_size];
            omp_set_lock(lock_ptr);
            is_join[start_node] = true;
            omp_unset_lock(lock_ptr);
            // debug_join_ids[i] = start_node;

            ll end_node = graph.kmer_to_node(transformed_sequence.substr(
                    transformed_sequence.length() - graph.get_k()));
            assert(end_node);
            lock_ptr = &node_locks[end_node / chunk_size];
            omp_set_lock(lock_ptr);
            is_split[end_node] = true;
            omp_unset_lock(lock_ptr);
            // debug_split_ids[i] = end_node;
        }
        for (uint64_t i = 0; i < node_locks.size(); i++) {
            omp_destroy_lock(&node_locks[i]);
        }
        statistics["additional_bifurcations_time"]
                = additional_bifurcations_timer.finished();
        statistics["additional_bifurcations_ram"] = get_used_memory();

//#pragma omp parallel for reduction(append : debug_join_ids, debug_split_ids)
#pragma omp parallel for
        for (uint64_t id = 0; id <= graph.num_nodes(); id += 64) {
            for (uint64_t node = id; node < id + 64 && node <= graph.num_nodes(); ++node) {
                if (!node)
                    continue;
                auto outdegree = graph.outdegree(node);
                is_split[node] = is_split[node] or outdegree > 1;
                //            if (outdegree > 1) {
                //                debug_split_ids.push_back(node);
                //            }
                auto indegree = graph.indegree(node);
                is_join[node] = is_join[node] or indegree > 1;
                //            if (indegree > 1) {
                //                debug_join_ids.push_back(node);
                //            }
                is_bifurcation[node] = is_split[node] || is_join[node];
            }
        }
        rank_is_split = decltype(rank_is_split)(&is_split);
        rank_is_join = decltype(rank_is_join)(&is_join);
        rank_is_bifurcation = decltype(rank_is_bifurcation)(&is_bifurcation);

        statistics["bifurcations_time"] = bifurcation_timer.finished();
        statistics["bifurcations_ram"] = get_used_memory();
    }


    bool node_is_join(node_index node) const { return is_join[node]; }

    bool node_is_split(node_index node) const { return is_split[node]; }

    void compress() {}

    size_t num_paths() const { return encoded_paths; };

    json get_statistics([[maybe_unused]] uint64_t verbosity = ~0u) const {
        // auto result = PathDatabase<pair<node_index,int>,GraphT>::get_statistics(verbosity);
        // result.update(statistics);
        // return result;
        return statistics;
    }

    string transform_sequence(const string &sequence) const {
        string result = sequence;
        transform_sequence_inplace(result);
        return result;
    }

    void transform_sequence_inplace(string &sequence) const {
        size_t kmer_end = graph.get_k();
        static_assert(has_member_transformations<RoutingTableT>::value
                      == std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,
                                         RoutingTableT>::value);
        if constexpr (std::is_base_of<TransformationsEnabler<DynamicRoutingTableCore<>>,
                                      RoutingTableT>::value) {
            graph.map_to_nodes_sequentially(all(sequence), [&](node_index node) {
                if (kmer_end < sequence.size()) {
                    sequence[kmer_end] = routing_table.transform(node, sequence[kmer_end]);
                }
                kmer_end++;
            });
        }
    }

    void serialize(const fs::path &) const {};

    // protected:
    json statistics;

    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int64_t encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)


    sdsl::bit_vector is_join;
    sdsl::bit_vector is_split;
    sdsl::bit_vector is_bifurcation;
    typename decltype(is_join)::rank_1_type rank_is_join;
    typename decltype(is_split)::rank_1_type rank_is_split;
    typename decltype(is_bifurcation)::rank_1_type rank_is_bifurcation;

    std::shared_ptr<const GraphT> graph_;
    const GraphT &graph;
    IncomingTableT incoming_table;
    RoutingTableT routing_table;
    int64_t chunks;

    static DBGSuccinct *
    buildGraph(DynamicPathDatabaseCore *self, vector<string> reads, int64_t kmer_length) {
        VerboseTimer building_graph_timer("building the graph");
        auto graph = new DBGSuccinct(dbg_succ_graph_constructor(reads, kmer_length));
#if defined(MASK_DUMMY_KMERS)
        graph->mask_dummy_kmers(1, false);
#endif
        auto elapsed = building_graph_timer.finished();
        self->statistics["graph_build_time"] = elapsed;
        return graph;
    }
};

template <typename RoutingTableT = DefaultDynamicRoutingTable, typename IncomingTableT = DefaultDynamicIncomingTable>
class PathDatabaseDynamic
    : public DecodeEnabler<DynamicPathDatabaseCore<RoutingTableT, IncomingTableT>> {
    using DecodeEnabler<DynamicPathDatabaseCore<RoutingTableT, IncomingTableT>>::DecodeEnabler;
};


#endif /* __PATH_DATABASE_DYNAMIC_HPP__ */
