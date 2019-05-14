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
#include <progress_bar.hpp>

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

template<typename GraphT=DBGSuccinct>
class PathDatabaseDynamic : public PathDatabase<pair<node_index,int>,GraphT> {
public:
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    explicit PathDatabaseDynamic(std::shared_ptr<const GraphT> graph) :
            PathDatabase<pair<node_index,int>,GraphT>(graph),
            graph(*(this->graph_)),
            incoming_table(*(this->graph_)), routing_table(nullptr) {}

    explicit PathDatabaseDynamic(const vector<string> &filenames,
                 size_t k_kmer = 21 /* default kmer */) :
            PathDatabase<pair<node_index,int>,GraphT>(filenames, k_kmer),
            graph(*(this->graph_)),
            incoming_table(*(this->graph_)),
            //routing_table((new GraphPreprocessor(*(this->graph_)))->find_weak_splits())
            routing_table()
            {}

    virtual ~PathDatabaseDynamic() {}

    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {
        Timer timer;

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
            // add additional bifurcation
            additional_joins_vec[i] = graph.kmer_to_node(
                sequence.substr(0, graph.get_k())
            );
            additional_splits_vec[i] = graph.kmer_to_node(
                sequence.substr(sequence.length() - graph.get_k())
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
        is_split = sdsl::bit_vector(graph.num_nodes() + 1);
        is_join = sdsl::bit_vector(graph.num_nodes() + 1);

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

        vector<path_id> encoded(sequences.size());

        vector<string> path_for_sequences(sequences.size());

        std::cout << "Finished preprocessing in " << timer.elapsed() << " sec" << std::endl;
        //ProgressBar progress_bar(sequences.size(), "Building dRT and dEM");

        const size_t batch_size = std::max(static_cast<size_t>(1),
                                           static_cast<size_t>(sequences.size() / 200));

        double join_time = 0;
        double split_time = 0;

#ifdef USE_LOCKS
        vector<omp_lock_t> locks(graph.num_nodes()+1);
        for(int i=0;i<locks.size();i++) {
            omp_init_lock(&locks[i]);
        }
#endif

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < sequences.size(); i++) {
            std::vector<std::tuple<node_index, char, char>> bifurcations;

            const auto &sequence = sequences[i];
            auto &path_for_sequence = path_for_sequences[i];
            path_for_sequence = "$" + sequence + "$";


            size_t kmer_begin = 0;
            size_t kmer_end = graph.get_k();

            // TODO: use map_to_nodes_sequentially when implemented
            graph.map_to_nodes(sequence, [&](node_index node) {
                path_for_sequence[kmer_end + 1] = routing_table.transform(node, path_for_sequence[kmer_end + 1]);
                char join_char = node_is_join(node)
                                 ? path_for_sequence[kmer_begin]
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
            for (const auto &[node, join_symbol, split_symbol] : bifurcations) {
#ifdef USE_LOCKS
                omp_set_lock(&locks[node]);
#endif
                if (join_symbol) {
                    if (join_symbol == '$') {
                        // always putting new read above all other reads
                        relative_position = incoming_table.branch_size(node, '$');
                        encoded[i] = {node, relative_position};
                    }
                    relative_position += incoming_table.branch_offset_and_increment(node, join_symbol);
                }

                if (split_symbol) {
                    routing_table.insert(node, relative_position, split_symbol);
                    relative_position = routing_table.new_relative_position(node, relative_position);
                    if (split_symbol == '$') {
                        //++progress_bar;
                    }
                }
#ifdef USE_LOCKS
                omp_unset_lock(&locks[node]);
#endif
            }
        }

#ifdef USE_LOCKS
        for(int i=0;i<locks.size();i++) {
            omp_destroy_lock(&locks[i]);
        }
#endif
        encoded_paths += encoded.size();
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

    size_t num_paths() const override {
        return encoded_paths;
    };

    std::string decode(path_id path) const override {
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
            node = graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,encoded_base); // 1 times base
            sequence_path.append(1,base); // 1 times base

            if (node_is_join(node)) {
                auto join_symbol = sequence_path[kmer_position-1];
                relative_position += incoming_table.branch_offset(node,join_symbol);
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
    // denote how many reads are joining from every branch ($ATCGN) ($ denotes start of a new read)
    int encoded_paths = 0;
    // denote where the reads should go ($ATCGN) ($ denodes the end of particular read)

    DynamicRoutingTable routing_table;
    DynamicIncomingTable<> incoming_table;

    tsl::hopscotch_set<node_index> additional_joins;
    tsl::hopscotch_set<node_index> additional_splits;

    const GraphT & graph;

#ifdef MEMOIZE
    sdsl::bit_vector is_join;
    sdsl::bit_vector is_split;
#endif
};





#endif /* path_database_baseline_hpp */
