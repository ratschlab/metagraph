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
#include "threading.hpp"


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
                incoming_table(*(this->graph_)) {}

    explicit PathDatabaseDynamic(const vector<string> &filenames,
                 size_t k_kmer = 21 /* default kmer */) :
                         PathDatabase<pair<node_index,int>,GraphT>(filenames, k_kmer),
                         graph(*(this->graph_)),
                         incoming_table(*(this->graph_)) {}

    virtual ~PathDatabaseDynamic() {}

    node_index starting_node(const string& sequence) {
        return graph.kmer_to_node(sequence.substr(0, graph.get_k()));
    }

    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {
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
            additional_joins_vec[i] = starting_node(sequence);
            additional_splits_vec[i] = graph.kmer_to_node(
                sequence.substr(sequence.length() - graph.get_k())
            );
        }

        additional_joins.insert(additional_joins_vec.begin(), additional_joins_vec.end());
        additional_splits.insert(additional_splits_vec.begin(), additional_splits_vec.end());

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

        vector<path_id> encoded;

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < sequences.size(); ++i) {
            const auto &sequence = sequences[i];

            // always putting new read above all other reads
            int relative_position = 0;
            int relative_starting_position = -1;

            // TODO: use map_to_nodes_sequentially when implemented
            std::vector<std::tuple<size_t, node_index, bool, bool>> bifurcations;
            size_t i = 0;
            graph.map_to_nodes(sequence, [&](node_index node) {
                bool is_join = node_is_join(node);
                bool is_split = node_is_split(node);
                if (is_join || is_split) {
                    bifurcations.emplace_back(i, node, is_join, is_split);
                }
                i++;
            });

            #pragma omp critical
            {
                for (const auto &[kmer_begin, node, is_join, is_split] : bifurcations) {

                    auto kmer_end = kmer_begin + graph.get_k();

                    if (is_join) {
                        auto join_symbol = kmer_begin ? sequence[kmer_begin - 1] : '$';

                        if (join_symbol == '$') {
                            relative_starting_position = incoming_table.branch_size(node, '$');
                            relative_position += relative_starting_position;
                        }

                        relative_position += incoming_table.branch_offset_and_increment(node, join_symbol);
                    }

                    if (is_split) {
                        auto split_symbol = kmer_end < sequence.size() ? sequence[kmer_end] : '$';
                        routing_table.insert(node, relative_position, split_symbol);
                        relative_position = routing_table.rank(node, split_symbol, relative_position);
                    }
                }

                encoded.emplace_back(starting_node(sequence), relative_starting_position);
                encoded_paths++;
            }
        }

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
    sdsl::bit_vector is_join;
    sdsl::bit_vector is_split;
#endif
};

#endif /* path_database_baseline_hpp */
