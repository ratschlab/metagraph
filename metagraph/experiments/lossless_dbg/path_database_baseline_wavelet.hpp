//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_wavelet_hpp
#define path_database_baseline_wavelet_hpp

#include "graph_patch.hpp"
#include "path_database.hpp"
#include "path_database_baseline.hpp"
#include "cxx-prettyprint.hpp"

#include "utils.hpp"
#include "utilities.hpp"
#include <iostream>
#include <set>
#include <functional>
#include <map>
#include "alphabets.hpp"
#include "routing_table.hpp"
#include "incoming_table.hpp"
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
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
class PathDatabaseBaselineWavelet : public PathDatabaseBaseline<DBGSuccinct> {
public:
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaselineWavelet(std::shared_ptr<const DBGSuccinct> graph) : PathDatabaseBaseline(graph),
                                                                              incoming_table(graph)
                                                                              {}
                                                                              
    PathDatabaseBaselineWavelet(const vector<string> &raw_reads,
                                size_t k_kmer = 21 /* default */) : PathDatabaseBaseline(raw_reads,k_kmer),
                                                                    incoming_table(graph)
                                                                    {}


    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {

        vector<path_id> encoded = PathDatabaseBaseline::encode(sequences);
            // add additional bifurcation
        construct_routing_table();
        construct_edge_multiplicity_table();
        return encoded;
    }

    void construct_routing_table() {
        vector<char> routing_table_array;
        for(int node=1;node<=graph.num_nodes();node++) {
            routing_table_array.push_back('#');// to always start a block with #
            if (PathDatabaseBaseline::node_is_split(node)) {
                auto& dynamic_table = PathDatabaseBaseline::routing_table;
                for(int i=0;i<dynamic_table.size(node);i++) {
                    routing_table_array.push_back(dynamic_table.get(node,i));
                }
            }
        }
        routing_table_array.push_back('#'); // to also always end a block with #
        //routing_table_array.push_back('\0'); // end of sequence

        routing_table = RoutingTable(routing_table_array);
    }

    void construct_edge_multiplicity_table() {
        vector<int> edge_multiplicity_table_builder;
        vector<bool> is_join_node;
        for(int node=1;node<=graph.num_nodes();node++) {
            is_join_node.push_back(1);
            if (PathDatabaseBaseline::node_is_join(node)) {
                auto new_reads = PathDatabaseBaseline::incoming_table.branch_size(node,'$');
                if (new_reads) {
                    is_join_node.push_back(0);
                    edge_multiplicity_table_builder.push_back(new_reads);
                }
#ifdef ALL_EDGES_COVERED
                for(auto& base : "ACGTN") {
                    auto branch_size = PathDatabaseBaseline::incoming_table.branch_size(node,base);
                    if (branch_size) {
                        // so it is an actual edge in a graph (because all edges are covered)
                        is_join_node.push_back(0);
                        edge_multiplicity_table_builder.push_back(branch_size);
                    }
                }
#else
                graph.call_incoming_kmers_mine(node,[&node,&edge_multiplicity_table_builder,
                        &is_join_node,this](node_index prev_node,char c) {
                    auto branch_size = PathDatabaseBaseline::incoming_table.branch_size(node,c);
                    is_join_node.push_back(0);
                    edge_multiplicity_table_builder.push_back(branch_size);
                });
#endif
            }
        }
        is_join_node.push_back(1); // to also always end a block with 1
        incoming_table.edge_multiplicity_table = sdsl::enc_vector<>(edge_multiplicity_table_builder);
        sdsl::bit_vector temporary_representation(is_join_node.size());
        for(int i=0;i<is_join_node.size();i++) {
            temporary_representation[i] = is_join_node[i];
        }
        incoming_table.joins = BitVector(temporary_representation);
    }

    using range_t = pair<int,int>;
    using history_t = vector<range_t>;
    using score_t = int;
    using next_nodes_with_extended_info_t =  map<node_index,pair<score_t,history_t>>;

    history_t get_initial_history(node_index node) const {
        history_t history = {{0, get_coverage(node)}};
        return history;
    }

    static bool range_is_empty(range_t range) {
        return range.first >= range.second;
    }

    history_t gather_history(const string& str_history) {
        // todo: add function to gather only consistent history (with strong aka first class, 0-th order support),
        //       now we are gathering reads that are inconsistent with longest path
        // todo: separate in two
        node_index node = graph.kmer_to_node(str_history.substr(0,graph.get_k()));
        assert(node);
        auto history = get_initial_history(node);
        for(auto& c : str_history.substr(graph.get_k())) {
            auto nodes_with_support = get_next_nodes_with_support(node,history);
            node = graph.traverse(node,c);
            history = nodes_with_support[node].second;
            if (history.empty()) {
                break;
            }
        }
        return history;
    }

    node_index get_next_consistent_node(node_index node) {
        history_t history = get_initial_history(node);
        return get_next_consistent_node(node,history);
    }

    node_index get_next_consistent_node(node_index node,const string& str_history) {
        // consistent node should have score 0 and be the only node to go to
        assert(node == graph.kmer_to_node(str_history.substr(str_history.size()-graph.get_k())));
        auto history = gather_history(str_history);
        auto support = get_next_nodes_with_support(node,history);
        node_index consistent_node = 0;
        for(auto& [next_node,info] : support) {
            auto& [score,history] = info;
            if (score == 0) {
                if (!consistent_node) {
                    consistent_node = next_node;
                }
                else {
                    consistent_node = 0;
                    break;
                }
            }
            else {
                // histories are sorted
                break;
            }
        }
        return consistent_node;
    }



    next_nodes_with_extended_info_t get_next_nodes_with_support(node_index node,history_t& history) {

        // todo: merge ranges for successive joining reads
        next_nodes_with_extended_info_t result;
        if (node_is_split(node)) {
            int range_score = 0;
            for(auto& range : history) {
                for(auto& c : "ACGT"s) {
                    range_t new_range;
                    new_range.first = routing_table.rank(node,range.first,c);
                    new_range.second = routing_table.rank(node,range.second,c);
                    if (!range_is_empty(new_range)) {
                        auto new_node = graph.traverse(node,c);
                        if (!result.count(new_node)) {
                            result[new_node].first = range_score;
                        }
                        result[new_node].second.push_back(new_range);
                    }
                }
                range_score++;
            }
        }
        else {
            auto base = graph.get_outgoing_base(node);
            auto new_node = graph.traverse(node,base);
            result[new_node].first = 0;
            result[new_node].second = history;
        }
        for(auto& [new_node,additional_info] : result) {
            auto& [score,result_history] = additional_info;
            if (node_is_join(new_node)) {
                auto offset = incoming_table.branch_offset(new_node,node);
                for(auto& range : result_history) {
                    range.first += offset;
                    range.second += offset;
                }
                if (number_of_reads_starting_at_node(new_node)) {
                    result_history.push_back({0,number_of_reads_starting_at_node(new_node)});
                }
            }
        }
        return result;
    }

    vector<path_id> get_paths_going_through(node_index node) const override  {
        //catch: relative indices in node can be from the same sequence if the read is going there multiple times
        //use set to filter out duplicates or be smarter and stop when arrived to the starting node again
        //todo: improve speed by not getting global path for the same reads
        auto coverage = get_coverage(node);
        set<path_id> out;
        for(int position=0;position<coverage;position++) {
            out.insert(get_global_path_id(node,position));
        }
        return vector<path_id>(all(out));
    }

    std::string decode(path_id path) const override {
        auto node = path.first;
        auto kmer = graph.get_node_sequence(node);
        string sequence = kmer;

        int relative_starting_position = path.second;
        int relative_position = incoming_table.branch_offset(node,0) + relative_starting_position;



        int kmer_position = 0;
        node_index prev_node;
        char base;
        while(true) {
            if (node_is_split(node)) {
                base = routing_table.get(node,relative_position);

#if defined(CHECK_CORECTNESS)
                auto& routing_table_naive = splits.at(node);
                auto rt_index = routing_table_naive.begin();
                advance(rt_index,relative_position);
                char base_check = *rt_index;
                auto new_relative_position = rank(routing_table_naive,base,relative_position)-1;
#endif
                relative_position = routing_table.rank(node,relative_position,base);

#if defined(CHECK_CORECTNESS)
                assert(base_check == base);
                assert(relative_position == new_relative_position);
#endif

            }
            else {
                assert(graph.outdegree(node) == 1);
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
            }
            if (base == '$') break;
            prev_node = node;
            node = graph.traverse(node,base);
            assert(node);
            kmer_position++;
            sequence.append(1,base); // 1 times base
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];

#if defined(CHECK_CORECTNESS)
                auto tv = PathDatabaseBaseline::branch_starting_offset(node,join_symbol);
                auto cv = incoming_table.branch_offset(node,prev_node);
                assert(tv==cv);
#endif

                relative_position += incoming_table.branch_offset(node,prev_node);
            }
        }

        return sequence;
    }

    char node_get_last_char(node_index node) const {
        auto kmer = graph.get_node_sequence(node);
        return kmer[kmer.size()-1];
    }

    char node_get_first_char(node_index node) const {
        auto kmer = graph.get_node_sequence(node);
        return kmer.front();

    }


    int get_coverage(node_index node) const {
        // one can also traverse backwards
        if (node_is_split(node)) {
            return routing_table.size(node);
        }

        char base;
        graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});

        auto new_node = graph.traverse(node,base);

        if (node_is_join(new_node)) {
            return incoming_table.branch_size(new_node,node);
        }

        return get_coverage(new_node);
    }

    vector<string> decode_all_reads() const {
        auto reads = vector<string>();
        for(node_index node=1;node<=graph.num_nodes();node++) {
            auto count = number_of_reads_starting_at_node(node);
            for(auto relative_index=0;relative_index<count;relative_index++) {
                reads.push_back(decode({node,relative_index}));
            }
        }
        return reads;
    }

    vector<string> decode_all_reads_inverse() const {
        auto reads = vector<string>();
        for(node_index node=1;node<=graph.num_nodes();node++) {
            auto count = number_of_reads_ending_at_node(node);
            for(auto read_index=0;read_index<count;read_index++) {
                auto relative_index = routing_table.select(node,read_index+1,'$');
                auto id = get_global_path_id(node,relative_index);
                reads.push_back(decode(id));
            }
        }
        return reads;
    }

    bool is_valid_path_id(path_id path_id) const {
        return node_is_join(path_id.first) && path_id.second < incoming_table.branch_size(path_id.first,0);
    }

    int number_of_reads_starting_at_node(node_index node) const {
        int result = 0;
        if (node_is_join(node)) {
            result = incoming_table.branch_size(node,0);
        }
        return result;
    }

    int number_of_reads_ending_at_node(node_index node) const {
        int result = 0;
        if (node_is_split(node)) {
            auto size = routing_table.select(node,1,'#');
            result = routing_table.rank(node,size,'$');
        }
        return result;
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



    void serialize(const fs::path& folder) const override {
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

    static PathDatabaseBaselineWavelet deserialize(const fs::path& folder) {
        ifstream edge_multiplicity_file(folder / "edge_multiplicity.bin");
        ifstream routing_table_file(folder / "routing_table.bin");
        ifstream joins_file(folder / "joins.bin");
        string graph_filename = folder / "graph.bin";

        auto graph = std::shared_ptr<DeBruijnGraph>{
                new DBGSuccinct(21)
                };
        graph->load(graph_filename);
        auto db = PathDatabaseBaselineWavelet(graph);
        db.incoming_table.edge_multiplicity_table.load(edge_multiplicity_file);
        db.routing_table.load(routing_table_file);
        db.incoming_table.joins.load(joins_file);
        return db;
    }

    json get_statistics(unsigned int verbosity = ~0u) const {
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
        json result = {{"true_joins", true_joins},
                       {"added_joins", added_joins},
                       {"true_splits", true_splits},
                       {"added_splits", added_splits},
                       {"num_of_nodes", graph.num_nodes()}
                      };
        if (verbosity & STATS_SPLITS_HISTOGRAM) {
            result["splits_diff_symbols_histogram"] = splits_diff_symbols_histogram;
            result["splits_size_histogram"] = splits_size_histogram;
        }
        if (verbosity & STATS_JOINS_HISTOGRAM) {
            result["joins_diff_symbols_histogram"] = joins_diff_symbols_histogram;
        }
        cerr << result.dump(4) << endl;
        return result;
    }

    path_id get_global_path_id(node_index node, int relative_position) const {
        node_index prev_node = 0;
        int prev_offset = 0;
        if (node_is_join(node)) {
            graph.call_incoming_kmers_mine(node,[&](node_index possible_node,char c) {
                auto offset = incoming_table.branch_offset(node,possible_node);
                if (offset <= relative_position && offset >= prev_offset) {
                    prev_node = possible_node;
                    prev_offset = offset;
                }
            });
            relative_position -= prev_offset;
        }
        else {
            assert(graph.indegree(node) == 1);
            graph.call_incoming_kmers_mine(node,[&](node_index possible_node,char c) {
                prev_node = possible_node;
            });
        }
        if (!prev_node) {
            assert(is_valid_path_id({node,relative_position}));
            return {node,relative_position};
        }
        assert(prev_node);
        if (node_is_split(prev_node)) {
            relative_position = routing_table.select(prev_node,relative_position+1,node_get_last_char(node));// +1 as relative_position is 0-based
        }
        return get_global_path_id(prev_node,relative_position);
    }

private:
    RoutingTable<Wavelet> routing_table;
    IncomingTable<BitVector> incoming_table;
};

#endif /* path_database_baseline_hpp */
