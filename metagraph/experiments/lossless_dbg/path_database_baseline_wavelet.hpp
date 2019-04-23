//
//  path_database_baseline.hpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#ifndef path_database_baseline_wavelet_hpp
#define path_database_baseline_wavelet_hpp

#include "path_database.hpp"
#include "path_database_baseline.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include <iostream>
#include <set>
#include <map>
#include "alphabets.hpp"
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>


template<typename POD>
std::istream &deserialize(std::istream &is, vector<POD> &v) {
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}

template<typename POD>
std::ostream &serialize(std::ostream &os, const vector<POD> &v) {
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}

#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"

#define STATS_JOINS_HISTOGRAM (1u << 0)
#define STATS_SPLITS_HISTOGRAM (1u << 1)

using namespace std;
using alphabets::log2;

// todo find a tool that removes this relative namespacing issue
// say to Mikhail that "de_bruijn_graph" instead of "metagraph/de_bruijn_graph" is the same violation as this


using routing_character_t = int;
const char RoutingTableAlphabet[] = {'$','A','C','G','T','N','#','?'};
// improvement (constexpr use https://github.com/serge-sans-paille/frozen)
const map<char,int> RoutingTableInverseAlphabet = {{'$',0},{'A',1},{'C',2},{'G',3},{'T',4},{'N',5},{'#',6},{'?',7}};
const auto& rte2int = RoutingTableInverseAlphabet;

int operator""_rc(char c) {
    return RoutingTableInverseAlphabet.at(c);
}
int rc(char c) {
    return RoutingTableInverseAlphabet.at(c);
}

char tochar(routing_character_t rc) {
    return RoutingTableAlphabet[rc];
}

//template <class Wavelet = sdsl::wt_rlmn<sdsl::sd_vector<>>>
template<class Wavelet = sdsl::wt_rlmn<>>
class PathDatabaseBaselineWavelet : public PathDatabaseBaseline {
public:
    using routing_table_t = vector<char>;
    // implicit assumptions
    // graph contains all reads
    // sequences are of size at least k
    PathDatabaseBaselineWavelet(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabaseBaseline(graph)
                                                                              {}

    PathDatabaseBaselineWavelet(const vector<string> &raw_reads,
                                size_t k_kmer = 21 /* default */) : PathDatabaseBaseline(raw_reads,k_kmer)
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
            if (splits.count(node)) {
                for(auto& choice : splits[node]) {
                    routing_table_array.push_back(choice);
                }
            }
        }
        routing_table_array.push_back('#'); // to also always end a block with #
        //routing_table_array.push_back('\0'); // end of sequence
        sdsl::int_vector<0> routing_table_array_encoded(routing_table_array.size());
        for(int i=0;i<routing_table_array.size();i++) {
            routing_table_array_encoded[i] = RoutingTableInverseAlphabet.at(routing_table_array[i]);
        }
        construct_im(routing_table,routing_table_array_encoded,0);
    }

    void construct_edge_multiplicity_table() {
        vector<bool> is_join_node(graph.num_nodes());
        for(int node=1;node<=graph.num_nodes();node++) {
            is_join_node[node-1] = PathDatabaseBaseline::node_is_join(node);
            if (PathDatabaseBaseline::node_is_join(node)) {
                for(int rc=0;rc<'N'_rc;rc++) { // don't need to store last branch as we only compute prefix sum excluding
                                               // the branch which we came from (N in this case)
                    edge_multiplicity_table.push_back(PathDatabaseBaseline::joins[node][tochar(rc)]);
                }
            }
        }
        joins = bit_vector_stat(is_join_node);

    }

    int routing_table_offset(node_index node) const { return routing_table.select(node, '#'_rc) + 1; }

    int routing_table_select(node_index node, int occurrence, char symbol) const {
        auto routing_table_block = routing_table_offset(node);
        auto occurrences_of_symbol_before_block = routing_table.rank(routing_table_block,rc(symbol));
        return routing_table.select(occurrences_of_symbol_before_block+occurrence,rc(symbol)) - routing_table_block;
    }

    int routing_table_rank(node_index node, int position, char symbol) const {
        auto routing_table_block = routing_table_offset(node);
        auto absolute_position = routing_table_block+position;
        auto occurrences_of_base_before_block = routing_table.rank(routing_table_block,rc(symbol));
        return routing_table.rank(absolute_position,rc(symbol)) - occurrences_of_base_before_block;
    }

    char routing_table_get(node_index node, int position) const {
        auto routing_table_block = routing_table_offset(node);
        return tochar(routing_table[routing_table_block+position]);
    }

    int routing_table_size(node_index node) const {
        return routing_table_select(node,1,'#');
    }

    void routing_table_print_content(node_index node) const {
        auto size = routing_table_size(node);
        for (int i=0;i<size;i++) {
            cout << routing_table_get(node,i);
        }
        cout << endl;
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
        int relative_position = branch_starting_offset(node,'$') + relative_starting_position;



        int kmer_position = 0;
        int base;
        while(true) {
            if (node_is_split(node)) {
                auto routing_table_block = routing_table_offset(node);
                auto absolute_position = routing_table_block+relative_position;
                base = routing_table[absolute_position];
                auto occurrences_of_base_before_block = routing_table.rank(routing_table_block,base);

//                //checkers
//                auto& routing_table_naive = splits.at(node);
//                auto rt_index = routing_table_naive.begin();
//                advance(rt_index,relative_position);
//                char base_check = *rt_index;
//                auto new_relative_position = rank(routing_table_naive,tochar(base),relative_position)-1;

                auto rank_of_base = routing_table.rank(absolute_position,base) - occurrences_of_base_before_block;

//                //checkers
//                assert(base_check == tochar(base));
//                assert(rank_of_base == new_relative_position);

                relative_position = rank_of_base;
            }
            else {
                assert(graph.outdegree(node) == 1);
                graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = rc(edge_label);});
            }
            if (base == '$'_rc) break;
            node = graph.traverse(node,tochar(base));
            assert(node);
            kmer_position++;
            sequence.append(1,tochar(base)); // 1 times base
            if (node_is_join(node)) {
                // todo better name (it is a symbol that determines from which branch we came)
                auto join_symbol = sequence[kmer_position-1];
//                auto tv = PathDatabaseBaseline::branch_starting_offset(node,join_symbol);
//                auto cv = branch_starting_offset(node,join_symbol);
//                assert(tv==cv);
                relative_position += branch_starting_offset(node,join_symbol);
            }
        }

        return sequence;
    }


    int get_coverage(node_index node, char inverse_edge_label='?') const {
        if (node_is_split(node)) {
            return routing_table_size(node);
        }
        if (node_is_join(node)) {
            if (inverse_edge_label == '?') {
                return joins_size(node);
            }
            else {
                return branch_size(node,inverse_edge_label);
            }
        }
        char base;
        // one can also traverse backwards
        graph.call_outgoing_kmers(node,[&base](node_index node,char edge_label ) { base = edge_label;});
        auto kmer = graph.get_node_sequence(node);
        auto last_base = kmer[kmer.size()-1];
        return get_coverage(graph.traverse(node,base),last_base);
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
                auto relative_index = routing_table_select(node,read_index+1,'$');
                auto id = get_global_path_id(node,relative_index);
                reads.push_back(decode(id));
            }
        }
        return reads;
    }

    bool is_valid_path_id(path_id path_id) const {
        return node_is_join(path_id.first) && path_id.second < branch_size(path_id.first,'$');
    }

    int number_of_reads_starting_at_node(node_index node) const {
        int result = 0;
        if (node_is_join(node)) {
            int starting_offset = (joins.rank1(node-1)-1)*'N'_rc;
            result = edge_multiplicity_table[starting_offset + '$'_rc];
        }
        return result;
    }

    int number_of_reads_ending_at_node(node_index node) const {
        int result = 0;
        if (node_is_split(node)) {
            auto size = routing_table_select(node,1,'#');
            result = routing_table_rank(node,size,'$');
        }
        return result;
    }

    // is join or start of the read
    bool node_is_join(node_index node) const {
        return joins[node-1]; // 0 based
    }

    // is split or end of the read
    bool node_is_split(node_index node) const {
        auto offset = routing_table.select(node,'#'_rc);
        // is not the last element of routing table and next character is not starting of new node
        return routing_table.size() != (offset + 1) and routing_table[offset+1] != '#'_rc;
    }

    int branch_starting_offset(node_index node,char branch_label) const {
        //node-1 as we are indexing from 0
        //rank1 - 1 because the rank is inclusive
        // * 'N' as we write multiple values
        int starting_offset = (joins.rank1(node-1)-1)*'N'_rc;
        int result = 0;
        assert(rc(branch_label) <= 'N'_rc); // no information for other symbols
        for(int previous_branch=0;previous_branch<rc(branch_label);previous_branch++) {
            result += edge_multiplicity_table[starting_offset+previous_branch];
        }
        return result;
    }

    int branch_size(node_index node, char branch_label) const {
        // todo: merge with branch_starting_offset
        int starting_offset = (joins.rank1(node-1)-1)*'N'_rc;
        int result = 0;
        assert(rc(branch_label) <= 'T'_rc); // no information for other symbols
        return edge_multiplicity_table[starting_offset+rc(branch_label)];
    }

    int joins_size(node_index node) const {
        // warning: correct only when N not used
        // todo: decide on adding also the last symbol
        // todo: rename
        return branch_starting_offset(node,'N');
    }

    void serialize(const fs::path& folder) const override {
        fs::create_directories(folder / "xxx.bin");
        ofstream edge_multiplicity_file(folder / "edge_multiplicity.bin", ios_base::trunc | ios_base::out);
        ofstream routing_table_file(folder / "routing_table.bin", ios_base::trunc | ios_base::out);
        ofstream joins_file(folder / "joins.bin", ios_base::trunc | ios_base::out);
        string graph_filename = folder / "graph.bin";

        ::serialize(edge_multiplicity_file,edge_multiplicity_table);
        routing_table.serialize(routing_table_file);
        joins.serialize(joins_file);
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
        ::deserialize(edge_multiplicity_file,db.edge_multiplicity_table);
        db.routing_table.load(routing_table_file);
        db.joins.load(joins_file);
        return db;
    }

    json get_statistics(unsigned int verbosity = 0) const {
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
                    int cardinality = 0;
                    for (char c : {'$','A','C','G','T','N'}) {
                        int cur = branch_starting_offset(node,c);
                        if (cur != prev) {
                            cardinality++;
                        }
                    }
                    joins_diff_symbols_histogram[cardinality]++; // size histogram doesn't have infromation whether N was present
                                                         // ToDo: fix this
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
                    int start = routing_table_offset(node);
                    int i = start;
                    for (; routing_table[i] != '#'_rc; i++) {
                        diff_symbols.insert(routing_table[i]);
                    }
                    splits_diff_symbols_histogram[diff_symbols.size()]++;
                    splits_size_histogram[i - start]++;
                }
            }
        }
        json result = {{"true_joins", true_joins},
                       {"added_joins", added_joins},
                       {"true_splits", true_splits},
                       {"added_splits", added_splits}
                      };
        if (verbosity & STATS_SPLITS_HISTOGRAM) {
            result["splits_diff_symbols_histogram"] = splits_diff_symbols_histogram;
            result["splits_size_histogram"] = splits_size_histogram;
        }
        if (verbosity & STATS_JOINS_HISTOGRAM) {
            result["splits_diff_symbols_histogram"] = splits_diff_symbols_histogram;
        }
        cerr << result.dump(4) << endl;
        return result;
    }

    path_id get_global_path_id(node_index node, int relative_position) const {
        char first_base = '\0';
        if (node_is_join(node)) {
            for(auto c = 'N'_rc; c >= 0;c--) {
                auto offset = branch_starting_offset(node,tochar(c));
                if (offset <= relative_position) {
                    first_base = tochar(c);
                    relative_position -= offset;
                    break;
                }
            }
        }
        else {
            assert(graph.indegree(node) == 1);
            // alternative for
            // graph.call_incoming_kmers(node,[&first_base](node_index node,char edge_label ) { first_base = rc(edge_label);});
            for(auto c : "ACGTN"s) {
                if (graph.traverse_back(node,c)) {
                    first_base = c;
                    break;
                }
            }

        }
        if (first_base == '$') {
            assert(is_valid_path_id({node,relative_position}));
            return {node,relative_position};
        }
        assert(first_base);
        // ToDo: get faster last character of a kmer
        auto kmer = graph.get_node_sequence(node);
        auto last_base = kmer[kmer.size()-1];
        node = graph.traverse_back(node,first_base);
        if (node_is_split(node)) {
            relative_position = routing_table_select(node,relative_position+1,last_base);// +1 as relative_position is 0-based
        }
        return get_global_path_id(node,relative_position);
    }

private:
    Wavelet routing_table;
    bit_vector_stat joins;
    vector<int> edge_multiplicity_table;
};

#endif /* path_database_baseline_hpp */
