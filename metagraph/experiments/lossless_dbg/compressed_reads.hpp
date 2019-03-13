//
// Created by Jan Studený on 2019-03-08.
//

#ifndef METAGRAPH_COMPRESSED_READS_HPP
#define METAGRAPH_COMPRESSED_READS_HPP


#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <vector>
#include <nlohmann/json.hpp>
#include "utilities.hpp"

using json = nlohmann::json;


using namespace std;
using namespace std::string_literals;
using node_index = SequenceGraph::node_index;
const int DEFAULT_K_KMER = 21;


class CompressedReads {
    using bifurcation_choices_t = vector<char>;
    using kmer_t = string;
    using compressed_read_t = pair<kmer_t, bifurcation_choices_t>;

public:
    // TODO: read about lvalues, rvalues, etc
    // http://thbecker.net/articles/rvalue_references/section_01.html
    // const value
    // CompressedReads(const vector<string> &raw_reads)
    // non-const pointer to modify
    // CompressedReads(vector<string> *raw_reads)
    CompressedReads(const vector<string>& raw_reads, int k_kmer = DEFAULT_K_KMER)
            : k_kmer(k_kmer),
              read_length(raw_reads[0].length()),
              graph(dbg_succ_graph_constructor(raw_reads,k_kmer)) {
        for(auto & read : raw_reads) {
            compressed_reads.insert(align_read(read));
        }
    }

    int compressed_size_without_reference() {
        // returns size in bits
        int size = 0;
        for(auto &read : compressed_reads) {
            // *2 for two bit encoding
            size += read.first.size() * 2;
            size += read.second.size() * 2;
        }
        // to be really fair
        // size += sizeof(read_length);
        return size;
    }

    double bits_per_symbol(bool include_reference = false) {
        assert(!include_reference);
        return static_cast<double>(compressed_size_without_reference())
               / (compressed_reads.size()*read_length);
    }

    vector<string> get_reads() {
        vector<string> reads;
        for(auto& compressed_read : compressed_reads) {
            reads.push_back(decompress_read(compressed_read));
        }
        return reads;
    }

    json get_statistics() {
        map<string,int> bifurcation_size_histogram;
        for(auto& read : compressed_reads) {
            bifurcation_size_histogram[to_string(read.second.size())]++;
        }
        json result = {{"bifurcation_size",bifurcation_size_histogram},
                       {"total_size",compressed_size_without_reference()},
                       {"number_of_reads",compressed_reads.size()}};
        cerr << result.dump(4) << endl;
        return result;
    }

private:
    compressed_read_t align_read(const string &read) {
        auto kmer = read.substr(0,k_kmer);
        auto bifurcation_choices = bifurcation_choices_t();
        auto node = graph.kmer_to_node(kmer);
        // for all other characters
        for (auto &character : read.substr(k_kmer)) {
            vector<node_index> outnodes;
            graph.adjacent_outgoing_nodes(node, &outnodes);
            if (outnodes.size() > 1) {
                bifurcation_choices.push_back(character);
            }
            node = graph.traverse(node,character);
            assert(!outnodes.empty());
        }
        return {kmer,bifurcation_choices};
    }
    static DBG_succ* dbg_succ_graph_constructor(const vector<string> &raw_reads, int k_kmer) {
        auto graph_constructor = DBGSuccConstructor(k_kmer - 1);// because DBG_succ has smaller kmers
        cerr << "Starting building the graph" << endl;
        for(auto &read : raw_reads) {
            assert(read.size() >= k_kmer);
            graph_constructor.add_sequence(read);
        }
        return new DBG_succ(&graph_constructor);
        //graph = DBGSuccinct(stupid_old_representation);
    }

    string decompress_read(const compressed_read_t &compressed_read) {
        auto& [starting_kmer,bifurcation_choices] = compressed_read;
        auto current_bifurcation_choice = bifurcation_choices.begin();
        auto node = graph.kmer_to_node(starting_kmer);
        string read = starting_kmer;
        char next_char;
        while(read.size() < read_length) {
            int outgoing_degree = 0;
            graph.call_outgoing_kmers(node,[&](node_index next_node, char character) {
                outgoing_degree++;
                if (outgoing_degree>1) {
                    if (character == *current_bifurcation_choice) {
                        //initial guess was wrong
                        node = next_node;
                        next_char = character;
                    }
                }
                else {
                    //initial guess
                    node = next_node;
                    next_char = character;
                }
            });
            if (outgoing_degree>1) {
                current_bifurcation_choice++;//prepare next choice
            }
            read += next_char;
        }
        return read;
    }


    multiset<compressed_read_t> compressed_reads;
    int read_length;
    int k_kmer;
    static const int defult_k_k_mer = 21;
    DBGSuccinct graph;
};

#endif //METAGRAPH_COMPRESSED_READS_HPP
