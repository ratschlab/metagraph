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

#include "path_encoder.hpp"
#include "utilities.hpp"

using json = nlohmann::json;


using namespace std;
using namespace std::string_literals;
using node_index = SequenceGraph::node_index;
const int DEFAULT_K_KMER = 21;


// Stores reads in a compressed format
class CompressedReads : public PathDatabase {
  public:
    // TODO: read about lvalues, rvalues, etc
    // http://thbecker.net/articles/rvalue_references/section_01.html
    // const value
    // CompressedReads(const vector<string> &raw_reads)
    // non-const pointer to modify
    // CompressedReads(vector<string> *raw_reads)
    //
    // Graph |graph| must contain all k-mers from the sequences passed
    CompressedReads(DBGSuccinct *graph,
                    const vector<string> &raw_reads,
                    size_t k_kmer = DEFAULT_K_KMER)
          : PathDatabase(std::shared_ptr<const DeBruijnGraph> { graph }),
            k_kmer_(k_kmer),
            read_length(raw_reads[0].length()) {
        for (const auto &read : raw_reads) {
            compressed_reads_.insert(encode_read(read));
        }
    }

    CompressedReads(const vector<string> &raw_reads,
                    size_t k_kmer = DEFAULT_K_KMER)
          : CompressedReads(new DBGSuccinct(dbg_succ_graph_constructor(raw_reads, k_kmer)),
                            raw_reads,
                            k_kmer) {}

    std::vector<string> get_reads() const {
        vector<string> reads;
        for (const auto &compressed_read : compressed_reads_) {
            reads.push_back(decode_read(compressed_read));
        }
        return reads;
    }

    size_t num_paths() const { return compressed_reads_.size(); }

    json get_statistics() const {
        std::map<std::string, int> bifurcation_size_histogram;
        for (const auto &read : compressed_reads_) {
            bifurcation_size_histogram[to_string(read.second.size())]++;
        }
        json result = {{"bifurcation_size",bifurcation_size_histogram},
                       {"total_size",compressed_size_without_reference()},
                       {"number_of_reads",compressed_reads_.size()}};
        cerr << result.dump(4) << endl;
        return result;
    }

    int compressed_size_without_reference() const {
        // returns size in bits
        int size = 0;
        for(auto &read : compressed_reads_) {
            // *2 for two bit encoding
            size += read.first.size() * 2;
            size += read.second.size() * 2;
        }
        // to be really fair
        // size += sizeof(read_length);
        return size;
    }

    double bits_per_symbol(bool include_reference = false) const {
        assert(!include_reference);
        return static_cast<double>(compressed_size_without_reference())
               / (compressed_reads_.size()*read_length);
    }

    std::vector<path_id> encode(const std::vector<std::string> &sequences) override {
        std::vector<path_id> ids;
        for (const auto &read : sequences) {
            compressed_reads_.insert(encode_read(read));
            ids.push_back(ids.size());
        }
        return ids;
    }

    node_index get_first_node(path_id path) const {}
    node_index get_last_node(path_id path) const {}

    // returns ids of all paths that go through sequence |str|
    std::vector<path_id> get_paths_going_through(const std::string &str) const {}
    std::vector<path_id> get_paths_going_through(node_index node) const {}

    // make one traversal step through the selected path
    node_index get_next_node(node_index node, path_id path) const {}

    // transition to the next node consistent with the history
    // return npos if there is no transition consistent with the history
    node_index get_next_consistent_node(const std::string &history) const {}

    std::string decode(path_id path) const {
        auto it = compressed_reads_.begin();
        for (size_t i = 0; i < path; ++i) {
            ++it;
        }
        return decode_read(*it);
    }

  private:
    using kmer_t = std::string;
    using edge_choices_t = std::vector<char>;
    using compressed_read_t = std::pair<kmer_t, edge_choices_t>;

    compressed_read_t encode_read(const string &read) const {
        auto kmer = read.substr(0, k_kmer_);
        auto node = graph_->kmer_to_node(kmer);

        edge_choices_t edge_choices;

        // for all other characters
        for (char character : read.substr(k_kmer_)) {
            vector<node_index> outnodes;
            graph_->adjacent_outgoing_nodes(node, &outnodes);
            if (outnodes.size() > 1) {
                edge_choices.push_back(character);
            }
            node = graph_->traverse(node, character);
            assert(!outnodes.empty());
        }
        return { kmer, edge_choices };
    }

    std::string decode_read(const compressed_read_t &compressed_read) const {
        auto& [starting_kmer, edge_choices] = compressed_read;
        auto current_bifurcation_choice = edge_choices.begin();
        auto node = graph_->kmer_to_node(starting_kmer);
        string read = starting_kmer;
        char next_char;
        while (read.size() < read_length) {
            int outgoing_degree = 0;
            graph_->call_outgoing_kmers(node,[&](node_index next_node, char character) {
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
            if (outgoing_degree > 1) {
                current_bifurcation_choice++;//prepare next choice
            }
            read += next_char;
        }
        return read;
    }

    static DBG_succ* dbg_succ_graph_constructor(const vector<string> &raw_reads,
                                                size_t k_kmer) {
        auto graph_constructor = DBGSuccConstructor(k_kmer - 1);// because DBG_succ has smaller kmers
        cerr << "Starting building the graph" << endl;
        for(auto &read : raw_reads) {
            assert(read.size() >= k_kmer);
            graph_constructor.add_sequence(read);
        }
        return new DBG_succ(&graph_constructor);
        //graph = DBGSuccinct(stupid_old_representation);
    }


    std::multiset<compressed_read_t> compressed_reads_;
    const int read_length;
    const int k_kmer_;
};

#endif //METAGRAPH_COMPRESSED_READS_HPP
