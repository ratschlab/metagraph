//
// Created by Jan Studený on 2019-03-08.
//

#ifndef __PATH_DATABASE_LIST_OF_BIFURCATION_CHOICES_HPP__
#define __PATH_DATABASE_LIST_OF_BIFURCATION_CHOICES_HPP__


#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <vector>
#include <nlohmann/json.hpp>

#include "dbg_succinct.hpp"
#include "path_database.hpp"
#include "utilities.hpp"

#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wreturn-type"




using namespace std;
using namespace std::string_literals;
using node_index = SequenceGraph::node_index;
const int64_t DEFAULT_K_KMER = 21;

// Stores reads in a compressed format
class PathDatabaseListBC : public PathDatabase<int64_t> {
  public:
    // const value
    // PathDatabaseListBC(const vector<string> &raw_reads)
    // non-const pointer to modify
    // PathDatabaseListBC(vector<string> *raw_reads)
    //
    // Graph |graph| must contain all k-mers from the sequences passed
    PathDatabaseListBC(DBGSuccinct *graph,
                    const vector<string> &raw_reads,
                    size_t kmer_length = DEFAULT_K_KMER)
          : PathDatabase(std::shared_ptr<const DBGSuccinct> { graph }),
            kmer_length_(kmer_length),
            read_length(raw_reads[0].length()) {
    }

    PathDatabaseListBC(const vector<string> &raw_reads,
                    size_t kmer_length = DEFAULT_K_KMER)
          : PathDatabase(raw_reads,kmer_length), kmer_length_(kmer_length),read_length(raw_reads[0].length()) {}

    std::vector<string> get_all_reads() const {
        vector<string> reads;
        for (const auto &compressed_read : compressed_reads_) {
            reads.push_back(decode_read(compressed_read));
        }
        return reads;
    }

    size_t num_paths() const override { return compressed_reads_.size(); }

    virtual json get_statistics(__attribute__((unused)) int64_t verbosity = 0) const override {
        std::map<std::string, int64_t> bifurcation_size_histogram;
        for (const auto &read : compressed_reads_) {
            bifurcation_size_histogram[to_string(read.second.size())]++;
        }
        json result = {{"bifurcation_histogram",bifurcation_size_histogram},
                       {"total_size",compressed_size_without_reference()},
                       {"bifurcation_size", compressed_size_without_reference()-2*num_paths()*graph_->get_k()},
                       {"number_of_reads",compressed_reads_.size()}};
        cerr << result.dump(4) << endl;
        return result;
    }

    int64_t compressed_size_without_reference() const {
        // returns size in bits
        int64_t size = 0;
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
            compressed_reads_.push_back(encode_read(read));
            ids.push_back(ids.size());
        }
        return ids;
    }

    // returns ids of all paths that go through sequence |str|
    std::vector<path_id> get_paths_going_through(const std::string &) const override {}
    std::vector<path_id> get_paths_going_through(node_index) const override {}

    // make one traversal step through the selected path
    node_index get_next_node(node_index, path_id) const override {}

    // transition to the next node consistent with the history
    // return npos if there is no transition consistent with the history
    node_index get_next_consistent_node(const std::string &) const override {}

    std::string decode(path_id path) const override {
        return decode_read(compressed_reads_[path]);
    }

  private:
    using kmer_t = std::string;
    using edge_choices_t = std::vector<char>;
    using compressed_read_t = std::pair<kmer_t, edge_choices_t>;

    compressed_read_t encode_read(const string &read) const {
        auto kmer = read.substr(0, kmer_length_);
        auto node = graph_->kmer_to_node(kmer);

        edge_choices_t edge_choices;

        // for all other characters
        for (char character : read.substr(kmer_length_)) {
            int outdegree = graph_->outdegree(node);
            if (outdegree > 1) {
                edge_choices.push_back(character);
            }
            node = graph_->traverse(node, character);
            assert(outdegree > 0);
        }
        return { kmer, edge_choices };
    }

    std::string decode_read(const compressed_read_t &compressed_read) const {
        auto& [starting_kmer, edge_choices] = compressed_read;
        auto current_bifurcation_choice = edge_choices.begin();
        auto node = graph_->kmer_to_node(starting_kmer);
        string read = starting_kmer;
        char next_char;
        while ((int64_t )read.size() < read_length) {
            int64_t outgoing_degree = 0;
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

    void serialize(const fs::path&) const override {};


    std::vector<compressed_read_t> compressed_reads_;
    const int64_t kmer_length_;
    const int64_t read_length;
};

#endif // __PATH_DATABASE_LIST_OF_BIFURCATION_CHOICES_HPP__
