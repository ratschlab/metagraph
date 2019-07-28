#ifndef __PATH_ENCODER_HPP__
#define __PATH_ENCODER_HPP__

//
// Created by Jan Studený on 2019-03-08.
//


#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <vector>
#include <nlohmann/json.hpp>

#include "sequence_graph.hpp"
#include "boss_construct.hpp"
#include "utilities.hpp"
#include "configuration.hpp"

#include "graph_patch.hpp"

using json = nlohmann::json;

template <typename PathID,
          typename GraphT = DBGSuccinct>
class PathDatabase {
    // This is an abstract class, for implementation, use PathDatabaseWavelet.
  public:
    // convenience constructor
    explicit PathDatabase(const vector<string> &reads, size_t kmer_length = 21) {
        Timer timer;
        cerr << "Started building the graph" << endl;
        auto graph = std::make_unique<GraphT>(dbg_succ_graph_constructor(reads, kmer_length));
#ifdef MASK_DUMMY_KMERS
        graph->mask_dummy_kmers(1, false);
#endif
        graph_.reset(graph.release());
        auto elapsed = timer.elapsed();
        cerr << "Building finished in " << elapsed << " sec." << endl;
        statistics["graph_build_time"] = elapsed;
    }

    PathDatabase(std::shared_ptr<const GraphT> graph) : graph_(graph) {}

    virtual ~PathDatabase() = default;

    using node_index = DeBruijnGraph::node_index;
    using path_id = PathID;

    // compress a batch of sequences
    virtual std::vector<path_id>
    encode(const std::vector<std::string> &sequences) = 0;

    // // All k-mers from the sequence must be represented in graph
    // path_id encode(const std::string &sequence);
    // void compress();

    virtual size_t num_paths() const = 0;

    virtual std::string decode(path_id path) const = 0;

    virtual node_index get_first_node(path_id path) const {
        auto decoded_path = decode(path);
        return graph_->kmer_to_node(decoded_path.substr(0,graph_->get_k()));
    };

    // returns ids of all paths that go through sequence |str|
    virtual std::vector<path_id> get_paths_going_through(const std::string &str) const = 0;
    virtual std::vector<path_id> get_paths_going_through(node_index node) const = 0;

    // make one traversal step through the selected path
    // TODO: Figure out what to do if the node is visited multiple times in the path.
    virtual node_index get_next_node(node_index node, path_id path) const = 0;

    // transition to the next node consistent with the history
    // return npos if there is no transition consistent with the history
    virtual node_index get_next_consistent_node(const std::string &history) const = 0;

    virtual void serialize(const fs::path& folder) const = 0;

    virtual json get_statistics(int64_t verbosity = 0) const {
        return statistics;
    }

  protected:
    std::shared_ptr<const GraphT> graph_;

    json statistics;
};


#endif // __PATH_ENCODER_HPP__
