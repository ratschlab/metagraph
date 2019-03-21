#ifndef __DBG_HASH_HPP__
#define __DBG_HASH_HPP__

#include <fstream>
#include <tsl/hopscotch_map.h>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"


class DBGHash : public DeBruijnGraph {
  public:
    explicit DBGHash(size_t k) : k_(k) {}

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Map k-mers from sequence to nodes of the graph similarly to map_to_nodes
    // Guarantees that the k-mers from sequence are called in their natural order
    void map_kmers_sequentially(std::string::const_iterator,
                                std::string::const_iterator,
                                const std::function<void(node_index)> &,
                                const std::function<bool()> &) const {
        // TODO: Complete map_sequence_sequentially for DBGHash.
        throw std::runtime_error("Not implemented");
    }

    void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const {
        // TODO: Complete call_outgoing_kmers for DBGHash.
        throw std::runtime_error("Not implemented");
    }

    void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const {
        // TODO: Complete call_incoming_kmers for DBGHash.
        throw std::runtime_error("Not implemented");
    }

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const;

    size_t outdegree(node_index node) const;

    node_index kmer_to_node(const std::string &kmer) const;
    std::string node_to_kmer(node_index node) const;

    std::string get_node_sequence(node_index node) const {
        return node_to_kmer(node);
    }

    size_t get_k() const { return k_; }
    uint64_t num_nodes() const { return kmers_.size(); }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

  private:
    std::string encode_sequence(const std::string &sequence) const;

    size_t k_;
    tsl::hopscotch_map<std::string, uint64_t> indices_;
    std::vector<std::string> kmers_;
    KmerExtractor2Bit seq_encoder_;

    static constexpr auto kExtension = ".hashdbg";
};

#endif // __DBG_HASH_HPP__
