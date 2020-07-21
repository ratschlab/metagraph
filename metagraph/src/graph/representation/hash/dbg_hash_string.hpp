#ifndef __DBG_HASH_STRING_HPP__
#define __DBG_HASH_STRING_HPP__

#include <iostream>
#include <tsl/ordered_set.h>

#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class DBGHashString : public DeBruijnGraph {
  public:
    explicit DBGHashString(size_t k) : k_(k) {}

    // Insert sequence to graph and invoke callback |on_insertion| for each new
    // node created in the graph.
    void add_sequence(std::string_view sequence,
                      const std::function<void(node_index)> &on_insertion = [](node_index) {});

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string_view sequence,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(std::string_view sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    void call_outgoing_kmers(node_index node, const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index node, const IncomingEdgeCallback &callback) const;

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;
    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const;

    size_t outdegree(node_index) const;
    bool has_single_outgoing(node_index) const;
    bool has_multiple_outgoing(node_index) const;

    size_t indegree(node_index) const;
    bool has_no_incoming(node_index) const;
    bool has_single_incoming(node_index) const;

    void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const;

    node_index kmer_to_node(std::string_view kmer) const;

    std::string get_node_sequence(node_index node) const;

    size_t get_k() const { return k_; }
    uint64_t num_nodes() const { return kmers_.size(); }

    // TODO: add the support for the canonical mode
    bool is_canonical_mode() const { return false; }

    void serialize(std::ostream &out) const;
    void serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const;

    const std::string& alphabet() const;

    static constexpr auto kExtension = ".hashstrdbg";

  private:
    std::vector<std::string> encode_sequence(std::string_view sequence) const;

    size_t k_;

    using KmerIndex = tsl::ordered_set<std::string,
                                       std::hash<std::string>,
                                       std::equal_to<std::string>,
                                       std::allocator<std::string>,
                                       std::deque<std::string>,
                                       std::uint64_t>;
    KmerIndex kmers_;

    static const std::string alphabet_;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_HASH_STRING_HPP__
