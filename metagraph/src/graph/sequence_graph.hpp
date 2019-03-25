#ifndef __SEQUENCE_GRAPH_HPP__
#define __SEQUENCE_GRAPH_HPP__

#include <vector>
#include <string>
#include <functional>

class bit_vector_dyn;


class SequenceGraph {
  public:
    // node indexes [1,...,num_nodes]
    typedef uint64_t node_index;
    static constexpr uint64_t npos = 0;

    virtual ~SequenceGraph() {}

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    virtual void add_sequence(const std::string &sequence,
                              bit_vector_dyn *nodes_inserted = NULL) = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         std::vector<node_index> *target_nodes) const = 0;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    virtual void adjacent_incoming_nodes(node_index node,
                                         std::vector<node_index> *source_nodes) const = 0;

    virtual uint64_t num_nodes() const = 0;

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index node_index) const = 0;
};


class DeBruijnGraph : public SequenceGraph {
  public:
    virtual ~DeBruijnGraph() {}

    virtual size_t get_k() const = 0;

    virtual bool is_canonical_mode() const { return false; }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const = 0;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are not mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                                        = [](){ return false; }) const = 0;

    virtual size_t outdegree(node_index) const = 0;
    virtual size_t indegree(node_index) const = 0;

    virtual node_index kmer_to_node(const char *begin) const;
    virtual node_index kmer_to_node(const std::string &kmer) const;

    using OutgoingEdgeCallback = std::function<void(node_index /* target_kmer */,
                                                    char /* last_target_char */)>;
    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const = 0;

    using IncomingEdgeCallback = std::function<void(node_index /* source_kmer */,
                                                    char /* first_source_char */)>;
    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const = 0;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(const std::string &sequence, double discovery_fraction = 1) const;
};

#endif // __SEQUENCE_GRAPH_HPP__
