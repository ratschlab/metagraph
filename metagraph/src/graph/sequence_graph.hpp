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

    // Maps k-mers from sequence to nodes of the graph similarly to map_to_nodes
    // Guarantees that the k-mers from sequence are called in their natural order
    virtual void sequence_map(const std::string::const_iterator &begin,
                              const std::string::const_iterator &end,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    virtual node_index kmer_to_node(const char *begin) const {
        node_index node = npos;
        map_to_nodes(std::string(begin, get_k()),
            [&node](node_index i) { node = i; },
        );
        return result;
    }

    virtual node_index kmer_to_node(const std::string &kmer) const {
        assert(kmer.size() == get_k());
        return kmer_to_node(kmer.data());
    }

    using OutgoingEdgeCallback = std::function<void(node_index /* target_kmer */,
                                                    char /* last_target_char */)>;
    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const = 0;

    using IncomingEdgeCallback = std::function<void(node_index /* source_kmer */,
                                                    char /* first_source_char */)>;
    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const = 0;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(const std::string &sequence,
                      double discovery_fraction = 1) const {
        if (sequence.length() < get_k())
            return false;

        const size_t num_kmers = sequence.length() - get_k() + 1;
        const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
        const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
        size_t num_kmers_discovered = 0;
        size_t num_kmers_missing = 0;

        map_to_nodes(sequence,
            [&](node_index node) {
                if (node) {
                    num_kmers_discovered++;
                } else {
                    num_kmers_missing++;
                }
            },
            [&]() { return num_kmers_missing > max_kmers_missing
                            || num_kmers_discovered >= min_kmers_discovered; }
        );

        return num_kmers_missing <= max_kmers_missing;
    }
};


#endif // __SEQUENCE_GRAPH_HPP__
