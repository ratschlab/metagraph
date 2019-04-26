#ifndef __DBG_HASH_ORDERED_HPP__
#define __DBG_HASH_ORDERED_HPP__

#include <fstream>
#include <tsl/ordered_set.h>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"
#include "utils.hpp"


class DBGHashOrdered : public DeBruijnGraph {
  public:
    explicit DBGHashOrdered(size_t k, bool canonical_mode = false);

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL) {
        hash_dbg_->add_sequence(sequence, nodes_inserted);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes(sequence, callback, terminate);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->map_to_nodes_sequentially(begin, end, callback, terminate);
    }

    void call_outgoing_kmers(node_index node,
                             const OutgoingEdgeCallback &callback) const {
        hash_dbg_->call_outgoing_kmers(node, callback);
    }

    void call_incoming_kmers(node_index node,
                             const IncomingEdgeCallback &callback) const {
        hash_dbg_->call_incoming_kmers(node, callback);
    }

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const {
        return hash_dbg_->traverse(node, next_char);
    }

    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const {
        return hash_dbg_->traverse_back(node, prev_char);
    }

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const {
        hash_dbg_->adjacent_outgoing_nodes(node, target_nodes);
    }

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const {
        hash_dbg_->adjacent_incoming_nodes(node, source_nodes);
    }


    size_t outdegree(node_index node) const { return hash_dbg_->outdegree(node); }
    size_t indegree(node_index node) const { return hash_dbg_->indegree(node); }

    node_index kmer_to_node(const std::string &kmer) const {
        return hash_dbg_->kmer_to_node(kmer);
    }

    std::string get_node_sequence(node_index node) const {
        return hash_dbg_->get_node_sequence(node);
    }

    size_t get_k() const { return hash_dbg_->get_k(); }
    bool is_canonical_mode() const { return hash_dbg_->is_canonical_mode(); }

    uint64_t num_nodes() const { return hash_dbg_->num_nodes(); }

    void serialize(std::ostream &out) const { hash_dbg_->serialize(out); }
    void serialize(const std::string &filename) const { hash_dbg_->serialize(filename); }

    bool load(std::istream &in);
    bool load(const std::string &filename);

    bool operator==(const DeBruijnGraph &other) const {
        if (!dynamic_cast<const DBGHashOrdered*>(&other)) {
            throw std::runtime_error("Not implemented");
            return false;
        }

        return *hash_dbg_ == *dynamic_cast<const DBGHashOrdered*>(&other)->hash_dbg_;
    }

    static constexpr auto kExtension = ".orhashdbg";

    class DBGHashOrderedInterface : public DeBruijnGraph {
      public:
        virtual ~DBGHashOrderedInterface() {}
        virtual void serialize(std::ostream &out) const = 0;
        virtual void serialize(const std::string &filename) const = 0;
        virtual bool load(std::istream &in) = 0;
        virtual bool load(const std::string &filename) = 0;
    };

  private:
    static std::unique_ptr<DBGHashOrderedInterface>
    initialize_graph(size_t k, bool canonical_mode);

    std::unique_ptr<DBGHashOrderedInterface> hash_dbg_;
};

#endif // __DBG_HASH_ORDERED_HPP__
