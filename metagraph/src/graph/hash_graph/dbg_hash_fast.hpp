#ifndef __DBG_HASH_FAST_HPP__
#define __DBG_HASH_FAST_HPP__

#include <fstream>

#include "sequence_graph.hpp"
#include "kmer_extractor.hpp"


class DBGHashFast : public DeBruijnGraph {
  public:
    explicit DBGHashFast(size_t k,
                            bool canonical_mode = false,
                            bool packed_serialization = true);

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

    void call_nodes(const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &stop_early) const {
        hash_dbg_->call_nodes(callback, stop_early);
    }

    // Given a starting node, traverse the graph forward following the edge
    // sequence delimited by begin and end. Terminate the traversal if terminate()
    // returns true, or if the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    void traverse(node_index start,
                  const char* begin,
                  const char* end,
                  const std::function<void(node_index)> &callback,
                  const std::function<bool()> &terminate = [](){ return false; }) const {
        hash_dbg_->traverse(start, begin, end, callback, terminate);
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

    // Given a node index, call the target nodes of all edges outgoing from it.
    void adjacent_outgoing_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        hash_dbg_->adjacent_outgoing_nodes(node, callback);
    }

    // Given a node index, call the source nodes of all edges incoming to it.
    void adjacent_incoming_nodes(node_index node,
                                 const std::function<void(node_index)> &callback) const {
        hash_dbg_->adjacent_incoming_nodes(node, callback);
    }

    size_t outdegree(node_index node) const { return hash_dbg_->outdegree(node); }
    bool has_single_outgoing(node_index node) const { return hash_dbg_->has_single_outgoing(node); }
    bool has_multiple_outgoing(node_index node) const { return hash_dbg_->has_multiple_outgoing(node); }

    size_t indegree(node_index node) const { return hash_dbg_->indegree(node); }
    bool has_no_incoming(node_index node) const { return hash_dbg_->has_no_incoming(node); }
    bool has_single_incoming(node_index node) const { return hash_dbg_->has_single_incoming(node); }

    node_index kmer_to_node(const std::string &kmer) const {
        return hash_dbg_->kmer_to_node(kmer);
    }

    std::string get_node_sequence(node_index node) const {
        return hash_dbg_->get_node_sequence(node);
    }

    size_t get_k() const { return hash_dbg_->get_k(); }
    bool is_canonical_mode() const { return hash_dbg_->is_canonical_mode(); }

    uint64_t num_nodes() const { return hash_dbg_->num_nodes(); }
    uint64_t max_index() const { return hash_dbg_->max_index(); }

    void serialize(std::ostream &out) const { hash_dbg_->serialize(out); }
    void serialize(const std::string &filename) const { hash_dbg_->serialize(filename); }

    bool load(std::istream &in);
    bool load(const std::string &filename);

    std::string file_extension() const { return kExtension; }

    bool operator==(const DeBruijnGraph &other) const {
        if (this == &other)
            return true;

        return other == *hash_dbg_;
    }

    const std::string& alphabet() const { return hash_dbg_->alphabet(); }

    bool in_graph(node_index node) const { return hash_dbg_->in_graph(node); }

    static constexpr auto kExtension = ".hashfastdbg";

    class DBGHashFastInterface : public DeBruijnGraph {
      public:
        virtual ~DBGHashFastInterface() {}
        virtual void serialize(std::ostream &out) const = 0;
        virtual void serialize(const std::string &filename) const = 0;
        virtual bool load(std::istream &in) = 0;
        virtual bool load(const std::string &filename) = 0;
        const std::string& alphabet() const = 0;
    };

  private:
    static std::unique_ptr<DBGHashFastInterface>
    initialize_graph(size_t k, bool canonical_mode, bool packed_serialization);

    std::unique_ptr<DBGHashFastInterface> hash_dbg_;
};

#endif // __DBG_HASH_FAST_HPP__
