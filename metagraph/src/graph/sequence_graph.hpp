#ifndef __SEQUENCE_GRAPH_HPP__
#define __SEQUENCE_GRAPH_HPP__

#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <memory>

class bit_vector_dyn;

namespace utils {
    std::string remove_suffix(const std::string &str, const std::string &suffix);
}


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

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const = 0;
    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const = 0;

    virtual uint64_t num_nodes() const = 0;

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;
    virtual std::string file_extension() const = 0;

    // Get string corresponding to |node|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index node) const = 0;

    // Check if the node index is a valid node in the graph
    virtual bool in_graph(node_index node) const = 0;

    /********************************************************/
    /******************* graph extensions *******************/
    /********************************************************/

    class GraphExtension {
      public:
        virtual ~GraphExtension() {}
        virtual bool load(const std::string &filename_base) = 0;
        virtual void serialize(const std::string &filename_base) const = 0;
        virtual bool is_compatible(const SequenceGraph &graph, bool verbose = true) const = 0;
    };

    // TODO: improve interface: either prohibit or support
    //       properly multiple extensions of the same type.
    //       Use GraphExtension::file_extension() or enum.
    void add_extension(std::shared_ptr<GraphExtension> extension);

    template <class ExtensionSubtype>
    std::shared_ptr<ExtensionSubtype> get_extension() const {
        static_assert(std::is_base_of<GraphExtension, ExtensionSubtype>::value);
        for (auto extension : extensions_) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(extension))
                return match;
        }
        return nullptr;
    }

    template <class ExtensionSubtype>
    void remove_extension() {
        static_assert(std::is_base_of<GraphExtension, ExtensionSubtype>::value);
        for (auto it = extensions_.begin(); it != extensions_.end(); ++it) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(*it)) {
                extensions_.erase(it);
                return;
            }
        }
    }

    template <class ExtensionSubtype>
    void for_each(std::function<void(ExtensionSubtype &extension)> callback) {
        static_assert(std::is_base_of<GraphExtension, ExtensionSubtype>::value);
        for (auto extension : extensions_) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(extension))
                callback(*match);
        }
    };

    template <class ExtensionSubtype>
    std::shared_ptr<ExtensionSubtype> load_extension(const std::string &filename) {
        static_assert(std::is_base_of<GraphExtension, ExtensionSubtype>::value);
        remove_extension<ExtensionSubtype>();
        auto extension = std::make_shared<ExtensionSubtype>();

        auto filename_base = utils::remove_suffix(filename, file_extension());

        if (!extension->load(filename_base + file_extension()))
            return nullptr;

        add_extension(extension);
        return extension;
    }

    void serialize_extensions(const std::string &filename_base) const;

  private:
    std::vector<std::shared_ptr<GraphExtension>> extensions_;
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
                                           const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Given a starting node, traverse the graph forward following the edge
    // sequence delimited by begin and end. Terminate the traversal if terminate()
    // returns true, or if the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void traverse(node_index start,
                          const char* begin,
                          const char* end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const;

    // TODO: move to graph_algorithm.hpp
    virtual void call_sequences(const std::function<void(const std::string&)> &callback) const;
    /**
     * Call all unitigs except short tips.
     * Def. Tips are the unitigs with InDegree(first) + OutDegree(last) < 2.
     */
    virtual void call_unitigs(const std::function<void(const std::string&)> &callback,
                              size_t min_tip_size = 1) const;
    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const;
    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const;

    virtual size_t outdegree(node_index) const = 0;
    virtual bool has_single_outgoing(node_index node) const { return outdegree(node) == 1; }
    virtual bool has_multiple_outgoing(node_index node) const { return outdegree(node) > 1; }

    virtual size_t indegree(node_index) const = 0;
    virtual bool has_no_incoming(node_index node) const { return indegree(node) == 0; }
    virtual bool has_single_incoming(node_index node) const { return indegree(node) == 1; }

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

    virtual bool operator==(const DeBruijnGraph &other) const;
    virtual bool operator!=(const DeBruijnGraph &other) const { return !operator==(other); }

    virtual const std::string& alphabet() const = 0;

    virtual void print(std::ostream &out) const;

    friend std::ostream& operator<<(std::ostream &out, const DeBruijnGraph &graph);

    // Call all nodes that have no incoming edges
    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const;

    // Check if the node index is a valid node in the graph
    virtual bool in_graph(node_index node) const = 0;
};


// returns the edge rank, starting from zero
size_t incoming_edge_rank(const DeBruijnGraph &graph,
                          DeBruijnGraph::node_index source,
                          DeBruijnGraph::node_index target);

#endif // __SEQUENCE_GRAPH_HPP__
