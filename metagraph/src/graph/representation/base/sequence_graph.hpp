#ifndef __SEQUENCE_GRAPH_HPP__
#define __SEQUENCE_GRAPH_HPP__

#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <memory>


namespace utils {
    std::string make_suffix(const std::string &str, const std::string &suffix);
} // namespace utils


namespace mtg {
namespace graph {

class SequenceGraph {
  public:
    // Node indexes [1,...,max_index], but only num_nodes of them are real.
    // For iteration, call `call_nodes`
    typedef uint64_t node_index;
    static constexpr uint64_t npos = 0;

    virtual ~SequenceGraph() {}

    // Insert sequence to graph and invoke callback |on_insertion| for each new
    // node index augmenting the range [1,...,max_index], including those not
    // pointing to any real node in graph. That is, the callback is invoked for
    // all new real nodes and all new dummy node indexes allocated in graph.
    // In short: max_index[after] = max_index[before] + {num_invocations}.
    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion = [](uint64_t) {}) = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const = 0;
    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const = 0;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const;

    virtual uint64_t num_nodes() const = 0;
    virtual uint64_t max_index() const { return num_nodes(); };

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;
    virtual std::string file_extension() const = 0;

    // Get string corresponding to |node|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index node) const = 0;

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

        if (!extension->load(utils::make_suffix(filename, file_extension())))
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
    enum Mode { BASIC = 0, CANONICAL, PRIMARY };

    virtual ~DeBruijnGraph() {}

    virtual size_t get_k() const = 0;

    virtual Mode get_mode() const = 0;

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const = 0;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are not mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Given a starting node and a sequence of edge labels, traverse the graph
    // forward. The traversal is terminated once terminate() returns true or
    // when the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones.
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const;

    typedef std::function<void(const std::string&, const std::vector<node_index>&)> CallPath;

    /**
     * Call contigs, a set of sequences covering each node in the graph exactly once.
     * @param num_threads number of threads to use for graph traversal
     * @param kmers_in_single_form if true, output each k-mer only in one of its forms
     * (canonical/non-canonical). That is, skip a k-mer if its reverse-complement has been
     * extracted.
     */
    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const;
    /**
     * Call all unitigs except short tips, where tips are
     * the unitigs with InDegree(first) + OutDegree(last) < 2.
     * If |kmers_in_single_form| is true, output each k-mer only in one of its
     * forms (canonical/non-canonical). That is, skip a k-mer if its
     * reverse-complement has been extracted.
     */
    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const;

    virtual size_t outdegree(node_index) const = 0;
    virtual bool has_single_outgoing(node_index node) const { return outdegree(node) == 1; }
    virtual bool has_multiple_outgoing(node_index node) const { return outdegree(node) > 1; }

    virtual size_t indegree(node_index) const = 0;
    virtual bool has_no_incoming(node_index node) const { return indegree(node) == 0; }
    virtual bool has_single_incoming(node_index node) const { return indegree(node) == 1; }

    virtual node_index kmer_to_node(std::string_view kmer) const;

    using OutgoingEdgeCallback = std::function<void(node_index /* target_kmer */,
                                                    char /* last_target_char */)>;
    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const = 0;

    using IncomingEdgeCallback = std::function<void(node_index /* source_kmer */,
                                                    char /* first_source_char */)>;
    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const = 0;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const;

    virtual bool operator==(const DeBruijnGraph &other) const;
    virtual bool operator!=(const DeBruijnGraph &other) const { return !operator==(other); }

    virtual const std::string& alphabet() const = 0;

    virtual void print(std::ostream &out) const;

    friend std::ostream& operator<<(std::ostream &out, const DeBruijnGraph &graph);

    // Call all nodes that have no incoming edges
    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const;

    virtual const DeBruijnGraph& get_base_graph() const { return *this; }
};


// returns the edge rank, starting from zero
size_t incoming_edge_rank(const SequenceGraph &graph,
                          SequenceGraph::node_index source,
                          SequenceGraph::node_index target);

std::vector<SequenceGraph::node_index>
map_sequence_to_nodes(const SequenceGraph &graph, std::string_view sequence);

void reverse_complement_seq_path(const SequenceGraph &graph,
                                 std::string &seq,
                                 std::vector<SequenceGraph::node_index> &path);

} // namespace graph
} // namespace mtg

#endif // __SEQUENCE_GRAPH_HPP__
