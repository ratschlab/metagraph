#ifndef __DBG_WRAPPER__
#define __DBG_WRAPPER__

#include "sequence_graph.hpp"


namespace mtg {
namespace graph {

/**
 * This abstract class stores a graph internally and transfers all calls to it.
 * It can be used as a parent to other wrappers which add additional functionality.
 * All children of this class must implement add_sequence, load, and operator==
 */
template <class Graph = DeBruijnGraph>
class DBGWrapper : public DeBruijnGraph {
  public:
    template <class InGraph>
    explicit DBGWrapper(std::shared_ptr<InGraph> graph)
          : graph_ptr_(std::dynamic_pointer_cast<Graph>(graph)), graph_(graph_ptr_) {
        assert(graph_);
        assert(graph_ptr_);
    }

    template <class InGraph>
    explicit DBGWrapper(std::shared_ptr<const InGraph> graph)
          : graph_(std::dynamic_pointer_cast<const Graph>(graph)) {
        assert(graph_);
        assert(!graph_ptr_);
    }

    // aliasing constructor
    template <class InGraph>
    explicit DBGWrapper(const InGraph &graph)
          : graph_(std::shared_ptr<const Graph>{}, dynamic_cast<const Graph*>(&graph)) {
        assert(graph_);
        assert(!graph_ptr_);
    }

    virtual ~DBGWrapper() {}

    /**
     * Added methods
     */
    virtual const Graph& get_graph() const { return *graph_; }
    virtual std::shared_ptr<const Graph> get_graph_ptr() const { return graph_; }
    virtual std::shared_ptr<Graph> get_mutable_graph_ptr() const { return graph_ptr_; }

    /**
     * Methods from SequenceGraph
     */
    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion
                                  = [](node_index) {}) override = 0;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override {
        graph_->map_to_nodes(sequence, callback, terminate);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override {
        graph_->map_to_nodes_sequentially(sequence, callback, terminate);
    }

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override {
        graph_->adjacent_outgoing_nodes(node, callback);
    }

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override {
        graph_->adjacent_incoming_nodes(node, callback);
    }

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override {
        graph_->call_nodes(callback, stop_early);
    }

    virtual uint64_t num_nodes() const override { return graph_->num_nodes(); }
    virtual uint64_t max_index() const override { return graph_->max_index(); }

    virtual bool load(const std::string &filename) override = 0;

    virtual void serialize(const std::string &filename) const override {
        graph_->serialize(filename);
    }

    virtual std::string file_extension() const override { return graph_->file_extension(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override {
        return graph_->get_node_sequence(index);
    }


    /**
     * Methods from DeBruijnGraph
     */
    virtual size_t get_k() const override { return graph_->get_k(); }
    virtual Mode get_mode() const override { return graph_->get_mode(); }

    virtual node_index traverse(node_index node, char next_char) const override {
        return graph_->traverse(node, next_char);
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override {
        return graph_->traverse_back(node, prev_char);
    }

    // Given a starting node and a sequence of edge labels, traverse the graph
    // forward. The traversal is terminated once terminate() returns true or
    // when the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones.
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const override {
        graph_->traverse(start, begin, end, callback, terminate);
    }

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override {
        graph_->call_sequences(callback, num_threads, kmers_in_single_form);
    }

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override {
        graph_->call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override {
        graph_->call_kmers(callback);
    }

    virtual size_t outdegree(node_index node) const override {
        return graph_->outdegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override {
        return graph_->has_single_outgoing(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override {
        return graph_->has_multiple_outgoing(node);
    }

    virtual size_t indegree(node_index node) const override {
        return graph_->indegree(node);
    }

    virtual bool has_no_incoming(node_index node) const override {
        return graph_->has_no_incoming(node);
    }

    virtual bool has_single_incoming(node_index node) const override {
        return graph_->has_single_incoming(node);
    }

    virtual node_index kmer_to_node(std::string_view kmer) const override {
        return graph_->kmer_to_node(kmer);
    }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override {
        graph_->call_outgoing_kmers(kmer, callback);
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override {
        graph_->call_incoming_kmers(kmer, callback);
    }

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override {
        return graph_->find(sequence, discovery_fraction);
    }

    virtual bool operator==(const DeBruijnGraph &other) const override = 0;

    virtual const std::string& alphabet() const override { return graph_->alphabet(); }
    virtual void print(std::ostream &out) const override { graph_->print(out); }

    // Call all nodes that have no incoming edges
    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override {
        graph_->call_source_nodes(callback);
    }

    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_->get_base_graph();
    }

  protected:
    std::shared_ptr<Graph> graph_ptr_;
    std::shared_ptr<const Graph> graph_;
};

template <typename Graph = DeBruijnGraph>
class DBGNodeModifyingWrapper : public DBGWrapper<Graph> {
  public:
    typedef typename Graph::node_index node_index;
    typedef typename Graph::CallPath CallPath;
    typedef typename Graph::OutgoingEdgeCallback OutgoingEdgeCallback;
    typedef typename Graph::IncomingEdgeCallback IncomingEdgeCallback;

    template <typename... Args>
    explicit DBGNodeModifyingWrapper(Args&&... args)
          : DBGWrapper<Graph>(std::forward<Args>(args)...) {}

    virtual ~DBGNodeModifyingWrapper() {}

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override = 0;

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override = 0;

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override = 0;

    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override = 0;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override = 0;

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override = 0;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override = 0;
    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override = 0;

    virtual uint64_t num_nodes() const override = 0;
    virtual uint64_t max_index() const override = 0;

    virtual std::string get_node_sequence(node_index index) const override = 0;

    virtual node_index traverse(node_index node, char next_char) const override = 0;
    virtual node_index traverse_back(node_index node, char prev_char) const override = 0;

    virtual size_t outdegree(node_index node) const override = 0;
    virtual size_t indegree(node_index node) const override = 0;

    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const override {
        DeBruijnGraph::traverse(start, begin, end, callback, terminate);
    }

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override {
        DeBruijnGraph::call_sequences(callback, num_threads, kmers_in_single_form);
    }

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override {
        DeBruijnGraph::call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual bool has_single_outgoing(node_index node) const override {
        return DeBruijnGraph::has_single_outgoing(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override {
        return DeBruijnGraph::has_multiple_outgoing(node);
    }

    virtual bool has_no_incoming(node_index node) const override {
        return DeBruijnGraph::has_no_incoming(node);
    }

    virtual bool has_single_incoming(node_index node) const override {
        return DeBruijnGraph::has_single_incoming(node);
    }

    virtual node_index kmer_to_node(std::string_view kmer) const override {
        return DeBruijnGraph::kmer_to_node(kmer);
    }

    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override {
        return DeBruijnGraph::find(sequence, discovery_fraction);
    }

    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override {
        return DeBruijnGraph::call_source_nodes(callback);
    }
};

} // namespace graph
} // namespace mtg

#endif // __DBG_WRAPPER__
