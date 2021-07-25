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
class IDBGWrapper : public DeBruijnGraph {
  public:
    virtual ~IDBGWrapper() {}

    /**
     * Added methods
     */
    virtual const Graph& get_graph() const { return *graph_; }
    virtual std::shared_ptr<const Graph> get_graph_ptr() const { return graph_; }
    virtual std::shared_ptr<Graph> get_mutable_graph_ptr() const { return graph_ptr_; }

    /**
     * Methods shared by all wrappers
     */
    virtual std::string file_extension() const override final { return graph_->file_extension(); }
    virtual size_t get_k() const override final { return graph_->get_k(); }
    virtual const std::string& alphabet() const override final { return graph_->alphabet(); }
    virtual void print(std::ostream &out) const override { graph_->print(out); }
    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_->get_base_graph();
    }
    virtual void serialize(const std::string &filename) const override {
        graph_->serialize(filename);
    }

    /**
     * Methods to implement
     */
    virtual bool operator==(const DeBruijnGraph &other) const override = 0;

  protected:
    std::shared_ptr<Graph> graph_ptr_;
    std::shared_ptr<const Graph> graph_;

    template <class InGraph>
    explicit IDBGWrapper(std::shared_ptr<InGraph> graph)
          : graph_ptr_(std::dynamic_pointer_cast<Graph>(graph)), graph_(graph_ptr_) {
        assert(graph_);
        assert(graph_ptr_);
    }

    template <class InGraph>
    explicit IDBGWrapper(std::shared_ptr<const InGraph> graph)
          : graph_(std::dynamic_pointer_cast<const Graph>(graph)) {
        assert(graph_);
        assert(!graph_ptr_);
    }

    // aliasing constructor
    template <class InGraph>
    explicit IDBGWrapper(const InGraph &graph)
          : graph_(std::shared_ptr<const Graph>{}, dynamic_cast<const Graph*>(&graph)) {
        assert(graph_);
        assert(!graph_ptr_);
    }
};

/**
 * This wrapper uses the default methods from Graph when available. This may be
 * used when the nodes and edges of the wrapped graph differ from the underlying
 * graph.
 */
template <typename Graph = DeBruijnGraph>
class DBGNodeModifyingWrapper : public IDBGWrapper<Graph> {
  public:
    typedef typename Graph::node_index node_index;

    template <typename... Args>
    explicit DBGNodeModifyingWrapper(Args&&... args)
          : IDBGWrapper<Graph>(std::forward<Args>(args)...) {}

    virtual ~DBGNodeModifyingWrapper() {}

    /**
     * These methods should be implemented. The Graph defaults cannot be used.
     */
    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override = 0;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override = 0;

    virtual uint64_t max_index() const override = 0;
};

/**
 * This wrapper uses the methods directly from the underlying graph by default.
 * This may be used when the nodes and edges of the wrapped graph are the
 * same as the underlying graph.
 */
template <class Graph = DeBruijnGraph>
class DBGWrapper : public IDBGWrapper<Graph> {
    typedef typename Graph::node_index node_index;
    typedef typename Graph::Mode Mode;
    typedef typename Graph::CallPath CallPath;
    typedef typename Graph::OutgoingEdgeCallback OutgoingEdgeCallback;
    typedef typename Graph::IncomingEdgeCallback IncomingEdgeCallback;

    template <typename... Args>
    explicit DBGWrapper(Args&&... args) : IDBGWrapper<Graph>(std::forward<Args>(args)...) {}

    virtual ~DBGWrapper() {}

    /**
     * Methods from SequenceGraph
     */
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override {
        this->graph_->map_to_nodes(sequence, callback, terminate);
    }

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate
                                               = [](){ return false; }) const override {
        this->graph_->map_to_nodes_sequentially(sequence, callback, terminate);
    }

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override {
        this->graph_->adjacent_outgoing_nodes(node, callback);
    }

    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override {
        this->graph_->adjacent_incoming_nodes(node, callback);
    }

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override {
        this->graph_->call_nodes(callback, stop_early);
    }

    virtual uint64_t num_nodes() const override { return this->graph_->num_nodes(); }
    virtual uint64_t max_index() const override { return this->graph_->max_index(); }

    virtual void serialize(const std::string &filename) const override {
        this->graph_->serialize(filename);
    }

    virtual std::string get_node_sequence(node_index index) const override {
        return this->graph_->get_node_sequence(index);
    }

    /**
     * Methods from DeBruijnGraph
     */
    virtual Mode get_mode() const override { return this->graph_->get_mode(); }

    virtual node_index traverse(node_index node, char next_char) const override {
        return this->graph_->traverse(node, next_char);
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override {
        return this->graph_->traverse_back(node, prev_char);
    }

    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const override {
        this->graph_->traverse(start, begin, end, callback, terminate);
    }

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override {
        this->graph_->call_sequences(callback, num_threads, kmers_in_single_form);
    }

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override {
        this->graph_->call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override {
        this->graph_->call_kmers(callback);
    }

    virtual size_t outdegree(node_index node) const override {
        return this->graph_->outdegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override {
        return this->graph_->has_single_outgoing(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override {
        return this->graph_->has_multiple_outgoing(node);
    }

    virtual size_t indegree(node_index node) const override {
        return this->graph_->indegree(node);
    }

    virtual bool has_no_incoming(node_index node) const override {
        return this->graph_->has_no_incoming(node);
    }

    virtual bool has_single_incoming(node_index node) const override {
        return this->graph_->has_single_incoming(node);
    }

    virtual node_index kmer_to_node(std::string_view kmer) const override {
        return this->graph_->kmer_to_node(kmer);
    }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override {
        this->graph_->call_outgoing_kmers(kmer, callback);
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override {
        this->graph_->call_incoming_kmers(kmer, callback);
    }

    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override {
        return this->graph_->find(sequence, discovery_fraction);
    }

    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override {
        this->graph_->call_source_nodes(callback);
    }
};

} // namespace graph
} // namespace mtg

#endif // __DBG_WRAPPER__
