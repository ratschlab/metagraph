#ifndef __DBG_WRAPPER__
#define __DBG_WRAPPER__

#include "sequence_graph.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {

/**
 * This class stores a graph internally and transfers all calls to it.
 * It can be used as a parent to other wrappers which add additional functionality.
 */
template <class Graph = DeBruijnGraph>
class DBGWrapper : public DeBruijnGraph {
  public:
    explicit DBGWrapper(std::shared_ptr<const Graph> graph) : graph_(graph) {
        assert(graph_);
    }

    explicit DBGWrapper(std::shared_ptr<Graph> graph) : graph_(graph), graph_ptr_(graph) {
        assert(graph_);
        assert(graph_ptr_);
    }

    virtual ~DBGWrapper() {}

    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_->get_base_graph();
    }

    virtual const Graph& get_graph() const { return *graph_; }
    virtual std::shared_ptr<const Graph> get_graph_ptr() const { return graph_; }
    virtual std::shared_ptr<Graph> get_mutable_graph_ptr() const { return graph_ptr_; }

    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion
                                  = [](node_index) {}) override {
        if (!graph_ptr_) {
            common::logger->error("add_sequence only supported for non-const graphs");
            exit(1);
        }

        graph_ptr_->add_sequence(sequence, on_insertion);
        flush();
    }

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
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override {
        graph_->adjacent_outgoing_nodes(node, callback);
    }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override {
        graph_->call_outgoing_kmers(kmer, callback);
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override {
        graph_->call_incoming_kmers(kmer, callback);
    }

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override {
        graph_->adjacent_incoming_nodes(node, callback);
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

    virtual uint64_t num_nodes() const override { return graph_->num_nodes(); }
    virtual uint64_t max_index() const override { return graph_->max_index(); }

    virtual bool load(const std::string &filename) override {
        if (!graph_ptr_) {
            common::logger->error("load only supported for non-const graphs");
            exit(1);
        }

        if (!graph_ptr_->load(filename))
            return false;

        flush();
        return true;
    }

    virtual void serialize(const std::string &filename) const override {
        graph_->serialize(filename);
    }

    virtual std::string file_extension() const override { return graph_->file_extension(); }

    virtual const std::string& alphabet() const override { return graph_->alphabet(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override {
        return graph_->get_node_sequence(index);
    }

    virtual size_t get_k() const override { return graph_->get_k(); }

    virtual Mode get_mode() const override { return graph_->get_mode(); }

    virtual node_index traverse(node_index node, char next_char) const override {
        return graph_->traverse(node, next_char);
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override {
        return graph_->traverse_back(node, prev_char);
    }

    virtual size_t outdegree(node_index node) const override { return graph_->outdegree(node); }
    virtual size_t indegree(node_index node) const override { return graph_->indegree(node); }

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override {
        graph_->call_nodes(callback, stop_early);
    }

  protected:
    std::shared_ptr<const Graph> graph_;
    std::shared_ptr<Graph> graph_ptr_;

    // clear any internal storage the wrapper may have
    virtual void flush() = 0;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_WRAPPER__
