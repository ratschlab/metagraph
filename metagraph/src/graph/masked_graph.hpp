#ifndef __MASKED_GRAPH_HPP__
#define __MASKED_GRAPH_HPP__


#include <functional>
#include <vector>

#include "sequence_graph.hpp"
#include "bit_vector.hpp"


class MaskedDeBruijnGraph : public DeBruijnGraph {
  public:
    MaskedDeBruijnGraph(std::shared_ptr<DeBruijnGraph> graph, bit_vector *mask = nullptr);
    MaskedDeBruijnGraph(std::shared_ptr<DeBruijnGraph> graph, const std::function<bool(const node_index&)> &is_target_callback);

    void add_sequence(const std::string &sequence,
                      bit_vector_dyn *nodes_inserted = NULL);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    void map_to_nodes_sequentially(std::string::const_iterator begin,
                                   std::string::const_iterator end,
                                   const std::function<void(node_index)> &callback,
                                   const std::function<bool()> &terminate = [](){ return false; }) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    void adjacent_outgoing_nodes(node_index node,
                                 std::vector<node_index> *target_nodes) const;

    void call_outgoing_kmers(node_index kmer,
                             const OutgoingEdgeCallback &callback) const;

    void call_incoming_kmers(node_index kmer,
                             const IncomingEdgeCallback &callback) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    void adjacent_incoming_nodes(node_index node,
                                 std::vector<node_index> *source_nodes) const;

    uint64_t num_nodes() const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    std::string get_node_sequence(node_index index) const;

    size_t get_k() const { return graph_->get_k(); }

    // Traverse the outgoing edge
    node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    node_index traverse_back(node_index node, char prev_char) const;

    size_t outdegree(node_index node) const;
    size_t indegree(node_index) const;

    void call_outgoing_paths(node_index node,
                             const std::function<void(const std::string&)> &callback,
                             size_t max_num_steps) const;

    virtual void
    call_nodes(const std::function<void(const node_index&)> &callback,
               const std::function<bool()> &stop_early = [](){ return false; }) const override final;

    const DeBruijnGraph& get_graph() const { return *graph_; }
    std::shared_ptr<DeBruijnGraph> get_graph_ptr() const { return graph_; }

    uint64_t unmasked_outdegree(node_index node) const;

    inline bool in_graph(node_index node) const {
        return node == DeBruijnGraph::npos ? false : (is_target_callback_.get()
            ? (*is_target_callback_)(node)
            : (is_target_mask_.get()
                ? (*is_target_mask_)[node]
                : true));
    }

    bool operator==(const MaskedDeBruijnGraph &other) const;
    bool operator==(const DeBruijnGraph &other) const;

  private:
    std::shared_ptr<DeBruijnGraph> graph_;
    std::unique_ptr<bit_vector> is_target_mask_;
    std::unique_ptr<std::function<bool(const node_index&)>> is_target_callback_;
};


#endif // __MASKED_GRAPH_HPP__
