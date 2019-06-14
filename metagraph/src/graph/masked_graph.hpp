#ifndef __MASKED_GRAPH_HPP__
#define __MASKED_GRAPH_HPP__


#include <functional>
#include <vector>

#include "sequence_graph.hpp"
#include "bitmap.hpp"


class MaskedDeBruijnGraph : public DeBruijnGraph {
  public:
    MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph, bitmap *mask = nullptr);

    MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                        std::function<bool(const DeBruijnGraph::node_index&)>&& callback,
                        size_t num_set_bits = -1);

    virtual void add_sequence(const std::string &sequence,
                              bit_vector_dyn *nodes_inserted = NULL);

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         std::vector<node_index> *target_nodes) const;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const;

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    virtual void adjacent_incoming_nodes(node_index node,
                                         std::vector<node_index> *source_nodes) const;

    virtual uint64_t num_nodes() const;

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;
    virtual std::string file_extension() const { return graph_->file_extension(); }

    virtual const std::string& alphabet() const { return graph_->alphabet(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const;

    virtual size_t get_k() const { return graph_->get_k(); }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const;

    virtual size_t outdegree(node_index node) const;
    virtual size_t indegree(node_index) const;

    virtual void
    call_nodes(const std::function<void(const node_index&)> &callback,
               const std::function<bool()> &stop_early = [](){ return false; }) const final;

    virtual const DeBruijnGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const DeBruijnGraph> get_graph_ptr() const { return graph_; }

    virtual uint64_t unmasked_outdegree(node_index node) const;
    virtual uint64_t unmasked_indegree(node_index node) const;

    virtual inline bool in_graph(node_index node) const {
        return node == DeBruijnGraph::npos || !is_target_mask_.get()
            ? false
            : (*is_target_mask_)[node];
    }

    virtual bool operator==(const MaskedDeBruijnGraph &other) const;
    virtual bool operator==(const DeBruijnGraph &other) const;

    virtual void set_mask(bitmap *mask) { is_target_mask_.reset(mask); }

    virtual const bitmap* get_mask() const { return is_target_mask_.get(); }

  private:
    std::shared_ptr<const DeBruijnGraph> graph_;
    std::unique_ptr<bitmap> is_target_mask_;
};


#endif // __MASKED_GRAPH_HPP__
