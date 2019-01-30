#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include "sequence_graph.hpp"


template <typename Weights = sdsl::int_vector<>>
class WeightedDBG : public DeBruijnGraph {
  public:
    using weight = typename Weights::value_type;
    using node_index = typename DeBruijnGraph::node_index;

    WeightedDBG(std::shared_ptr<DeBruijnGraph> graph, Weights&& weights);

    virtual ~WeightedDBG() {}

    /**
     * Weight functions
     */

    // virtual void set_weights() { weights_ = std::move(weights); }
    virtual weight get_weight(node_index i) const { return weights_[i]; }

    /**
     * Graph functions
     */

    virtual size_t get_k() const { return graph_->get_k(); }

    virtual bool is_canonical_mode() const { return graph_->is_canonical_mode(); }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const {
        return graph_->traverse(node, next_char);
    }
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const {
        return graph_->traverse_back(node, prev_char);
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are not mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate) const {
        graph_->map_to_nodes_sequentially(begin, end, callback, terminate);
    }

    // Perform extension on a provided seed based on the string iterator.
    // If seed is npos, perform seeding automatically.
    // Extend until the termination condition is satisfied or reached the end of the query.
    // In canonical mode, non-canonical k-mers are not mapped to canonical ones
    virtual void extend_from_seed(std::string::const_iterator,
                                           std::string::const_iterator,
                                           const std::function<void(node_index)>&,
                                           const std::function<bool()>&,
                                           node_index) const {
        // TODO: Complete extend_from_seed.
        throw std::runtime_error("Not implemented");
    }

    virtual size_t outdegree(node_index node) const { return graph_->outdegree(node); }
    virtual size_t indegree(node_index node) const { return graph_->indegree(node); }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const {
        graph_->call_outgoing_kmers(kmer, callback);
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const {
        graph_->call_incoming_kmers(kmer, callback);
    }

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    virtual void add_sequence(const std::string & /* sequence */,
                              bit_vector_dyn * /* nodes_inserted */) {
        throw std::runtime_error("Not implemented");
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const {
        graph_->map_to_nodes(sequence, callback, terminate);
    }

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         std::vector<node_index> *target_nodes) const {
        graph_->adjacent_outgoing_nodes(node, target_nodes);
    }
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    virtual void adjacent_incoming_nodes(node_index node,
                                         std::vector<node_index> *source_nodes) const {
        graph_->adjacent_incoming_nodes(node, source_nodes);
    }

    virtual uint64_t num_nodes() const { return graph_->num_nodes(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index node) const {
        return graph_->get_node_sequence(node);
    }

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;
    virtual std::string file_extension() const { return graph_->file_extension(); }

    virtual const std::string& alphabet() const { return graph_->alphabet(); }

  private:
    std::shared_ptr<DeBruijnGraph> graph_;
    Weights weights_;

    static constexpr auto kWeightsExtension = ".weights";
};

#endif // __WEIGHTED_GRAPH_HPP__
