#ifndef __MASKED_GRAPH_HPP__
#define __MASKED_GRAPH_HPP__

#include <functional>
#include <vector>

#include "common/vectors/bitmap.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {

class MaskedDeBruijnGraph : public DeBruijnGraph {
  public:
    MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                        std::unique_ptr<bitmap>&& kmers_in_graph,
                        bool only_valid_nodes_in_mask = false,
                        Mode mode = BASIC);

    MaskedDeBruijnGraph(std::shared_ptr<const DeBruijnGraph> graph,
                        std::function<bool(node_index)>&& callback,
                        bool only_valid_nodes_in_mask = false,
                        Mode mode = BASIC);

    MaskedDeBruijnGraph(MaskedDeBruijnGraph&&) = default;
    MaskedDeBruijnGraph& operator=(MaskedDeBruijnGraph&&) = default;

    virtual ~MaskedDeBruijnGraph() {}

    virtual const DeBruijnGraph& get_base_graph() const override final { return graph_->get_base_graph(); }

    virtual void add_sequence(std::string_view,
                              const std::function<void(node_index)> &) override {
        throw std::runtime_error("Not implemented");
    }

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override;

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override;

    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override;

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override;

    virtual uint64_t num_nodes() const override { return kmers_in_graph_->num_set_bits(); }
    virtual uint64_t max_index() const override { return graph_->max_index(); }

    virtual bool load(const std::string &) override {
        throw std::runtime_error("Not implemented");
    }

    virtual void serialize(const std::string &) const override {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string file_extension() const override { return graph_->file_extension(); }

    virtual const std::string& alphabet() const override { return graph_->alphabet(); }

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override;

    virtual size_t get_k() const override { return graph_->get_k(); }

    virtual Mode get_mode() const override { return mode_; }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const override;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const override;

    virtual size_t outdegree(node_index) const override;
    virtual size_t indegree(node_index) const override;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const override;

    virtual const DeBruijnGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const DeBruijnGraph> get_graph_ptr() const { return graph_; }

    virtual inline bool in_subgraph(node_index node) const {
        assert(node > 0 && node <= graph_->max_index());
        assert(kmers_in_graph_.get());

        return (*kmers_in_graph_)[node];
    }

    virtual bool operator==(const MaskedDeBruijnGraph &other) const;
    virtual bool operator==(const DeBruijnGraph &other) const override;

    virtual void set_mask(bitmap *mask) { kmers_in_graph_.reset(mask); }

    virtual const bitmap& get_mask() const { return *kmers_in_graph_; }

  private:
    std::shared_ptr<const DeBruijnGraph> graph_;
    std::unique_ptr<bitmap> kmers_in_graph_;
    bool only_valid_nodes_in_mask_;
    Mode mode_;
};

} // namespace graph
} // namespace mtg

#endif // __MASKED_GRAPH_HPP__
