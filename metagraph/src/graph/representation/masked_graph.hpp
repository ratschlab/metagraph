#ifndef __MASKED_GRAPH_HPP__
#define __MASKED_GRAPH_HPP__

#include <functional>
#include <vector>

#include "common/vectors/bitmap.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"


namespace mtg {
namespace graph {

class MaskedDeBruijnGraph : public DBGWrapper<DeBruijnGraph> {
  public:
    template <class Graph>
    MaskedDeBruijnGraph(Graph&& graph,
                        std::unique_ptr<bitmap>&& kmers_in_graph,
                        bool only_valid_nodes_in_mask = false,
                        Mode mode = BASIC)
          : DBGWrapper<DeBruijnGraph>(std::forward<Graph>(graph)),
            kmers_in_graph_(std::move(kmers_in_graph)),
            only_valid_nodes_in_mask_(only_valid_nodes_in_mask),
            mode_(mode) {
        assert(kmers_in_graph_.get());
        assert(kmers_in_graph_->size() == graph_->max_index() + 1);

        if (graph_->get_mode() == PRIMARY && mode_ != PRIMARY) {
            throw std::runtime_error("Any subgraph of a primary graph is primary, thus"
                                     " the mode of the subgraph must be set to PRIMARY");
        }

        if (graph_->get_mode() != CANONICAL && mode_ == CANONICAL)
            throw std::runtime_error("Canonical subgraph requires canonical base graph");
    }

    template <class Graph>
    MaskedDeBruijnGraph(Graph&& graph,
                        std::function<bool(node_index)>&& callback,
                        bool only_valid_nodes_in_mask = false,
                        Mode mode = BASIC)
          : MaskedDeBruijnGraph(std::forward<Graph>(graph),
                                std::make_unique<bitmap_lazy>(
                                    std::move(callback), graph->max_index() + 1),
                                only_valid_nodes_in_mask, mode) {}

    virtual ~MaskedDeBruijnGraph() {}

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

    // Get string corresponding to |node_index|.
    // Note: Not efficient if sequences in nodes overlap. Use sparingly.
    virtual std::string get_node_sequence(node_index index) const override;

    virtual Mode get_mode() const override { return mode_; }

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const override;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const override;

    virtual size_t outdegree(node_index) const override;
    virtual size_t indegree(node_index) const override;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early = [](){ return false; }) const override;

    virtual inline bool in_subgraph(node_index node) const {
        assert(node > 0 && node <= max_index());
        assert(kmers_in_graph_.get());

        return (*kmers_in_graph_)[node];
    }

    virtual bool operator==(const MaskedDeBruijnGraph &other) const;
    virtual bool operator==(const DeBruijnGraph &other) const override;

    virtual void set_mask(bitmap *mask) { kmers_in_graph_.reset(mask); }

    virtual const bitmap& get_mask() const { return *kmers_in_graph_; }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override;

    virtual node_index kmer_to_node(std::string_view kmer) const override {
        node_index node = graph_->kmer_to_node(kmer);
        return (*kmers_in_graph_)[node] ? node : npos;
    }

    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override {
        graph_->call_source_nodes([&](node_index node) {
            if ((*kmers_in_graph_)[node])
                callback(node);
        });
    }

    virtual void serialize(const std::string &) const override final {
        throw std::runtime_error("serialize not supported for MaskedDeBruijnGraph");
    }

    virtual bool load(const std::string &) override final {
        throw std::runtime_error("koad not supported for MaskedDeBruijnGraph");
    }

    virtual void add_sequence(std::string_view,
                              const std::function<void(node_index)>&) override final {
        throw std::runtime_error("add_sequence not supported for MaskedDeBruijnGraph");
    }

  private:
    std::unique_ptr<bitmap> kmers_in_graph_;
    bool only_valid_nodes_in_mask_;
    Mode mode_;
};

} // namespace graph
} // namespace mtg

#endif // __MASKED_GRAPH_HPP__
