#ifndef __RC_DBG_HPP__
#define __RC_DBG_HPP__

#include "graph/representation/base/dbg_wrapper.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {

/**
 * RCDBG is a wrapper which represents the reverse complement of the underlying Graph.
 * e.g., get_node_sequence(n) := reverse_complement(graph.get_node_sequence(n))
 * e.g., traverse(n, c) := graph.traverse_back(n, complement(c))
 */
class RCDBG : public DBGNodeModifyingWrapper<DeBruijnGraph> {
  public:
    template <typename... Args>
    RCDBG(Args&&... args)
          : DBGNodeModifyingWrapper<DeBruijnGraph>(std::forward<Args>(args)...) {}

    virtual uint64_t num_nodes() const override final { return graph_->num_nodes(); }
    virtual uint64_t max_index() const override final { return graph_->max_index(); }

    virtual node_index traverse(node_index node, char next_char) const override final {
        return graph_->traverse_back(node, complement(next_char));
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        return graph_->traverse(node, complement(prev_char));
    }

    virtual void
    map_to_nodes_sequentially(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override final {
        if (terminate())
            return;

        std::string rc(sequence);
        ::reverse_complement(rc.begin(), rc.end());
        std::vector<node_index> nodes = map_sequence_to_nodes(*graph_, rc);
        for (auto it = nodes.rbegin(); it != nodes.rend() && !terminate(); ++it) {
            callback(*it);
        }
    }

    virtual size_t outdegree(node_index node) const override final {
        return graph_->indegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override final {
        return graph_->has_single_incoming(node);
    }

    virtual size_t indegree(node_index node) const override final {
        return graph_->outdegree(node);
    }

    virtual bool has_single_incoming(node_index node) const override final {
        return graph_->has_single_outgoing(node);
    }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override final {
        graph_->call_incoming_kmers(kmer, [&](node_index prev, char c) {
            callback(prev, complement(c));
        });
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override final {
        graph_->call_outgoing_kmers(kmer, [&](node_index next, char c) {
            callback(next, complement(c));
        });
    }

    virtual void map_to_nodes(std::string_view /* sequence */,
                              const std::function<void(node_index)> & /* callback */,
                              const std::function<bool()> & /* terminate */
                                  = [](){ return false; }) const override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_incoming_nodes(node, callback);
    }
    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_outgoing_nodes(node, callback);
    }

    virtual std::string get_node_sequence(node_index node) const override final {
        std::string rc = graph_->get_node_sequence(node);
        ::reverse_complement(rc.begin(), rc.end());
        return rc;
    }

    virtual bool operator==(const DeBruijnGraph &other) const override final {
        if (const auto *other_rc = dynamic_cast<const RCDBG*>(&other)) {
            return *graph_ == *other_rc->graph_;
        } else {
            return DeBruijnGraph::operator==(other);
        }
    }

    virtual bool load(const std::string &filename) override final {
        if (!graph_ptr_)
            throw std::runtime_error("load only supported for non-const graphs");

        return graph_ptr_->load(filename);
    }

    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion
                                  = [](node_index) {}) override final {
        if (!graph_ptr_)
            throw std::runtime_error("load only supported for non-const graphs");

        graph_ptr_->add_sequence(sequence, on_insertion);
    }

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override {
        // all node IDs are the same
        graph_->call_nodes(callback, stop_early);
    }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override {
        // all node IDs are the same, but represent the reverse complement sequences
        graph_->call_kmers([&](node_index node, const std::string &seq) {
            std::string rseq(seq);
            ::reverse_complement(rseq.begin(), rseq.end());
            callback(node, rseq);
        });
    }
};

} // namespace mtg
} // namespace graph

#endif // __RC_DBG_HPP__
