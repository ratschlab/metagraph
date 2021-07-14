#ifndef __RC_DBG_HPP__
#define __RC_DBG_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {

/**
 * RCDBG is a wrapper which represents the reverse complement of the underlying
 * DeBruijnGraph.
 * e.g., get_node_sequence(n) := reverse_complement(graph.get_node_sequence(n))
 * e.g., traverse(n, c) := graph.traverse_back(n, complement(c))
 */
class RCDBG : public DeBruijnGraph {
  public:
    RCDBG(const DeBruijnGraph &graph) : graph_(graph) {}

    virtual const DeBruijnGraph& get_base_graph() const override final { return graph_.get_base_graph(); }

    virtual size_t get_k() const override final { return graph_.get_k(); }
    virtual Mode get_mode() const override final { return graph_.get_mode(); }

    virtual node_index traverse(node_index node, char next_char) const override final {
        return graph_.traverse_back(node, complement(next_char));
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        return graph_.traverse(node, complement(prev_char));
    }

    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override final {
        if (terminate())
            return;

        std::string rc(sequence);
        ::reverse_complement(rc.begin(), rc.end());
        std::vector<node_index> nodes = map_sequence_to_nodes(graph_, rc);
        for (auto it = nodes.rbegin(); it != nodes.rend() && !terminate(); ++it) {
            callback(*it);
        }
    }

    virtual size_t outdegree(node_index node) const override final {
        return graph_.indegree(node);
    }

    virtual bool has_single_outgoing(node_index node) const override final {
        return graph_.has_single_incoming(node);
    }

    virtual bool has_multiple_outgoing(node_index node) const override final {
        return graph_.indegree(node) > 1;
    }

    virtual size_t indegree(node_index node) const override final {
        return graph_.outdegree(node);
    }

    virtual bool has_no_incoming(node_index node) const override final {
        return graph_.outdegree(node) == 0;
    }

    virtual bool has_single_incoming(node_index node) const override final {
        return graph_.has_single_outgoing(node);
    }

    virtual void call_outgoing_kmers(node_index kmer,
                                     const OutgoingEdgeCallback &callback) const override final {
        graph_.call_incoming_kmers(kmer, [&](node_index prev, char c) {
            callback(prev, complement(c));
        });
    }

    virtual void call_incoming_kmers(node_index kmer,
                                     const IncomingEdgeCallback &callback) const override final {
        graph_.call_outgoing_kmers(kmer, [&](node_index next, char c) {
            callback(next, complement(c));
        });
    }

    virtual const std::string& alphabet() const override final { return graph_.alphabet(); }


    virtual void add_sequence(std::string_view /* sequence */,
                              const std::function<void(node_index)> & /* on_insertion */ = [](uint64_t) {}) override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void map_to_nodes(std::string_view /* sequence */,
                              const std::function<void(node_index)> & /* callback */,
                              const std::function<bool()> & /* terminate */ = [](){ return false; }) const override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final {
        graph_.adjacent_incoming_nodes(node, callback);
    }
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final {
        graph_.adjacent_outgoing_nodes(node, callback);
    }

    virtual uint64_t num_nodes() const override final { return graph_.num_nodes(); }

    virtual bool load(const std::string &) override final {
        throw std::runtime_error("Not implemented");
    }

    virtual void serialize(const std::string &) const override final {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string file_extension() const override final {
        throw std::runtime_error("Not implemented");
    }

    virtual std::string get_node_sequence(node_index node) const override final {
        std::string rc = graph_.get_node_sequence(node);
        ::reverse_complement(rc.begin(), rc.end());
        return rc;
    }

  private:
    const DeBruijnGraph &graph_;
};

} // namespace mtg
} // namespace graph

#endif // __RC_DBG_HPP__
