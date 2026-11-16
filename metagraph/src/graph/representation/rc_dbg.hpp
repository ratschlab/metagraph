#ifndef __RC_DBG_HPP__
#define __RC_DBG_HPP__

#include <algorithm>

#include "graph/representation/base/dbg_wrapper.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"

namespace mtg {
namespace graph {

class DBGSuccinct;

/**
 * RCDBG is a wrapper which represents the reverse complement of the underlying Graph.
 * e.g., get_node_sequence(n) := reverse_complement(graph.get_node_sequence(n))
 * e.g., traverse(n, c) := graph.traverse_back(n, complement(c))
 */
class RCDBG : public DBGWrapper<DeBruijnGraph> {
  public:
    template <typename... Args>
    RCDBG(Args&&... args) : DBGWrapper<DeBruijnGraph>(std::forward<Args>(args)...) {
        if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph_.get()))
            add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));
    }

    virtual node_index traverse(node_index node, char next_char) const override final {
        return graph_->traverse_back(node, complement(next_char));
    }

    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        return graph_->traverse(node, complement(prev_char));
    }

    virtual std::unique_ptr<NodeStream>
    node_stream(std::string_view sequence, bool sequential) const override final {
        std::string rc(sequence);
        ::reverse_complement(rc.begin(), rc.end());
        std::vector<node_index> nodes = sequential
                ? graph::map_to_nodes_sequentially(*graph_, rc)
                : graph::map_to_nodes(*graph_, rc);
        std::reverse(nodes.begin(), nodes.end());
        return std::make_unique<BufferedNodeStream>(std::move(nodes));
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

    virtual Edges outgoing_kmers(node_index kmer) const override final {
        const auto cache = get_extension<NodeFirstCache>();
        Edges edges = cache ? cache->incoming_kmers(kmer) : graph_->incoming_kmers(kmer);
        for (auto &[prev, c] : edges) {
            c = complement(c);
        }
        return edges;
    }

    virtual Edges incoming_kmers(node_index kmer) const override final {
        Edges edges = graph_->outgoing_kmers(kmer);
        for (auto &[next, c] : edges) {
            c = complement(c);
        }
        return edges;
    }

    virtual Nodes adjacent_outgoing_nodes(node_index node) const override final {
        return graph_->adjacent_incoming_nodes(node);
    }
    virtual Nodes adjacent_incoming_nodes(node_index node) const override final {
        return graph_->adjacent_outgoing_nodes(node);
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

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; },
                            size_t num_threads = 1,
                            size_t batch_size = 1'000'000) const override final {
        // all node IDs are the same
        graph_->call_nodes(callback, stop_early, num_threads, batch_size);
    }

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override {
        // all node IDs are the same, but represent the reverse complement sequences
        graph_->call_kmers([&](node_index node, const std::string &seq) {
            std::string rseq(seq);
            ::reverse_complement(rseq.begin(), rseq.end());
            callback(node, rseq);
        }, stop_early);
    }
};

} // namespace mtg
} // namespace graph

#endif // __RC_DBG_HPP__
