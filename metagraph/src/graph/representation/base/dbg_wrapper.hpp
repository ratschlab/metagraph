#ifndef __DBG_WRAPPER__
#define __DBG_WRAPPER__

#include "sequence_graph.hpp"


namespace mtg {
namespace graph {

/**
 * This abstract class stores a graph internally and transfers all calls to it.
 * This wrapper uses the default methods from Graph when available. This may be
 * used when the nodes and edges of the wrapped graph differ from the underlying
 * graph.
 */
template <class Graph = DeBruijnGraph>
class DBGWrapper : public DeBruijnGraph {
  public:
    template <class InGraph>
    DBGWrapper(std::shared_ptr<const InGraph> graph)
          : graph_(std::dynamic_pointer_cast<const Graph>(graph)) { assert(graph_); }

    template <class InGraph>
    DBGWrapper(std::shared_ptr<InGraph> graph)
          : graph_(std::dynamic_pointer_cast<const Graph>(graph)) { assert(graph_); }

    // aliasing constructors
    template <class InGraph>
    explicit DBGWrapper(const InGraph &graph)
          : graph_(std::shared_ptr<const Graph>{}, dynamic_cast<const Graph*>(&graph)) {
        assert(graph_);
    }

    template <class InGraph>
    explicit DBGWrapper(InGraph &graph)
          : graph_(std::shared_ptr<const Graph>{}, dynamic_cast<const Graph*>(&graph)) {
        assert(graph_);
    }

    virtual ~DBGWrapper() {}

    /**
     * Added methods
     */
    virtual const Graph& get_graph() const { return *graph_; }
    virtual std::shared_ptr<const Graph> get_graph_ptr() const { return graph_; }
    virtual void set_graph(std::shared_ptr<const Graph> graph) = 0;

    /**
     * Methods shared by all wrappers
     */
    virtual std::string file_extension() const override final { return graph_->file_extension(); }
    virtual size_t get_k() const override final { return graph_->get_k(); }
    virtual const std::string& alphabet() const override final { return graph_->alphabet(); }
    virtual void print(std::ostream &out) const override { graph_->print(out); }
    virtual Mode get_mode() const override { return graph_->get_mode(); }
    virtual const DeBruijnGraph& get_base_graph() const override final {
        return graph_->get_base_graph();
    }

    /**
     * Override these if the wrapper changes the graph's indexing
     */
    virtual uint64_t num_nodes() const override { return graph_->num_nodes(); }
    virtual uint64_t max_index() const override { return graph_->max_index(); }

    /**
     * The Graph defaults of these are likely to break in a wrapped graph,
     * so these should be implemented.
     */
    virtual bool operator==(const DeBruijnGraph &other) const override = 0;

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override = 0;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override = 0;

  protected:
    std::shared_ptr<const Graph> graph_;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_WRAPPER__
