#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include "sequence_graph.hpp"
#include "bit_vector.hpp"
#include "config.hpp"
#include "boss.hpp"


class DBGSuccinct : public DeBruijnGraph {
  public:
    explicit DBGSuccinct(size_t k, bool canonical_mode = false);
    explicit DBGSuccinct(BOSS *boss_graph, bool canonical_mode = false);

    virtual ~DBGSuccinct() {}

    virtual size_t get_k() const override final;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(const std::string &sequence,
                      double discovery_fraction = 1) const override final;

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const override final;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const override final;

    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the outgoing edges and pushes back indices of their target nodes.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         std::vector<node_index> *target_nodes) const override final;
    // Given a node index and a pointer to a vector of node indices, iterates
    // over all the incoming edges and pushes back indices of their source nodes.
    virtual void adjacent_incoming_nodes(node_index node,
                                         std::vector<node_index> *source_nodes) const override final;

    // Insert sequence to graph and mask the inserted nodes if |nodes_inserted|
    // is passed. If passed, |nodes_inserted| must have length equal
    // to the number of nodes in graph.
    virtual void add_sequence(const std::string &sequence,
                              bit_vector_dyn *nodes_inserted = NULL) override final;

    virtual std::string get_node_sequence(node_index) const override final;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string::const_iterator begin,
                                           std::string::const_iterator end,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override final;

    // TODO: compare these to DeBruijnGraph implementations
    virtual void call_sequences(const std::function<void(const std::string&)> &callback) const override final;

    virtual void call_unitigs(const std::function<void(const std::string&)> &callback,
                              size_t min_tip_size = 1) const override final;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override final;

    virtual void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const override final;

    virtual void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const override final {
        // TODO: Implement
        throw std::runtime_error("Not implemented");
    }

    virtual size_t outdegree(node_index node) const override final;
    virtual size_t indegree(node_index node) const override final;
    virtual uint64_t num_nodes() const override final;

    virtual void mask_dummy_kmers(size_t num_threads, bool with_pruning) final;

    virtual bool load_without_mask(const std::string &filename_base) final;
    virtual bool load(const std::string &filename_base) override;
    virtual void serialize(const std::string &filename_base) const override;
    virtual std::string file_extension() const override final { return kExtension; }

    virtual void switch_state(Config::StateType new_state) final;
    virtual Config::StateType get_state() const final;

    virtual bool is_canonical_mode() const override final { return canonical_mode_; }

    virtual const BOSS& get_boss() const final { return *boss_graph_; }
    virtual BOSS& get_boss() final { return *boss_graph_; }
    virtual BOSS* release_boss() final { return boss_graph_.release(); }

    virtual bool operator==(const DeBruijnGraph &other) const override final;

    virtual const std::string& alphabet() const override final;

    virtual void print(std::ostream &out) const override final;


  private:
    void add_seq(const std::string &sequence, bit_vector_dyn *nodes_inserted);
    uint64_t kmer_to_boss_index(node_index kmer_index) const;
    node_index boss_to_kmer_index(uint64_t boss_index) const;

    std::unique_ptr<BOSS> boss_graph_;
    // all edges in boss except dummy
    std::unique_ptr<bit_vector> valid_edges_;

    bool canonical_mode_;

    static constexpr auto kExtension = ".dbg";
    static constexpr auto kDummyMaskExtension = ".edgemask";
};


#endif // __DBG_SUCCINCT_HPP__
