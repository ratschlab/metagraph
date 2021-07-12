#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include "common/vectors/bit_vector.hpp"
#include "kmer/kmer_bloom_filter.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "boss.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct : public DeBruijnGraph {
  public:
    friend class MaskedDeBruijnGraph;

    explicit DBGSuccinct(size_t k, Mode mode = BASIC);
    explicit DBGSuccinct(boss::BOSS *boss_graph, Mode mode = BASIC);

    virtual ~DBGSuccinct() {}

    virtual size_t get_k() const override final;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(std::string_view sequence,
                      double discovery_fraction = 1) const override final;

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char next_char) const override final;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char prev_char) const override final;

    // Given a node index, call the target nodes of all edges outgoing from it.
    virtual void adjacent_outgoing_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;
    // Given a node index, call the source nodes of all edges incoming to it.
    virtual void adjacent_incoming_nodes(node_index node,
                                         const std::function<void(node_index)> &callback) const override final;

    // Insert sequence to graph and invoke callback |on_insertion| for each new
    // node index augmenting the range [1,...,max_index], including those not
    // pointing to any real node in graph. That is, the callback is invoked for
    // all new real nodes and all new dummy node indexes allocated in graph.
    // In short: max_index[after] = max_index[before] + {num_invocations}.
    virtual void add_sequence(std::string_view sequence,
                              const std::function<void(node_index)> &on_insertion = [](node_index) {}) override final;

    virtual std::string get_node_sequence(node_index node) const override final;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const override;

    // Traverse graph mapping sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; }) const override final;

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override final;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback) const override final;

    // Find nodes with a common suffix matching the maximal prefix of the string |str|,
    // and call these nodes. If more than |max_num_allowed_matches| are found,
    // or if the maximal prefix is shorter than |min_match_length|, return
    // without calling.
    void call_nodes_with_suffix_matching_longest_prefix(
            std::string_view str,
            std::function<void(node_index, uint64_t /* match length */)> callback,
            size_t min_match_length = 1,
            size_t max_num_allowed_matches = std::numeric_limits<size_t>::max()) const;

    // Given a starting node, traverse the graph forward following the edge
    // sequence delimited by begin and end. Terminate the traversal if terminate()
    // returns true, or if the sequence is exhausted.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate = [](){ return false; }) const override final;

    virtual void call_outgoing_kmers(node_index, const OutgoingEdgeCallback&) const override final;

    virtual void call_incoming_kmers(node_index, const IncomingEdgeCallback&) const override final;

    virtual size_t outdegree(node_index) const override final;
    virtual bool has_single_outgoing(node_index) const override final;
    virtual bool has_multiple_outgoing(node_index) const override final;
    virtual size_t indegree(node_index) const override final;
    virtual bool has_no_incoming(node_index) const override final;
    virtual bool has_single_incoming(node_index) const override final;

    /**
     * Returns the number of nodes (k-mers) in the graph, which is equal to the number of
     * edges in the BOSS graph (because an edge in the BOSS graph represents a k-mer).
     */
    virtual uint64_t num_nodes() const override final;

    virtual void mask_dummy_kmers(size_t num_threads, bool with_pruning) final;

    // Return a pointer to the mask, or NULL if not initialized
    virtual const bit_vector* get_mask() const final { return valid_edges_.get(); }

    virtual void reset_mask() final { valid_edges_.reset(); }
    virtual bit_vector* release_mask() final { return valid_edges_.release(); }

    virtual bool load_without_mask(const std::string &filename_base) final;
    virtual bool load(const std::string &filename_base) override;
    virtual void serialize(const std::string &filename_base) const override;
    virtual std::string file_extension() const override final { return kExtension; }
    std::string bloom_filter_file_extension() const { return kBloomFilterExtension; }

    /*
     * Available representations:
     *  STAT: provides the best space/time trade-off
     *      Representation:
     *            BOSS::last -- bit_vector_stat
     *               BOSS::W -- wavelet_tree_stat
     *           valid_edges -- bit_vector_small
     *
     *  SMALL: is the smallest, useful for storage or when RAM is limited
     *      Representation:
     *            BOSS::last -- bit_vector_small
     *               BOSS::W -- wavelet_tree_small
     *           valid_edges -- bit_vector_small
     *
     *  FAST: is the fastest but large
     *      Representation:
     *            BOSS::last -- bit_vector_stat
     *               BOSS::W -- wavelet_tree_fast
     *           valid_edges -- bit_vector_stat
     *
     *  DYN: is a dynamic representation supporting insert and delete
     *      Representation:
     *            BOSS::last -- bit_vector_dyn
     *               BOSS::W -- wavelet_tree_dyn
     *           valid_edges -- bit_vector_dyn
     */
    virtual void switch_state(boss::BOSS::State new_state) final;
    virtual boss::BOSS::State get_state() const final;

    virtual Mode get_mode() const override final { return mode_; }

    virtual const boss::BOSS& get_boss() const final { return *boss_graph_; }
    virtual boss::BOSS& get_boss() final { return *boss_graph_; }
    virtual boss::BOSS* release_boss() final { return boss_graph_.release(); }

    virtual bool operator==(const DeBruijnGraph &other) const override final;

    virtual const std::string& alphabet() const override final;

    virtual void print(std::ostream &out) const override final;

    virtual void call_source_nodes(const std::function<void(node_index)> &callback) const override final;

    node_index row_diff_successor(node_index node, const bit_vector &rd_succ) const;
    template <class Callback>
    void call_row_diff_successors(node_index node, const bit_vector &rd_succ, const Callback &callback) const;

    void add_rd_successors_at_forks(size_t num_threads,
                                    const sdsl::bit_vector &anchors,
                                    sdsl::bit_vector *rd_succ,
                                    size_t max_length) const;

    uint64_t kmer_to_boss_index(node_index kmer_index) const;
    node_index boss_to_kmer_index(uint64_t boss_index) const;

    void initialize_bloom_filter_from_fpr(double false_positive_rate,
                                          uint32_t max_num_hash_functions = -1);

    void initialize_bloom_filter(double bits_per_kmer,
                                 uint32_t max_num_hash_functions = -1);

    const mtg::kmer::KmerBloomFilter<>* get_bloom_filter() const { return bloom_filter_.get(); }

    static constexpr auto kExtension = ".dbg";
    static constexpr auto kDummyMaskExtension = ".edgemask";
    static constexpr auto kBloomFilterExtension = ".bloom";

  private:
    std::unique_ptr<boss::BOSS> boss_graph_;
    // all edges in boss except dummy
    std::unique_ptr<bit_vector> valid_edges_;

    Mode mode_;

    std::unique_ptr<mtg::kmer::KmerBloomFilter<>> bloom_filter_;
};

inline DBGSuccinct::node_index
DBGSuccinct::row_diff_successor(node_index node, const bit_vector &rd_succ) const {
    const boss::BOSS &boss = *boss_graph_;
    boss::BOSS::edge_index edge = kmer_to_boss_index(node);
    boss::BOSS::TAlphabet d = boss.get_W(edge) % boss.alph_size;
    assert(d != boss::BOSS::kSentinelCode && "sinks have no row-diff successors");
    // make one traversal step
    edge = boss.fwd(edge, d);
    node = boss_to_kmer_index(edge);

    if (!rd_succ.size() || boss.get_last(edge - 1))
        return node;

    // pick the row-diff successor
    while (!rd_succ[node]) {
        node--;
        edge--;
        assert(!boss.get_last(edge) && "a row-diff successor must exist");
    }
    return node;
}

template <class Callback>
inline void DBGSuccinct::call_row_diff_successors(node_index node,
                                                  const bit_vector &rd_succ,
                                                  const Callback &callback) const {
    const boss::BOSS &boss = *boss_graph_;
    boss::BOSS::edge_index edge = kmer_to_boss_index(node);
    boss::BOSS::TAlphabet d = boss.get_W(edge) % boss.alph_size;
    assert(d != boss::BOSS::kSentinelCode && "sinks have no row-diff successors");
    // make one traversal step
    edge = boss.fwd(edge, d);
    node = boss_to_kmer_index(edge);

    if (!rd_succ.size() || boss.get_last(edge - 1)) {
        callback(node);
        return;
    }

    // pick the row-diff successor
    do {
        if (rd_succ[node])
            callback(node);
        node--;
    } while (!boss.get_last(--edge));
}

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_HPP__
