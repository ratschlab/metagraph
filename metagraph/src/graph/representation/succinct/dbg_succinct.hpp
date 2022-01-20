#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include "common/vectors/bit_vector.hpp"
#include "common/caches.hpp"
#include "kmer/kmer_bloom_filter.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/base/dbg_wrapper.hpp"
#include "boss.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct : public DeBruijnGraph {
  public:
    friend class MaskedDeBruijnGraph;

    /**
     * A wrapper which caches computed node sequences for DBGSuccinct graphs.
     * This allows for faster get_node_sequence and call_incoming_kmers calls.
     *
     * In particular, for each cached node, it stores its decoded k-mer sequence
     * and the BOSS node whose last character is the first character of the k-mer.
     *
     * call_incoming_kmers, call_outgoing_kmers, and map_to_nodes_sequentially
     * are overridden to fill the cache alongside their normal functions.
     *
     * Bytes stored per k-mer: ~56 (k<=21), ~64 (k<=42), ~80 (k<=85)
     */
    class CachedView;

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
    // For each k-mer satisfying the skip condition, run the callback on npos.
    // Guarantees that nodes are called in the same order as the input sequence.
    // In canonical mode, non-canonical k-mers are NOT mapped to canonical ones
    virtual void map_to_nodes_sequentially(std::string_view sequence,
                                           const std::function<void(node_index)> &callback,
                                           const std::function<bool()> &terminate = [](){ return false; },
                                           const std::function<bool()> &skip = [](){ return false; }) const override final;

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
    // Initialize DBGSuccinct and dump to disk without loading to RAM.
    // FYI: Note that suffix ranges will not be indexed.
    static void serialize(boss::BOSS::Chunk&& chunk,
                          const std::string &filename_base,
                          Mode mode,
                          boss::BOSS::State state = boss::BOSS::State::STAT);
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

    uint64_t kmer_to_boss_index(node_index kmer_index) const;
    node_index boss_to_kmer_index(uint64_t boss_index) const;

    void initialize_bloom_filter_from_fpr(double false_positive_rate,
                                          uint32_t max_num_hash_functions = -1);

    void initialize_bloom_filter(double bits_per_kmer,
                                 uint32_t max_num_hash_functions = -1);

    const mtg::kmer::KmerBloomFilter<>* get_bloom_filter() const { return bloom_filter_.get(); }

    std::unique_ptr<CachedView> get_cached_view(size_t cache_size = 100'000) const;

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


class DBGSuccinct::CachedView : public DBGWrapper<DBGSuccinct> {
  public:
    typedef boss::BOSS::edge_index edge_index;
    typedef boss::BOSS::TAlphabet TAlphabet;

    virtual ~CachedView() {}

    // get the encoding of the first character of this node's sequence
    virtual TAlphabet get_first_value(edge_index i) const = 0;

    /**
     * Methods from DBGWrapper
     */
    virtual node_index get_base_node(node_index node) const override final { return node; }
    virtual bool is_base_node(node_index node) const override final { return graph_->is_base_node(node); }

    virtual std::pair<std::vector<node_index>, bool /* is reversed */>
    get_base_path(const std::vector<node_index> &path,
                  const std::string&) const override final {
        return std::make_pair(path, false);
    }

    /**
     * Methods from DeBruijnGraph
     */
    virtual bool operator==(const DeBruijnGraph &other) const override final;

    #define DELEGATE_METHOD(RETURN_TYPE, METHOD, ARG_TYPE, ARG_NAME) \
    virtual RETURN_TYPE METHOD(ARG_TYPE ARG_NAME) const override final { \
        return graph_->METHOD(ARG_NAME); \
    }

    DELEGATE_METHOD(void, call_kmers, const std::function<void(node_index, const std::string&)> &, callback)
    DELEGATE_METHOD(void, call_source_nodes, const std::function<void(node_index)> &, callback)

    DELEGATE_METHOD(uint64_t, num_nodes, , )
    DELEGATE_METHOD(uint64_t, max_index, , )
    DELEGATE_METHOD(size_t, outdegree, node_index, node)
    DELEGATE_METHOD(bool, has_single_outgoing, node_index, node)
    DELEGATE_METHOD(bool, has_multiple_outgoing, node_index, node)
    DELEGATE_METHOD(size_t, indegree, node_index, node)
    DELEGATE_METHOD(bool, has_no_incoming, node_index, node)
    DELEGATE_METHOD(bool, has_single_incoming, node_index, node)
    DELEGATE_METHOD(node_index, kmer_to_node, std::string_view, kmer)

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override final {
        graph_->call_nodes(callback, stop_early);
    }

    virtual void call_sequences(const CallPath &callback,
                                size_t num_threads = 1,
                                bool kmers_in_single_form = false) const override final {
        graph_->call_sequences(callback, num_threads, kmers_in_single_form);
    }

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final {
        graph_->call_unitigs(callback, num_threads, min_tip_size, kmers_in_single_form);
    }

    virtual bool find(std::string_view sequence, double discovery_fraction = 1) const override final {
        return graph_->find(sequence, discovery_fraction);
    }

    virtual void map_to_nodes(std::string_view sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate
                                  = [](){ return false; }) const override final {
        graph_->map_to_nodes(sequence, callback, terminate);
    }

    virtual void
    adjacent_outgoing_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_outgoing_nodes(node, callback);
    }

    virtual void
    adjacent_incoming_nodes(node_index node,
                            const std::function<void(node_index)> &callback) const override final {
        graph_->adjacent_incoming_nodes(node, callback);
    }

    virtual void serialize(const std::string &filename_base) const override final {
        graph_->serialize(filename_base);
    }

    // TODO: these can be overloaded to cache values, but this functionality is not
    //       needed now.
    virtual node_index traverse(node_index node, char next_char) const override final {
        return graph_->traverse(node, next_char);
    }
    virtual node_index traverse_back(node_index node, char prev_char) const override final {
        return graph_->traverse_back(node, prev_char);
    }
    virtual void traverse(node_index start,
                          const char *begin,
                          const char *end,
                          const std::function<void(node_index)> &callback,
                          const std::function<bool()> &terminate
                              = [](){ return false; }) const override final {
        graph_->traverse(start, begin, end, callback, terminate);
    }

  protected:
    template <typename Graph>
    explicit CachedView(Graph&& graph) : DBGWrapper(std::forward<Graph>(graph)) {}
};


} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_HPP__
