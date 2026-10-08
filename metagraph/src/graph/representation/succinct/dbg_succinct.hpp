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

    virtual void call_nodes(const std::function<void(node_index)> &callback,
                            const std::function<bool()> &terminate = [](){ return false; },
                            size_t num_threads = 1,
                            size_t batch_size = 1'000'000) const override final;

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
                                bool kmers_in_single_form = false,
                                bool verbose = common::get_verbose()) const override final;

    virtual void call_unitigs(const CallPath &callback,
                              size_t num_threads = 1,
                              size_t min_tip_size = 1,
                              bool kmers_in_single_form = false) const override final;

    virtual void call_kmers(const std::function<void(node_index, const std::string&)> &callback,
                            const std::function<bool()> &stop_early
                                = [](){ return false; }) const override final;

    // Find nodes with a common suffix matching the maximal prefix of the string |str|,
    // and call these nodes. If more than |max_num_allowed_matches| are found,
    // or if the maximal prefix is shorter than |min_match_length|, return
    // without calling.
    void call_nodes_with_suffix_matching_longest_prefix(
            std::string_view str,
            std::function<void(node_index, uint64_t /* match length */)> callback,
            size_t min_match_length = 1,
            size_t max_num_allowed_matches = std::numeric_limits<size_t>::max()) const;

    /**
     * Pattern-search primitives (docs/DESIGN-pattern-search.md §4.1), beside
     * call_nodes_with_suffix_matching_longest_prefix, which they leave as it is: they
     * count first and call only what the caller asks for (its TODO), and they never count
     * a dummy or pruned edge as a k-mer (LastSymbolEdges::candidates is an upper bound that
     * includes them; the caller resolves it with the scans when invalid_non_sentinel > 0).
     * Each works on a normalised BOSS edge range [first, last] (1 <= first <= last <=
     * max_index(); whole node groups, as BOSS::tighten_range returns them, or [1,
     * max_index()] for the empty suffix), and each requires the valid-edge mask
     * (get_mask() != NULL): without it a dummy edge is indistinguishable from a k-mer.
     * The mask is trusted as written: the primitives assume what mask_dummy_kmers
     * guarantees (build/transform --mask-dummy, --pattern-build-mask), that every dummy
     * edge, every edge with W = $ included, is 0 in it. A mask that marks a dummy valid (as
     * DBGSuccinct::add_sequence on a masked graph writes for the dummies it inserts, see its
     * TODO) makes them count it as a k-mer. Of that premise the pattern search checks the
     * W = $ half once per load (count_valid_sentinel_edges; refused as mask_invalid); a
     * source dummy marked valid is not detected (`metagraph extend` re-masks its output).
     */

    // The valid edges of the range: the k-mers leaving its nodes (two ranks).
    uint64_t count_valid_edges_in_range(node_index first, node_index last) const;

    // The edges of a range that can end with one symbol c (see count_edges_with_last_symbol)
    struct LastSymbolEdges {
        // edges with W in {c, c + alph_size} (plain + marked): the k-mers u.c of the
        // range's nodes u, valid or not; an upper bound of the valid ones
        uint64_t candidates = 0;
        // invalid edges of the range
        uint64_t invalid = 0;
        // invalid edges whose W is not the sentinel ($ or its marked form): the only invalid
        // edges that can be candidates, since every edge with W = $ is a dummy sink, so
        // candidates - invalid_non_sentinel is a lower bound of the valid candidates
        uint64_t invalid_non_sentinel = 0;
    };

    /**
     * Counts the edges of the range whose W is c (plain) or c + alph_size (marked: another
     * edge into the same target carries c), with O(1) ranks. The valid ones among them
     * (the k-mers of the range's nodes ending with c) number exactly |candidates| when
     * |invalid_non_sentinel| is 0; otherwise the caller resolves them with
     * next_invalid_edge or next_edge_with_last_symbol, a scan it can charge and interrupt.
     * 1 <= c < alph_size.
     */
    LastSymbolEdges count_edges_with_last_symbol(node_index first, node_index last,
                                                 boss::BOSS::TAlphabet c) const;

    // The first edge e in [from, last] with W[e] in {c, c + alph_size}, valid or not;
    // npos if there is none (from > last included).
    node_index next_edge_with_last_symbol(node_index from, node_index last,
                                          boss::BOSS::TAlphabet c) const;

    // The first valid edge in [from, last]; npos if there is none.
    node_index next_valid_edge(node_index from, node_index last) const;

    // The first invalid (dummy or pruned) edge in [from, last]; npos if there is none.
    node_index next_invalid_edge(node_index from, node_index last) const;

    /**
     * The edges with W = $ (plain or marked: the sink dummies and the main dummy edge 1) that
     * the mask marks valid, which the primitives above assume none is (see their comment).
     * 0 for every mask mask_dummy_kmers builds; positive for a mask that add_sequence updated
     * (`metagraph extend` on a masked graph before it re-masked, or a stale or foreign
     * .edgemask of the right size). O(number of W = $ edges) selects on W, never a pass over
     * all edges: meant to run once per load (the pattern search refuses such a mask,
     * mask_invalid). 0 when there is no mask.
     */
    uint64_t count_valid_sentinel_edges() const;

    /**
     * The same primitives for a graph without the valid-edge mask (owner decision #16 of
     * 2026-10-08: the pattern search counts upper bounds there). They read W only, never the
     * mask, and count every edge that can carry a pattern base: a k-mer, or a source dummy
     * (a k-mer starting with '$', which only the range's unspelled node symbols can hold;
     * BOSS::node_has_sentinel tells one). A sink dummy (W = $) never carries a base and is
     * never counted. Same range conventions as above.
     */

    // The edges of the range whose W is not $ (plain or marked): the k-mers leaving its
    // nodes and the source dummies among them (four ranks).
    uint64_t count_non_sink_edges_in_range(node_index first, node_index last) const;

    // The edges of the range with W in {c, c + alph_size}: LastSymbolEdges::candidates,
    // without the mask's invalid counts (two ranks per form). 1 <= c < alph_size.
    uint64_t count_edges_with_symbol(node_index first, node_index last,
                                     boss::BOSS::TAlphabet c) const;

    // The first edge in [from, last] whose W is not $; npos if there is none.
    node_index next_non_sink_edge(node_index from, node_index last) const;

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
    virtual uint64_t max_index() const override final;

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

    virtual bool in_graph(node_index node) const override final;
    node_index validate_edge(node_index node) const;
    node_index select_node(uint64_t rank) const;
    uint64_t rank_node(node_index node) const;

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

    std::unique_ptr<bit_vector> generate_valid_kmer_mask(size_t num_threads, bool with_pruning) const;
};

} // namespace graph
} // namespace mtg

#endif // __DBG_SUCCINCT_HPP__
