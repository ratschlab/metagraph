#ifndef __BOSS_HPP__
#define __BOSS_HPP__

#include <type_traits>

#include "common/vectors/bit_vector.hpp"
#include "common/vectors/wavelet_tree.hpp"
#include "kmer/kmer_extractor.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace boss {

class BOSSConstructor;

/**
 * This class implements the BOSS table, a succinct representation
 * of the de Bruijn graph following ideas and suggestions presented here:
 * http://link.springer.com/chapter/10.1007/978-3-642-33122-0_18
 *
 * There is also conceptual code available at
 * https://code.google.com/p/csalib/downloads/list
 * that has been used as a reference for this implementation.
 */
class BOSS {
  public:
    /**
     * Constant representing an invalid node or edge index.
     */
    static constexpr uint64_t npos = 0;
    // in [1,...,num_nodes], 0 = npos (invalid index)
    typedef uint64_t node_index;
    // in [1,...,num_edges], 0 = npos (invalid index)
    typedef uint64_t edge_index;

    typedef uint8_t TAlphabet;

    static constexpr auto ALWAYS_FALSE = []() { return false; };

    /**
     * Construct a BOSS table with the given node (k-mer) size.
     */
    explicit BOSS(size_t k = 1);
    explicit BOSS(BOSSConstructor *builder);
    ~BOSS();

    explicit BOSS(const BOSS &other) = delete;
    bool operator=(const BOSS &other) = delete;

    /**
     * Check whether the BOSS tables store identical k-mer sets.
     * FYI: this function reconstructs all the k-mers in O(k x n) time.
     */
    bool operator==(const BOSS &other) const;

    /**
     * Perform an element-wise comparison of the arrays W, last and
     * F and will only check for identity. If any element differs, return
     * false and true otherwise.
     */
    bool equals_internally(const BOSS &other, bool verbose = false) const;

    // Return the k-mer length
    size_t get_k() const { return k_; }
    uint64_t num_nodes() const;
    uint64_t num_edges() const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    bool load(std::ifstream &instream);
    void serialize(std::ofstream &outstream) const;

    /**
     * Load the index of node ranges constructed with index_suffix_ranges()
     * to speed up the search in the BOSS table.
     */
    bool load_suffix_ranges(std::ifstream &instream);
    void serialize_suffix_ranges(std::ofstream &outstream) const;

    // Traverse graph mapping k-mers from sequence to the graph edges
    // and run callback for each edge until the termination condition is satisfied
    // Invokes #callback with npos if a k-mer can't be mapped to the graph edges
    void map_to_edges(std::string_view sequence,
                      const std::function<void(edge_index)> &callback,
                      const std::function<bool()> &terminate = ALWAYS_FALSE,
                      const std::function<bool()> &skip = ALWAYS_FALSE) const;

    std::vector<edge_index> map_to_edges(std::string_view sequence) const;

    // |seq_encoded| must have no sentinels (zeros)
    void map_to_edges(const std::vector<TAlphabet> &seq_encoded,
                      const std::function<void(edge_index)> &callback,
                      const std::function<bool()> &terminate = ALWAYS_FALSE,
                      const std::function<bool()> &skip = ALWAYS_FALSE) const;

    // |seq_encoded| must have no sentinels (zeros)
    std::vector<edge_index>
    map_to_edges(const std::vector<TAlphabet> &seq_encoded) const;

    template <class... T>
    using Call = typename std::function<void(T...)>;

    /**
     * Given the last outgoing edge of a given node, call all edges from that node
     */
    void call_outgoing(edge_index edge, const Call<edge_index> &callback) const;

    /**
     * Traverse the boss graph and call all its edges
     * except for the dummy source nodes and the dummy sink nodes
     */
    void call_kmers(Call<edge_index, const std::string&> callback) const;

    // call all non-dummy edges without other adjacent incoming non-dummy edges
    void call_start_edges(Call<edge_index> callback) const;

    /**
     * Call contigs (or unitigs if |unitigs| is true) covering
     * all the edges of the BOSS graph (or its subgraph, if
     * |subgraph_mask| is specified), including the dummy ones.
     * Fetch only the primary sequences if |kmers_in_single_form|
     * is true. All the sentinel characters are omitted in this case.
     * If |kmers_in_single_form| is false, set |trim_dummy| to true
     * to fetch paths without sentinels.
     */
    void call_paths(Call<std::vector<edge_index> &&, std::vector<TAlphabet> &&> callback,
                    size_t num_threads = 1,
                    bool unitigs = false,
                    bool kmers_in_single_form = false,
                    const bitmap *subgraph_mask = NULL,
                    bool trim_sentinels = false) const;

    /**
     * Call contigs (dummy edges are skipped).
     */
    void call_sequences(Call<std::string&&, std::vector<edge_index>&&> callback,
                        size_t num_threads = 1,
                        bool kmers_in_single_form = false,
                        const bitmap *subgraph_mask = NULL) const;

    /**
     * Traversal starts at dummy source edges, then at forks and in the end at cycles.
     * A contig is terminated when we reach dead end or a fork where the rd successor
     * is visited, but not marked as being near a terminal edge.
     *
     * @param num_threads parallelize the graph traversal on this many threads
     * @param max_length maximum distance between two terminal nodes; this is a soft
     *        limit - in the worst case the distance between to terminal nodes can
     *        be 2*max_length
     */
    void row_diff_traverse(size_t num_threads,
                           size_t max_length,
                           const bit_vector &rd_succ,
                           sdsl::bit_vector *terminal) const;

    edge_index row_diff_successor(edge_index edge, const bit_vector &rd_succ) const {
        TAlphabet d = get_W(edge) % alph_size;
        assert(d != kSentinelCode && "sinks have no row-diff successors");
        // make one traversal step
        edge = fwd(edge, d);
        // pick the row-diff successor
        if (!get_last(edge - 1)) {
            while (!rd_succ[edge]) {
                edge--;
                assert(!get_last(edge) && "a row-diff successor must exist");
            }
        }
        return edge;
    }

    /**
     * Call unitigs (dummy edges are skipped).
     */
    void call_unitigs(Call<std::string&&, std::vector<edge_index>&&> callback,
                      size_t num_threads = 1,
                      size_t max_pruned_dead_end_size = 0,
                      bool kmers_in_single_form = false,
                      const bitmap *subgraph_mask = NULL) const;

    // |edge| must be the first incoming edge
    void call_incoming_to_target(edge_index edge, TAlphabet w,
                                 const Call<edge_index> &callback) const;

    /**
     * Add a full sequence to the graph.
     * If |try_extend| is true, search for the first k-mer in the graph
     * and extend it from that point. If the search fails, start from the dummy source.
     */
    void add_sequence(std::string_view seq,
                      bool try_extend = false,
                      std::vector<edge_index> *edges_inserted = NULL);

    // Given an edge list, remove them from the graph.
    // TODO: fix the implementation (anchoring the isolated nodes)
    void erase_edges_dyn(const std::set<edge_index> &edges);

    /**
     * Traverse the entire dummy subgraph (which is a tree)
     * and erase all redundant dummy edges.
     * If passed, mark |source_dummy_edges| with positions
     * of non-redundant dummy source edges.
     * Return value: edges removed from the initial graph.
     */
    sdsl::bit_vector
    erase_redundant_dummy_edges(sdsl::bit_vector *source_dummy_edges = NULL,
                                size_t num_threads = 0,
                                bool verbose = false);

    /**
     * Mark source dummy edges into mask by traversing the entire dummy subgraph (which
     * is a tree).
     * @param[out] mask a bit mask where sink dummy edges are marked. Must have the same
     * size as #W_ or must be nullptr.
     * @param[in] num_threads number of threads to use in the traversal (1 thread if <=1).
     * @return the number of source dummy edges
     */
    uint64_t mark_source_dummy_edges(sdsl::bit_vector *mask = NULL,
                                     size_t num_threads = 0) const;

    /**
     * Mark sink dummy edges into mask. Does not include the main dummy edge (with
     * edge_index = 1). A sink dummy edge is an outgoing edge labeled with $.
     * @param[out] mask a bit mask where sink dummy edges are marked. Must have the same
     * size as #W_ or must be nullptr.
     * @return the number of sink dummy edges (i.e. the number of bits that would be set
     * in #mask, if provided)
     */
    uint64_t mark_sink_dummy_edges(sdsl::bit_vector *mask = NULL) const;

    // Mark npos, i.e. 0, as well as all the source
    // and all the sink dummy edges in graph.
    sdsl::bit_vector mark_all_dummy_edges(size_t num_threads) const;

    // Prune redundant dummy edges in graph
    // and mark all dummy edges that cannot be removed.
    sdsl::bit_vector prune_and_mark_all_dummy_edges(size_t num_threads);

    /**
     * Depth first edge traversal.
     * Traverse all edges reachable from the given one.
     */
    void edge_DFT(edge_index start,
                  Call<edge_index> pre_visit,
                  Call<edge_index> post_visit,
                  std::function<bool(edge_index)> end_branch) const;

    /**
    * Heavily borrowing from the graph sequence traversal, this function gets
    * a graph `other` and merges its nodes into the target graph object.
    * The edges of `other` are fully traversed and nodes are added if not existing yet.
    * This function is well suited to merge small graphs into large ones.
    */
    void merge(const BOSS &other, size_t num_threads = 1);

    /**
     * Using the offset structure F this function returns the value of the last
     * position of the source node for edge i.
     */
    TAlphabet get_node_last_value(edge_index i) const;

    /**
     * Given an edge index i, this function returns the k-mer sequence of its
     * source node as a string.
     */
    std::string get_node_str(edge_index i) const;

    /**
     * Given an edge index i, this function returns the k-mer sequence of its
     * source node.
     */
    std::vector<TAlphabet> get_node_seq(edge_index i) const;

    /**
     * Given an edge index i, this function returns the k-mer sequence of its
     * source node, and the node whose last character corresponds to the first
     * character of the sequence.
     */
    std::pair<std::vector<TAlphabet>, edge_index> get_node_seq_with_end_node(edge_index i) const;

    /**
     * Given index i of an edge and a value k, this function
     * returns the k-th last character of the source node for edge i.
     */
    std::pair<TAlphabet, edge_index> get_minus_k_value(edge_index i, size_t k) const;

    /**
     * Index ranges of nodes for all possible suffixes (alph_size-1)^t
     * Suffixes with sentinel characters are not indexed.
     * After the index is constructed, it speeds up search in the BOSS table
     * by narrowing down the initial node range and skipping several fwd calls.
     */
    void index_suffix_ranges(size_t suffix_length);

    size_t get_indexed_suffix_length() const { return indexed_suffix_length_; }

    /**
     * Print current representation of the graph to stream.
     */
    void print(std::ostream &os = std::cout) const;

    /**
     * Print vectors F, L, and W to stream.
     */
    void print_internal_representation(std::ostream &os = std::cout) const;

    /**
     * Alternative representation states of the BOSS table.
     *
     * STAT: provides the best space/time trade-off
     *      Representation:
     *          last -- bit_vector_stat
     *             W -- wavelet_tree_stat
     *
     * SMALL: is the smallest, useful for storage or when RAM is limited
     *      Representation:
     *          last -- bit_vector_small
     *             W -- wavelet_tree_small
     *
     * FAST: is the fastest but large
     *      Representation:
     *          last -- bit_vector_stat
     *             W -- wavelet_tree_fast
     *
     * DYN: is a dynamic representation supporting insert and delete
     *      Representation:
     *          last -- bit_vector_dyn
     *             W -- wavelet_tree_dyn
     */
    enum State { SMALL = 1, DYN, STAT, FAST };

    State get_state() const { return state; }
    void switch_state(State state);

    /**
     * Write the adjacency list to file |filename| or
     * print it to stdout of the filename is not provided.
     */
    void print_adj_list(std::ostream &os = std::cout) const;


    /**
     * Given an edge index i, this function returns true if that is
     * the only outgoing edge from its source node.
     */
    bool is_single_outgoing(edge_index i) const;

    /**
     * Given an edge index i (first incoming) and label d, this function returns
     * the number of edges incoming to its target node.
     */
    size_t num_incoming_to_target(edge_index i, TAlphabet d) const;

    /**
     * Given an edge index i and label w, this function returns
     * true if that is the only edge incoming to its target node.
     */
    bool is_single_incoming(edge_index i, TAlphabet w) const;

    /**
     * Given an edge index |i| and a character |c|, get the index of an adjacent
     * incoming edge with the first character c if such exists and npos otherwise.
     */
    edge_index pick_incoming_edge(edge_index i, TAlphabet c) const;

    /**
     * Given an edge index |i| and a character |c|, get the index of the edge with
     * label c outgoing from the same source node if such exists and npos otherwise.
     */
    edge_index pick_edge(edge_index i, TAlphabet c) const;

    /**
     * Given a node label kmer, this function returns the index
     * of the corresponding node or the closest predecessor, if no node
     * with the sequence is not found.
     */
    node_index pred_kmer(const std::vector<TAlphabet> &kmer) const;

    /**
     * Return value of last at position i.
     */
    bool get_last(edge_index i) const { return (*last_)[i]; }
    const bit_vector& get_last() const { return *last_; }

    /**
     * Transforms a boss edge index to the index of its source node.
     * Uses the object's array last and a position and
     * returns the number of set bits up to that postion.
     */
    node_index rank_last(edge_index i) const;

    /**
     * Transforms a boss node index to the index of its last outgoing edge.
     * Uses the object's array last and a given position i and
     * returns the position of the i-th set bit in last[1..i].
     */
    edge_index select_last(node_index i) const;

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the last set bit in last[1..i].
     */
    edge_index pred_last(edge_index i) const;

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the first set bit in last[i..N].
     */
    edge_index succ_last(edge_index i) const;

    /**
     * Return value of W at position i.
     */
    TAlphabet get_W(edge_index i) const { return (*W_)[i]; }
    const wavelet_tree& get_W() const { return *W_; }

    /**
     * For the given position i in W and a character c from the alphabet,
     * return the number of occurrences of c in W up to (including) position i.
     */
    uint64_t rank_W(edge_index i, TAlphabet c) const;

    /**
     * Return the position of the last occurrence of |c| in W[1..i] or zero
     * if such does not exist.
     */
    edge_index pred_W(edge_index i, TAlphabet c) const;

    /**
     * For characters |first| and |second|, return the last occurrence
     * of them in W[1..i], i.e. max(pred_W(i, first), pred_W(i, second)).
     */
    edge_index pred_W(edge_index i, TAlphabet first, TAlphabet second) const;

    /**
     * Return position of the first occurrence of |c| in W[i..N].
     */
    edge_index succ_W(edge_index i, TAlphabet c) const;

    /**
     * For characters |first| and |second|, return the first occurrence
     * of them in W[i..N], i.e. min(succ_W(i, first), succ_W(i, second)).
     */
    std::pair<edge_index, TAlphabet>
    succ_W(edge_index i, TAlphabet first, TAlphabet second) const;

    /**
     * Return the offset for the edges with the penultimate character c.
     */
    edge_index get_F(TAlphabet c) const { return F_.at(c); }

    /**
     * This functions gets a position i reflecting the r-th occurrence of the corresponding
     * character c in W and returns the position of the r-th occurrence of c in last.
     */
    edge_index fwd(edge_index i, TAlphabet c) const;

    /**
     * This function gets a position i that reflects the i-th node and returns the
     * position in W that corresponds to the i-th node's last character.
     */
    edge_index bwd(edge_index i) const;

    // Given the alphabet index return the corresponding symbol
    char decode(TAlphabet s) const;
    std::string decode(const std::vector<TAlphabet> &seq_encoded) const;
    // Given the alphabet character return its corresponding number
    TAlphabet encode(char s) const;
    std::vector<TAlphabet> encode(std::string_view sequence) const;

    /**
     * Given iterators to an input sequence, this function finds the index range
     * of nodes with the maximal length suffix matching a prefix of the sequence.
     * The tuple includes the indices of the boundary nodes (inclusive) and an
     * iterator to the first character in the input not matched (end if the full
     * sequence matched). If a match is not found, it returns std::make_tuple(0, 0, begin)
     */
    template <typename RandomAccessIt>
    inline std::tuple<edge_index, edge_index, RandomAccessIt>
    index_range(RandomAccessIt begin, RandomAccessIt end) const;

    inline bool tighten_range(edge_index *rl, edge_index *ru, TAlphabet s) const;

    /**
     * The size of the alphabet for kmers that this graph encodes.
     * For DNA, this value is 5 ($,A,C,G,T).
     */
    const TAlphabet alph_size;
    /**
     * The symbols of the alphabet, for DNA this will be "$ACGT".
     */
    const std::string &alphabet;

    static constexpr size_t kSentinelCode = 0;
    static constexpr size_t kSentinel = '$';

  private:
    // file dump extension
    static constexpr auto kExtension = ".dbg";

    const mtg::kmer::KmerExtractorBOSS kmer_extractor_;
    const size_t bits_per_char_W_;

    // k-mer size
    size_t k_;
    // the bit array indicating the last outgoing edge of a node
    bit_vector *last_;
    // offsets for the boss edges with a given penultimate character
    // F_[c] points to the last boss edge with its penultimate character < c
    std::vector<edge_index> F_;
    // node offsets: NF_[c] = rank_last(F_[c]), a helper array
    std::vector<node_index> NF_;
    // the array containing the edge labels
    wavelet_tree *W_;

    State state = State::DYN;

    size_t indexed_suffix_length_ = 0;
    std::vector<std::pair<edge_index, edge_index>> indexed_suffix_ranges_;

    /**
     * This function gets a character c and updates the edge offsets F_
     * by incrementing them with +1 (for edge insertion) or decrementing
     * (for edge delition) and updates the node offsets NF_ accordingly.
     */
    void update_F(TAlphabet c, int delta, bool is_representative);

    /**
     * Recompute the node offsets NF_ from F_. Call after changes in F_.
     */
    void recompute_NF();

    /**
     * Given a character c and an edge index, this function
     * creates an outgoing edge from the same source node with
     * label c if it is not a part of the graph yet.
     */
    edge_index append_pos(TAlphabet c, edge_index source_node,
                          const TAlphabet *source_node_kmer,
                          std::vector<edge_index> *edges_inserted = NULL);

    /**
     * Helper function used by the append_pos function
     */
    edge_index insert_edge(TAlphabet c, edge_index begin, edge_index end);

    /**
     * Erase exactly all the masked edges from the graph,
     * may invalidate the graph (if leaves nodes with no incoming edges).
     * Returns the number of edges erased.
     */
    uint64_t erase_edges(const sdsl::bit_vector &edges_to_remove_mask);

    /**
     * This function gets two edge indices and returns if their source
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(edge_index first, edge_index second) const;
    /**
     * This function gets an edge indix and checks if its source
     * node has the same k-1 suffix as k-mer |second|.
     */
    bool compare_node_suffix(edge_index first, const TAlphabet *second) const;

    /**
     * Given an edge index i, this function returns the k-mer sequence of its
     * source node, and the node whose last character corresponds to the first
     * character of the sequence. If the graph is suffix indexed and the third
     * returned value is true, then the returned node is the last node visited after
     * k - indexed_suffix_length_ bwd steps.
     */
    std::tuple<std::vector<TAlphabet>, edge_index, bool>
    get_node_seq_with_end_node_indexed(edge_index i) const;

    /**
     * Given a (k+1)-mer, this function returns the index
     * of the corresponding edge, if such exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    inline edge_index map_to_edge(RandomAccessIt begin, RandomAccessIt end) const;

    /**
     * Maps a prefix of the query to a range of nodes in the BOSS table
     * with the same suffixes.
     * If suffix ranges are indexed in the BOSS table and the query is longer
     * than |indexed_suffix_length_| and contains no sentinels, take a prefix
     * of that size.
     * Otherwise, match only the first character *begin.
     * Returns: range and length of the prefix matched.
     */
    template <typename RandomAccessIt>
    inline std::tuple<edge_index, edge_index, size_t>
    get_initial_range(RandomAccessIt begin, RandomAccessIt end) const;

    /**
     * Given a k-mer, this function returns the index of last edge going out
     * from the k-mer's corresponding node, if such a node exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    inline edge_index index(RandomAccessIt begin, RandomAccessIt end) const;

    void verbose_cout() const {}

#ifdef DBGDEBUG
    template <typename T, typename... Targs>
    void verbose_cout(const T &arg, Targs ...rest) const {
        std::cout << arg << std::flush;
        verbose_cout(rest...);
    }
#else
    template <typename T, typename... Targs>
    void verbose_cout(const T&, Targs...) const {}
#endif

    bool is_valid() const;

  public:
    class Chunk;
    void initialize(Chunk *chunk);
};

std::ostream& operator<<(std::ostream &os, const BOSS &graph);


template <typename RandomAccessIt>
std::tuple<BOSS::edge_index, BOSS::edge_index, size_t>
BOSS::get_initial_range(RandomAccessIt begin, RandomAccessIt end) const {
    assert(std::all_of(begin, end, [&](TAlphabet c) { return c < alph_size; }));

    // node range
    edge_index rl, ru;
    size_t offset;

    if (indexed_suffix_length_
            && begin + indexed_suffix_length_ <= end
            // only node suffixes without sentinels are indexed
            && !std::count(begin, begin + indexed_suffix_length_, kSentinelCode)) {
        // find range of nodes with suffixes matching
        // prefix [begin, begin + indexed_suffix_length_)
        uint64_t index = 0;
        // use the co-lex order to assign close indexes to suffixes of
        // nodes close in the boss table
        for (auto it = begin + indexed_suffix_length_ - 1, end = begin - 1; it != end; --it) {
            assert(*it && *it < alph_size);
            // shift the alphabet: suffixes with sentinels '$' are not indexed
            index = index * (alph_size - 1) + (*it - 1);
        }

        // query range
        std::tie(rl, ru) = indexed_suffix_ranges_[index];
        offset = indexed_suffix_length_;

    } else {
        // get first
        TAlphabet s = *begin;
        // use the initial range
        rl = F_.at(s) + 1 < W_->size()
             ? F_.at(s) + 1
             : W_->size(); // lower bound
        ru = s + 1 < alph_size
             ? F_[s + 1]
             : W_->size() - 1; // upper bound

        offset = 1;
    }

    return std::make_tuple(rl, ru, offset);
}

inline bool BOSS::tighten_range(edge_index *rl, edge_index *ru, TAlphabet s) const {
    // Tighten the range of edges including only those ending with |s| and do fwd.
    uint64_t rk_rl = rank_W(*rl - 1, s) + 1;
    uint64_t rk_ru = rank_W(*ru, s);
    if (rk_rl > rk_ru)
        return false;

    // select the index of the position in last that is rank many positions after offset
    *rl = select_last(NF_[s] + rk_rl - 1) + 1;
    *ru = select_last(NF_[s] + rk_ru);
    return true;
}

template <typename RandomAccessIt>
BOSS::edge_index BOSS::index(RandomAccessIt begin, RandomAccessIt end) const {
    static_assert(std::is_same_v<std::decay_t<decltype(*begin)>, TAlphabet>,
                  "Only encoded sequences can be queried");
    assert(begin + k_ == end);
    assert(std::all_of(begin, end, [&](TAlphabet c) { return c <= alph_size; }));

    // return npos if invalid characters are found
    if (std::find(begin, end, alph_size) != end)
        return npos;

    // get the initial node range
    auto [rl, ru, offset] = get_initial_range(begin, end);
    if (rl > ru)
        return npos;

    // update range iteratively while scanning through s
    for (auto it = begin + offset; it != end; ++it) {
        if (!tighten_range(&rl, &ru, *it))
            return npos;
    }
    assert(succ_last(rl) <= ru);
    return ru;
}

template <typename RandomAccessIt>
std::tuple<BOSS::edge_index, BOSS::edge_index, RandomAccessIt>
BOSS::index_range(RandomAccessIt begin, RandomAccessIt end) const {
    static_assert(std::is_same_v<std::decay_t<decltype(*begin)>, TAlphabet>,
                  "Only encoded sequences can be queried");

    assert(end >= begin);
    assert(end <= begin + k_);
    assert(std::all_of(begin, end, [&](TAlphabet c) { return c <= alph_size; }));

    if (begin == end)
        return std::make_tuple((edge_index)1, (edge_index)1, begin);

    // check if all characters belong to the alphabet
    if (std::find(begin, end, alph_size) != end)
        return std::make_tuple((edge_index)0, (edge_index)0, begin);

    // get the initial node range
    auto [rl, ru, offset] = get_initial_range(begin, end);
    if (rl > ru) {
        // start search from scratch
        TAlphabet s = *begin;

        rl = F_.at(s) + 1 < W_->size()
                ? F_.at(s) + 1
                : W_->size(); // lower bound
        ru = s + 1 < alph_size
                ? F_[s + 1]
                : W_->size() - 1; // upper bound

        if (rl > ru)
            return std::make_tuple((edge_index)0, (edge_index)0, begin);

        offset = 1;
    }

    auto it = begin + offset;
    // update range iteratively while scanning through s
    for (; it != end; ++it) {
        if (!tighten_range(&rl, &ru, *it))
            return std::make_tuple(succ_last(rl), ru, it);
    }
    assert(succ_last(rl) <= ru);
    return std::make_tuple(succ_last(rl), ru, it);
}

template <typename RandomAccessIt>
BOSS::edge_index
BOSS::map_to_edge(RandomAccessIt begin, RandomAccessIt end) const {
    assert(begin + k_ + 1 == end);
    assert(std::all_of(begin, end, [&](TAlphabet c) { return c <= alph_size; }));

    edge_index edge = index(begin, end - 1);

    return edge && *(end - 1) < alph_size
            ? pick_edge(edge, *(end - 1))
            : npos;
}

inline void BOSS::call_outgoing(edge_index edge, const Call<edge_index> &callback) const {
    assert(get_last(edge));
    do {
        callback(edge);
    } while (--edge && !get_last(edge));
}

} // namespace boss
} // namespace graph
} // namespace mtg

#endif // __BOSS_HPP__
