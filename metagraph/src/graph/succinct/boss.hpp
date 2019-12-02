#ifndef __BOSS_HPP__
#define __BOSS_HPP__

#include <type_traits>

#include "common/config.hpp"
#include "graph/base/sequence_graph.hpp"
#include "kmer/kmer_extractor.hpp"
#include "utils/bit_vectors/bit_vector.hpp"
#include "utils/bit_vectors/wavelet_tree.hpp"

class BOSSConstructor;

auto ALWAYS_FALSE = []() { return false; };

/**
 * This class implements the BOSS table, a succinct representation
 * of the de Bruijn graph following ideas and suggestions presented here:
 * https://www.researchgate.net/publication/262273032_Succinct_de_Bruijn_Graphs
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

    /**
     * Construct a boos table with the given kmer size.
     */
    explicit BOSS(size_t k = 1);
    explicit BOSS(BOSSConstructor *builder);
    ~BOSS();

    explicit BOSS(const BOSS &other) = delete;
    bool operator=(const BOSS &other) = delete;

    /**
     * Check whether BOSS tables store the same data.
     * FYI: this function reconstructs all the kmers, and
     * the complexity is at least O(k x n).
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

    // Traverse graph mapping k-mers from sequence to the graph edges
    // and run callback for each edge until the termination condition is satisfied
    // Invokes #callback with npos if a k-mer can't be mapped to the graph edges
    void map_to_edges(const std::string &sequence,
                      const std::function<void(edge_index)> &callback,
                      const std::function<bool()> &terminate = ALWAYS_FALSE,
                      const std::function<bool()> &skip = ALWAYS_FALSE) const;

    std::vector<edge_index> map_to_edges(const std::string &sequence) const;

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
     * Traverse boss graph and call all its edges
     * except for the dummy source of sink ones
     */
    void call_kmers(Call<edge_index, const std::string&> callback) const;

    // call all non-dummy edges without other adjacent incoming non-dummy edges
    void call_start_edges(Call<edge_index> callback) const;

    // call contigs (or unitigs if |unitigs| is true) that cover
    // exactly all edges in graph (or subgraph, if |subgraph_mask| is passed)
    void call_paths(Call<std::vector<edge_index>&&,
                         std::vector<TAlphabet>&&> callback,
                    bool unitigs = false,
                    bool kmers_in_single_form = false,
                    const bitmap *subgraph_mask = NULL) const;

    void call_sequences(Call<std::string&&, std::vector<edge_index>&&> callback,
                        bool kmers_in_single_form = false,
                        const bitmap *subgraph_mask = NULL) const;

    void call_unitigs(Call<std::string&&, std::vector<edge_index>&&> callback,
                      size_t max_pruned_dead_end_size = 0,
                      bool kmers_in_single_form = false,
                      const bitmap *subgraph_mask = NULL) const;

    // |edge| must be the first incoming edge
    void call_incoming_to_target(edge_index edge,
                                 std::function<void(edge_index)> callback) const;

    /**
     * Add a full sequence to the graph.
     * If |try_extend| is true, search for the first k-mer in the graph
     * and extend it from that point. If the search fails, start from the dummy source.
     */
    void add_sequence(const std::string &seq,
                      bool try_extend = false,
                      std::vector<uint64_t> *edges_inserted = NULL);

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
     * Does not include the main dummy edge (with edge_index = 1).
     * @param[out] mask a bit mask where sink dummy edges are marked. Must have the same
     * size as #W_ or must be nullptr.
     * @param[in] num_threads number of threads to use in the traversal (1 thread if <=1).
     * @param[in] verbose logging verbosity
     * @return the number of sink dummy edges
     */
    uint64_t mark_source_dummy_edges(sdsl::bit_vector *mask = NULL,
                                     size_t num_threads = 0,
                                     bool verbose = false) const;

    /**
     * Mark sink dummy edges into mask. Does not include the main dummy edge (with
     * edge_index = 1). A sink dummy edge is an outgoing edge labeled with $.
     * @param[out] mask a bit mask where sink dummy edges are marked. Must have the same
     * size as #W_ or must be nullptr.
     * @return the number of sink dummy edges
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
    void merge(const BOSS &other);

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
     * Given index i of an edge and a value k, this function
     * returns the k-th last character of the source node for edge i.
     */
    std::pair<TAlphabet, edge_index> get_minus_k_value(edge_index i, size_t k) const;

    /**
     * Print current representation of the graph to stream.
     */
    void print(std::ostream &os = std::cout) const;

    /**
     * Print vectors F, L, and W to stream.
     */
    void print_internal_representation(std::ostream &os = std::cout) const;

    Config::StateType get_state() const { return state; }
    void switch_state(Config::StateType state);

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
     * Given an edge index i (first incoming), this function returns
     * the number of edges incoming to its target node.
     */
    size_t num_incoming_to_target(edge_index i) const;

    /**
     * Given an edge index i, this function returns true if that is
     * the only edge incoming to its target node.
     */
    bool is_single_incoming(edge_index i) const;

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
    bool get_last(uint64_t i) const { return (*last_)[i]; }
    const bit_vector& get_last() const { return *last_; }

    /**
     * Uses the object's array last and a position and
     * returns the number of set bits up to that postion.
     */
    uint64_t rank_last(uint64_t i) const;

    /**
     * Uses the object's array last and a given position i and
     * returns the position of the i-th set bit in last[1..i].
     */
    uint64_t select_last(uint64_t i) const;

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the last set bit in last[1..i].
     */
    uint64_t pred_last(uint64_t i) const;

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the first set bit in last[i..N].
     */
    uint64_t succ_last(uint64_t i) const;

    /**
     * Return value of W at position i.
     */
    TAlphabet get_W(uint64_t i) const { return (*W_)[i]; }
    const wavelet_tree& get_W() const { return *W_; }

    /**
     * Uses the object's array W, a given position i in W and a character c
     * from the alphabet and returns the number of occurences of c in W up to
     * position i.
     */
    uint64_t rank_W(uint64_t i, TAlphabet c) const;

    /**
     * Uses the array W and gets a count i and a character c from
     * the alphabet and returns the positions of the i-th occurence of c in W.
     */
    uint64_t select_W(uint64_t i, TAlphabet c) const;

    /**
     * For characters |first| and |second|, return the last occurrence
     * of them in W[1..i], i.e. max(pred_W(i, first), pred_W(i, second)).
     */
    uint64_t pred_W(uint64_t i, TAlphabet first, TAlphabet second) const;

    /**
     * Return position of the first occurrence of |c| in W[i..N].
     */
    uint64_t succ_W(uint64_t i, TAlphabet c) const;

    /**
     * For characters |first| and |second|, return the first occurrence
     * of them in W[i..N], i.e. min(succ_W(i, first), succ_W(i, second)).
     */
    std::pair<uint64_t, TAlphabet> succ_W(uint64_t i, TAlphabet first, TAlphabet second) const;

    /**
     * Return value of F vector at index k.
     * The index is over the alphabet!
     */
    uint64_t get_F(TAlphabet k) const { return F_.at(k); }
    const std::vector<uint64_t>& get_F() const { return F_; }

    /**
     * This functions gets a position i reflecting the r-th occurence of the corresponding
     * character c in W and returns the position of the r-th occurence of c in last.
     */
    uint64_t fwd(uint64_t i) const;

    /**
     * This function gets a position i that reflects the i-th node and returns the
     * position in W that corresponds to the i-th node's last character.
     */
    uint64_t bwd(uint64_t i) const;

    // Given the alphabet index return the corresponding symbol
    char decode(TAlphabet s) const;
    std::string decode(const std::vector<TAlphabet> &seq_encoded) const;
    // Given the alphabet character return its corresponding number
    TAlphabet encode(char s) const;
    std::vector<TAlphabet> encode(const std::string &sequence) const;

    /**
     * Given iterators to an input sequence, this function finds the index range
     * of nodes with the maximal length suffix matching a prefix of the sequence.
     * The tuple includes the indices of the boundary nodes (inclusive) and an
     * iterator to the first character in the input not matched (end if the full
     * sequence matched). If a match is not found, it returns std::make_tuple(0, 0, begin)
     */
    template <typename RandomAccessIt>
    std::tuple<edge_index, edge_index, RandomAccessIt>
    index_range(RandomAccessIt begin, RandomAccessIt end) const {
        static_assert(std::is_same_v<TAlphabet&, decltype(*begin)>
                        || std::is_same_v<const TAlphabet&, decltype(*begin)>,
                      "Only encoded sequences can be queried");

        assert(end >= begin);
        assert(end <= begin + k_);

        if (begin == end)
            return std::make_tuple(edge_index(1), edge_index(1), begin);

        // check if all characters belong to the alphabet
        if (std::any_of(begin, end, [&](TAlphabet c) { return c >= alph_size; }))
            return std::make_tuple(edge_index(0), edge_index(0), begin);

        // get first
        TAlphabet s = *begin;

        // initial range
        edge_index rl = F_.at(s) + 1 < W_->size()
                        ? succ_last(F_.at(s) + 1)
                        : W_->size(); // lower bound
        edge_index ru = s < F_.size() - 1
                        ? F_.at(s + 1)
                        : W_->size() - 1; // upper bound
        if (rl > ru)
            return std::make_tuple(edge_index(0), edge_index(0), begin);

        auto it = begin + 1;
        edge_index rl_h;
        // update range iteratively while scanning through s
        for (; it != end; ++it) {
            s = *it;

            // Include the head of the first node with the given suffix.
            rl_h = pred_last(rl - 1) + 1;

            // Tighten the range including all edges where
            // the source nodes have the given suffix.
            uint64_t rk_rl = rank_W(rl_h - 1, s) + 1;
            uint64_t rk_ru = rank_W(ru, s);
            if (rk_rl > rk_ru)
                return std::make_tuple(rl, ru, it);

            uint64_t offset = rank_last(F_[s]);

            // select the index of the position in last that is rank many positions after offset
            ru = select_last(offset + rk_ru);
            rl = select_last(offset + rk_rl);
        }
        assert(rl <= ru);
        return std::make_tuple(rl, ru, it);
    }

    /**
     * The size of the alphabet for kmers that this graph encodes. For DNA, this value is
     * 5 (A,C,G,T and $)
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

    const KmerExtractorBOSS kmer_extractor_;
    const size_t bits_per_char_W_;

    // k-mer size
    size_t k_;
    // the bit array indicating the last outgoing edge of a node
    bit_vector *last_;
    // the offset array to mark the offsets for the last column in the implicit node list
    std::vector<uint64_t> F_;
    // the array containing the edge labels
    wavelet_tree *W_;

    Config::StateType state = Config::DYN;

    /**
     * This function gets a value of the alphabet c and updates the offset of
     * all following values by adding |value|.
     */
    void update_F(TAlphabet c, int value);

    /**
     * Given a character c and an edge index, this function
     * creates an outgoing edge from the same source node with
     * label c if it is not a part of the graph yet.
     */
    edge_index append_pos(TAlphabet c, edge_index source_node,
                          const TAlphabet *source_node_kmer,
                          std::vector<uint64_t> *edges_inserted = NULL);

    /**
     * Helper function used by the append_pos function
     */
    uint64_t insert_edge(TAlphabet c, uint64_t begin, uint64_t end);

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
     * Given a (k+1)-mer, this function returns the index
     * of the corresponding edge, if such exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    edge_index map_to_edge(RandomAccessIt begin, RandomAccessIt end) const {
        assert(begin + k_ + 1 == end);

        uint64_t edge = index(begin, end - 1);

        return edge ? pick_edge(edge, *(end - 1)) : npos;
    }

    /**
     * Given a k-mer, this function returns the index of last edge going out
     * from the k-mer's corresponding node, if such a node exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    uint64_t index(RandomAccessIt begin, RandomAccessIt end) const {
        static_assert(std::is_same_v<TAlphabet&, decltype(*begin)>
                        || std::is_same_v<const TAlphabet&, decltype(*begin)>,
                      "Only encoded sequences can be queried");
        assert(begin + k_ == end);

        auto match = index_range(begin, end);

        return std::get<2>(match) == end ? std::get<1>(match) : 0;
    }

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
};

std::ostream& operator<<(std::ostream &os, const BOSS &graph);


#endif // __BOSS_HPP__
