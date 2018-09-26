#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include <type_traits>

#include "datatypes.hpp"
#include "config.hpp"
#include "wavelet_tree.hpp"
#include "bit_vector.hpp"

typedef uint8_t TAlphabet;


class SequenceGraph {
  public:
    typedef uint64_t node_index;
    typedef uint64_t edge_index;
    static const uint64_t npos;

    virtual ~SequenceGraph() {};

    virtual void add_sequence(const std::string &sequence,
                              bool try_extend = false,
                              bit_vector_dyn *edges_inserted = NULL) = 0;

    // Traverse graph mapping k-mers from sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    virtual void map_to_nodes(const std::string &sequence,
                              const std::function<void(node_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Traverse graph mapping k-mers from sequence to the graph edges
    // and run callback for each edge until the termination condition is satisfied
    virtual void map_to_edges(const std::string &sequence,
                              const std::function<void(edge_index)> &callback,
                              const std::function<bool()> &terminate = [](){ return false; }) const = 0;

    // Check whether graph contains fraction of nodes from the sequence
    virtual bool find(const std::string &sequence,
                      double discovery_fraction = 1) const = 0;

    // Traverse the outgoing edge
    virtual node_index traverse(node_index node, char edge_label) const = 0;
    // Traverse the incoming edge
    virtual node_index traverse_back(node_index node, char edge_label) const = 0;

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;
};


/*
class DBGCondensed : public SequenceGraph {

};
*/


class DBGSuccConstructor;

class DBG_succ : public SequenceGraph {
  public:
    // in [1,...,num_nodes], 0 = npos (invalid index)
    typedef uint64_t node_index;
    // in [1,...,num_edges], 0 = npos (invalid index)
    typedef uint64_t edge_index;

    explicit DBG_succ(size_t k = 1);
    explicit DBG_succ(DBGSuccConstructor *builder);
    ~DBG_succ();

    /**
     * Check whether graphs store the same data.
     * FYI: this function reconstructs all the kmers, and
     * the complexity is at least O(k x n).
     */
    bool operator==(const DBG_succ &other) const;

    /**
     * Perform an element wise comparison of the arrays W, last and
     * F and will only check for identity. If any element differs, return
     * false and true otherwise.
     */
    bool equals_internally(const DBG_succ &other, bool verbose = false) const;

    // Return the k-mer length
    size_t get_k() const { return k_; }
    uint64_t num_nodes() const;
    uint64_t num_edges() const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    // Traverse graph mapping k-mers from sequence to the graph nodes
    // and run callback for each node until the termination condition is satisfied
    // Call npos if a k-mer can't be mapped to the graph nodes
    void map_to_nodes(const std::string &sequence,
                      const std::function<void(node_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    /**
     * Breaks the sequence into k-mers and searches for the index of each
     * k-mer in the graph. Returns these indices.
     * Default: kmer_size = k
     */
    // TODO: remove parameter kmer_size and implement the arbitrary k-mer size
    // functionality in align_fuzzy()
    std::vector<node_index> map_to_nodes(const std::string &sequence,
                                         size_t kmer_size = 0) const;

    // Traverse graph mapping k-mers from sequence to the graph edges
    // and run callback for each edge until the termination condition is satisfied
    // Call npos if a k-mer can't be mapped to the graph edges
    void map_to_edges(const std::string &sequence,
                      const std::function<void(edge_index)> &callback,
                      const std::function<bool()> &terminate = [](){ return false; }) const;

    std::vector<edge_index> map_to_edges(const std::string &sequence) const;

    // Check whether the graph contains a fraction of k-mers from the sequence
    bool find(const std::string &sequence,
              double kmer_discovery_fraction = 1) const;

    bool find(const std::string &sequence,
              double kmer_discovery_fraction,
              size_t kmer_mapping_mode) const;

    // TODO: revise the implementation, write unit tests
    std::vector<std::vector<HitInfo>> align_fuzzy(const std::string &sequence,
                                                  size_t max_distance = 0,
                                                  size_t alignment_length = 0) const;

    template <class... T>
    using Call = typename std::function<void(T...)>;

    // traverse all nodes in graph except for the dummy source of sink ones
    void call_kmers(Call<node_index, const std::string&> callback) const;
    void call_edges(Call<edge_index, const std::vector<TAlphabet>&> callback) const;
    // call paths (or simple paths if |split_to_contigs| is true) that cover
    // exactly all edges in graph
    void call_paths(Call<const std::vector<edge_index>,
                         const std::vector<TAlphabet>&> callback,
                    bool split_to_contigs = false) const;
    void call_sequences(Call<const std::string&> callback,
                        bool split_to_contigs = false) const;

    node_index traverse(node_index node, char edge_label) const;
    node_index traverse_back(node_index node, char edge_label) const;

    /**
     * Add a full sequence to the graph.
     * If |try_extend| is true, search for the first k-mer in the graph
     * and extend it from that point. If the search fails, start from the dummy source.
     */
    void add_sequence(const std::string &seq,
                      bool try_extend = false,
                      bit_vector_dyn *edges_inserted = NULL);

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
    std::vector<bool>
    erase_redundant_dummy_edges(std::vector<bool> *source_dummy_edges = NULL,
                                size_t num_threads = 0,
                                bool verbose = false);

    uint64_t mark_source_dummy_edges(std::vector<bool> *mask,
                                     size_t num_threads = 0,
                                     bool verbose = false) const;

    uint64_t mark_sink_dummy_edges(std::vector<bool> *mask) const;

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
    void merge(const DBG_succ &other);

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
     * Given a node index i, this function returns the number of outgoing
     * edges from the node i.
     */
    size_t outdegree(node_index i) const;

    /**
     * Given an edge index i, this function returns true if that is
     * the only outgoing edge from its source node.
     */
    bool is_single_outgoing(edge_index i) const;

    /**
     * Given a node index i, this function returns the number of incoming
     * edges to node i.
     */
    size_t indegree(node_index i) const;

    /**
     * Given an edge index i, this function returns true if that is
     * the only edge incoming to its target node.
     */
    bool is_single_incoming(edge_index i) const;

    /**
     * Given a node index i and an edge label c, this function returns the
     * index of the node the edge is pointing to.
     */
    node_index outgoing(node_index i, TAlphabet c) const;

    /**
     * Given a node index i and an edge label c, this function returns the
     * index of the node the incoming edge belongs to.
     */
    node_index incoming(node_index i, TAlphabet c) const;

    /**
     * Given a node index i and an edge label c, this function returns the
     * index of the outgoing edge with label c if it exists and npos otherwise.
     */
    edge_index outgoing_edge_idx(node_index i, TAlphabet c) const;

    /**
     * Given an edge index i and a character c, get the index of the edge with
     * label c outgoing from the same source node if such exists and npos otherwise.
     */
    edge_index pick_edge(edge_index edge, node_index node, TAlphabet c) const;

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
     * This is a convenience function that returns for array W, a position i and
     * a character c the last index of a character c preceding in W[1..i].
     */
    uint64_t pred_W(uint64_t i, TAlphabet c) const;

    /**
     * This is a convenience function that returns for array W, a position i and
     * a character c the first index of a character c in W[i..N].
     */
    uint64_t succ_W(uint64_t i, TAlphabet c) const;

    /**
     * Return value of F vector at index k.
     * The index is over the alphabet!
     */
    uint64_t get_F(TAlphabet k) const { return F_.at(k); }

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
    static char decode(TAlphabet s);
    static std::string decode(const std::vector<TAlphabet> &sequence);
    // Given the alphabet character return its corresponding number
    static TAlphabet encode(char s);
    static std::vector<TAlphabet> encode(const std::string &sequence);

    static const std::string alphabet;
    static const TAlphabet alph_size;
    static const size_t kLogSigma;

    static const size_t kSentinelCode = 0;
    static const size_t kSentinel = '$';

  private:
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
                          bit_vector_dyn *edges_inserted = NULL);

    /**
     * Helper function used by the append_pos function
     */
    bool insert_edge(TAlphabet c, uint64_t begin, uint64_t end,
                     bit_vector_dyn *edges_inserted = NULL);

    /**
     * Erase exactly all the masked edges from the graph,
     * may invalidate the graph (if leaves nodes with no incoming edges).
     * Returns the number of edges erased.
     */
    uint64_t erase_edges(const std::vector<bool> &edges_to_remove_mask);

    /**
     * This function gets two edge indices and returns if their source
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(edge_index first, edge_index second) const;
    bool compare_node_suffix(edge_index first, const TAlphabet *second) const;

    inline node_index get_source_node(edge_index i) const;

    /**
     * Given a k-mer, this function returns the index
     * of the corresponding node, if such exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    node_index map_to_node(RandomAccessIt begin, RandomAccessIt end) const {
        uint64_t edge = index(begin, end);
        return edge ? get_source_node(edge) : npos;
    }

    /**
     * Given a (k+1)-mer, this function returns the index
     * of the corresponding edge, if such exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    node_index map_to_edge(RandomAccessIt begin, RandomAccessIt end) const {
        assert(begin + k_ + 1 == end);

        uint64_t edge = index(begin, end - 1);

        return edge ? pick_edge(edge, get_source_node(edge), *(end - 1))
                    : npos;
    }

    /**
     * Given a k-mer, this function returns the index of last edge going out
     * from the k-mer's corresponding node, if such a node exists and 0 otherwise.
     */
    template <typename RandomAccessIt>
    uint64_t index(RandomAccessIt begin, RandomAccessIt end) const {
        static_assert(std::is_same<TAlphabet&, decltype(*begin)>::value,
                      "Only encoded sequences can be queried");

        assert(begin + k_ == end);

        return index_range(begin, end).second;
    }

    /**
     * Given a node label str, this function returns the index range
     * of nodes sharing the suffix str, if no such range exists, the pair
     * (0, 0) is returnd.
     */
    template <typename RandomAccessIt>
    std::pair<uint64_t, uint64_t> index_range(RandomAccessIt begin,
                                              RandomAccessIt end) const {
        static_assert(std::is_same<TAlphabet&, decltype(*begin)>::value,
                      "Only encoded sequences can be queried");

        assert(end > begin);
        assert(end <= begin + k_);

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
            return std::make_pair(0, 0);

        // update range iteratively while scanning through s
        for (auto it = ++begin; it != end; ++it) {
            s = *it;

            // Include the head of the first node with the given suffix.
            rl = pred_last(rl - 1) + 1;

            // Tighten the range including all edges where
            // the source nodes have the given suffix.
            uint64_t rk_rl = rank_W(rl - 1, s) + 1;
            uint64_t rk_ru = rank_W(ru, s);
            if (rk_rl > rk_ru)
                return std::make_pair(0, 0);

            uint64_t offset = rank_last(F_[s]);

            // select the index of the position in last that is rank many positions after offset
            ru = select_last(offset + rk_ru);
            rl = select_last(offset + rk_rl);
        }
        assert(rl <= ru);
        return std::make_pair(rl, ru);
    }

    // TODO: revise the implementation, write unit tests
    std::vector<HitInfo> index_fuzzy(const std::string &str,
                                     size_t max_distance) const;

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

std::ostream& operator<<(std::ostream &os, const DBG_succ &graph);

#endif // __DBG_SUCCINCT_HPP__
