#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include <type_traits>

#include "datatypes.hpp"
#include "config.hpp"
#include "wavelet_tree.hpp"
#include "bit_vector.hpp"

class DBGSuccConstructor;


typedef uint8_t TAlphabet;


class SequenceGraph {
  public:
    typedef uint64_t node_index;
    typedef uint64_t edge_index;

    static const node_index npos;

    virtual ~SequenceGraph() {};

    virtual void add_sequence(const std::string &sequence,
                              bool try_extend = false,
                              bit_vector_dyn *edges_inserted = NULL) = 0;

    // Traverse graph aligning the sequence
    // and run callback until the termination condition is satisfied
    virtual void align(const std::string &sequence,
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


class DBG_succ : public SequenceGraph {
  public:
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
    bool equals_internally(const DBG_succ &other) const;

    // Return the k-mer length
    uint64_t get_k() const { return this->k_; }
    uint64_t num_nodes() const;
    uint64_t num_edges() const;

    node_index traverse(node_index node, char edge_label) const;
    node_index traverse_back(node_index node, char edge_label) const;

    // Traverse graph aligning k-mers to the sequence
    // and run callback until the termination condition is satisfied
    void align(const std::string &sequence,
               const std::function<void(edge_index)> &callback,
               const std::function<bool()> &terminate = [](){ return false; }) const;

    // Check whether graph contains fraction of k-mers from the sequence
    bool find(const std::string &sequence,
              double kmer_discovery_fraction = 1) const;

    bool find(const std::string &sequence,
              double kmer_discovery_fraction,
              size_t kmer_mapping_mode) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    /**
     * Add a full sequence to the graph.
     * If |try_extend| is true, search for the first k-mer in the graph
     * and extend it from that point. If the search fails, start from the dummy source.
     */
    void add_sequence(const std::string &seq,
                      bool try_extend = false,
                      bit_vector_dyn *edges_inserted = NULL);

    // Given an edge list, remove them from the graph.
    // TODO: fix the implementation
    void remove_edges(const std::set<edge_index> &edges);

    /**
    * Heavily borrowing from the graph sequence traversal, this function gets 
    * a graph `other` and merges its nodes into the target graph object.
    * The edges of `other` are fully traversed and nodes are added if not existing yet.
    * This function is well suited to merge small graphs into large ones.
    */
    void merge(const DBG_succ &other);

    /**
     * Return value of W at position k.
     */
    TAlphabet get_W(uint64_t k) const { return (*W_)[k]; }
    const wavelet_tree& get_W() const { return *W_; }

    /**
     * Return value of last at position k.
     */
    bool get_last(uint64_t k) const { return (*last_)[k]; }

    /**
     * Using the offset structure F this function returns the value of the last
     * position of node i.
     */
    TAlphabet get_node_last_value(uint64_t i) const;

    /**
     * Given a node index k, this function returns the k-mer sequence of the
     * node as a string.
     */
    std::string get_node_str(uint64_t k) const;

    /**
     * Given index i of a node and a value k, this function
     * will return the k-th last character of node i.
     */
    std::pair<TAlphabet, uint64_t> get_minus_k_value(uint64_t i, uint64_t k) const;

    /**
     * Given a node index k_node, this function returns the k-mer sequence of the
     * node.
     */
    std::vector<TAlphabet> get_node_seq(uint64_t k_node) const;

    /**
     * Breaks the seq into k-mers and searches for the index of each
     * k-mer in the graph. Returns these indices.
     * Default: kmer_size = k
     */
    std::vector<uint64_t> map_kmers(const std::string &sequence,
                                    size_t kmer_size = 0) const;

    std::vector<std::vector<HitInfo>> align_fuzzy(const std::string &sequence,
                                                  size_t max_distance = 0,
                                                  size_t alignment_length = 0) const;

    typedef std::function<void(const std::vector<uint64_t>,
                               const std::vector<TAlphabet>&)> PathCallback;

    void call_paths(const PathCallback &callback) const;

    typedef std::function<void(const std::string&)> SequenceCallback;

    void call_sequences(const SequenceCallback &callback) const;

    typedef std::function<void(edge_index, const std::vector<TAlphabet>&)> EdgeCallback;

    void call_edges(const EdgeCallback &callback) const;

    typedef std::function<void(node_index, const std::string&)> KmerCallback;
    // traverse all nodes in graph except for the dummy source of sink ones
    void call_kmers(const KmerCallback &callback) const;

    /**
     * This is a debug function that prints the current representation of the graph to
     * the screen.
     */
    void print_state(std::ostream &os = std::cout) const;

    Config::StateType get_state() const { return state; }

    /**
     * Write the adjacency list to file |filename| or
     * print it to stdout of the filename is not provided.
     */
    void print_adj_list(std::ostream &os = std::cout) const;

    void switch_state(Config::StateType state);

    /**
     * Given a node index i, this function returns the number of outgoing
     * edges from node i.
     */
    uint64_t outdegree(uint64_t i) const;

    /**
     * Given a node index i, this function returns true if there is
     * only one outgoing edges from node i.
     */
    bool is_single_outgoing(uint64_t i) const;

    /**
     * Given a node index i, this function returns the number of incoming
     * edges to node i.
     */
    uint64_t indegree(uint64_t i) const;

    /**
     * Given a position i in W and an edge label c, this function returns the
     * index of the node the edge is pointing to.
     */
    uint64_t outgoing(uint64_t i, TAlphabet c) const;

    /**
     * Given a node index i and an edge label c, this function returns the
     * index of the node the incoming edge belongs to.
     */
    uint64_t incoming(uint64_t i, TAlphabet c) const;

    /**
     * Given a node label kmer, this function returns the index
     * of the corresponding node or the closest predecessor, if no node
     * with the sequence is not found.
     */
    uint64_t pred_kmer(const std::vector<TAlphabet> &kmer) const;

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
     * This function takes a character c and appends it to the end of the graph
     * sequence given that the corresponding note is not part of the graph yet.
     */
    uint64_t append_pos(TAlphabet c, uint64_t source_node, TAlphabet *ckmer = NULL,
                        bit_vector_dyn *edges_inserted = NULL);

    /**
     * Helper function used by the append_pos function
     */
    bool insert_edge(TAlphabet c, uint64_t begin, uint64_t end,
                     bit_vector_dyn *edges_inserted = NULL);

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
     * This function gets a local range in W from lower bound l
     * to upper bound u and swaps the inserted element to the
     * righ location.
     */
    void fix_order_in_W_local(uint64_t l, uint64_t u);

    /**
     * Return value of F vector at index k.
     * The index is over the alphabet!
     */
    uint64_t get_F(TAlphabet k) const { return F_.at(k); }

    /**
     * This function gets two node indices and returns if the
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(uint64_t first, uint64_t second) const;
    bool compare_node_suffix(TAlphabet *ref, uint64_t j) const;

    /**
     * Given a node label s, this function returns the index
     * of the corresponding node, if this node exists and 0 otherwise.
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

        // get first
        TAlphabet s = *begin;

        // initial range
        uint64_t rl = F_.at(s) + 1 < W_->size()
                      ? succ_last(F_.at(s) + 1)
                      : W_->size();
        uint64_t ru = s < F_.size() - 1
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
            rl = std::min(succ_W(rl, s),
                          succ_W(rl, s + alph_size));
            if (rl >= W_->size())
                return std::make_pair(0, 0);

            ru = std::max(pred_W(ru, s),
                          pred_W(ru, s + alph_size));
            if (rl > ru)
                return std::make_pair(0, 0);

            // Translate the node indices from the sources to the targets.
            rl = outgoing(rl, s);
            ru = outgoing(ru, s);
        }
        assert(rl <= ru);
        return std::make_pair(rl, ru);
    }

    std::vector<HitInfo> index_fuzzy(const std::string &str, uint64_t eops) const;

    /**
     * Given a position i in W and an edge label c, this function returns the
     * index of the outgoing edge with label c.
     */
    uint64_t outgoing_edge_idx(uint64_t i, TAlphabet c) const;

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
