#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include <type_traits>

#include "datatypes.hpp"
#include "kmer.hpp"
#include "config.hpp"
#include "dbg_succinct_merge.hpp"
#include "wavelet_tree.hpp"
#include "bit_vector.hpp"


/*
GenomeGraph {
  public:
    template <class T>
    get_index(const T &k_mer)
    get_indices(const std::string &sequence)

    get_incoming_edges(node_it node, char edge_label)
    get_outgoing_edges(node_it node, char edge_label)

    void add_sequence(const std::string &sequence)
};


class DBGCondensed : public GenomeGraph {

};
*/


class DBG_succ { //: public GenomeGraph{
    friend DBG_succ* merge::build_chunk(const std::vector<const DBG_succ*> &graphs,
                                        Config *config);

    friend DBG_succ* merge::merge_chunks(const std::string &filenamebase,
                                         size_t num_chunks);

    friend void merge::merge(DBG_succ *Gt, const DBG_succ &Gm);

    friend DBG_succ* merge::merge(const std::vector<const DBG_succ*> &Gv,
                                  std::vector<uint64_t> kv,
                                  std::vector<uint64_t> nv);
  public:
    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef uint64_t TAlphabet;

    DBG_succ(size_t k = 1, bool sentinel = true);
    ~DBG_succ();

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
     * This functions gets a position i reflecting the r-th occurence of the corresponding
     * character c in W and returns the position of the r-th occurence of c in last.
     */
    uint64_t fwd(uint64_t i) const;

    /**
    * Perform an element wise comparison of the arrays W, last and
    * F and will only check for identity. If any element differs, return
    * false and true otherwise.
    */
    bool operator==(const DBG_succ &other) const;

    bool load(const std::string &filename_base);
    void serialize(const std::string &filename_base) const;

    /**
     * Add a full sequence to the graph.
     * If |try_extend| is true, search for the first k-mer in the graph
     * and extend it from that point. If the search fails, start from the dummy source.
     */
    void add_sequence(const std::string &seq, bool try_extend = false);

    void add_sequence_fast(const std::string &seq,
                           bool add_bridge = true,
                           unsigned int parallel = 1,
                           std::string suffix = "");

    void construct_succ(unsigned int parallel = 1);

    void remove_edges(const std::set<uint64_t> &edges);

    /**
     * Given index i of a node and a value k, this function
     * will return the k-th last character of node i.
     */
    std::pair<TAlphabet, uint64_t> get_minus_k_value(uint64_t i, uint64_t k) const;

    // the bit array indicating the last outgoing edge of a node
    bit_vector *last = new bit_vector_dyn();

    // the array containing the edge labels
    wavelet_tree *W = new wavelet_tree_dyn(4); // 4 is log (sigma)

    //Temporary storage for kmers before succinct representation construction
    //the second element stores an ID for each kmer indicating which sequence it came from
    //an even ID indicates that it's a normal kmer, an odd ID indicates that it's a bridge
    std::vector<KMer> kmers;

    static const std::string alphabet;
    static const size_t alph_size;

    Config::StateType state = Config::DYN;

#ifdef DBGDEBUG
    bool verbose = true;
#else
    bool verbose = false;
#endif

    /**
     * Using the offset structure F this function returns the value of the last
     * position of node i.
     */
    TAlphabet get_node_last_char(uint64_t i) const;

    /**
     * Given a node index i, this function returns the number of outgoing
     * edges from node i.
     */
    uint64_t outdegree(uint64_t i) const;

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
     * Given a node label s, this function returns the index
     * of the corresponding node or the closest predecessor, if no node
     * with the sequence is not found.
     */
    uint64_t colex_upper_bound(std::deque<TAlphabet> str) const;

    /**
    * Given a node index k_node, this function returns the k-mer sequence of the
    * node in a deque data structure.
    */
    std::deque<TAlphabet> get_node_seq(uint64_t k_node) const;

    /**
     * Return k-mer length of current graph.
     */
    uint64_t get_k() const { return this->k_; }

    /**
     * Return value of W at position k.
     */
    TAlphabet get_W(uint64_t k) const { return (*W)[k]; }
    const wavelet_tree& get_W() const { return *W; }

    // Given the alphabet index return the corresponding symbol
    static char decode(TAlphabet s);
    // Given the alphabet character return its corresponding number
    static TAlphabet encode(char s);

    /**
     * Breaks the seq into k-mers and searches for the index of each
     * k-mer in the graph. Returns these indices.
     */
    std::vector<uint64_t> align(const std::string &sequence,
                                uint64_t alignment_length = 0) const;

    std::vector<std::vector<HitInfo>> align_fuzzy(const std::string &sequence,
                                                  uint64_t max_distance = 0,
                                                  uint64_t alignment_length = 0) const;

    uint64_t num_nodes() const;

    uint64_t num_edges() const;

    void switch_state(Config::StateType state);

    /**
     * This is a debug function that prints the current representation of the graph to
     * the screen.
     */
    void print_state(std::ostream &os) const;

    /**
     * Returns the sequence stored in W and prints the node
     * information in an overview.
     * Useful for debugging purposes.
     */
    void print_seq() const;

    /**
     * Write the adjacency list to file |filename| or
     * print it to stdout of the filename is not provided.
     */
    void print_adj_list(const std::string &filename = "") const;

  private:
    // k-mer size
    size_t k_;

    // index of position that marks the source of the graph
    const uint64_t kDummySource = 1;

    // the offset array to mark the offsets for the last column in the implicit node list
    std::vector<uint64_t> F;

    /**
     * This function gets a value of the alphabet c and updates the offset of
     * all following values by adding |value|.
     */
    void update_F(TAlphabet c, int value);

    /**
     * This is a convenience function that returns for array W, a position i and
     * a character c the first index of a character c in W[i..N].
     */
    uint64_t succ_W(uint64_t i, TAlphabet c) const;

    /**
     * This function takes a character c and appends it to the end of the graph
     * sequence given that the corresponding note is not part of the graph yet.
     */
    uint64_t append_pos(uint64_t c, uint64_t source_node, TAlphabet *ckmer = NULL);

    /**
     * Helper function used by the append_pos function
     */
    void insert_edge(TAlphabet c, uint64_t begin, uint64_t end);

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
    void sort_W_locally(uint64_t l, uint64_t u);

    /**
     * Return value of last at position k.
     */
    bool get_last(uint64_t k) const { return (*last)[k]; }

    /**
     * Return value of F vector at index k.
     * The index is over the alphabet!
     */
    TAlphabet get_F(size_t k) const { return F.at(k); }

    /**
    * Given a node index k, this function returns the k-mer sequence of the
    * node as a string.
    */
    std::string get_node_str(uint64_t k) const;

    /**
     * This function returns true if node i is a terminal node.
     */
    bool is_terminal_node(uint64_t i) const;

    /**
     * This function gets two node indices and returns if the
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(uint64_t i1, uint64_t i2) const;
    bool compare_node_suffix(TAlphabet *ref, uint64_t i2) const;

    /**
     * Given a node label s, this function returns the index
     * of the corresponding node, if this node exists and 0 otherwise.
     */
    template <typename T>
    uint64_t index(const T &k_mer) const {
        static_assert(std::is_base_of<std::string, T>::value
                        || std::is_base_of<std::deque<TAlphabet>, T>::value
                        || std::is_base_of<std::vector<TAlphabet>, T>::value);
        assert(k_mer.size() >= k_);

        auto range = index_range(k_mer, k_);
        return std::max(range.first, range.second);
    }

    /**
     * Given a node label str, this function returns the index range
     * of nodes sharing the suffix str, if no such range exists, the pair
     (0, 0) is returnd.
     */
    template <typename T>
    std::pair<uint64_t, uint64_t> index_range(const T &str, uint64_t length) const {
        static_assert(std::is_base_of<std::string, T>::value
                        || std::is_base_of<std::deque<TAlphabet>, T>::value
                        || std::is_base_of<std::vector<TAlphabet>, T>::value);
        // get first
        TAlphabet s = std::is_same<T, std::string>::value
                        ? encode(str[0])
                        : str[0];

        // initial range
        uint64_t rl = succ_last(F.at(s) + 1);
        uint64_t ru = s < F.size() - 1
                            ? F.at(s + 1)
                            : W->size() - 1; // upper bound
        if (rl > ru)
            return std::make_pair(0, 0);

        // update range iteratively while scanning through s
        for (uint64_t i = 1; i < length; ++i) {
            s = std::is_same<T, std::string>::value
                        ? encode(str[i])
                        : str[i];

            // Include the head of the first node with the given suffix.
            rl = pred_last(rl - 1) + 1;

            // Tighten the range including all edges where
            // the source nodes have the given suffix.
            rl = std::min(succ_W(rl, s),
                          succ_W(rl, s + alph_size));
            if (rl >= W->size())
                return std::make_pair(0, 0);

            ru = std::max(pred_W(ru, s),
                          pred_W(ru, s + alph_size));
            if (rl > ru)
                return std::make_pair(0, 0);

            // Translate the node indices from the sources to the targets.
            rl = outgoing(rl, s);
            ru = outgoing(ru, s);
        }
        return std::make_pair(rl, ru);
    }

    std::vector<HitInfo> index_fuzzy(const std::string &str, uint64_t eops) const;

    /**
     * Given a position i in W and an edge label c, this function returns the
     * index of the outgoing edge with label c.
     */
    uint64_t outgoing_edge_idx(uint64_t i, TAlphabet c) const;

    /**
     * Given index of node i, the function returns the
     * first character of the node.
     */
    TAlphabet get_node_begin_value(uint64_t i) const;

    /**
     * This function gets a position i that reflects the i-th node and returns the
     * position in W that corresponds to the i-th node's last character.
     */
    uint64_t bwd(uint64_t i) const;

    /**
     * This is a convenience function that returns for array W, a position i and
     * a character c the last index of a character c preceding in W[1..i].
     */
    uint64_t pred_W(uint64_t i, TAlphabet c) const;

    void verbose_cout() const {}

    template <typename T, typename... Targs>
    void verbose_cout(const T &arg, Targs ...rest) const {
        if (verbose) {
            std::cout << arg;
            verbose_cout(rest...);
        }
    }
};

inline std::ostream& operator<<(std::ostream &os,
                                const DBG_succ &graph) {
    graph.print_state(os);
    return os;
}

#endif // __DBG_SUCCINCT_HPP__
