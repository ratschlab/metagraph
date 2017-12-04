#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

#include <zlib.h>
#include <unordered_map>
#include <type_traits>

#include "datatypes.hpp"
#include "kmer.hpp"
#include "config.hpp"
#include "dbg_succinct_merge.hpp"


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
    friend void merge::merge(DBG_succ *Gt, DBG_succ *Gm);

    friend void merge::merge(DBG_succ *Gt,
                             std::vector<DBG_succ*> Gv,
                             std::vector<uint64_t> kv,
                             std::vector<uint64_t> nv);
  public:
    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef uint64_t TAlphabet;

    DBG_succ(size_t k = 1, bool sentinel = true);
    ~DBG_succ();

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
     */
    void add_sequence(const std::string &seq, bool append = true);

    void add_sequence_fast(const std::string &seq,
                           bool add_bridge = true,
                           unsigned int parallel = 1,
                           std::string suffix = "");

    void construct_succ(unsigned int parallel = 1);

    /** This function takes a pointer to a graph structure and concatenates the arrays W, last
     * and F to this graph's arrays. In almost all cases this will not produce a valid graph and
     * should only be used as a helper in the parallel merge procedure.
     */
    void append_graph(DBG_succ *other);

    /**
     * This function takes a pointer to a graph structure and concatenates the arrays W, last
     * and F to this graph's static containers last_stat and W_stat. In almost all cases
     * this will not produce a valid graph and should only be used as a helper in the
     * parallel merge procedure.
     */
    void append_graph_static(DBG_succ *other);

    uint64_t remove_edges(std::set<uint64_t> &edges, uint64_t ref_point = 0);

    void add_sink(unsigned int parallel = 1, std::string suffix = "");

    /**
     * This is a debug function that prints the current representation of the graph to
     * the screen.
     */
    void print_state() const;

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

    /**
     * Given index i of a node and a value k, this function
     * will return the k-th last character of node i.
     */
    std::pair<TAlphabet, uint64_t> get_minus_k_value(uint64_t i, uint64_t k) const;

    //Temporary storage for kmers before succinct representation construction
    //the second element stores an ID for each kmer indicating which sequence it came from
    //an even ID indicates that it's a normal kmer, an odd ID indicates that it's a bridge
    std::vector<KMer> kmers;

    // the bit array indicating the last outgoing edge of a node
    bit_vector *last = new bit_vector_dyn();

    // the array containing the edge labels
    wavelet_tree *W = new wavelet_tree_dyn(4); // 4 is log (sigma)

    // the bit array indicating the last outgoing edge of a node (static container for full init)
    std::vector<bool> last_stat;
    std::vector<uint8_t> last_stat_safe; // need uint8 here for thread safety

    // the array containing the edge labels
    std::vector<uint8_t> W_stat;

    //read coverage
    std::vector<size_t> coverage;

    //read bridge indicator
    //std::vector<uint8_t> bridge_stat;
    //bit_vector *bridge = new bit_vector_dyn();

    // the offset array to mark the offsets for the last column in the implicit node list
    std::vector<TAlphabet> F;

    // k-mer size
    size_t k_;

    // index of position that marks end in graph
    uint64_t p_;

    static const std::string alphabet;
    static const size_t alph_size;

    // state of graph
    Config::StateType state = Config::STAT;

    //default values for the sink node
    std::string start;
    std::string sink;

#ifdef DBGDEBUG
    bool verbose = true;
#else
    bool verbose = false;
#endif

    /**
     * This is a convenience function that returns for array W, a position i and
     * a character c the first index of a character c in W[i..N].
     */
    uint64_t succ_W(uint64_t i, TAlphabet c) const;

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
     * Using the offset structure F this function returns the value of the last
     * position of node i.
     */
    TAlphabet get_node_end_value(uint64_t i) const;

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

    /**
     * This function gets a value of the alphabet c and updates the offset of
     * all following values by +1 is positive is true and by -1 otherwise.
     */
    void update_F(TAlphabet c, bool positive);

    /**
     * This is a convenience function to replace the value at
     * position i in W with val.
     */
    void W_set_value(size_t i, TAlphabet val);

    void switch_state(Config::StateType state);

  private:
    /**
     * This function takes a character c and appends it to the end of the graph
     * sequence given that the corresponding note is not part of the graph yet.
     */
    uint64_t append_pos(uint64_t c, uint64_t *ckmer = NULL, uint64_t i = 0);

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
     * Given a position i, this function returns the boundaries of the interval
     * of nodes identical to node i (ignoring the values in W).
     */
    std::pair<uint64_t, uint64_t> get_equal_node_range(uint64_t i) const;

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
    template <class T>
    uint64_t index(const T &str, uint64_t length) const {
        auto range = index_range(str, length);
        if (range.first == 0 && range.second == 0) {
            return 0;
        } else {
            return range.second > range.first
                                    ? range.second
                                    : range.first;
        }
    }

    std::vector<HitInfo> index_fuzzy(const std::string &str, uint64_t eops) const;

    /**
     * Given a node label str, this function returns the index range
     * of nodes sharing the suffix str, if no such range exists, the pair
     (0, 0) is returnd.
     */
    template <class T>
    std::pair<uint64_t, uint64_t> index_range(const T &str, uint64_t length) const {

        // get first
        TAlphabet s;
        if (std::is_same<T, std::string>::value) {
            s = encode(str[0]);
        } else {
            s = str[0];
        }
        // init range
        uint64_t rl = succ_last(F.at(s) + 1);
        uint64_t ru = (s < F.size() - 1) ? F.at(s + 1) : (W->size() - 1); // upper bound
        uint64_t pl;
        //fprintf(stderr, "char: %i rl: %i ru: %i\n", (int) s, (int) rl, (int) ru);
        // update range iteratively while scanning through s
        for (uint64_t i = 1; i < length; ++i) {
            if (std::is_same<T, std::string>::value) {
                s = encode(str[i]);
            } else {
                s = str[i];
            }
            pl = pred_last(rl - 1) + 1;
            rl = std::min(succ_W(pl, s), succ_W(pl, s + alph_size));
            if (rl >= W->size())
                return std::make_pair(0, 0);
            ru = std::max(pred_W(ru, s), pred_W(ru, s + alph_size));
            if (ru >= W->size())
                return std::make_pair(0, 0);
            if (rl > ru)
                return std::make_pair(0, 0);

            rl = outgoing(rl, s);
            ru = outgoing(ru, s);
        }
        return std::make_pair(rl, ru);
    }

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

};

#endif // __DBG_SUCCINCT_HPP__
