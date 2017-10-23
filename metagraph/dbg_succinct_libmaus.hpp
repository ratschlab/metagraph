#ifndef __DBG_SUCCINCT_LIBM_HPP__
#define __DBG_SUCCINCT_LIBM_HPP__

#include <zlib.h>
#include <unordered_map>

#include "config.hpp"
#include "datatypes.hpp"
#include "dbg_succinct_boost.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include "dbg_succinct_boost.hpp"

#include <type_traits>

typedef boost::multiprecision::uint256_t ui256;

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef uint64_t TAlphabet;

    public:
    //Temporary storage for kmers before succinct representation construction
    //the second element stores an ID for each kmer indicating which sequence it came from
    //an even ID indicates that it's a normal kmer, an odd ID indicates that it's a bridge
    std::vector<kmer_boost::KMer> kmers;
    //std::vector<std::pair<ui256, size_t> > kmers;

    //Annotation hash generator
    AnnotationHash hasher;

    //store the position of the last character position modified in F
    size_t lastlet = 0;

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
    std::vector<uint8_t> bridge_stat;
    bit_vector *bridge = new bit_vector_dyn();

    // the offset array to mark the offsets for the last column in the implicit node list
    std::vector<TAlphabet> F; 

    // k-mer size
    size_t k;
    // index of position that marks end in graph
    uint64_t p;

    // alphabet size
    size_t alph_size = 7;

    // infile base when loaded from file
    std::string infbase;

    // config object
    Config* config;

    // alphabet
    const std::string alphabet; //("$ACGTNX$ACGTNXn"); 

    // state of graph
    enum state_type {stat = 1,
                     dyn,
                     cstr};
    state_type state = stat;

    // annotation containers
    std::deque<uint32_t> annotation; // list that associates each node in the graph with an annotation hash
    std::vector<std::string> id_to_label; // maps the label ID back to the original string
    std::unordered_map<std::string, uint32_t> label_to_id_map; // maps each label string to an integer ID
    std::map<uint32_t, uint32_t> annotation_map; // maps the hash of a combination to the position in the combination vector
    std::vector<uint32_t> combination_vector; // contains all known combinations
    uint64_t combination_count = 0;

    //std::vector<sdsl::rrr_vector<63>* > annotation_full;
    std::vector<sdsl::sd_vector<>* > annotation_full;
    sdsl::bit_vector* annotation_curr = NULL;
    std::vector<sdsl::rrr_vector<63> > colors_to_bits;
    std::vector<std::string> bits_to_labels;
    size_t anno_labels;

#ifdef DBGDEBUG
    bool debug = true;
#else
    bool debug = false;
#endif 

    // construct empty graph instance
    DBG_succ(size_t k_,
             Config* config_, 
             bool sentinel = true);

    // load graph instance from a provided file name base
    DBG_succ(std::string infbase_, 
             Config* config_);

    // destructor
    ~DBG_succ();
    
    //
    //
    // QUERY FUNCTIONS
    //
    //
    
    /** 
     * Uses the object's array W, a given position i in W and a character c
     * from the alphabet and returns the number of occurences of c in W up to
     * position i.
     */
    uint64_t rank_W(uint64_t i, TAlphabet c);

    /**
     * Uses the array W and gets a count i and a character c from 
     * the alphabet and returns the positions of the i-th occurence of c 
     * in W.
     */
    uint64_t select_W(uint64_t i, TAlphabet c);

    /**
     * This is a convenience function that returns for array W, a position i and 
     * a character c the last index of a character c preceding in W[1..i].
     */
    uint64_t pred_W(uint64_t i, TAlphabet c);

    /**
     * This is a convenience function that returns for array W, a position i and 
     * a character c the first index of a character c in W[i..N].
     */
    uint64_t succ_W(uint64_t i, TAlphabet c);

    /** 
     * Uses the object's array last and a position and
     * returns the number of set bits up to that postion.
     */
    uint64_t rank_last(uint64_t i);

    /**
     * Uses the object's array last and a given position i and
     * returns the position of the i-th set bit in last[1..i].
     */
    uint64_t select_last(uint64_t i);

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the last set bit in last[1..i].
     */
    uint64_t pred_last(uint64_t i);

    /**
     * This is a convenience function that returns for the object's array last
     * and a given position i the position of the first set bit in last[i..N].
     */
    uint64_t succ_last(uint64_t i);

    /**
     * This function gets a position i that reflects the i-th node and returns the
     * position in W that corresponds to the i-th node's last character. 
     */
    uint64_t bwd(uint64_t i);

    /**
     * This functions gets a position i reflecting the r-th occurence of the corresponding
     * character c in W and returns the position of the r-th occurence of c in last.
     */
    uint64_t fwd(uint64_t i);

    /**
     * Using the offset structure F this function returns the value of the last 
     * position of node i.
     */
    TAlphabet get_node_end_value(uint64_t i);
    
    /**
     * Given index of node i, the function returns the 
     * first character of the node.
     */
    TAlphabet get_node_begin_value(uint64_t i);

    /**
     * Given a position i in W and an edge label c, this function returns the
     * index of the outgoing edge with label c.
     */
    uint64_t outgoing_edge_idx(uint64_t i, TAlphabet c);

    /**
     * Given a position i in W and an edge label c, this function returns the
     * index of the node the edge is pointing to.
     */
    uint64_t outgoing(uint64_t i, TAlphabet c);

    /**
     * Given a node index i and an edge label c, this function returns the
     * index of the node the incoming edge belongs to.
     */
    uint64_t incoming(uint64_t i, TAlphabet c);

    /**
     * Given a node index i, this function returns the number of outgoing
     * edges from node i.
     */
    uint64_t outdegree(uint64_t i);

    /**
     * Given a node index i, this function returns the number of incoming
     * edges to node i.
     */
    uint64_t indegree(uint64_t i);

    /**
     * Given a node label s, this function returns the index
     * of the corresponding node, if this node exists and 0 otherwise.
     */
    template <class T> uint64_t index(T &s_, uint64_t length) {
        TAlphabet s;
        if (std::is_same<T, std::string>::value) {
            s = get_alphabet_number(s_[0]);
        } else {
            s = s_[0];
        }    
        // init range
        uint64_t rl = succ_last(F[s] + 1);
        uint64_t ru = s+1 < alph_size ? F[s + 1] : last->size()-1; // upper bound
        // update range iteratively while scanning through s
        for (uint64_t i = 1; i < length; i++) {
            if (std::is_same<T, std::string>::value) {
                s = get_alphabet_number(s_[i]);
            } else {
                s = s_[i];
            }    
            rl = std::min(succ_W(pred_last(rl - 1) + 1, s), succ_W(pred_last(rl - 1) + 1, s + alph_size));
            if (rl >= W->size())
                return 0;
            ru = std::max(pred_W(ru, s), pred_W(ru, s + alph_size));
            if (ru >= W->size())
                return 0;
            if (rl > ru)
                return 0;
            rl = outgoing(rl, s);
            ru = outgoing(ru, s);
        }    
        return (ru > rl) ? ru : rl;

    }

    uint64_t index(std::string &s_) {
        return index(s_, s_.length());
    }

    uint64_t index(std::deque<TAlphabet> str);

    std::vector<HitInfo> index_fuzzy(std::string &str, uint64_t eops);

    std::pair<uint64_t, uint64_t> index_range(std::deque<TAlphabet> str);

    uint64_t index_predecessor(std::deque<TAlphabet> str);

    /**
     * Given a position i, this function returns the boundaries of the interval
     * of nodes identical to node i (ignoring the values in W).
     */
    std::pair<uint64_t, uint64_t> get_equal_node_range(uint64_t i);

    /**
     * Given index i of a node and a value k, this function 
     * will return the k-th last character of node i.
     */
    std::pair<TAlphabet, uint64_t> get_minus_k_value(uint64_t i, uint64_t k);

    /** 
     * This function gets two node indices and returns if the
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(uint64_t i1, uint64_t i2);
    bool compare_node_suffix(TAlphabet *ref, uint64_t i2);

    /**
     * This function returns true if node i is a terminal node.
     */
    bool is_terminal_node(uint64_t i);

    /**
    * Given a node index k_node, this function returns the k-mer sequence of the
    * node in a deque data structure.
    */
    std::deque<TAlphabet> get_node_seq(uint64_t k_node); 

    /**
    * Given a node index k_node, this function returns the k-mer sequence of the 
    * node as a string.
    */
    std::string get_node_str(uint64_t k_node); 

    /**
     * Return number of edges in the current graph.
     * TODO: delete
     */
    uint64_t get_size();

    /**
     * Return number of nodes in the current graph.
     */
    uint64_t get_nodes();

    /**
     * Return k-mer length of current graph.
     * TODO: delete
     */
    uint64_t get_k();

    /**
     * Return value of W at position k.
     */
    TAlphabet get_W(uint64_t k);

    /** 
     * Return value of F vector at index k.
     * The index is over the alphabet!
     */
    TAlphabet get_F(uint64_t k);

    /**
     * Return value of last at position k.
     */
    bool get_last(uint64_t k);

    // Given the alphabet index return the corresponding symbol
    char get_alphabet_symbol(uint64_t s);

    // Given the alphabet character return its corresponding number
    TAlphabet get_alphabet_number(char s);

    /**
     * Breaks the seq into k-mers and searches for the index of each
     * k-mer in the graph. Returns these indices.
     */
    std::vector<uint64_t> align(kstring_t seq);

    std::vector<std::vector<HitInfo> > align_fuzzy(kstring_t seq, uint64_t max_distance = 0, uint64_t alignment_length = 0);

    uint64_t get_node_count();

    uint64_t get_edge_count();


    //
    //
    // APPEND
    // 
    //
    
    /**
     * This function gets a value of the alphabet c and updates the offset of 
     * all following values by +1 is positive is true and by -1 otherwise.
     */
    void update_F(TAlphabet c, bool positive);

    /**
     * This function gets a local range in W from lower bound l
     * to upper bound u and swaps the inserted element to the
     * righ location.
     */
    void sort_W_locally(uint64_t l, uint64_t u);

    /** 
     * This is a convenience function to replace the value at
     * position i in W with val.
     */
    void replaceW(size_t i, TAlphabet val);

    void switch_state(state_type state);
    
    //
    //
    // MERGE
    //
    //

   
    uint64_t next_non_zero(std::vector<uint64_t> v, uint64_t pos);
    uint64_t next_non_zero(std::vector<std::pair<uint64_t, std::deque<TAlphabet> > > v, uint64_t pos);

    /*
     * Helper function that will split up a given range in the graph
     * into bins, one for each character in the alphabet. The split is performed based
     * on the k - d column of the node label. It is assumed that the all nodes in the
     * given range share a common suffix of length d.
     */
    std::vector<uint64_t> split_range(uint64_t start, uint64_t end, uint64_t d /*depth*/);
    void split_range(std::deque<TAlphabet>* str, std::pair<uint64_t, uint64_t>& range);

    /*
     * Helper function to generate the prefix corresponding to 
     * a given bin ID.
     */
    std::deque<TAlphabet> bin_id_to_string(uint64_t bin_id, uint64_t binlen);

    //
    //
    // SERIALIZE
    //
    //
    
    /**
     * This is a debug function that prints the current state of the graph arrays to
     * the screen.
     */
    void print_state();
    
    
    /**
     * This is a debug function that prints the current representation of the graph to
     * the screen.
     */
    void print_state_str();

    /*
     * Returns the sequence stored in W and prints the node
     * information in an overview. 
     * Useful for debugging purposes.
     */
    void print_seq();

    void print_adj_list();

    /**
     * Take the current graph content and store in a file.
     *
     */
    void toFile(unsigned int total = 1, unsigned int idx = 0); 

    /**
     * Visualization, Serialization and Deserialization of annotation content.
     */
    void annotationToFile();
    void annotationFromFile();

};
#endif
