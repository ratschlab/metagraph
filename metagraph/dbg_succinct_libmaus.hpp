#ifndef __DBG_SUCCINCT_LIBM_HPP__
#define __DBG_SUCCINCT_LIBM_HPP__

#include <libmaus2/bitbtree/bitbtree.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>

#include <config.hpp>
#include <datatypes.hpp>

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'X'> > Dna5F; 
    typedef uint64_t TAlphabet;

    public:

    // the bit array indicating the last outgoing edge of a node
    libmaus2::bitbtree::BitBTree<6, 64> *last = new libmaus2::bitbtree::BitBTree<6, 64>();

    // the array containing the edge labels
    libmaus2::wavelet::DynamicWaveletTree<6, 64> *W = new libmaus2::wavelet::DynamicWaveletTree<6, 64>(4); // 4 is log (sigma)

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
    CFG config;

#ifdef DBGDEBUG
    bool debug = true;
#else
    bool debug = false;
#endif 

    /**
     * This object collects information about branches during graph traversal, so 
     * we know where to jump back to when we reached a dead end.
     */
    struct BranchInfo;

    /**
     * This object collects information about branches during graph traversal for the
     * purpose of merging, so we know where to jump back to when we reached a dead end.
     */
    struct BranchInfoMerge;

    /**
     * This will hold the graph edges that will be written to the SQL graph output.
     */
    struct JoinInfo;

    // construct empty graph instance
    DBG_succ(size_t k_,
             CFG config_, 
             bool sentinel = true);

    // load graph instance from a provided file name base
    DBG_succ(std::string infbase_, 
             CFG config_);
    
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
    uint64_t index(seqan::String<Dna5F> &s_);

    uint64_t index(std::deque<TAlphabet> str);

    std::vector<HitInfo> index_fuzzy(seqan::String<seqan::Dna5> &str, uint64_t eops);

    std::pair<uint64_t, uint64_t> index_range(std::deque<TAlphabet> str);

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
    * This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
    * as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
    * returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
    * second value set to true if G2(k2_node) < G1(k1_node).
    */
    std::pair<bool, bool> compare_nodes(DBG_succ *G1, uint64_t k1_node, DBG_succ *G2, uint64_t k2_node); 

    /** 
     * This function gets two node indices and returns if the
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(uint64_t i1, uint64_t i2);

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

    /**
     * Breaks the seq into k-mers and searches for the index of each
     * k-mer in the graph. Returns these indices.
     */
    std::vector<uint64_t> align(seqan::String<seqan::Dna5> seq);

    std::vector<std::vector<HitInfo> > align_fuzzy(seqan::String<seqan::Dna5> seq, uint64_t max_distance = 0, uint64_t alignment_length = 0);



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

    // add a full sequence to the graph
    void add_seq (seqan::String<Dna5F> seq);

    /** This function takes a character c and appends it to the end of the graph sequence
     * given that the corresponding note is not part of the graph yet.
     */
    void append_pos(TAlphabet c);

    //
    //
    // TRAVERSAL
    //
    //
    
    /**
     * This is a convenience function that pops the last branch and updates the traversal state.
     */
    BranchInfo pop_branch(std::stack<BranchInfo> &branchnodes, uint64_t &seqId, uint64_t &seqPos, uint64_t &nodeId, uint64_t &lastEdge, bool &isFirst);

    BranchInfoMerge pop_branch(std::stack<BranchInfoMerge> &branchnodes, uint64_t &nodeId, uint64_t &lastEdge, std::deque<TAlphabet> &last_k);

    bool finish_sequence(seqan::String<Dna5F> &sequence, uint64_t seqId, std::ofstream &SQLstream); 

    size_t traverseGraph(std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream); 

    void allelesFromSeq(seqan::String<Dna5F>seq, unsigned int f, std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream, bool isRefRun = false, size_t seqNum = 0);


    //
    //
    // MERGE
    //
    //

    /**
    * Heavily borrowing from the graph sequence traversal, this function gets a graph pointer G and merges its
    * nodes into the current graph object. The edges of the graph G are fully traversed and nodes are added to
    * the object graph if not existing yet. This function is well suited to merge small graphs into large ones.
    */
    void merge(DBG_succ* G); 

    /*
     * Given two other graph structures G1 and G2, this function 
     * integrate both into a new graph G.
     */
    void merge(DBG_succ* G1, DBG_succ* G2, uint64_t k1 = 1, uint64_t k2 = 1, uint64_t n1 = 0, uint64_t n2 = 0); 

    /**
    * Given a pointer to a graph structure G, the function compares its elements to the
    * current graph. It will perform an element wise comparison of the arrays W, last and
    * F and will only check for identity. If any element differs, the function will return 
    * false and true otherwise.
    */
    bool compare(DBG_succ* G); 

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

    /*
     * Returns the sequence stored in W and prints the node
     * information in an overview. 
     * Useful for debugging purposes.
     */
    void print_seq();

    /**
     * Take the current graph content and return it in SQL
     * format (GA4GH Spec).
     *
     * We will perform one depth first search traversal of the graph. While we will record
     * one long reference string, we will output all sidepaths on the way.
     */
    void toSQL(); 

    /**
     * Take the current graph content and store in a file.
     *
     */
    void toFile(); 

};
#endif
