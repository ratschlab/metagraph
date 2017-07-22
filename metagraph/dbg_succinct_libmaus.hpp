#ifndef __DBG_SUCCINCT_LIBM_HPP__
#define __DBG_SUCCINCT_LIBM_HPP__

#include <zlib.h>
#include <unordered_map>

#include <libmaus2/bitbtree/bitbtree.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>

#include "config.hpp"
#include "datatypes.hpp"

#include "kseq.h"
#include <sdsl/wavelet_trees.hpp>

class rs_bit_vector: public sdsl::bit_vector {
    private:
        sdsl::rank_support_v5<> rk;
        sdsl::select_support_mcl<> slct;
        bool update_rs = true;
        void init_rs() {
            rk = sdsl::rank_support_v5<>(this);
            slct = sdsl::select_support_mcl<>(this);
            update_rs = false;
        }
    public:
        rs_bit_vector(size_t size, bool def) : sdsl::bit_vector(size, def) {
        }
        rs_bit_vector() : sdsl::bit_vector() {
        }
        void set(size_t id, bool val) {
            this->operator[](id) = val;
            update_rs = true;
        }
        void setBitQuick(size_t id, bool val) {
            this->operator[](id) = val;
            update_rs = true;
        }
        void insertBit(size_t id, bool val) {
            this->resize(this->size()+1);
            if (this->size() > 1)
                std::copy_backward(this->begin()+id,(this->end())-1,this->end());
            set(id, val);
            update_rs = true;
        }
        void deleteBit(size_t id) {
            if (this->size() > 1)
                std::copy(this->begin()+id+1,this->end(),this->begin()+id);
            this->resize(this->size()-1);
            update_rs = true;
        }
        void deserialise(std::istream &in) {
            this->load(in);
        }
        void serialise(std::ostream &out) {
            this->serialize(out);
        }
        uint64_t select1(size_t id) {
            if (update_rs)
                init_rs();
            //compensating for libmaus weirdness
            id++;
            size_t maxrank = rk(this->size());
            if (id > maxrank) {
                //TODO: should this line ever be reached?
                return this->size();
            }
            return slct(id);
        }
        uint64_t rank1(size_t id) {
            if (update_rs)
                init_rs();
            //the rank method in SDSL does not include id in the count
            return rk(id >= this->size() ? this->size() : id+1);
        }
};

class dyn_wavelet: public sdsl::int_vector<> {
    private:
        sdsl::wt_int<> wwt;
        size_t logsigma;
        bool update_rs;
        void init_wt() {
            this->resize(n);
            sdsl::construct_im(wwt, *this);
            update_rs=false;
        }
    public:
        size_t n;
        dyn_wavelet(size_t logsigma, size_t size, uint64_t def) 
            : sdsl::int_vector<>(2 * size + 1, def, 1<< logsigma) {
            n = size;
            update_rs = true;
        }
        dyn_wavelet(size_t logsigma)
            : dyn_wavelet(logsigma, 0, 0) {
        }
        void deserialise(std::istream &in) {
            wwt.load(in);
            this->load(in);
            n = this->size();
        }
        dyn_wavelet(std::istream &in) {
            this->deserialise(in);
        }
        void insert(uint64_t val, size_t id) {
            if (n == this->size()) {
                this->resize(2*n+1);
            }
            n++;
            if (this->size() > 1)
                std::copy_backward(this->begin()+id,this->begin()+n-1,this->begin()+n);
            this->operator[](id) = val;
            update_rs = true;
        }
        void remove(size_t id) {
            if (this->size() > 1)
                std::copy(this->begin()+id+1,this->begin()+n,this->begin()+id);
            n--;
            update_rs = true;
        }
        uint64_t rank(uint64_t c, uint64_t i) {
            if (update_rs)
                init_wt();
            return wwt.rank(i >= wwt.size() ? wwt.size() : i+1, c);
        }
        uint64_t select(uint64_t c, uint64_t i) {
            if (update_rs)
                init_wt();
            i++;
            uint64_t maxrank = wwt.rank(wwt.size(), c);
            if (i > maxrank) {
                //TODO: should this line ever be reached?
                return wwt.size();
            }
            return wwt.select(i, c);
        }
        void serialise(std::ostream &out) {
            this->resize(n);
            wwt.serialize(out);
            this->serialize(out);
        }
        sdsl::int_vector<>::reference operator[](const size_t id) {
            update_rs = true;
            return sdsl::int_vector<>::operator[](id);
        }
};

//Child of DynamicWaveletTree allowing for construction from an int vector
class dyn_wavelet2 : public libmaus2::wavelet::DynamicWaveletTree<6, 64> {
    public:
        dyn_wavelet2(size_t b)
            : libmaus2::wavelet::DynamicWaveletTree<6, 64>(b) {
        }
        dyn_wavelet2(std::istream &in)
            : libmaus2::wavelet::DynamicWaveletTree<6, 64>(in) {
        }
        libmaus2::bitbtree::BitBTree<6, 64>* makeTree(std::vector<uint8_t> *W_stat, size_t b) {
            size_t n = W_stat->size();
            std::vector<uint64_t> offsets((1ull << (b-1)) - 1, 0);
            uint64_t v, m, o, p;
            for (size_t i=0;i<n;++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat->at(i);
                o = 0;
                for (size_t ib = 1; ib < b; ++ib) {
                    bool const bit = m & v;
                    if (!bit)
                        offsets.at(o) += 1;
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
            }
            libmaus2::bitbtree::BitBTree<6, 64> *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n*b, false);
            uint64_t co;
            bool bit;
            std::vector<uint64_t> upto_offsets ((1ull << (b - 1)) - 1, 0);
            for (size_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat->at(i);
                o = 0;
                p = i;
                co = 0;
                for (size_t ib = 0; ib < b - 1; ++ib) {
                    bit = m & v;
                    if (bit) {
                        tmp->setBitQuick(ib * n + p + co, true);
                        co += offsets.at(o);
                        p -= upto_offsets.at(o);
                    } else {
                        p -= (p - upto_offsets.at(o)); 
                        upto_offsets.at(o) += 1;
                    }
                    //dtd::cerr << "o: " << o << " offset[o]: " << offsets.at(o) << std::endl;
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
                bit = m & v;
                if (bit) {
                   // std::cerr << "b - 1: " << b - 1 << " n: " << n << " p: " << p << " co: " << co << std::endl;
                    tmp->setBitQuick((b - 1) * n + p + co, true); 
                }
            }
            return tmp;
        }

        dyn_wavelet2(std::vector<uint8_t> *W_stat, size_t b)
            : libmaus2::wavelet::DynamicWaveletTree<6, 64>(makeTree(W_stat, b), b, W_stat->size()) {
        }
};

//libmaus2 structures
typedef libmaus2::bitbtree::BitBTree<6, 64> BitBTree;
//typedef libmaus2::wavelet::DynamicWaveletTree<6, 64> WaveletTree;

//SDSL-based structures
//typedef rs_bit_vector BitBTree;
//typedef dyn_wavelet WaveletTree;
typedef dyn_wavelet2 WaveletTree;

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef uint64_t TAlphabet;

    public:

    // the bit array indicating the last outgoing edge of a node
    BitBTree *last = new BitBTree();
    //libmaus2::bitbtree::BitBTree<6, 64> *last = new libmaus2::bitbtree::BitBTree<6, 64>();

    // the array containing the edge labels
    //libmaus2::wavelet::DynamicWaveletTree<6, 64> *W = new libmaus2::wavelet::DynamicWaveletTree<6, 64>(4); // 4 is log (sigma)
    WaveletTree *W = new WaveletTree(4);

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

    // annotation containers
    std::deque<uint32_t> annotation; // list that associates each node in the graph with an annotation hash
    std::vector<std::string> id_to_label; // maps the label ID back to the original string
    std::unordered_map<std::string, uint32_t> label_to_id_map; // maps each label string to an integer ID
    std::map<uint32_t, uint32_t> annotation_map; // maps the hash of a combination to the position in the combination vector
    std::vector<uint32_t> combination_vector; // contains all known combinations
    uint64_t combination_count = 0;

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
    uint64_t index(std::string &s_);

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
    * This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
    * as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
    * returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
    * second value set to true if G2(k2_node) < G1(k1_node).
    */
    std::pair<bool, bool> compare_nodes(DBG_succ *G1, uint64_t k1_node, DBG_succ *G2, uint64_t k2_node); 

    std::pair<std::vector<bool>, uint64_t> compare_nodes(std::vector<DBG_succ*> G, std::vector<uint64_t> k, std::vector<uint64_t> n, size_t &cnt);

    /** 
     * This function gets two node indices and returns if the
     * node labels share a k-1 suffix.
     */
    bool compare_node_suffix(uint64_t i1, uint64_t i2);

    /**
     *  This function checks whether two given strings given as deques are 
     *  identical.
     */
    bool compare_seq(std::deque<TAlphabet> s1, std::deque<TAlphabet> s2, size_t start = 0);

    /**
     *  This function checks whether string s1 is lexicographically inverse 
     *  greater than s2.
     */
    bool seq_is_greater(std::deque<TAlphabet> s1, std::deque<TAlphabet> s2);

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

    // add a full sequence to the graph
    void add_seq (kstring_t &seq);
    void add_seq_alt (kstring_t &seq);

    /** This function takes a character c and appends it to the end of the graph sequence
     * given that the corresponding note is not part of the graph yet.
     */
    void append_pos(TAlphabet c);

    /** This function takes a pointer to a graph structure and concatenates the arrays W, last 
     * and F to this graph's arrays. In almost all cases this will not produce a valid graph and 
     * should only be used as a helper in the parallel merge procedure.
     */
    void append_graph(DBG_succ *g);

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

    bool finish_sequence(std::string &sequence, uint64_t seqId, std::ofstream &SQLstream); 

    size_t traverseGraph(std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream); 

    void allelesFromSeq(kstring_t &seq, unsigned int f, std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream, bool isRefRun = false, size_t seqNum = 0);

    void traversalHash();

    //
    //
    // ANNOTATE
    //
    //

    void annotate_seq(kstring_t &seq, kstring_t &label, uint64_t start=0, uint64_t end=0, pthread_mutex_t* anno_mutex=NULL);

    void annotate_kmer(std::string &kmer, uint32_t &label, uint64_t &previous, pthread_mutex_t* anno_mutex, bool ignore=false);

    std::vector<uint32_t> classify_path(std::vector<uint64_t> path);

    std::set<uint32_t> classify_read(kstring_t &read, uint64_t max_distance);


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
    
    uint64_t next_non_zero(std::vector<uint64_t> v, uint64_t pos);
    uint64_t next_non_zero(std::vector<std::pair<uint64_t, std::deque<TAlphabet> > > v, uint64_t pos);

    void merge_bins(DBG_succ* G1, DBG_succ* G2, std::deque<TAlphabet>* curr_range, std::pair<uint64_t, uint64_t>& r1, std::pair<uint64_t, uint64_t>& r2);
    
    void merge_fast(DBG_succ* G1, DBG_succ* G2, uint64_t k1 = 1, uint64_t k2 = 1, uint64_t n1 = 0, uint64_t n2 = 0, bool is_parallel = false);
    /*
     * Given two other graph structures G1 and G2, this function 
     * integrate both into a new graph G.
     */
    void merge(DBG_succ* G1, DBG_succ* G2, uint64_t k1 = 1, uint64_t k2 = 1, uint64_t n1 = 0, uint64_t n2 = 0, bool is_parallel = false); 
    void merge2(DBG_succ* G1, DBG_succ* G2, uint64_t k1 = 1, uint64_t k2 = 1, uint64_t n1 = 0, uint64_t n2 = 0, bool is_parallel = false); 
    void merge3(std::vector<DBG_succ*> Gv, std::vector<uint64_t> kv, std::vector<uint64_t> nv, bool is_parallel = false);

    /**
    * Given a pointer to a graph structure G, the function compares its elements to the
    * current graph. It will perform an element wise comparison of the arrays W, last and
    * F and will only check for identity. If any element differs, the function will return 
    * false and true otherwise.
    */
    bool compare(DBG_succ* G); 

    /*
     * Helper function that will split up a given range in the graph
     * into bins, one for each character in the alphabet. The split is performed based
     * on the k - d column of the node label. It is assumed that the all nodes in the
     * given range share a common suffix of length d.
     */
    std::vector<uint64_t> split_range(uint64_t start, uint64_t end, uint64_t d /*depth*/);
    void split_range(std::deque<TAlphabet>* str, std::pair<uint64_t, uint64_t>& range);

    /* 
     * Helper function to determine the bin boundaries, given 
     * a number of threads.
     */
    std::vector<std::pair<uint64_t, uint64_t> > get_bins(uint64_t threads, uint64_t bins_per_thread, DBG_succ* G);
    std::vector<std::pair<uint64_t, uint64_t> > get_bins(uint64_t bins);

    std::vector<std::pair<uint64_t, uint64_t> > get_bins_relative(DBG_succ* G, std::vector<std::pair<uint64_t, uint64_t> > ref_bins, uint64_t first_pos, uint64_t last_pos);

    /*
     * Helper function to generate the prefix corresponding to 
     * a given bin ID.
     */
    std::deque<TAlphabet> bin_id_to_string(uint64_t bin_id, uint64_t binlen);

    /*
     * Distribute the merging of two graph structures G1 and G2 over
     * bins, such that n parallel threads are used. The number of bins
     * is determined dynamically.
     */
    //void merge_parallel(DBG_succ* G1, DBG_succ* G2, uint64_t k1, uint64_t k2, uint64_t n1, uint64_t n2, uint64_t threads);


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
    void toFile(unsigned int total = 1, unsigned int idx = 0); 

    /**
     * Visualization, Serialization and Deserialization of annotation content.
     */
    void annotationToScreen();
    void annotationToFile();
    void annotationFromFile();

};
#endif
