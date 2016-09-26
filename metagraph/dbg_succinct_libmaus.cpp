/**
 * This class contains a succinct representation of the de bruijn graph
 * following ideas and suggestions presented here:
 * http://link.springer.com/chapter/10.1007/978-3-642-33122-0_18
 *
 * There is also conceptual code available at 
 * https://code.google.com/p/csalib/downloads/list
 * that has been used as a reference for this implementation.
 */

#include <vector>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <cmath>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include <libmaus2/bitbtree/bitbtree.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>
#include <libmaus2/digest/md5.hpp>
#include <libmaus2/util/NumberSerialisation.hpp>

#include "config.hpp"
#include "datatypes.hpp"
#include "serialization.hpp"
#include "dbg_succinct_libmaus.hpp"
// use Heng Li's kseq structure for string IO
#include "kseq.h"

// define an extended alphabet for W --> somehow this does not work properly as expected
typedef uint64_t TAlphabet;
typedef DBG_succ::BranchInfo BranchInfo;
typedef DBG_succ::BranchInfoMerge BranchInfoMerge;
typedef DBG_succ::JoinInfo JoinInfo;

/**
 * This object collects information about branches during graph traversal, so 
 * we know where to jump back to when we reached a dead end.
 */
struct DBG_succ::BranchInfo {
    uint64_t nodeId;
    uint64_t seqId;
    uint64_t seqPos;
    TAlphabet lastEdge;

    BranchInfo(uint64_t nodeId_ = 0, uint64_t seqId_ = 0, uint64_t seqPos_ = 0, TAlphabet lastEdge_ = 0):
        nodeId(nodeId_),
        seqId(seqId_),
        seqPos(seqPos_),
        lastEdge(lastEdge_) {}
};

/**
 * This object collects information about branches during graph traversal for the
 * purpose of merging, so we know where to jump back to when we reached a dead end.
 */
struct DBG_succ::BranchInfoMerge {
    uint64_t nodeId;
    TAlphabet lastEdge;
    std::deque<TAlphabet> last_k;

    BranchInfoMerge() {}

    BranchInfoMerge(uint64_t nodeId_, TAlphabet lastEdge_, std::deque<TAlphabet> last_k_):
        nodeId(nodeId_),
        lastEdge(lastEdge_),
        last_k(last_k_) {}
};


/**
 * This will hold the graph edges that will be written to the SQL graph output.
 */
struct DBG_succ::JoinInfo {
    uint64_t seqId1;
    uint64_t seqPos1;
    uint64_t seqId2;
    uint64_t seqPos2;

    JoinInfo(uint64_t seqId1_ = 0, uint64_t seqPos1_ = 0, uint64_t seqId2_ = 0, uint64_t seqPos2_ = 0):
        seqId1(seqId1_),
        seqPos1(seqPos1_),
        seqId2(seqId2_),
        seqPos2(seqPos2_) {}
};

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
// alphabet
const char alphabet[] = "$ACGTNX$ACGTNXn";

// infile base when loaded from file
std::string infbase;

// config object
//CFG config;


#ifdef DBGDEBUG
    bool debug = true;
#else
    bool debug = false;
#endif 

//
//
// CONSTRUCTORS
//
//
DBG_succ::DBG_succ(size_t k_, Config* config_, bool sentinel) : 
    k(k_),
    config(config_) {

    last->insertBit(0, false);
    if (sentinel)
        last->insertBit(1, true);

    W->insert(0, 0);
    if (sentinel)
        W->insert(0, 0);

    F.push_back(0);
    if (sentinel) {
        for (size_t j = 1; j < alph_size; j++)
            F.push_back(1);
        p = 1;
    } else {
        for (size_t j = 1; j < alph_size; j++)
            F.push_back(0);
        p = 0;
    }
}

DBG_succ::DBG_succ(std::string infbase_, Config* config_) : 
    infbase(infbase_),
    config(config_) {

    // load last array
    std::ifstream instream((infbase + ".l.dbg").c_str());
    last->deserialise(instream);
    instream.close();

    // load W array
    //delete W;
    instream.open((infbase + ".W.dbg").c_str());
    W = new libmaus2::wavelet::DynamicWaveletTree<6, 64>(instream);
    instream.close();

    // load F and k and p
    instream.open((infbase + ".F.dbg").c_str());
    std::string line;
    size_t mode = 0;
    while (std::getline(instream, line)) {
        if (strcmp(line.c_str(), ">F") == 0) {
            mode = 1;
        } else if (strcmp(line.c_str(), ">k") == 0) {
            mode = 2;
        } else if (strcmp(line.c_str(), ">p") == 0) {
            mode = 3;
        } else {
            if (mode == 1) {
                F.push_back(std::strtoul(line.c_str(), NULL, 10));
            } else if (mode == 2) {
                k = strtoul(line.c_str(), NULL, 10);
            } else if (mode == 3) {
                p = strtoul(line.c_str(), NULL, 10);
            } else {
                fprintf(stderr, "ERROR: input file corrupted\n");
                exit(1);
            }
        }
    }
    instream.close();
}

DBG_succ::~DBG_succ() {
    delete W;
    delete last;
}

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
uint64_t DBG_succ::rank_W(uint64_t i, TAlphabet c) {

    // deal with  border conditions
    if (i <= 0)
        return 0;
    return W->rank(c, std::min(i, W->n - 1)) - (c == 0);
}

/**
 * Uses the array W and gets a count i and a character c from 
 * the alphabet and returns the positions of the i-th occurence of c 
 * in W.
 */
uint64_t DBG_succ::select_W(uint64_t i, TAlphabet c) {
    
    // deal with  border conditions
    if (i <= 0)
        return 0;

    // count occurences of c and store them in cnt
    //fprintf(stderr, "query select W -- c: %i i: %lu return: %lu \n", c, i-1+(c==0), W->select(c, i-1+(c==0)));
    return std::min(W->select(c, i - 1 + (c == 0)), W->n);
}

/**
 * This is a convenience function that returns for array W, a position i and 
 * a character c the last index of a character c preceding in W[1..i].
 */
uint64_t DBG_succ::pred_W(uint64_t i, TAlphabet c) {
    return select_W(rank_W(i, c), c);
}

/**
 * This is a convenience function that returns for array W, a position i and 
 * a character c the first index of a character c in W[i..N].
 */
uint64_t DBG_succ::succ_W(uint64_t i, TAlphabet c) {
    return select_W(rank_W(i - 1, c) + 1, c);
}

/** 
 * Uses the object's array last and a position and
 * returns the number of set bits up to that postion.
 */
uint64_t DBG_succ::rank_last(uint64_t i) {
    // deal with  border conditions
    if (i <= 0)
        return 0;
    return last->rank1(i);
}

/**
 * Uses the object's array last and a given position i and
 * returns the position of the i-th set bit in last[1..i].
 */
uint64_t DBG_succ::select_last(uint64_t i) {
    // deal with  border conditions
    if (i <= 0)
        return 0;
    // for some reason the libmaus2 select is 0 based ...
    return std::min(last->select1(i - 1), last->size());
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the last set bit in last[1..i].
 */
uint64_t DBG_succ::pred_last(uint64_t i) {
    return select_last(rank_last(i));
}

/**
 * This is a convenience function that returns for the object's array last
 * and a given position i the position of the first set bit in last[i..N].
 */
uint64_t DBG_succ::succ_last(uint64_t i) {
    return select_last(rank_last(i - 1) + 1);
}

/**
 * This function gets a position i that reflects the i-th node and returns the
 * position in W that corresponds to the i-th node's last character. 
 */
uint64_t DBG_succ::bwd(uint64_t i) {
    // get value of last position in node i
    TAlphabet c = get_node_end_value(i);
    // get the offset for the last position in node i
    uint64_t o = F[c];
    //fprintf(stdout, "i %lu c %i o %lu rank(i) %lu rank(o) %lu difference %lu\n", i, (int) c, o, rank_last(i), rank_last(o), rank_last(i) - rank_last(o));
    // compute the offset for this position in W and select it
    return select_W(rank_last(i) - rank_last(o), c);
}

/**
 * This functions gets a position i reflecting the r-th occurence of the corresponding
 * character c in W and returns the position of the r-th occurence of c in last.
 */
uint64_t DBG_succ::fwd(uint64_t i) {
    // get value of W at position i
    TAlphabet c = (*W)[i] % alph_size; 
    // get the offset for position c
    uint64_t o = F[c];
    // get the rank of c in W at position i
    uint64_t r = rank_W(i, c);
    // select the index of the position in last that is rank many positions after offset
    return select_last(rank_last(o) + r);
}


/**
 * Using the offset structure F this function returns the value of the last 
 * position of node i.
 */
TAlphabet DBG_succ::get_node_end_value(uint64_t i) {
    if (i == 0)
        return 0;
    for (size_t j = 0; j < F.size(); j++) {
        if (F[j] >= i)
            return j - 1;
    }
    return F.size() - 1;
}


/**
 * Given index of node i, the function returns the 
 * first character of the node.
 */
TAlphabet DBG_succ::get_node_begin_value(uint64_t i) {
    if (i == 1)
        return 0;

    for (size_t j = 0; j < k-1; ++j) {
        i = bwd(succ_last(i));
        if (i == 1)
            return 0;
    }
    return get_node_end_value(i);
}


/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the node the edge is pointing to.
 */
uint64_t DBG_succ::outgoing(uint64_t i, TAlphabet c) {
    if (i > W->n)
        return 0;
    if (i == 0)
        return 0;
    std::pair<uint64_t, uint64_t> R = get_equal_node_range(i);

    uint64_t j1 = pred_W(R.second, c);
    uint64_t j2 = pred_W(R.second, c + alph_size);
    uint64_t j = (j1 < j2) ? j2 : j1;
    if (j < R.first || j >= W->n)
        return 0;
    j = fwd(j);
    if (j == 0 || j == W->n)
        return 0;
    return j;
}

/**
 * Given a node index i and an edge label c, this function returns the
 * index of the node the incoming edge belongs to.
 */
uint64_t DBG_succ::incoming(uint64_t i, TAlphabet c) {
    if (i == 1)
        return 0;
    c %= alph_size;
    TAlphabet d = get_node_end_value(i);
    uint64_t x = bwd(i);
    if (get_node_begin_value(x) == c) {
        return succ_last(x);
    }
    uint64_t y = succ_W(x + 1, d);
    while (true) {
        x = succ_W(x+1, d + alph_size);
        if (x >= y) {
            break;
        }
        if (get_node_begin_value(x) == c) {
            return succ_last(x);
        }
    }
    return 0;
}

/**
 * Given a node index i, this function returns the number of outgoing
 * edges from node i.
 */
uint64_t DBG_succ::outdegree(uint64_t i) {
    return (i < W->n) ? succ_last(i) - pred_last(i - 1) : 0;
}


/**
 * Given a node index i, this function returns the number of incoming
 * edges to node i.
 */
uint64_t DBG_succ::indegree(uint64_t i) {
    if (i < 2)
        return 0;
    uint64_t x = bwd(succ_last(i));
    TAlphabet d = get_node_end_value(i);
    uint64_t y = succ_W(x + 1, d);
    return 1 + rank_W(y, d + alph_size) - rank_W(x, d + alph_size);
}

/**
 * Given a node label s, this function returns the index
 * of the corresponding node, if this node exists and 0 otherwise.
 */
uint64_t DBG_succ::index(std::string &s_) {
    TAlphabet s = get_alphabet_number(s_[0]);
    // init range
    uint64_t rl = succ_last(F[s] + 1);
    uint64_t ru = F[s + 1]; // upper bound
    // update range iteratively while scanning through s
    for (uint64_t i = 1; i < s_.length(); i++) {
        s = get_alphabet_number(s_[i]);
        rl = std::min(succ_W(pred_last(rl - 1) + 1, s), succ_W(pred_last(rl - 1) + 1, s + alph_size));
        if (rl >= W->n)
            return 0;
        ru = std::max(pred_W(ru, s), pred_W(ru, s + alph_size));
        if (ru >= W->n)
            return 0;
        if (rl > ru)
            return 0;
        rl = outgoing(rl, s);
        ru = outgoing(ru, s);
    }
    return (ru > rl) ? ru : rl;
}

/** 
 * Given a string str and a maximal number of edit operations
 * max_distance, this function returns all nodes with labels at most
 * max_distance many edits away from str.
 */
std::vector<HitInfo> DBG_succ::index_fuzzy(std::string &str, uint64_t max_distance) {
    
    std::vector<HitInfo> result; 
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits;
    std::priority_queue<HitInfo, std::vector<HitInfo>, HitInfoCompare> hits2;
    uint64_t rl;
    uint64_t ru;

    // walk through pattern, thereby collecting possible partial matches
    // once the end of the pattern is readched, add match to results
    
    // init match/mismatch to first pattern position
    TAlphabet s = get_alphabet_number(str[0]);
    for (TAlphabet b = 1; b < 5; ++b) {
        rl = succ_last(F[b] + 1);
        ru = F[b + 1];
        //std::cout << "pushing: rl " << rl << " ru " << ru << " str_pos 1 max_distance " << (uint64_t) (b != s) << std::endl;
        //std::cout << "s " << s << " b " << b << std::endl;
        std::vector<uint64_t> tmp;
        hits.push(HitInfo(rl, ru, 1, 1, (uint64_t) (b != s), std::string(1, get_alphabet_symbol(b)), tmp));

        // opening/extending a gap in the pattern starting with the first position
        for (size_t p = 1; p < str.length() - 1; ++p) {
            TAlphabet ss = get_alphabet_number(str[p]);
            if ((p + (b != ss)) > max_distance)
                break;
            hits.push(HitInfo(rl, ru, p + 1, 1, p + (b != ss), std::string(p, 'd') + std::string(1, get_alphabet_symbol(b)), tmp));
            //std::cout << "a) adding '-'" << std::endl;
        }
    }

    // walk through pattern thereby extending all partial hits
    while (hits.size() > 0) {
        while (hits.size() > 0) {
            HitInfo curr_hit = hits.top();
            hits.pop();
            //std::cout << "loaded: rl " << curr_hit.rl << " ru " << curr_hit.ru << " dist " << curr_hit.distance << std::endl;

            if (curr_hit.str_pos < str.length()) {
            
                // opening/extending a gap in the graph, leaving current pattern position unmatched
                if (curr_hit.distance < max_distance) {
                    hits2.push(HitInfo(curr_hit.rl, curr_hit.ru, curr_hit.str_pos + 1, curr_hit.graph_pos, curr_hit.distance + 1, curr_hit.cigar + 'd', curr_hit.path));
                    //std::cout << "b) " << curr_hit.cigar << " adding '-'" << std::endl;
                }

                s = get_alphabet_number(str[curr_hit.str_pos]);

                // has the number of matches exceeded the node length?
                // there are three possible scenarios for extension of the path:
                //  1) pattern is shorter than the node length --> get an interval of matching nodes
                //  2) pattern length exactly mathces the node length --> there is one correponding node
                //  3) pattern is longer than the node length --> we append to a path 
                if (curr_hit.graph_pos >= k) {
                //    std::cout << "push back tp path " << curr_hit.rl << std::endl;
                    curr_hit.path.push_back(curr_hit.rl);
                }

                // iterate through all possible extensions of current position
                for (TAlphabet b = 1; b < 5; ++b) {
                    if (curr_hit.distance <= max_distance) {

                        // we cannot afford any more mismatches
                        if ((curr_hit.distance + (b != s)) > max_distance)
                            continue;

                        // re-define range of nodes to check for outgoing nodes
                        rl = std::min(succ_W(pred_last(curr_hit.rl - 1) + 1, b), succ_W(pred_last(curr_hit.rl - 1) + 1, b + alph_size));
                        ru = std::max(pred_W(curr_hit.ru, b), pred_W(curr_hit.ru, b + alph_size));

                        // the current range in W does not contain our next symbol
                        if ((rl >= W->n) || (ru >= W->n) || (rl > ru))
                            continue;

                        // update the SA range with the current symbol b
                        rl = outgoing(rl, b);
                        ru = outgoing(ru, b);

                        // range is empty
                        if ((rl == 0) && (ru == 0))
                            continue;

                        // add hit for extension in next step
                        hits2.push(HitInfo(rl, ru, curr_hit.str_pos + 1, curr_hit.graph_pos + 1, curr_hit.distance + (b != s), curr_hit.cigar + get_alphabet_symbol(b), curr_hit.path));
                        //std::cout << "curr rl " << curr_hit.rl << " ru " << curr_hit.ru << std::endl;
                        //std::cout << "c) " << curr_hit.cigar << " adding '" << get_alphabet_symbol(b) << "' - graph pos " << curr_hit.graph_pos + 1 << " rl " << rl << " ru " << ru << std::endl;
                        
                        // opening/extending a gap in the pattern, leaving current graph position unmatched 
                        // --> choose any available mismatching next edge 
                        if (b != s) {
                            hits2.push(HitInfo(rl, ru, curr_hit.str_pos, curr_hit.graph_pos + 1, curr_hit.distance + 1, curr_hit.cigar + 'i', curr_hit.path));
                        }
                    }
                }
            } else {
                // collect results
                //std::cout << "pushing " << curr_hit.cigar << "  " << curr_hit.distance << " rl " << curr_hit.rl << " ru " << curr_hit.ru << std::endl;
                result.push_back(curr_hit); //std::make_pair(curr_hit.rl < curr_hit.ru ? curr_hit.ru : curr_hit.rl, curr_hit.cigar)); 
            }
        }
        hits.swap(hits2);
    }

    return result;
}


uint64_t DBG_succ::index(std::deque<TAlphabet> str) {
    std::pair<uint64_t, uint64_t> tmp = index_range(str);
    return (tmp.second > tmp.first) ? tmp.second : tmp.first;
}

std::pair<uint64_t, uint64_t> DBG_succ::index_range(std::deque<TAlphabet> str) {
    // get first 
    TAlphabet s = str.front() % alph_size;
    str.pop_front();
    // init range
    uint64_t rl = succ_last(F[s] + 1);
    uint64_t ru = F[s + 1];                // upper bound
    //fprintf(stderr, "char: %i rl: %i ru: %i\n", (int) s, (int) rl, (int) ru);
    // update range iteratively while scanning through s
    while (!str.empty()) {
        s = str.front() % alph_size;
        str.pop_front();
        rl = std::min(succ_W(pred_last(rl - 1) + 1, s), succ_W(pred_last(rl - 1) + 1, s + alph_size));
        if (rl >= W->n)
            return std::make_pair(0, 0);
        ru = std::max(pred_W(ru, s), pred_W(ru, s + alph_size));
        if (ru >= W->n)
            return std::make_pair(0, 0);
        if (rl > ru)
            return std::make_pair(0, 0);
        rl = outgoing(rl, s);
        ru = outgoing(ru, s);
    }
    return std::make_pair(rl, ru);
}

/**
 * Given a position i, this function returns the boundaries of the interval
 * of nodes identical to node i (ignoring the values in W).
 */
std::pair<uint64_t, uint64_t> DBG_succ::get_equal_node_range(uint64_t i) {
    return std::make_pair(pred_last(i - 1) + 1, succ_last(i));
}

/**
 * Given index i of a node and a value k, this function 
 * will return the k-th last character of node i.
 */
std::pair<TAlphabet, uint64_t> DBG_succ::get_minus_k_value(uint64_t i, uint64_t k) {
    for (; k > 0; --k)
        i = bwd(succ_last(i));
    return std::make_pair(get_node_end_value(i), bwd(succ_last(i)));
}

/**
* This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
* as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
* returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
* second value set to true if G2(k2_node) < G1(k1_node).
*/
std::pair<bool, bool> DBG_succ::compare_nodes(DBG_succ *G1, uint64_t k1_node, DBG_succ *G2, uint64_t k2_node) {

    std::pair<TAlphabet, uint64_t> k1_val;
    std::pair<TAlphabet, uint64_t> k2_val;
    uint64_t curr_k = 0;
    while (curr_k < this->k) {
        k1_val = G1->get_minus_k_value(k1_node, 0);
        k2_val = G2->get_minus_k_value(k2_node, 0);
        if (k1_val.first == k2_val.first) {
            ++curr_k;
            k1_node = k1_val.second;
            k2_node = k2_val.second;
        } else {
            break;
        }
    }
    //std::cerr << "k1_val: " << k1_val.first << " k2_val: " << k2_val.first << " curr_k:" << curr_k << std::endl;
    return std::make_pair(k1_val.first < k2_val.first, k2_val.first < k1_val.first);
}

/** 
 * This function gets two node indices and returns whether the
 * node labels share a k-1 suffix.
 */
bool DBG_succ::compare_node_suffix(uint64_t i1, uint64_t i2) {
    for (size_t ii = 0; ii < k-1; ii++) {
        if (get_node_end_value(i1) != get_node_end_value(i2)) {
            return false;
        }
        i1 = bwd(succ_last(i1));
        i2 = bwd(succ_last(i2));
    }
    return true;
}

/**
 * This function returns true if node i is a terminal node.
 */
bool DBG_succ::is_terminal_node(uint64_t i) {
    for (size_t ii = 0; ii < k-1; ii++) {
        if (get_node_end_value(i) % alph_size != 0) {
            return false;
        }
        i = bwd(i);
    }
    return true;
}

/**
* Given a node index k_node, this function returns the k-mer sequence of the
* node in a deque data structure.
*/
std::deque<TAlphabet> DBG_succ::get_node_seq(uint64_t k_node) {
    std::deque<TAlphabet> ret;
    std::pair<TAlphabet, uint64_t> k_val;
    for (uint64_t curr_k = 0; curr_k < this->k; ++curr_k) {
        k_val = get_minus_k_value(k_node, 0);
        ret.push_front(k_val.first);
        k_node = k_val.second;
    }
    
    return ret;
}

/**
* Given a node index k_node, this function returns the k-mer sequence of the 
* node as a string.
*/
std::string DBG_succ::get_node_str(uint64_t k_node) {
    std::stringstream ret;
    std::pair<TAlphabet, uint64_t> k_val;
    for (uint64_t curr_k = 0; curr_k < this->k; ++curr_k) {
        k_val = get_minus_k_value(k_node, 0);
        ret << k_val.first;
        k_node = k_val.second;
    }
    return ret.str();
}

/**
 * Return number of edges in the current graph.
 */
uint64_t DBG_succ::get_size() {
    return W->n;
}

/**
 * Return k-mer length of current graph.
 */
uint64_t DBG_succ::get_k() {
    return this->k;
}

/**
 * Return value of W at position k.
 */
TAlphabet DBG_succ::get_W(uint64_t k) {
    return (*W)[k];
}

/** 
 * Return value of F vector at index k.
 * The index is over the alphabet!
 */
TAlphabet DBG_succ::get_F(uint64_t k) {
    return F.at(k);
}

/**
 * Return value of last at position k.
 */
bool DBG_succ::get_last(uint64_t k) {
    return (*last)[k];
}

char DBG_succ::get_alphabet_symbol(uint64_t s) {
    return s < strlen(alphabet) ? alphabet[s] : alphabet[14];
}

std::vector<uint64_t> DBG_succ::align(kstring_t seq) {

  uint64_t kmer_length = this->get_k();
  uint64_t no_kmers_in_seq = seq.l - kmer_length + 1;

  std::vector<uint64_t> indices (no_kmers_in_seq, 0);

  for (uint64_t i = 0; i < no_kmers_in_seq; ++i) {
  // TODO make sure that kmer is shorter than seq.l
    std::string kmer(seq.s + i, seq.s + i + kmer_length);
    indices[i] = this->index(kmer);
  }
  
  return indices;
}

std::vector<std::vector<HitInfo> > DBG_succ::align_fuzzy(kstring_t seq, uint64_t alignment_length, uint64_t max_distance) {

    size_t seq_length = seq.l;
    std::vector<std::vector<HitInfo> > hit_list;

    if (alignment_length == 0) {

    } else {
        alignment_length = alignment_length < 2 ? 2 : alignment_length;
        alignment_length = alignment_length < seq_length ? alignment_length : seq_length;

        for (uint64_t i = 0; i < seq_length - alignment_length + 1; ++i) {
            std::string kmer(seq.s + i, seq.s + i + alignment_length);
            hit_list.push_back(this->index_fuzzy(kmer, max_distance));

        }
    }
    return hit_list;
}

/*
 * Returns the number of nodes on the current graph.
 */
uint64_t DBG_succ::get_node_count() {
    return rank_last(last->size() - 1);
}

/*
 * Return the number of edges in the current graph.
 */
uint64_t DBG_succ::get_edge_count() {
    return W->n - 1;
}

//
//
// APPEND
// 
//

/**
 * This function gets a value of the alphabet c and updates the offset of 
 * all following values by +1 is positive is true and by -1 otherwise.
 */
void DBG_succ::update_F(TAlphabet c, bool positive) {
    for (TAlphabet i = c+1; i < F.size(); i++)
        F[i] += (2 * (TAlphabet) positive - 1);
}

/**
 * This function gets a local range in W from lower bound l
 * to upper bound u and swaps the inserted element to the
 * righ location.
 */
void DBG_succ::sort_W_locally(uint64_t l, uint64_t u) {
    for (uint64_t s = u; s > l; s--) {
        TAlphabet tmp;
        if (((*W)[s] % alph_size) < ((*W)[s-1] % alph_size)) {
            tmp = (*W)[s-1];
            replaceW(s-1, (*W)[s]);
            replaceW(s, tmp);
        }
    }
    for (uint64_t s = l; s < u; ++s) {
        TAlphabet tmp;
        if (((*W)[s] % alph_size) > ((*W)[s+1] % alph_size)) {
            tmp = (*W)[s+1];
            replaceW(s+1, (*W)[s]);
            replaceW(s, tmp);
        }
    }
}


/** 
 * This is a convenience function to replace the value at
 * position i in W with val.
 */
void DBG_succ::replaceW(size_t i, TAlphabet val) {
    W->remove(i);
    W->insert(val, i);
}


TAlphabet DBG_succ::get_alphabet_number(char s) {

     TAlphabet nt_lookup[128] = {
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
        5, 5, 5, 5,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5, 
        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5, 
        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5, 
        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5 
    };

    return nt_lookup[(int) s];
}


// add a full sequence to the graph
void DBG_succ::add_seq (kstring_t &seq) {

    if (debug) {
        print_seq();
        print_state();
        std::cout << "======================================" << std::endl;
    }

    // Padding of the input genome / read
    if (W->n == 2) {
        for (size_t j = 0; j < k; j++) {
            append_pos(6);
            if (debug) {
                print_seq();
                print_state();
                std::cout << "======================================" << std::endl;
            }
        }
    }

    /** Iterate over input sequence and enumerae all
     * k-mers.
     */
    for (uint64_t i = 0; i < seq.l; ++i) {
        if (i > 0 && i % 1000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10000 == 0) {
                fprintf(stdout, "%lu - edges %lu / nodes %lu\n", i, get_edge_count(), get_node_count());
            }
        }
        // if (debug) {
        //    fprintf(stdout, "appending %i\n", (int) ordValue(seq[i]));
        //    cerr << "seq[i] " << seq[i] << std::endl;
        //    cerr << "seq[i] ord " << ordValue(seq[i]) + 1 << std::endl;
        // }

        append_pos(get_alphabet_number(seq.s[i]));

        if (debug) {
            print_seq();
            print_state();
            std::cout << "======================================" << std::endl;
        }
    }

    // Padding after sequence to get back into default state.
    for (size_t j = 0; j < k; j++) {
        append_pos(6);
        if (debug) {
            print_seq();
            print_state();
            std::cout << "======================================" << std::endl;
        }
    }

    fprintf(stdout, "edges %lu / nodes %lu\n", get_edge_count(), get_node_count());

    //toSQL();
    //String<Dna5F> test = "CCT";
    //fprintf(stdout, "\nindex of CCT: %i\n", (int) index(test));
}

/** This function takes a character c and appends it to the end of the graph sequence
 * given that the corresponding note is not part of the graph yet.
 */
void DBG_succ::append_pos(TAlphabet c) {

    // check that the last position of the graph is indeed a terminal
    assert((*W)[p] == 0);
    TAlphabet c_p = get_node_end_value(p);
    // get range of identical nodes (without W) pos current end position
    std::pair<uint64_t, uint64_t> R = get_equal_node_range(this->p);
    //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

    // get position of first occurence of c in W after p
    uint64_t next_c = succ_W(p, c);
    // check if c is part of range
    bool exist_c = (next_c <= R.second);
    if (!exist_c) {
        // get position of first occurence of c- in W after p
        next_c = succ_W(p, c + alph_size);
        // check if c- is part of range
        exist_c = (next_c <= R.second);
    }

    /**
     * if the character already exists in the range, we delete the terminal symbol
     * at p, insert c at fwd(next_c) and update p.
     */
    if (exist_c) {
        uint64_t p_new = fwd(next_c);
        // remove old terminal symbol
        last->deleteBit(p);
        W->remove(p);
        // adapt position if altered by previous deletion
        p_new -= (p < p_new);
        // insert new terminal symbol 
        // we have to insert 0 into last as the node already existed in the range 
        // and the terminal symbol is always first
        last->insertBit(p_new, false);
        W->insert(0, p_new);
        // update new terminal position
        p = p_new;
        // take care of updating the offset array F
        update_F(c_p, false);
        //assert(get_node_end_value(p) == c);
        update_F(c, true);
    } else {
        /**
         * We found that c does not yet exist in the current range and now have to
         * figure out if we need to add c or c- to the range.
         * To do this, we check if there is a previous position j1 with W[j1] == c
         * whose node shares a k-1 suffix with the current node. If yes, we add c- 
         * instead of c.
         */
        // get position of last occurence of c before p (including p - 1)
        uint64_t last_c = pred_W(p - 1, c);
        // if this position exists
        if (last_c > 0) {
            uint64_t x = fwd(last_c);
            assert((*last)[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

            // check, if there are any c or c- symbols following after position p
            uint64_t next_c = succ_W(p + 1, c);
            uint64_t next_cm = succ_W(p + 1, c + alph_size);
            // there is no c between p and next_cm and next_cm is a c- ==> we should add a c- 
            // all nodes between W[i] = c and W[j] = c- share a common suffix of length k-1
            bool minus1 = (next_cm < next_c);
            // check, if we share a k-1 suffix with last_c
            if (!minus1) {
                minus1 = compare_node_suffix(p, last_c);
            }

            // adding a new node can influence following nodes that share a k-1 suffix with the
            // new node -> need to adapt the respektive cc to a cc-
            bool minus2 = false;
            if (next_c < W->n) {
                minus2 = compare_node_suffix(p, next_c);
                if (minus2) {
                    replaceW(next_c, (*W)[next_c] + alph_size);
                }
            }

            replaceW(p, minus1 ? c + alph_size : c);
            // after we are done, assert that the order within the range we created 
            // is still valid within W
            if (p - R.second > 0) {
                sort_W_locally(p, R.second);
            }

            // if minus1 is true, we share a k-1 suffix with the node at 
            // last_c and thus need to adapt our insertion position by -1, 
            // as we would like to insert before it. Otherwise we insert directly after
            // it as we are now sorted after it. 
            if (minus1) {
                p = x;
                last->insertBit(x, false);
                W->insert(0, x);
            } else if (minus2) {
                p = x + 1;
                last->insertBit(x + 1, false);
                W->insert(0, x + 1);
            // no node shares a k-1 suffix with last_c and thus the new node comes after
            // the forward of last_c (as the current node came after last_c as well)
            } else {
                p = x + 1;
                last->insertBit(x + 1, true);
                W->insert(0, x + 1);
            }
        } else {
            uint64_t x = F[c] + 1;
            uint64_t next_c = succ_W(p + 1, c);
            bool minus = false;
            if (next_c < W->n) {
                minus = compare_node_suffix(p, next_c);
            }
            replaceW(p, c);
            if (p - R.second > 0) {
                sort_W_locally(p, R.second);
            }
            p = x;
            if (minus) {
                replaceW(next_c, (*W)[next_c] + alph_size);
                last->insertBit(x, false);
            } else {
                last->insertBit(x, true);
            }
            W->insert(0, x);
        }
        update_F(c, true);
    }
    // update sorting at new location of p
    // with this we assert that $ is always inserted at the first position 
    // of a range of equal nodes --> this will help us to prevent multiple insertions
    // of already existing nodes
    R = get_equal_node_range(this->p);
    if (R.second - R.first > 0) {
        sort_W_locally(R.first, R.second);
        while ((*W)[p] != 0)
            p--;
        assert((*W)[p] == 0);
    }
}


/** This function takes a pointer to a graph structure and concatenates the arrays W, last 
 * and F to this graph's arrays. In almost all cases this will not produce a valid graph and 
 * should only be used as a helper in the parallel merge procedure.
 */
void DBG_succ::append_graph(DBG_succ *g) {

    size_t curr_pos = this->get_size();

    // handle last and W
    for (size_t j = 1; j < g->get_size(); ++j) {
        this->last->insertBit(curr_pos, g->get_last(j));
        this->W->insert(g->get_W(j), curr_pos);
        ++curr_pos;
    }

    // handle F
    assert(this->F.size() == g->F.size());
    for (size_t j = 0; j < this->F.size(); ++j) {
        this->F.at(j) += g->F.at(j);
    }
}

//
//
// TRAVERSAL
//
//
    

/**
 * This is a convenience function that pops the last branch and updates the traversal state.
 */
BranchInfo DBG_succ::pop_branch(std::stack<BranchInfo> &branchnodes, uint64_t &seqId, uint64_t &seqPos, uint64_t &nodeId, uint64_t &lastEdge, bool &isFirst) {
    BranchInfo branch = branchnodes.top();
    branchnodes.pop();
    isFirst = true;
    seqPos = branch.seqPos;
    seqId = branch.seqId;
    lastEdge = branch.lastEdge;
    nodeId = branch.nodeId;

    return branch;
}

BranchInfoMerge DBG_succ::pop_branch(std::stack<BranchInfoMerge> &branchnodes, uint64_t &nodeId, uint64_t &lastEdge, std::deque<TAlphabet> &last_k) {
    BranchInfoMerge branch = branchnodes.top();
    branchnodes.pop();
    lastEdge = branch.lastEdge;
    nodeId = branch.nodeId;
    last_k = branch.last_k;

    return branch;
}


bool DBG_succ::finish_sequence(std::string &sequence, uint64_t seqId, std::ofstream &SQLstream) {
    if (sequence.length() > 0) {
        if (seqId == 1)
            SQLstream << "INSERT INTO FASTA VALUES (1, '" << config->sqlfbase << ".fa');" << std::endl;
        std::ofstream stream;
        if (seqId == 1)
            stream.open((config->sqlfbase + ".fa").c_str());
        else
            stream.open((config->sqlfbase + ".fa").c_str(), std::ofstream::app);
        stream << ">seq" << seqId << std::endl;
        uint64_t i = 0;
        while ((i + 80) < sequence.length()) {
            stream << sequence.substr(i, 80) << std::endl;
            i += 80;
        }
        if (i != sequence.length())
            stream << sequence.substr(i) << std::endl;
        stream.close();

        if (debug)
            std::cout << sequence << std::endl;

        std::string md5;
        std::ostringstream test;
        test << sequence;
        libmaus2::util::MD5::md5(test.str(), md5);
        SQLstream << "INSERT INTO Sequence VALUES (" << seqId << ", 1, 'seq" << seqId << "', '" << md5 << "', " << sequence.length() << ");" << std::endl;
        sequence.clear();
        return true;
    } else {
        return false;
    }
}

size_t DBG_succ::traverseGraph(std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream) {
    // store all branch nodes on the way
    std::stack<BranchInfo> branchnodes;
    std::map<uint64_t, std::pair<uint64_t, uint64_t> > nodeId2seqId;
    // bool vector that keeps track of visited nodes
    std::vector<bool> visited(last->size());
    for (std::vector<bool>::iterator it = visited.begin(); it != visited.end(); ++it) {
        *it = false;
    }
    std::string sequence;
    // for nodes with indegree > 1 we store sequence and index of the 
    // sequence that visited them, so we know where to anchor branches into it
    std::map<uint64_t, std::pair<uint64_t, uint64_t> > nodeId2seqPos; 
    // some initializations
    uint64_t nodeId = 1; // start at source node
    uint64_t seqPos = 0; // position in currently traversed sequence, will increase with every visited node and be reset upon new branch
    uint64_t seqId = 1;  // first sequence ID is 1
    size_t seqCnt = 1; // number of total sequences
    bool isFirst = true; // is this the first node in a new sequence?
    size_t out = outdegree(nodeId);
    std::pair<uint64_t, uint64_t> branchPos;
    BranchInfo branch;
    TAlphabet val;
    TAlphabet lastEdge = 0;
    TAlphabet joinEdge = 0;
    bool joinOpen = false;
    JoinInfo currJoin;

    while (out > 0 || branchnodes.size() > 0) {

        //fprintf(stderr, "1: nodeId %lu out %lu\n", nodeId, out);
        // we have reached the sink but there are unvisited nodes left on the stack
        if (out == 0) {
            if (branchnodes.size() == 0)
                break;
            // get new branch
            branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
            out = outdegree(nodeId);
            if (debug)
                fprintf(stderr, " -- popped %lu -- ", nodeId); 
            joinOpen = true;
            joinEdge = lastEdge + 1;
        }

        // we have not visited that node before
        if (!visited.at(nodeId)) {
            visited.at(nodeId) = true;
            seqPos += isFirst ? 0 : 1;
            isFirst = false;
            val = get_node_end_value(nodeId);
            sequence.append(1, get_alphabet_symbol(val % alph_size));
            // store seq position of this node (we will join to it later)
            if (indegree(nodeId) > 1) {
                nodeId2seqPos.insert(std::make_pair(nodeId, std::make_pair(seqId, seqPos)));
            }
        }

        // there is only one child
        if (out == 1) {
            uint64_t next = fwd(nodeId);
            // the next node is new
            if (!visited.at(next)) {
                if (joinOpen) {
                    if (sequence.length() > 0) {
                        finish_sequence(sequence, seqCnt++, SQLstream);
                    }
                    seqId = seqCnt;
                    seqPos = 0;
                    branchPos = nodeId2seqId[nodeId];
                    joins.push_back(JoinInfo(branchPos.first, branchPos.second, seqId, seqPos));
                    joinOpen = false;
                    branchMap.insert(std::make_pair(std::make_pair(nodeId, joinEdge), joins.size() - 1));
                }
                nodeId = next;
                lastEdge = 0;
            // we have seen the next node before
            } else {
                // look up the sequence info of that node
                joins.push_back(JoinInfo(seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                branchMap.insert(std::make_pair(std::make_pair(nodeId, 1ul), joins.size() - 1));
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
                out = outdegree(nodeId);
                if (debug)
                    fprintf(stderr, " -- popped %lu -- ", nodeId); 
                joinOpen = true;
                joinEdge = lastEdge + 1;
            }
            if (debug)
                fprintf(stderr, " new nodeId: %lu\n", nodeId);
        // there are several children
        } else {
            size_t cnt = 0;
            bool updated = false;
            for (TAlphabet c = 1; c < alph_size; ++c) {
                uint64_t next = outgoing(nodeId, c);
                if (next > 0) {
                    cnt++;
                    // we already handled this edge erlier
                    if (cnt <= lastEdge)
                        continue;

                    lastEdge++;
                    if (!visited.at(next)) {
                        // there are remaining branches - push node to stack
                        if (cnt < out && next != nodeId) {
                            // each node can branch from exactly one sequence ID
                            if (nodeId2seqId.find(nodeId) == nodeId2seqId.end())
                                nodeId2seqId.insert(std::make_pair(nodeId, std::make_pair(seqId, seqPos)));
                            // push node information to stack
                            branchnodes.push(BranchInfo(nodeId, seqId, seqPos, lastEdge));
                            if (debug)
                                fprintf(stderr, " -- pushed %lu : seqId %lu seqPos %lu lastEdge %lu -- ", nodeId, seqId, seqPos, lastEdge); 
                        }
                        if (joinOpen) {
                            if (sequence.length() > 0) {
                                finish_sequence(sequence, seqCnt++, SQLstream);
                            }
                            seqId = seqCnt;
                            seqPos = 0;
                            branchPos = nodeId2seqId[nodeId];
                            joins.push_back(JoinInfo(branchPos.first, branchPos.second, seqId, seqPos));
                            joinOpen = false;
                            branchMap.insert(std::make_pair(std::make_pair(nodeId, lastEdge), joins.size() - 1));
                        }

                        nodeId = next;
                        updated = true;
                        lastEdge = 0;
                        break;
                    } else {
                        // look up the sequence info of that node
                        if (nodeId == next) { 
                            joins.push_back(JoinInfo(nodeId2seqPos[next].first, nodeId2seqPos[next].second, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                        } else {
                            joins.push_back(JoinInfo(seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                        }
                        branchMap.insert(std::make_pair(std::make_pair(nodeId, lastEdge), joins.size() - 1));
                    }
                }
            }
            // we are done with this branch
            // we should end up here, when nodes branch to themselves with their last edge
            if (!updated) {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, seqId, seqPos, nodeId, lastEdge, isFirst);
                out = outdegree(nodeId);
                if (debug)
                    fprintf(stderr, " -- popped %lu -- ", nodeId); 
                joinOpen = true;
                joinEdge = lastEdge + 1;
            }
            //if (debug)
            //    fprintf(stderr, " new nodeId: %lu\n", nodeId);
        }
        out = outdegree(nodeId);
    }
    // for completeness
    if (sequence.length() > 0)
        finish_sequence(sequence, seqCnt++, SQLstream);
    else
        seqCnt--;

    if (debug) {
        // output joins
        for (size_t i = 0; i < joins.size(); ++i) {
            std::cout << "(" << joins.at(i).seqId1 << ":" << joins.at(i).seqPos1 << "--" << joins.at(i).seqId2 << ":" << joins.at(i).seqPos2 << ")" << std::endl;
        }
        // output branches for tracking
        for (std::map<std::pair<uint64_t, uint64_t>, uint64_t>::iterator it = branchMap.begin(); it != branchMap.end(); ++it) {
            std::cout << "[" << (*it).first.first << ", " << (*it).first.second << "] -- " << (*it).second << std::endl;
        }
    }

    return seqCnt;
}

void DBG_succ::allelesFromSeq(kstring_t &seq, unsigned int f, std::vector<JoinInfo> &joins, std::map<std::pair<uint64_t, TAlphabet>, uint64_t> &branchMap, std::ofstream &SQLstream, bool isRefRun, size_t seqNum) {
    
    uint64_t nodeId = 1;
    uint64_t seqId = 1;
    uint64_t seqPos = 0;
    uint64_t alleleSeqPos = 0;
    // iterate over nodes
    uint64_t out = 0;
    uint64_t edge = 0;
    uint64_t currStart = 0;
    unsigned int alleleCnt = 1;
    TAlphabet seqVal, seqValNext;
    TAlphabet nodeVal;
    JoinInfo currJoin;
    std::vector<bool> isRef;

    if (!isRefRun) {
        SQLstream << "INSERT INTO Allele VALUES (" << f+1 << ", 1, NULL);" << std::endl;
    } else {
        for (size_t i = 0; i < seqNum; ++i)
            isRef.push_back(false);
    }
    if (config->verbose)
        fprintf(stderr, "processing alleles for file %u\n", f);

    while (true) {
        nodeVal = get_node_end_value(nodeId);
        //fprintf(stderr, "nodeId %lu nodeVal %lu seqVal %u\n", nodeId, nodeVal, ordValue(seq[seqPos]) + 1);
        if (nodeVal == 0) {
            nodeId = fwd(nodeId);
            alleleSeqPos++;
            currStart++;
            continue;
        }
        seqVal = get_alphabet_number(seq.s[seqPos]);
            
        //fprintf(stderr, "nodeVal %lu seqVal %lu nodeId %lu seqId %lu seqPos %lu alleleSeqPos %lu\n", nodeVal, seqVal, nodeId, seqId, seqPos, alleleSeqPos);
        assert(nodeVal % alph_size == 6 || nodeVal % alph_size == seqVal);
        if (seqPos + 1 == seq.l)
            break;
        if (nodeVal % alph_size != 6) {
            seqPos++;
        } else {
            currStart++;
        }
        seqValNext = get_alphabet_number(seq.s[seqPos]) + 1;

        // find edge to next node
        out = outdegree(nodeId);
        //fprintf(stderr, "out %lu nodeId %lu\n", out, nodeId);
        if (out == 1) {
            if (branchMap.find(std::make_pair(nodeId, out)) != branchMap.end()) {
                currJoin = joins.at(branchMap[std::make_pair(nodeId, out)]);
                //fprintf(stderr, "cseqId1 %lu cseqId2 %lu cseqPos1 %lu cseqPos2 %lu seqId %lu alleleSeqPos %lu\n", currJoin.seqId1, currJoin.seqId2, currJoin.seqPos1, currJoin.seqPos2, seqId, alleleSeqPos);
                assert(currJoin.seqId1 == seqId);
                assert(currJoin.seqPos1 == alleleSeqPos);
                if (!isRefRun) {
                    if (alleleSeqPos - currStart + 1 > 0) {
                        SQLstream << "INSERT INTO AllelePathItem VALUES (" << f+1 << ", " << alleleCnt++ << ", " << seqId << ", ";
                        SQLstream << currStart << ", " << alleleSeqPos - currStart + 1 << ", 'TRUE')" << std::endl;
                    }
                } else {
                    isRef[seqId - 1] = true;                 
                }
                seqId = currJoin.seqId2;
                currStart = alleleSeqPos = currJoin.seqPos2;
                nodeId = outgoing(nodeId, seqValNext);
            } else {
                nodeId = fwd(nodeId);
                alleleSeqPos++;
            }
        } else {
            // find edge corresponding to the proceding seqPos
            uint64_t start = pred_last(nodeId - 1);
            uint64_t stop = succ_last(nodeId);
            assert(stop - start == out);
            edge = 0;
            size_t k;
            for (k = start + 1; k <= stop; ++k) {
                if ((*W)[k] % alph_size == seqValNext) {
                    edge = k - start;
                    break;
                }
            }
            if (nodeVal % alph_size == 6)
                edge--;
            //fprintf(stderr, "seqValnext %lu edge %lu\n", seqValNext, edge);
            assert(edge > 0);
            if (branchMap.find(std::make_pair(nodeId, edge)) != branchMap.end()) {
                currJoin = joins.at(branchMap[std::make_pair(nodeId, edge)]);
                //fprintf(stderr, "cseqId1 %lu cseqId2 %lu cseqPos1 %lu cseqPos2 %lu seqId %lu seqPos %lu\n", currJoin.seqId1, currJoin.seqId2, currJoin.seqPos1, currJoin.seqPos2, seqId, alleleSeqPos);
                assert(currJoin.seqId1 == seqId);
                assert(currJoin.seqPos1 == alleleSeqPos);
                if (!isRefRun) {
                    if (alleleSeqPos - currStart + 1 > 0) {
                        SQLstream << "INSERT INTO AllelePathItem VALUES (" << f+1 << ", " << alleleCnt++ << ", " << seqId << ", ";
                        SQLstream << currStart << ", " << alleleSeqPos - currStart + 1 << ", 'TRUE')" << std::endl;
                    }
                } else {
                    isRef[seqId - 1] = true;
                }
                seqId = currJoin.seqId2;
                currStart = alleleSeqPos = currJoin.seqPos2;
                nodeId = outgoing(nodeId, seqValNext);
            } else {
                nodeId = fwd(k);
                alleleSeqPos++;
            }
        }
    }
    if (!isRefRun) {
        if (alleleSeqPos - currStart + 1 > 0) {
            SQLstream << "INSERT INTO AllelePathItem VALUES (" << f+1 << ", " << alleleCnt++ << ", " << seqId << ", ";
            SQLstream << currStart << ", " << alleleSeqPos - currStart + 1 << ", 'TRUE')" << std::endl;
        }
    } else {
        isRef[seqId - 1] = true;
        // write reference information to SQL stream 
        for (size_t i = 0; i < isRef.size(); ++i) {
            if (isRef[i])
                SQLstream << "INSERT INTO Reference VALUES (" << i + 1 << ", 'seq" << i + 1 << "', date('now'), " << i + 1 << ", 0, NULL, NULL, NULL, NULL, NULL, 'TRUE');" << std::endl;
            else
                SQLstream << "INSERT INTO Reference VALUES (" << i + 1 << ", 'seq" << i + 1 << "', date('now'), " << i + 1 << ", 0, NULL, NULL, NULL, NULL, NULL, 'FALSE');" << std::endl;
        }
        SQLstream << "INSERT INTO ReferenceSet VALUES (1, NULL, NULL, 'normal', 'FALSE');" << std::endl;
        for (size_t i = 0; i < isRef.size(); ++i)
            SQLstream << "INSERT INTO Reference_ReferenceSet_Join VALUES (" << i + 1 << ", 1);" << std::endl;
        for (size_t i = 0; i < joins.size(); ++i)
            SQLstream << "INSERT INTO GraphJoin_ReferenceSet_Join VALUES (" << i + 1 << ", 1);" << std::endl;
    }

}


//
//
// ANNOTATE
//
//

void DBG_succ::annotate_kmer(std::string &kmer, std::string &label) {

    // get index of the k-mer we would like to annotate
    uint64_t idx = this->index(kmer);
    if (idx == 0)
        return;

    // get annotation of current kmer
    size_t curr_hash = this->annotation[idx];
    std::set<std::string> curr_anno;
    if (curr_hash > 0) {
        curr_anno = annotation_map.at(curr_hash);    
    }

    // check whether current label is already part of current annotation
    // and add it if not
    std::set<std::string>::const_iterator it = curr_anno.find(label); 
    if (it == curr_anno.end()) {
        curr_anno.insert(label);
        curr_hash = AnnotationHash{}(curr_anno);
        annotation_map[curr_hash] = curr_anno;
        annotation[idx] = curr_hash;
    }
}

//void annotate_kmers(rocksdb::DB* db, kstring_t &seq, kstring_t &label) {
void DBG_succ::annotate_kmers(kstring_t &seq, kstring_t &label) {

    std::string curr_kmer;
    std::string label_str = std::string(label.s);
    //rocksdb::Status s;

    for (size_t i = 0; i < seq.l; ++i) {

        if (config->verbose && i > 0 && i % 1000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10000 == 0)
                std::cout << i << " kmers added" << std::endl;
        }

        if (i < k) {
            curr_kmer.push_back(seq.s[i]);
            continue;
        }

        annotate_kmer(curr_kmer, label_str);
        //db->Put(rocksdb::WriteOptions(), curr_kmer, label_str); 
        
        //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
        curr_kmer.push_back(seq.s[i]);
        curr_kmer = curr_kmer.substr(1, config->k);
    }
    // add last kmer and label to database
    //db->Put(rocksdb::WriteOptions(), curr_kmer, label_str); 
    annotate_kmer(curr_kmer, label_str);
    //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
}   



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
void DBG_succ::merge(DBG_succ* G) {

    // store all branch nodes on the way
    std::stack<BranchInfoMerge> branchnodes;
    // bool vector that keeps track of visited nodes
    std::vector<bool> visited(G->get_size());
    for (std::vector<bool>::iterator it = visited.begin(); it != visited.end(); ++it) {
        *it = false;
    }

    // some initializations
    uint64_t nodeId = 1; // start at source node
    size_t out = G->outdegree(nodeId);
    BranchInfoMerge branch;
    TAlphabet val;
    TAlphabet lastEdge = 0;
    // keep a running list of the last k-1 characters we have seen
    std::deque<TAlphabet> last_k;
    bool new_branch = false;
    bool old_last = (*last)[p];
    bool initial_k = true;
    uint64_t added = 0;

    // keep traversing until we reach the think and have worked off all branches from the stack
    while (out > 0 || branchnodes.size() > 0) {

        // verbose output
        if (added > 0 && added % 1000 == 0) {
            std::cout << "." << std::flush;
            if (added % 10000 == 0) {
                fprintf(stdout, "merged %lu / %lu - edges %lu / nodes %lu\n", added, G->get_size(), W->n - 1, rank_last((last->size() - 1)));
            }
        }

        // we have reached the sink but there are unvisited nodes left on the stack
        if (out == 0) {
            if (branchnodes.size() == 0)
                break;
            // get new branch
            branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
            out = G->outdegree(nodeId);
            new_branch = true;
        }
        //std::cerr << "starting loop with nodID " << nodeId << std::endl;

        if (new_branch) {
            // find node where to restart insertion
            uint64_t ridx = this->index(last_k);
            // put at the beginning of equal node range
            ridx = this->pred_last(ridx - 1) + 1;
            ridx -= (p < ridx);
            //std::cerr << "ridx: " << ridx << std::endl;
            //std::cerr << "bef move " << std::endl;
            //this->print_seq();
            assert(!(*last)[p]); // must not remove a dangling end
            // move p to last position
            W->remove(p);
            last->deleteBit(p);
            update_F(this->get_node_end_value(p), false);
            //std::cerr << "moving p from " << p << " to " << ridx << std::endl;
            p = ridx;
            W->insert(0, p);
            last->insertBit(p, 0); // it's at the beginning of a branch node
            update_F(this->get_node_end_value(p), true);

            new_branch = false;
            //std::cerr << "aft move " << std::endl;
            //this->print_seq();
        }

        // we have not visited that node before
        if (!visited.at(nodeId)) {
            visited.at(nodeId) = true;
            val = G->get_W(nodeId) % alph_size;
            //std::cerr << "current val " << val % alph_size << " nodeID: " << nodeId << std::endl;
            //G->print_seq();
            last_k.push_back(G->get_node_end_value(nodeId));
            if (last_k.size() > k)
                last_k.pop_front();
        }

        // there is only one child
        if (out == 1) {
            uint64_t next = G->fwd(nodeId);
            val = G->get_W(nodeId) % alph_size;
            if ((val != 6 || !initial_k) && val != 0) {
                initial_k = false;
                this->append_pos(val % alph_size);
                added++;
                //std::cerr << "append " << val % alph_size << " nodeID: " << nodeId << std::endl;
                //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                //this->print_seq();
            }

            // the next node is new
            if (!visited.at(next)) {
                nodeId = next;
                lastEdge = 0;
            // we have seen the next node before
            } else {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // append next node
                if ((*last)[p]) {
                    val = G->get_W(next) % alph_size;
                    if ((val != 6 || !initial_k) && val != 0) {
                        initial_k = false;
                        this->append_pos(val % alph_size);
                        added++;
                        //std::cerr << "..append " << val % alph_size << " nodeID: " << nodeId << std::endl;
                        //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                    }
                }
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = G->outdegree(nodeId);
                //std::cerr << "new branch 1 - nodeID: " << nodeId << " next is " << next << std::endl;
                new_branch = true;
            }
        // there are several children
        } else {
            size_t cnt = 0;
            bool updated = false;
            // account for sentinel symbol as possible outgoing edge
            if (G->get_W(nodeId) == 0)
                cnt++;
            // loop over outgoing edges
            for (TAlphabet c = 1; c < alph_size; ++c) {
                uint64_t next = G->outgoing(nodeId, c);
                if (next > 0) {
                    cnt++;
                    // we already handled this edge erlier
                    if (cnt <= lastEdge)
                        continue;
                    lastEdge++;

                    //std::cerr << "mult edge - val is now " << c << std::endl;
                    uint64_t curr_p = p;
                    if ((c != 6 || !initial_k) && c != 0) {
                        initial_k = false;
                        //std::cerr << "p (before): " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                        this->append_pos(c);
                        added++;
                        curr_p += (p <= curr_p);
                        //std::cerr << "append " << c % alph_size << " nodeID: " << nodeId << std::endl;
                        //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                        //this->print_seq();
                    }

                    if (!visited.at(next)) {
                        // there are remaining branches - push node to stack
                        if (cnt < out && next != nodeId) {
                            // push node information to stack
                            branchnodes.push(BranchInfoMerge(nodeId, lastEdge, last_k));
                            //std::cerr << "pushing nodeID " << nodeId << " onto stack" << std::endl;
                        }
                        nodeId = next;
                        updated = true;
                        lastEdge = 0;
                        break;
                    } else {
                        //std::cerr << "visited next before: " << next <<std::endl;
                        // append next node
                        if ((*last)[p]) {
                            c = G->get_W(next) % alph_size;
                            if ((c != 6 || !initial_k) && c != 0) {
                                initial_k = false;
                                this->append_pos(c);
                                added++;
                                //std::cerr << "...append " << c % alph_size << " nodeID: " << nodeId << std::endl;
                                //std::cerr << "p: " << p << " W size: " << W->n << std::endl;
                                //this->print_seq();
                            }
                        }
                        // reset to previous position
                        if (nodeId != next) {
                            W->remove(p);
                            last->deleteBit(p);
                            update_F(this->get_node_end_value(p), false);
                            curr_p -= (p < curr_p);
                            //std::cerr << ".moving p from " << p << " to " << curr_p << std::endl;
                            p = curr_p;
                            W->insert(0, p);
                            last->insertBit(p, 0); // it's at the beginning of a branch node
                            update_F(this->get_node_end_value(p), true);
                            //std::cerr << ".aft move " << std::endl;
                            //this->print_seq();
                        }
                    }
                }
            }
            // we are done with this branch
            // we should end up here, when nodes branch to themselves with their last edge
            if (!updated) {
                // there are no branches left
                if (branchnodes.size() == 0)
                    break;
                // otherwise go back to last branch
                branch = pop_branch(branchnodes, nodeId, lastEdge, last_k);
                out = G->outdegree(nodeId);
                new_branch = true;
                //std::cerr << "new branch 2" << std::endl;
            }
        }
        out = G->outdegree(nodeId);
    }
    // bring graph into default state
    std::deque<TAlphabet> tmp_p;
    for (size_t t = 0; t < k; t++)
        tmp_p.push_back(6);
    uint64_t old_p = pred_last(this->index(tmp_p) - 1) + 1;

    old_p -= (p < old_p);
    W->remove(p);
    last->deleteBit(p);
    update_F(this->get_node_end_value(p), false);
    p = old_p;
    W->insert(0, p);
    last->insertBit(p, old_last);

    // locally update sorting
    std::pair<uint64_t, uint64_t> R = get_equal_node_range(this->p);
    if (R.second - R.first > 0) {
        sort_W_locally(R.first, R.second);
        while ((*W)[p] != 0)
            p--;
        assert((*W)[p] == 0);
    }
    update_F(this->get_node_end_value(p), true);
}

/* 
 * Helper function to determine the bin boundaries, given 
 * a number of threads.
 */
std::vector<std::pair<uint64_t, uint64_t> > DBG_succ::get_bins(uint64_t threads, uint64_t bins_per_thread, DBG_succ* G) {

    uint64_t binlen;
    uint64_t bins = threads * bins_per_thread;

    // depending on the number of threads, we will use a different length prefix
    // to compute the bin boundaries
    if (bins < alph_size)
        binlen = 1;
    else if (bins < (alph_size * alph_size) || (k <= 2))
        binlen = 2;
    else if (bins < (alph_size * alph_size * alph_size) || (k <= 3))
        binlen = 3;
    else if (bins < (alph_size * alph_size * alph_size * alph_size) || (k <= 4))
        binlen = 4;
    else if (bins < (alph_size * alph_size * alph_size * alph_size * alph_size) || (k <= 5))
        binlen = 5;
    else
        binlen = 6;

    std::vector<std::pair<uint64_t, uint64_t> > tmp;
    for (uint64_t i = 0; i < std::pow(alph_size, binlen); ++i) {
        //std::deque<TAlphabet> ttt = bin_id_to_string(i, binlen);
        //std::cerr << std::endl;
        //for (size_t ii = 0; ii < ttt.size(); ++ii) {
        //    std::cerr << ttt[ii] << get_alphabet_symbol(ttt[ii] % alph_size);
        //}
        //std::cerr << std::endl;
        std::pair<uint64_t, uint64_t> idx = G->index_range(bin_id_to_string(i, binlen));
        if (idx.first > idx.second)
            idx.first = idx.second = 0;
        if (idx.second >= G->get_size())
            idx.second = G->get_size() - 1;
        if (idx.first > 0)
            idx.first = G->pred_last(idx.first - 1) + 1;
        tmp.push_back(idx);
    }
    return tmp;
}

std::deque<TAlphabet> DBG_succ::bin_id_to_string(uint64_t bin_id, uint64_t binlen) {
    std::deque<TAlphabet> str;
    while (str.size() < binlen) {
        str.push_back(bin_id % alph_size);
        bin_id -= (bin_id % alph_size);
        bin_id /= alph_size;
    }
    return str;
}


/*
 * Given two pointers to graph structures G1 and G2, this function 
 * integrate both into a new graph G.
 */
void DBG_succ::merge(DBG_succ* G1, DBG_succ* G2, uint64_t k1, uint64_t k2, uint64_t n1, uint64_t n2, bool is_parallel) {

    // check whether we can merge the given graphs
    if (G1->get_k() != G2->get_k()) {
        fprintf(stderr, "Graphs have different k-mer lengths - cannot be merged!\n");
        exit(1);
    }
    
    // positions in the graph for respective traversal
    n1 = (n1 == 0) ? G1->get_size() : n1;
    n2 = (n2 == 0) ? G2->get_size() : n2;

    // handle special cases where one or both input graphs are empty
    k1 = (k1 == 0) ? G1->get_size() : k1;
    k2 = (k2 == 0) ? G2->get_size() : k2;

    // keep track of how many nodes we added and from which graph the 
    // last node originated
    uint64_t added = 0;
    uint64_t last_added_k = 0;
    DBG_succ* last_added_G = NULL;

    bool k1_smaller, k2_smaller, advance_k1, advance_k2, insert_k1, insert_k2, identical;

    std::string k1_str, k2_str;
    
    // Send two pointers running through each of the two graphs, where k1 runs through 
    // G1 and k2 through G2. At each step, compare the two graph nodes at positions 
    // k1 and k2 with each other. Insert the lexicographically smaller one into the 
    // common merge graph G. 
    while (k1 < n1 || k2 < n2) {

        if (!is_parallel && config->verbose && added > 0 && added % 1000 == 0) {
            std::cout << "." << std::flush;
            if (added % 10000 == 0) {
                fprintf(stdout, "added %lu - G1: edge %lu/%lu - G2: edge %lu/%lu\n", added, k1, G1->get_size(), k2, G2->get_size());
            }
        }

        k1_smaller = false;
        k2_smaller = false;
        advance_k1 = false;
        advance_k2 = false;
        insert_k1 = false;
        insert_k2 = false;
        identical = false;

        //G1->print_seq();
        //G2->print_seq();
        if (k1 < n1 && k2 < n2) {
            std::pair<bool, bool> tmp = compare_nodes(G1, k1, G2, k2);
            k1_smaller = tmp.first;
            k2_smaller = tmp.second;
        } else {
            k1_smaller = k1 < n1;
            k2_smaller = k2 < n2;
        }

        // the node sequences are identical
        if (!k1_smaller && !k2_smaller) {
            // insert G1 and advance both pointers
            if ((G1->get_W(k1) % alph_size) == (G2->get_W(k2) % alph_size)) {
                advance_k1 = true;
                advance_k2 = true;
                insert_k1 = true;
                identical = true;
                //std::cerr << "a) identical - advance both" << std::endl;
            // insert G1 and advance k1
            } else if ((G1->get_W(k1) % alph_size) < (G2->get_W(k2) % alph_size)) {
                advance_k1 = true;
                insert_k1 = true;
                //std::cerr << "b) inserting " << G1->get_W(k1) << " from position " << k1 << std::endl;
            // insert G2 and advance k2
            } else {
                advance_k2 = true;
                insert_k2 = true;
                //std::cerr << "b) inserting " << G2->get_W(k2) << " from position " << k2 << std::endl;
            }
        // the node in graph 1 is smaller
        } else if (k1_smaller) {
            // insert G1 and advance k1
            advance_k1 = true;
            insert_k1 = true;
            //std::cerr << "d) inserting " << G1->get_W(k1) << " from position " << k1 << std::endl;
        // the node in graph 2 is smaller
        } else if (k2_smaller) {
            // insert G2 and advance k2
            advance_k2 = true;
            insert_k2 = true;
            //std::cerr << "e) inserting " << G2->get_W(k2) << " from position " << k2 << std::endl;
        } else {
            std::cerr << "This should never happen." << std::endl;
            std::exit(1);
        }

        // assert XOR of insert; we take a node from either G1 or G2
        assert(insert_k1 ^ insert_k2);

        // insert node from G1
        if (insert_k1) {
            TAlphabet val = G1->get_W(k1);
            if (identical) {
                W->insert(std::max(val, G2->get_W(k2)), W->n);
            } else if (val < alph_size) {
                std::deque<TAlphabet> seq1 = G1->get_node_seq(k1);
                seq1.pop_front();
                seq1.push_back(val);
                uint64_t kidx2 = G2->index(seq1);
                if (kidx2 > 0 && kidx2 < G2->get_size()) {
                    uint64_t kidx3 = G2->bwd(kidx2);
                    //std::cerr << "kidx3 " << kidx3 << std::endl;
                    if (G2->get_minus_k_value(kidx3, k-1) < G1->get_minus_k_value(k1, k-1)) {
                        W->insert(val + alph_size, W->n);
                    } else {
                        W->insert(val, W->n);
                    }
                    //std::cerr << "c" << std::endl;
                } else {
                    W->insert(val, W->n);
                    //std::cerr << "d" << std::endl;
                }
            } else {
                W->insert(val, W->n);
                //std::cerr << "e" << std::endl;
            }
            update_F(G1->get_node_end_value(k1), true);
        } 

        // insert node from G2
        if (insert_k2) {
            TAlphabet val = G2->get_W(k2);
            if (val % alph_size == val) {
                std::deque<TAlphabet> seq1 = G2->get_node_seq(k2);
                seq1.pop_front();
                seq1.push_back(val);
                uint64_t kidx2 = G1->index(seq1);
                //std::cerr << " kidx2 " << kidx2 << std::endl; 
                if ((kidx2 > 0) && (kidx2 < G1->get_size())) {
                    uint64_t kidx3 = G1->bwd(kidx2);
                    //std::cerr << "kidx3 " << kidx3 << std::endl;
                    if (G1->get_minus_k_value(kidx3, k-1) < G2->get_minus_k_value(k2, k-1)) {
                        W->insert(val + alph_size, W->n);
                    } else {
                        W->insert(val, W->n);
                    }
                    //std::cerr << "h" << std::endl;
                } else {
                    W->insert(val, W->n);
                    //std::cerr << "i" << std::endl;
                }
            } else {
                W->insert(val, W->n);
                //std::cerr << "j" << std::endl;
            }
            update_F(G2->get_node_end_value(k2), true);
        }
        if ((insert_k1 && !G1->get_last(k1)) || (insert_k2 && !G2->get_last(k2))) {
            //std::cerr << "insert_k1 " << insert_k1 << " last_k1 " << G1->get_last(k1) << "insert_k2 " << insert_k2 << " last_k2 " << G2->get_last(k2) << std::endl;
            last->insertBit(W->n - 1, false);
        } else {
            last->insertBit(W->n - 1, true);
        }

        if (last_added_k > 0 && W->n > 2 && (*last)[W->n-2]) {
            // compare the last two added nodes
            std::pair<bool, bool> tmp = insert_k1 ? compare_nodes(G1, k1, last_added_G, last_added_k) : compare_nodes(G2, k2, last_added_G, last_added_k); 
            if (!tmp.first && !tmp.second) {
                last->set(W->n - 2, false);
            }
        }
        last_added_k = insert_k1 ? k1 : k2;
        last_added_G = insert_k1 ? G1 : G2;
        k1 += advance_k1;
        k2 += advance_k2;

        ++added;
    }
    p = succ_W(1, 0);
}


/**
* Given a pointer to a graph structure G, the function compares its elements to the
* current graph. It will perform an element wise comparison of the arrays W, last and
* F and will only check for identity. If any element differs, the function will return 
* false and true otherwise.
*/
bool DBG_succ::compare(DBG_succ* G) {
    
    // compare size
    if (W->n != G->get_size())
        return false;
    
    // compare W
    for (size_t i = 0; i < W->n; ++i) {
        if ((*W)[i] != G->get_W(i)) {
            std::cerr << "W differs at position " << i << std::endl;
            return false;
        }
    }

    // compare last
    for (size_t i = 0; i < W->n; ++i) {
        if ((*last)[i] != G->get_last(i)) {
            std::cerr << "last differs at position " << i << std::endl;
            return false;
        }
    }

    // compare F
    for (size_t i = 0; i < F.size(); ++i) {
        if (F.at(i) != G->get_F(i)) {
            std::cerr << "F differs at position " << i << std::endl;
            return false;
        }
    }

    return true;

}

//
//
// SERIALIZE
//
//
    
/**
 * This is a debug function that prints the current state of the graph arrays to
 * the screen.
 */
void DBG_succ::print_state() {

    fprintf(stderr, "W:\n");
    for (uint64_t i = 0; i < W->n; i++) {
        fprintf(stderr, "\t%lu", (*W)[i]);
        if (i == p)
            fprintf(stderr, "*");
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "last:\n");
    for (uint64_t i = 0; i < last->size(); i++)
        fprintf(stderr, "\t%i", (int) (*last)[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "F:\n");
    for (uint64_t i = 0; i < F.size(); i++)
        fprintf(stderr, "\t%i", (int) F[i]);
    fprintf(stderr, "\n");

}


/*
 * Returns the sequence stored in W and prints the node
 * information in an overview. 
 * Useful for debugging purposes.
 */
void DBG_succ::print_seq() {

    uint64_t linelen = 80;
    uint64_t start = 1;
    uint64_t end = start + linelen < W->n ? start + linelen : W->n;

    while (start < W->n) {
        for (uint64_t i = start; i < end; i++) {
            if (i % 10 == 0)
                fprintf(stdout, "%lu", (i / 10) % 10);
            else
                fprintf(stdout, " ");
        }
        std::cout << std::endl;

        for (uint64_t i = start; i < end; i++) {
            if ((*W)[i] >= alph_size)
                fprintf(stdout, "-");
            else
                fprintf(stdout, " ");
        }
        std::cout << std::endl;

        for (uint64_t i = start; i < end; i++) {
            if ((*W)[i] % alph_size == 0)
                fprintf(stdout, "$");
            else
                std::cout << get_alphabet_symbol((*W)[i] % alph_size);
        }
        std::cout << std::endl;

        for (uint64_t i = start; i < end; i++) {
            if (p == i)
                fprintf(stdout, "*");
            else
                fprintf(stdout, " ");
        }
        std::cout << std::endl;

        size_t j;
        for (size_t l = 0; l < k; l++) {
            for (uint64_t i = start; i < end; i++) {
                j = get_minus_k_value(i, l).first;
                if (j % alph_size == 0)
                    std::cout << "$";
                else
                    std::cout << get_alphabet_symbol(j % alph_size);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        for (uint64_t i = start; i < end; i++) {
            fprintf(stdout, "%i", (int) (*last)[i]);
        }
        std::cout << std::endl;
        std::cout << std::endl;

        for (uint64_t i = start; i < end; ++i) {
            std::cout << indegree(i);  
        }
        std::cout << std::endl;
        for (uint64_t i = start; i < end; ++i) {
            std::cout << outdegree(i);  
        }
        std::cout << std::endl;
        std::cout << std::endl;

        start += linelen;
        end = start + linelen < W->n ? start + linelen : W->n;
    }
}

/**
 * Take the current graph content and return it in SQL
 * format (GA4GH Spec).
 *
 * We will perform one depth first search traversal of the graph. While we will record
 * one long reference string, we will output all sidepaths on the way.
 */
void DBG_succ::toSQL() {
    
    // this vector stores the joins between the sequence objects we wrote
    std::vector<JoinInfo> joins;
    // we also store for each branching edge the join it creates.
    // we will use this for allele traversal
    std::map<std::pair<uint64_t, TAlphabet>, uint64_t> branchMap;

    // open sql filestream
    std::ofstream SQLstream;
    SQLstream.open((config->sqlfbase + ".sql").c_str());

    // traverse the graph, thereby filling joins vector, branchMap and 
    // writing the sequences to individual fasta files
    size_t seqNum = traverseGraph(joins, branchMap, SQLstream); 
    
    // write graph joins to SQL file
    for (size_t i = 0; i < joins.size(); ++i) {
        if (joins.at(i).seqId1 < joins.at(i).seqId2 || (joins.at(i).seqId1 == joins.at(i).seqId2 && joins.at(i).seqPos1 < joins.at(i).seqPos2)) {
            SQLstream << "INSERT INTO GraphJoin VALUES (" << i + 1 << ", " << joins.at(i).seqId1 << ", " << joins.at(i).seqPos1 << ", 'FALSE', " 
                                                          << joins.at(i).seqId2 << ", " << joins.at(i).seqPos2 << ", 'TRUE');" << std::endl;
        } else {
            SQLstream << "INSERT INTO GraphJoin VALUES (" << i + 1 << ", " << joins.at(i).seqId2 << ", " << joins.at(i).seqPos2 << ", 'TRUE', " 
                                                          << joins.at(i).seqId1 << ", " << joins.at(i).seqPos1 << ", 'FALSE');" << std::endl;
        }
    }

    // for each input sequence traverse the graph once more and
    // collect allele path information 
    for (unsigned int f = 0; f < config->fname.size(); ++f) {

        // first traversal is for reference info
        if (f == 0) {
            // open stream to fasta file
            gzFile input_p = gzopen(config->fname.at(f).c_str(), "r");
            kseq_t *stream = kseq_init(input_p);

            if (stream != NULL)
                std::cerr << "ERROR while opening input file " << config->fname.at(f) << std::endl;

            while (kseq_read(stream) >= 0) {
                allelesFromSeq(stream->seq, f, joins, branchMap, SQLstream, true, seqNum);
            }
            kseq_destroy(stream);
            gzclose(input_p);

            // open variant set
            SQLstream << "INSERT INTO VariantSet VALUES (1, 1, 'deBruijnGraph');" << std::endl;
        }
        // open stream to fasta file
        gzFile input_p = gzopen(config->fname.at(f).c_str(), "r");
        kseq_t *stream = kseq_init(input_p);
        if (stream != NULL)
            std::cerr << "ERROR while opening input file " << config->fname.at(f) << std::endl;

        while (kseq_read(stream) >= 0) 
            allelesFromSeq(stream->seq, f, joins, branchMap, SQLstream);

        kseq_destroy(stream);
        gzclose(input_p);
    }

    // write call set (one per input file)
    for (unsigned int f = 0; f < config->fname.size(); ++f)
        SQLstream << "INSERT INTO CallSet VALUES (" << f+1 <<", '" << config->fname.at(f) << "', 'DBG" << f + 1 << "');" << std::endl;
    for (unsigned int f = 0; f < config->fname.size(); ++f)
        SQLstream << "INSERT INTO VariantSet_CallSet_Join VALUES (1, " << f + 1 << ");" << std::endl;
    for (unsigned int f = 0; f < config->fname.size(); ++f) {
        for (unsigned int ff = 0; ff < config->fname.size(); ++ff) {
            if (f == ff)
                SQLstream << "INSERT INTO AlleleCall VALUES (" << ff + 1 << ", " << f + 1 << ", 1);" << std::endl;
            else
                SQLstream << "INSERT INTO AlleleCall VALUES (" << ff + 1 << ", " << f + 1 << ", 0);" << std::endl;

        }
    }
}


/**
 * Take the current graph content and store in a file.
 *
 */
void DBG_succ::toFile() {

    // write Wavelet Tree
    std::ofstream outstream((config->outfbase + ".W.dbg").c_str());
    W->serialise(outstream);
    outstream.close();

    // write last array
    outstream.open((config->outfbase + ".l.dbg").c_str());
    last->serialise(outstream);
    outstream.close();

    // write F values and k
    outstream.open((config->outfbase + ".F.dbg").c_str());
    outstream << ">F" << std::endl;
    for (size_t i = 0; i < F.size(); ++i)
        outstream << F.at(i) << std::endl;
    outstream << ">k" << std::endl;
    outstream << k << std::endl;
    outstream << ">p" << std::endl;
    outstream << p << std::endl;
    outstream.close();
}

// write annotation to disk
void DBG_succ::annotationToFile() {
    std::ofstream outstream((config->infbase + ".anno.dbg").c_str());
    libmaus2::util::NumberSerialisation::serialiseNumberShortDeque<uint16_t>(outstream, annotation);
    serialize_annotation_map(outstream, annotation_map);
    outstream.close();
}

// read annotation from disk
void DBG_succ::annotationFromFile() {
    std::ifstream instream((config->infbase + ".anno.dbg").c_str());
    if (instream.good()) {
        annotation = libmaus2::util::NumberSerialisation::deserialiseNumberShortDeque<uint16_t>(instream);
        annotation_map = deserialize_annotation_map(instream);
    } else {
        annotation.resize(get_size(), 0);
    }
    instream.close();
}
