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
#include <pthread.h>
#include <ctime>
#include <parallel/algorithm>

// use Heng Li's kseq structure for string IO
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
#include "dbg_succinct_boost.hpp"
#include "annotation.hpp"

// define an extended alphabet for W --> somehow this does not work properly as expected
typedef uint64_t TAlphabet;
typedef DBG_succ::BranchInfo BranchInfo;
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
//libmaus2::bitbtree::BitBTree<6, 64> *last = new libmaus2::bitbtree::BitBTree<6, 64>();

// the array containing the edge labels
//libmaus2::wavelet::DynamicWaveletTree<6, 64> *W = new libmaus2::wavelet::DynamicWaveletTree<6, 64>(4); // 4 is log (sigma)

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
    id_to_label.push_back("");
    combination_vector.push_back(0);
}



DBG_succ::DBG_succ(std::string infbase_, Config* config_) : 
    infbase(infbase_),
    config(config_) {

    // load last array
    std::ifstream instream((infbase + ".l.dbg").c_str());
    last->deserialise(instream);
    instream.close();

    // load W array
    delete W;
    instream.open((infbase + ".W.dbg").c_str());
    W = new WaveletTree(instream);
    //W = new libmaus2::wavelet::DynamicWaveletTree<6, 64>(instream);
    instream.close();

    // load F and k and p
    for (size_t j = 0; j < alph_size; j++)
        F.push_back(0);
    instream.open((infbase + ".F.dbg").c_str());
    std::string line;
    size_t mode = 0;
    size_t fidx = 0;
    while (std::getline(instream, line)) {
        if (strcmp(line.c_str(), ">F") == 0) {
            mode = 1;
        } else if (strcmp(line.c_str(), ">k") == 0) {
            mode = 2;
        } else if (strcmp(line.c_str(), ">p") == 0) {
            mode = 3;
        } else {
            if (mode == 1) {
                F.at(fidx) += std::strtoul(line.c_str(), NULL, 10);
                fidx++;
                //F.push_back(std::strtoul(line.c_str(), NULL, 10));
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
    id_to_label.push_back("");
    combination_vector.push_back(0);
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
    // once the end of the pattern is reached, add match to results
    
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
    std::deque<TAlphabet>::iterator it = str.begin();
    TAlphabet s = *it % alph_size;
    it++;
    // init range
    uint64_t rl = succ_last(F.at(s) + 1);
    uint64_t ru = (s < F.size() - 1) ? F.at(s + 1) : (W->n - 1);                // upper bound
    uint64_t pl;
    //fprintf(stderr, "char: %i rl: %i ru: %i\n", (int) s, (int) rl, (int) ru);
    // update range iteratively while scanning through s
    for (; it != str.end(); it++) {
        s = *it % alph_size;
        pl = pred_last(rl - 1) + 1;
        rl = std::min(succ_W(pl, s), succ_W(pl, s + alph_size));
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
 * Given a node label s, this function returns the index
 * of the corresponding node or the closest predecessor, if no node 
 * with the sequence is not found.
 */
uint64_t DBG_succ::index_predecessor(std::deque<TAlphabet> str) {
    // get first 
    std::deque<TAlphabet>::iterator it = str.begin();
    TAlphabet s1 = *it % alph_size;
    // init range
    uint64_t rl = succ_last(F.at(s1) + 1);                           // lower bound
    uint64_t ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);    // upper bound
    while (rl > ru && s1 > 0) {
        s1--;
        rl = succ_last(F.at(s1) + 1);
        ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);
    }
    /*if (s1 == 0) {
        s1 = *it % alph_size + 1;
        rl = succ_last(F.at(s1) + 1);
        ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);
        while (rl > ru && s1 < alph_size) {
            s1++;
            rl = succ_last(F.at(s1) + 1);
            ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->n - 1);
        }
    }*/
    //std::cerr << "s: " << s1 << " rl: " << rl << " ru: " << ru << std::endl;

    it++;
    uint64_t pll, puu;
    bool before = false;
    //fprintf(stderr, "char: %i rl: %i ru: %i\n", (int) s1, (int) rl, (int) ru);
    // update range iteratively while scanning through s
    for (; it != str.end(); it++) {
        s1 = *it % alph_size;
        pll = this->pred_last(rl - 1) + 1;
        puu = this->succ_last(ru);
        before = false;

        //std::cerr << "s: " << s1 << " rl: " << rl << " ru: " << ru << " pll: " << pll << std::endl;

        rl = std::min(succ_W(pll, s1), succ_W(pll, s1 + alph_size)); 
        ru = std::max(pred_W(puu, s1), pred_W(puu, s1 + alph_size));
        if (rl > puu) {
            rl = std::max(pred_W(pll, s1), pred_W(pll, s1 + alph_size));
            if (rl == 0) {
                rl = std::min(succ_W(pll, s1), succ_W(pll, s1 + alph_size));
                if (rl >= W->n) {
                    s1--;
                    while (s1 > 0) {
                        rl = std::max(pred_W(W->n - 1, s1), pred_W(W->n - 1, s1));
                        if (rl < W->n)
                            break;
                        s1--;
                    }
                    if (s1 == 0) {
                        s1 = (*it % alph_size) + 1;
                        before = true;
                        while (s1 < alph_size) {
                            rl = std::min(succ_W(1, s1), succ_W(1 + alph_size, s1));
                            if (rl < W->n)
                                break;
                            s1++;
                        }
                    }
                } else {
                    s1 = (*W)[rl];
                    before = true;
                }
            } else {
                s1 = (*W)[rl];
            }
        }

        if (ru == 0)
            ru = rl;

        //std::cerr << "endloop - rl: " << rl << " s1: " << s1 << " ru: " << ru << std::endl;
        assert(rl <= ru);

        rl = outgoing(rl, s1);
        ru = outgoing(ru, s1);
    }
    return pred_last(rl) - before;
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
 * Return number nodes in the current graph.
 */
uint64_t DBG_succ::get_nodes() {
    return this->rank_last(this->get_size() - 1);
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

void DBG_succ::add_seq_alt (kstring_t &seq, bool bridge, unsigned int parallel, std::string suffix) {

    if (debug) {
        print_seq();
        print_state();
        std::cout << "======================================" << std::endl;
    }

    char *nt_lookup = (char*)malloc(128);
    uint8_t def = 5;
    memset(nt_lookup, def, 128);
    for (size_t i=0;i<alph_size;++i) {
        nt_lookup[(uint8_t)alphabet[i]]=i;
        nt_lookup[(uint8_t)tolower(alphabet[i])]=i;
    }

    //clock_t start = clock();
    //std::cerr << "Loading kmers\n";
    if (!kmers.size()) {
        seqtokmer(kmers, "$", 1, k, nt_lookup);
        kmers.push_back(stokmer(std::string(k-1,'X')+std::string("$$"), nt_lookup));
    }
    seqtokmer(kmers, seq.s, seq.l, k, nt_lookup, bridge, parallel, suffix);
    //std::cerr << (clock()-start)/CLOCKS_PER_SEC << "\n";
    free(nt_lookup);
}

void DBG_succ::construct_succ(unsigned int parallel) {
    clock_t start=clock();
    std::cerr << "Sorting kmers\t";
    //omp_set_num_threads(std::max((int)parallel-1,1));
    //__gnu_parallel::sort(kmers.begin(),kmers.end());
    std::sort(kmers.begin(),kmers.end());
    kmers.erase(std::unique(kmers.begin(), kmers.end() ), kmers.end() );
    std::cerr << (clock()-start)/CLOCKS_PER_SEC << "\n";
    start=clock();
    std::cerr << "Constructing succinct representation\t";
    
    delete last;
    delete W;
    last = new BitBTree(kmers.size(), true);
    //sdsl::int_vector<> wbv(kmers.size(), 0, strlen(alphabet));
    
    std::vector<uint8_t> Wvec(kmers.size());
    size_t lastlet=0;
    std::cerr << "\n";
    //std::cerr << getPos(kmers[10],k-1,alphabet,alph_size) << "-" << getPos(kmers[10],k,alphabet,alph_size) << "\n";
    //assert(getPos(kmers[10],k,alphabet,alph_size) == 0);
    #pragma omp parallel num_threads(parallel)
    {
        #pragma omp for
        for (size_t i=0;i<kmers.size();++i) {
            //set last
            if (i+1 < kmers.size()) {
                bool dup = compare_kmer_suffix(kmers[i], kmers[i+1]);
                if (dup) {
                    #pragma omp critical
                    //last->setBitQuick(i, !compare_kmer_suffix(kmers[i], kmers[i+1]));
                    last->setBitQuick(i, false);
                }
            }
            //set F
            //set W
            uint64_t curW = getW(kmers[i]);
            Wvec[i] = curW;
            if (!curW && i)
                p=i;
            if (i) {
                //uint64_t cF=get_alphabet_number(getPos(kmers[i], k-1, alphabet, alph_size));
                //#pragma omp critical
                //F[cF] = std::min(F[cF],i-1);
                for (int j=i-1;j>=0 && compare_kmer_suffix(kmers[j], kmers[i], 1);--j) {
                    if ((Wvec[j] % alph_size) == curW) {
                        Wvec[i] += alph_size;
                        break;
                    }
             
                }
            }
        }
    }
    F[0] = 0;
//    uint64_t cF = get_alphabet_number(getPos(kmers[0], k-1, alphabet, alph_size));
    for (size_t i=0;i<kmers.size();++i) {
        uint64_t cF=getPos(kmers[i], k-1, alphabet, alph_size);
        /*
        char *curseq = kmertos(kmers[i], alphabet, alph_size);
        std::cerr << "-" << (*last)[i] << " " << curseq << " " << get_alphabet_symbol(Wvec[i]) << (Wvec[i]>=alph_size ? "-" : "") << "\n"; 
        free(curseq);
        */
        if (cF != alphabet[lastlet]) {
            //std::cerr << cF << " " << lastlet << " " << i-1 << "\n";
            for (lastlet++;lastlet<alph_size;lastlet++) {
                F[lastlet]=i-1;
                if (alphabet[lastlet]==cF) {
                    break;
                }
            }
        }
    }
    std::cerr << (clock()-start)/CLOCKS_PER_SEC << "\n";
    start=clock();

    std::cerr << "Building wavelet tree\t";
    W = new WaveletTree(Wvec, 4, parallel);
    assert(W->size() == Wvec.size());
    assert((*W)[3] == Wvec[3]);
    Wvec.clear();
    kmers.clear();
    /*
    for (size_t i=0;i<W->size();++i) {
        char *curseq = kmertos(kmers[i], alphabet, alph_size);
        std::cout << i << " " << (*last)[i] << " " << curseq+1 << " " << curseq[0] << ((*Wvec)[i] >= alph_size ? "-":"") << " " << get_alphabet_symbol((*W)[i]) << "\n";
        free(curseq);
    }
    */
    std::cerr << (clock()-start)/CLOCKS_PER_SEC << "\n";
    start=clock();
    FILE* sfile = fopen("/proc/self/status","r");
    char line[128];
    while (fgets(line, 128, sfile) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            break;
        }
    }
    fclose(sfile);
    fprintf(stdout, "edges %lu / nodes %lu/ %s\n", get_edge_count(), get_node_count(), line);
    //std::cerr << kmers.size() << " " << W->n << " " << last->size() << "\n";
    //assert(W->n == last->size());
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

    if (config->verbose)
        std::cout << "    adding " << g->get_size() << " edges" << std::endl;
    // handle last and W
    for (size_t j = 1; j < g->get_size(); ++j) {
        this->last->insertBit(curr_pos, g->get_last(j));
        this->W->insert(g->get_W(j), curr_pos);
        ++curr_pos;
    }
    if (config->verbose)
        std::cout << "new total edges: " << g->W->n << std::endl;

    // handle F
    assert(this->F.size() == g->F.size());
    for (size_t j = 0; j < this->F.size(); ++j) {
        this->F.at(j) += g->F.at(j);
    }
}

/** 
 * This function takes a pointer to a graph structure and concatenates the arrays W, last 
 * and F to this graph's static containers last_stat and W_stat. In almost all cases 
 * this will not produce a valid graph and should only be used as a helper in the 
 * parallel merge procedure.
 */
void DBG_succ::append_graph_static(DBG_succ *g) {

    size_t n = g->get_size();
    if (config->verbose)
        std::cout << "    adding " << n << " edges" << std::endl;

    //size_t n_old = this->last_stat.size();
    //this->last_stat.resize(n_old + n);
    //this->W_stat.resize(n_old + n);

    size_t const b = 4;
    std::vector<uint64_t> offsets (1ull << b, 0);
    std::queue<uint64_t> blocks;
    std::queue<uint64_t> new_blocks;
    blocks.push(n);
    size_t pos = 0;
    uint64_t o = 0;
    for (size_t ib = 0; ib < b; ++ib) {
        while (!blocks.empty()) {
            uint64_t cnt = blocks.front();
            blocks.pop();
            uint64_t epos = pos + cnt;
            for ( ; pos < epos; ++pos) {
                offsets.at(o) += !(*(g->W->R))[pos];
            }
            if (ib < b - 1) {
                new_blocks.push(offsets.at(o));
                new_blocks.push(cnt - offsets.at(o));
            }
            o++; 
        }
        if (ib < b - 1)
            blocks.swap(new_blocks);
    }

    //std::cerr << "R size: " << g->W->R->size() << std::endl;

    bool bit;
    std::vector<uint64_t> upto_offsets ((1ull << (b - 1)) - 1, 0);
    uint64_t p, co, v, m;
    for (size_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        //v = (uint64_t) W_stat.at(i);
        v = 0;
        o = 0;
        p = i;
        co = 0;
        for (size_t ib = 0; ib < b - 1; ++ib) {
            bit = (*(g->W->R))[ib * n + p + co];
            if (bit) {
                v |= m;
                co += offsets.at(o);
                p -= upto_offsets.at(o);
            } else {
                p -= (p - upto_offsets.at(o)); 
                upto_offsets.at(o) += 1;
            }
            o = 2*o + 1 + bit;
            m >>= 1;
        }
        bit = (*(g->W->R))[(b - 1) * n + p + co];
        if (bit) {
            v |= m;
        }
        if (i == 0)
            continue;
        this->W_stat.push_back(v);
        this->last_stat.push_back(g->get_last(i));
    }

    // handle F
    assert(this->F.size() == g->F.size());
    for (size_t j = 0; j < this->F.size(); ++j) {
        this->F.at(j) += g->F.at(j);
    }
}

void DBG_succ::toDynamic() {

    size_t const b = 4;
    size_t const n = W_stat.size();

    // compute total offsets for the individual bins
    std::vector<uint64_t> offsets ((1ull << (b - 1)) - 1, 0);
    uint64_t v, m, o, p;
    for (size_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = (uint64_t) W_stat.at(i);
        o = 0;
        for (size_t ib = 1; ib < b; ++ib) {
            bool const bit  = m & v;
            if (!bit)
                offsets.at(o) += 1;
            o = 2*o + 1 + bit;
            m >>= 1;
        }
    }

    libmaus2::bitbtree::BitBTree<6, 64> *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);  

    uint64_t co;
    bool bit;
    std::vector<uint64_t> upto_offsets ((1ull << (b - 1)) - 1, 0);
    for (size_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = (uint64_t) W_stat.at(i);
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
    W_stat.clear();
    delete W;
    //W = new libmaus2::wavelet::DynamicWaveletTree<6, 64> (tmp, b, n);
    W = new WaveletTree(tmp, b, n);

    libmaus2::bitbtree::BitBTree<6, 64> *last_new = new libmaus2::bitbtree::BitBTree<6, 64>(last_stat.size(), false);
    for (size_t i = 0; i < last_stat.size(); ++i)
        if (last_stat.at(i))
            last_new->setBitQuick(i, true);
    last_stat.clear();
    delete last;
    last = last_new;
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

void DBG_succ::annotate_kmer(std::string &kmer, uint32_t &label_id, uint64_t &idx, pthread_mutex_t* anno_mutex, bool ignore) {

    // we just need to walk one step in the path
    if (idx > 0) {
        TAlphabet s = get_alphabet_number(kmer[kmer.length() - 1]);
        uint64_t rl = std::min(succ_W(pred_last(idx - 1) + 1, s), succ_W(pred_last(idx - 1) + 1, s + alph_size));
        uint64_t ru = std::max(pred_W(idx, s), pred_W(idx, s + alph_size));
        rl = outgoing(rl, s);
        ru = outgoing(ru, s);
        idx = (ru > rl) ? ru : rl;
    // we need to look up the full kmer
    } else {
        idx = this->index(kmer);
    }
    //std::cerr << "kmer: " << kmer << " idx: " << idx << " ignore: " << (ignore?"yes":"no") << " label_id: " << label_id << std::endl;
    //assert(idx > 0);
    if ((idx == 0) || ignore)
        return;

    std::vector<uint32_t> curr_combo;

    if (anno_mutex)
        pthread_mutex_lock(anno_mutex);

    // get annotation of current kmer
    size_t curr_hash = this->annotation[idx];
    if (curr_hash > 0)
        curr_combo = get_curr_combination(combination_vector, annotation_map[curr_hash]);
   
    //if (anno_mutex)
    //    pthread_mutex_unlock(anno_mutex);

    // check whether current label is already part of current annotation
    // and add it if not
    if (!std::binary_search(curr_combo.begin(), curr_combo.end(), label_id)) {

        curr_combo = add_to_combination(curr_combo, label_id); 
        curr_hash = AnnotationHash{}(curr_combo);

       // if (anno_mutex)
       //     pthread_mutex_lock(anno_mutex);

        annotation[idx] = curr_hash;
        //annotation_map[curr_hash] = curr_anno;
        if (annotation_map.find(curr_hash) == annotation_map.end()) {
            annotation_map[curr_hash] = insert_new_combination(combination_vector, curr_combo);
            combination_count++;
        }
    }
    if (anno_mutex)
        pthread_mutex_unlock(anno_mutex);
    //}
}

void DBG_succ::annotate_seq(kstring_t &seq, kstring_t &label, uint64_t start, uint64_t end, pthread_mutex_t* anno_mutex) {

    std::string curr_kmer;
    std::string label_str = std::string(label.s);
    uint32_t label_id;

    if (anno_mutex)
        pthread_mutex_lock(anno_mutex);

     // does the current label already have an ID?
    std::unordered_map<std::string, uint32_t>::iterator id_it = label_to_id_map.find(label_str);
    if (id_it == label_to_id_map.end()) {
        label_id = (uint32_t) id_to_label.size();
        id_to_label.push_back(label_str);
        label_to_id_map[label_str] = label_id;
        if (config->verbose)
            std::cout << "added label ID " << label_id << " for label string " << label_str << std::endl;
    } else { 
        label_id = id_it->second;
    }

    if (anno_mutex)
        pthread_mutex_unlock(anno_mutex);

    end = (end == 0) ? seq.l : end;

    uint64_t previous_idx = 0;
    size_t i;
    for (i = start; i < end; ++i) {

        if (config->verbose && i > 0 && i % 1000 == 0) {
            std::cout << "." << std::flush;
            if (!anno_mutex && (i % 10000 == 0))
                std::cout << i << " kmers added" << std::endl;
        }

        if (curr_kmer.size() < k) {
            curr_kmer.push_back(seq.s[i]);
            continue;
        }
        assert(curr_kmer.size() == k);
        annotate_kmer(curr_kmer, label_id, previous_idx, anno_mutex, (i % config->frequency) > 0);
        
        //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
        curr_kmer.push_back(seq.s[i]);
        curr_kmer = curr_kmer.substr(1, k);
    }
    // add last kmer and label to database
    if (curr_kmer.size() == k)
        annotate_kmer(curr_kmer, label_id, previous_idx, anno_mutex, (i % config->frequency) > 0);
    //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
}   


std::vector<uint32_t> DBG_succ::classify_path(std::vector<uint64_t> path) {
    
    uint32_t curr_anno;
    std::vector<uint32_t> labels;
    std::vector<uint32_t> current_combination;
    std::map<uint32_t, uint64_t> label_counter;

    // collect all annotated combinations for the path
    // take majority vote as consensus for now
    for (size_t i = 0; i < path.size(); i++) {
        curr_anno = annotation.at(path.at(i));
        if (curr_anno > 0) {
            current_combination = get_curr_combination(combination_vector, annotation_map[curr_anno]);
            for (std::vector<uint32_t>::iterator c = current_combination.begin(); c != current_combination.end(); c++) {
                if (label_counter.find(*c) != label_counter.end()) {
                    label_counter[*c] += 1;
                } else {
                    label_counter[*c] = 1;
                }
            }
        }
    }

    // take majority vote as consensus for now
    if (label_counter.size() == 1) {
        labels.push_back(label_counter.begin()->first);
    } else if (label_counter.size() > 0) {
        uint32_t curr_max = 0;
        for (std::map<uint32_t, uint64_t>::iterator c = label_counter.begin(); c != label_counter.end(); c++) {
            if (c->second > curr_max) {
                labels.clear();
                curr_max = c->second;
            } 
            if (c->second == curr_max) {
                labels.push_back(c->first);
            }
        }
    }

    return labels;
}


std::set<uint32_t> DBG_succ::classify_read(kstring_t &read, uint64_t max_distance) {

    // containers for label information
    std::vector<uint32_t> path_labels;
    std::set<uint32_t> all_labels;

    // get alignment of the read to the graph
    std::string read_str = std::string(read.s);
    std::vector<HitInfo> alignment = index_fuzzy(read_str, max_distance);

    // classify hits
    for (size_t i = 0; i < alignment.size(); i++) {
        path_labels = classify_path(alignment.at(i).path);
        for (size_t j = 0; j < path_labels.size(); j++) {
            all_labels.insert(path_labels.at(j));
        }
    }

    return all_labels;
}


//
//
// MERGE
//
//


void DBG_succ::split_range(std::deque<TAlphabet>* str, std::pair<uint64_t, uint64_t> &range) {
    range = index_range(*str);
    if (range.first > range.second) {
        range.first = 0;
        range.second = 0;
    } else if (range.first > 0) {
        range.first = pred_last(range.first - 1) + 1;
    }
}

/*
 * Helper function that will split up a given range in the graph
 * into bins, one for each character in the alphabet. The split is performed based
 * on the k - d column of the node label. It is assumed that the all nodes in the
 * given range share a common suffix of length d.
 */
std::vector<uint64_t> DBG_succ::split_range(uint64_t start, uint64_t end, uint64_t d /*depth*/) {

    std::vector<uint64_t> result;
    uint64_t _end = end;

    // special case d == 0
    if (d == 0) {
        for (size_t i = 1; i < F.size(); ++i) {
            result.push_back((F.at(i) != F.at(i-1))?F.at(i-1)+1:0);  
        }
        result.push_back(F.back() + 1);
    // walk d-1 steps backwards
    } else{
        std::stack<uint64_t> backtrace;
        end--;
        while (d > 0) {
            std::cerr << "BT push " << get_node_end_value(start) << " start " << start << " --> " << bwd(succ_last(start)) << " end " << end << " --> " << bwd(succ_last(end)) << std::endl;
            backtrace.push(get_node_end_value(start));
            start = bwd(succ_last(start));
            end = bwd(succ_last(end));
            d--;
        }
        uint64_t c;
        uint64_t cc = backtrace.top();
        uint64_t v = start;
        for (size_t i = 0; i < alph_size; ++i) {
            c = succ_W(pred_last(v - 1) + 1, cc);
            //std::cerr << "v: " << v << " cc: " << cc << " c: " << c << std::endl;
            if (i < (alph_size - 1)) {
                std::cerr << "v is: " << v << " push back to result: " << (c < F.at(i + 1) ? c : 0) << std::endl;
                result.push_back((c < F.at(i + 1)) ? c : 0);
                v = std::max(v, F.at(i + 1) + 1);
            } else {
                result.push_back(c < W->n ? c : 0);
            }
        }
        // translate into final coordinates
        while (backtrace.size() > 0) {
            for (size_t i = 0; i < result.size(); ++i) {
                //std::cerr << "res[i] " << result.at(i) << " pred_last(res[i] - 1) " << pred_last(result.at(i) - 1) << " bt " << backtrace.top() << std::endl;
                if (result.at(i) > 0) {
                    uint64_t tmp = result.at(i);
                    uint64_t nextc = succ_W(pred_last(result.at(i) - 1) + 1, backtrace.top());
                    uint64_t maxc = succ_last(result.at(i));
                    std::cerr << "nextc " << nextc << " maxc " << maxc << std::endl;
                    result.at(i) = (nextc <= maxc) ? (pred_last(fwd(nextc) - 1) + 1) : 0;
                    std::cerr << "tracing back on " << backtrace.top() << " index " << i << " position " << tmp << " resulting in " << result.at(i) << std::endl;
                }
            }
            backtrace.pop();
        }
        // shift to beginning of equal node range
        for (size_t i = 0; i < result.size(); ++i) {
            if (result.at(i) > 0)
                result.at(i) = pred_last(result.at(i) - 1) + 1;
        }
    }
    result.push_back(_end);

    return result;

}

/* 
 * Helper function to determine the bin boundaries, given 
 * a number of threads.
 */
std::vector<std::pair<uint64_t, uint64_t> > DBG_succ::get_bins(uint64_t threads, uint64_t bins_per_thread, DBG_succ* G) {

    uint64_t binlen = 0;
    uint64_t bins = threads * bins_per_thread * 10;

    // depending on the number of threads, we will use a different length prefix
    // to compute the bin boundaries
    /*if (bins < alph_size)
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
        binlen = std::min(6lu, k);
    */
    uint64_t exp_binlen = ((uint64_t) std::log((double) bins) / std::log((double) alph_size)) + 1;
    binlen = std::min(k, exp_binlen);
    std::cerr << "target binlen is " << binlen << ", generating " << std::pow(binlen, alph_size) << " possible bins" << std::endl;

    std::vector<std::pair<uint64_t, uint64_t> > tmp;
    for (uint64_t i = 0; i < std::pow(alph_size, binlen); ++i) {
        std::pair<uint64_t, uint64_t> idx = G->index_range(bin_id_to_string(i, binlen));
        if (idx.first > idx.second)
            idx.first = idx.second = 0;
        if (idx.second >= G->get_size())
            idx.second = G->get_size() - 1;
        if (idx.first > 0)
            idx.first = G->pred_last(idx.first - 1) + 1;
        tmp.push_back(idx);
        /*std::deque<TAlphabet> ttt = bin_id_to_string(i, binlen);
        std::cerr << std::endl << "id: " << i << " ";
        for (size_t ii = 0; ii < ttt.size(); ++ii) {
            //std::cerr << ttt[ii] << get_alphabet_symbol(ttt[ii] % alph_size);
            std::cerr << get_alphabet_symbol(ttt[ii] % alph_size);
        }
        std::cerr << " - " << idx.first << ":" << idx.second;
        */
    }
    return tmp;
}

std::vector<std::pair<uint64_t, uint64_t> > DBG_succ::get_bins(uint64_t bins) {

    uint64_t nodes = this->rank_last(this->get_size() - 1);
    uint64_t orig_bins = bins;
    std::cerr << "working with " << orig_bins << " orig bins; " << nodes << " nodes" <<  std::endl;
    if (bins > nodes) {
        std::cerr << "[WARNING] There are max " << nodes << " slots available for binning. Your current choice is " << bins << " which will create " << bins - nodes << " empty slots." << std::endl;
        bins = nodes;
    }

    std::vector<std::pair<uint64_t, uint64_t> > result;
    uint64_t binsize = (nodes + bins - 1) / bins;
    uint64_t thresh = (nodes - (bins * (nodes / bins))) * binsize;
    uint64_t pos = 1;
    for (uint64_t i = 0; i < nodes;) {
        if (i >= thresh) {
            binsize = nodes / bins;
        }
        //std::cerr << "push " << pos << " - " << this->select_last(std::min(nodes, i + binsize)) << std::endl;
        result.push_back(std::make_pair(pos, this->select_last(std::min(nodes, i + binsize))));
        pos = result.back().second + 1;
        i += binsize;
    }

    for (uint64_t i = bins; i < orig_bins; i++) {
        //result.push_back(std::make_pair(pos, pos));
        result.push_back(std::make_pair(1, 0));
    }
    
    std::cerr << "created " << result.size() << " bins" << std::endl;
    return result;
}

std::vector<std::pair<uint64_t, uint64_t> > DBG_succ::get_bins_relative(DBG_succ* G, std::vector<std::pair<uint64_t, uint64_t> > ref_bins, uint64_t first_pos, uint64_t last_pos) {
    
    std::vector<std::pair<uint64_t, uint64_t> > result;
    uint64_t pos = (first_pos == 0) ? 1 : this->index_predecessor(G->get_node_seq(first_pos)) + 1;
    uint64_t upper;
    for (size_t i = 0; i < ref_bins.size(); i++) {
        if (ref_bins.at(i).second == 0) { // this happens if we have more bins than nodes
            result.push_back(std::make_pair(0, 0));
        } else {
            upper = this->index_predecessor(G->get_node_seq(ref_bins.at(i).second));
            std::cerr << "ref bin " << ref_bins.at(i).second << " rel upper " << upper << std::endl;
            result.push_back(std::make_pair(pos, upper));
            pos = upper + 1;
        }
    }
    result.back().second = (last_pos == 0) ? this->get_size() - 1 : result.back().second;
    return result;
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

uint64_t DBG_succ::next_non_zero(std::vector<std::pair<uint64_t, std::deque<TAlphabet> > > v, uint64_t pos) {

    uint64_t val = (pos < v.size()) ? v.at(pos).first : 0;

    while ((pos < v.size() - 1) && (v.at(pos).first == 0)) {
        val = v.at(pos + 1).first;
        pos++;
    }
    return val;
}


uint64_t DBG_succ::next_non_zero(std::vector<uint64_t> v, uint64_t pos) {

    uint64_t val = (pos < v.size()) ? v.at(pos) : 0;

    while ((pos < v.size() - 1) && (v.at(pos) == 0)) {
        val = v.at(pos + 1);
        pos++;
    }
    return val;
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

void DBG_succ::print_adj_list() {
    std::pair<uint64_t, uint64_t> R;
    uint64_t i = 1;
    R.first = i;
    R.second = succ_last(R.first);
    uint64_t n = rank_last(R.second);
    while (R.first < W->n) {
        printf("%lu\t", n);
        for (uint64_t j = R.first; j<=R.second; ++j) {
            if (j > R.first)
                fprintf(stdout, ",");
            fprintf(stdout, "%lu", rank_last(fwd(j)));
        }
        if (this->annotation.size() > 0)
            fprintf(stdout, "\t%u", this->annotation.at(n));
        printf("\n");
        R.first = R.second+1;
        if (R.first >= W->n)
            break;
        R.second = succ_last(R.first);
        n++;
    }
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
void DBG_succ::toFile(unsigned int total, unsigned int idx) {

    // adapt outbase depending on merging strategy
    std::string outbase = config->outfbase;
    if (total > 1)
        outbase += "." + std::to_string(idx) + "_" + std::to_string(total);

    // write Wavelet Tree
    std::ofstream outstream((outbase + ".W.dbg").c_str());
    W->serialise(outstream);
    outstream.close();

    // write last array
    outstream.open((outbase + ".l.dbg").c_str());
    last->serialise(outstream);
    outstream.close();

    // write F values and k
    outstream.open((outbase + ".F.dbg").c_str());
    outstream << ">F" << std::endl;
    for (size_t i = 0; i < F.size(); ++i)
        outstream << F.at(i) << std::endl;
    outstream << ">k" << std::endl;
    outstream << k << std::endl;
    outstream << ">p" << std::endl;
    outstream << p << std::endl;
    outstream.close();
}


// write annotation to screen
void DBG_succ::annotationToScreen() {
    std::deque<uint32_t>::iterator ait = annotation.begin();
    //std::set<uint32_t>::iterator sit;
    std::vector<uint32_t>::iterator vit;
    std::vector<uint32_t> curr_comb;
    for (; ait != annotation.end(); ++ait) {
        std::cout << *ait << " : " << annotation_map[*ait] << " : ";
        if (*ait > 0) {
            curr_comb = get_curr_combination(combination_vector, annotation_map[*ait]);
            //sit = annotation_map[*ait].begin();
            for (vit = curr_comb.begin(); vit != curr_comb.end(); ++vit) {
                std::cout << id_to_label.at(*vit) << " ";
            }
        }
        std::cout << std::endl;
    }
    //std::cerr << "ANNO MAP " << std::endl;
    //for (std::unordered_map<uint32_t, uint32_t>::iterator itt = annotation_map.begin(); itt != annotation_map.end(); itt++)
    //    std::cerr << itt->first << ":" << itt->second << "; ";
    //std::cerr << std::endl;
}

// write annotation to disk
void DBG_succ::annotationToFile() {
    std::ofstream outstream((config->infbase + ".anno.dbg").c_str());
    //annotation
    libmaus2::util::NumberSerialisation::serialiseNumber32Deque<uint32_t>(outstream, annotation);
    // annotation_map
    serialize_annotation_map(outstream, annotation_map);
    //combination vector
    serialize_combination_vector(outstream, combination_vector);
    // id_to_label
    serialize_annotation_id_vector(outstream, id_to_label);
    // label_to_id_map
    serialize_label_to_id_map(outstream, label_to_id_map);
    // combination_count
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, combination_count);
    outstream.close();
}

// read annotation from disk
void DBG_succ::annotationFromFile() {
    std::ifstream instream((config->infbase + ".anno.dbg").c_str());
    if (instream.good()) {
        if (config->verbose)
            std::cerr << "get deque from disk" << std::endl;
        // annotation
        annotation = libmaus2::util::NumberSerialisation::deserialiseNumber32Deque<uint32_t>(instream);
        if (config->verbose)
            std::cerr << "get map from disk" << std::endl;
        // annotation_map
        deserialize_annotation_map(instream, annotation_map);
        // combination_vector
        combination_vector = deserialize_combination_vector(instream);
        // id_to_label
        id_to_label = deserialize_annotation_id_vector(instream);
        // label_to_id_map
        deserialize_label_to_id_map(instream, label_to_id_map);
        combination_count = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
    } else {
        annotation.resize(get_size(), 0);
    }
    instream.close();
}

