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

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include <libmaus2/util/NumberSerialisation.hpp>

#include "config.hpp"
#include "datatypes.hpp"
#include "serialization.hpp"
#include "dbg_succinct_libmaus.hpp"

// define an extended alphabet for W --> somehow this does not work properly as expected
typedef uint64_t TAlphabet;

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
    config(config_),
    alphabet("$ACGTNX$ACGTNXn") {

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
    label_to_id_map[""]=0;
    combination_vector.push_back(0);
    state = Config::dyn;
}



DBG_succ::DBG_succ(std::string infbase_, Config* config_) : 
    infbase(infbase_),
    config(config_),
    alphabet("$ACGTNX$ACGTNXn") {

    std::ifstream instream;
    // if not specified in the file, the default for loading is dynamic
    state = Config::dyn;

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
        } else if (strcmp(line.c_str(), ">s") == 0) {
            mode = 4;
        } else {
            if (mode == 1) {
                F.at(fidx) += std::strtoul(line.c_str(), NULL, 10);
                fidx++;
            } else if (mode == 2) {
                k = strtoul(line.c_str(), NULL, 10);
            } else if (mode == 3) {
                p = strtoul(line.c_str(), NULL, 10);
            } else if (mode == 4) {
                state = (Config::state_type) strtoul(line.c_str(), NULL, 10);
            } else {
                fprintf(stderr, "ERROR: input file corrupted\n");
                exit(1);
            }
        }
    }
    instream.close();
    id_to_label.push_back("");
    label_to_id_map[""]=0;
    combination_vector.push_back(0);

    // load last array
    //std::ifstream instream((infbase + ".l.dbg").c_str());
    //last->deserialise(instream);
    //instream.close();

    // load W and last arrays
    delete W;
    delete last;
    std::ifstream instream_W((infbase + ".W.dbg").c_str());
    std::ifstream instream_l((infbase + ".l.dbg").c_str());
    switch (state) {
        case Config::dyn: {
            W = new wavelet_tree_dyn(instream_W);
            last = new bit_vector_dyn(instream_l);
        } break;
        case Config::stat: {
            W = new wavelet_tree_stat(instream_W);
            last = new bit_vector_stat(instream_l);
        } break;
    }
    instream_W.close();
    instream_l.close();
}

DBG_succ::~DBG_succ() {

    delete W;
    delete last;
    if (bridge != NULL)
        delete bridge;
    for (auto it = annotation_full.begin(); it != annotation_full.end(); ++it) {
        if (*it != NULL) {
            delete *it;
        }
    }
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
    return W->rank(c, std::min(i, W->size() - 1)) - (c == 0);
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
    return std::min(W->select(c, i - 1 + (c == 0)), W->size());
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
    //usleep(500000);
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
 * index of the outgoing edge with label c.
 */
uint64_t DBG_succ::outgoing_edge_idx(uint64_t i, TAlphabet c) {

    if (i > W->size())
        return 0;
    if (i == 0)
        return 0;
    std::pair<uint64_t, uint64_t> R = get_equal_node_range(i);

    uint64_t j1 = pred_W(R.second, c);
    uint64_t j2 = pred_W(R.second, c + alph_size);
    uint64_t j = (j1 < j2) ? j2 : j1;
    if (j < R.first || j >= W->size())
        return 0;

    return j;
}



/**
 * Given a position i in W and an edge label c, this function returns the
 * index of the node the edge is pointing to.
 */
uint64_t DBG_succ::outgoing(uint64_t i, TAlphabet c) {
    if (i > W->size())
        return 0;
    if (i == 0)
        return 0;

    uint64_t j = outgoing_edge_idx(i, c);
    if (j == 0)
        return 0;
    j = fwd(j);
    if (j == 0 || j == W->size())
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
    return (i < W->size()) ? succ_last(i) - pred_last(i - 1) : 0;
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
        if (max_distance > 0) {
            for (size_t p = 1; p < str.length() - 1; ++p) {
                TAlphabet ss = get_alphabet_number(str[p]);
                if ((p + (b != ss)) > max_distance)
                    break;
                hits.push(HitInfo(rl, ru, p + 1, 1, p + (b != ss), std::string(p, 'd') + std::string(1, get_alphabet_symbol(b)), tmp));
                //std::cout << "a) adding '-'" << std::endl;
            }
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
                        if ((rl >= W->size()) || (ru >= W->size()) || (rl > ru))
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
    uint64_t ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->size() - 1);    // upper bound
    while (rl > ru && s1 > 0) {
        s1--;
        rl = succ_last(F.at(s1) + 1);
        ru = (s1 < F.size() - 1) ? F.at(s1 + 1) : (W->size() - 1);
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
                if (rl >= W->size()) {
                    s1--;
                    while (s1 > 0) {
                        rl = std::max(pred_W(W->size() - 1, s1), pred_W(W->size() - 1, s1));
                        if (rl < W->size())
                            break;
                        s1--;
                    }
                    if (s1 == 0) {
                        s1 = (*it % alph_size) + 1;
                        before = true;
                        while (s1 < alph_size) {
                            rl = std::min(succ_W(1, s1), succ_W(1 + alph_size, s1));
                            if (rl < W->size())
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

bool DBG_succ::compare_node_suffix(TAlphabet *ref, uint64_t i2) {
    TAlphabet *i1 = &ref[k-1];
    for (size_t ii=0; ii < k-1;ii++) {
        if (*i1 != get_node_end_value(i2)) {
            return false;
        }
        i1 = &ref[k-2-ii];
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
        ret << get_alphabet_symbol(k_val.first);
        k_node = k_val.second;
    }
    std::string ret_str = ret.str();
    return std::string(ret_str.rbegin(), ret_str.rend());
}

/**
 * Return number of edges in the current graph.
 */
uint64_t DBG_succ::get_size() {
    return W->size();
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
    return s < alphabet.size() ? alphabet[s] : alphabet[14];
}

std::vector<uint64_t> DBG_succ::align(kstring_t seq, uint64_t alignment_length) {

    size_t seq_length = seq.l;
    //std::vector<std::vector<HitInfo> > hit_list;
    std::vector<uint64_t> indices;

    if (alignment_length == 0)
        alignment_length = this->get_k();

  std::vector<HitInfo> curr_result;
  for (uint64_t i = 0; i < seq_length - alignment_length + 1; ++i) {
    std::string kmer(seq.s + i, seq.s + i + alignment_length);
    indices.push_back(this->index(kmer, kmer.size()));
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
 //TODO: Check this!!!
uint64_t DBG_succ::get_edge_count() {
    return W->size() - 1;
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


void DBG_succ::switch_state(Config::state_type state) {

    //std::cerr << "switching state from " << this->state << " to " << state << std::endl;
    if (this->state == state)
        return;
    
    switch (state) {
        case Config::cstr: {
            this->state = Config::cstr;
        } break;

        case Config::dyn: {
            if (this->state == Config::cstr) {
                delete W;
                W = new wavelet_tree_dyn(W_stat, 4);
                W_stat.clear();

                delete last;
                if (last_stat.size()) {
                    last = new bit_vector_dyn(last_stat);
                    last_stat.clear();
                } else {
                    last = new bit_vector_dyn(last_stat_safe); 
                    last_stat_safe.clear();
                }

                delete bridge;
                if (bridge_stat.size()) {
                    bridge = new bit_vector_dyn(bridge_stat);
                    bridge_stat.clear();
                } else {
                    bridge = NULL;
                }
            } else {
                wavelet_tree* W_new = new wavelet_tree_dyn(W, 4);
                delete W;
                W = W_new;

                bit_vector* last_new = new bit_vector_dyn(last);
                delete last;
                last = last_new;
            }
            this->state = Config::dyn;
        } break;

        case Config::stat: {
            wavelet_tree* W_new = new wavelet_tree_stat(W, 4);
            delete W;
            W = W_new;

            bit_vector* last_new = new bit_vector_stat(last);
            delete last;
            last = last_new;

            this->state = Config::stat;
        } break;
    }
}



//
//
// MERGE
//
//


void DBG_succ::split_range(std::deque<TAlphabet>* str, std::pair<uint64_t, uint64_t> &range) {
    range = index_range(*str, str->size());
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
                result.push_back(c < W->size() ? c : 0);
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
    for (uint64_t i = 0; i < W->size(); i++) {
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

//TODO: assume that a sentinel was used during construction
void DBG_succ::print_state_str() {
    for (uint64_t i = 1; i < W->size(); i++) {
        std::cout << i << "\t" << get_last(i) << "\t" << get_node_str(i) << "\t" << get_alphabet_symbol((*W)[i]) << ((*W)[i] > alph_size ? "-" : "") << (i==this->p ? "<" : "") << std::endl;
    }    
}


void DBG_succ::print_adj_list() {
    for (uint64_t edge = 1; edge < W->size(); ++edge) {
            fprintf(stdout, "%lu\t%lu\t", rank_last(succ_last(edge)), rank_last(outgoing(edge, (*W)[edge])));
            bool is_first = true;
            for (uint64_t k = 0; k < this->annotation_full.size(); ++k) {
                if ((*(this->annotation_full.at(k)))[edge] == 1) {
                    if (!is_first)
                        fprintf(stdout, ",");
                    is_first = false;
                    fprintf(stdout, "%lu", k+1);
                }
            }
            if (is_first)
                fprintf(stdout, "0");
            printf("\n");
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
    uint64_t end = start + linelen < W->size() ? start + linelen : W->size();

    while (start < W->size()) {
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
        end = start + linelen < W->size() ? start + linelen : W->size();
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
    outstream << ">s" << std::endl;
    outstream << state << std::endl;
    outstream.close();
}


// write annotation to disk
void DBG_succ::annotationToFile() {
    std::ofstream outstream((config->infbase + ".anno.dbg").c_str());
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, annotation_full.size());
    for (size_t i = 0; i < annotation_full.size(); ++i) {
        annotation_full.at(i)->serialize(outstream);
    }

    // id_to_label
    serialize_annotation_id_vector(outstream, id_to_label);
    // label_to_id_map
    serialize_label_to_id_map(outstream, label_to_id_map);
/*
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
*/
    outstream.close();
}

// read annotation from disk
void DBG_succ::annotationFromFile() {
    // generate annotation object
    // populate it with existing annotation if available
    std::ifstream instream((config->infbase + ".anno.dbg").c_str());
    if (instream.good()) {
        if (config->verbose)
            std::cerr << "get annotation from disk" << std::endl;
        size_t anno_size = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
        for (size_t i = 0; i < anno_size; ++i) {
            //annotation_full.push_back(new sdsl::rrr_vector<63>());
            annotation_full.push_back(new sdsl::sd_vector<>());
            annotation_full.back()->load(instream);
        }

        // id_to_label
        id_to_label = deserialize_annotation_id_vector(instream);
        // label_to_id_map
        deserialize_label_to_id_map(instream, label_to_id_map);

        //annotation_full = new sdsl::bit_vector(get_size() * 100, 0);

        /*
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
        */
   // } else {
        //annotation.resize(get_size(), 0);
    }
    instream.close();
}

