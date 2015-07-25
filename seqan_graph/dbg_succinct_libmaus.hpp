#ifndef __DBG_SUCCINCT_LIBM_HPP__
#define __DBG_SUCCINCT_LIBM_HPP__

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
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>

#include <libmaus/bitbtree/bitbtree.hpp>
#include <libmaus/wavelet/DynamicWaveletTree.hpp>

/** 
 * We use seqan for an efficient representation of our alphabet.
 */
#include <seqan/sequence.h>
#include <seqan/basic.h>

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'X'> > Dna5F; 
    typedef uint64_t TAlphabet;

    private:
        // the bit array indicating the last outgoing edge of a node
        libmaus::bitbtree::BitBTree<6, 64> *last = new libmaus::bitbtree::BitBTree<6, 64>();

        // the array containing the edge labels
        libmaus::wavelet::DynamicWaveletTree<6, 64> *W = new libmaus::wavelet::DynamicWaveletTree<6, 64>(4); // 4 is log (sigma)

        // the offset array to mark the offsets for the last column in the implicit node list
        std::vector<TAlphabet> F; 

        // k-mer size
        size_t k;
        // index of position that marks end in graph
        uint64_t p;
        // number of edges in the graph
        uint64_t m;
        // alphabet size
        size_t alph_size = 7;

        uint64_t W_size = 0;

#ifdef DBGDEBUG
        bool debug = true;
#else
        bool debug = false;
#endif 
    
    public:
        DBG_succ(size_t k) : k(k) {

            last->insertBit(0, true);
            last->insertBit(0, false);

            W->insert(0, 0);
            W->insert(0, 0);
            W_size += 2;

            F.push_back(0);
            for (size_t j = 1; j < alph_size; j++)
                F.push_back(1);
            
            m = 1;
            p = 1;
        }

        void
        add_seq (
            String<Dna5F> seq
        ) {
            if (debug) {
                print_seq();
                print_state();
                std::cout << "======================================" << std::endl;
            }

            if (W_size == 2) {
                for (size_t j = 0; j < k - 1; j++) {
                    append_pos(6);
                    if (debug) {
                        print_seq();
                        print_state();
                        std::cout << "======================================" << std::endl;
                    }
                }
            }

            for (uint64_t i = 0; i < length(seq); ++i) {
                //if (i > 0 && i % 100000 == 0) {
                if (i > 0 && i % 1000 == 0) {
                    std::cout << "." << std::flush;
                    if (i % 10000 == 0) {
                        fprintf(stdout, "%lu - edges %lu / nodes %lu\n", i, W_size - 1, rank_last((last->size() - 1)));
                    }
                }
                //fprintf(stdout, "appending %i\n", (int) ordValue(seq[i]));
                //cerr << "seq[i] " << seq[i] << std::endl;
                //cerr << "seq[i] ord " << ordValue(seq[i]) + 1 << std::endl;
                append_pos((TAlphabet) ordValue(seq[i]) + 1);
                if (debug) {
                    print_seq();
                    print_state();
                    std::cout << "======================================" << std::endl;
                }
            }

            for (size_t j = 0; j < k - 1; j++) {
                append_pos(6);
                if (debug) {
                    print_seq();
                    print_state();
                    std::cout << "======================================" << std::endl;
                }
            }

            fprintf(stdout, "edges %lu / nodes %lu\n", W_size - 1, rank_last((last->size() - 1)));

            //toSQL();
            //String<Dna5F> test = "CCT";
            //fprintf(stdout, "\nindex of CCT: %i\n", (int) index(test));
        }

    private:
    
        /** 
         * Uses the object's array W, a given position i in W and a character c
         * from the alphabet and returns the number of occurences of c in W up to
         * position i.
         */
        uint64_t rank_W(uint64_t i, TAlphabet c) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            return W->rank(c, std::min(i, W_size - 1));
        }

        /**
         * Uses the array W and gets a count i and a character c from 
         * the alphabet and returns the positions of the i-th occurence of c 
         * in W.
         */
        uint64_t select_W(uint64_t i, TAlphabet c) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            //fprintf(stderr, "query select W -- c: %i i: %lu return: %lu \n", c, i-1+(c==0), W->select(c, i-1+(c==0)));
            return std::min(W->select(c, i - 1 + (c == 0)), W_size);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the last index of a character c preceding in W[1..i].
         */
        uint64_t pred_W(uint64_t i, TAlphabet c) {
            return select_W(rank_W(i, c), c);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the first index of a character c in W[i..N].
         */
        uint64_t succ_W(uint64_t i, TAlphabet c) {
            return select_W(rank_W(i - 1, c) + 1, c);
        }

        /** 
         * Uses the object's array last and a position and
         * returns the number of set bits up to that postion.
         */
        uint64_t rank_last(uint64_t i) {
            // deal with  border conditions
            if (i <= 0)
                return 0;
            //if (i > last->size() - 1) {
            //    fprintf(stderr, "i %lu size %lu\n", i, last->size());
            //    return last->size() - 1;
            //}
            return last->rank1(i);
        }

        /**
         * Uses the object's array last and a given position i and
         * returns the position of the i-th set bit in last[1..i].
         */
        uint64_t select_last(uint64_t i) {
            // deal with  border conditions
            if (i <= 0)
                return 0;
            //fprintf(stderr, "i %lu size %lu\n", i, last->size());
            //if (i >= last->size())
            //    return last->size();
            // for some reason the libmaus select is 0 based ...
            return std::min(last->select1(i - 1), last->size());
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the last set bit in last[1..i].
         */
        uint64_t pred_last(uint64_t i) {
            return select_last(rank_last(i));
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the first set bit in last[i..N].
         */
        uint64_t succ_last(uint64_t i) {
            return select_last(rank_last(i - 1) + 1);
        }

        /**
         * Given a position i in W and an edge label c, this function returns the
         * index of the node the edge is pointing to.
         */
        uint64_t outgoing(uint64_t i, TAlphabet c) {
            if (i >= W_size - 1)
                return 0;
            std::pair<uint64_t, uint64_t> R = get_equal_node_range(i);

            uint64_t j1 = pred_W(R.second, c);
            uint64_t j2 = pred_W(R.second, c + alph_size);
            uint64_t j = (j1 < j2) ? j2 : j1;
            if (j < R.first || j >= W_size - 1)
                return 0;
            j = fwd(j);
            if (j == 0 || j == W_size)
                return 0;
            //if (j > 36000)
            //    fprintf(stdout, "i %i, c %i, j1 %i, j2 %i size W %i\n", (int) i, (int) c, (int) j1, (int) j2, W_size);
            return j;
        }


        /**
         * Given a node index i and an edge label c, this function returns the
         * index of the node the incoming edge belongs to.
         */
        uint64_t incoming(uint64_t i, TAlphabet c) {
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
        uint64_t outdegree(uint64_t i) {
            return (i < W_size - 1) ? succ_last(i) - pred_last(i - 1) : 0;
        }


        /**
         * Given a node index i, this function returns the number of incoming
         * edges to node i.
         */
        uint64_t indegree(uint64_t i) {
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
        uint64_t index(String<Dna5F> &s) {
            // init range
            uint64_t rl = succ_last(F[(TAlphabet) ordValue(s[0]) + 1] + 1);
            uint64_t ru = F[(TAlphabet) ordValue(s[0]) + 2];                // upper bound
            // update range iteratively while scanning through s
            for (uint64_t i = 1; i < length(s); i++) {
                rl = outgoing(rl, (TAlphabet) ordValue(s[i]) + 1);
                ru = outgoing(ru, (TAlphabet) ordValue(s[i]) + 1);
               // fprintf(stdout, "char: %i rl: %i ru: %i\n", (int) ordValue(s[i]) + 1, (int) rl, (int) ru);
            }
            return (ru > rl) ? ru : rl;
        }


        /**
         * Using the offset structure F this function returns the value of the last 
         * position of node i.
         */
        TAlphabet get_node_end_value(uint64_t i) {
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
        TAlphabet get_node_begin_value(uint64_t i) {
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
         * This function gets a position i that reflects the i-th node and returns the
         * position in W that corresponds to the i-th node's last character. 
         */
        uint64_t bwd(uint64_t i) {
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
        uint64_t fwd(uint64_t i) {
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
         * Given a position i, this function returns the boundaries of the interval
         * of nodes identical to node i (ignoring the values in W).
         */
        std::pair<uint64_t, uint64_t> get_equal_node_range(uint64_t i) {
            return std::make_pair(pred_last(i - 1) + 1, succ_last(i));
        }

        /**
         * This is a debug function that prints the current state of the graph arrays to
         * the screen.
         */
        void print_state() {

            fprintf(stderr, "W:\n");
            for (uint64_t i = 0; i < W_size; i++) {
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

        /** This function takes a character c and appends it to the end of the graph sequence
         * given that the corresponding note is not part of the graph yet.
         */
        void append_pos(TAlphabet c) {

            // check that the last position of the graph is indeed a terminal
            assert((*W)[p] == 0);
            TAlphabet c_p = get_node_end_value(p);
            // get range of identical nodes (without W) pos current end position
            std::pair<uint64_t, uint64_t> R = get_equal_node_range(this->p);
            //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

            // get position of first occurence of c in W after p
            uint64_t next_c = succ_W(p, c);
            //fprintf(stdout, "p %i, c %i, next_c %i\n", (int) p, (int) c, (int) next_c);
            // check if c is part of range
            bool exist_c = (next_c <= R.second); // && W[next_c] != 0 && !is_terminal_node(next_c));
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
                // we need an addtional pred_W here as fwd is not defined if next_c is element of the minus set
                //uint64_t p_new = fwd(pred_W(next_c, c));
                uint64_t p_new = fwd(next_c);
                // remove old terminal symbol
                last->deleteBit(p);
                W->remove(p);
                W_size--;
                // adapt position if altered by previous deletion
                p_new -= (p < p_new);
                // insert new terminal symbol 
                // we have to insert 0 into last as the node already existed in the range 
                // and the terminal symbol is always first
                last->insertBit(p_new, false);
                W->insert(0, p_new);
                W_size++;
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
                    if (next_c < W_size) {
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
                    W_size++;
                } else {
                    uint64_t x = F[c] + 1;
                    uint64_t next_c = succ_W(p + 1, c);
                    bool minus = false;
                    if (next_c < W_size) {
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
                    W_size++;
                }
                m++;
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


        /** 
         * This function gets two node indices and returns if the
         * node labels share a k-1 suffix.
         */
        bool compare_node_suffix(uint64_t i1, uint64_t i2) {
            for (size_t ii = 0; ii < k-1; ii++) {
                //std::cout << "node1 - " << i1 << ": " << Dna5F(get_node_end_value(i1) % alph_size - 1) << std::endl; 
                //std::cout << "node2 - " << i2 << ": " << Dna5F(get_node_end_value(i2) % alph_size - 1) << std::endl; 
                if (get_node_end_value(i1) != get_node_end_value(i2)) {
                    return false;
                }
                //std::cout << "succ (i1): " << succ_last(i1) << " succ (i2): " << succ_last(i2) << std::endl;
                i1 = bwd(succ_last(i1));
                i2 = bwd(succ_last(i2));
            }
            return true;
        }

        bool is_terminal_node(uint64_t i) {
            for (size_t ii = 0; ii < k-1; ii++) {
                if (get_node_end_value(i) % alph_size != 0) {
                    return false;
                }
                i = bwd(i);
            }
            return true;
        }


        /**
         * This function gets a value of the alphabet c and updates the offset of 
         * all following values by +1 is positive is true and by -1 otherwise.
         */
        void update_F(TAlphabet c, bool positive) {
            for (TAlphabet i = c+1; i < F.size(); i++)
                F[i] += (2 * (TAlphabet) positive - 1);
        }


        /**
         * Given index i of a node and a value k, this function 
         * will return the k-th last character of node i.
         */
        TAlphabet get_minus_k_value(uint64_t i, uint64_t k) {
            for (; k > 0; --k)
                i = bwd(succ_last(i));
            return get_node_end_value(i);
        }


        /*
         * Returns the sequence stored in W and prints the node
         * information in an overview. 
         * Useful for debugging purposes.
         */
        void print_seq() {

            for (uint64_t i = 1; i < W_size; i++) {
                if ((*W)[i] % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    std::cout << Dna5F(((*W)[i] % alph_size) - 1);
            }
            std::cout << std::endl;

            for (uint64_t i = 1; i < W_size; i++) {
                if (p == i)
                    fprintf(stdout, "*");
                else
                    fprintf(stdout, " ");
            }
            std::cout << std::endl;

            size_t j;
            for (size_t l = 0; l < k; l++) {
                for (uint64_t i = 1; i < W_size; i++) {
                    j = get_minus_k_value(i, l);
                    if (j % alph_size == 0)
                        std::cout << "$";
                    else
                        std::cout << Dna5F((j % alph_size) - 1);
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (uint64_t i = 1; i < last->size(); i++) {
                fprintf(stdout, "%i", (int) (*last)[i]);
            }
            std::cout << std::endl;
            std::cout << std::endl;

            for (uint64_t i = 1; i < W_size; ++i) {
                std::cout << indegree(i);  
            }
            std::cout << std::endl;
            for (uint64_t i = 1; i < W_size; ++i) {
                std::cout << outdegree(i);  
            }
            std::cout << std::endl;
            std::cout << std::endl;
            /*for (TAlphabet c = 0; c <= 5; ++c) {
                if (c == 0)
                    std::cout << "$";
                else
                    std::cout << Dna5F(c - 1);
                for (uint64_t i = 1; i < W_size; ++i) {
                    std::cout << " " << incoming(i, c);  
                }
                std::cout << std::endl;
            }*/
            /*for (TAlphabet c = 1; c <= 6; ++c) {
                if (c == 0)
                    std::cout << "$";
                else
                    std::cout << Dna5F(c - 1);
                for (uint64_t i = 1; i < W_size; ++i) {
                    std::cout << " " << outgoing(i, c);  
                }
                std::cout << std::endl;
            }*/


        }


        /**
         * This function gets a local range in W from lower bound l
         * to upper bound u and swaps the inserted element to the
         * righ location.
         */
        void sort_W_locally(uint64_t l, uint64_t u) {
            for (uint64_t s = u; s > l; s--) {
                TAlphabet tmp;
                if (((*W)[s] % alph_size) < ((*W)[s-1] % alph_size)) {
                    tmp = (*W)[s-1];
                    replaceW(s-1, (*W)[s]);
                    replaceW(s, tmp);
                }
            }
        }


        /** 
         * This is a convenience function to replace the value at
         * position i in W with val.
         */
        void replaceW(size_t i, TAlphabet val) {
            W->remove(i);
            W->insert(val, i);
        }


        /**
         * This object collects information about branches during graph traversal, so 
         * we know where to jump back to when we reached a dead end.
         */
        struct BranchInfo {
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
        struct JoinInfo {
            uint64_t seqId1;
            uint64_t seqPos1;
            uint64_t seqId2;
            uint64_t seqPos2;

            JoinInfo(uint64_t seqId1_, uint64_t seqPos1_, uint64_t seqId2_, uint64_t seqPos2_):
                seqId1(seqId1_),
                seqPos1(seqPos1_),
                seqId2(seqId2_),
                seqPos2(seqPos2_) {}
        };


        /**
         * This is a convenience function that pops the last branch and updates the traversal state.
         */
        BranchInfo pop_branch(std::stack<BranchInfo> &branchnodes, uint64_t &seqPos, uint64_t &nodeId, uint64_t &lastEdge, bool &isFirst) {
            BranchInfo branch = branchnodes.top();
            branchnodes.pop();
            isFirst = true;
            seqPos = 0;
            lastEdge = branch.lastEdge;
            nodeId = branch.nodeId;

            return branch;
        }

        /*std::ofstream get_fasta_stream(uint64_t seqId) {
            std::ofstream stream(sprintf("sequence_%lu.fa", seqId));
            stream << ">seq" << seqId << std::endl;
            return stream
        }*/

        bool finish_sequence(seqan::String<Dna5F> &sequence) {
            if (seqan::length(sequence) > 0) {
                std::cout << sequence << std::endl;
                seqan::clear(sequence);
                return true;
            } else {
                return false;
            }
        }

        /**
         * Take the current graph content and return it in SQL
         * format (GA4GH Spec).
         *
         * We will perform one depth first search traversal of the graph. While we will record
         * one long reference string, we will output all sidepaths on the way.
         */
        public:
        void toSQL() {
            
            // store all branch nodes on the way
            std::stack<BranchInfo> branchnodes;
            std::vector<bool> visited(last->size());
            for (std::vector<bool>::iterator it = visited.begin(); it != visited.end(); ++it) {
                *it = false;
            }
            std::vector<JoinInfo> joins;
            seqan::String<Dna5F> sequence;

            // for nodes with indegree > 1 we store sequence and index of the 
            // sequence that visited them, so we know where to anchor branches into it
            std::map<uint64_t, std::pair<uint64_t, uint64_t> > nodeId2seqPos; 
            uint64_t nodeId = 1; // start at source node
            uint64_t seqPos = 0; // position in currently traversed sequence, will increase with every visited node and be reset upon new branch
            uint64_t seqId = 1;  // first sequence ID is 1
            uint64_t seqCnt = 1; // number of total sequences
            bool isFirst = true;
            size_t out = outdegree(nodeId);
            BranchInfo branch;
            TAlphabet val;
            TAlphabet lastEdge = 0;
           // std::cout << " (" << seqId << ":" << seqPos << ") ";
           //
            while (out > 0 || branchnodes.size() > 0) {

                // we have reached the sink but there are unvisited nodes left on the stack
                if (out == 0) {
                    //std::cout << " (sink)" << std::endl;
                    if (branchnodes.size() == 0)
                        break;
                    // get new branch
                    branch = pop_branch(branchnodes, seqPos, nodeId, lastEdge, isFirst);
                    if (finish_sequence(sequence))
                        seqId = seqCnt += 1;
                    out = outdegree(nodeId);
                    if (debug)
                        fprintf(stderr, " -- popped %lu -- ", nodeId); 
                    joins.push_back(JoinInfo(branch.seqId, branch.seqPos, seqId, seqPos));
                    //fprintf(stderr, "1: join %lu %lu %lu %lu -- lastEdge %lu\n", branch.seqId, branch.seqPos, seqId, seqPos, lastEdge);
                    //std::cout << " (" << branch.seqId << ":" << branch.seqPos << ") ";
                }

                // we have not visited that node before
                if (!visited.at(nodeId)) {
                    visited.at(nodeId) = true;
                    seqPos += isFirst ? 0 : 1;
                    isFirst = false;
                    val = get_node_end_value(nodeId);
                    if (val % alph_size == 0) {
                        seqan::append(sequence, "$");
                        //std::cout << "$";
                    } else {
                        seqan::append(sequence, Dna5F((val % alph_size) - 1));
                        //std::cout << Dna5F((val % alph_size) - 1);
                    }
                    // store seq position of this node (we will join to it later)
                    if (indegree(nodeId) > 1) {
                        nodeId2seqPos.insert(std::make_pair(nodeId, std::make_pair(seqId, seqPos)));
                        //fprintf(stderr, "nodeId %lu seqId %lu seqPos %lu val %lu\n", nodeId, seqId, seqPos, val);
                    }
                }

                // there is only one child
                if (out == 1) {
                    uint64_t next = fwd(nodeId);
                    // the next node is new
                    if (!visited.at(next)) {
                        nodeId = next;
                        lastEdge = 0;
                    // we have seen the next node before
                    } else {
                        // look up the sequence info of that node
                        //std::cout << " (" << nodeId2seqPos[next].first << ":" << std::max(nodeId2seqPos[next].second, 1ul) - 1 << ") " << std::endl;
                        joins.push_back(JoinInfo(seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                        //fprintf(stderr, "2: join %lu %lu %lu %lu\n", seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second);
                        // there are no branches left
                        if (branchnodes.size() == 0)
                            break;
                        // otherwise go back to last branch
                        branch = pop_branch(branchnodes, seqPos, nodeId, lastEdge, isFirst);
                        if (finish_sequence(sequence))
                            seqId = seqCnt += 1;
                        out = outdegree(nodeId);
                        if (debug)
                            fprintf(stderr, " -- popped %lu -- ", nodeId); 
                        //std::cout << " (" << branch.seqId << ":" << branch.seqPos << ") ";
                        joins.push_back(JoinInfo(branch.seqId, branch.seqPos, seqId, seqPos));
                        //fprintf(stderr, "3: join %lu %lu %lu %lu\n", branch.seqId, branch.seqPos, seqId, seqPos);
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
                                    branchnodes.push(BranchInfo(nodeId, seqId, seqPos, lastEdge));
                                    if (debug)
                                        fprintf(stderr, " -- pushed %lu -- ", nodeId); 
                                }
                                nodeId = next;
                                updated = true;
                                lastEdge = 0;
                                break;
                            } else {
                                // look up the sequence info of that node
                                //std::cout << " -- (" << nodeId2seqPos[next].first << ":" << std::max(nodeId2seqPos[next].second, 1ul) - 1 << ") " << std::endl;
                                if (nodeId == next) 
                                    joins.push_back(JoinInfo(nodeId2seqPos[next].first, nodeId2seqPos[next].second, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                                else
                                    joins.push_back(JoinInfo(seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second));
                                //fprintf(stderr, "4: join %lu %lu %lu %lu next %lu\n", seqId, seqPos, nodeId2seqPos[next].first, nodeId2seqPos[next].second, next);
                            }
                        }
                    }
                    // we are done with this branch
                    // we should end up here, when nodes branch to themselves
                    if (!updated) {
                        // there are no branches left
                        if (branchnodes.size() == 0)
                            break;
                        // otherwise go back to last branch
                        branch = pop_branch(branchnodes, seqPos, nodeId, lastEdge, isFirst);
                        if (finish_sequence(sequence))
                            seqId = seqCnt += 1;
                        out = outdegree(nodeId);
                        if (debug)
                            fprintf(stderr, " -- popped %lu -- ", nodeId); 
                        //std::cout << " (" << branch.seqId << ":" << branch.seqPos << ") ";
                        joins.push_back(JoinInfo(branch.seqId, branch.seqPos, seqId, seqPos));
                        //fprintf(stderr, "5: join %lu %lu %lu %lu\n", branch.seqId, branch.seqPos, seqId, seqPos);
                    }
                    if (debug)
                        fprintf(stderr, " new nodeId: %lu\n", nodeId);
                }
                out = outdegree(nodeId);
            }
            // for completeness
            if (seqan::length(sequence) > 0)
                finish_sequence(sequence);

            for (size_t i = 0; i < joins.size(); ++i) {
                std::cout << "(" << joins.at(i).seqId1 << ":" << joins.at(i).seqPos1 << "--" << joins.at(i).seqId2 << ":" << joins.at(i).seqPos2 << ")" << std::endl;
            }
            /*if (out == 0 && branchnodes.size() == 0) {
                std::cout << " (sink)" << std::endl;
            } else {
                std::cout << std::endl;
            }*/

        }


};
#endif
