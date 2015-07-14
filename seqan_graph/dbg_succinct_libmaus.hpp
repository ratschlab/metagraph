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
#include <assert.h>
#include <iostream>

#include <libmaus/bitbtree/bitbtree.hpp>
#include <libmaus/wavelet/DynamicWaveletTree.hpp>

/** 
 * We use seqan for an efficient representation of our alphabet.
 */
#include <seqan/sequence.h>
#include <seqan/basic.h>
using namespace seqan;

using namespace std;

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef ModifiedAlphabet<Dna5, ModExpand<'X'> > Dna5F; 
   typedef uint64_t TAlphabet;

    private:
        // the bit array indicating the last outgoing edge of a node
        libmaus::bitbtree::BitBTree<6, 64> *last2 = new libmaus::bitbtree::BitBTree<6, 64>();


        // the array containing the edge labels
        libmaus::wavelet::DynamicWaveletTree<6, 64> *W2 = new libmaus::wavelet::DynamicWaveletTree<6, 64>(4); // 4 is log (sigma)

        // the offset array to mark the offsets for the last column in the implicit node list
        vector<TAlphabet> F; 

        // k-mer size
        size_t k;
        // index of position that marks end in graph
        uint64_t p;
        // number of edges in the graph
        uint64_t m;
        // alphabet size
        size_t alph_size = 7;

        uint64_t W2_size = 0;

        bool debug = false; //true;
    
    public:
        DBG_succ(size_t k) : k(k) {

            last2->insertBit(0, true);
            last2->insertBit(0, false);

            W2->insert(0, 0);
            W2->insert(0, 0);
            W2_size += 2;

            //F = vector<unsigned int>(5, 0);
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
                cout << "======================================" << endl;
            }

            for (uint64_t i = 0; i < length(seq); ++i) {
                //if (i > 0 && i % 100000 == 0) {
                if (i > 0 && i % 1000 == 0) {
                    std::cout << "." << std::flush;
                    if (i % 10000 == 0) {
                        fprintf(stdout, "%lu - edges %lu / nodes %lu\n", i, W2_size - 1, rank_last2((last2->size() - 1)));
                    }
                }
                //fprintf(stdout, "appending %i\n", (int) ordValue(seq[i]));
                //cerr << "seq[i] " << seq[i] << endl;
                //cerr << "seq[i] ord " << ordValue(seq[i]) + 1 << endl;
                append_pos((TAlphabet) ordValue(seq[i]) + 1);
                if (debug) {
                    print_seq();
                    print_state();
                    cout << "======================================" << endl;
                }
            }

            for (size_t j = 0; j < k - 1; j++) {
                append_pos(6);
                if (debug) {
                    print_seq();
                    print_state();
                    cout << "======================================" << endl;
                }
            }

            fprintf(stdout, "edges %lu / nodes %lu\n", W2_size - 1, rank_last2((last2->size() - 1)));
            //String<Dna5F> test = "CCT";
            //fprintf(stdout, "\nindex of CCT: %i\n", (int) index(test));
        }

    private:
    
        /** 
         * Uses the object's array W, a given position i in W and a character c
         * from the alphabet and returns the number of occurences of c in W up to
         * position i.
         */
        uint64_t rank_W2(uint64_t i, TAlphabet c) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            if (i > W2_size - 1)
                return W2_size - 1;
            return W2->rank(c, i);
        }

        /**
         * Uses the array W and gets a count i and a character c from 
         * the alphabet and returns the positions of the i-th occurence of c 
         * in W.
         */
        uint64_t select_W2(uint64_t i, TAlphabet c) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            //fprintf(stderr, "query select W2 -- c: %i i: %lu return: %lu \n", c, i-1+(c==0), W2->select(c, i-1+(c==0)));
            return min(W2->select(c, i - 1 + (c == 0)), W2_size);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the last index of a character c preceding in W[1..i].
         */
        uint64_t pred_W2(uint64_t i, TAlphabet c) {
            return select_W2(rank_W2(i, c), c);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the first index of a character c in W[i..N].
         */
        uint64_t succ_W2(uint64_t i, TAlphabet c) {
            return select_W2(rank_W2(i - 1, c) + 1, c);
        }

        /** 
         * Uses the object's array last and a position and
         * returns the number of set bits up to that postion.
         */
        uint64_t rank_last2(uint64_t i) {
            // deal with  border conditions
            if (i <= 0)
                return 0;
            //if (i > last2->size() - 1) {
            //    fprintf(stderr, "i %lu size %lu\n", i, last2->size());
            //    return last2->size() - 1;
            //}
            return last2->rank1(i);
        }

        /**
         * Uses the object's array last and a given position i and
         * returns the position of the i-th set bit in last[1..i].
         */
        uint64_t select_last2(uint64_t i) {
            // deal with  border conditions
            if (i <= 0)
                return 0;
            //fprintf(stderr, "i %lu size %lu\n", i, last2->size());
            //if (i >= last2->size())
            //    return last2->size();
            // for some reason the libmaus select is 0 based ...
            return last2->select1(i - 1);
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the last set bit in last[1..i].
         */
        uint64_t pred_last2(uint64_t i) {
            return select_last2(rank_last2(i));
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the first set bit in last[i..N].
         */
        uint64_t succ_last2(uint64_t i) {
            //if (i > 65000) {
            //    fprintf(stderr, "i: %lu ", i-1);
            //    fprintf(stderr, "rank: %lu\n", rank_last2(i - 1));
            //}
            return select_last2(rank_last2(i - 1) + 1);
        }

        /**
         * Given a position i in W and an edge label c, this function returns the
         * index of the node the edge is pointing to.
         */
        uint64_t outgoing(uint64_t i, TAlphabet c) {
            uint64_t j1 = pred_W2(i, c);
            uint64_t j2 = pred_W2(i, c + alph_size);
            //fprintf(stdout, "i %i, c %i, j1 %i, j2 %i\n", (int) i, (int) c, (int) j1, (int) j2);
            uint64_t j = (j1 < j2) ? j2 : j1;
            if (j == 0 || j == W2_size)
                return 0;
            j = fwd(j);
            if (j == 0 || j == W2_size)
                return 0;
            return j;
        }

        uint64_t index(String<Dna5F> &s) {
            // init range
            uint64_t rl = succ_last2(F[(TAlphabet) ordValue(s[0]) + 1] + 1);
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
            //fprintf(stderr, "This does not make sense: %i\n", (int) i);
            return F.size() - 1;
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
            if (select_W2(rank_last2(i) - rank_last2(o), c) > 1000000) {
                uint64_t tmp = rank_last2(i) - rank_last2(o);
                fprintf(stdout, "%lu %lu\n", select_W2(tmp, c), select_W2(rank_last2(i) - rank_last2(o), c));
                fprintf(stdout, "i %lu c %i o %lu rank(i) %lu rank(o) %lu difference %lu\n", i, (int) c, o, rank_last2(i), rank_last2(o), rank_last2(i) - rank_last2(o));
            }
            // compute the offset for this position in W and select it
            return select_W2(rank_last2(i) - rank_last2(o), c);
        }


        /**
         * This functions gets a position i reflecting the r-th occurence of the corresponding
         * character c in W and returns the position of the r-th occurence of c in last.
         */
        uint64_t fwd(uint64_t i) {
            // get value of W at position i
            TAlphabet c = (*W2)[i]; 
            // get the offset for position c
            uint64_t o = F[c];
            // get the rank of c in W at position i
            uint64_t r = rank_W2(i, c);
            // select the index of the position in last that is rank many positions after offset
            return select_last2(rank_last2(o) + r);
        }

        /**
         * Given a position i, this function returns the boundaries of the interval
         * of nodes identical to node i (ignoring the values in W).
         */
        pair<uint64_t, uint64_t> get_equal_node_range(uint64_t i) {
            return make_pair(pred_last2(i - 1) + 1, succ_last2(i));
        }

        /**
         * This is a debug function that prints the current state of the graph arrays to
         * the screen.
         */
        void print_state() {

            fprintf(stderr, "W2:\n");
            for (uint64_t i = 0; i < W2_size; i++) {
                fprintf(stderr, "\t%lu", (*W2)[i]);
                if (i == p)
                    fprintf(stderr, "*");
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "last2:\n");
            for (uint64_t i = 0; i < last2->size(); i++)
                fprintf(stderr, "\t%i", (int) (*last2)[i]);
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
            assert((*W2)[p] == 0);
            TAlphabet c_p = get_node_end_value(p);
            // get range of identical nodes (without W) pos current end position
            pair<uint64_t, uint64_t> R = get_equal_node_range(this->p);
            //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

            // get position of first occurence of c in W after p
            uint64_t next_c = succ_W2(p, c);
            //fprintf(stdout, "p %i, c %i, next_c %i\n", (int) p, (int) c, (int) next_c);
            // check if c is part of range
            bool exist_c = (next_c <= R.second); // && W[next_c] != 0 && !is_terminal_node(next_c));
            if (!exist_c) {
                // get position of first occurence of c- in W after p
                next_c = succ_W2(p, c + alph_size);
                // check if c- is part of range
                exist_c = (next_c <= R.second);
            }

            /**
             * if the character already exists in the range, we delete the terminal symbol
             * at p, insert c at fwd(next_c) and update p.
             */
            //fprintf(stdout, "exist_c %i\n", (int) exist_c);
            if (exist_c) {
                // we need an addtional pred_W here as fwd is not defined if next_c is element of the minus set
                uint64_t p_new = fwd(pred_W2(next_c, c));
                // remove old terminal symbol
                last2->deleteBit(p);
                W2->remove(p);
                W2_size--;
                // adapt position if altered by previous deletion
                p_new -= (p < p_new);
                // insert new terminal symbol 
                // we have to insert 0 into last as the node already existed in the range 
                // and the terminal symbol is always first
                last2->insertBit(p_new, false);
                W2->insert(0, p_new);
                W2_size++;
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
                uint64_t last_c = pred_W2(p - 1, c);
                // if this position exists
                if (last_c > 0) {
                    uint64_t x = fwd(last_c);
                    assert((*last2)[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

                    // check, if there are any c or c- symbols following after position p
                    uint64_t next_c = succ_W2(p + 1, c);
                    uint64_t next_cm = succ_W2(p + 1, c + alph_size);
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
                    if (next_c < W2_size) {
                        minus2 = compare_node_suffix(p, next_c);
                        if (minus2) {
                            replaceW(next_c, (*W2)[next_c] + alph_size);
                        }
                    }

                    replaceW(p, minus1 ? c + alph_size : c);
                    // after we are done, assert that the order within the range we created 
                    // is still valid within W
                    if (p - R.second > 0) {
                        sort_W_locally(p, R.second);
                    }

                    // if one of the minuses is true, at least one node shares a k-1 suffix with last_c
                    // and has the same character in W --> this the fwd(last_c) creates a node that is not unique
                    // as we insert at the beginning of the range, last is false. We insert at position x
                    // as this is the beginning of the range (last_c)
                    if (minus1 || minus2) {
                        p = x;
                        last2->insertBit(x, false);
                        W2->insert(0, x);
                        W2_size++;
                    // no node shares a k-1 suffix with last_c and thus the new node comes after
                    // the forward of last_c (as the current node came after last_c as well)
                    } else {
                        p = x + 1;
                        last2->insertBit(x + 1, true);
                        W2->insert(0, x + 1);
                        W2_size++;
                    }
                } else {
                    uint64_t x = F[c];
                    uint64_t next_c = succ_W2(p + 1, c);
                    bool minus = false;
                    if (next_c < W2_size) {
                        minus = compare_node_suffix(p, next_c);
                    }
                    replaceW(p, c);
                    if (p - R.second > 0) {
                        sort_W_locally(p, R.second);
                    }
                    p = x + 1;
                    if (minus) {
                        replaceW(next_c, (*W2)[next_c] + alph_size);
                        last2->insertBit(x + 1, false);
                    } else {
                        last2->insertBit(x + 1, true);
                    }
                    W2->insert(0, x + 1);
                    W2_size++;
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
                while ((*W2)[p] != 0)
                    p--;
                assert((*W2)[p] == 0);
            }
        }


        /** 
         * This function gets two node indices and returns if the
         * node labels share a k-1 suffix.
         */
        bool compare_node_suffix(uint64_t i1, uint64_t i2) {
            for (size_t ii = 0; ii < k-1; ii++) {
                //cout << "node1 - " << i1 << ": " << Dna5F(get_node_end_value(i1) % alph_size - 1) << endl; 
                //cout << "node2 - " << i2 << ": " << Dna5F(get_node_end_value(i2) % alph_size - 1) << endl; 
                if (get_node_end_value(i1) != get_node_end_value(i2)) {
                    return false;
                }
                //cout << "succ (i1): " << succ_last(i1) << " succ (i2): " << succ_last(i2) << endl;
                i1 = bwd(succ_last2(i1));
                i2 = bwd(succ_last2(i2));
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

        /*
         * Returns the sequence stored in W
         */
        void print_seq() {

            for (uint64_t i = 1; i < W2_size; i++) {
                if ((*W2)[i] % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F(((*W2)[i] % alph_size) - 1);
            }
            cout << endl;

            for (uint64_t i = 1; i < W2_size; i++) {
                if (p == i)
                    fprintf(stdout, "*");
                else
                    fprintf(stdout, " ");
            }
            cout << endl;

            size_t j;
            for (uint64_t i = 1; i < W2_size; i++) {
                j = get_node_end_value(i);
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (uint64_t i = 1; i < W2_size; i++) {
                j = get_node_end_value(bwd(succ_last2(i)));
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (uint64_t i = 1; i < W2_size; i++) {
                j = get_node_end_value(bwd(succ_last2(bwd(succ_last2(i)))));
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (uint64_t i = 1; i < last2->size(); i++) {
                fprintf(stdout, "%i", (int) (*last2)[i]);
            }
            cout << endl;
        }

        /**
         * This function gets a local range in W from lower bound l
         * to upper bound u and swaps the inserted element to the
         * righ location.
         */
        void sort_W_locally(uint64_t l, uint64_t u) {
            for (uint64_t s = u; s > l; s--) {
                TAlphabet tmp;
                if (((*W2)[s] % alph_size) < ((*W2)[s-1] % alph_size)) {
                    tmp = (*W2)[s-1];
                    replaceW(s-1, (*W2)[s]);
                    replaceW(s, tmp);
                }
            }
        }

        /** This is a convenience function to replace value at
         * position i in W with val.
         */
        void replaceW(size_t i, TAlphabet val) {
            W2->remove(i);
            W2->insert(val, i);
        }
};
#endif
