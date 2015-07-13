#ifndef __DBG_SUCCINCT_HPP__
#define __DBG_SUCCINCT_HPP__

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

/** 
 * We use seqan for an efficient representation of our alphabet.
 */
#include <seqan/sequence.h>
#include <seqan/basic.h>
using namespace seqan;

/** In future we would like to use some efficient low level implementation 
 * of rank and select on the arrays but for now we go the straighforward way to
 * see if the concept works for us. 
 *
 * In case we want to use libraries, we could use the sdsl-lite library to take 
 * care of the low level implementation of rank and select structures.
 * https://github.com/simongog/sdsl-lite
 */

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
using namespace std;
using namespace sdsl;

class DBG_succ {

    // define an extended alphabet for W --> somehow this does not work properly as expected
    typedef ModifiedAlphabet<Dna5, ModExpand<'X'> > Dna5F; 
   // typedef ModifiedAlphabet<Dna5F, ModExpand<'B'> > tmp_alph2; //A- 
   // typedef ModifiedAlphabet<tmp_alph2, ModExpand<'D'> > tmp_alph3; //C-
   // typedef ModifiedAlphabet<tmp_alph3, ModExpand<'H'> > tmp_alph4; //G-
   // typedef ModifiedAlphabet<tmp_alph4, ModExpand<'U'> > TAlphabet; //T-
   typedef unsigned short TAlphabet;

    private:
        // the bit array indicating the last outgoing edge of a node
        String<bool> last; // can also be sdsl bit_vector

        // the array containing the edge labels
        String<unsigned short> W;
        // convert our alphabet to int and use a compile time fixed integer vector
        //int_vector<> W;

        // the offset array to mark the offsets for the last column in the implicit node list
        vector<TAlphabet> F; 

        // k-mer size
        size_t k;
        // index of position that marks end in graph
        size_t p;
        // number of edges in the graph
        size_t m;
        // alphabet size
        size_t alph_size = 7;

        bool debug = false; //true;
    
    public:
        DBG_succ(unsigned k) : k(k) {

            append(last, false); //bit_vector(0, 0);
            append(last, true);
            //rank_support_v<1> last_rs(&last);

            append(W, (unsigned short) 0);
            append(W, (unsigned short) 0);
            //W = int_vector<>(0, 0, 3);
            //int_vector<0>::rank_1_type W_rs(&W);

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

            for (unsigned i = 0; i < length(seq); ++i) {
                //if (i > 0 && i % 100000 == 0) {
                if (i > 0 && i % 100 == 0) {
                    std::cout << "." << std::flush;
                    if (i % 1000 == 0)
                        fprintf(stdout, "%i - edges %i / nodes %i\n", i, (int) length(W) - 1, (int) rank_last((length(last) - 1)));
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

            fprintf(stdout, "edges %i / nodes %i\n", (int) length(W) - 1, (int) rank_last((length(last) - 1)));
            //String<Dna5F> test = "CCT";
            //fprintf(stdout, "\nindex of CCT: %i\n", (int) index(test));
        }

    private:
    
        /** 
         * Uses the object's array W, a given position i in W and a character c
         * from the alphabet and returns the number of occurences of c in W up to
         * position i.
         */
         size_t rank_W(size_t i, TAlphabet c) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            if (i > length(this->W) - 1)
                return length(this->W) - 1;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            for (size_t j = 1; j <= i; j++)
                cnt += (this->W[j] == c); 
            // return occurences
            return cnt;
        }

        /**
         * Uses the array W and gets a count i and a character c from 
         * the alphabet and returns the positions of the i-th occurence of c 
         * in W.
         */
        size_t select_W(size_t i, TAlphabet c) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            size_t j;
            for (j = 1; j < length(this->W); j++) {
                cnt += (this->W[j] == c); 
                if (cnt == i)
                    return j;
            }
            /**
             * If i is greater than the rank of c at end of W, return length(W). 
             * The length of W is one more than we have elements - this is to offset counting 
             * from 0.
             */
            return length(this->W);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the last index of a character c preceding in W[1..i].
         */
        size_t pred_W(size_t i, TAlphabet c) {
            return select_W(rank_W(i, c), c);
        }

        /**
         * This is a convenience function that returns for array W, a position i and 
         * a character c the first index of a character c in W[i..N].
         */
        size_t succ_W(size_t i, TAlphabet c) {
            return select_W(rank_W(i - 1, c) + 1, c);
        }

        /** 
         * Uses the object's array last and a position and
         * returns the number of set bits up to that postion.
         */
         size_t rank_last(size_t i) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            if (i > length(this->last) - 1)
                return length(this->last) - 1;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            for (size_t j = 1; j <= i; j++)
                cnt += this->last[j];
            // return occurences
            return cnt;
        }

        /**
         * Uses the object's array last and a given position i and
         * returns the position of the i-th set bit in last[1..i].
         */
        size_t select_last(size_t i) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            size_t cnt = 0, j = 0;
            for (j = 1; j < length(this->last); j++) {
                cnt += (this->last[j]);
                if (cnt == i)
                    return j;
            }
            // if i is greater than the rank of c at end of last, return length(last)
            return length(this->last);
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the last set bit in last[1..i].
         */
        size_t pred_last(size_t i) {
            return select_last(rank_last(i));
        }

        /**
         * This is a convenience function that returns for the object's array last
         * and a given position i the position of the first set bit in last[i..N].
         */
        size_t succ_last(size_t i) {
            return select_last(rank_last(i - 1) + 1);
        }


        /**
         * Given a position i in W and an edge label c, this function returns the
         * index of the node the edge is pointing to.
         */
        size_t outgoing(size_t i, TAlphabet c) {
            size_t j1 = pred_W(i, c);
            size_t j2 = pred_W(i, c + alph_size);
            //fprintf(stdout, "i %i, c %i, j1 %i, j2 %i\n", (int) i, (int) c, (int) j1, (int) j2);
            size_t j = (j1 < j2) ? j2 : j1;
            if (j == 0 || j == length(W))
                return 0;
            j = fwd(j);
            if (j == 0 || j == length(W))
                return 0;
            return j;
        }

        size_t index(String<Dna5F> &s) {
            // init range
            size_t rl = succ_last(F[(TAlphabet) ordValue(s[0]) + 1] + 1); // lower bound
            size_t ru = F[(TAlphabet) ordValue(s[0]) + 2];                // upper bound
            // update range iteratively while scanning through s
            for (size_t i = 1; i < length(s); i++) {
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
        TAlphabet get_node_end_value(size_t i) {
            if (i == 0)
                return 0;
            for (TAlphabet j = 0; j < F.size(); j++) {
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
        size_t bwd(size_t i) {
            // get value of last position in node i
            TAlphabet c = get_node_end_value(i);
            // get the offset for the last position in node i
            size_t o = F[c];
            //fprintf(stdout, "i %i c %i o %i rank(i) %i rank(o) %i\n", (int) i, (int) c, (int) o, (int) rank_last(i), (int) rank_last(o));
            // compute the offset for this position in W and select it
            return select_W(rank_last(i) - rank_last(o), c);
        }


        /**
         * This functions gets a position i reflecting the r-th occurence of the corresponding
         * character c in W and returns the position of the r-th occurence of c in last.
         */
        size_t fwd(size_t i) {
            // get value of W at position i
            TAlphabet c = W[i]; 
            // get the offset for position c
            size_t o = F[c];
            // get the rank of c in W at position i
            size_t r = rank_W(i, c);
            // select the index of the position in last that is rank many positions after offset
            return select_last(rank_last(o) + r);
        }

        /**
         * Given a position i, this function returns the boundaries of the interval
         * of nodes identical to node i (ignoring the values in W).
         */
        pair<size_t, size_t> get_equal_node_range(size_t i) {
            return make_pair(pred_last(i - 1) + 1, succ_last(i));
        }

        /**
         * This is a debug function that prints the current state of the graph arrays to
         * the screen.
         */
        void print_state() {

            fprintf(stderr, "W:\n");
            for (size_t i = 0; i < length(W); i++) {
                fprintf(stderr, "\t%i", (int) W[i]);
                if (i == p)
                    fprintf(stderr, "*");
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "last:\n");
            for (size_t i = 0; i < length(last); i++)
                fprintf(stderr, "\t%i", (int) last[i]);
            fprintf(stderr, "\n");

            fprintf(stderr, "F:\n");
            for (size_t i = 0; i < F.size(); i++)
                fprintf(stderr, "\t%i", (int) F[i]);
            fprintf(stderr, "\n");

        }

        /** This function takes a character c and appends it to the end of the graph sequence
         * given that the corresponding note is not part of the graph yet.
         */
        void append_pos(TAlphabet c) {

            // check that the last position of the graph is indeed a terminal
            assert(W[p] == 0);
            TAlphabet c_p = get_node_end_value(p);
            // get range of identical nodes (without W) pos current end position
            pair<size_t, size_t> R = get_equal_node_range(this->p);
            //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

            // get position of first occurence of c in W after p
            size_t next_c = succ_W(p, c);
            //fprintf(stdout, "p %i, c %i, next_c %i\n", (int) p, (int) c, (int) next_c);
            // check if c is part of range
            bool exist_c = (next_c <= R.second); // && W[next_c] != 0 && !is_terminal_node(next_c));
            bool skip = false;
            if (!exist_c) {
                // get position of first occurence of c- in W after p
                next_c = succ_W(p, c + alph_size);
                // check if c- is part of range
                exist_c = (next_c <= R.second);
                //skip = (next_c <= R.second); // --> node already exists - we are in a self-loop
            }

            /**
             * if the character already exists in the range, we delete the terminal symbol
             * at p, insert c at fwd(next_c) and update p.
             */
            //fprintf(stdout, "exist_c %i\n", (int) exist_c);
            if (exist_c) {
                // we need an addtional pred_W here as fwd is not defined if next_c is element of the minus set
                size_t p_new = fwd(pred_W(next_c, c));
                // remove old terminal symbol
                erase(last, p);
                erase(W, p);
                // adapt position if altered by previous deletion
                p_new -= (p < p_new);
                // insert new terminal symbol 
                // we have to insert 0 into last as the node already existed in the range 
                // and the terminal symbol is always first
                insert(last, p_new, false);//p == p_new);
                insert(W, p_new, 0);
                // update new terminal position
                p = p_new;
                // take care of updating the offset array F
                update_F(c_p, false);
                //assert(get_node_end_value(p) == c);
                update_F(c, true);
            } else if (!skip) {
                /**
                 * We found that c does not yet exist in the current range and now have to
                 * figure out if we need to add c or c- to the range.
                 * To do this, we check if there is a previous position j1 with W[j1] == c
                 * whose node shares a k-1 suffix with the current node. If yes, we add c- 
                 * instead of c.
                 */
                // get position of last occurence of c before p (including p - 1)
                size_t last_c = pred_W(p - 1, c);
                // if this position exists
                if (last_c > 0) {
                    size_t x = fwd(last_c);
                    assert(last[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

                    // check, if there are any c or c- symbols following after position p
                    size_t next_c = succ_W(p + 1, c);
                    size_t next_cm = succ_W(p + 1, c + alph_size);
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
                    if (next_c < length(W)) {
                        minus2 = compare_node_suffix(p, next_c);
                        if (minus2)
                            W[next_c] += alph_size;
                    }

                    W[p] = minus1 ? c + alph_size : c;
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
                        insert(last, x, false, Generous());
                        insert(W, x, 0, Generous());
                    // no node shares a k-1 suffix with last_c and thus the new node comes after
                    // the forward of last_c (as the current node came after last_c as well)
                    } else {
                        p = x + 1;
                        insert(last, x + 1, true, Generous());
                        insert(W, x + 1, 0, Generous());
                    }
                } else {
                    size_t x = F[c];
                    size_t next_c = succ_W(p + 1, c);
                    bool minus = false;
                    if (next_c < length(W)) {
                        minus = compare_node_suffix(p, next_c);
                    }
                    W[p] = c;
                    if (p - R.second > 0) {
                        sort_W_locally(p, R.second);
                    }
                    p = x + 1;
                    if (minus) {
                        W[next_c] += alph_size;
                        insert(last, x + 1, 0, Generous());
                    } else {
                        insert(last, x + 1, 1, Generous());
                    }
                    insert(W, x + 1, 0, Generous());
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
                while (W[p] != 0)
                    p--;
                assert(W[p] == 0);
            }
        }


        /** 
         * This function gets two node indices and returns if the
         * node labels share a k-1 suffix.
         */
        bool compare_node_suffix(size_t i1, size_t i2) {
            for (size_t ii = 0; ii < k-1; ii++) {
                //cout << "node1 - " << i1 << ": " << Dna5F(get_node_end_value(i1) % alph_size - 1) << endl; 
                //cout << "node2 - " << i2 << ": " << Dna5F(get_node_end_value(i2) % alph_size - 1) << endl; 
                if (get_node_end_value(i1) != get_node_end_value(i2)) {
                    return false;
                }
                //cout << "succ (i1): " << succ_last(i1) << " succ (i2): " << succ_last(i2) << endl;
                i1 = bwd(succ_last(i1));
                i2 = bwd(succ_last(i2));
            }
            return true;
        }

        bool is_terminal_node(size_t i) {
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
            for (size_t i = 1; i < length(W); i++) {
            //while (i > 0 && get_node_end_value(i) > 0) {
                if (W[i] % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((W[i] % alph_size) - 1);
            }
            cout << endl;

            for (size_t i = 1; i < length(W); i++) {
                if (p == i)
                    fprintf(stdout, "*");
                else
                    fprintf(stdout, " ");
            }
            cout << endl;

            size_t j;
            for (size_t i = 1; i < length(W); i++) {
                j = get_node_end_value(i);
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (size_t i = 1; i < length(W); i++) {
                j = get_node_end_value(bwd(succ_last(i)));
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (size_t i = 1; i < length(W); i++) {
                j = get_node_end_value(bwd(succ_last(bwd(succ_last(i)))));
                if (j % alph_size == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5F((j % alph_size) - 1);
            }
            cout << endl;
            for (size_t i = 1; i < length(last); i++) {
                fprintf(stdout, "%i", (int) last[i]);
            }
            cout << endl;
        }

        /**
         * This function gets a local range in W from lower bound l
         * to upper bound u and swaps the inserted element to the
         * righ location.
         */
        void sort_W_locally(size_t l, size_t u) {
            for (size_t s = u; s > l; s--) {
                TAlphabet tmp;
                if ((W[s] % alph_size) < (W[s-1] % alph_size)) {
                    tmp = W[s-1];
                    W[s-1] = W[s];
                    W[s] = tmp;
                }
            }
        }
};
#endif
