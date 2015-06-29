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

/** In future we would like to loose some efficient low level implementation 
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
   // typedef ModifiedAlphabet<Dna5, ModExpand<'$'> > Dna5F; 
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
        size_t alph_size = 12;
    
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
            String<Dna5> seq
        ) {
            for (size_t j = 0; j < k - 1; j++)
                append_pos(0);

            for (unsigned i = 0; i < length(seq); ++i) {
                //if (i > 0 && i % 100000 == 0) {
                if (i > 0 && i % 1000 == 0) {
                    std::cout << "." << std::flush;
                    if (i % 1000000 == 0)
                        fprintf(stdout, "%i\n", i);
                }
                //fprintf(stdout, "appending %i\n", (int) ordValue(seq[i]));
                //cerr << "seq[i] " << seq[i] << endl;
                //cerr << "seq[i] ord " << ordValue(seq[i]) + 1 << endl;
                append_pos((TAlphabet) ordValue(seq[i]) + 1);
            }
            print_seq();

            for (size_t j = 0; j < k - 1; j++)
                append_pos(0);


        }

    private:
    
        /** 
         * Gets a reference to the array W, a position i in W and a character c
         * from the alphabet and returns the number of occurences of c in W up to
         * position i.
         */
         //size_t rank_W(String<unsigned short> &W, size_t i, unsigned short &c) {
         size_t rank_W(String<TAlphabet> &W, size_t i, TAlphabet c) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            if (i > length(W) - 1)
                return length(W) - 1;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            for (size_t j = 1; j <= i; j++)
                cnt += (W[j] == c); 
            // return occurences
            return cnt;
        }

        /**
         * Gets a reference to the array W, a count i and a character c from 
         * the alphabet and returns the positions of the i-th occurence of c 
         * in W.
         */
        size_t select_W(String<TAlphabet> &W, size_t i, TAlphabet c) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            size_t j;
            for (j = 1; j < length(W); j++) {
                cnt += (W[j] == c); 
                if (cnt == i)
                    return j;
            }
            /**
             * If i is greater than the rank of c at end of W, return length(W). 
             * The length of W is one more than we have elements - this is to offset counting 
             * from 0.
             */
            return j;
        }

        /**
         * This is a convenience function that returns for a given array W, a position i and 
         * a character c the last index of a character c preceding in W[1..i].
         */
        size_t pred_W(String<TAlphabet> &W, size_t i, TAlphabet c) {
            return select_W(W, rank_W(W, i, c), c);
        }

        /**
         * This is a convenience function that returns for a given array W, a position i and 
         * a character c the first index of a character c in W[i..N].
         */
        size_t succ_W(String<TAlphabet> &W, size_t i, TAlphabet c) {
            return select_W(W, rank_W(W, i - 1, c) + 1, c);
        }

        /** 
         * Gets a reference to the array last and a position and
         * returns the number of set bits up to that postion.
         */
         size_t rank_last(String<bool> &last, size_t i) {

            // deal with  border conditions
            if (i <= 0)
                return 0;
            if (i > length(last) - 1)
                return length(last) - 1;

            // count occurences of c and store them in cnt
            size_t cnt = 0;
            for (size_t j = 1; j <= i; j++)
                cnt += last[j];
            // return occurences
            return cnt;
        }

        /**
         * Gets a reference to the array last and a position i and
         * returns the position of the i-th set bit in last[1..i].
         */
        size_t select_last(String<bool> &last, size_t i) {
            
            // deal with  border conditions
            if (i <= 0)
                return 0;

            // count occurences of c and store them in cnt
            size_t cnt = 0, j = 0;
            for (j = 1; j < length(last); j++) {
                cnt += (last[j]);
                if (cnt == i)
                    return j;
            }
            // if i is greater than the rank of c at end of W, return length(W)
            return j;
        }

        /**
         * This is a convenience function that returns for a given binary array last
         * and a position i the position of the last set bit in last[1..i].
         */
        size_t pred_last(String<bool> &last, size_t i) {
            return select_last(last, rank_last(last, i));
        }

        /**
         * This is a convenience function that returns for a given binary array last
         * and a position i the position of the first set bit in last[i..N].
         */
        size_t succ_last(String<bool> &last, size_t i) {
            return select_last(last, rank_last(last, i - 1) + 1);
        }


        /**
         * Using the offset structure F this function returns the value of the last 
         * position of node i.
         */
        TAlphabet get_node_end_value(size_t i) {
            for (TAlphabet j = 0; j < F.size(); j++)
                if (F[j] >= i)
                    return j - 1;
            return F.size();
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
            // compute the offset for this position in W and select it
            return select_W(W, rank_last(last, i) - rank_last(last, o), c);
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
            size_t r = rank_W(W, i, c);
            //fprintf(stderr, "i %i, c %i, o %i, r %i\n", (int) i, (int) c, (int) o, (int) r);
            // select the index of the position in last that is rank many positions after offset
            return select_last(last, rank_last(last, o) + r);
        }

        /**
         * Given a position i, this function returns the boundaries of the interval
         * of nodes identical to node i (ignoring the values in W).
         */
        pair<size_t, size_t> get_equal_node_range(size_t i) {
            return make_pair(pred_last(last, i-1) + 1, succ_last(last, i));
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
            //print_state();

            // get position of first occurence of c in W after p
            size_t pc = succ_W(W, p, c);
            // check if c is part of range
            bool exist_c = (pc <= R.second);
            if (~exist_c) {
                // get position of first occurence of c- in W after p
                pc = succ_W(W, p, c + alph_size);
                // check if c- is part of range
                exist_c = (pc <= R.second);
            }

            /**
             * if the character already exists in the range, we delete the terminal symbol
             * at p, insert c at fwd(pc) and update p.
             */
            if (exist_c) {
                size_t p_new = fwd(pc);
                // remove old terminal symbol
                erase(last, p);
                erase(W, p);
                // adapt position if altered by previous deletion
                p_new -= (p < p_new);
                // insert new terminal symbal 
                insert(last, p_new, false);
                insert(W, p_new, 0);
                // update new terminal position
                p = p_new;
                // take care of updating the offset array F
                update_F(c_p, false);
                update_F(c, true);
            } else {
                /**
                 * We found that c does not yet exist in the current range and now have to
                 * figure out if we need to add c or c- to the range.
                 * To do this, we check if there is a previous position j1 with W[j1] == c
                 * having the same node sequence. If yes, we add c- instead of c.
                 */
                size_t j1 = pred_W(W, p - 1, c);
                if (j1 > 0) {
                    size_t x = fwd(j1);
                    assert(last[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

                    // check, if there are any c or c- symbols following after position p
                    size_t j2 = succ_W(W, p + 1, c);
                    size_t j3 = succ_W(W, p + 1, c + 6);
                    // there is no c between p and j3 and j3 is a c- ==> we should minus a c- 
                    bool minus1 = (j3 < j2);
                    // check, if we share a k-1 suffix with j1
                    if (~minus1) {
                        minus1 = compare_node_suffix(p, j1);
                    }
                    if (minus1)
                        W[p] = c + 6;
                    else
                        W[p] = c;

                    bool minus2 = compare_node_suffix(p, j2);
                    if ((j2 < m) && minus2) {
                        W[j2] = (W[j2] % 6) + 6;
                    }
                    if (minus1 || minus2) {
                        insert(last, x, false, Generous());
                        insert(W, x, 0, Generous());
                        p = x;
                    } else {
                        insert(last, x + 1, true, Generous());
                        insert(W, x + 1, 0, Generous());
                        p = x + 1;
                    }
                } else {
                    size_t x = F[c];
                    size_t j2 = succ_W(W, p + 1, c);
                    bool minus = false;
                    if (j2 < length(W)) {
                        minus = compare_node_suffix(p, j2);
                    }
                    W[p] = c;
                    //fprintf(stderr, "c %i\n", (int) c);
                    //fprintf(stderr, "j2 %i\n", (int) j2);
                    if (minus) {
                        W[j2] = c + 6;
                        insert(last, x + 1, 0, Generous());
                    } else {
                        insert(last, x + 1, 1, Generous());
                    }
                    insert(W, x + 1, 0, Generous());
                    p = x + 1;
                }
                m++;
                update_F(c, true);
            }
        }


        /** 
         * This function gets two node indices and returns if the
         * node labels share a k-1 suffix.
         */
        bool compare_node_suffix(size_t i1, size_t i2) {
            bool share_suffix = true;
            for (size_t ii = 0; ii < k-1; ii++) {
                if (get_node_end_value(i1) != get_node_end_value(i2)) {
                    share_suffix = false;
                    break;
                }
                i1 = bwd(i1);
                i2 = bwd(i2);
            }
            return share_suffix;
        }


        /**
         * This function gets a value of the alphabet c and updates the offset of 
         * all following values by +1 is positive is true and by -1 otherwise.
         */
        void update_F(TAlphabet c, bool positive) {
            for (TAlphabet i = c+1; i < (TAlphabet) alph_size; i++)
                F[i] += (2 * (TAlphabet) positive - 1);
        }

        /*
         * Returns the sequence stored in W
         */
        void print_seq() {
            for (size_t i = 1; i < length(W); i++) {
            //while (i > 0 && get_node_end_value(i) > 0) {
                if (W[i] == 0)
                    fprintf(stdout, "$");
                else
                    cout << Dna5((W[i] % 6) - 1);
            }
            cout << endl;
        }
};
#endif
