#include "construct.hpp"

#include "kseq.h"
#include "dbg_succinct_boost.hpp"
#include "dbg_succinct_libmaus.hpp"

namespace construct {

#ifdef DBGDEBUG
    bool debug = true;
#else
    bool debug = false;
#endif 

    // add a full sequence to the graph
    void add_seq(DBG_succ* G, kstring_t &seq) {

        if (debug) {
            G->print_seq();
            G->print_state();
            std::cout << "======================================" << std::endl;
        }

        // Padding of the input genome / read
        if (G->W->n == 2) {
            for (size_t j = 0; j < G->k; j++) {
                append_pos(G, 6);
                if (debug) {
                    G->print_seq();
                    G->print_state();
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
                    fprintf(stdout, "%lu - edges %lu / nodes %lu\n", i, G->get_edge_count(), G->get_node_count());
                }
            }
            // if (debug) {
            //    fprintf(stdout, "appending %i\n", (int) ordValue(seq[i]));
            //    cerr << "seq[i] " << seq[i] << std::endl;
            //    cerr << "seq[i] ord " << ordValue(seq[i]) + 1 << std::endl;
            // }

            append_pos(G, G->get_alphabet_number(seq.s[i]));

            if (debug) {
                G->print_seq();
                G->print_state();
                std::cout << "======================================" << std::endl;
            }
        }

        // Padding after sequence to get back into default state.
        for (size_t j = 0; j < G->k; j++) {
            append_pos(G, 6);
            if (debug) {
                G->print_seq();
                G->print_state();
                std::cout << "======================================" << std::endl;
            }
        }

        fprintf(stdout, "edges %lu / nodes %lu\n", G->get_edge_count(), G->get_node_count());
    }


    void add_seq_alt(DBG_succ* G, kstring_t &seq, bool bridge, unsigned int parallel, std::string suffix) {

        if (debug) {
            G->print_seq();
            G->print_state();
            std::cout << "======================================" << std::endl;
        }

        /*
        char *nt_lookup = (char*)malloc(128);
        uint8_t def = 5;
        memset(nt_lookup, def, 128);
        for (size_t i=0;i<alph_size;++i) {
            nt_lookup[(uint8_t)alphabet[i]]=i;
            nt_lookup[(uint8_t)tolower(alphabet[i])]=i;
        }
        */
        const char nt_lookup[128] = {
                        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
                        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
                        5, 5, 5, 5,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
                        5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5, 
                        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5, 
                        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5, 
                        5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5, 
                        5, 5, 5, 5,  4, 4, 5, 5,  6, 5, 5, 5,  5, 5, 5, 5 
        };   

        seqtokmer(G->kmers, seq.s, seq.l, G->k, nt_lookup, G->alphabet, bridge, parallel, suffix);
        //free(nt_lookup);
    }


    void construct_succ(DBG_succ* G, unsigned int parallel) {
        omp_set_num_threads(std::max((int)parallel,1));
        __gnu_parallel::sort(G->kmers.begin(),G->kmers.end());
        //std::sort(G->kmers.begin(), G->kmers.end());
        G->kmers.erase(std::unique(G->kmers.begin(), G->kmers.end() ), G->kmers.end() );

        //DEBUG: output kmers in current bin
        /*
        for (size_t i=0;i<G->kmers.size();++i) {
            char* curseq = kmertos(G->kmers[i], G->alphabet, G->alph_size);
            std::cerr << G->kmers[i] << "\t" << curseq+1 << " " << curseq[0] << "\n";
            free(curseq);
        }
        */

        size_t curpos = G->W_stat.size();
        G->W_stat.resize(G->W_stat.size()+G->kmers.size());
        G->last_stat_safe.resize(G->last_stat_safe.size()+G->kmers.size(), true);
        
        #pragma omp parallel num_threads(parallel)
        {    
            #pragma omp for nowait
            for (size_t i=0;i<G->kmers.size();++i) {
                //set last
                if (i+1 < G->kmers.size()) {
                    bool dup = compare_kmer_suffix(G->kmers[i], G->kmers[i+1]);
                    if (dup) {
                        G->last_stat_safe[curpos+i] = false;
                    }    
                }    
                //set W
                uint8_t curW = getW(G->kmers[i]);
                if (curW == 127) {
                    char* curseq = kmertos(G->kmers[i], G->alphabet, G->alph_size);
                    std::cerr << "Failure decoding kmer " << i << "\n" << G->kmers[i] << "\n" << curseq << "\n";
                    free(curseq);
                    exit(1);
                }    
                if (!curW && curpos+i)
                    G->p=curpos+i;
                if (i) {
                    for (size_t j=i-1;compare_kmer_suffix(G->kmers[j], G->kmers[i], 1);--j) {
                        //TODO: recalculating W is probably faster than doing a pragma for ordered
                        if (getW(G->kmers[j]) == curW) {
                            curW += G->alph_size;
                            break;
                        }    
                        if (!j) 
                            break;
                    }    
                }    
                G->W_stat[curpos+i] = curW;
            }    
        }    
        for (size_t i=0;i<G->kmers.size();++i) {
            char cF=getPos(G->kmers[i], G->k-1, G->alphabet, G->alph_size);
            if (cF != G->alphabet[G->lastlet]) {
                for ((G->lastlet)++; G->lastlet<G->alph_size; (G->lastlet)++) {
                    G->F[G->lastlet]=curpos+i-1;
                    if (G->alphabet[G->lastlet]==cF) {
                        break;
                    }    
                }    
            }    
        }    
        G->kmers.clear();
    }



    /** This function takes a character c and appends it to the end of the graph sequence
     * given that the corresponding note is not part of the graph yet.
     */
    void append_pos(DBG_succ* G, TAlphabet c) {

        // check that the last position of the graph is indeed a terminal
        assert((*(G->W))[G->p] == 0);
        TAlphabet c_p = G->get_node_end_value(G->p);
        // get range of identical nodes (without W) pos current end position
        std::pair<uint64_t, uint64_t> R = G->get_equal_node_range(G->p);
        //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

        // get position of first occurence of c in W after p
        uint64_t next_c = G->succ_W(G->p, c);
        // check if c is part of range
        bool exist_c = (next_c <= R.second);
        if (!exist_c) {
            // get position of first occurence of c- in W after p
            next_c = G->succ_W(G->p, c + G->alph_size);
            // check if c- is part of range
            exist_c = (next_c <= R.second);
        }

        /**
         * if the character already exists in the range, we delete the terminal symbol
         * at p, insert c at fwd(next_c) and update p.
         */
        if (exist_c) {
            uint64_t p_new = G->fwd(next_c);
            // remove old terminal symbol
            G->last->deleteBit(G->p);
            G->W->remove(G->p);
            // adapt position if altered by previous deletion
            p_new -= (G->p < p_new);
            // insert new terminal symbol 
            // we have to insert 0 into last as the node already existed in the range 
            // and the terminal symbol is always first
            G->last->insertBit(p_new, false);
            G->W->insert(0, p_new);
            // update new terminal position
            G->p = p_new;
            // take care of updating the offset array F
            G->update_F(c_p, false);
            //assert(get_node_end_value(p) == c);
            G->update_F(c, true);
        } else {
            /**
             * We found that c does not yet exist in the current range and now have to
             * figure out if we need to add c or c- to the range.
             * To do this, we check if there is a previous position j1 with W[j1] == c
             * whose node shares a k-1 suffix with the current node. If yes, we add c- 
             * instead of c.
             */
            // get position of last occurence of c before p (including p - 1)
            uint64_t last_c = G->pred_W(G->p - 1, c);
            // if this position exists
            if (last_c > 0) {
                uint64_t x = G->fwd(last_c);
                assert((*(G->last))[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

                // check, if there are any c or c- symbols following after position p
                uint64_t next_c = G->succ_W(G->p + 1, c);
                uint64_t next_cm = G->succ_W(G->p + 1, c + G->alph_size);
                // there is no c between p and next_cm and next_cm is a c- ==> we should add a c- 
                // all nodes between W[i] = c and W[j] = c- share a common suffix of length k-1
                bool minus1 = (next_cm < next_c);
                // check, if we share a k-1 suffix with last_c
                if (!minus1) {
                    minus1 = G->compare_node_suffix(G->p, last_c);
                }

                // adding a new node can influence following nodes that share a k-1 suffix with the
                // new node -> need to adapt the respektive cc to a cc-
                bool minus2 = false;
                if (next_c < G->W->n) {
                    minus2 = G->compare_node_suffix(G->p, next_c);
                    if (minus2) {
                        G->replaceW(next_c, (*(G->W))[next_c] + G->alph_size);
                    }
                }

                G->replaceW(G->p, minus1 ? c + G->alph_size : c);
                // after we are done, assert that the order within the range we created 
                // is still valid within W
                if (G->p - R.second > 0) {
                    G->sort_W_locally(G->p, R.second);
                }

                // if minus1 is true, we share a k-1 suffix with the node at 
                // last_c and thus need to adapt our insertion position by -1, 
                // as we would like to insert before it. Otherwise we insert directly after
                // it as we are now sorted after it. 
                if (minus1) {
                    G->p = x;
                    G->last->insertBit(x, false);
                    G->W->insert(0, x);
                } else if (minus2) {
                    G->p = x + 1;
                    G->last->insertBit(x + 1, false);
                    G->W->insert(0, x + 1);
                // no node shares a k-1 suffix with last_c and thus the new node comes after
                // the forward of last_c (as the current node came after last_c as well)
                } else {
                    G->p = x + 1;
                    G->last->insertBit(x + 1, true);
                    G->W->insert(0, x + 1);
                }
            } else {
                uint64_t x = G->F[c] + 1;
                uint64_t next_c = G->succ_W(G->p + 1, c);
                bool minus = false;
                if (next_c < G->W->n) {
                    minus = G->compare_node_suffix(G->p, next_c);
                }
                G->replaceW(G->p, c);
                if (G->p - R.second > 0) {
                    G->sort_W_locally(G->p, R.second);
                }
                G->p = x;
                if (minus) {
                    G->replaceW(next_c, (*(G->W))[next_c] + G->alph_size);
                    G->last->insertBit(x, false);
                } else {
                    G->last->insertBit(x, true);
                }
                G->W->insert(0, x);
            }
            G->update_F(c, true);
        }
        // update sorting at new location of p
        // with this we assert that $ is always inserted at the first position 
        // of a range of equal nodes --> this will help us to prevent multiple insertions
        // of already existing nodes
        R = G->get_equal_node_range(G->p);
        if (R.second - R.first > 0) {
            G->sort_W_locally(R.first, R.second);
            while ((*(G->W))[G->p] != 0)
                G->p--;
            assert((*(G->W))[G->p] == 0);
        }
    }


    /** This function takes a pointer to a graph structure and concatenates the arrays W, last 
     * and F to this graph's arrays. In almost all cases this will not produce a valid graph and 
     * should only be used as a helper in the parallel merge procedure.
     */
    void append_graph(DBG_succ* G_s, DBG_succ* G_t) {

        size_t curr_pos = G_t->get_size();

        if (G_t->config->verbose)
            std::cout << "    adding " << G_s->get_size() << " edges" << std::endl;
        // handle last and W
        for (size_t j = 1; j < G_s->get_size(); ++j) {
            G_t->last->insertBit(curr_pos, G_s->get_last(j));
            G_t->W->insert(G_s->get_W(j), curr_pos);
            ++curr_pos;
        }
        if (G_t->config->verbose)
            std::cout << "new total edges: " << G_s->W->n << std::endl;

        // handle F
        assert(G_t->F.size() == G_s->F.size());
        for (size_t j = 0; j < G_t->F.size(); ++j) {
            G_t->F.at(j) += G_s->F.at(j);
        }
    }


    /** 
     * This function takes a pointer to a graph structure and concatenates the arrays W, last 
     * and F to this graph's static containers last_stat and W_stat. In almost all cases 
     * this will not produce a valid graph and should only be used as a helper in the 
     * parallel merge procedure.
     */
    void append_graph_static(DBG_succ* G_t, DBG_succ* G_s) {

        size_t n = G_s->get_size();
        if (G_t->config->verbose)
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
                    offsets.at(o) += !(*(G_s->W->R))[pos];
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

        //std::cerr << "R size: " << G_s->W->R->size() << std::endl;

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
                bit = (*(G_s->W->R))[ib * n + p + co];
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
            bit = (*(G_s->W->R))[(b - 1) * n + p + co];
            if (bit) {
                v |= m;
            }
            if (i == 0)
                continue;
            G_t->W_stat.push_back(v);
            G_t->last_stat.push_back(G_s->get_last(i));
        }

        // handle F
        assert(G_t->F.size() == G_s->F.size());
        for (size_t j = 0; j < G_t->F.size(); ++j) {
            G_t->F.at(j) += G_s->F.at(j);
        }
    }
}
