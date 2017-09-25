#include "construct.hpp"

#include <unordered_set>
#include <unordered_map>


#include "kseq.h"
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

    bool check_suffix(DBG_succ* G, char* target, std::string& suffix, const char* nt_lookup) {
        std::string cursuff = std::string(target+(G->k)-suffix.length(), target+G->k);

        for (std::string::iterator it = cursuff.begin(); it!=cursuff.end(); ++it)
            *it = G->alphabet[(uint8_t)nt_lookup[(uint8_t)*it]];
        return cursuff == suffix;
    }


    void add_seq_fast(DBG_succ* G, kstring_t &seq, kstring_t &name, bool add_bridge, unsigned int parallel, std::string suffix, size_t cid) {

        if (debug) {
            G->print_seq();
            G->print_state();
            std::cout << "======================================" << std::endl;
        }

        std::string label_str = std::string(name.s);
        uint32_t label_id,max_label_id=0;
        std::unordered_map<std::string, uint32_t>::iterator id_it = G->label_to_id_map.find(label_str);
        if (id_it == G->label_to_id_map.end()) {
            label_id = (uint32_t) G->id_to_label.size();
            G->id_to_label.push_back(label_str);
            G->label_to_id_map[label_str] = label_id;
            max_label_id = std::max(max_label_id, label_id);
        } else {
            label_id = id_it->second;
        }
        if (max_label_id >= G->annotation_full.size()) {
            G->annotation_full.resize(max_label_id+1);
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
		char *bridge = (char*)malloc(G->k+2);
        memset(bridge, 'X', G->k); 
        if (!seq.l)
            return;
        bridge[G->k] = seq.s[0];
        bridge[G->k+1] = 0;
        size_t i=0;
        //std::cout << "Loading next sequence with " << parallel << " threads\n";
        if (add_bridge) {
            for (i=0;i<std::min(G->k,seq.l); ++i) {
                if (check_suffix(G, bridge, suffix, nt_lookup)) {
                    //add 1 to the CID if the kmer is a bridge
                    //assumes that CID is an even number
                    //kmers.push_back(std::make_pair(stokmer(bridge, k+1, nt_lookup), cid+1));
                    G->kmers.push_back(kmer_boost::KMer{kmer_boost::stokmer(bridge, G->k+1, nt_lookup), cid+1, 0});
                }   
                memmove(bridge, bridge+1, G->k); 
                bridge[G->k]= (i+1 < seq.l) ? seq.s[i+1] : 'X';
            }   
        }   
        if (G->k<seq.l) {
            #pragma omp parallel num_threads(parallel)
            {   
                std::vector<kmer_boost::KMer> kmer_priv;
                //std::vector<std::pair<ui256, size_t> > kmer_priv;
                #pragma omp for nowait
                for (i=0; i < seq.l-G->k; ++i) {
                    if (check_suffix(G, seq.s+i, suffix, nt_lookup)) {
                        kmer_priv.push_back(kmer_boost::KMer{kmer_boost::stokmer(seq.s+i, G->k+1, nt_lookup), cid, label_id});
                        //kmer_priv.push_back(std::make_pair(stokmer(seq+i, k+1, nt_lookup), cid));
                    }   
                }   
                #pragma omp critical
                G->kmers.insert(G->kmers.end(), std::make_move_iterator(kmer_priv.begin()), std::make_move_iterator(kmer_priv.end()));
            }   
            memcpy(bridge, seq.s+seq.l-G->k, G->k); 
            bridge[G->k]='X';
        }   
        if (add_bridge) {
            for (i=0;i<G->k;++i) {
                if (check_suffix(G, bridge, suffix, nt_lookup)) {
                    //add 1 to the CID if the kmer is a bridge
                    G->kmers.push_back(kmer_boost::KMer{kmer_boost::stokmer(bridge, G->k+1, nt_lookup), cid+1, 0});
                    //kmers.push_back(std::make_pair(stokmer(bridge, k+1, nt_lookup), cid+1));
                }   
                memmove(bridge, bridge+1, G->k); 
                bridge[G->k]='X';
            }   
        }   
        free(bridge); 
        //free(nt_lookup);
    }

    //removes duplicate kmers, counts coverage, and assigns labels
    std::vector<kmer_boost::KMer>::iterator unique_count(std::vector<kmer_boost::KMer>::iterator first, std::vector<kmer_boost::KMer>::iterator last) {
    //std::vector<std::pair<ui256, size_t> >::iterator unique_count(std::vector<std::pair<ui256, size_t> >::iterator first, std::vector<std::pair<ui256, size_t> >::iterator last) {
        //Adapted from http://en.cppreference.com/w/cpp/algorithm/unique
        std::unordered_set<size_t> cids;
        if (first == last)
            return last;
        std::vector<kmer_boost::KMer>::iterator result = first, start = first;
        cids.insert(result->second);
        cpp_int annot;
        if (result->annot != 0) {
            annot = (cpp_int(1) << static_cast<size_t>(result->annot-1));
        } else {
            annot = 0;
        }
        while (++first != last) {
            if (!(result->first == first->first)) {
                //set a specific code for the coverage
                if (cids.find(0) != cids.end()) {
                    //0 if in reference and not a bridge
                    result->second = 0;
                } else if (cids.find(1) != cids.end()) {
                    //1 if in reference and not a bridge
                    result->second = 1;
                } else {
                    //2*coverage + is_bridge
                    result->second = cids.size()*2 + (*cids.begin() % 2);
                }
                result->second = cids.size();
                result->annot = annot;
                if (++result != first) {
                    *result = std::move(*first);
                }
                cids.clear();
                if (result->annot != 0) {
                    annot = (cpp_int(1) << static_cast<size_t>(result->annot-1));
                } else {
                    annot = 0;
                }
            }
            if (first->annot) {
                //num_annots = std::max(num_annots, static_cast<size_t>(first->annot));
                //TODO: add check for first->annot being too big
                annot |= (cpp_int(1) << static_cast<size_t>(first->annot-1));
            }
            cids.insert(first->second);
        }
        result->second = cids.size();
        result->annot = annot;
        return ++result;
    }


    void construct_succ(DBG_succ* G, unsigned int parallel) {
        omp_set_num_threads(std::max((int)parallel,1));
        __gnu_parallel::sort(G->kmers.begin(),G->kmers.end());
        std::vector<kmer_boost::KMer>::iterator uniq_count = unique_count(G->kmers.begin(), G->kmers.end());
        G->kmers.erase(uniq_count, G->kmers.end() );
        std::vector<std::vector<uint8_t> > annots(G->annotation_full.size(), std::vector<uint8_t>(G->kmers.size()));
        //TODO: set the bit vector, then store it
        //TODO: when merging into G, keep track of the fact that this may have a different number of annotations that the rest of hte graph

        //DEBUG: output kmers in current bin
        /*
        std::cerr << "\n";
        for (size_t i=0;i<G->kmers.size();++i) {
            char* curseq = kmertos(G->kmers[i].first, G->alphabet, G->alph_size);
            std::cerr << G->kmers[i].first << "\t" << curseq+1 << " " << curseq[0] << " " << G->kmers[i].second << "\n";
            free(curseq);
        }
        */

        size_t curpos = G->W_stat.size();
        G->W_stat.resize(G->W_stat.size()+G->kmers.size());
        G->last_stat_safe.resize(G->last_stat_safe.size()+G->kmers.size(), true);
        G->coverage.resize(G->coverage.size()+G->kmers.size(),0);
        G->bridge_stat.resize(G->bridge_stat.size()+G->kmers.size(),false);
        
        #pragma omp parallel num_threads(parallel)
        {    
            #pragma omp for nowait
            for (size_t i=0;i<G->kmers.size();++i) {
                //set last
                if (i+1 < G->kmers.size()) {
                    bool dup = kmer_boost::compare_kmer_suffix(G->kmers[i].first, G->kmers[i+1].first);
                    if (dup) {
                        G->last_stat_safe[curpos+i] = false;
                    }    
                }
                //set bridge and coverage if from a read
                if (G->kmers[i].second > 1) {
                    G->coverage[curpos+i] = G->kmers[i].second / 2;
                    if (G->kmers[i].second % 2) {
                        G->bridge_stat[curpos+i] = true;
                    }
                }
                //set W
                uint8_t curW = kmer_boost::getW(G->kmers[i].first);
                if (curW == 127) {
                    char* curseq = kmer_boost::kmertos(G->kmers[i].first, G->alphabet, G->alph_size);
                    std::cerr << "Failure decoding kmer " << i << "\n" << G->kmers[i].first << "\n" << curseq << "\n";
                    free(curseq);
                    exit(1);
                }    
                if (!curW && curpos+i)
                    G->p=curpos+i;
                if (i) {
                    for (size_t j=i-1;kmer_boost::compare_kmer_suffix(G->kmers[j].first, G->kmers[i].first, 1);--j) {
                        //TODO: recalculating W is probably faster than doing a pragma for ordered
                        if (kmer_boost::getW(G->kmers[j].first) == curW) {
                            curW += G->alph_size;
                            break;
                        }    
                        if (!j) 
                            break;
                    }    
                }    
                G->W_stat[curpos+i] = curW;
                //set  bitvector
                for (size_t j=0;(G->kmers[i].annot) >> j; ++j) {
                    assert(j < annots.size());
                    annots[j][i] = (uint8_t)(((G->kmers[i].annot) >> j) % 2);
                }
            }    
        }    
        for (size_t i=0;i<G->kmers.size();++i) {
            char cF=kmer_boost::getPos(G->kmers[i].first, G->k-1, G->alphabet, G->alph_size);
            if (cF != G->alphabet[G->lastlet]) {
                for ((G->lastlet)++; G->lastlet<G->alph_size; (G->lastlet)++) {
                    G->F[G->lastlet]=curpos+i-1;
                    if (G->alphabet[G->lastlet]==cF) {
                        break;
                    }    
                }    
            }    
        }    
        sdsl::bit_vector bv(G->W_stat.size());
        for (size_t i=0;i<annots.size();++i) {
            size_t j=0;
            for (;i < G->annotation_full.size() && G->annotation_full.at(i) != NULL && j<G->annotation_full.at(i)->size();++j) {
                bv[j] = G->annotation_full.at(i)->operator[](j);
            }
            for (;j<curpos;++j)
                bv[j]=0;
            for (size_t k=0;j<bv.size();++j,++k) {
                bv[j] = annots.at(i).at(k);
            }
            if (G->annotation_full.at(i) != NULL)
                delete G->annotation_full[i];
            G->annotation_full[i] = new sdsl::sd_vector<>(bv);
        }
        G->kmers.clear();
    }

    void clean_bridges(DBG_succ* G) {
        //TODO: MERGE FROM SEQMERGE
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
