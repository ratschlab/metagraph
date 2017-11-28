#include "dbg_succinct.hpp"

#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <cstdio>
#include <parallel/algorithm>


// add a full sequence to the graph
void DBG_succ::add_seq(kstring_t &seq, bool append) {

#ifdef DBGDEBUG
    print_seq();
    print_state();
    std::cout << "======================================" << std::endl;
#endif

    // Padding of the input genome / read
    if (W->size() == 2) {
        for (size_t j = 0; j < k; j++) {
            append_pos(6);
#ifdef DBGDEBUG
            print_seq();
            print_state();
            std::cout << "======================================" << std::endl;
#endif
        }
    }

    /** Iterate over input sequence and enumerae all
     * k-mers.
     */

    size_t i;
    uint64_t *ckmer = new uint64_t[k+1];
    for (i = 0; i < std::min(seq.l, k); ++i) {
        ckmer[i] = 6;
    }
    ckmer[i] = 0;
    uint64_t c;
    uint64_t ind;
    if (!append) {
        ind = index(ckmer, i);
        if (!ind) {
            ind = p;
            i = 0;
        } else {
            if (k >= seq.l) {
                delete[] ckmer;
                return;
            }
            i = k;
        }
    } else {
        i = 0;
        ind = 0;
    }
    for (; i < seq.l; ++i) {
        if (i > 0 && i % 1'000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10'000 == 0) {
                fprintf(stdout, "%lu - edges %" PRIu64 " / nodes %" PRIu64 "\n",
                                i, get_edge_count(), get_node_count());
            }
        }

        c = get_alphabet_number(seq.s[i]);
        ind = append_pos(c, ckmer, ind) * !append;
        memmove(ckmer, &ckmer[1], sizeof(uint64_t)*(k-1));
        ckmer[k-1] = c;

#ifdef DBGDEBUG
        print_seq();
        print_state();
        std::cout << "======================================" << std::endl;
#endif
    }
    assert(!append || W->operator[](ind)==0);
    // Padding after sequence to get back into default state.
    if (W->operator[](ind)==0) {
        for (size_t j = 0; j < k; j++) {
            ind = append_pos(6, ckmer, ind) * !append;
            memmove(ckmer, &ckmer[1], sizeof(uint64_t)*(k-1));
            ckmer[k-1] = 6;
#ifdef DBGDEBUG
            print_seq();
            print_state();
            std::cout << "======================================" << std::endl;
#endif
        }
        if (!append)
            ind = append_pos(0, ckmer, ind);
    }

    fprintf(stdout, "edges %" PRIu64 " / nodes %" PRIu64 "\n",
                    get_edge_count(), get_node_count());
    delete[] ckmer;
}


bool check_suffix(DBG_succ *G, const char *target, std::string& suffix) {

    std::string cursuff = std::string(target + G->k - suffix.length(), target + G->k);

    for (auto it = cursuff.begin(); it != cursuff.end(); ++it) {
        *it = G->alphabet[static_cast<uint8_t>(DBG_succ::get_alphabet_number(*it))];
    }
    return cursuff == suffix;
}


void DBG_succ::add_seq_fast(const std::string &seq,
                            bool add_bridge, unsigned int parallel,
                            std::string suffix) {

#ifdef DBGDEBUG
        print_seq();
        print_state();
        std::cout << "======================================" << std::endl;
#endif

    std::vector<uint32_t> label_id;

    // ther is nothing to parse
    if (!seq.size()) {
        return;
    }

	char *bridge = (char*) malloc(k+2);
    memset(bridge, 'X', k);
    bridge[k] = seq[0];
    bridge[k+1] = 0;

    size_t i = 0;
    //std::cout << "Loading next sequence with " << parallel << " threads\n";
    if (add_bridge) {
        for (i = 0; i < std::min(k, seq.length()); ++i) {
            if (check_suffix(this, bridge, suffix)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k + 1), DBG_succ::get_alphabet_number));
            }
            memmove(bridge, bridge+1, k);
            bridge[k] = (i+1 < seq.length()) ? seq[i+1] : 'X';
        }
    }
    if (k < seq.length()) {
        #pragma omp parallel num_threads(parallel)
        {
            std::vector<KMer> kmer_priv;
            #pragma omp for nowait
            for (i = 0; i < seq.length() - k; ++i) {
                if (check_suffix(this, seq.c_str() + i, suffix)) {
                    kmer_priv.push_back(KMer::from_string(
                        std::string(seq.c_str() + i, k + 1),
                        DBG_succ::get_alphabet_number
                    ));
                }
            }
            #pragma omp critical
            kmers.insert(kmers.end(), std::make_move_iterator(kmer_priv.begin()), std::make_move_iterator(kmer_priv.end()));
        }
        memcpy(bridge, seq.c_str() + seq.length() - k, k);
        bridge[k]='X';
    }
    if (add_bridge) {
        for (i = 0; i < k; ++i) {
            if (check_suffix(this, bridge, suffix)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k + 1), DBG_succ::get_alphabet_number));
            }
            memmove(bridge, bridge+1, k);
            bridge[k] = 'X';
        }
    }
    free(bridge);
}

void DBG_succ::construct_succ(unsigned int parallel) {

    // parallel sort of all kmers
    omp_set_num_threads(std::max((int)parallel,1));
    __gnu_parallel::sort(this->kmers.begin(),this->kmers.end());

    auto last = std::unique(this->kmers.begin(), this->kmers.end());
    this->kmers.erase(last, this->kmers.end()); 

    //DEBUG: output kmers in current bin
    /*
    std::cerr << "\n";
    for (size_t i=0;i<this->kmers.size();++i) {
        char* curseq = kmer_to_s(this->kmers[i].first, this->alphabet, this->alph_size);
        std::cerr << this->kmers[i].first << "\t" << curseq+1 << " " << curseq[0] << " " << this->kmers[i].second << "\n";
        free(curseq);
    }
    */

    size_t curpos = this->W_stat.size();
    this->W_stat.resize(this->W_stat.size() + this->kmers.size());
    this->last_stat_safe.resize(this->last_stat_safe.size() + this->kmers.size(), true);
    this->coverage.resize(this->coverage.size() + this->kmers.size(),0);
    //this->bridge_stat.resize(this->bridge_stat.size() + this->kmers.size(), false);

    #pragma omp parallel num_threads(parallel)
    {
        #pragma omp for nowait
        for (size_t i = 0; i < this->kmers.size(); ++i) {
            //set last
            if (i + 1 < this->kmers.size()) {
                bool dup = KMer::compare_kmer_suffix(this->kmers[i], this->kmers[i + 1]);
                if (dup) {
                    this->last_stat_safe[curpos + i] = false;
                }
            }
            //set W
            uint8_t curW = this->kmers[i][0];
            if (curW == 127) {
                std::string curseq = this->kmers[i].to_string(this->alphabet);
                std::cerr << "Failure decoding kmer " << i << "\n" << this->kmers[i] << "\n" << curseq << "\n";
                exit(1);
            }
            if (!curW && curpos+i)
                this->p = curpos + i;
            if (i) {
                for (size_t j = i - 1; KMer::compare_kmer_suffix(this->kmers[j], this->kmers[i], 1); --j) {
                    //TODO: recalculating W is probably faster than doing a pragma for ordered
                    if (this->kmers[j][0] == curW) {
                        curW += this->alph_size;
                        break;
                    }
                    if (!j)
                        break;
                }
            }
            this->W_stat[curpos+i] = curW;
        }
    }
    for (size_t i = 0; i < this->kmers.size(); ++i) {
        char cF = this->alphabet[this->kmers[i][this->k]];
        if (cF != this->alphabet[this->lastlet]) {
            for ((this->lastlet)++; this->lastlet<this->alph_size; (this->lastlet)++) {
                this->F[this->lastlet] = curpos + i - 1;
                if (this->alphabet[this->lastlet] == cF) {
                    break;
                }
            }
        }
    }
    this->kmers.clear();
}

/** This function takes a character c and appends it to the end of the graph sequence
 * given that the corresponding note is not part of the graph yet.
 */
uint64_t DBG_succ::append_pos(uint64_t c, uint64_t *ckmer, uint64_t i) {

    // check that the last position of the graph is indeed a terminal
    assert((*(this->W))[this->p] == 0);
    uint64_t *p;
    bool append = false;
    uint64_t c_p = this->get_node_end_value(this->p);
    // get range of identical nodes (without W) pos current end position
    std::pair<uint64_t, uint64_t> R;
    //std::pair<uint64_t, uint64_t> R = this->get_equal_node_range(this->p);
    //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

    if (!i) {
        p = &(this->p);
        append = true;
        R.second = this->succ_last(*p);
    } else {
        append = false;
        R = this->get_equal_node_range(this->p);
        p = &(R.first);
    }

    // get position of first occurence of c in W after p
    uint64_t next_c = this->succ_W(*p, c);
    // check if c is part of range
    bool exist_c = (next_c <= R.second);
    if (!exist_c) {
        // get position of first occurence of c- in W after p
        next_c = this->succ_W(*p, c + this->alph_size);
        // check if c- is part of range
        exist_c = (next_c <= R.second);
    }

    /**
     * if the character already exists in the range, we delete the terminal symbol
     * at p, insert c at fwd(next_c) and update p.
     */
    if (exist_c) {
        if (!append) {
            //if not appending, then this has already been observed and we only need to delete the extra pointer
            for (i=0;this->W->operator[](*p)==0;++i,++(*p)) {
                if (*p==(this->p))
                    continue;
                if (this->debug) {
                    assert(i<2);
                    assert(*p < R.second);
                }
                this->W->remove(*p);
                this->last->deleteBit(*p);
                this->update_F(c_p, false);
                if (*p <= this->p) {
                    (this->p)--;
                    assert(this->W->operator[](*p) == 0);
                }
                if (*p <= next_c)
                    next_c--;
            }
            return this->fwd(next_c);
        }
        uint64_t p_new = this->fwd(next_c);
        // remove old terminal symbol
        this->last->deleteBit(*p);
        this->W->remove(*p);
        // adapt position if altered by previous deletion
        p_new -= (*p < p_new);
        // insert new terminal symbol
        // we have to insert 0 into last as the node already existed in the range
        // and the terminal symbol is always first
        this->last->insertBit(p_new, false);
        this->W->insert(0, p_new);
        // update new terminal position
        *p = p_new;
        // take care of updating the offset array F
        this->update_F(c_p, false);
        //assert(get_node_end_value(p) == c);
        this->update_F(c, true);
    } else {
        /**
         * We found that c does not yet exist in the current range and now have to
         * figure out if we need to add c or c- to the range.
         * To do this, we check if there is a previous position j1 with W[j1] == c
         * whose node shares a k-1 suffix with the current node. If yes, we add c-
         * instead of c.
         */
        if (!append) {
            //we need to insert a new pointer
            if (*p == this->p)
                (*p)++;
            if (this->W->operator[](*p)) {
                //if no placeholder exists
                this->W->insert(0, *p);
                this->update_F(c_p, true);
                this->last->insertBit(*p, false);
                R.second++;
                if (*p <= this->p) {
                    (this->p)++;
                    assert(this->W->operator[](this->p) == 0);
                }
            }
        }
        // get position of last occurence of c before p (including p - 1)
        uint64_t last_c = this->pred_W(*p - 1, c);
        // if this position exists
        if (last_c > 0) {
            uint64_t x = this->fwd(last_c);
            assert((*(this->last))[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

            // check, if there are any c or c- symbols following after position p
            uint64_t next_c = this->succ_W(*p + 1, c);
            uint64_t next_cm = this->succ_W(*p + 1, c + this->alph_size);
            // there is no c between p and next_cm and next_cm is a c- ==> we should add a c-
            // all nodes between W[i] = c and W[j] = c- share a common suffix of length k-1
            bool minus1 = (next_cm < next_c);
            // check, if we share a k-1 suffix with last_c
            if (!minus1) {
                if (ckmer) {
                    minus1 = this->compare_node_suffix(ckmer, last_c);
                } else {
                    minus1 = this->compare_node_suffix(*p, last_c);
                }
            }

            // adding a new node can influence following nodes that share a k-1 suffix with the
            // new node -> need to adapt the respektive cc to a cc-
            bool minus2 = false;
            if (next_c < this->W->size()) {
                if (ckmer) {
                    minus2 = this->compare_node_suffix(ckmer, next_c);
                } else {
                    minus2 = this->compare_node_suffix(*p, next_c);
                }
                if (minus2) {
                    this->replaceW(next_c, (*(this->W))[next_c] + this->alph_size);
                }
            }

            this->replaceW(*p, minus1 ? c + this->alph_size : c);
            // after we are done, assert that the order within the range we created
            // is still valid within W
            if (*p - R.second > 0) {
                this->sort_W_locally(*p, R.second);
            }

            // if minus1 is true, we share a k-1 suffix with the node at
            // last_c and thus need to adapt our insertion position by -1,
            // as we would like to insert before it. Otherwise we insert directly after
            // it as we are now sorted after it.
            if (minus1) {
                *p = x;
                this->last->insertBit(x, false);
                this->W->insert(0, x);
            } else if (minus2) {
                *p = x + 1;
                this->last->insertBit(x + 1, false);
                this->W->insert(0, x + 1);
            // no node shares a k-1 suffix with last_c and thus the new node comes after
            // the forward of last_c (as the current node came after last_c as well)
            } else {
                *p = x + 1;
                this->last->insertBit(x + 1, true);
                this->W->insert(0, x + 1);
            }
        } else {
            uint64_t x = this->F[c] + 1;
            uint64_t next_c = this->succ_W(*p + 1, c);
            bool minus = false;
            if (next_c < this->W->size()) {
                if (ckmer) {
                    minus = this->compare_node_suffix(ckmer, next_c);
                } else {
                    minus = this->compare_node_suffix(*p, next_c);
                }
            }
            this->replaceW(*p, c);
            if (*p - R.second > 0) {
                this->sort_W_locally(*p, R.second);
            }
            *p = x;
            if (minus) {
                this->replaceW(next_c, (*(this->W))[next_c] + this->alph_size);
                this->last->insertBit(x, false);
            } else {
                this->last->insertBit(x, true);
            }
            this->W->insert(0, x);
        }
        if (*p < this->p || (!append && *p == this->p)) {
            (this->p)++;
            assert(this->W->operator[](this->p) == 0);
        }
        this->update_F(c, true);
    }
    // update sorting at new location of p
    // with this we assert that $ is always inserted at the first position
    // of a range of equal nodes --> this will help us to prevent multiple insertions
    // of already existing nodes
    R = this->get_equal_node_range(*p);
    if (R.second - R.first > 0) {
        this->sort_W_locally(R.first, R.second);
        *p = R.first;
        if (!append && *p == this->p)
            (*p)++;
        /*
        while ((*(this->W))[*p] != 0)
            (*p)--;
        */
        assert(this->W->operator[](*p)==0);
    }
    return *p;
}

/** This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's arrays. In almost all cases this will not produce a valid graph and
 * should only be used as a helper in the parallel merge procedure.
 */
void DBG_succ::append_graph(DBG_succ *G, bool verbose) {

    size_t curr_pos = this->get_size();

    if (verbose)
        std::cout << "    adding " << G->get_size() << " edges" << std::endl;
    // handle last and W
    for (size_t j = 1; j < G->get_size(); ++j) {
        this->last->insertBit(curr_pos, G->get_last(j));
        this->W->insert(G->get_W(j), curr_pos);
        ++curr_pos;
    }
    if (verbose)
        std::cout << "new total edges: " << G->W->size() << std::endl;

    // handle F
    assert(this->F.size() == G->F.size());
    for (size_t j = 0; j < this->F.size(); ++j) {
        this->F.at(j) += G->F.at(j);
    }
}

/**
 * This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's static containers last_stat and W_stat. In almost all cases
 * this will not produce a valid graph and should only be used as a helper in the
 * parallel merge procedure.
 */
void DBG_succ::append_graph_static(DBG_succ *G, bool verbose) {

    size_t n = G->get_size();
    if (verbose)
        std::cout << "    adding " << n << " edges" << std::endl;

    //size_t n_old = this->last_stat.size();
    //this->last_stat.resize(n_old + n);
    //this->W_stat.resize(n_old + n);

    const size_t b = 4;
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
                offsets.at(o) += !G->W->get_bit_raw(pos);
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

    //std::cerr << "R size: " << G->W->R->size() << std::endl;

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
            bit = G->W->get_bit_raw(ib * n + p + co);
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
        bit = G->W->get_bit_raw((b - 1) * n + p + co);
        if (bit) {
            v |= m;
        }
        if (i == 0)
            continue;
        this->W_stat.push_back(v);
        this->last_stat.push_back(G->get_last(i));
    }

    // handle F
    assert(this->F.size() == G->F.size());
    for (size_t j = 0; j < this->F.size(); ++j) {
        this->F.at(j) += G->F.at(j);
    }
}

// Given an edge list, remove them from the graph.
// If a ref_point is given, it's updated index is return, otherwise return 0
uint64_t DBG_succ::remove_edges(std::set<uint64_t> &edges, uint64_t ref_point) {
    uint64_t offset = 0;
    uint64_t i, j, d;

    for (auto x = edges.begin(); x != edges.end(); ++x) {
        assert(*x >= offset);
        d = this->W->operator[](*x - offset);
        if (d < this->alph_size) {
            //fix W array
            j = this->succ_W(*x-offset + 1, d);
            for (i = *x - offset + 1; i < j; ++i) {
                if (this->W->operator[](i) == d + this->alph_size) {
                    this->replaceW(i, d);
                    break;
                }
            }
        }
        this->W->remove(*x - offset);
        this->update_F(this->get_node_end_value(*x - offset), false);
        //if the current node has multiple outgoing edges, remove one of the 0s from
        //last instead of 1
        if (this->get_last(*x-offset) && (*x >= offset + 1) && !this->get_last(*x - offset - 1)) {
            this->last->deleteBit(*x - offset - 1);
        } else {
            this->last->deleteBit(*x - offset);
        }
        //fix pointers
        if (ref_point && *x-offset <= ref_point)
            ref_point--;
        if (*x-offset <= this->p) {
            (this->p)--;
            assert(this->W->operator[](this->p) == 0);
        }
        offset++;
    }
    return ref_point;
}

//Given a graph and a minimum number of splits, generate a list of suffices from the alphabet
std::deque<std::string> DBG_succ::generate_suffices(unsigned int nsplits) {
    unsigned int suffix_len = (unsigned int) ceil(log2(nsplits) / log2(alph_size - 1));

    //should be set to at most k-1 so that W calculation is correct
    suffix_len = std::min(suffix_len, (unsigned int) k - 1);
    std::deque<std::string> suffices = {""};
    for (size_t i = 0; i < suffix_len; ++i) {
         while (suffices[0].length() < suffix_len) {
             for (size_t j = 0; j < alph_size; ++j) {
                  suffices.push_back(alphabet[j] + suffices[0]);
             }
             suffices.pop_front();
         }
    }
    assert(suffices.size() == pow(alph_size, suffix_len));
    return suffices;
}

void DBG_succ::add_sink(unsigned int parallel, std::string suffix) {
    add_seq_fast(std::string(start.s, start.l), false, parallel, suffix);
    add_seq_fast(std::string(graphsink.s, graphsink.l), true, parallel, suffix);
}
