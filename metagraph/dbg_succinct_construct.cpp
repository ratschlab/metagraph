#include "dbg_succinct.hpp"

#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <cstdio>
#include <parallel/algorithm>


// add a full sequence to the graph
void DBG_succ::add_sequence(const std::string &seq, bool append) {

    // Padding of the input genome / read
    if (W->size() == 2) {
        for (size_t j = 0; j < k_; j++) {
            append_pos(6);
            // append_pos(encode('X'));
        }
    }

    // Iterate over input sequence and enumerae all k-mers.
    size_t i;
    std::vector<uint64_t> ckmer(k_ + 1, encode('X'));
    for (i = 0; i < std::min(seq.length(), k_); ++i) {
        ckmer[i] = 6;
    }
    ckmer[i] = 0;

    uint64_t ind;
    if (!append) {
        ind = index(ckmer, i);
        if (!ind) {
            ind = p_;
            i = 0;
        } else {
            if (k_ >= seq.length())
                return;
            i = k_;
        }
    } else {
        i = 0;
        ind = 0;
    }
    for (; i < seq.length(); ++i) {
        if (i > 0 && i % 1'000 == 0) {
            std::cout << "." << std::flush;
            if (i % 10'000 == 0) {
                verbose_cout(i, " - edges ", num_edges(), " / nodes ", num_nodes(), "\n");
            }
        }

        uint64_t c = encode(seq[i]);
        ind = append_pos(c, ckmer.data(), ind) * !append;
        memmove(ckmer.data(), &ckmer[1], sizeof(uint64_t) * (k_ - 1));
        ckmer[k_ - 1] = c;

    }
    assert(!append || W->operator[](ind) == 0);

    // Padding after sequence to get back into default state.
    if (W->operator[](ind)==0) {
        for (size_t j = 0; j < k_; j++) {
            ind = append_pos(6, ckmer.data(), ind) * !append;
            memmove(ckmer.data(), &ckmer[1], sizeof(uint64_t) * (k_ - 1));
            ckmer[k_ - 1] = 6;
        }
        if (!append)
            ind = append_pos(0, ckmer.data(), ind);
    }

    verbose_cout("edges ", num_edges(), " / nodes ", num_nodes(), "\n");
}


bool equal_encodings(const char first, const char second) {
    return DBG_succ::encode(first) == DBG_succ::encode(second);
}

void DBG_succ::add_sequence_fast(const std::string &seq,
                                 bool add_bridge, unsigned int parallel,
                                 std::string suffix) {
    std::vector<uint32_t> label_id;

    // ther is nothing to parse
    if (!seq.size()) {
        return;
    }

	char *bridge = (char*)malloc(k_ + 2);
    memset(bridge, 'X', k_);
    bridge[k_] = seq[0];
    bridge[k_ + 1] = 0;

    size_t i = 0;
    //std::cout << "Loading next sequence with " << parallel << " threads\n";
    if (add_bridge) {
        for (i = 0; i < std::min(k_, seq.length()); ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           bridge + k_ - suffix.length(),
                           equal_encodings)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k_ + 1), DBG_succ::encode));
            }
            memmove(bridge, bridge + 1, k_);
            bridge[k_] = (i + 1 < seq.length()) ? seq[i + 1] : 'X';
        }
    }
    if (k_ < seq.length()) {
        #pragma omp parallel num_threads(parallel)
        {
            std::vector<KMer> kmer_priv;
            #pragma omp for nowait
            for (i = 0; i < seq.length() - k_; ++i) {
                if (std::equal(suffix.begin(), suffix.end(),
                               seq.c_str() + i + k_ - suffix.length(),
                               equal_encodings)) {
                    kmer_priv.push_back(KMer::from_string(
                        std::string(seq.c_str() + i, k_ + 1),
                        DBG_succ::encode
                    ));
                }
            }
            #pragma omp critical
            kmers.insert(kmers.end(),
                std::make_move_iterator(kmer_priv.begin()),
                std::make_move_iterator(kmer_priv.end())
            );
        }
        memcpy(bridge, seq.c_str() + seq.length() - k_, k_);
        bridge[k_] = 'X';
    }
    if (add_bridge) {
        for (i = 0; i < k_; ++i) {
            if (std::equal(suffix.begin(), suffix.end(),
                           bridge + k_ - suffix.length(),
                           equal_encodings)) {
                kmers.push_back(KMer::from_string(std::string(bridge, k_ + 1), DBG_succ::encode));
            }
            memmove(bridge, bridge + 1, k_);
            bridge[k_] = 'X';
        }
    }
    free(bridge);
}

void DBG_succ::construct_succ(unsigned int parallel) {

    // parallel sort of all kmers
    omp_set_num_threads(std::max(static_cast<int>(parallel), 1));
    __gnu_parallel::sort(kmers.begin(),kmers.end());

    auto last = std::unique(kmers.begin(), kmers.end());
    kmers.erase(last, kmers.end()); 

    //DEBUG: output kmers in current bin
    /*
    std::cerr << "\n";
    for (size_t i=0;i<kmers.size();++i) {
        char* curseq = kmer_to_s(kmers[i].first, alphabet, alph_size);
        std::cerr << kmers[i].first << "\t" << curseq+1 << " " << curseq[0] << " " << kmers[i].second << "\n";
        free(curseq);
    }
    */

    size_t curpos = W_stat.size();
    W_stat.resize(W_stat.size() + kmers.size());
    last_stat_safe.resize(last_stat_safe.size() + kmers.size(), true);
    coverage.resize(coverage.size() + kmers.size(),0);
    //bridge_stat.resize(bridge_stat.size() + kmers.size(), false);

    #pragma omp parallel num_threads(parallel)
    {
        #pragma omp for nowait
        for (size_t i = 0; i < kmers.size(); ++i) {
            //set last
            if (i + 1 < kmers.size()) {
                bool dup = KMer::compare_kmer_suffix(kmers[i], kmers[i + 1]);
                if (dup) {
                    last_stat_safe[curpos + i] = false;
                }
            }
            //set W
            uint8_t curW = kmers[i][0];
            if (curW == 127) {
                std::string curseq = kmers[i].to_string(alphabet);
                std::cerr << "Failure decoding kmer " << i << "\n" << kmers[i] << "\n" << curseq << "\n";
                exit(1);
            }
            if (!curW && curpos+i)
                p_ = curpos + i;
            if (i) {
                for (size_t j = i - 1; KMer::compare_kmer_suffix(kmers[j], kmers[i], 1); --j) {
                    //TODO: recalculating W is probably faster than doing a pragma for ordered
                    if (kmers[j][0] == curW) {
                        curW += alph_size;
                        break;
                    }
                    if (!j)
                        break;
                }
            }
            W_stat[curpos+i] = curW;
        }
    }
    size_t lastlet = 0;

    for (size_t i = 0; i < kmers.size(); ++i) {
        char cF = alphabet[kmers[i][k_]];
        if (cF != alphabet[lastlet]) {
            for (lastlet++; lastlet < alph_size; lastlet++) {
                F[lastlet] = curpos + i - 1;
                if (alphabet[lastlet] == cF) {
                    break;
                }
            }
        }
    }
    kmers.clear();
}

/** This function takes a character c and appends it to the end of the graph sequence
 * given that the corresponding note is not part of the graph yet.
 */
uint64_t DBG_succ::append_pos(uint64_t c, uint64_t *ckmer, uint64_t i) {

    // check that the last position of the graph is indeed a terminal
    assert((*W)[p_] == 0);
    uint64_t *p;
    bool append = false;
    uint64_t c_p = get_node_end_value(p_);
    // get range of identical nodes (without W) pos current end position
    std::pair<uint64_t, uint64_t> R;
    //std::pair<uint64_t, uint64_t> R = get_equal_node_range(p_);
    //fprintf(stdout, "range [%i %i]\n", (int) R.first, (int) R.second);

    if (!i) {
        p = &(p_);
        append = true;
        R.second = succ_last(*p);
    } else {
        append = false;
        R = get_equal_node_range(p_);
        p = &(R.first);
    }

    // get position of first occurence of c in W after p
    uint64_t next_c = succ_W(*p, c);
    // check if c is part of range
    bool exist_c = (next_c <= R.second);
    if (!exist_c) {
        // get position of first occurence of c- in W after p
        next_c = succ_W(*p, c + alph_size);
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
            for (i=0; (*W)[*p] == 0; ++i, ++(*p)) {
                if (*p == p_)
                    continue;

                assert(i < 2);
                assert(*p < R.second);

                W->remove(*p);
                last->deleteBit(*p);
                update_F(c_p, false);
                if (*p <= p_) {
                    p_--;
                    assert((*W)[*p] == 0);
                }
                if (*p <= next_c)
                    next_c--;
            }
            return fwd(next_c);
        }
        uint64_t p_new = fwd(next_c);
        // remove old terminal symbol
        last->deleteBit(*p);
        W->remove(*p);
        // adapt position if altered by previous deletion
        p_new -= (*p < p_new);
        // insert new terminal symbol
        // we have to insert 0 into last as the node already existed in the range
        // and the terminal symbol is always first
        last->insertBit(p_new, false);
        W->insert(p_new, 0);
        // update new terminal position
        *p = p_new;
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
        if (!append) {
            //we need to insert a new pointer
            if (*p == p_)
                (*p)++;
            if (W->operator[](*p)) {
                //if no placeholder exists
                W->insert(*p, 0);
                update_F(c_p, true);
                last->insertBit(*p, false);
                R.second++;
                if (*p <= p_) {
                    (p_)++;
                    assert(W->operator[](p_) == 0);
                }
            }
        }
        // get position of last occurence of c before p (including p - 1)
        uint64_t last_c = pred_W(*p - 1, c);
        // if this position exists
        if (last_c > 0) {
            uint64_t x = fwd(last_c);
            assert((*(last))[x]); // this should always be true - unless x is 0 - I do not get the logic in the reference implementation

            // check, if there are any c or c- symbols following after position p
            uint64_t next_c = succ_W(*p + 1, c);
            uint64_t next_cm = succ_W(*p + 1, c + alph_size);
            // there is no c between p and next_cm and next_cm is a c- ==> we should add a c-
            // all nodes between W[i] = c and W[j] = c- share a common suffix of length k-1
            bool minus1 = (next_cm < next_c);
            // check, if we share a k-1 suffix with last_c
            if (!minus1) {
                if (ckmer) {
                    minus1 = compare_node_suffix(ckmer, last_c);
                } else {
                    minus1 = compare_node_suffix(*p, last_c);
                }
            }

            // adding a new node can influence following nodes that share a k-1 suffix with the
            // new node -> need to adapt the respektive cc to a cc-
            bool minus2 = false;
            if (next_c < W->size()) {
                if (ckmer) {
                    minus2 = compare_node_suffix(ckmer, next_c);
                } else {
                    minus2 = compare_node_suffix(*p, next_c);
                }
                if (minus2) {
                    W_set_value(next_c, (*(W))[next_c] + alph_size);
                }
            }

            W_set_value(*p, minus1 ? c + alph_size : c);
            // after we are done, assert that the order within the range we created
            // is still valid within W
            if (*p - R.second > 0) {
                sort_W_locally(*p, R.second);
            }

            // if minus1 is true, we share a k-1 suffix with the node at
            // last_c and thus need to adapt our insertion position by -1,
            // as we would like to insert before it. Otherwise we insert directly after
            // it as we are now sorted after it.
            if (minus1) {
                *p = x;
                last->insertBit(x, false);
                W->insert(x, 0);
            } else if (minus2) {
                *p = x + 1;
                last->insertBit(x + 1, false);
                W->insert(x + 1, 0);
            // no node shares a k-1 suffix with last_c and thus the new node comes after
            // the forward of last_c (as the current node came after last_c as well)
            } else {
                *p = x + 1;
                last->insertBit(x + 1, true);
                W->insert(x + 1, 0);
            }
        } else {
            uint64_t x = F[c] + 1;
            uint64_t next_c = succ_W(*p + 1, c);
            bool minus = false;
            if (next_c < W->size()) {
                if (ckmer) {
                    minus = compare_node_suffix(ckmer, next_c);
                } else {
                    minus = compare_node_suffix(*p, next_c);
                }
            }
            W_set_value(*p, c);
            if (*p - R.second > 0) {
                sort_W_locally(*p, R.second);
            }
            *p = x;
            if (minus) {
                W_set_value(next_c, (*(W))[next_c] + alph_size);
                last->insertBit(x, false);
            } else {
                last->insertBit(x, true);
            }
            W->insert(x, 0);
        }
        if (*p < p_ || (!append && *p == p_)) {
            (p_)++;
            assert(W->operator[](p_) == 0);
        }
        update_F(c, true);
    }
    // update sorting at new location of p
    // with this we assert that $ is always inserted at the first position
    // of a range of equal nodes --> this will help us to prevent multiple insertions
    // of already existing nodes
    R = get_equal_node_range(*p);
    if (R.second - R.first > 0) {
        sort_W_locally(R.first, R.second);
        *p = R.first;
        if (!append && *p == p_)
            (*p)++;
        /*
        while ((*(W))[*p] != 0)
            (*p)--;
        */
        assert(W->operator[](*p)==0);
    }
    return *p;
}

/** This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's arrays. In almost all cases this will not produce a valid graph and
 * should only be used as a helper in the parallel merge procedure.
 */
void DBG_succ::append_graph(const DBG_succ &G) {
    verbose_cout("    adding ", G.W->size(), " edges\n");

    // handle last and W
    for (size_t j = 1, curr_pos = W->size(); j < G.W->size(); ++j, ++curr_pos) {
        W->insert(curr_pos, G.get_W(j));
        last->insertBit(curr_pos, G.get_last(j));
    }

    verbose_cout("new total edges: ", W->size(), "\n");

    // handle F
    assert(F.size() == G.F.size());
    for (size_t j = 0; j < F.size(); ++j) {
        F.at(j) += G.F.at(j);
    }
}

/**
 * This function takes a pointer to a graph structure and concatenates the arrays W, last
 * and F to this graph's static containers last_stat and W_stat. In almost all cases
 * this will not produce a valid graph and should only be used as a helper in the
 * parallel merge procedure.
 */
void DBG_succ::append_graph_static(const DBG_succ &G) {
    verbose_cout("    adding ", G.W->size(), " edges\n");

    assert(dynamic_cast<wavelet_tree_dyn*>(G.W));
    auto G_W_stat = dynamic_cast<wavelet_tree_dyn*>(G.W)->to_vector();

    W_stat.insert(W_stat.end(), G_W_stat.begin() + 1, G_W_stat.end());

    for (size_t i = 1; i < G.W->size(); ++i) {
        last_stat.push_back(G.get_last(i));
    }

    verbose_cout("new total edges: ", W->size(), "\n");

    // handle F
    assert(F.size() == G.F.size());
    for (size_t j = 0; j < F.size(); ++j) {
        F.at(j) += G.F.at(j);
    }
}

// Given an edge list, remove them from the graph.
void DBG_succ::remove_edges(const std::set<uint64_t> &edges) {
    uint64_t shift = 0;

    for (const auto &edge : edges) {
        assert(edge >= shift);
        uint64_t edge_id = edge - shift;

        uint64_t d = W->operator[](edge_id);
        if (d < alph_size) {
            //fix W array
            uint64_t next = edge_id + 1;
            uint64_t j = succ_W(next, d);
            for (uint64_t i = next; i < j; ++i) {
                if (W->operator[](i) == d + alph_size) {
                    W_set_value(i, d);
                    break;
                }
            }
        }
        W->remove(edge_id);
        update_F(get_node_end_value(edge_id), false);
        // If the current node has multiple outgoing edges,
        // remove one of the 0s from last instead of 1.
        if (get_last(edge_id) && (edge >= shift + 1)
                              && !get_last(edge_id - 1)) {
            last->deleteBit(edge_id - 1);
        } else {
            last->deleteBit(edge_id);
        }
        if (edge_id <= p_) {
            p_--;
            assert(W->operator[](p_) == 0);
        }
        shift++;
    }
}

void DBG_succ::add_sink(unsigned int parallel, std::string suffix) {
    add_sequence_fast(start, false, parallel, suffix);
    add_sequence_fast(sink, true, parallel, suffix);
}
