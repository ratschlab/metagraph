#include "annotate.hpp"

#include <string>
#include <algorithm>
#include <pthread.h>

#include "dbg_succinct_libmaus.hpp"
#include "datatypes.hpp"

namespace annotate {

    typedef uint64_t TAlphabet;

    void annotate_kmer(DBG_succ* G, std::string &kmer, uint32_t &label_id, uint64_t &idx, pthread_mutex_t* anno_mutex, bool ignore) {

        // we just need to walk one step in the path
        if (idx > 0) {
            TAlphabet s = G->get_alphabet_number(kmer[kmer.length() - 2]);
            uint64_t rl = std::min(G->succ_W(G->pred_last(idx - 1) + 1, s), G->succ_W(G->pred_last(idx - 1) + 1, s + G->alph_size));
            uint64_t ru = std::max(G->pred_W(idx, s), G->pred_W(idx, s + G->alph_size));
            rl = G->outgoing(rl, s);
            ru = G->outgoing(ru, s);
            idx = (ru > rl) ? ru : rl;
        // we need to look up the full kmer
        } else {
            idx = G->index(kmer, G->k);
        }
        //std::cerr << "kmer: " << kmer << " idx: " << idx << " ignore: " << (ignore?"yes":"no") << " label_id: " << label_id << std::endl;
        //assert(idx > 0);
        if ((idx == 0) || ignore)
            return;
        idx = G->outgoing_edge_idx(idx, G->get_alphabet_number(kmer[kmer.length() - 1]));
        if (idx == 0)
            return;

        //std::vector<uint32_t> curr_combo;
        //sdsl::rrr_vector<64> curr_color;

        if (anno_mutex)
            pthread_mutex_lock(anno_mutex);

        (*(G->annotation_curr))[idx] = 1;

/*
        // get annotation of current kmer
       // size_t curr_hash = G->annotation[idx];
        if (curr_hash > 0)
            curr_color = get_curr_color(G, idx);
            //curr_combo = get_curr_combination(G->combination_vector, G->annotation_map[curr_hash]);
       
        //if (anno_mutex)
        //    pthread_mutex_unlock(anno_mutex);

        // check whether current label is already part of current annotation
        // and add it if not
        if (!std::binary_search(curr_combo.begin(), curr_combo.end(), label_id)) {

            curr_combo = add_to_combination(curr_combo, label_id); 
            curr_hash = G->hasher.bithash(curr_combo);

           // if (anno_mutex)
           //     pthread_mutex_lock(anno_mutex);

            G->annotation[idx] = curr_hash;
            //annotation_map[curr_hash] = curr_anno;
            if (G->annotation_map.find(curr_hash) == G->annotation_map.end()) {
                G->annotation_map[curr_hash] = insert_new_combination(G->combination_vector, curr_combo);
                G->combination_count++;
            }
        }
*/
        if (anno_mutex)
            pthread_mutex_unlock(anno_mutex);
        //}
    }

    void annotate_seq(DBG_succ* G, kstring_t &seq, kstring_t &label, uint64_t start, uint64_t end, pthread_mutex_t* anno_mutex) {

        std::string curr_kmer;
        std::string label_str = std::string(label.s);
        uint32_t label_id;

        if (anno_mutex)
            pthread_mutex_lock(anno_mutex);

        if (G->annotation_curr == NULL) {

             // does the current label already have an ID?
            std::unordered_map<std::string, uint32_t>::iterator id_it = G->label_to_id_map.find(label_str);
            if (id_it == G->label_to_id_map.end()) {
                label_id = (uint32_t) G->id_to_label.size();
                G->id_to_label.push_back(label_str);
                G->label_to_id_map[label_str] = label_id;
                G->annotation_curr = new sdsl::bit_vector(G->get_size(), 0);
                if (G->config->verbose)
                    std::cout << "added label ID " << label_id << " for label string " << label_str << std::endl;
            } else { 
                label_id = id_it->second;
                sdsl::select_support_sd<> slct = sdsl::select_support_sd<>(G->annotation_full.at(label_id - 1));
                sdsl::rank_support_sd<> rank = sdsl::rank_support_sd<>(G->annotation_full.at(label_id - 1));
                size_t maxrank = rank(G->annotation_full.at(label_id - 1)->size());
                G->annotation_curr = new sdsl::bit_vector(G->annotation_full.at(label_id - 1)->size(), 0);
                size_t idx;
                for (size_t i = 1; i <= maxrank; ++i) {
                    idx = slct(i);
                    if (i < G->annotation_full.at(label_id - 1)->size())
                        (*(G->annotation_curr))[i] = (*(G->annotation_full.at(label_id - 1)))[i];
                }
            }
        }

        if (anno_mutex)
            pthread_mutex_unlock(anno_mutex);

        end = (end == 0) ? seq.l : end;

        uint64_t previous_idx = 0;
        size_t i;
        for (i = start; i < end; ++i) {

            if (G->config->verbose && i > 0 && i % 1000 == 0) {
                std::cout << "." << std::flush;
                if (!anno_mutex && (i % 10000 == 0))
                    std::cout << i << " kmers added" << std::endl;
            }

            if (curr_kmer.size() <= G->k) {
                curr_kmer.push_back(seq.s[i]);
                continue;
            }
            assert(curr_kmer.size() == G->k + 1);
            annotate_kmer(G, curr_kmer, label_id, previous_idx, anno_mutex, (i % G->config->frequency) > 0);
            
            //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
            curr_kmer.push_back(seq.s[i]);
            curr_kmer = curr_kmer.substr(1, G->k + 1);
        }
        // add last kmer and label to database
        if (curr_kmer.size() == G->k + 1)
            annotate_kmer(G, curr_kmer, label_id, previous_idx, anno_mutex, (i % G->config->frequency) > 0);
        //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
        
        if (label_id - 1 == G->annotation_full.size()) {
            //G->annotation_full.push_back(new sdsl::rrr_vector<63>(*(G->annotation_curr)));
            G->annotation_full.push_back(new sdsl::sd_vector<>(*(G->annotation_curr)));
        } else {
            delete G->annotation_full.at(label_id - 1);
            //(G->annotation_full.at(label_id - 1)) = new sdsl::rrr_vector<63>(*(G->annotation_curr));
            (G->annotation_full.at(label_id - 1)) = new sdsl::sd_vector<>(*(G->annotation_curr));
        }
        delete G->annotation_curr;
        G->annotation_curr = NULL;
    }   

    std::vector<uint32_t> classify_path(DBG_succ* G, std::vector<uint64_t> path) {
        
        //uint32_t curr_anno;
        std::vector<uint32_t> labels(G->annotation_full.size(), 0);
        //std::vector<uint32_t> current_combination;
        //std::map<uint32_t, uint64_t> label_counter;

        // collect all annotated combinations for the path
        // take majority vote as consensus for now
        for (size_t i = 0; i < path.size(); i++) {
            for (size_t j = 0; j < G->annotation_full.size(); ++j) {
                labels.at(j) = (*(G->annotation_full.at(j)))[path.at(i)] ? 1 : std::max(labels.at(j), 0u);
            }
            /*curr_anno = G->annotation.at(path.at(i));
            if (curr_anno > 0) {
                current_combination = get_curr_combination(G->combination_vector, G->annotation_map[curr_anno]);
                for (std::vector<uint32_t>::iterator c = current_combination.begin(); c != current_combination.end(); c++) {
                    if (label_counter.find(*c) != label_counter.end()) {
                        label_counter[*c] += 1;
                    } else {
                        label_counter[*c] = 1;
                    }
                }
            }*/
        }

        /*
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
        }*/

        return labels;
    }


    std::set<uint32_t> classify_read(DBG_succ* G, kstring_t &read, uint64_t max_distance) {

        // containers for label information
        std::vector<uint32_t> path_labels;
        std::set<uint32_t> all_labels;

        // get alignment of the read to the graph
        std::string read_str = std::string(read.s);
        std::vector<HitInfo> alignment = G->index_fuzzy(read_str, max_distance);

        // classify hits
        for (size_t i = 0; i < alignment.size(); i++) {
            path_labels = classify_path(G, alignment.at(i).path);
            for (size_t j = 0; j < path_labels.size(); j++) {
                if (path_labels.at(j) > 0)
                    all_labels.insert(j);
            }
        }

        return all_labels;
    }

    std::vector<uint32_t> get_curr_combination(std::vector<uint32_t> &combinations, uint32_t idx) {
        std::vector<uint32_t>::iterator start = combinations.begin() + idx;
        uint32_t length = *start;
        std::vector<uint32_t>::iterator end = combinations.begin() + idx + length + 1;
        start++;

        std::vector<uint32_t> curr_combo(length);
        std::copy(start, end, curr_combo.begin());

        return curr_combo;
    }

/*    sdsl::rrr_vector<63> get_curr_color(DBG_succ* G, uint64_t idx) {
        uint64_t color_begin = G->annotation_new->select1(idx);
        uint64_t color_end = G->annotation_new->select1(idx + 1);
        uint64_t color_idx = 0;
        uint64_t mask = 1;
        for (size_t i = color_begin; i < color_end; ++i) {
            if (G->annotation_colors[i])
                color_idx |= mask;
            mask <<= 1;
        }
        return G->colors_to_bits.at(color_idx);
    }
*/
    /*bool in_combo(std::vector<uint32_t> &combo, uint32_t value) {
        
        std::cerr << "search for " << value << " at idx " << idx << std::endl;
        std::vector<uint32_t>::iterator start = anno.begin() + idx;
        std::cerr << "start val " << *start << std::endl;
        std::vector<uint32_t>::iterator end = anno.begin() + idx + *start;
        start++;

        return std::binary_search(start, end, value);
    };*/

    std::vector<uint32_t> add_to_combination(std::vector<uint32_t> &old_combination, uint32_t value) {

        std::vector<uint32_t> new_combination;

        if (old_combination.size() > 0) {
            std::vector<uint32_t>::iterator insert_idx = std::upper_bound(old_combination.begin(), old_combination.end(), value);
            std::vector<uint32_t>::iterator it;
            for (it = old_combination.begin(); it != insert_idx; ++it)
                new_combination.push_back(*it);
            new_combination.push_back(value);
            for (it = insert_idx; it != old_combination.end(); ++it)
                new_combination.push_back(*it);
        } else {
            new_combination.push_back(value);
        }

        return new_combination;
    }

    uint32_t insert_new_combination(std::vector<uint32_t> &combination_vector, std::vector<uint32_t> &new_combination) {

        uint32_t insert_pos = combination_vector.size();
        //std::cerr << "ins pos " << insert_pos << std::endl;

        // are we close to the next billion? -> pre-allocate another billion entries
        // this prevents the vector size from doubling with each re-allocation
        // we give up on the constant amortized cost, but this is ok for here

        if ((insert_pos + 1000) % 1000000000 == 0)
            combination_vector.reserve(insert_pos + 1000 + 1000000000);

        combination_vector.push_back(new_combination.size());
        for (std::vector<uint32_t>::iterator it = new_combination.begin(); it != new_combination.end(); ++it)
            combination_vector.push_back(*it);
        return insert_pos;
    }


    // write annotation to screen
    void annotationToScreen(DBG_succ* G) {
        for (size_t i = 1; i < G->get_size(); ++i) {
            std::cout << i;
            for (size_t j = 0; j < G->annotation_full.size(); ++j) {
                if ((*(G->annotation_full.at(j)))[i])
                    std::cout << " : " << G->id_to_label.at(j + 1);
            }
            std::cout << std::endl;
        }
        /*std::deque<uint32_t>::iterator ait = G->annotation.begin();
        std::vector<uint32_t>::iterator vit;
        std::vector<uint32_t> curr_comb;
        for (; ait != G->annotation.end(); ++ait) {
            std::cout << *ait << " : " << G->annotation_map[*ait] << " : ";
            if (*ait > 0) {
                curr_comb = get_curr_combination(G->combination_vector, G->annotation_map[*ait]);
                //sit = annotation_map[*ait].begin();
                for (vit = curr_comb.begin(); vit != curr_comb.end(); ++vit) {
                    std::cout << G->id_to_label.at(*vit) << " ";
                }
            }
            std::cout << std::endl;
        }*/
        //std::cerr << "ANNO MAP " << std::endl;
        //for (std::unordered_map<uint32_t, uint32_t>::iterator itt = annotation_map.begin(); itt != annotation_map.end(); itt++)
        //    std::cerr << itt->first << ":" << itt->second << "; ";
        //std::cerr << std::endl;
    }

                
}
