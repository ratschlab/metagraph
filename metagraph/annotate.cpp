#include "annotate.hpp"

#include <string>
#include <algorithm>
#include <pthread.h>

#include "datatypes.hpp"


namespace annotate {

void annotate_kmer(DBG_succ *G, sdsl::bit_vector *annotation_curr,
                   std::string &kmer, uint64_t &idx, bool ignore) {

    // we just need to walk one step in the path
    if (idx > 0) {
        uint64_t s = G->get_alphabet_number(kmer[kmer.length() - 2]);
        uint64_t pl = G->pred_last(idx - 1);
        uint64_t rl = std::min(G->succ_W(pl + 1, s), G->succ_W(pl + 1, s + G->alph_size));
        uint64_t ru = std::max(G->pred_W(idx, s), G->pred_W(idx, s + G->alph_size));
        rl = G->outgoing(rl, s);
        ru = G->outgoing(ru, s);
        idx = (ru > rl) ? ru : rl;
    // we need to look up the full kmer
    } else {
        idx = G->index(kmer, G->k);
    }
    //std::cerr << "kmer: " << kmer << " idx: " << idx << " ignore: " << (ignore?"yes":"no") << " label_id: " << label_id << std::endl;
    if ((idx == 0) || ignore)
        return;
    idx = G->outgoing_edge_idx(idx, G->get_alphabet_number(kmer[kmer.length() - 1]));
    if (idx == 0)
        return;

    annotation_curr->operator[](idx) = 1;
}

sdsl::bit_vector* inflate_annotation(DBG_succ* G, uint64_t id) {
    sdsl::select_support_sd<> slct = sdsl::select_support_sd<>(G->annotation_full.at(id));
    sdsl::rank_support_sd<> rank = sdsl::rank_support_sd<>(G->annotation_full.at(id));
    size_t maxrank = rank(G->annotation_full.at(id)->size());
    sdsl::bit_vector* result = new sdsl::bit_vector(G->annotation_full.at(id)->size(), 0);
    size_t idx;
    for (size_t i = 1; i <= maxrank; ++i) {
        idx = slct(i);
        if (idx < G->annotation_full.at(id)->size()) {
            result->operator[](idx) = G->annotation_full.at(id)->operator[](idx);
        }
    }
    return result;
}

void annotate_seq(DBG_succ *G, Config *config, kstring_t &seq, kstring_t &label,
                  uint64_t start, uint64_t end,
                  pthread_mutex_t *anno_mutex) {

    std::string curr_kmer;
    std::string label_str = std::string(label.s);
    uint32_t label_id;

    if (anno_mutex)
        pthread_mutex_lock(anno_mutex);

     // does the current label already have an ID?
    std::unordered_map<std::string, uint32_t>::iterator id_it = G->label_to_id_map.find(label_str);
    sdsl::bit_vector* annotation_curr = NULL;
    if (id_it == G->label_to_id_map.end()) {
        label_id = (uint32_t) G->id_to_label.size();
        G->id_to_label.push_back(label_str);
        G->label_to_id_map[label_str] = label_id;
        G->annotation_full.push_back(NULL);
        annotation_curr = new sdsl::bit_vector(G->get_size(), 0);
        if (config->verbose)
            std::cout << "added label ID " << label_id
                      << " for label string " << label_str << std::endl;
    } else {
        label_id = id_it->second;
        if (G->annotation_full.at(label_id - 1) != NULL) {
            annotation_curr = inflate_annotation(G, label_id - 1);
        } else {
            annotation_curr = new sdsl::bit_vector(G->get_size(), 0);
        }
    }

    if (anno_mutex)
        pthread_mutex_unlock(anno_mutex);

    end = (end == 0) ? seq.l : end;

    uint64_t previous_idx = 0;
    size_t i;
    for (i = start; i < end; ++i) {

        if (config->verbose && i > 0 && i % 1'000 == 0) {
            std::cout << "." << std::flush;
            if (!anno_mutex && (i % 10'000 == 0))
                std::cout << i << " kmers added" << std::endl;
        }

        if (curr_kmer.size() <= G->k) {
            curr_kmer.push_back(seq.s[i]);
            continue;
        }
        assert(curr_kmer.size() == G->k + 1);
        annotate_kmer(G, annotation_curr, curr_kmer, previous_idx, (i % config->frequency) > 0);

        //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
        curr_kmer.push_back(seq.s[i]);
        curr_kmer = curr_kmer.substr(1, G->k + 1);
    }
    // add last kmer and label to database
    if (curr_kmer.size() == G->k + 1)
        annotate_kmer(G, annotation_curr, curr_kmer, previous_idx, (i % config->frequency) > 0);

    if (anno_mutex)
        pthread_mutex_lock(anno_mutex);

    while (G->annotation_full.size() < label_id)
        G->annotation_full.push_back(NULL);

    if (G->annotation_full.at(label_id - 1) != NULL) {
        sdsl::bit_vector* tmp = inflate_annotation(G, label_id - 1);
        *annotation_curr |= *tmp;
        delete tmp;
        delete G->annotation_full.at(label_id - 1);
    }
    //(G->annotation_full.at(label_id - 1)) = new sdsl::rrr_vector<63>(*annotation_curr);
    G->annotation_full.at(label_id - 1) = new sdsl::sd_vector<>(*annotation_curr);

    if (anno_mutex)
        pthread_mutex_unlock(anno_mutex);

    delete annotation_curr;
    annotation_curr = NULL;
}

std::vector<uint32_t> classify_path(DBG_succ* G, std::vector<uint64_t> path) {

    //uint32_t curr_anno;
    std::vector<uint32_t> labels(G->annotation_full.size(), 0);
    //std::vector<uint32_t> current_combination;
    //std::map<uint32_t, uint64_t> label_counter;

    // collect all annotated combinations for the path
    // take majority vote as consensus for now
    for (size_t i = 0; i < path.size(); i += 1) {
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


std::set<uint32_t> classify_read(DBG_succ *G, kstring_t &read, uint64_t max_distance) {

    // containers for label information
    std::vector<uint32_t> path_labels;
    // TODO: unordered_set
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

// write annotation to screen
void annotationToScreen(DBG_succ *G) {
    for (size_t i = 1; i < G->get_size(); ++i) {
        std::cout << i;
        for (size_t j = 0; j < G->annotation_full.size(); ++j) {
            if ((*(G->annotation_full.at(j)))[i])
                std::cout << " : " << G->id_to_label.at(j + 1);
        }
        std::cout << std::endl;
    }
}

} // namespace annotate
