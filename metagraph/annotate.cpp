#include "annotate.hpp"

#include <string>
#include <algorithm>
#include <pthread.h>

#include "datatypes.hpp"


namespace annotate {

template<typename LabelType>
ColorCompressed<LabelType>::ColorCompressed(const std::vector<ColorCompressed<LabelType>> &categories,
                                            const std::vector<std::set<size_t>> &merge_plan) :
                                             annotation_curr_(NULL) {
                                            }

template<typename LabelType>
void ColorCompressed<LabelType>::set(Index i, const LabelType &label) {
    
    auto id_it = label_to_id_.find(label);
    uint64_t id = 0;

    // current label does not exist yet -> assign new column
    if (id_it == label_to_id_.end()) {
        id = id_to_label_.size();
        id_to_label_.emplace_back(label);
        label_to_id_[label] = id;
        bitmatrix_.push_back(NULL);

        if (annotation_curr_) {
            flush();
            delete annotation_curr_;
        }
        annotation_curr_ = new sdsl::bit_vector(graph_size_, 0);
        label_curr_ = label;

    // decompress existing column
    } else if (!annotation_curr_) {
        id = id_it->second;
        if (bitmatrix_.at(id - 1) != NULL) {
            annotation_curr_ = inflate_column(id - 1);
        } else {
            annotation_curr_ = new sdsl::bit_vector(graph_size_, 0);
        }
    }

    annotation_curr_->operator[](i) = 1;
}


template<typename LabelType>
void ColorCompressed<LabelType>::flush() {
    
    if (annotation_curr_) {
        auto id_it = label_to_id_.find(label_curr_);
        if (id_it != label_to_id_.end()) {
            uint64_t id = id_it->second;
            while(bitmatrix_.size() < id)
                bitmatrix_.push_back(NULL);

            if (bitmatrix_.at(id - 1) != NULL) {
                sdsl::bit_vector* tmp = inflate_column(id - 1);
                *annotation_curr_ |= *tmp;
                delete tmp;
                delete bitmatrix_.at(id - 1);
            }

            bitmatrix_.at(id - 1) = new sdsl::sd_vector<>(*annotation_curr_);
            delete annotation_curr_;
            annotation_curr_ = NULL;
        }
    }
}


template<typename LabelType>
sdsl::bit_vector* ColorCompressed<LabelType>::inflate_column(const uint64_t id) const { 

   sdsl::sd_vector<>* col = bitmatrix_.at(id);
   sdsl::select_support_sd<> slct = sdsl::select_support_sd<>(col);
   sdsl::rank_support_sd<> rank = sdsl::rank_support_sd<>(col);

   size_t maxrank = rank(col->size());
   sdsl::bit_vector* result = new sdsl::bit_vector(col->size(), 0);
   size_t idx = 0;
   for (size_t i = 1; i <= maxrank; ++i) {
       idx = slct(i);
       if (idx < col->size()) {
           result->operator[](idx) = col->operator[](idx);
           continue;
       }
       break;
   }
        
   return result;
}

// void annotate_kmer(DBG_succ *G, sdsl::bit_vector *annotation_curr,
//                    std::string &kmer, uint64_t &idx, bool ignore) {

//     // we just need to walk one step in the path
//     if (idx > 0) {
//         uint64_t s = G->get_alphabet_number(kmer[kmer.length() - 2]);
//         uint64_t pl = G->pred_last(idx - 1);
//         uint64_t rl = std::min(G->succ_W(pl + 1, s), G->succ_W(pl + 1, s + G->alph_size));
//         uint64_t ru = std::max(G->pred_W(idx, s), G->pred_W(idx, s + G->alph_size));
//         rl = G->outgoing(rl, s);
//         ru = G->outgoing(ru, s);
//         idx = (ru > rl) ? ru : rl;
//     // we need to look up the full kmer
//     } else {
//         idx = G->index(kmer, G->k);
//     }
//     //std::cerr << "kmer: " << kmer << " idx: " << idx << " ignore: " << (ignore?"yes":"no") << " label_id: " << label_id << std::endl;
//     if ((idx == 0) || ignore)
//         return;
//     idx = G->outgoing_edge_idx(idx, G->get_alphabet_number(kmer[kmer.length() - 1]));
//     if (idx == 0)
//         return;

//     annotation_curr->operator[](idx) = 1;
// }

// void annotate_seq(DBG_succ *G, Config *config, kstring_t &seq, kstring_t &label,
//                   uint64_t start, uint64_t end,
//                   pthread_mutex_t *anno_mutex) {

//     std::string curr_kmer;
//     std::string label_str = std::string(label.s);
//     uint32_t label_id;

//     if (anno_mutex)
//         pthread_mutex_lock(anno_mutex);

//      // does the current label already have an ID?
//     std::unordered_map<std::string, uint32_t>::iterator id_it = G->label_to_id_map.find(label_str);
//     sdsl::bit_vector* annotation_curr = NULL;
//     if (id_it == G->label_to_id_map.end()) {
//         label_id = (uint32_t) G->id_to_label.size();
//         G->id_to_label.push_back(label_str);
//         G->label_to_id_map[label_str] = label_id;
//         G->annotation_full.push_back(NULL);
//         annotation_curr = new sdsl::bit_vector(G->get_size(), 0);
//         if (config->verbose)
//             std::cout << "added label ID " << label_id
//                       << " for label string " << label_str << std::endl;
//     } else {
//         label_id = id_it->second;
//         if (G->annotation_full.at(label_id - 1) != NULL) {
//             annotation_curr = inflate_annotation(G, label_id - 1);
//         } else {
//             annotation_curr = new sdsl::bit_vector(G->get_size(), 0);
//         }
//     }

//     if (anno_mutex)
//         pthread_mutex_unlock(anno_mutex);

//     end = (end == 0) ? seq.l : end;

//     uint64_t previous_idx = 0;
//     size_t i;
//     for (i = start; i < end; ++i) {

//         if (config->verbose && i > 0 && i % 1'000 == 0) {
//             std::cout << "." << std::flush;
//             if (!anno_mutex && (i % 10'000 == 0))
//                 std::cout << i << " kmers added" << std::endl;
//         }

//         if (curr_kmer.size() <= G->k) {
//             curr_kmer.push_back(seq.s[i]);
//             continue;
//         }
//         assert(curr_kmer.size() == G->k + 1);
//         annotate_kmer(G, annotation_curr, curr_kmer, previous_idx, (i % config->frequency) > 0);

//         //std::cerr << curr_kmer << ":" << std::string(label.s) << std::endl;
//         curr_kmer.push_back(seq.s[i]);
//         curr_kmer = curr_kmer.substr(1, G->k + 1);
//     }
//     // add last kmer and label to database
//     if (curr_kmer.size() == G->k + 1)
//         annotate_kmer(G, annotation_curr, curr_kmer, previous_idx, (i % config->frequency) > 0);

//     if (anno_mutex)
//         pthread_mutex_lock(anno_mutex);

//     while (G->annotation_full.size() < label_id)
//         G->annotation_full.push_back(NULL);

//     if (G->annotation_full.at(label_id - 1) != NULL) {
//         sdsl::bit_vector* tmp = inflate_annotation(G, label_id - 1);
//         *annotation_curr |= *tmp;
//         delete tmp;
//         delete G->annotation_full.at(label_id - 1);
//     }
//     //(G->annotation_full.at(label_id - 1)) = new sdsl::rrr_vector<63>(*annotation_curr);
//     G->annotation_full.at(label_id - 1) = new sdsl::sd_vector<>(*annotation_curr);

//     if (anno_mutex)
//         pthread_mutex_unlock(anno_mutex);

//     delete annotation_curr;
//     annotation_curr = NULL;
// }

// std::vector<uint32_t> classify_path(DBG_succ* G, std::vector<uint64_t> path) {

//     //uint32_t curr_anno;
//     std::vector<uint32_t> labels(G->annotation_full.size(), 0);
//     //std::vector<uint32_t> current_combination;
//     //std::map<uint32_t, uint64_t> label_counter;

//     // collect all annotated combinations for the path
//     // take majority vote as consensus for now
//     for (size_t i = 0; i < path.size(); i += 1) {
//         for (size_t j = 0; j < G->annotation_full.size(); ++j) {
//             labels.at(j) = (*(G->annotation_full.at(j)))[path.at(i)] ? 1 : std::max(labels.at(j), 0u);
//         }
//     }
//     return labels;
// }


// std::set<uint32_t> classify_read(DBG_succ *G, kstring_t &read, uint64_t max_distance) {

//     // containers for label information
//     std::vector<uint32_t> path_labels;
//     // TODO: unordered_set
//     std::set<uint32_t> all_labels;

//     // get alignment of the read to the graph
//     std::string read_str = std::string(read.s);
//     std::vector<HitInfo> alignment = G->index_fuzzy(read_str, max_distance);

//     // classify hits
//     for (size_t i = 0; i < alignment.size(); i++) {
//         path_labels = classify_path(G, alignment.at(i).path);
//         for (size_t j = 0; j < path_labels.size(); j++) {
//             if (path_labels.at(j) > 0)
//                 all_labels.insert(j);
//         }
//     }

//     return all_labels;
// }

// // write annotation to screen
// void annotationToScreen(DBG_succ *G) {
//     for (size_t i = 1; i < G->get_size(); ++i) {
//         std::cout << i;
//         for (size_t j = 0; j < G->annotation_full.size(); ++j) {
//             if ((*(G->annotation_full.at(j)))[i])
//                 std::cout << " : " << G->id_to_label.at(j + 1);
//         }
//         std::cout << std::endl;
//     }
// }

} // namespace annotate
