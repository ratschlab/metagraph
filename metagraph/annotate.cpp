#include "annotate.hpp"

#include <pthread.h>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "datatypes.hpp"
#include "serialization.hpp"


namespace annotate {

// // Merge constructor
// ColorCompressed::ColorCompressed(const std::vector<ColorCompressed> &categories,
//                                  const std::vector<std::vector<size_t>> &merge_plan)
//       : annotation_curr_(NULL) {
//     std::cout << categories.size();
//     std::cout << merge_plan.size();
// }

ColorCompressed::SetStr ColorCompressed::get(Index i) const {
    assert(i < graph_size_);

    const_cast<ColorCompressed*>(this)->flush();

    SetStr label;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j]);
        if (bitmatrix_[j]->operator[](i)) {
            label.insert(id_to_label_[j]);
        }
    }
    return label;
}

std::vector<uint64_t> ColorCompressed::get_row(Index i) const {
    assert(i < graph_size_);

    std::vector<uint64_t> row((bitmatrix_.size() + 63) >> 6);

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j]);
        if (bitmatrix_[j]->operator[](i)) {
            row[j >> 6] |= (1llu << (j % 64));
        }
    }
    return row;
}

void ColorCompressed::serialize_uncompressed_rows(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();
    std::ofstream outstream(filename);
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, graph_size_);

    std::cout << graph_size_ << " " << bitmatrix_.size() << "\n";

    for (size_t i = 0; i < graph_size_; ++i) {
        auto row = get_row(i);
        libmaus2::util::NumberSerialisation::serialiseNumberVector(outstream, row);
    }
    outstream.close();
}

bool ColorCompressed::has_label(Index i, const SetStr &label) const {
    assert(i < graph_size_);

    for (auto &column_name : label) {
        if (!has_label(i, column_name))
            return false;
    }
    return true;
}

bool ColorCompressed::has_label(Index i, const std::string &label) const {
    assert(i < graph_size_);

    const_cast<ColorCompressed*>(this)->flush();

    auto it = label_to_id_.find(label);
    if (it == label_to_id_.end())
        return false;

    return bitmatrix_[it->second]->operator[](i);
}

void ColorCompressed::set_label(Index i, const SetStr &label) {
    assert(i < graph_size_);

    for (const auto &value : label) {
        add_label(i, value);
    }
}

void ColorCompressed::add_labels(const std::string &sequence, const SetStr &labels) {
    if (!graph_) {
        throw std::runtime_error("Please initialize with a graph\n");
    }
    graph_->align(sequence, [&](uint64_t i) {
        if (i > 0) {
            for (const auto &label : labels) {
                add_label(i, label);
            }
        }
    });
}

void ColorCompressed::add_label(Index i, const std::string &label) {
    assert(i < graph_size_);

    auto id_it = label_to_id_.find(label);

    if (id_it == label_to_id_.end()) {
        // current label does not exist yet -> assign new column
        if (annotation_curr_) {
            flush();
            delete annotation_curr_;
        }

        curr_id_ = id_to_label_.size();
        annotation_curr_ = new sdsl::bit_vector(graph_size_, 0);

        label_to_id_[label] = curr_id_;
        id_to_label_.push_back(label);
    } else if (!annotation_curr_ || curr_id_ != id_it->second) {
        if (annotation_curr_) {
            flush();
            delete annotation_curr_;
        }

        // decompress existing column
        curr_id_ = id_it->second;
        assert(bitmatrix_[curr_id_]);
        annotation_curr_ = inflate_column(curr_id_);
    }

    (*annotation_curr_)[i] = 1;
}

// write annotation to disk
void ColorCompressed::serialize(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();

    std::ofstream outstream(filename);
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, graph_size_);
    serialize_string_number_map(outstream, label_to_id_);
    libmaus2::util::StringSerialisation::serialiseStringVector(outstream, id_to_label_);
    for (auto *column : bitmatrix_) {
        assert(column);
        column->serialize(outstream);
    }
}

// read annotation from disk
bool ColorCompressed::load(const std::string &filename) {
    release();

    std::ifstream instream(filename);
    if (!instream.good())
        return false;

    graph_size_ = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
    label_to_id_ = load_string_number_map(instream);
    id_to_label_ = libmaus2::util::StringSerialisation::deserialiseStringVector(instream);
    bitmatrix_.resize(id_to_label_.size(), NULL);
    for (auto &column : bitmatrix_) {
        column = new sdsl::sd_vector<>();
        column->load(instream);
    }
    return true;
}

void ColorCompressed::release() {
    label_to_id_.clear();
    id_to_label_.clear();
    for (auto *column : bitmatrix_) {
        if (column)
            delete column;
    }
    bitmatrix_.clear();
    if (annotation_curr_) {
        delete annotation_curr_;
        annotation_curr_ = NULL;
    }
}

void ColorCompressed::flush() {
    if (!annotation_curr_)
        return;

    if (bitmatrix_.size() <= curr_id_) {
        bitmatrix_.push_back(NULL);
    }
    assert(curr_id_ < bitmatrix_.size());

    if (bitmatrix_[curr_id_]) {
        sdsl::bit_vector *tmp = inflate_column(curr_id_);
        *annotation_curr_ |= *tmp;
        delete tmp;
        delete bitmatrix_[curr_id_];
    }

    bitmatrix_[curr_id_] = new sdsl::sd_vector<>(*annotation_curr_);
}

sdsl::bit_vector* ColorCompressed::inflate_column(uint32_t id) const {
    auto *col = bitmatrix_.at(id);

    assert(col);

    sdsl::bit_vector *result = new sdsl::bit_vector(col->size(), 0);

    auto slct = sdsl::select_support_sd<>(col);
    auto rank = sdsl::rank_support_sd<>(col);
    size_t num_set_bits = rank(col->size());

    for (size_t i = 1; i <= num_set_bits; ++i) {
        size_t idx = slct(i);
        if (idx >= col->size())
            break;
        (*result)[idx] = 1;
    }

    return result;
}


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


AnnotationCategoryBloom::SetStr
AnnotationCategoryBloom::get(Index i) const {
    auto annotation = hash_annotate::BloomAnnotator::unpack(annotator_.get_annotation(i));
    SetStr result;
    for (size_t value : annotation) {
        result.insert(column_to_label_[value]);
    }
    return result;
}

void AnnotationCategoryBloom::set_label(Index i, const SetStr &label) {
    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);

    for (const auto &value : label) {
        add_label(kmer_edge, value);
    }
}

void AnnotationCategoryBloom::add_labels(const std::string &sequence,
                                         const SetStr &labels,
                                         size_t num_elements) {
    for (const auto &label : labels) {
        add_label(sequence, label, num_elements);
    }
}

void AnnotationCategoryBloom::add_label(const std::string &sequence,
                                        const std::string &label,
                                        size_t num_elements) {
    //TODO: set size of the Bloom filter based on the number of edges in graph
    if (label_to_column_.find(label) == label_to_column_.end()) {
        label_to_column_[label] = column_to_label_.size();
        column_to_label_.push_back(label);
    }
    annotator_.add_sequence(sequence, label_to_column_[label], num_elements);
}

bool AnnotationCategoryBloom::has_label(Index i, const SetStr &label) const {
    std::set<size_t> sorted_labels;
    for (const auto &value : label) {
        auto it = label_to_column_.find(value);
        if (it == label_to_column_.end())
            return false;
        sorted_labels.insert(it->second);
    }

    auto annotation = hash_annotate::BloomAnnotator::unpack(annotator_.get_annotation(i));
    return std::equal(sorted_labels.begin(), sorted_labels.end(),
                      annotation.begin(), annotation.end());
}

bool AnnotationCategoryBloom::load(const std::string &filename) {
    std::ifstream instream(filename);
    if (!instream.good())
        return false;

    label_to_column_ = load_string_number_map(instream);
    column_to_label_ = libmaus2::util::StringSerialisation::deserialiseStringVector(instream);
    annotator_.load(instream);
    return true;
}

void AnnotationCategoryBloom::serialize(const std::string &filename) const {
    std::ofstream outstream(filename);
    serialize_string_number_map(outstream, label_to_column_);
    libmaus2::util::StringSerialisation::serialiseStringVector(outstream, column_to_label_);
    annotator_.serialize(outstream);
}


AnnotationCategoryHash::AnnotationCategoryHash(const DBG_succ &graph)
      : graph_(graph), annotator_(graph_) {}

AnnotationCategoryHash::SetStr
AnnotationCategoryHash::get(Index i) const {
    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);

    auto annotation = annotator_.annotation_from_kmer(kmer_edge);
    SetStr result;
    for (size_t value : annotation) {
        result.insert(column_to_label_[value]);
    }
    return result;
}

void AnnotationCategoryHash::set_label(Index i, const SetStr &label) {
    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);

    for (const auto &value : label) {
        add_label(kmer_edge, value);
    }
}

void AnnotationCategoryHash::add_label(const std::string &sequence,
                                        const std::string &label) {
    if (label_to_column_.find(label) == label_to_column_.end()) {
        label_to_column_[label] = column_to_label_.size();
        column_to_label_.push_back(label);
    }
    annotator_.add_sequence(sequence, label_to_column_[label]);
}

bool AnnotationCategoryHash::has_label(Index i, const SetStr &label) const {
    std::set<size_t> sorted_labels;
    for (const auto &value : label) {
        auto it = label_to_column_.find(value);
        if (it == label_to_column_.end())
            return false;
        sorted_labels.insert(it->second);
    }

    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);
    auto annotation = hash_annotate::BloomAnnotator::unpack(
        annotator_.annotation_from_kmer(kmer_edge)
    );
    return std::equal(sorted_labels.begin(), sorted_labels.end(),
                      annotation.begin(), annotation.end());
}

void AnnotationCategoryHash::compare_annotations(const AnnotationCategoryBloom &bloom,
                                                 size_t step) const {
    bloom.compare_annotations(*this, step);
}

void AnnotationCategoryHash::compare_annotations(const hash_annotate::BloomAnnotator &bloom_annotator,
                                                 size_t step) const {
    bloom_annotator.test_fp_all(annotator_, step);
}

} // namespace annotate
