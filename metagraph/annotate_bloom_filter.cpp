#include "annotate_bloom_filter.hpp"

#include "serialization.hpp"


namespace annotate {

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
