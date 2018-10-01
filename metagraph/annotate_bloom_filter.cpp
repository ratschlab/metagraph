#include "annotate_bloom_filter.hpp"

#include "serialization.hpp"

using libmaus2::util::StringSerialisation;
using hash_annotate::BloomAnnotator;


namespace annotate {

void AnnotationCategoryBloom::set_labels(Index, const VLabels &) {
    throw std::runtime_error(
        "ERROR: can't erase from Bloom filter. Use add_labels instead."
    );
}

AnnotationCategoryBloom::VLabels
AnnotationCategoryBloom::get_labels(Index i) const {
    auto annotation = BloomAnnotator::unpack(annotator_.get_annotation(i));
    VLabels result;
    for (auto value : annotation) {
        result.push_back(column_to_label_[value]);
    }
    return result;
}

void AnnotationCategoryBloom::add_label(Index i, const Label &label) {
    add_label(graph_.get_node_kmer(i) + graph_.get_edge_label(i), label);
}

void AnnotationCategoryBloom::add_labels(Index i, const VLabels &labels) {
    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);

    for (const auto &label : labels) {
        add_label(kmer_edge, label);
    }
}

void AnnotationCategoryBloom::add_labels(const std::vector<Index> &indices,
                                         const VLabels &labels) {
    for (Index i : indices) {
        add_labels(i, labels);
    }
}

void AnnotationCategoryBloom::add_labels(const std::string &sequence,
                                         const VLabels &labels,
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

bool AnnotationCategoryBloom::has_label(Index i, const Label &label) const {
    auto it = label_to_column_.find(label);
    if (it == label_to_column_.end())
        return false;

    auto annotation = BloomAnnotator::unpack(annotator_.get_annotation(i));
    return std::find(annotation.begin(), annotation.end(), it->second)
            != annotation.end();
}

bool AnnotationCategoryBloom::has_labels(Index i, const VLabels &labels) const {
    std::set<size_t> indices;
    for (const auto &label : labels) {
        auto it = label_to_column_.find(label);
        if (it == label_to_column_.end())
            return false;
        indices.insert(it->second);
    }

    auto annotation = BloomAnnotator::unpack(annotator_.get_annotation(i));
    std::set<size_t> sorted_annotation(annotation.begin(), annotation.end());

    return std::includes(indices.begin(), indices.end(),
                         sorted_annotation.begin(), sorted_annotation.end());
}

bool AnnotationCategoryBloom::merge_load(const std::vector<std::string> &filenames) {
    std::ifstream instream(filenames.at(0) + ".bloom.annodbg");
    if (!instream.good())
        return false;

    label_to_column_ = load_string_number_map(instream);
    column_to_label_ = StringSerialisation::deserialiseStringVector(instream);
    annotator_.load(instream);
    return true;
}

void AnnotationCategoryBloom::serialize(const std::string &filename) const {
    std::ofstream outstream(filename + ".bloom.annodbg");
    serialize_string_number_map(outstream, label_to_column_);
    StringSerialisation::serialiseStringVector(outstream, column_to_label_);
    annotator_.serialize(outstream);
}

AnnotationCategoryBloom::VLabels
AnnotationCategoryBloom::get_labels(const std::vector<Index> &indices,
                                    double presence_ratio) const {
    assert(presence_ratio >= 0 && presence_ratio <= 1);

    const size_t min_labels_discovered =
                        presence_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * presence_ratio);
    // const size_t max_labels_missing = indices.size() - min_labels_discovered;

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        const auto &labels_encoded =
            BloomAnnotator::unpack(annotator_.get_annotation(i));

        for (auto code : labels_encoded) {
            encoded_counter[code]++;
        }
    }

    VLabels filtered_labels;

    for (auto it = encoded_counter.begin(); it != encoded_counter.end(); ++it) {
        if (it->second >= min_labels_discovered)
            filtered_labels.push_back(column_to_label_[it->first]);
    }

    return filtered_labels;
}

std::vector<std::pair<AnnotationCategoryBloom::Label, size_t>>
AnnotationCategoryBloom
::get_top_labels(const std::vector<Index> &indices, size_t num_top) const {
    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        const auto &labels_encoded =
            BloomAnnotator::unpack(annotator_.get_annotation(i));

        for (auto code : labels_encoded) {
            encoded_counter[code]++;
        }
    }

    std::vector<std::pair<size_t, size_t>> counts(
        encoded_counter.begin(), encoded_counter.end()
    );
    // sort in decreasing order
    std::sort(counts.begin(), counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    counts.resize(std::min(counts.size(), num_top));

    std::vector<std::pair<Label, size_t>> top_counts;
    for (const auto &encoded_pair : counts) {
        top_counts.emplace_back(column_to_label_[encoded_pair.first],
                                encoded_pair.second);
    }

    return top_counts;
}

uint64_t AnnotationCategoryBloom::num_objects() const {
    return graph_.get_num_edges();
}

size_t AnnotationCategoryBloom::num_labels() const {
    return column_to_label_.size();
}

double AnnotationCategoryBloom::sparsity() const {
    throw std::runtime_error("Unknown sparsity, method not implemented");
}

} // namespace annotate
