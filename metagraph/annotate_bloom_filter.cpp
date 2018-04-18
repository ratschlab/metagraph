#include "annotate_bloom_filter.hpp"

#include "serialization.hpp"

using libmaus2::util::StringSerialisation;
using hash_annotate::BloomAnnotator;


namespace annotate {

void AnnotationCategoryBloom::set_coloring(Index i, const Coloring &coloring) {
    std::ignore = i;
    std::ignore = coloring;
    throw std::runtime_error(
        "ERROR: can't erase from Bloom filter. Use add_colors instead."
    );
}

AnnotationCategoryBloom::Coloring
AnnotationCategoryBloom::get_coloring(Index i) const {
    auto annotation = BloomAnnotator::unpack(annotator_.get_annotation(i));
    Coloring result;
    for (auto value : annotation) {
        result.push_back(column_to_label_[value]);
    }
    return result;
}

void AnnotationCategoryBloom::add_color(Index i, const Color &color) {
    add_color(graph_.get_node_kmer(i) + graph_.get_edge_label(i), color);
}

void AnnotationCategoryBloom::add_colors(Index i, const Coloring &colors) {
    auto kmer_edge = graph_.get_node_kmer(i) + graph_.get_edge_label(i);

    for (const auto &color : colors) {
        add_color(kmer_edge, color);
    }
}

void AnnotationCategoryBloom::add_colors(const std::vector<Index> &indices,
                                             const Coloring &colors) {
    for (Index i : indices) {
        add_colors(i, colors);
    }
}

void AnnotationCategoryBloom::add_colors(const std::string &sequence,
                                             const Coloring &colors,
                                             size_t num_elements) {
    for (const auto &color : colors) {
        add_color(sequence, color, num_elements);
    }
}

void AnnotationCategoryBloom::add_color(const std::string &sequence,
                                            const std::string &color,
                                            size_t num_elements) {
    //TODO: set size of the Bloom filter based on the number of edges in graph
    if (label_to_column_.find(color) == label_to_column_.end()) {
        label_to_column_[color] = column_to_label_.size();
        column_to_label_.push_back(color);
    }
    annotator_.add_sequence(sequence, label_to_column_[color], num_elements);
}

bool AnnotationCategoryBloom::has_color(Index i, const Color &color) const {
    auto it = label_to_column_.find(color);
    if (it == label_to_column_.end())
        return false;

    auto annotation = BloomAnnotator::unpack(annotator_.get_annotation(i));
    return std::find(annotation.begin(), annotation.end(), it->second)
            != annotation.end();
}

bool AnnotationCategoryBloom::has_colors(Index i, const Coloring &coloring) const {
    std::set<size_t> indices;
    for (const auto &color : coloring) {
        auto it = label_to_column_.find(color);
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

AnnotationCategoryBloom::Coloring
AnnotationCategoryBloom::aggregate_colors(const std::vector<Index> &indices,
                                          double discovery_ratio) const {
    assert(discovery_ratio >= 0 && discovery_ratio <= 1);

    const size_t min_colors_discovered =
                        discovery_ratio == 0
                            ? 1
                            : std::ceil(indices.size() * discovery_ratio);
    // const size_t max_colors_missing = indices.size() - min_colors_discovered;

    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        const auto &coloring_encoded =
            BloomAnnotator::unpack(annotator_.get_annotation(i));

        for (auto code : coloring_encoded) {
            encoded_counter[code]++;
        }
    }

    Coloring filtered_colors;

    for (auto it = encoded_counter.begin(); it != encoded_counter.end(); ++it) {
        if (it->second >= min_colors_discovered)
            filtered_colors.push_back(column_to_label_[it->first]);
    }

    return filtered_colors;
}

std::vector<std::pair<AnnotationCategoryBloom::Color, size_t>>
AnnotationCategoryBloom
::get_most_frequent_colors(const std::vector<Index> &indices,
                           size_t num_top) const {
    std::unordered_map<size_t, size_t> encoded_counter;

    for (Index i : indices) {
        const auto &coloring_encoded =
            BloomAnnotator::unpack(annotator_.get_annotation(i));

        for (auto code : coloring_encoded) {
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

    std::vector<std::pair<Color, size_t>> top_counts;
    for (const auto &encoded_pair : counts) {
        top_counts.emplace_back(column_to_label_[encoded_pair.first],
                                encoded_pair.second);
    }

    return top_counts;
}

} // namespace annotate
