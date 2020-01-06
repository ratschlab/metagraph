#include "annotated_dbg.hpp"

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"

typedef std::pair<std::string, size_t> StringCountPair;


AnnotatedDBG::AnnotatedDBG(std::shared_ptr<SequenceGraph> dbg,
                           std::unique_ptr<Annotator>&& annotation,
                           size_t num_threads,
                           bool force_fast)
      : graph_(dbg), annotator_(std::move(annotation)),
        thread_pool_(num_threads > 1 ? num_threads : 0),
        force_fast_(force_fast) {
    assert(graph_.get());
    assert(annotator_.get());
}

void AnnotatedDBG::annotate_sequence_thread_safe(const std::string &sequence,
                                                 const std::vector<std::string> &labels) {
    std::vector<uint64_t> indices;
    indices.reserve(sequence.size());

    graph_->map_to_nodes(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(mutex_);

    if (force_fast_) {
        auto row_major = dynamic_cast<annotate::RowCompressed<std::string>*>(annotator_.get());
        if (row_major) {
            row_major->add_labels_fast(indices, labels);
            return;
        }
    }

    annotator_->add_labels(indices, labels);
}

void AnnotatedDBG::annotate_sequence(std::string&& sequence,
                                     const std::vector<std::string> &labels) {
    assert(check_compatibility());

    thread_pool_.enqueue(
        [this](const auto&... args) {
            this->annotate_sequence_thread_safe(args...);
        },
        std::move(sequence), labels
    );
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    tsl::hopscotch_map<uint64_t, size_t> index_counts;
    size_t num_present_kmers = 0;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](uint64_t i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        } else {
            num_missing_kmers++;
        }
    });

    size_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                 * (num_present_kmers
                                                     + num_missing_kmers)));

    if (num_present_kmers < min_count)
        return {};

    return get_labels(index_counts, min_count);
}

std::vector<std::string> AnnotatedDBG
::get_labels(const std::vector<std::string> &sequences,
             const std::vector<double> &weights,
             double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());
    assert(sequences.size() == weights.size());

    tsl::hopscotch_map<uint64_t, double> index_weights;
    double weighted_num_missing_kmers = 0;
    double scale = std::accumulate(weights.begin(), weights.end(), 0.0);

    for (size_t j = 0; j < sequences.size(); ++j) {
        graph_->map_to_nodes(sequences[j], [&](uint64_t i) {
            if (i > 0) {
                index_weights[graph_to_anno_index(i)] += weights[j];
            } else {
                weighted_num_missing_kmers += weights[j];
            }
        });
    }

    tsl::hopscotch_map<uint64_t, size_t> index_counts;

    size_t num_missing_kmers = std::floor(weighted_num_missing_kmers / scale);
    size_t num_present_kmers = 0;

    for (const auto &pair : index_weights) {
        double count = std::floor(pair.second / scale);
        index_counts[pair.first] = count;
        num_present_kmers += count;
    }

    size_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                 * (num_present_kmers
                                                     + num_missing_kmers)));

    if (num_present_kmers < min_count)
        return {};

    return get_labels(index_counts, min_count);
}

std::vector<std::string>
AnnotatedDBG::get_labels(const tsl::hopscotch_map<row_index, size_t> &index_counts,
                         size_t min_count) const {
    assert(check_compatibility());

    auto code_counts = annotator_->count_labels(index_counts, min_count, min_count);

    std::vector<std::string> labels;
    labels.reserve(code_counts.size());

    const auto &label_encoder = annotator_->get_label_encoder();

    for (const auto &pair : code_counts) {
        assert(pair.second >= min_count);
        labels.push_back(label_encoder.decode(pair.first));
    }

    return labels;
}

std::vector<std::string> AnnotatedDBG::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get(graph_to_anno_index(index));
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels,
                             double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    tsl::hopscotch_map<uint64_t, size_t> index_counts;
    size_t num_present_kmers = 0;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](uint64_t i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        } else {
            num_missing_kmers++;
        }
    });

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                    * (num_present_kmers
                                                        + num_missing_kmers)));
    if (num_present_kmers < min_count)
        return {};

    return get_top_labels(index_counts, num_top_labels, min_count);
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::vector<std::string> &sequences,
                             const std::vector<double> &weights,
                             size_t num_top_labels,
                             double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());
    assert(sequences.size() == weights.size());

    tsl::hopscotch_map<uint64_t, double> index_weights;
    double weighted_num_missing_kmers = 0;
    double scale = std::accumulate(weights.begin(), weights.end(), 0.0);

    for (size_t j = 0; j < sequences.size(); ++j) {
        graph_->map_to_nodes(sequences[j], [&](uint64_t i) {
            if (i > 0) {
                index_weights[graph_to_anno_index(i)] += weights[j];
            } else {
                weighted_num_missing_kmers += weights[j];
            }
        });
    }

    tsl::hopscotch_map<uint64_t, size_t> index_counts;

    size_t num_missing_kmers = std::floor(weighted_num_missing_kmers / scale);
    size_t num_present_kmers = 0;

    for (const auto &pair : index_weights) {
        double count = std::floor(pair.second / scale);
        index_counts[pair.first] = count;
        num_present_kmers += count;
    }

    size_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                 * (num_present_kmers
                                                     + num_missing_kmers)));

    if (num_present_kmers < min_count)
        return {};

    return get_top_labels(index_counts, num_top_labels, min_count);
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const tsl::hopscotch_map<node_index, size_t> &index_counts,
                             size_t num_top_labels,
                             size_t min_count) const {
    assert(check_compatibility());

    auto code_counts = annotator_->count_labels(index_counts, min_count);

    assert(std::all_of(
        code_counts.begin(), code_counts.end(),
        [&](const auto &code_count) { return code_count.second >= min_count; }
    ));

    std::sort(code_counts.begin(), code_counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    // leave only first |num_top_labels| top labels
    if (code_counts.size() > num_top_labels)
        code_counts.erase(code_counts.begin() + num_top_labels,
                           code_counts.end());

    const auto &label_encoder = annotator_->get_label_encoder();

    std::vector<StringCountPair> label_counts(code_counts.size());

    for (size_t i = 0; i < code_counts.size(); ++i) {
        label_counts[i].first = label_encoder.decode(code_counts[i].first);
        label_counts[i].second = code_counts[i].second;
    }

    return label_counts;
}

bool AnnotatedDBG::label_exists(const std::string &label) const {
    return annotator_->label_exists(label);
}

bool AnnotatedDBG::has_label(node_index index, const std::string &label) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->has_label(graph_to_anno_index(index), label);
}

void AnnotatedDBG
::call_annotated_nodes(const std::string &label,
                       std::function<void(node_index)> callback) const {
    assert(check_compatibility());

    annotator_->call_objects(
        label,
        [&](row_index index) { callback(anno_to_graph_index(index)); }
    );
}

AnnotatedDBG::row_index
AnnotatedDBG::graph_to_anno_index(node_index kmer_index) {
    assert(kmer_index);
    return kmer_index - 1;
}

AnnotatedDBG::node_index
AnnotatedDBG::anno_to_graph_index(row_index anno_index) {
    return anno_index + 1;
}

bool AnnotatedDBG::check_compatibility() const {
    return graph_->max_index() == annotator_->num_objects();
}
