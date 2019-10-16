#include "annotated_dbg.hpp"

#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"

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

void AnnotatedDBG::insert_zero_rows(Annotator *annotator,
                                    const bit_vector_dyn &inserted_edges) {
    assert(annotator);

    std::vector<uint64_t> inserted_edge_idx;

    // transform indexes of the inserved k-mers to the annotation format
    inserted_edges.call_ones([&](auto i) {
        inserted_edge_idx.push_back(graph_to_anno_index(i));
    });

    annotator->insert_rows(inserted_edge_idx);
}

void AnnotatedDBG::annotate_sequence_thread_safe(const std::string &sequence,
                                                 const std::vector<std::string> &labels) {
    std::vector<uint64_t> indices;

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

void AnnotatedDBG::annotate_sequence(const std::string &sequence,
                                     const std::vector<std::string> &labels) {
    assert(check_compatibility());

    thread_pool_.enqueue(
        [this](auto... args) {
            this->annotate_sequence_thread_safe(args...);
        },
        sequence, labels
    );
}

std::vector<StringCountPair>
count_labels(const annotate::ColumnCompressed<> &annotation,
             const std::unordered_map<uint64_t, size_t> &index_counts,
             std::function<bool(size_t /* checked */,
                                size_t /* matched */)> stop_counting_for_label) {

    std::vector<uint64_t> indices(index_counts.size());
    std::transform(index_counts.begin(), index_counts.end(), indices.begin(),
                   [](const auto &pair) { return pair.first; });

    std::vector<StringCountPair> label_counts;
    label_counts.reserve(annotation.num_labels());

    for (const auto &label : annotation.get_all_labels()) {
        uint64_t total_checked = 0;
        uint64_t total_matched = 0;
        annotation.call_relations(
            indices,
            label,
            [&](auto i, bool matched) {
                auto count = index_counts.at(i);
                total_checked += count;
                total_matched += count * matched;
            },
            [&]() { return stop_counting_for_label(total_checked, total_matched); }
        );

        if (total_matched)
            label_counts.emplace_back(label, total_matched);
    }

    return label_counts;
}

std::vector<StringCountPair>
count_labels(const AnnotatedDBG::Annotator &annotation,
             const std::unordered_map<uint64_t, size_t> &index_counts,
             std::function<bool(size_t /* checked */,
                                size_t /* min_matched */,
                                size_t /* max_matched */)> /* stop_counting_labels */) {

    std::vector<uint64_t> indices(index_counts.size());
    std::transform(index_counts.begin(), index_counts.end(), indices.begin(),
                   [](const auto &pair) { return pair.first; });

    std::vector<size_t> code_counts(annotation.num_labels(), 0);

    auto it = index_counts.begin();
    annotation.call_rows(
        indices,
        [&](const auto &row) {
            assert(it != index_counts.end());

            for (size_t label_code : row) {
                assert(label_code < code_counts.size());

                code_counts[label_code] += it->second;
            }

            ++it;
        }
        //TODO: implement the support for calling the smallest and the largest label count
        // [&]() { return stop_counting_labels(checked, min_matched, max_matched); }
        //
        // Use unordered map: count -> number of labels with that count (histogram)
        // and maintain a pointer to the smallest count with positive number of labels.
        // The order of pointers is a non decreasing sequence, so the complexity is O(\sum_i count_i)
    );

    std::vector<StringCountPair> label_counts;
    label_counts.reserve(annotation.num_labels());

    for (size_t label_code = 0; label_code < code_counts.size(); ++label_code) {
        if (code_counts[label_code]) {
            label_counts.emplace_back(annotation.get_label_encoder().decode(label_code),
                                      code_counts[label_code]);
        }
    }

    return label_counts;
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    std::unordered_map<uint64_t, size_t> index_counts;
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
::get_labels(const std::unordered_map<row_index, size_t> &index_counts,
             size_t min_count) const {
    assert(check_compatibility());

    uint64_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<StringCountPair> label_counts
        = dynamic_cast<const annotate::ColumnCompressed<>*>(annotator_.get())
            // Iterate by column instead of by row for column-major annotators
            ? count_labels(dynamic_cast<const annotate::ColumnCompressed<>&>(*annotator_),
                           index_counts,
                           [&](size_t checked, size_t matched) {
                               return matched >= min_count
                                        || matched + (total_sum_count - checked) < min_count;
                           })
            : count_labels(*annotator_,
                           index_counts,
                           [&](size_t checked, size_t min_matched, size_t max_matched) {
                               return min_matched >= min_count
                                        || max_matched + (total_sum_count - checked) < min_count;
                           });

    std::vector<std::string> labels;
    for (auto&& pair : label_counts) {
        if (pair.second >= min_count)
            labels.emplace_back(std::move(pair.first));
    }

    return labels;
}

std::vector<std::string> AnnotatedDBG::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get_labels(graph_to_anno_index(index));
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels,
                             double min_label_frequency) const {
    assert(min_label_frequency >= 0.);
    assert(min_label_frequency <= 1.);
    assert(check_compatibility());

    std::unordered_map<uint64_t, size_t> index_counts;
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

    uint64_t min_count = std::max(1.0, std::ceil(min_label_frequency
                                                    * (num_present_kmers
                                                        + num_missing_kmers)));
    if (num_present_kmers < min_count)
        return {};

    return get_top_labels(index_counts, num_top_labels, min_count);
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::unordered_map<node_index, size_t> &index_counts,
                             size_t num_top_labels,
                             size_t min_count) const {
    assert(check_compatibility());

    uint64_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<StringCountPair> label_counts
        = dynamic_cast<const annotate::ColumnCompressed<>*>(annotator_.get())
            // Iterate by column instead of by row for column-major annotators
            ? count_labels(dynamic_cast<const annotate::ColumnCompressed<>&>(*annotator_),
                           index_counts,
                           [&](size_t checked, size_t matched) {
                               return matched + (total_sum_count - checked) < min_count;
                           })
            : count_labels(*annotator_,
                           index_counts,
                           [&](size_t checked, size_t /* min_matched */, size_t max_matched) {
                               return max_matched + (total_sum_count - checked) < min_count;
                           });

    // remove labels which don't meet |min_label_frequency| criterion
    label_counts.erase(
        std::remove_if(label_counts.begin(), label_counts.end(),
                       [&](const auto &pair) { return pair.second < min_count; }),
        label_counts.end()
    );

    std::sort(label_counts.begin(), label_counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    // leave only first |num_top_labels| top labels
    if (label_counts.size() > num_top_labels)
        label_counts.erase(label_counts.begin() + num_top_labels,
                           label_counts.end());
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
    return graph_->num_nodes() == annotator_->num_objects();
}
