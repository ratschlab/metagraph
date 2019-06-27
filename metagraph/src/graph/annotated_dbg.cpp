#include "annotated_dbg.hpp"

#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"


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

void AnnotatedDBG::annotate_sequence_thread_safe(std::string sequence,
                                                 std::vector<std::string> labels) {
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

std::vector<std::string> AnnotatedDBG::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get_labels(graph_to_anno_index(index));
}

std::vector<std::pair<std::string, size_t>>
AnnotatedDBG
::get_labels(const std::string &sequence,
             std::function<bool(uint64_t, uint64_t)> start,
             std::function<bool(const std::vector<std::pair<std::string,
                                                            size_t>>&)> terminate) const {
    assert(check_compatibility());

    std::unordered_map<uint64_t, uint64_t> index_counts;
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

    if (!start(num_present_kmers, num_missing_kmers))
        return {};

    std::vector<uint64_t> indices;
    indices.reserve(index_counts.size());
    std::transform(index_counts.begin(), index_counts.end(),
                   std::back_inserter(indices),
                   [](const auto &pair) { return pair.first; });

    std::vector<std::pair<std::string, size_t>> label_counts;
    label_counts.reserve(annotator_->num_labels());

    if (dynamic_cast<annotate::ColumnCompressed<>*>(annotator_.get())) {
        // Iterate by column instead of by row for column-major annotators
        annotator_->call_labels([&](const auto &label) {
            label_counts.emplace_back(label, 0);
            annotator_->call_indices_until(
                indices,
                label,
                [&](auto i) {
                    auto it = index_counts.find(i);
                    assert(it != index_counts.end());
                    label_counts.back().second += it->second;
                },
                [&]() { return terminate(label_counts); }
            );
        });
    } else {
        annotator_->call_labels([&](const auto &label) {
            label_counts.emplace_back(label, 0);
        });

        auto it = index_counts.begin();
        annotator_->call_rows(
            indices,
            [&](auto&& label_indices) {
                assert(it != index_counts.end());

                for (auto j : label_indices) {
                    label_counts[j].second += it->second;
                }

                ++it;
            },
            [&]() { return terminate(label_counts); }
        );
    }

    return label_counts;
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);

    uint64_t min_count = 0;

    std::function<bool(const std::vector<std::pair<std::string, size_t>>&)> terminate;
    if (dynamic_cast<annotate::ColumnCompressed<>*>(annotator_.get())) {
        terminate = [&](const auto &label_counts) {
            return label_counts.back().second >= min_count;
        };
    } else {
        terminate = [&](const auto &label_counts) {
            assert(label_counts.size() == annotator_->num_labels());
            return std::all_of(label_counts.begin(), label_counts.end(),
                               [&](const auto &pair) {
                                   return pair.second >= min_count;
                               });
        };
    }

    auto label_counts = get_labels(
        sequence,
        [&](auto num_present_kmers, auto num_missing_kmers) {
            min_count = std::max(1.0,
                                 std::ceil(presence_ratio
                                     * (num_present_kmers + num_missing_kmers)));
            return num_present_kmers >= min_count;
        },
        terminate
    );

    std::vector<std::string> labels;
    for (auto&& pair : label_counts) {
        if (pair.second >= min_count)
            labels.emplace_back(std::move(pair.first));
    }

    return labels;
}

std::vector<std::pair<std::string, size_t>>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels,
                             double min_label_frequency) const {
    assert(min_label_frequency >= 0.);
    assert(min_label_frequency <= 1.);

    if (!num_top_labels)
        return {};

    uint64_t min_count = 0;
    auto label_counts = get_labels(
        sequence,
        [&](auto num_present_kmers, auto num_missing_kmers) {
            min_count = std::max(1.0,
                                 std::ceil(min_label_frequency
                                     * (num_present_kmers + num_missing_kmers)));
            return num_present_kmers >= min_count;
        }
    );

    std::sort(label_counts.begin(), label_counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    // remove labels which don't meet num_top_labels and min_label_frequency criteria
    label_counts.erase(std::find_if(label_counts.begin(),
                                    num_top_labels < label_counts.size()
                                        ? label_counts.begin() + num_top_labels
                                        : label_counts.end(),
                                    [&](const auto &pair) {
                                        return pair.second < min_count;
                                    }),
                       label_counts.end());
    label_counts.shrink_to_fit();

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
