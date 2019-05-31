#include "annotated_dbg.hpp"

#include "annotate_row_compressed.hpp"


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

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(presence_ratio >= 0);
    assert(check_compatibility());

    std::vector<uint64_t> indices;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](uint64_t i) {
        if (i > 0) {
            indices.push_back(graph_to_anno_index(i));
        } else {
            num_missing_kmers++;
        }
    });

    // presence_ratio = num_present / all
    // new_presence_ratio = num_present / (all - missing)
    if (indices.size() > 0
        && indices.size()
            >= presence_ratio * (indices.size() + num_missing_kmers)) {
        return annotator_->get_labels(
            indices,
            presence_ratio
                * (indices.size() + num_missing_kmers)
                / indices.size()
        );
    } else {
        return {};
    }
}

std::vector<std::string> AnnotatedDBG::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get_labels(graph_to_anno_index(index));
}

std::vector<std::pair<std::string, size_t>>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels,
                             double min_label_frequency) const {
    assert(min_label_frequency >= 0);
    assert(check_compatibility());

    std::vector<uint64_t> indices;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](uint64_t i) {
        if (i > 0) {
            indices.push_back(graph_to_anno_index(i));
        } else {
            num_missing_kmers++;
        }
    });

    if (indices.size() > 0
        && indices.size()
            >= min_label_frequency * (indices.size() + num_missing_kmers)) {

        return annotator_->get_top_labels(indices, num_top_labels,
            min_label_frequency
                * (indices.size() + num_missing_kmers)
                / indices.size()
        );
    } else {
        return {};
    }
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
