#include "annotated_dbg.hpp"

#include "annotate_column_compressed.hpp"
#include "annotate_column_compressed_fast.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_bloom_filter.hpp"

using annotate::AnnotationCategoryBloom;


AnnotatedDBG::AnnotatedDBG(DBG_succ *dbg, Annotator *annotation, size_t num_threads)
      : graph_(dbg), annotator_(annotation),
        thread_pool_(num_threads > 1 ? num_threads : 0) {}

AnnotatedDBG::AnnotatedDBG(Annotator *annotation, DBG_succ *dbg, size_t num_threads)
      : AnnotatedDBG(dbg, annotation, num_threads) {}

void AnnotatedDBG::adjust_annotation(const bit_vector_dyn &inserted_edges) {
    std::vector<uint64_t> inserted_edge_idx;

    for (uint64_t j = 1; j <= inserted_edges.num_set_bits(); ++j) {
        inserted_edge_idx.push_back(inserted_edges.select1(j));
    }

    annotator_->insert_rows(inserted_edge_idx);
}

void annotate_sequence_thread_safe(std::string sequence,
                                   std::vector<std::string> labels,
                                   const DBG_succ &graph,
                                   AnnotatedDBG::Annotator *annotator,
                                   std::mutex &annotation_mutex) {
    std::vector<uint64_t> indices;

    graph.map_to_edges(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(i);
    });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(annotation_mutex);
    annotator->add_labels(indices, labels);
}

void AnnotatedDBG::annotate_sequence(const std::string &sequence,
                                     const std::vector<std::string> &labels) {
    assert(graph_.get() && annotator_.get());
    assert(graph_->num_edges() + 1 == annotator_->num_objects());

    if (sequence.size() <= graph_->get_k())
        return;

    try {
        dynamic_cast<AnnotationCategoryBloom&>(*annotator_).add_labels(
            sequence, labels, graph_->num_edges()
        );
    } catch (const std::bad_cast &) {
        thread_pool_.enqueue(
            annotate_sequence_thread_safe,
            sequence, labels,
            std::ref(*graph_), annotator_.get(), std::ref(mutex_)
        );
    }
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(graph_.get() && annotator_.get());
    assert(graph_->num_edges() + 1 == annotator_->num_objects());

    std::vector<uint64_t> indices;
    size_t num_missing_kmers = 0;

    graph_->map_to_edges(sequence, [&](uint64_t i) {
        if (i > 0) {
            indices.push_back(i);
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

std::vector<std::pair<std::string, size_t>>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels) const {
    assert(graph_.get() && annotator_.get());
    assert(graph_->num_edges() + 1 == annotator_->num_objects());

    return annotator_->get_top_labels(
        graph_->map_to_edges(sequence),
        num_top_labels
    );
}
