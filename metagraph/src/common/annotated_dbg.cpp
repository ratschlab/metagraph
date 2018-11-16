#include "annotated_dbg.hpp"

#include "annotate_column_compressed.hpp"
#include "annotate_column_compressed_fast.hpp"
#include "annotate_row_compressed.hpp"


AnnotatedDBG::AnnotatedDBG(DBG_succ *dbg, Annotator *annotation, size_t num_threads)
      : graph_(dbg), annotator_(annotation),
        thread_pool_(num_threads > 1 ? num_threads : 0) {}

AnnotatedDBG::AnnotatedDBG(DBG_succ *dbg, size_t num_threads)
      : AnnotatedDBG(dbg, NULL, num_threads) {}

AnnotatedDBG::AnnotatedDBG(Annotator *annotation, size_t num_threads)
      : AnnotatedDBG(NULL, annotation, num_threads) {}

AnnotatedDBG::AnnotatedDBG(size_t num_threads)
      : AnnotatedDBG(NULL, NULL, num_threads) {}

void AnnotatedDBG::adjust_annotation(const bit_vector_dyn &inserted_edges) {
    std::vector<uint64_t> inserted_edge_idx;

    for (uint64_t j = 1; j <= inserted_edges.num_set_bits(); ++j) {
        inserted_edge_idx.push_back(inserted_edges.select1(j));
    }

    // update the list of edges that can be annotated
    if (annotation_mask_.get()) {
        assert(!annotator_.get()
                || annotator_->num_objects() == annotation_mask_->num_set_bits());

        auto edge_mask = annotation_mask_->convert_to<sdsl::bit_vector>();

        utils::insert_default_values(inserted_edge_idx, &edge_mask);

        for (auto index : inserted_edge_idx) {
            assert(!edge_mask[index]);
            edge_mask[index] = true;
        }
        annotation_mask_.reset(new bit_vector_stat(std::move(edge_mask)));
    }

    if (annotator_) {
        // pretend that we don't have annotator to pass the compatibility check
        std::unique_ptr<Annotator> annotator { annotator_.release() };
        // transform indexes of the inserved k-mers to the annotation format
        for (auto &value : inserted_edge_idx) {
            value = graph_to_anno_index(value);
        }
        annotator->insert_rows(inserted_edge_idx);
        // move the annotation back
        annotator_.swap(annotator);
    }

    assert(check_compatibility());
}

void AnnotatedDBG::annotate_sequence_thread_safe(std::string sequence,
                                                 std::vector<std::string> labels) {
    std::vector<uint64_t> indices;

    graph_->map_to_edges(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(mutex_);
    annotator_->add_labels(indices, labels);
}

void AnnotatedDBG::annotate_sequence(const std::string &sequence,
                                     const std::vector<std::string> &labels) {
    assert(graph_.get() && annotator_.get());
    assert(check_compatibility(true));

    if (sequence.size() <= graph_->get_k())
        return;

    thread_pool_.enqueue(
        [this](auto... args) {
            this->annotate_sequence_thread_safe(args...);
        },
        sequence, labels
    );
}

std::vector<std::string> AnnotatedDBG::get_labels(const std::string &sequence,
                                                  double presence_ratio) const {
    assert(graph_.get() && annotator_.get());
    assert(check_compatibility(true));

    std::vector<uint64_t> indices;
    size_t num_missing_kmers = 0;

    graph_->map_to_edges(sequence, [&](uint64_t i) {
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

std::vector<std::pair<std::string, size_t>>
AnnotatedDBG::get_top_labels(const std::string &sequence,
                             size_t num_top_labels) const {
    assert(graph_.get() && annotator_.get());
    assert(check_compatibility(true));

    std::vector<uint64_t> indices;

    graph_->map_to_edges(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    });

    return annotator_->get_top_labels(indices, num_top_labels);
}

uint64_t AnnotatedDBG::num_anno_rows() const {
    return annotation_mask_.get()
            ? annotation_mask_->num_set_bits()
            : graph_->num_edges();
}

uint64_t AnnotatedDBG::graph_to_anno_index(uint64_t kmer_index) const {
    assert(kmer_index);
    assert(check_compatibility(true));

    if (!annotation_mask_.get())
        return kmer_index - 1;

    assert(kmer_index < annotation_mask_->size());

    if (!(*annotation_mask_)[kmer_index])
        throw std::runtime_error("Error: trying to annotate k-mer"
                                 " that cannot be annotated");

    return annotation_mask_->rank1(kmer_index) - 1;
}

bool AnnotatedDBG::check_compatibility(bool verbose) const {
    if (annotation_mask_.get()) {
        if (graph_.get() && annotation_mask_->size()
                            != graph_->num_edges() + 1) {
            if (verbose)
                std::cerr << "Error: edge mask is not compatible with graph."
                          << std::endl;
            return false;
        }
        if (annotator_.get() && annotation_mask_->num_set_bits()
                                != annotator_->num_objects()) {
            if (verbose)
                std::cerr << "Error: edge mask is not compatible with annotation."
                          << std::endl;
            return false;
        }
    } else if (graph_.get() && annotator_.get()
                            && graph_->num_edges()
                                != annotator_->num_objects()) {
        if (verbose)
            std::cerr << "Error: graph and annotation are not compatible."
                      << std::endl;
        return false;
    }

    return true;
}

void AnnotatedDBG::initialize_annotation_mask(size_t num_threads,
                                              bool prune_redundant_dummy) {
    assert(graph_.get());

    annotation_mask_.reset();

    std::vector<bool> edge_mask(graph_->num_edges() + 1, 0);

    if (prune_redundant_dummy) {
        graph_->erase_redundant_dummy_edges(&edge_mask, num_threads);
    } else {
        graph_->mark_source_dummy_edges(&edge_mask, num_threads);
    }
    // exclude 0 as the dummy index that denotes not existing k-mers
    edge_mask[0] = true;

    graph_->mark_sink_dummy_edges(&edge_mask);

    for (size_t i = 0; i < edge_mask.size(); ++i) {
        edge_mask[i] = !edge_mask[i];
    }

    annotation_mask_.reset(new bit_vector_stat(std::move(edge_mask)));

    assert(check_compatibility(true));
}

bool AnnotatedDBG::load_annotation_mask(const std::string &filename_base) {
    std::ifstream instream(filename_base + kAnnotationMaskExtension);

    if (!instream.good())
        return false;

    // release the old vector
    annotation_mask_.reset();

    // load a new vector
    std::unique_ptr<bit_vector_stat> vector { new bit_vector_stat() };
    if (!vector->load(instream))
        return false;

    annotation_mask_.swap(vector);

    if (check_compatibility(true))
        return true;

    // clear the loaded mask if it's incompatible with graph or annotation
    annotation_mask_.reset();
    return false;
}

void AnnotatedDBG::serialize_annotation_mask(const std::string &filename_base) const {
    if (!annotation_mask_.get())
        throw std::runtime_error("Trying to serialize uninitialized annotation mask");

    std::ofstream outstream(filename_base + kAnnotationMaskExtension);
    if (!outstream.good())
        throw std::ios_base::failure("Can't write to file "
                                        + filename_base + kAnnotationMaskExtension);

    annotation_mask_->serialize(outstream);
}
