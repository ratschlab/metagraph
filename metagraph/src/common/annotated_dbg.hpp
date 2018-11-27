#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include "sequence_graph.hpp"
#include "annotate.hpp"
#include "utils.hpp"


class AnnotatedDBG {
  public:
    typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

    explicit AnnotatedDBG(size_t num_threads = 0);
    AnnotatedDBG(SequenceGraph *dbg, size_t num_threads = 0);
    AnnotatedDBG(Annotator *annotation, size_t num_threads = 0);
    AnnotatedDBG(SequenceGraph *dbg, Annotator *annotation, size_t num_threads = 0);

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string>
    get_labels(const std::string &sequence, double presence_ratio) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence, size_t num_top_labels) const;

    void adjust_annotation(const bit_vector_dyn &inserted_edges);

    void annotate_sequence(const std::string &sequence,
                           const std::vector<std::string> &labels);

    void initialize_annotation_mask(std::vector<bool>&& skipped_nodes);
    bool load_annotation_mask(const std::string &filename_base);
    void serialize_annotation_mask(const std::string &filename_base) const;

    void join() { thread_pool_.join(); }

    void set_annotation(Annotator *annotator) { annotator_.reset(annotator); }
    Annotator& get_annotation() { return *annotator_; }
    const Annotator& get_annotation() const { return *annotator_; }
    uint64_t num_anno_rows() const;

    SequenceGraph* release_graph() { return graph_.release(); }
    SequenceGraph& get_graph() { return *graph_; }
    const SequenceGraph& get_graph() const { return *graph_; }

    bool check_compatibility(bool verbose = false) const;

    uint64_t graph_to_anno_index(uint64_t kmer_index) const;

  private:
    void annotate_sequence_thread_safe(std::string sequence,
                                       std::vector<std::string> labels);

    std::unique_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;
    // marks graph edges that can be annotated
    std::unique_ptr<bit_vector_stat> annotation_mask_;

    utils::ThreadPool thread_pool_;
    std::mutex mutex_;

    static constexpr auto kAnnotationMaskExtension = ".edgemask";
};

#endif // __ANNOTATED_DBG_HPP__
