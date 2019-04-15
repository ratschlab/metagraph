#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <memory>

#include "sequence_graph.hpp"
#include "annotate.hpp"
#include "threading.hpp"


class AnnotatedDBG {
  public:
    typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

    AnnotatedDBG(std::shared_ptr<SequenceGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 size_t num_threads = 0,
                 bool force_fast = false);

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string>
    get_labels(const std::string &sequence, double presence_ratio = 0) const;

    std::vector<std::string> get_labels(const SequenceGraph::node_index &index) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence,
                   size_t num_top_labels,
                   double min_label_frequency = 0.0) const;

    bool label_exists(const std::string &label) const;

    bool has_label(const SequenceGraph::node_index &index, const std::string &label) const;

    void call_indices(const std::string &label,
                      const std::function<void(const SequenceGraph::node_index&)> callback) const;

    uint64_t count_labels(SequenceGraph::node_index index, const std::vector<std::string> &labels_to_match) const;

    void annotate_sequence(const std::string &sequence,
                           const std::vector<std::string> &labels);

    void join() { thread_pool_.join(); }

    uint64_t num_anno_rows() const;

    uint64_t num_nodes() const { return graph_->num_nodes(); }

    bool check_compatibility() const;

    const SequenceGraph& get_graph() const { return *graph_; }
    const Annotator& get_annotation() const { return *annotator_; }

    static void insert_zero_rows(Annotator *annotator,
                                 const bit_vector_dyn &inserted_edges);

  private:
    static Annotator::Index
    graph_to_anno_index(SequenceGraph::node_index kmer_index);

    static SequenceGraph::node_index
    anno_to_graph_index(Annotator::Index anno_index);

    void annotate_sequence_thread_safe(std::string sequence,
                                       std::vector<std::string> labels);

    std::shared_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    ThreadPool thread_pool_;
    std::mutex mutex_;
    bool force_fast_;
};

#endif // __ANNOTATED_DBG_HPP__
