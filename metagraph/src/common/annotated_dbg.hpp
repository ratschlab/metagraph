#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include "sequence_graph.hpp"
#include "annotate.hpp"
#include "utils.hpp"


class AnnotatedDBG {
  public:
    typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

    AnnotatedDBG(SequenceGraph *dbg,
                 Annotator *annotation,
                 size_t num_threads = 0);

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string>
    get_labels(const std::string &sequence, double presence_ratio) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence, size_t num_top_labels) const;

    void annotate_sequence(const std::string &sequence,
                           const std::vector<std::string> &labels);

    void join() { thread_pool_.join(); }

    uint64_t num_anno_rows() const;

    bool check_compatibility() const;

    const SequenceGraph& get_graph() const { return *graph_; }
    const Annotator& get_annotation() const { return *annotator_; }

    static void insert_zero_rows(Annotator *annotator,
                                 const bit_vector_dyn &inserted_edges);
    static Annotator::Index
    graph_to_anno_index(SequenceGraph::node_index kmer_index);

  private:
    void annotate_sequence_thread_safe(std::string sequence,
                                       std::vector<std::string> labels);

    std::unique_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    utils::ThreadPool thread_pool_;
    std::mutex mutex_;
};

#endif // __ANNOTATED_DBG_HPP__
