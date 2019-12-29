#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <memory>

#include "sequence_graph.hpp"
#include "annotate.hpp"
#include "threading.hpp"


// TODO: rename to AnnotatedSequenceGraph
class AnnotatedDBG {
  public:
    typedef annotate::MultiLabelEncoded<uint64_t, std::string> Annotator;
    using node_index = SequenceGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedDBG(std::shared_ptr<SequenceGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 size_t num_threads = 0,
                 bool force_fast = false);

    std::vector<std::string> get_labels(node_index index) const;

    bool has_label(node_index index,
                   const std::string &label) const;

    void annotate_sequence(std::string&& sequence,
                           const std::vector<std::string> &labels);

    void call_annotated_nodes(const std::string &label,
                              std::function<void(node_index)> callback) const;

    void join() { thread_pool_.join(); }

    bool label_exists(const std::string &label) const;

    bool check_compatibility() const;

    const SequenceGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const SequenceGraph> get_graph_ptr() const { return graph_; }

    const Annotator& get_annotation() const { return *annotator_; }

    /*********************** Special queries **********************/

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string> get_labels(const std::string &sequence,
                                        double presence_ratio) const;
    // Weights don't need to sum to 1
    std::vector<std::string> get_labels(const std::vector<std::string> &sequences,
                                        const std::vector<double> &weights,
                                        double presence_ratio) const;
    std::vector<std::string> get_labels(const std::unordered_map<row_index, size_t> &index_counts,
                                        size_t min_count) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence,
                   size_t num_top_labels,
                   double presence_ratio = 0.0) const;

    // Weights don't need to sum to 1
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::vector<std::string> &sequences,
                   const std::vector<double> &weights,
                   size_t num_top_labels,
                   double presence_ratio = 0.0) const;

    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::unordered_map<row_index, size_t> &index_counts,
                   size_t num_top_labels,
                   size_t min_count = 0) const;

    static row_index graph_to_anno_index(node_index kmer_index);
    static node_index anno_to_graph_index(row_index anno_index);

  private:
    void annotate_sequence_thread_safe(const std::string &sequence,
                                       const std::vector<std::string> &labels);

    std::shared_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    ThreadPool thread_pool_;
    std::mutex mutex_;
    bool force_fast_;
};

#endif // __ANNOTATED_DBG_HPP__
