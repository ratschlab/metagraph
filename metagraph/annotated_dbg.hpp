#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include "dbg_succinct.hpp"
#include "annotate.hpp"
#include "utils.hpp"


class AnnotatedDBG {
  public:
    typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

    AnnotatedDBG(DBG_succ *dbg, Annotator *annotation = NULL, size_t num_threads = 1);
    AnnotatedDBG(Annotator *annotation, DBG_succ *dbg = NULL, size_t num_threads = 1);

    void adjust_annotation(const bit_vector_dyn &inserted_edges);

    void annotate_sequence(const std::string &sequence,
                           const std::vector<std::string> &labels);

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string>
    get_labels(const std::string &sequence, double presence_ratio) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence, size_t num_top_labels) const;

    Annotator& get_annotation() { return *annotator_; }
    const Annotator& get_annotation() const { return *annotator_; }

    DBG_succ& get_graph() { return *graph_; }
    const DBG_succ& get_graph() const { return *graph_; }

    void join() { thread_pool_.join(); }

  private:
    std::unique_ptr<DBG_succ> graph_;
    std::unique_ptr<Annotator> annotator_;

    utils::ThreadPool thread_pool_;
    std::mutex mutex_;
};

#endif // __ANNOTATED_DBG_HPP__