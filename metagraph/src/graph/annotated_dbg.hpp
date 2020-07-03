#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <cassert>
#include <memory>
#include <mutex>

#include <sdsl/int_vector.hpp>

#include "representation/base/sequence_graph.hpp"
#include "annotation/representation/base/annotation.hpp"

using namespace mtg;


class AnnotatedSequenceGraph {
  public:
    typedef anno::MultiLabelEncoded<std::string> Annotator;
    using node_index = SequenceGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graphh,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast = false);

    virtual ~AnnotatedSequenceGraph() {}

    virtual std::vector<std::string> get_labels(node_index index) const;

    virtual bool has_label(node_index index, const std::string &label) const;

    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequence(const std::string &sequence,
                                   const std::vector<std::string> &labels);

    virtual void call_annotated_nodes(const std::string &label,
                                      std::function<void(node_index)> callback) const;

    virtual bool label_exists(const std::string &label) const;

    virtual bool check_compatibility() const;

    virtual const SequenceGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const SequenceGraph> get_graph_ptr() const { return graph_; }

    virtual const Annotator& get_annotation() const { return *annotator_; }

    static row_index graph_to_anno_index(node_index kmer_index) {
        assert(kmer_index);
        return kmer_index - 1;
    }
    static node_index anno_to_graph_index(row_index anno_index) {
        return anno_index + 1;
    }

  protected:
    std::shared_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    std::mutex mutex_;
    bool force_fast_;
};


class AnnotatedDBG : public AnnotatedSequenceGraph {
  public:
    typedef anno::MultiLabelEncoded<std::string> Annotator;
    using node_index = DeBruijnGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 bool force_fast = false);

    using AnnotatedSequenceGraph::get_labels;

    const DeBruijnGraph& get_graph() const { return dbg_; }

    /*********************** Special queries **********************/

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<std::string> get_labels(const std::string &sequence,
                                        double presence_ratio) const;

    std::vector<std::string> get_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                                        size_t min_count) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::string &sequence,
                   size_t num_top_labels,
                   double presence_ratio = 0.0) const;

    std::vector<std::pair<std::string, size_t>>
    get_top_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                   size_t num_top_labels,
                   size_t min_count = 0) const;

    std::vector<std::pair<std::string, sdsl::bit_vector>>
    get_top_label_signatures(const std::string &sequence,
                             size_t num_top_labels,
                             double presence_ratio = 0.0) const;

    int32_t score_kmer_presence_mask(const sdsl::bit_vector &kmer_presence_mask,
                                     int32_t match_score = 1,
                                     int32_t mismatch_score = 2) const;

  private:
    DeBruijnGraph &dbg_;
};

#endif // __ANNOTATED_DBG_HPP__
