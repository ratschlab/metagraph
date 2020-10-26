#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <cassert>
#include <memory>
#include <mutex>

#include <sdsl/int_vector.hpp>

#include "representation/base/sequence_graph.hpp"
#include "annotation/representation/base/annotation.hpp"


namespace mtg {
namespace graph {

class AnnotatedSequenceGraph {
  public:
    typedef std::string Label;
    typedef annot::MultiLabelEncoded<Label> Annotator;
    using node_index = SequenceGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graphh,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast = false);

    virtual ~AnnotatedSequenceGraph() {}

    virtual std::vector<Label> get_labels(node_index index) const;

    virtual bool has_label(node_index index, const Label &label) const;

    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequence(std::string_view sequence,
                                   const std::vector<Label> &labels);

    virtual void call_annotated_nodes(const Label &label,
                                      std::function<void(node_index)> callback) const;

    virtual bool label_exists(const Label &label) const;

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
    AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 bool force_fast = false);

    using AnnotatedSequenceGraph::get_labels;

    const DeBruijnGraph& get_graph() const { return dbg_; }

    /*********************** Special queries **********************/

    // return labels that occur at least in |presence_ratio| k-mers
    std::vector<Label> get_labels(std::string_view sequence,
                                  double presence_ratio) const;

    std::vector<Label> get_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                                  size_t min_count) const;

    // return top |num_top_labels| labels with their counts
    std::vector<std::pair<Label, size_t>>
    get_top_labels(std::string_view sequence,
                   size_t num_top_labels,
                   double presence_ratio = 0.0) const;

    std::vector<std::pair<Label, size_t>>
    get_top_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                   size_t num_top_labels,
                   size_t min_count = 0) const;

    std::vector<std::pair<Label, sdsl::bit_vector>>
    get_top_label_signatures(std::string_view sequence,
                             size_t num_top_labels,
                             double presence_ratio = 0.0) const;

    int32_t score_kmer_presence_mask(const sdsl::bit_vector &kmer_presence_mask,
                                     int32_t match_score = 1,
                                     int32_t mismatch_score = 2) const;

  private:
    DeBruijnGraph &dbg_;
};

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_DBG_HPP__
