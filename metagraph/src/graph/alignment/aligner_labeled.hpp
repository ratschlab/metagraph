#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include "dbg_aligner.hpp"
#include "annotation_buffer.hpp"
#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace graph {
namespace align {


class LabeledExtender : public DefaultColumnExtender {
  public:
    typedef AnnotationBuffer::Columns Columns;
    typedef AnnotationBuffer::CoordinateSet CoordinateSet;

    LabeledExtender(const DeBruijnGraph &graph,
                    const DBGAlignerConfig &config,
                    AnnotationBuffer &annotation_buffer,
                    std::string_view query)
          : DefaultColumnExtender(graph, config, query),
            annotation_buffer_(annotation_buffer) {}

    // |aligner| must be an instance of LabeledAligner<>
    LabeledExtender(const IDBGAligner &aligner, std::string_view query);

    virtual ~LabeledExtender() {}

  private:
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window,
                                             score_t right_end_bonus,
                                             const std::vector<size_t> &tips,
                                             node_index target_node = DeBruijnGraph::npos) override final {
        // extract all annotations for explored nodes
        flush();

        // run backtracking
        auto alignments = DefaultColumnExtender::backtrack(
            min_path_score, window, right_end_bonus, tips, target_node
        );

        for (Alignment &alignment : alignments) {
            alignment.label_encoder = &annotation_buffer_.get_annotator().get_label_encoder();
        }

        return alignments;
    }

    virtual bool set_seed(const Alignment &seed) override final;

    // overrides for backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &) const override final {
        // we are done with backtracking if all seed labels have been accounted for
        return !remaining_labels_i_;
    }

    virtual bool skip_backtrack_start(size_t i) override final;

    // this override ensures that outgoing nodes are label- and coordinate-consistent
    // (when applicable)
    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char, score_t)> &callback,
                               size_t table_i,
                               bool force_fixed_seed = false) override final;

    // this method calls multiple label-consistent alignments by backtracking
    virtual void call_alignments(score_t end_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> &trace,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 score_t extra_score,
                                 const std::function<void(Alignment&&)> &callback) override final;

    // ensure that when terminating a path early, the per-node label storage is also
    // correctly handled
    virtual void pop(size_t i) override final {
        assert(i < node_labels_.size());
        DefaultColumnExtender::pop(i);
        node_labels_.erase(node_labels_.begin() + i);
    }

    // this override flushes the AnnotationBuffer, and checks elements in the
    // dynamic programming table for label- (and coordinate-)consistency
    void flush(size_t table_i = 0, const std::vector<node_index> &outgoing = {});

    // stores annotations for nodes
    AnnotationBuffer &annotation_buffer_;

    // map each table element to a corresponding label set index
    std::vector<size_t> node_labels_;

    // if the seed has coordinates, then the coordinates of the initial node in the
    // extension is stored here
    CoordinateSet base_coords_;

    // index of the set of labels that still have not been seen during backtracking
    size_t remaining_labels_i_;

    // the label set intersection and difference between the current backtracking
    // operation and the ones still left to be observed
    Columns label_intersection_;
    Columns label_diff_;
};

class ILabeledAligner {
  public:
    virtual AnnotationBuffer& get_annotation_buffer() const = 0;
};


template <class Seeder = SuffixSeeder<ExactSeeder>,
          class Extender = LabeledExtender,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public DBGAligner<Seeder, Extender, AlignmentCompare>, public ILabeledAligner {
    friend class LabeledExtender;
  public:
    typedef AnnotatedDBG::Annotator Annotator;

    LabeledAligner(const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config,
                   const Annotator &annotator);

    virtual ~LabeledAligner();

    virtual AnnotationBuffer& get_annotation_buffer() const override final {
        return annotation_buffer_;
    }

    virtual bool has_coordinates() const override final {
        return annotation_buffer_.has_coordinates();
    }

  private:
    mutable AnnotationBuffer annotation_buffer_;

    typedef typename DBGAligner<Seeder, Extender, AlignmentCompare>::BatchSeeders BatchSeeders;
    BatchSeeders
    virtual build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                          const std::vector<AlignmentResults> &wrapped_seqs) const override final;

    // helper for the build_seeders method
    size_t filter_seeds(std::vector<Seed> &seeds) const;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
