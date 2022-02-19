#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include "dbg_aligner.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/hashers/hash.hpp"
#include "common/vector_set.hpp"
#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace graph {
namespace align {

// caches queried annotations to speed up next queries (labels with or w/o coordinates)
class AnnotationBuffer {
  public:
    typedef AnnotatedDBG::Annotator Annotator;
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;
    typedef Alignment::Columns Columns;
    typedef Alignment::CoordinateSet CoordinateSet;

    AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator);

    void queue_path(std::vector<node_index>&& path) {
        queued_paths_.push_back(std::move(path));
    }

    // fetch annotations for the queued nodes from the buffer and reset the buffer
    void fetch_queued_annotations();

    bool has_coordinates() const { return multi_int_; }

    // Get the annotations and coordinates of a node if they have been fetched.
    // The returned pointers are valid until next fetch_queued_annotations().
    std::pair<const Columns*, const CoordinateSet*>
    get_labels_and_coords(node_index node) const;

    // get the labels of a node if they have been fetched
    inline const Columns* get_labels(node_index node) const {
        return get_labels_and_coords(node).first;
    }

    const Annotator& get_annotator() const { return annotator_; }

    size_t num_nodes_buffered() const { return node_to_cols_.size(); }
    size_t num_column_sets() const { return column_sets_.size(); }

    // This method lets the caller push additional column sets to dictionary
    // `column_sets_`. These column sets can later be fetched with `get_cached_column_set()`
    template <typename... Args>
    inline size_t cache_column_set(Args&&... args) {
        auto it = column_sets_.emplace(std::forward<Args>(args)...).first;
        assert(std::is_sorted(it->begin(), it->end()));
        return it - column_sets_.begin();
    }

    // Fetch a label set given its index returned by `cache_column_set()`
    inline const Columns& get_cached_column_set(size_t i) const {
        assert(i < column_sets_.size());
        return column_sets_.data()[i];
    }

  private:
    const DeBruijnGraph &graph_;
    const Annotator &annotator_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // keep a unique set of annotation rows
    // the first element is the empty label set
    VectorSet<Columns, utils::VectorHash> column_sets_;
    // map nodes to indexes in column_sets_
    VectorMap<node_index, size_t> node_to_cols_; // map node to index in |column_sets_|
    // coordinate sets for all nodes in |node_to_cols_| in the same order
    std::vector<CoordinateSet> label_coords_;
    // buffer of paths to later querying with fetch_queued_annotations()
    std::vector<std::vector<node_index>> queued_paths_;
};


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
                                             node_index target_node = DeBruijnGraph::npos) override final {
        // extract all annotations for explored nodes
        flush();

        // run backtracking
        return DefaultColumnExtender::backtrack(min_path_score, window,
                                                right_end_bonus, target_node);
    }

    virtual bool set_seed(const Alignment &seed) override final;

    // overrides for backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &) const override final {
        // we are done with backtracking if all seed labels have been accounted for
        return !remaining_labels_i_;
    }

    virtual bool skip_backtrack_start(size_t i) override final;

    // since multi-node seeds may span across different labels, we no longer
    // want the restriction that the seed must be a prefix of the extended alignment
    virtual bool fixed_seed() const override final { return false; }

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
                                 const std::vector<score_t> &score_trace,
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
        assert(node_labels_.size() == node_labels_switched_.size());
        assert(i < node_labels_.size());
        DefaultColumnExtender::pop(i);
        last_flushed_table_i_ = std::min(i, last_flushed_table_i_);
        node_labels_.erase(node_labels_.begin() + i);
        node_labels_switched_.erase(node_labels_switched_.begin() + i);
    }

    // this override flushes the AnnotationBuffer, and checks elements in the
    // dynamic programming table for label- (and coordinate-)consistency
    void flush();

    // stores annotations for nodes
    AnnotationBuffer &annotation_buffer_;

    // index of the last dynamic programming table element whose node was flushed
    size_t last_flushed_table_i_;

    // map each table element to a corresponding label set index
    std::vector<size_t> node_labels_;
    std::vector<bool> node_labels_switched_;

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


template <class Seeder = SuffixSeeder<UniMEMSeeder>,
          class Extender = LabeledExtender,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public DBGAligner<Seeder, Extender, AlignmentCompare> {
    friend class LabeledExtender;
  public:
    typedef AnnotatedDBG::Annotator Annotator;

    LabeledAligner(const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config,
                   const Annotator &annotator);

    virtual ~LabeledAligner();

  private:
    mutable AnnotationBuffer annotation_buffer_;
    const Annotator &annotator_;

    typedef typename DBGAligner<Seeder, Extender, AlignmentCompare>::BatchSeeders BatchSeeders;
    BatchSeeders
    virtual build_seeders(const std::vector<IDBGAligner::Query> &seq_batch,
                          const std::vector<AlignmentResults> &wrapped_seqs) const override;

    // helper for the build_seeders method
    size_t filter_seeds(std::vector<Alignment> &seeds) const;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
