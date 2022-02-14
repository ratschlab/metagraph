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

class AnnotationBuffer {
  public:
    typedef AnnotatedDBG::Annotator Annotator;
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;
    typedef Alignment::LabelSet LabelSet;
    typedef Alignment::CoordinateSet CoordinateSet;

    // placeholder index for an unfetched annotation
    static constexpr size_t nannot = std::numeric_limits<size_t>::max();

    AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator);

    const DeBruijnGraph& get_graph() const { return graph_; }
    const Annotator& get_annotator() const { return annotator_; }
    const annot::matrix::MultiIntMatrix* get_coordinate_matrix() const { return multi_int_; }

    // flush the buffer and fetch their annotations from the Annotator
    void flush();

    // push a node to the buffer
    node_index add_node(node_index node);

    // Push the nodes in a path to the buffer. Return the a pair with the
    // vector of corresponding labeled nodes, and a bool indicating if this vector
    // is reversed relative to the input path.
    std::pair<std::vector<node_index>, bool>
    add_path(const std::vector<node_index> &path, std::string&& sequence);

    // get the annotations and coordinates of a node if they have been fetched
    std::pair<const LabelSet*, const CoordinateSet*>
    get_labels_and_coordinates(node_index node) const {
        std::pair<const LabelSet*, const CoordinateSet*> ret_val { nullptr, nullptr };

        auto it = labels_.find(node);

        // if the node hasn't been seen before, or if its annotations haven't
        // been flushed, return nothing
        if (it == labels_.end() || it->second.second == nannot)
            return ret_val;

        ret_val.first = &labels_set_.data()[it->second.second];
        assert(std::is_sorted(ret_val.first->begin(), ret_val.first->end()));

        if (multi_int_) {
            assert(static_cast<size_t>(it - labels_.begin()) < label_coords_.size());
            ret_val.second = &label_coords_[it - labels_.begin()];
        }

        return ret_val;
    }

    // get the labels of a node if they have been fetched
    inline const LabelSet* get_labels(node_index node) const {
        return get_labels_and_coordinates(node).first;
    }

    // get the number of nodes that are waiting to have their annotations fetched
    size_t num_cached() const { return added_rows_.size(); }

    // add a label set (represented as a vector of sorted column indices) to the buffer
    template <typename... Args>
    size_t emplace_label_set(Args&&... args) {
        auto it = labels_set_.emplace(std::forward<Args>(args)...).first;
        assert(std::is_sorted(it->begin(), it->end()));
        return it - labels_set_.begin();
    }

    // get the index of the label set
    size_t get_index(const LabelSet &labels) const {
        assert(std::is_sorted(labels.begin(), labels.end()));
        auto find = labels_set_.find(labels);
        return find != labels_set_.end() ? find - labels_set_.begin() : nannot;
    }

    // fetch a label set given its index
    const LabelSet& get_labels_from_index(size_t i) const {
        assert(i != nannot);
        assert(i < labels_set_.size());
        assert(std::is_sorted(labels_set_.data()[i].begin(),
                              labels_set_.data()[i].end()));
        return labels_set_.data()[i];
    }

    size_t num_nodes_buffered() const { return labels_.size(); }
    size_t num_label_sets() const { return labels_set_.size(); }

  private:
    const DeBruijnGraph &graph_;
    const Annotator &annotator_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // keep a unique set of annotation rows
    // the first element is the empty label set
    VectorSet<LabelSet, utils::VectorHash> labels_set_;

    // map nodes to indexes in labels_set_
    VectorMap<node_index, std::pair<Row, size_t>> labels_;

    // map each key in labels_ to a set of coordinates
    std::vector<CoordinateSet> label_coords_;

    // buffer of accessed nodes and their corresponding annotation rows
    std::vector<node_index> added_nodes_;
    std::vector<Row> added_rows_;
};

class LabeledExtender : public DefaultColumnExtender {
  public:
    typedef AnnotationBuffer::Column Column;
    typedef AnnotationBuffer::Tuple Tuple;
    typedef AnnotationBuffer::LabelSet LabelSet;
    typedef AnnotationBuffer::CoordinateSet CoordinateSet;

    LabeledExtender(AnnotationBuffer &annotation_buffer,
                    const DBGAlignerConfig &config,
                    std::string_view query)
          : DefaultColumnExtender(annotation_buffer.get_graph(), config, query),
            annotation_buffer_(annotation_buffer) {}

    // |aligner| must be an instance of LabeledAligner<>
    static LabeledExtender make(const IDBGAligner &aligner,
                                const DBGAlignerConfig &config,
                                std::string_view query);

    virtual ~LabeledExtender() {}

  private:
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window,
                                             node_index target_node = DeBruijnGraph::npos) override final {
        // extract all annotations for explored nodes
        flush();

        // run backtracking
        return DefaultColumnExtender::backtrack(min_path_score, window, target_node);
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
    virtual void call_alignments(score_t cur_cell_score,
                                 score_t end_score,
                                 score_t min_path_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> &trace,
                                 size_t table_i,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 score_t extra_penalty,
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
    void flush();

    // stores annotations for nodes
    AnnotationBuffer &annotation_buffer_;

    // index of the last dynamic programming table element whose node was flushed
    size_t last_flushed_table_i_;

    // map each table element to a corresponding label set index
    std::vector<size_t> node_labels_;

    // if the seed has coordinates, then the coordinates of the initial node in the
    // extension is stored here
    CoordinateSet base_coords_;

    // index of the set of labels that still have not been seen during backtracking
    size_t remaining_labels_i_;

    // the label set intersection and difference between the current backtracking
    // operation and the ones still left to be observed
    LabelSet label_intersection_;
    LabelSet label_diff_;
};

template <class Seeder = SuffixSeeder<UniMEMSeeder>,
          class Extender = LabeledExtender,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public DBGAligner<Seeder, Extender, AlignmentCompare> {
  public:
    typedef AnnotationBuffer::Annotator Annotator;
    typedef AnnotationBuffer::LabelSet LabelSet;
    typedef Alignment::node_index node_index;
    typedef Alignment::Column Column;

    LabeledAligner(const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config,
                   const Annotator &annotator);

    virtual ~LabeledAligner();

    auto& get_annotation_buffer() const { return annotation_buffer_; }

  private:
    typedef typename DBGAligner<Seeder, Extender, AlignmentCompare>::BatchSeeders BatchSeeders;
    mutable AnnotationBuffer annotation_buffer_;

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
