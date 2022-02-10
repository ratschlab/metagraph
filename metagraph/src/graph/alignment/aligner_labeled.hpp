#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <optional>

#include <tsl/hopscotch_map.h>
#include <tsl/ordered_set.h>

#include "dbg_aligner.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace graph {
namespace align {

template <typename Key, class Hash = std::hash<Key>, typename IndexType = uint64_t,
          class EqualTo = std::equal_to<Key>, class Allocator = std::allocator<Key>,
          class Container = std::vector<Key, Allocator>>
using VectorSet = tsl::ordered_set<Key, Hash, EqualTo, Allocator, Container, IndexType>;


class AnnotationBuffer {
  public:
    typedef AnnotatedDBG::Annotator Annotator;
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;

    typedef std::reference_wrapper<const Alignment::LabelSet> LabelSet;
    typedef std::reference_wrapper<const Alignment::CoordinateSet> CoordsSet;

    static constexpr Row nrow = std::numeric_limits<Row>::max();

    // placeholder index for an unfetched annotation
    static constexpr size_t nannot = std::numeric_limits<size_t>::max();

    AnnotationBuffer(const DeBruijnGraph &graph, const Annotator &annotator);

    const DeBruijnGraph& get_graph() const { return graph_; }
    const Annotator& get_annotator() const { return annotator_; }
    const annot::matrix::MultiIntMatrix* get_coordinate_matrix() const { return multi_int_; }

    // flush the buffer and fetch their annotations from the AnnotatedDBG
    void flush();

    // push (a) node(s) to the buffer
    node_index add_node(node_index node);

    std::pair<std::vector<node_index>, bool>
    add_path(const std::vector<node_index> &path, std::string&& sequence);

    // get the annotations and coordinates of a node if they have been fetched
    std::pair<std::optional<LabelSet>, std::optional<CoordsSet>>
    get_labels_and_coordinates(node_index node) const {
        std::pair<std::optional<LabelSet>, std::optional<CoordsSet>> ret_val {
            std::nullopt, std::nullopt
        };

        auto it = labels_.find(node);

        // if the node hasn't been seen before, or if its annotations haven't
        // been flushed, return nothing
        if (it == labels_.end() || it->second.second == nannot)
            return ret_val;

        ret_val.first = std::cref(labels_set_.data()[it->second.second]);

        // if no coordinates are present, return just the labels
        if (!multi_int_)
            return ret_val;

        assert(static_cast<size_t>(it - labels_.begin()) < label_coords_.size());
        ret_val.second = std::cref(label_coords_[it - labels_.begin()]);
        return ret_val;
    }

    inline bool is_flushed(node_index node) const {
        auto it = labels_.find(node);

        // if the node hasn't been seen before, or if its annotations haven't
        // been flushed, return false
        return it != labels_.end() && it->second.second != nannot;
    }

    inline bool is_flushed(const std::vector<node_index> &nodes) const {
        for (node_index node : nodes) {
            if (!is_flushed(node))
                return false;
        }

        return true;
    }

    // get the annotations of a node if they have been fetched
    inline std::optional<LabelSet> get_labels(node_index node) const {
        return get_labels_and_coordinates(node).first;
    }

    size_t num_cached() const { return added_rows_.size(); }

    template <typename... Args>
    size_t emplace_label_set(Args&&... args) {
        auto it = labels_set_.emplace(std::forward<Args>(args)...).first;
        return it - labels_set_.begin();
    }

    size_t get_index(const Vector<Column> &labels) const {
        auto find = labels_set_.find(labels);
        if (find == labels_set_.end())
            return nannot;

        return find - labels_set_.begin();
    }

    const Vector<Column>& get_labels_from_index(size_t i) const {
        assert(i != nannot);
        assert(i < labels_set_.size());
        return labels_set_.data()[i];
    }

  private:
    const DeBruijnGraph &graph_;
    const Annotator &annotator_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // keep a unique set of annotation rows
    VectorSet<Vector<Column>, utils::VectorHash> labels_set_;

    // map nodes to indexes in labels_set_
    VectorMap<node_index, std::pair<Row, size_t>> labels_;

    // map each element in labels_ to a set of coordinates
    std::vector<Vector<Tuple>> label_coords_;

    // buffer of accessed nodes and their corresponding annotation rows
    std::vector<Row> added_rows_;
    std::vector<node_index> added_nodes_;
};

class LabeledExtender : public DefaultColumnExtender {
  public:
    typedef AnnotationBuffer::Column Column;
    typedef AnnotationBuffer::Tuple Tuple;

    LabeledExtender(AnnotationBuffer &labeled_graph,
                    const DBGAlignerConfig &config,
                    std::string_view query)
          : DefaultColumnExtender(labeled_graph.get_graph(), config, query),
            labeled_graph_(labeled_graph) {}

    virtual ~LabeledExtender() {}

  private:
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window) override final {
        // extract all labels for explored nodes
        flush();

        // run backtracking
        return DefaultColumnExtender::backtrack(min_path_score, window);
    }

    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed) override final {
        // the first node of the seed has already been buffered and flushed
        last_buffered_table_i_ = 1;
        last_flushed_table_i_ = 1;
        return DefaultColumnExtender::extend(min_path_score, force_fixed_seed);
    }

    virtual bool set_seed(const Alignment &seed) override final {
        if (DefaultColumnExtender::set_seed(seed)) {
            assert(labeled_graph_.is_flushed(seed.get_nodes()));
            fetched_label_i_ = labeled_graph_.emplace_label_set(seed.label_columns);
            assert(fetched_label_i_ != AnnotationBuffer::nannot);
            node_labels_.assign(1, fetched_label_i_);
            return true;
        }

        return false;
    }

    // overrides for backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &) const override final {
        return !fetched_label_i_;
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

    virtual void pop(size_t i) override final {
        assert(i < node_labels_.size());
        DefaultColumnExtender::pop(i);
        node_labels_.erase(node_labels_.begin() + i, node_labels_.begin() + i + 1);
    }

    void flush();

    AnnotationBuffer &labeled_graph_;
    size_t last_buffered_table_i_;
    size_t last_flushed_table_i_;
    std::vector<size_t> node_labels_;
    size_t fetched_label_i_;
    Vector<Column> label_intersection_;
};

template <class AlignmentCompare = LocalAlignmentLess>
class ILabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    typedef AnnotationBuffer::Annotator Annotator;
    typedef Alignment::score_t score_t;
    typedef Alignment::node_index node_index;
    typedef Alignment::Column Column;

  protected:
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::BatchSeeders BatchSeeders;
    mutable AnnotationBuffer labeled_graph_;

    ILabeledAligner(const DeBruijnGraph &graph,
                    const Annotator &annotator,
                    const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(graph, config),
            labeled_graph_(graph, annotator) {}

    virtual void filter_seeds(BatchSeeders &seeders) const override final;

  private:
    // find the most frequent labels among the seeds and restrict graph traversal
    // to those labeled paths during extension
    size_t filter_seeds(std::vector<Alignment> &seeds) const;
};

template <class AlignmentCompare = LocalAlignmentLess,
          class Extender = LabeledExtender,
          class Seeder = UniMEMSeeder>
class LabeledAligner : public ILabeledAligner<AlignmentCompare> {
  public:
    typedef typename ILabeledAligner<AlignmentCompare>::Annotator Annotator;

    LabeledAligner(const DeBruijnGraph &graph,
                   const Annotator &annotator,
                   const DBGAlignerConfig &config)
          : ILabeledAligner<AlignmentCompare>(graph, annotator, config) {
        if (this->labeled_graph_.get_coordinate_matrix()
                && std::is_same_v<Extender, LabeledExtender>) {
            // do not use a global xdrop cutoff since we need separate cutoffs
            // for each label
            this->config_.global_xdrop = false;
        }
    }

  protected:
    virtual std::shared_ptr<IExtender>
    build_extender(std::string_view query, const DBGAlignerConfig &config) const override final {
        return std::make_shared<Extender>(this->labeled_graph_, config, query);
    }

    virtual std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<IDBGAligner::node_index> &nodes) const override final {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
