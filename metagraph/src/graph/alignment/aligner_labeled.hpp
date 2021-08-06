#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <optional>

#include <tsl/hopscotch_map.h>
#include <tsl/ordered_set.h>

#include "dbg_aligner.hpp"
#include "common/hashers/hash.hpp"
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
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;

    typedef std::reference_wrapper<const Alignment::LabelSet> LabelSet;
    typedef std::reference_wrapper<const Alignment::CoordinateSet> CoordsSet;

    static constexpr Row nrow = std::numeric_limits<Row>::max();

    AnnotationBuffer(const AnnotatedDBG &anno_graph);

    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }
    const annot::matrix::MultiIntMatrix* get_coordinate_matrix() const { return multi_int_; }

    // flush the buffer and fetch their annotations from the AnnotatedDBG
    void flush();

    // push (a) node(s) to the buffer
    void add_node(node_index node);
    void add_path(const std::vector<node_index> &path, std::string sequence);

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

    // get the annotations of a node if they have been fetched
    inline std::optional<LabelSet> get_labels(node_index node) const {
        return get_labels_and_coordinates(node).first;
    }

  private:
    const AnnotatedDBG &anno_graph_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // placeholder index for an unfetched annotation
    static constexpr size_t nannot = std::numeric_limits<size_t>::max();

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


template <class AlignmentCompare = LocalAlignmentLess>
class ILabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    typedef IDBGAligner::node_index node_index;

    ILabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(anno_graph.get_graph(), config),
            labeled_graph_(anno_graph) {}

    virtual ~ILabeledAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const IDBGAligner::AlignmentCallback &callback) const override {
        ISeedAndExtendAligner<AlignmentCompare>::align_batch(seq_batch,
            [&](std::string_view header, QueryAlignment&& alignments) {
                alignments.data().erase(
                    std::remove_if(alignments.data().begin(), alignments.data().end(),
                                   [](const auto &a) { return a.label_columns.empty(); }),
                    alignments.data().end()
                );
                callback(header, std::move(alignments));
            }
        );
    }


  protected:
    mutable AnnotationBuffer labeled_graph_;
};


class LabeledBacktrackingExtender : public DefaultColumnExtender {
  public:
    typedef AnnotationBuffer::Column Column;
    typedef AnnotationBuffer::Tuple Tuple;
    typedef AlignmentAggregator<LocalAlignmentLess> Aggregator;

    LabeledBacktrackingExtender(AnnotationBuffer &labeled_graph,
                                const DBGAlignerConfig &config,
                                const Aggregator &aggregator,
                                std::string_view query)
          : DefaultColumnExtender(labeled_graph.get_anno_graph().get_graph(), config, query),
            labeled_graph_(labeled_graph),
            aggregator_(aggregator),
            no_chain_config_(disable_chaining(this->config_)),
            extensions_(labeled_graph.get_anno_graph().get_graph(),
                        aggregator_.get_query(false),
                        aggregator_.get_query(true), no_chain_config_) {}

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed) override final;

    // backtrack through the DP table to reconstruct alignments
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window) override final {
        // extract all labels for explored nodes
        labeled_graph_.flush();

        // reset the per-node temporary label storage
        diff_label_sets_.clear();

        // run backtracking
        return DefaultColumnExtender::backtrack(min_path_score, window);
    }

    // overrides for backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &) const override final { return false; }
    virtual bool terminate_backtrack() const override final { return label_intersection_.empty(); }
    virtual bool skip_backtrack_start(size_t i) override final;

    virtual bool update_seed_filter(node_index node,
                                    size_t query_start,
                                    const score_t *s_begin,
                                    const score_t *s_end) override final {
        if (SeedFilteringExtender::update_seed_filter(node, query_start, s_begin, s_end)) {
            // if this node was explored, add it to the buffer
            labeled_graph_.add_node(node);
            return true;
        } else {
            return false;
        }
    }

    // since multi-node seeds may span across different labels, we no longer
    // want the restriction that the seed must be a prefix of the extended alignment
    virtual bool fixed_seed() const override final { return false; }

    // this override ensures that outgoing nodes are label- and coordinate-consistent
    // (when applicable)
    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char)> &callback) override final;

    // this method calls multiple label-consistent alignments by backtracking
    virtual void call_alignments(score_t cur_cell_score,
                                 score_t end_score,
                                 score_t min_path_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> &trace,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 const std::function<void(Alignment&&)> &callback) override final;

  private:
    AnnotationBuffer &labeled_graph_;

    // global set of alignments
    const Aggregator &aggregator_;

    // local set of alignments
    DBGAlignerConfig no_chain_config_;
    Aggregator extensions_;

    // keep track of the label set for the current backtracking
    Vector<Column> label_intersection_;
    size_t last_path_size_;

    Vector<Tuple> label_intersection_coords_;

    // After a node has been visited during backtracking, we keep track of which
    // of its labels haven't been considered yet. This way, if backtracking is
    // called from this node, then we can restrict it to these labels.
    tsl::hopscotch_map<size_t, Vector<Column>> diff_label_sets_;

    // we don't want to chain alignments in the local set of alignments, so
    // this generates a modified config for the local aggregator
    static DBGAlignerConfig disable_chaining(DBGAlignerConfig config) {
        config.chain_alignments = false;
        return config;
    }
};


template <class Extender = LabeledBacktrackingExtender,
          class Seeder = UniMEMSeeder,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ILabeledAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    LabeledAligner(Args&&... args)
          : ILabeledAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  protected:
    std::shared_ptr<IExtender>
    build_extender(std::string_view query,
                   const typename Extender::Aggregator &aggregator) const override final {
        return std::make_shared<Extender>(this->labeled_graph_, this->get_config(), aggregator, query);
    }

    std::shared_ptr<ISeeder>
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
