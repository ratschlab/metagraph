#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <optional>

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_set.hpp"
#include "common/hashers/hash.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

namespace mtg {

namespace annot {
namespace matrix {
    class MultiIntMatrix;
}
}

namespace graph {
namespace align {


class DynamicLabeledGraph {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::binmat::BinaryMatrix::Row Row;

    static constexpr Row nrow = std::numeric_limits<Row>::max();

    DynamicLabeledGraph(const AnnotatedDBG &anno_graph);

    const AnnotatedDBG& get_anno_graph() const { return anno_graph_; }
    const annot::matrix::MultiIntMatrix* get_coordinate_matrix() const { return multi_int_; }

    // flush the buffer and fetch their annotations from the AnnotatedDBG
    void flush();

    // push (a) node(s) to the buffer
    void add_node(node_index node);
    void add_path(const std::vector<node_index> &path, std::string sequence);

    // Checks if the edge (node, next), with its corresponding sequence, is
    // consistent w.r.t. coordinates. Always returns true if a coordinate annotation
    // is not loaded
    bool is_coord_consistent(node_index node, node_index next, std::string sequence) const;

    // get the annotations of a node if they have been fetched
    std::optional<std::reference_wrapper<const Vector<Column>>>
    operator[](node_index node) const {
        auto it = targets_.find(node);
        if (it == targets_.end() || it->second.second == nannot) {
            // if the node hasn't been seen before, or if its annotations haven't
            // been flushed, return nothing
            return std::nullopt;
        } else {
            return std::cref(targets_set_.data()[it->second.second]);
        }
    }

    std::vector<Row> get_anno_rows(const std::vector<node_index> &path) const {
        std::vector<Row> rows;
        rows.reserve(path.size());
        for (node_index n : path) {
            auto find = targets_.find(n);
            rows.push_back(find != targets_.end() ? find->second.first : nrow);
        }
        return rows;
    }

  private:
    const AnnotatedDBG &anno_graph_;
    const annot::matrix::MultiIntMatrix *multi_int_;

    // placeholder index for an unfetched annotation
    static constexpr size_t nannot = std::numeric_limits<size_t>::max();

    // keep a unique set of annotation rows
    VectorSet<Vector<Column>, utils::VectorHash> targets_set_;

    // map nodes to indexes in targets_set_
    tsl::hopscotch_map<node_index, std::pair<Row, size_t>> targets_;

    // buffer of accessed nodes and their corresponding annotation rows
    std::vector<Row> added_rows_;
    std::vector<node_index> added_nodes_;
};


template <class AlignmentCompare = LocalAlignmentLess>
class ILabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    ILabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(anno_graph.get_graph(), config),
            labeled_graph_(anno_graph) {}

    virtual ~ILabeledAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const IDBGAligner::AlignmentCallback &callback) const override {
        ISeedAndExtendAligner<AlignmentCompare>::align_batch(seq_batch,
            [&](std::string_view header, IDBGAligner::DBGQueryAlignment&& alignments) {
                auto it = std::remove_if(
                    alignments.begin(), alignments.end(),
                    [](const auto &a) { return a.target_columns.empty(); }
                );
                alignments.erase(it, alignments.end());

                for (auto &alignment : alignments) {
                    set_target_coordinates(alignment);
                }

                callback(header, std::move(alignments));
            }
        );
    }


  protected:
    mutable DynamicLabeledGraph labeled_graph_;

    void set_target_coordinates(IDBGAligner::DBGAlignment &alignment) const;
};


template <typename NodeType = DeBruijnGraph::node_index>
class LabeledBacktrackingExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef DynamicLabeledGraph::Column Column;
    typedef DefaultColumnExtender<DeBruijnGraph::node_index> BaseExtender;
    typedef typename BaseExtender::score_t score_t;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;
    typedef AlignmentAggregator<node_index, LocalAlignmentLess> Aggregator;

    LabeledBacktrackingExtender(DynamicLabeledGraph &labeled_graph,
                                const DBGAlignerConfig &config,
                                const Aggregator &aggregator,
                                std::string_view query)
          : BaseExtender(labeled_graph.get_anno_graph().get_graph(), config, query),
            labeled_graph_(labeled_graph),
            aggregator_(aggregator),
            no_chain_config_(disable_chaining(this->config_)),
            extensions_(labeled_graph.get_anno_graph().get_graph(),
                        aggregator_.get_query(false),
                        aggregator_.get_query(true), no_chain_config_) {}

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual std::vector<DBGAlignment> extend(score_t min_path_score, bool fixed_seed) override;

    virtual void init_backtrack() override {
        labeled_graph_.flush();
        diff_target_sets_.clear();
    }

    virtual bool terminate_backtrack_start(const std::vector<DBGAlignment> &) const override final { return false; }

    virtual bool terminate_backtrack() const override final { return target_intersection_.empty(); }

    virtual bool skip_backtrack_start(size_t i) override final;

    virtual bool update_seed_filter(node_index node,
                                    size_t query_start,
                                    const score_t *s_begin,
                                    const score_t *s_end) override final;

    virtual bool fixed_seed() const override final { return false; }

    virtual void call_outgoing(NodeType node,
                               size_t max_prefetch_distance,
                               const std::function<void(NodeType, char)> &callback,
                               size_t table_idx) override final;

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
                                 const std::function<void(DBGAlignment&&)> &callback) override final;

  private:
    DynamicLabeledGraph &labeled_graph_;

    // global set of alignments
    const Aggregator &aggregator_;

    // local set of alignments
    DBGAlignerConfig no_chain_config_;
    Aggregator extensions_;

    // keep track of the label set for the current backtracking
    Vector<Column> target_intersection_;
    size_t last_path_size_;

    // After a node has been visited during backtracking, we keep track of which
    // of its labels haven't been considered yet. This way, if backtracking is
    // called from this node, then we can restrict it to these labels.
    tsl::hopscotch_map<size_t, Vector<Column>> diff_target_sets_;

    static DBGAlignerConfig disable_chaining(DBGAlignerConfig config) {
        config.chain_alignments = false;
        return config;
    }
};


template <class Extender = LabeledBacktrackingExtender<>,
          class Seeder = UniMEMSeeder<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ILabeledAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    LabeledAligner(Args&&... args)
          : ILabeledAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  protected:
    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query,
                   const typename Extender::Aggregator &aggregator) const override final {
        return std::make_shared<Extender>(this->labeled_graph_, this->get_config(), aggregator, query);
    }

    std::shared_ptr<ISeeder<DeBruijnGraph::node_index>>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<DeBruijnGraph::node_index> &nodes) const override final {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
