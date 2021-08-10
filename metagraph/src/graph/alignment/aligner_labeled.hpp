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
                                          bool force_fixed_seed) override final {
        last_buffered_table_i_ = 0;

        // the overridden backtrack populates extensions_, so this should return nothing
        DefaultColumnExtender::extend(min_path_score, force_fixed_seed);

        // fetch the alignments from extensions_
        return extensions_.get_alignments();
    }

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

    // since multi-node seeds may span across different labels, we no longer
    // want the restriction that the seed must be a prefix of the extended alignment
    virtual bool fixed_seed() const override final { return false; }

    // this override ensures that outgoing nodes are label- and coordinate-consistent
    // (when applicable)
    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char)> &callback,
                               size_t table_i) override final;

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
    size_t last_buffered_table_i_;

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
class LabeledAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::score_t score_t;

    LabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : ISeedAndExtendAligner<AlignmentCompare>(anno_graph.get_graph(), config),
            labeled_graph_(anno_graph) {
        if (labeled_graph_.get_coordinate_matrix()) {
            // do not use a global xdrop cutoff since we need separate cutoffs for each label
            this->config_.global_xdrop = false;
        }
    }

  protected:
    mutable AnnotationBuffer labeled_graph_;

    std::shared_ptr<IExtender>
    build_extender(std::string_view query,
                   const typename Extender::Aggregator &aggregator) const override final {
        return std::make_shared<Extender>(labeled_graph_, this->config_, aggregator, query);
    }

    std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<IDBGAligner::node_index> &nodes) const override final {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }

    // Generates seeds and extends them. If force_fixed_seed is true, then
    // all alignments must have the seed as a prefix. Otherwise, only the first
    // node of the seed is used as an alignment starting node.
    void align_core(const ISeeder &seeder,
                    IExtender &extender,
                    const std::function<void(Alignment&&)> &callback,
                    const std::function<score_t(const Alignment&)> &get_min_path_score,
                    bool force_fixed_seed) const override final {
        auto seeds = seeder.get_seeds();

        if (!force_fixed_seed && seeds.size())
            filter_seeds(seeds);

        for (const Alignment &seed : seeds) {
            if (seed.empty())
                continue;

            score_t min_path_score = get_min_path_score(seed);

            DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

            for (auto&& extension : extender.get_extensions(seed, min_path_score, force_fixed_seed)) {
                callback(std::move(extension));
            }
        }
    }

  private:
    // find the most frequent labels among the seeds and restrict graph traversal
    // to those labeled paths during extension
    void filter_seeds(std::vector<Alignment> &seeds) const {
        for (const Alignment &seed : seeds) {
            labeled_graph_.add_path(
                seed.get_nodes(),
                std::string(seed.get_nodes().size() + this->graph_.get_k() - 1, '#')
            );
        }

        labeled_graph_.flush();

        typedef AnnotationBuffer::Column Column;
        typedef AnnotationBuffer::node_index node_index;

        VectorMap<Column, uint64_t> label_counter;
        for (const Alignment &seed : seeds) {
            for (node_index node : seed.get_nodes()) {
                if (auto labels = labeled_graph_.get_labels(node)) {
                    for (uint64_t label : labels->get()) {
                        ++label_counter[label];
                    }
                }
            }
        }

        if (label_counter.empty())
            return;

        std::vector<std::pair<Column, uint64_t>> label_counts
            = const_cast<std::vector<std::pair<Column, uint64_t>>&&>(label_counter.values_container());

        std::sort(label_counts.begin(), label_counts.end(), utils::GreaterSecond());

        double cutoff = static_cast<double>(label_counts[0].second)
            * this->config_.rel_score_cutoff;
        auto it = std::find_if(label_counts.begin() + 1, label_counts.end(),
                               [cutoff](const auto &a) { return a.second < cutoff; });

        if (it == label_counts.end())
            return;

        label_counts.erase(it, label_counts.end());

        Vector<Column> labels;
        labels.reserve(label_counts.size());
        for (const auto &[label, count] : label_counts) {
            labels.push_back(label);
        }
        std::sort(labels.begin(), labels.end());

        for (Alignment &seed : seeds) {
            const std::vector<node_index> &nodes = seed.get_nodes();
            auto fetch_labels = labeled_graph_.get_labels(nodes[0]);
            if (!fetch_labels)
                continue;

            std::set_intersection(fetch_labels->get().begin(),
                                  fetch_labels->get().end(),
                                  labels.begin(), labels.end(),
                                  std::back_inserter(seed.label_columns));
            for (size_t i = 1; i < nodes.size() && seed.label_columns.size(); ++i) {
                if (auto next_fetch_labels = labeled_graph_.get_labels(nodes[i])) {
                    Vector<Column> temp;
                    std::set_intersection(next_fetch_labels->get().begin(),
                                          next_fetch_labels->get().end(),
                                          seed.label_columns.begin(),
                                          seed.label_columns.end(),
                                          std::back_inserter(temp));
                    std::swap(temp, seed.label_columns);
                } else {
                    seed.label_columns.clear();
                }
            }
        }

        auto seed_it = std::remove_if(seeds.begin(), seeds.end(), [](const auto &a) {
            return a.label_columns.empty();
        });

        seeds.erase(seed_it, seeds.end());
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
