#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_map.hpp"

namespace mtg {
namespace graph {

class AnnotatedDBG;

namespace align {


class ILabeledDBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::QueryGenerator QueryGenerator;

    // undefined target column
    static constexpr uint64_t kNTarget = std::numeric_limits<uint64_t>::max() - 1;

    ILabeledDBGAligner(const AnnotatedDBG &anno_graph,
                       const DBGAlignerConfig &config,
                       size_t num_top_labels = 1);

    virtual ~ILabeledDBGAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    typedef std::pair<std::vector<node_index> /* forward */,
                      std::vector<node_index> /* reverse complement */ > Mapping;
    typedef std::pair<sdsl::bit_vector /* forward */,
                      sdsl::bit_vector /* reverse complement */ > Signature;
    typedef std::vector<Mapping> BatchMapping;
    typedef std::vector<VectorMap<uint64_t, Signature>> BatchLabels;

    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
    size_t num_top_labels_;

    std::pair<BatchMapping, BatchLabels>
    map_and_label_query_batch(const QueryGenerator &generate_query) const;
};

template <class BaseSeeder>
class LabeledSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::node_index node_index;
    typedef typename BaseSeeder::Seed Seed;

    template <typename... Args>
    LabeledSeeder(uint64_t target_column, Args&&... args)
          : BaseSeeder(std::forward<Args>(args)...), target_column_(target_column) {}

    virtual ~LabeledSeeder() {}

    virtual std::vector<Seed> get_seeds() const override {
        std::vector<Seed> seeds = BaseSeeder::get_seeds();
        for (Seed &seed : seeds) {
            seed.target_column = !seed.get_offset()
                ? target_column_
                : ILabeledDBGAligner::kNTarget;
        }

        return seeds;
    }

  private:
    uint64_t target_column_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class BaseSeeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledDBGAligner : public ILabeledDBGAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    template <typename... Args>
    LabeledDBGAligner(Args&&... args) : ILabeledDBGAligner(std::forward<Args>(args)...) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

  protected:
    typedef LabeledSeeder<BaseSeeder> Seeder;

    Seeder build_seeder(uint64_t target_column,
                        std::string_view query,
                        bool is_reverse_complement,
                        std::vector<node_index>&& base_nodes,
                        const sdsl::bit_vector &signature) const;
};

template <typename NodeType>
class LabeledColumnExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::score_t score_t;
    typedef typename IExtender<NodeType>::node_index node_index;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~LabeledColumnExtender() {}

    virtual void initialize(const DBGAlignment &path) override;

  protected:
    typedef std::vector<std::pair<NodeType, char>> Edges;
    typedef typename DefaultColumnExtender<NodeType>::Column Column;
    typedef typename DefaultColumnExtender<NodeType>::Scores Scores;
    typedef typename DefaultColumnExtender<NodeType>::AlignNode AlignNode;
    typedef typename DefaultColumnExtender<NodeType>::AlignNodeHash AlignNodeHash;

    virtual Edges get_outgoing(const AlignNode &node) const override;

    virtual void add_scores_to_column(Column &column,
                                      Scores&& scores,
                                      const AlignNode &node) override {
        const AlignNode &prev = std::get<6>(scores);
        assert(align_node_to_target_.count(prev));

        if (!align_node_to_target_.count(node)) {
            align_node_to_target_[node] = align_node_to_target_[prev];
#ifndef NDEBUG
        } else {
            assert(align_node_to_target_[prev] == align_node_to_target_[node]
                || target_columns_.at(align_node_to_target_[prev])
                    == ILabeledDBGAligner::kNTarget);
#endif
        }

        DefaultColumnExtender<NodeType>::add_scores_to_column(
            column, std::move(scores), node
        );
    }

  private:
    const AnnotatedDBG &anno_graph_;
    mutable std::vector<uint64_t> target_columns_;
    mutable tsl::hopscotch_map<NodeType, std::vector<Edges>> cached_edge_sets_;
    mutable tsl::hopscotch_map<AlignNode, uint8_t, AlignNodeHash> align_node_to_target_;
};


template <class BaseSeeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    auto mapped_batch = map_and_label_query_batch(generate_query);

    size_t i = 0;
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        const auto &[query_nodes_pair, target_columns] = mapped_batch;

        SeedAndExtendAlignerCore<AlignmentCompare> aligner_core(graph_, config_);
        DBGQueryAlignment paths(query, is_reverse_complement);
        std::string_view this_query = paths.get_query(is_reverse_complement);
        std::string_view reverse = paths.get_query(true);
        assert(this_query == query);

        assert(config_.num_alternative_paths);
        assert(target_columns[i].size());

        const auto &[nodes, nodes_rc] = query_nodes_pair[i];

        Extender extender(anno_graph_, config_, this_query);
        Extender extender_rc(anno_graph_, config_, reverse);

        for (const auto &[target_column, signature_pair] : target_columns[i]) {
            const auto &[signature, signature_rc] = signature_pair;

            Seeder seeder = build_seeder(
                target_column,
                this_query, // use this_query since paths stores a copy
                is_reverse_complement,
                std::vector<node_index>(nodes),
                signature
            );

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
                assert(!is_reverse_complement);

                auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                         const auto &callback) {
                    callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
                };

                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                Seeder seeder_rc = build_seeder(target_column,
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_both_directions(paths, seeder, seeder_rc,
                                                   extender, extender_rc,
                                                   build_rev_comp_alignment_core);

            } else if (config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);

                Seeder seeder_rc = build_seeder(target_column,
                                                reverse,
                                                !is_reverse_complement,
                                                std::vector<node_index>(nodes_rc),
                                                signature_rc);

                aligner_core.align_best_direction(paths, seeder, seeder_rc,
                                                  extender, extender_rc);

            } else {
                aligner_core.align_one_direction(paths, is_reverse_complement,
                                                 seeder, extender);
            }

            callback(header, std::move(paths));
        }

        ++i;
    });
}

template <class BaseSeeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::build_seeder(uint64_t target_column,
               std::string_view query,
               bool is_reverse_complement,
               std::vector<node_index>&& base_nodes,
               const sdsl::bit_vector &signature) const -> Seeder {
    assert(base_nodes.size() == signature.size());

    // mask out nodes not in present in target column
    if (is_reverse_complement) {
        for (size_t i = 0; i < base_nodes.size(); ++i) {
            if (!signature[base_nodes.size() - i - 1])
                base_nodes[i] = DeBruijnGraph::npos;
        }
    } else {
        for (size_t i = 0; i < base_nodes.size(); ++i) {
            if (!signature[i])
                base_nodes[i] = DeBruijnGraph::npos;
        }
    }

    return Seeder(target_column, graph_, query, is_reverse_complement,
                  std::move(base_nodes), config_);
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
