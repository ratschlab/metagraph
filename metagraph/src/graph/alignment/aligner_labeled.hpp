#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"

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
    static constexpr uint64_t kNTarget = std::numeric_limits<uint64_t>::max();

    ILabeledDBGAligner(const AnnotatedDBG &anno_graph,
                       const DBGAlignerConfig &config,
                       size_t num_top_labels = 1);

    virtual ~ILabeledDBGAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    typedef std::vector<std::vector<node_index>> BatchMapping;
    typedef std::vector<tsl::hopscotch_map<uint64_t, sdsl::bit_vector>> BatchTargets;

    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
    size_t num_top_labels_;

    std::pair<BatchMapping, BatchTargets>
    map_query_batch(const QueryGenerator &generate_query) const;
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

    uint64_t get_target_column() const { return target_column_; }

  private:
    uint64_t target_column_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class LabeledColumnExtender;

template <class BaseSeeder = ExactSeeder<>,
          class Extender = LabeledColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess<>>
class LabeledDBGAligner : public ILabeledDBGAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    template <typename... Args>
    LabeledDBGAligner(Args&&... args)
          : ILabeledDBGAligner(std::forward<Args>(args)...),
            aligner_core_(graph_, config_) {}

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

  protected:
    typedef LabeledSeeder<BaseSeeder> Seeder;

    SeedAndExtendAlignerCore<AlignmentCompare> aligner_core_;

    Extender build_extender(std::string_view query,
                            const ISeeder<node_index> &seeder) const;

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
    typedef typename IExtender<NodeType>::ExtensionCallback ExtensionCallback;

    LabeledColumnExtender(const AnnotatedDBG &anno_graph,
                          const DBGAlignerConfig &config,
                          std::string_view query,
                          uint64_t target_column = ILabeledDBGAligner::kNTarget);

    virtual ~LabeledColumnExtender() {}

    virtual void
    operator()(ExtensionCallback callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &path) override;

    void set_target_column(uint64_t target_column) { target_column_ = target_column; }

  protected:
    typedef std::deque<std::pair<NodeType, char>> Edges;

    virtual Edges fork_extension(NodeType node,
                                 ExtensionCallback callback,
                                 score_t min_path_score) override;

  private:
    const AnnotatedDBG &anno_graph_;

    // default target column
    uint64_t main_target_column_;

    // Target column when extending the current seed. This is equal to main_target_column_
    // for all but suffix seeds
    uint64_t target_column_;

    // alternative seed used to replace a suffix seed
    DBGAlignment alt_seed_;

    // Extension to an input suffix seed s.t.
    // input_seed.append(seed_extension_) results in alt_seed_
    DBGAlignment seed_extension_;

    // the last node of the input suffix seed
    NodeType old_start_;

    tsl::hopscotch_map<NodeType, Edges> cached_edge_sets_;
};


template <class BaseSeeder, class Extender, class AlignmentCompare>
inline void LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    auto mapped_batch = map_query_batch(generate_query);

    size_t num_queries = 0;
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        const auto &[query_nodes, target_columns] = mapped_batch;
        DBGQueryAlignment paths(query, is_reverse_complement);
        std::string_view this_query = paths.get_query(is_reverse_complement);
        assert(this_query == query);

        assert(config_.num_alternative_paths);
        assert(target_columns[num_queries].size());

        for (const auto &[target_column, signature] : target_columns[num_queries]) {
            Seeder seeder = build_seeder(
                target_column,
                this_query, // use this_query since paths stores a copy
                is_reverse_complement,
                std::vector<node_index>(query_nodes[num_queries]),
                signature
            );

            Extender extender = build_extender(this_query, seeder);

            if (graph_.is_canonical_mode()) {
                assert(!is_reverse_complement);

                std::string_view reverse = paths.get_query(true);

                auto build_rev_comp_seeder = [&](const auto &forward_seeder,
                                                 auto&& rev_comp_seeds,
                                                 const auto &callback) {
                    callback(LabeledSeeder<ManualSeeder<node_index>>(
                        dynamic_cast<const Seeder&>(forward_seeder).get_target_column(),
                        std::move(rev_comp_seeds)
                    ));
                };

                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                aligner_core_.align_both_directions(paths, seeder, std::move(extender),
                                                    build_extender(reverse, seeder),
                                                    build_rev_comp_seeder);
            } else if (config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);

                std::string_view reverse = paths.get_query(true);

                Seeder seeder_rc = build_seeder(target_column,
                                                reverse,
                                                !is_reverse_complement,
                                                map_sequence_to_nodes(graph_, reverse),
                                                signature);

                aligner_core_.align_best_direction(paths, seeder, seeder_rc,
                                                   std::move(extender),
                                                   build_extender(reverse, seeder_rc));
            } else {
                aligner_core_.align_one_direction(paths, is_reverse_complement,
                                                  seeder, std::move(extender));
            }

            callback(header, std::move(paths));
        }

        ++num_queries;
    });

    assert(num_queries == mapped_batch.second.size());
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

template <class BaseSeeder, class Extender, class AlignmentCompare>
inline auto LabeledDBGAligner<BaseSeeder, Extender, AlignmentCompare>
::build_extender(std::string_view query,
                 const ISeeder<node_index> &seeder) const -> Extender {
    typedef LabeledSeeder<ManualSeeder<node_index>> LabeledManualSeeder;
    uint64_t target_column = kNTarget;

    if (const auto *labeled_seeder = dynamic_cast<const Seeder*>(&seeder)) {
        target_column = labeled_seeder->get_target_column();
    } else if (const auto *manual = dynamic_cast<const LabeledManualSeeder*>(&seeder)) {
        target_column = manual->get_target_column();
    }

    return { anno_graph_, config_, query, target_column };
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
