#ifndef __LABELED_ALIGNER_HPP__
#define __LABELED_ALIGNER_HPP__

#include <tsl/hopscotch_map.h>

#include "dbg_aligner.hpp"
#include "common/vector_set.hpp"
#include "common/utils/template_utils.hpp"
#include "common/hashers/hash.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace graph {
namespace align {

class ILabeledAligner : public ISeedAndExtendAligner {
  public:
    ILabeledAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : anno_graph_(anno_graph), graph_(anno_graph_.get_graph()), config_(config) {}

    virtual ~ILabeledAligner() {}

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
};

template <typename NodeType = DeBruijnGraph::node_index>
class LabeledBacktrackingExtender : public DefaultColumnExtender<NodeType> {
  public:
    typedef DefaultColumnExtender<DeBruijnGraph::node_index> BaseExtender;
    typedef typename BaseExtender::score_t score_t;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::AlignNode AlignNode;
    typedef typename BaseExtender::AlignNodeHash AlignNodeHash;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;

    LabeledBacktrackingExtender(const AnnotatedDBG &anno_graph,
                                const DBGAlignerConfig &config,
                                std::string_view query)
          : BaseExtender(anno_graph.get_graph(), config, query), anno_graph_(anno_graph) {
        targets_set_.emplace(Vector<uint64_t>{});
    }

    virtual ~LabeledBacktrackingExtender() {}

  protected:
    virtual std::vector<AlignNode>
    backtrack(score_t min_path_score,
              AlignNode best_node,
              tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
              std::vector<DBGAlignment> &extensions) const override;

    virtual void init_backtrack() const override;

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &,
                                      const AlignNode &) const override { return false; }

  private:
    const AnnotatedDBG &anno_graph_;
    mutable VectorSet<Vector<uint64_t>, utils::VectorHash> targets_set_;
    mutable tsl::hopscotch_map<node_index, size_t> targets_;
};

template <class Seeder = ExactSeeder<>,
          class Extender = LabeledBacktrackingExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class LabeledAligner : public ILabeledAligner {
  public:
    template <typename... Args>
    LabeledAligner(Args&&... args) : ILabeledAligner(std::forward<Args>(args)...) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final {
        generate_query([&](std::string_view header,
                           std::string_view query,
                           bool is_reverse_complement) {
            SeedFilter seed_filter(graph_.get_k());
            SeedAndExtendAlignerCore<AlignmentCompare> aligner_core(
                graph_, config_, seed_filter, query, is_reverse_complement
            );
            auto &paths = aligner_core.get_paths();
            std::string_view this_query = paths.get_query(is_reverse_complement);
            std::string_view reverse = paths.get_query(!is_reverse_complement);
            assert(this_query == query);

            std::vector<node_index> nodes = map_sequence_to_nodes(graph_, query);
            std::vector<node_index> nodes_rc;
            if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                    || config_.forward_and_reverse_complement) {
                assert(!is_reverse_complement);
                std::string dummy(query);
                nodes_rc = nodes;
                reverse_complement_seq_path(graph_, dummy, nodes_rc);
                assert(dummy == paths.get_query(true));
                assert(nodes_rc.size() == nodes.size());
            }

            std::shared_ptr<ISeeder<node_index>> seeder;
            std::shared_ptr<ISeeder<node_index>> seeder_rc;
            if (config_.min_seed_length < graph_.get_k()
                    && SuffixSeeder<Seeder>::get_base_dbg_succ(this->graph_)) {
                seeder = std::make_shared<SuffixSeeder<Seeder>>(
                    graph_, this_query, is_reverse_complement, std::move(nodes), config_
                );
            } else {
                seeder = std::make_shared<Seeder>(
                    graph_, this_query, is_reverse_complement, std::move(nodes), config_
                );
            }

            if (config_.forward_and_reverse_complement
                    || graph_.get_mode() == DeBruijnGraph::CANONICAL) {
                if (config_.min_seed_length < graph_.get_k()
                        && SuffixSeeder<Seeder>::get_base_dbg_succ(this->graph_)) {
                    seeder_rc = std::make_shared<SuffixSeeder<Seeder>>(
                        graph_, reverse, !is_reverse_complement, std::move(nodes_rc), config_
                    );
                } else {
                    seeder_rc = std::make_shared<Seeder>(
                        graph_, reverse, !is_reverse_complement, std::move(nodes_rc), config_
                    );
                }
            }

            auto extender = build_extender(this_query);

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
                auto extender_rc = build_extender(reverse);

                auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                         const auto &callback) {
                    callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
                };

                // From a given seed, align forwards, then reverse complement and
                // align backwards. The graph needs to be canonical to ensure that
                // all paths exist even when complementing.
                aligner_core.align_both_directions(*seeder, *seeder_rc,
                                                   *extender, *extender_rc,
                                                   build_rev_comp_alignment_core);

            } else if (config_.forward_and_reverse_complement) {
                auto extender_rc = build_extender(reverse);
                aligner_core.align_best_direction(*seeder, *seeder_rc, *extender, *extender_rc);

            } else {
                aligner_core.align_one_direction(is_reverse_complement, *seeder, *extender);
            }

            aligner_core.flush([&](const DBGAlignment &alignment) {
                return alignment.target_columns.empty();
            });

            callback(header, std::move(paths));
        });
    }

  protected:
    std::shared_ptr<IExtender<DeBruijnGraph::node_index>>
    build_extender(std::string_view query) const override {
        return std::make_shared<Extender>(anno_graph_, config_, query);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __LABELED_ALIGNER_HPP__
