#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include <tsl/hopscotch_map.h>

#include "aligner_alignment.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "aligner_aggregator.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {


class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    typedef std::function<void(std::string_view /* header */,
                               std::string_view /* seq */,
                               bool /* orientation of seq */)> QueryCallback;
    typedef std::function<void(const QueryCallback&)> QueryGenerator;
    typedef std::function<void(std::string_view /* header */,
                               DBGQueryAlignment&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    // Main aligner
    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const = 0;

    // Convenience methods
    DBGQueryAlignment align(std::string_view query,
                            bool is_reverse_complement = false) const;
    void align_batch(const std::vector<std::pair<std::string, std::string>> &seq_batch,
                     const AlignmentCallback &callback) const;
};

class ISeedAndExtendAligner : public IDBGAligner {
  public:
    virtual ~ISeedAndExtendAligner() {}
    virtual const DBGAlignerConfig& get_config() const = 0;
};


template <class AlignmentCompare = LocalAlignmentLess>
class SeedAndExtendAlignerCore;

template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores())
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
    }

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
};

template <class AlignmentCompare>
class SeedAndExtendAlignerCore {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;

    typedef std::function<void(DBGAlignment&&)> LocalAlignmentCallback;
    typedef std::function<score_t(const DBGAlignment&)> MinScoreComputer;
    typedef std::function<void(const LocalAlignmentCallback&,
                               const MinScoreComputer&)> AlignmentGenerator;

    typedef std::function<void(const ISeeder<node_index>&)> SeederCallback;
    typedef std::function<void(std::vector<DBGAlignment>&& /* rev_comp_seeds */,
                               const SeederCallback&)> SeederGenerator;

    template <typename... Args>
    SeedAndExtendAlignerCore(const DeBruijnGraph &graph,
                             const DBGAlignerConfig &config,
                             Args&&... args)
          : graph_(graph), config_(config),
            paths_(std::forward<Args>(args)...),
            aggregator_(paths_.get_query(false), paths_.get_query(true), config_) {}

    void flush(const std::function<bool()> &terminate = []() { return false; }) {
        aggregator_.call_alignments([&](auto&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            paths_.emplace_back(std::move(alignment));
        }, terminate);
    }

    DBGQueryAlignment& get_paths() { return paths_; }

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    void align_one_direction(bool orientation_to_align,
                             const ISeeder<node_index> &seeder,
                             IExtender<node_index> &extender);

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    void align_best_direction(const ISeeder<node_index> &seeder,
                              const ISeeder<node_index> &seeder_rc,
                              IExtender<node_index> &extender,
                              IExtender<node_index> &extender_rc);

    // Align the forward and reverse complement of the query sequence in both
    // directions and return the overall best alignment. e.g., for the forward query
    // 1. Find all seeds of its reverse complement
    // 2. Given a seed, extend forwards to get alignment A
    // 3. Reverse complement the alignment to get A', treat it like a new seed
    // 4. Extend A' forwards to get the final alignment A''
    void align_both_directions(const ISeeder<node_index> &forward_seeder,
                               const ISeeder<node_index> &reverse_seeder,
                               IExtender<node_index> &forward_extender,
                               IExtender<node_index> &reverse_extender,
                               const SeederGenerator &rev_comp_core_generator);

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    void align_aggregate(const AlignmentGenerator &alignment_generator);

    AlignmentAggregator<node_index, AlignmentCompare>& get_aggregator() {
        return aggregator_;
    }

  protected:
    // Generate seeds, then extend them
    void align_core(std::string_view query,
                    const ISeeder<node_index> &seeder,
                    IExtender<node_index> &extender,
                    const LocalAlignmentCallback &callback,
                    const MinScoreComputer &get_min_path_score);

    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;

    DBGQueryAlignment paths_;
    AlignmentAggregator<node_index, AlignmentCompare> aggregator_;
};


template <class Seeder, class Extender, class AlignmentCompare>
inline void DBGAligner<Seeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        SeedAndExtendAlignerCore<AlignmentCompare> aligner_core(
            graph_, config_, query, is_reverse_complement
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

        Extender extender(graph_, config_, this_query);

        if (graph_.get_mode() == DeBruijnGraph::CANONICAL) {
            Extender extender_rc(graph_, config_, reverse);

            auto build_rev_comp_alignment_core = [&](auto&& rev_comp_seeds,
                                                     const auto &callback) {
                callback(ManualSeeder<node_index>(std::move(rev_comp_seeds)));
            };

            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            aligner_core.align_both_directions(*seeder, *seeder_rc,
                                               extender, extender_rc,
                                               build_rev_comp_alignment_core);

        } else if (config_.forward_and_reverse_complement) {
            Extender extender_rc(graph_, config_, reverse);
            aligner_core.align_best_direction(*seeder, *seeder_rc, extender, extender_rc);

        } else {
            aligner_core.align_one_direction(is_reverse_complement, *seeder, extender);
        }

        aligner_core.flush([this,&paths]() {
            assert(paths.size() <= config_.num_alternative_paths);
            return paths.size() == config_.num_alternative_paths;
        });

        callback(header, std::move(paths));
    });
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
