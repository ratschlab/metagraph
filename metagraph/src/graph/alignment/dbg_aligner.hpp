#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include "aligner_helper.hpp"
#include "aligner_methods.hpp"
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
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment<node_index> DBGAlignment;
    typedef QueryAlignment<node_index> DBGQueryAlignment;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~ISeedAndExtendAligner() {}

    virtual const DBGAlignerConfig& get_config() const = 0;
};


template <class AlignmentCompare = LocalAlignmentLess<>>
class SeedAndExtendAlignerCore;

template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = LocalAlignmentLess<>>
class DBGAligner : public ISeedAndExtendAligner {
  public:
    typedef IDBGAligner::node_index node_index;
    typedef IDBGAligner::DBGAlignment DBGAlignment;
    typedef IDBGAligner::DBGQueryAlignment DBGQueryAlignment;
    typedef IDBGAligner::score_t score_t;
    typedef IDBGAligner::QueryGenerator QueryGenerator;
    typedef IDBGAligner::AlignmentCallback AlignmentCallback;

    DBGAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : config_(config), aligner_core_(graph, config_), graph_(graph) {}

    virtual void align_batch(const QueryGenerator &generate_query,
                             const AlignmentCallback &callback) const override final;

    virtual const DBGAlignerConfig& get_config() const override final { return config_; }

  protected:
    DBGAlignerConfig config_;
    SeedAndExtendAlignerCore<AlignmentCompare> aligner_core_;
    const DeBruijnGraph &graph_;
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
    typedef const std::function<void(const LocalAlignmentCallback&,
                                     const MinScoreComputer&)> AlignmentGenerator;

    typedef std::function<void(const ISeeder<node_index>&,
                               IExtender<node_index>&&)> AlignCoreCallback;
    typedef std::function<void(std::string_view /* reverse_query */,
                               const ISeeder<node_index> & /* forward seeder */,
                               std::vector<DBGAlignment>&& /* rev_comp_seeds */,
                               const AlignCoreCallback&)> AlignCoreGenerator;

    SeedAndExtendAlignerCore(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {
        assert(config_.num_alternative_paths);
        if (!config_.check_config_scores()) {
            throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
        }
    }

    // Align the query sequence in the given orientation (false is forward,
    // true is reverse complement)
    void align_one_direction(DBGQueryAlignment &paths,
                             bool orientation_to_align,
                             const ISeeder<node_index> &seeder,
                             IExtender<node_index>&& extender) const;

    // Align both the forward and reverse complement of the query sequence,
    // then report the best scoring alignment.
    void align_best_direction(DBGQueryAlignment &paths,
                              const ISeeder<node_index> &seeder,
                              const ISeeder<node_index> &seeder_rc,
                              IExtender<node_index>&& extender,
                              IExtender<node_index>&& extender_rc) const;

    // Align both forwards and backwards from a given seed. Procedure
    // 1. Given each seed, extend forward to produce an alignment A
    // 2. Reverse complement the alignment to get A', treated like a new seed
    // 3. Extend A' forwards
    // 4. Reverse complement A' to get the final alignment A''
    void align_both_directions(DBGQueryAlignment &paths,
                               const ISeeder<node_index> &forward_seeder,
                               IExtender<node_index>&& forward_extender,
                               const AlignCoreGenerator &rev_comp_core_generator) const;

  protected:
    // Generate seeds, then extend them
    void align_core(std::string_view query,
                    const ISeeder<node_index> &seeder,
                    IExtender<node_index>&& extender,
                    const LocalAlignmentCallback &callback,
                    const MinScoreComputer &get_min_path_score) const;

    // Given alignments generated by a generator, add them to a priority queue
    // and add the top ones to paths.
    virtual void
    align_aggregate(DBGQueryAlignment &paths,
                    const AlignmentGenerator &alignment_generator) const;

    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
};


template <class Seeder, class Extender, class AlignmentCompare>
inline void DBGAligner<Seeder, Extender, AlignmentCompare>
::align_batch(const QueryGenerator &generate_query,
              const AlignmentCallback &callback) const {
    generate_query([&](std::string_view header,
                       std::string_view query,
                       bool is_reverse_complement) {
        DBGQueryAlignment paths(query, is_reverse_complement);
        std::string_view this_query = paths.get_query(is_reverse_complement);
        assert(this_query == query);

        assert(graph_.num_alternative_paths);
        Seeder seeder(graph_, this_query, // use this_query since paths stores a copy
                      is_reverse_complement, map_sequence_to_nodes(graph_, query),
                      config_);

        Extender extender(graph_, config_, this_query);

        if (graph_.is_canonical_mode()) {
            assert(!is_reverse_complement);

            auto build_rev_comp_alignment_core = [&](std::string_view reverse,
                                                     const auto &,
                                                     auto&& rev_comp_seeds,
                                                     const auto &callback) {
                ManualSeeder<node_index> seeder_rc(std::move(rev_comp_seeds));
                callback(seeder_rc, Extender(graph_, config_, reverse));
            };

            // From a given seed, align forwards, then reverse complement and
            // align backwards. The graph needs to be canonical to ensure that
            // all paths exist even when complementing.
            aligner_core_.align_both_directions(paths, seeder, std::move(extender),
                                                build_rev_comp_alignment_core);
        } else if (config_.forward_and_reverse_complement) {
            assert(!is_reverse_complement);
            std::string_view reverse = paths.get_query(true);

            Seeder seeder_rc(graph_, reverse, !is_reverse_complement,
                             map_sequence_to_nodes(graph_, reverse), config_);

            aligner_core_.align_best_direction(paths, seeder, seeder_rc,
                                               std::move(extender),
                                               Extender(graph_, config_, reverse));
        } else {
            aligner_core_.align_one_direction(paths, is_reverse_complement, seeder,
                                              std::move(extender));
        }

        callback(header, std::move(paths));
    });
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
