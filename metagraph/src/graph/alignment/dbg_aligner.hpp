#ifndef __DBG_ALIGNER_HPP__
#define __DBG_ALIGNER_HPP__

#include <cassert>
#include <functional>

#include "aligner_alignment.hpp"
#include "aligner_aggregator.hpp"
#include "aligner_seeder_methods.hpp"
#include "aligner_extender_methods.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {


class IDBGAligner {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    typedef std::tuple<std::string /* header */,
                       std::string /* seq */,
                       bool /* orientation of seq */> Query;
    typedef std::function<void(std::string_view /* header */,
                               QueryAlignment&& /* alignments */)> AlignmentCallback;

    virtual ~IDBGAligner() {}

    // Main aligner
    virtual void align_batch(const std::vector<Query> &seq_batch,
                             const AlignmentCallback &callback) const = 0;

    // Convenience method
    QueryAlignment align(std::string_view query, bool is_reverse_complement = false) const;
};

template <class AlignmentCompare = LocalAlignmentLess>
class ISeedAndExtendAligner : public IDBGAligner {
  public:
    ISeedAndExtendAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    virtual ~ISeedAndExtendAligner() {}

    virtual void align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
                             const AlignmentCallback &callback) const override;

    const DBGAlignerConfig& get_config() const { return config_; }

  protected:
    typedef AlignmentAggregator<AlignmentCompare> Aggregator;
    const DeBruijnGraph &graph_;

    virtual std::shared_ptr<IExtender>
    build_extender(std::string_view query, const Aggregator &aggregator) const = 0;

    virtual std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<node_index> &nodes) const = 0;

    template <class Seeder>
    std::shared_ptr<ISeeder>
    build_seeder_impl(std::string_view query,
                      bool is_reverse_complement,
                      const std::vector<node_index> &nodes) const {
        return std::make_shared<SuffixSeeder<Seeder>>(
            graph_, query, is_reverse_complement, nodes, config_
        );
    }

  private:
    DBGAlignerConfig config_;

    // Align the forward and reverse complement of the query sequence in both
    // directions and return the overall best alignment. e.g., for the forward query
    // 1. Find all seeds of its reverse complement
    // 2. Given a seed, extend forwards to get alignment A
    // 3. Reverse complement the alignment to get A', treat it like a new seed
    // 4. Extend A' forwards to get the final alignment A''
    void align_both_directions(std::string_view forward,
                               std::string_view reverse,
                               const ISeeder &forward_seeder,
                               const ISeeder &reverse_seeder,
                               IExtender &forward_extender,
                               IExtender &reverse_extender,
                               const std::function<void(Alignment&&)> &callback,
                               const std::function<score_t(const Alignment&)> &get_min_path_score) const;

    // Generates seeds and extends them
    void align_core(std::string_view query,
                    const ISeeder &seeder,
                    IExtender &extender,
                    const std::function<void(Alignment&&)> &callback,
                    const std::function<score_t(const Alignment&)> &get_min_path_score) const;
};

template <class Extender = DefaultColumnExtender,
          class Seeder = UniMEMSeeder,
          class AlignmentCompare = LocalAlignmentLess>
class DBGAligner : public ISeedAndExtendAligner<AlignmentCompare> {
  public:
    template <typename... Args>
    DBGAligner(Args&&... args)
          : ISeedAndExtendAligner<AlignmentCompare>(std::forward<Args>(args)...) {}

  private:
    typedef typename ISeedAndExtendAligner<AlignmentCompare>::Aggregator Aggregator;
    std::shared_ptr<IExtender>
    build_extender(std::string_view query, const Aggregator&) const override {
        return std::make_shared<Extender>(this->graph_, this->get_config(), query);
    }

    std::shared_ptr<ISeeder>
    build_seeder(std::string_view query,
                 bool is_reverse_complement,
                 const std::vector<IDBGAligner::node_index> &nodes) const override {
        return this->template build_seeder_impl<Seeder>(query, is_reverse_complement, nodes);
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_HPP__
