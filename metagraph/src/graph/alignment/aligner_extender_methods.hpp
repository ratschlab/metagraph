#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__

#include <queue>

#include <priority_deque.hpp>

#include "aligner_alignment.hpp"
#include "aligner_dp_table.hpp"
#include "common/aligned_vector.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace align {

template <typename NodeType = uint64_t>
class IExtender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;
    typedef std::function<void(DBGAlignment&&, NodeType)> ExtensionCallback;

    virtual ~IExtender() {}

    virtual void
    operator()(ExtensionCallback callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &seed) = 0;

    virtual void call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &callback) const = 0;

  protected:
    virtual void reset() = 0;
    virtual const DBGAlignment& get_seed() const = 0;
};


template <typename NodeType = uint64_t>
class DefaultColumnExtender : public IExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;
    typedef typename IExtender<NodeType>::ExtensionCallback ExtensionCallback;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual void
    operator()(ExtensionCallback callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &seed) override;

    virtual void call_visited_nodes(const std::function<void(NodeType, size_t, size_t)> &callback) const override;

  protected:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;

    typedef std::pair<NodeType, size_t> AlignNode;
    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<AlignNode> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec,
                       OpVec, OpVec, OpVec,
                       PrevVec, PrevVec,
                       size_t /* offset */,
                       size_t /* max_pos */> Scores;
    typedef std::pair<std::vector<Scores>, bool> Column;

    tsl::hopscotch_map<NodeType, Column> table_;

    virtual void reset() override { table_.clear(); }

    virtual const DBGAlignment& get_seed() const override { return *seed_; }

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;

    // the initial seed
    const DBGAlignment *seed_;

    std::string_view extend_window_;

    // start of the partial sum table
    const score_t *match_score_begin_;

    static bool has_converged(const Column &column, size_t begin, size_t end);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_METHODS_HPP__
