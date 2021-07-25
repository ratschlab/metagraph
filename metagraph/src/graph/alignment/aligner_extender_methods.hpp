#ifndef __DBG_EXTENDER_METHODS_HPP__
#define __DBG_EXTENDER_METHODS_HPP__

#include <tsl/hopscotch_map.h>

#include "aligner_alignment.hpp"
#include "common/aligned_vector.hpp"


namespace mtg {
namespace graph {
namespace align {

class IExtender {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    virtual ~IExtender() {}

    virtual std::vector<Alignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const Alignment &seed) = 0;

    virtual void
    call_visited_nodes(const std::function<void(node_index,
                                                size_t /* range begin */,
                                                size_t /* range end */)> &callback) const = 0;

  protected:
    virtual void reset() = 0;
    virtual const Alignment& get_seed() const = 0;
};


class DefaultColumnExtender : public IExtender {
  public:
    enum NodeId : uint8_t {
        NONE,
        PREV,
        CUR
    };

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual std::vector<Alignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const Alignment &seed) override;

    virtual void
    call_visited_nodes(const std::function<void(node_index,
                                                size_t /* range begin */,
                                                size_t /* range end */)> &callback) const override;

  protected:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;

    typedef std::tuple<node_index,
                       char /* last character of the node label */,
                       size_t /* copy number */,
                       size_t /* distance from origin */> AlignNode;

    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<NodeId> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec,
                       OpVec, OpVec, OpVec, AlignNode,
                       PrevVec, PrevVec,
                       size_t /* offset */,
                       size_t /* max_pos */> Scores;
    typedef std::pair<std::vector<Scores>, bool> Column;

    tsl::hopscotch_map<node_index, Column> table_;

    virtual void reset() override { table_.clear(); }

    virtual const Alignment& get_seed() const override { return *seed_; }

    virtual std::vector<std::pair<node_index, char>> get_outgoing(const AlignNode &node) const;

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;

    // the initial seed
    const Alignment *seed_;

    std::string_view extend_window_;

    // start of the partial sum table
    const score_t *match_score_begin_;

    static bool has_converged(const Column &column, const Scores &next);

    static void sanitize(Scores &scores);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
