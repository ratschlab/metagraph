#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "aligner_alignment.hpp"
#include "common/aligned_vector.hpp"


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

    virtual ~IExtender() {}

    virtual std::vector<DBGAlignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &seed) = 0;

    virtual void
    call_visited_nodes(const std::function<void(NodeType,
                                                size_t /* range begin */,
                                                size_t /* range end */)> &callback) const = 0;

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

    enum NodeId : uint8_t {
        NONE,
        PREV,
        CUR
    };

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual std::vector<DBGAlignment>
    get_extensions(score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &seed) override;

    virtual void
    call_visited_nodes(const std::function<void(NodeType,
                                                size_t /* range begin */,
                                                size_t /* range end */)> &callback) const override;

  protected:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;
    std::string_view extend_window_;

    typedef std::tuple<NodeType,
                       char /* last character of the node label */,
                       size_t /* copy number */,
                       size_t /* distance from origin */> AlignNode;

    AlignNode start_node_;

    typedef AlignedVector<score_t> ScoreVec;
    typedef AlignedVector<NodeId> PrevVec;
    typedef AlignedVector<Cigar::Operator> OpVec;
    typedef std::tuple<ScoreVec, ScoreVec, ScoreVec,
                       OpVec, OpVec, OpVec, AlignNode,
                       PrevVec, PrevVec,
                       size_t /* offset */,
                       size_t /* max_pos */> Scores;
    typedef std::pair<std::vector<Scores>, bool> Column;

    struct AlignNodeHash {
        uint64_t operator()(const AlignNode &x) const {
            uint64_t seed = hasher1(std::get<0>(x));
            return seed ^ (hasher2(std::get<2>(x)) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
        }

        std::hash<NodeType> hasher1;
        std::hash<size_t> hasher2;
    };

    tsl::hopscotch_map<NodeType, Column> table_;

    // the initial seed
    const DBGAlignment *seed_;

    score_t xdrop_cutoff_;
    size_t start_;

    virtual void reset() override { table_.clear(); }

    virtual const DBGAlignment& get_seed() const override { return *seed_; }

    virtual std::vector<AlignNode> get_outgoing(const AlignNode &node) const;

    virtual void add_scores_to_column(Column &column, Scores&& scores, const AlignNode&) {
        column.first.emplace_back(std::move(scores));
    }

    virtual void backtrack(score_t min_path_score,
                           AlignNode best_node,
                           tsl::hopscotch_set<AlignNode, AlignNodeHash> &prev_starts,
                           std::vector<DBGAlignment> &extensions) const;

    virtual bool skip_backtrack_start(const std::vector<DBGAlignment> &extensions,
                                      const AlignNode &node) const;

    static std::pair<size_t, size_t> get_band(const AlignNode &prev,
                                              const Column &column_prev,
                                              score_t xdrop_cutoff);

    virtual void pop(const AlignNode &) {}
    virtual void init_backtrack() const {}

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;

    static bool has_converged(const Column &column, const Scores &next);

    static void sanitize(Scores &scores, score_t max_score);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_METHODS_HPP__
