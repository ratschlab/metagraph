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

    typedef std::tuple<NodeType, score_t, bool /* converged */> ColumnRef;
    typedef boost::container::priority_deque<ColumnRef,
                                             std::vector<ColumnRef>,
                                             utils::LessSecond> ColumnQueue;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual void
    operator()(ExtensionCallback callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &seed) override;

    const DPTable<NodeType>& get_dp_table() const { return dp_table; }

  protected:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query;

    // keep track of which columns to use next
    ColumnQueue columns_to_update;

    DPTable<NodeType> dp_table;

    virtual void reset() override { dp_table.clear(); }

    virtual std::pair<typename DPTable<NodeType>::iterator, bool>
    emplace_node(NodeType node,
                 NodeType incoming_node,
                 char c,
                 size_t size,
                 size_t best_pos = 0,
                 size_t last_priority_pos = 0,
                 size_t begin = 0,
                 size_t end = std::numeric_limits<size_t>::max());

    virtual bool add_seed(size_t clipping);

    virtual const DBGAlignment& get_seed() const override { return *path_; }

    virtual void check_and_push(ColumnRef&& next_column);

    void extend_main(ExtensionCallback callback, score_t min_path_score);

    void update_columns(NodeType incoming_node,
                        const std::deque<std::pair<NodeType, char>> &out_columns,
                        score_t min_path_score);

    virtual std::deque<std::pair<NodeType, char>>
    fork_extension(NodeType /* fork after this node */, ExtensionCallback, score_t);

    inline bool ram_limit_reached() const {
        return static_cast<double>(dp_table.num_bytes()) / 1024 / 1024
            > config_.max_ram_per_alignment;
    }

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op;

    // the initial seed
    const DBGAlignment *path_;

    std::string_view extend_window_;

    // max size of a column
    size_t size;

    // start of the partial sum table
    const score_t *match_score_begin;
    NodeType start_node;
    score_t start_score;
    score_t score_cutoff;
    size_t begin;
    size_t end;
    score_t xdrop_cutoff;
    bool overlapping_range_;
    size_t max_num_nodes;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_ALIGNER_METHODS_HPP__
