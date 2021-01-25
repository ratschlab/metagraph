#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__

#include <queue>

#include <priority_deque.hpp>

#include "aligner_helper.hpp"
#include "common/utils/template_utils.hpp"
#include "common/aligned_vector.hpp"
#include "common/vectors/bitmap.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace align {

/******************* Abstract Seeder classes *******************/

template <typename NodeType = typename DeBruijnGraph::node_index>
class ISeeder {
  public:
    typedef Alignment<NodeType> Seed;

    virtual ~ISeeder() {}

    virtual void call_seeds(std::function<void(Seed&&)> callback) const = 0;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactMapSeeder : public ISeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;

    ExactMapSeeder(const DeBruijnGraph &graph,
                   std::string_view query,
                   bool orientation,
                   std::vector<NodeType>&& nodes,
                   const DBGAlignerConfig &config);

    virtual ~ExactMapSeeder() {}

  protected:
    const DeBruijnGraph &graph_;
    std::string_view query_;
    bool orientation_;
    std::vector<NodeType> query_nodes_;
    const DBGAlignerConfig &config_;
    std::vector<DBGAlignerConfig::score_t> partial_sum_;
    size_t num_matching_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class MEMSeeder : public ExactMapSeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;

    template <typename... Args>
    MEMSeeder(Args&&... args) : ExactMapSeeder<NodeType>(std::forward<Args>(args)...) {}

    virtual ~MEMSeeder() {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;

    virtual const bitmap& get_mem_terminator() const = 0;
};


/************************ Seeder classes ***********************/

template <typename NodeType = typename DeBruijnGraph::node_index>
class ManualSeeder : public ISeeder<NodeType> {
  public:
    typedef NodeType node_index;
    typedef Alignment<NodeType> Seed;

    ManualSeeder(std::vector<Seed>&& seeds) : seeds_(std::move(seeds)) {}

    void call_seeds(std::function<void(Seed&&)> callback) const override {
        for (const Seed &seed : seeds_) {
            callback(Seed(seed));
        }
    }

  private:
    std::vector<Seed> seeds_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactSeeder : public ExactMapSeeder<NodeType> {
  public:
    typedef NodeType node_index;
    typedef typename ISeeder<NodeType>::Seed Seed;

    template <typename... Args>
    ExactSeeder(Args&&... args) : ExactMapSeeder<NodeType>(std::forward<Args>(args)...) {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class UniMEMSeeder : public MEMSeeder<NodeType> {
  public:
    typedef NodeType node_index;
    typedef typename ISeeder<NodeType>::Seed Seed;

    template <typename... Args>
    UniMEMSeeder(Args&&... args)
          : MEMSeeder<NodeType>(std::forward<Args>(args)...),
            is_mem_terminus_([&](auto i) {
                                 return this->graph_.has_multiple_outgoing(i)
                                     || this->graph_.indegree(i) > 1;
                             },
                             this->graph_.max_index() + 1) {
        assert(is_mem_terminus_.size() == this->graph_.max_index() + 1);
    }

    const bitmap& get_mem_terminator() const override { return is_mem_terminus_; }

  private:
    bitmap_lazy is_mem_terminus_;
};

template <class BaseSeeder>
class SuffixSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::node_index node_index;
    typedef typename BaseSeeder::Seed Seed;

    template <typename... Args>
    SuffixSeeder(Args&&... args)
          : BaseSeeder(std::forward<Args>(args)...),
            dbg_succ_(dynamic_cast<const DBGSuccinct&>(this->graph_)) {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;

  private:
    const DBGSuccinct &dbg_succ_;
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class IExtender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~IExtender() {}

    virtual void
    operator()(std::function<void(DBGAlignment&&, NodeType)> callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &seed) = 0;

  protected:
    virtual void reset() = 0;
    virtual const DBGAlignment& get_seed() const = 0;
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class DefaultColumnExtender : public IExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;

    typedef std::tuple<NodeType, score_t, bool /* converged */> ColumnRef;
    typedef boost::container::priority_deque<ColumnRef,
                                             std::vector<ColumnRef>,
                                             utils::LessSecond> ColumnQueue;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);
    DefaultColumnExtender(DefaultColumnExtender&&) = default;

    virtual ~DefaultColumnExtender() {}

    virtual void
    operator()(std::function<void(DBGAlignment&&, NodeType)> callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &path) override;

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
    emplace_node(NodeType node, NodeType incoming_node, char c, size_t size,
                 size_t best_pos = 0, size_t last_priority_pos = 0,
                 size_t begin = 0, size_t end = std::numeric_limits<size_t>::max());

    virtual bool add_seed(size_t clipping);

    virtual const DBGAlignment& get_seed() const override { return *path_; }

    virtual void check_and_push(ColumnRef&& next_column);

    void extend_main(std::function<void(DBGAlignment&&, NodeType)> callback,
                     score_t min_path_score);

    void update_columns(NodeType incoming_node,
                        const std::deque<std::pair<NodeType, char>> &out_columns,
                        score_t min_path_score);

    virtual std::deque<std::pair<NodeType, char>>
    fork_extension(NodeType /* fork after this node */,
                   std::function<void(DBGAlignment&&, NodeType)>,
                   score_t);

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op;

    // the initial seed
    const DBGAlignment *path_;

    // starting position of the alignment
    const char *align_start;

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
