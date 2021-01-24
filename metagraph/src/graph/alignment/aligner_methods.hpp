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
namespace align {

/******************* Abstract Seeder classes *******************/

template <typename NodeType = typename DeBruijnGraph::node_index>
class ISeeder {
  public:
    typedef Alignment<NodeType> Seed;

    virtual ~ISeeder() {}

    virtual void initialize(std::string_view /* query string */,
                            bool /* orientation */) {}

    virtual void call_seeds(std::function<void(Seed&&)> callback) const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactMapSeeder : public ISeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;
    typedef DBGAlignerConfig::score_t score_t;

    virtual ~ExactMapSeeder() {}

    virtual void initialize(std::string_view query, bool orientation) override;

    const std::string_view get_query() const { return query_; }
    const std::vector<NodeType>& get_query_nodes() const { return query_nodes_; }
    const std::vector<score_t>& get_partial_sums() const { return partial_sum_; }
    bool get_orientation() const { return orientation_; }

    uint64_t get_num_matching_kmers() const {
        return query_nodes_.size() - std::count(query_nodes_.begin(), query_nodes_.end(), NodeType());
    }

  private:
    std::string_view query_;
    bool orientation_;
    std::vector<score_t> partial_sum_;
    std::vector<NodeType> query_nodes_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class MEMSeeder : public ExactMapSeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;
    typedef DBGAlignerConfig::score_t score_t;

    void call_seeds(std::function<void(Seed&&)> callback) const;

    virtual const bitmap& get_mem_terminator() const = 0;
};


/************************ Seeder classes ***********************/

template <typename NodeType = typename DeBruijnGraph::node_index>
class ManualSeeder : public ISeeder<NodeType> {
  public:
    typedef Alignment<NodeType> Seed;

    ManualSeeder(const DeBruijnGraph &graph,
                 const DBGAlignerConfig &config,
                 std::vector<Seed>&& seeds)
          : graph_(graph), config_(config), seeds_(std::move(seeds)) {}

    void call_seeds(std::function<void(Seed&&)> callback) const override {
        for (const Seed &seed : seeds_) {
            callback(Seed(seed));
        }
    }

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::vector<Seed> seeds_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactSeeder : public ExactMapSeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;
    typedef DBGAlignerConfig::score_t score_t;

    ExactSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) { assert(config_.check_config_scores()); }

    void call_seeds(std::function<void(Seed&&)> callback) const override;

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class UniMEMSeeder : public MEMSeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;

    UniMEMSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config),
            is_mem_terminus_([&](auto i) {
                                 return graph_.has_multiple_outgoing(i)
                                     || graph_.indegree(i) > 1;
                             },
                             graph_.max_index() + 1) {
        assert(is_mem_terminus_.size() == graph_.max_index() + 1);
    }

    const DeBruijnGraph& get_graph() const override { return graph_; }
    const DBGAlignerConfig& get_config() const override { return config_; }

    const bitmap& get_mem_terminator() const override { return is_mem_terminus_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    bitmap_lazy is_mem_terminus_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class SuffixSeeder : public ISeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;

    SuffixSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    void initialize(std::string_view query, bool orientation) {
        base_seeder_->initialize(query, orientation);
    }

    void call_seeds(std::function<void(Seed&&)> callback) const;

    const DeBruijnGraph& get_graph() const { return base_seeder_->get_graph(); }
    const std::string_view get_query() const { return base_seeder_->get_query(); }
    const DBGAlignerConfig& get_config() const { return base_seeder_->get_config(); }
    const std::vector<NodeType>& get_query_nodes() const { return base_seeder_->get_query_nodes(); }
    bool get_orientation() const { return base_seeder_->get_orientation(); }

    uint64_t get_num_matching_kmers() const { return base_seeder_->get_num_matching_kmers(); }

  private:
    std::unique_ptr<ExactMapSeeder<NodeType>> base_seeder_;
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

    virtual void initialize(const DBGAlignment &path) = 0;
    virtual void initialize_query(const std::string_view query) = 0;
    virtual const std::string_view get_query() const = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;

  protected:
    virtual void set_graph(const DeBruijnGraph &graph) = 0;
    virtual void reset() = 0;
    virtual const DBGAlignment& get_seed() const = 0;
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class DefaultColumnExtender : public IExtender<NodeType> {
  public:
    typedef typename IExtender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename IExtender<NodeType>::node_index node_index;
    typedef typename IExtender<NodeType>::score_t score_t;
    typedef std::tuple<NodeType,
                       score_t,
                       bool /* converged */> ColumnRef;
    typedef boost::container::priority_deque<ColumnRef,
                                             std::vector<ColumnRef>,
                                             utils::LessSecond> ColumnQueue;

    DefaultColumnExtender(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(&graph), config_(config) { assert(config_.check_config_scores()); }

    DefaultColumnExtender(DefaultColumnExtender&&) = default;

    virtual ~DefaultColumnExtender() {}

    virtual void
    operator()(std::function<void(DBGAlignment&&, NodeType)> callback,
               score_t min_path_score = std::numeric_limits<score_t>::min()) override;

    virtual void initialize(const DBGAlignment &path) override;

    virtual void initialize_query(const std::string_view query) override;
    const std::string_view get_query() const override { return query; }

    const DeBruijnGraph& get_graph() const override {
        assert(graph_);
        return *graph_;
    }

    const DBGAlignerConfig& get_config() const override { return config_; }

    const DPTable<NodeType>& get_dp_table() const { return dp_table; }
    const ColumnQueue& get_column_queue() const { return columns_to_update; }

  protected:
    virtual void set_graph(const DeBruijnGraph &graph) override { graph_ = &graph; }
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

    virtual DPTable<NodeType>& get_dp_table() { return dp_table; }
    virtual ColumnQueue& get_column_queue() { return columns_to_update; }

  private:
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;

    // keep track of which columns to use next
    ColumnQueue columns_to_update;

    DPTable<NodeType> dp_table;

    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    tsl::hopscotch_map<char, AlignedVector<int8_t>> profile_score;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op;

    std::string_view query;

    const DBGAlignment *path_;
    const char *align_start;
    size_t size;
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
