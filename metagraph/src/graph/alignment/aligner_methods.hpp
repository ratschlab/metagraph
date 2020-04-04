#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__

#include "aligner_helper.hpp"
#include "common/vectors/bitmap.hpp"


template <typename NodeType = typename DeBruijnGraph::node_index>
class Seeder {
  public:
    virtual ~Seeder() {}

    virtual void initialize(std::string_view /* query string */,
                            bool /* orientation */) {}

    virtual void call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const = 0;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactSeeder : public Seeder<NodeType> {
  public:
    ExactSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) { assert(config_.check_config_scores()); }

    void initialize(std::string_view query, bool orientation);

    void call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const;

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }
    const std::string_view get_query() const { return query_; }
    const std::vector<NodeType>& get_query_nodes() const { return query_nodes_; }
    const std::vector<score_t>& get_partial_sums() const { return partial_sum; }
    bool get_orientation() const { return orientation_; }

  protected:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;
    bool orientation_;
    std::vector<score_t> partial_sum;
    std::vector<NodeType> query_nodes_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class SuffixSeeder : public Seeder<NodeType> {
  public:
    SuffixSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    void initialize(std::string_view query, bool orientation) {
        exact_seeder_.initialize(query, orientation);
    }

    void call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const;

    const DeBruijnGraph& get_graph() const { return exact_seeder_.get_graph(); }
    const DBGAlignerConfig& get_config() const { return exact_seeder_.get_config(); }

  private:
    ExactSeeder<NodeType> exact_seeder_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class MEMSeeder : public Seeder<NodeType> {
  public:
    virtual ~MEMSeeder() {}

    virtual void initialize(std::string_view query, bool orientation) override final;

    virtual void call_seeds(std::function<void(Alignment<NodeType>&&)> callback) const override final;

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }
    const bitmap& get_mem_terminator() const {
        assert(is_mem_terminus_.get());
        return *is_mem_terminus_;
    }

  protected:
    MEMSeeder(const DeBruijnGraph &graph,
              const DBGAlignerConfig &config,
              std::unique_ptr<bitmap>&& is_mem_terminus)
          : graph_(graph),
            config_(config),
            orientation_(false),
            is_mem_terminus_(std::move(is_mem_terminus)) {
        assert(is_mem_terminus_->size() == graph.max_index() + 1);
        assert(config_.check_config_scores());
    }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;
    std::vector<NodeType> query_nodes_;
    bool orientation_;
    std::vector<score_t> partial_sum;
    const std::unique_ptr<bitmap> is_mem_terminus_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class UniMEMSeeder : public MEMSeeder<NodeType> {
  public:
    UniMEMSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
        : MEMSeeder<NodeType>(graph, config,
                              std::make_unique<bitmap_lazy>(
                                      [&](auto i) {
                                          return graph.has_multiple_outgoing(i)
                                                  || graph.indegree(i) > 1;
                                      },
                                      graph.max_index() + 1)) {}
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class Extender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~Extender() {}

    virtual void
    operator()(const DBGAlignment &path,
               std::string_view query,
               std::function<void(DBGAlignment&&)> callback,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min()) = 0;

    virtual void initialize(const DBGAlignment &path) = 0;
    virtual void initialize_query(const std::string_view query) = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class DefaultColumnExtender : public Extender<NodeType> {
  public:
    typedef typename Extender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename Extender<NodeType>::node_index node_index;
    typedef typename Extender<NodeType>::score_t score_t;

    DefaultColumnExtender(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) { assert(config_.check_config_scores()); }

    void
    operator()(const DBGAlignment &path,
               std::string_view query,
               std::function<void(DBGAlignment&&)> callback,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min());

    void initialize(const DBGAlignment &) {}

    void initialize_query(const std::string_view query);

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;

    DPTable<NodeType> dp_table;

    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;
};


#endif // __DBG_ALIGNER_METHODS_HPP__
