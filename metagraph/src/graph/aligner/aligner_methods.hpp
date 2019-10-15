#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__


#include "aligner_helper.hpp"
#include "bitmap.hpp"


template <typename NodeType = typename DeBruijnGraph::node_index>
class Seeder {
  public:
    virtual ~Seeder() {}

    virtual std::vector<Alignment<NodeType>>
    operator()(const char *begin, const char *end,
               size_t clipping = 0,
               bool orientation = false) const = 0;

    virtual void initialize(const std::string& /* query string */,
                            bool /* orientation */) {}
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class ExactSeeder : public Seeder<NodeType> {
  public:
    ExactSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {}

    std::vector<Alignment<NodeType>> operator()(const char *begin, const char *end,
                                                size_t clipping = 0,
                                                bool orientation = false) const;

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class SuffixSeeder : public Seeder<NodeType> {
  public:
    SuffixSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : exact_seeder_(graph, config) {}

    std::vector<Alignment<NodeType>> operator()(const char *begin, const char *end,
                                                size_t clipping = 0,
                                                bool orientation = false) const;

    const DeBruijnGraph& get_graph() const { return exact_seeder_.get_graph(); }
    const DBGAlignerConfig& get_config() const { return exact_seeder_.get_config(); }

  private:
    const ExactSeeder<NodeType> exact_seeder_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class MEMSeeder : public Seeder<NodeType> {
  public:
    std::vector<Alignment<NodeType>> operator()(const char *begin, const char *end,
                                                size_t clipping = 0,
                                                bool orientation = false) const override;

    virtual void initialize(const std::string &query, bool orientation) override final {
        orientation_ = orientation;
        query_nodes_.clear();
        query_nodes_.reserve(query.size() - graph_.get_k() + 1);
        graph_.map_to_nodes_sequentially(query.begin(), query.end(),
                                         [&](auto i) { query_nodes_.push_back(i); });
    }

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
        assert(is_mem_terminus_->size() == graph.num_nodes() + 1);
    }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
    std::vector<NodeType> query_nodes_;
    bool orientation_;
    const std::unique_ptr<bitmap> is_mem_terminus_;
};

template <typename NodeType = typename DeBruijnGraph::node_index>
class UniMEMSeeder : public MEMSeeder<NodeType> {
  public:
    UniMEMSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : MEMSeeder<NodeType>(graph,
                                config,
                                std::make_unique<bitmap_lazy>([&](auto i) {
                                    return graph.outdegree(i) > 1 || graph.indegree(i) > 1;
                                },
                                graph.num_nodes() + 1)) {}
};


template <typename NodeType = typename DeBruijnGraph::node_index>
class Extender {
  public:
    typedef Alignment<NodeType> DBGAlignment;
    typedef typename DBGAlignment::node_index node_index;
    typedef typename DBGAlignment::score_t score_t;

    virtual ~Extender() {}

    virtual std::vector<DBGAlignment>
    operator()(const DBGAlignment &path,
               const char *sequence_end_iterator,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min()) const = 0;

    virtual void initialize(const DBGAlignment &path) = 0;

    virtual const DeBruijnGraph& get_graph() const = 0;
    virtual const DBGAlignerConfig& get_config() const = 0;
};


template <typename NodeType = typename DeBruijnGraph::node_index,
          class Compare = std::less<typename Alignment<NodeType>::Column>>
class DefaultColumnExtender : public Extender<NodeType> {
  public:
    typedef typename Extender<NodeType>::DBGAlignment DBGAlignment;
    typedef typename Extender<NodeType>::node_index node_index;
    typedef typename Extender<NodeType>::score_t score_t;
    typedef typename DBGAlignment::DPTable DPTable;
    typedef typename DBGAlignment::Step Step;
    typedef typename DBGAlignment::Column Column;

    DefaultColumnExtender(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : graph_(graph), config_(config) {}

    std::vector<DBGAlignment>
    operator()(const DBGAlignment &path,
               const char *sequence_end_iterator,
               bool orientation,
               score_t min_path_score = std::numeric_limits<score_t>::min()) const;

    void initialize(const DBGAlignment &) {}

    const DeBruijnGraph& get_graph() const { return graph_; }
    const DBGAlignerConfig& get_config() const { return config_; }

  private:
    const DeBruijnGraph &graph_;
    const DBGAlignerConfig &config_;
};

#endif // __DBG_ALIGNER_METHODS_HPP__
