#ifndef __ALIGNER_SEEDER_METHODS_HPP__
#define __ALIGNER_SEEDER_METHODS_HPP__

#include "aligner_alignment.hpp"
#include "common/vectors/bitmap.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace align {

class ISeeder {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment Seed;

    virtual ~ISeeder() {}

    virtual std::vector<Seed> get_seeds() const = 0;
};

class ManualSeeder : public ISeeder {
  public:
    ManualSeeder(std::vector<Seed>&& seeds) : seeds_(std::move(seeds)) {}

    virtual ~ManualSeeder() {}

    std::vector<Seed> get_seeds() const override { return seeds_; }

  private:
    std::vector<Seed> seeds_;
};

class ExactSeeder : public ISeeder {
  public:
    typedef DBGAlignerConfig::score_t score_t;

    ExactSeeder(const DeBruijnGraph &graph,
                std::string_view query,
                bool orientation,
                std::vector<node_index>&& nodes,
                const DBGAlignerConfig &config);

    virtual ~ExactSeeder() {}

    std::vector<Seed> get_seeds() const override;

  protected:
    const DeBruijnGraph &graph_;
    std::string_view query_;
    bool orientation_;
    std::vector<node_index> query_nodes_;
    const DBGAlignerConfig &config_;
    std::vector<score_t> partial_sum_;
    size_t num_matching_;
};

class MEMSeeder : public ExactSeeder {
  public:
    template <typename... Args>
    MEMSeeder(Args&&... args) : ExactSeeder(std::forward<Args>(args)...) {}

    virtual ~MEMSeeder() {}

    std::vector<Seed> get_seeds() const override;

    virtual const bitmap& get_mem_terminator() const = 0;
};

class UniMEMSeeder : public MEMSeeder {
  public:
    template <typename... Args>
    UniMEMSeeder(Args&&... args)
          : MEMSeeder(std::forward<Args>(args)...),
            is_mem_terminus_([&](auto i) {
                                 return graph_.has_multiple_outgoing(i)
                                     || graph_.indegree(i) > 1;
                             },
                             graph_.max_index() + 1) {
        assert(is_mem_terminus_.size() == graph_.max_index() + 1);
    }

    virtual ~UniMEMSeeder() {}

    const bitmap& get_mem_terminator() const override { return is_mem_terminus_; }

  private:
    bitmap_lazy is_mem_terminus_;
};

template <class BaseSeeder>
class SuffixSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::Seed Seed;
    typedef typename BaseSeeder::node_index node_index;
    typedef typename BaseSeeder::score_t score_t;

    template <typename... Args>
    SuffixSeeder(Args&&... args)
          : BaseSeeder(std::forward<Args>(args)...),
            dbg_succ_(get_base_dbg_succ(this->graph_)) {
        assert(this->config_.min_seed_length < this->graph_.get_k());
    }

    virtual ~SuffixSeeder() {}

    std::vector<Seed> get_seeds() const override;

    BaseSeeder& get_base_seeder() { return dynamic_cast<BaseSeeder&>(*this); }

  private:
    const DBGSuccinct &dbg_succ_;

    static const DBGSuccinct& get_base_dbg_succ(const DeBruijnGraph &graph);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_SEEDER_METHODS_HPP__
