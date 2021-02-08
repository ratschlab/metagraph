#ifndef __ALIGNER_SEEDER_METHODS_HPP__
#define __ALIGNER_SEEDER_METHODS_HPP__

#include "aligner_alignment.hpp"
#include "common/vectors/bitmap.hpp"


namespace mtg {
namespace graph {

class DBGSuccinct;

namespace align {

template <typename NodeType = uint64_t>
class ISeeder {
  public:
    typedef Alignment<NodeType> Seed;

    virtual ~ISeeder() {}

    virtual void call_seeds(std::function<void(Seed&&)> callback) const = 0;
};

template <typename NodeType = uint64_t>
class ManualSeeder : public ISeeder<NodeType> {
  public:
    typedef NodeType node_index;
    typedef Alignment<NodeType> Seed;

    ManualSeeder(std::vector<Seed>&& seeds) : seeds_(std::move(seeds)) {}

    virtual ~ManualSeeder() {}

    void call_seeds(std::function<void(Seed&&)> callback) const override {
        for (const Seed &seed : seeds_) {
            callback(Seed(seed));
        }
    }

  private:
    std::vector<Seed> seeds_;
};

template <typename NodeType = uint64_t>
class ExactSeeder : public ISeeder<NodeType> {
  public:
    typedef NodeType node_index;
    typedef typename ISeeder<NodeType>::Seed Seed;

    ExactSeeder(const DeBruijnGraph &graph,
                std::string_view query,
                bool orientation,
                std::vector<NodeType>&& nodes,
                const DBGAlignerConfig &config);

    virtual ~ExactSeeder() {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;

  protected:
    const DeBruijnGraph &graph_;
    std::string_view query_;
    bool orientation_;
    std::vector<NodeType> query_nodes_;
    const DBGAlignerConfig &config_;
    std::vector<DBGAlignerConfig::score_t> partial_sum_;
    size_t num_matching_;
};

template <typename NodeType = uint64_t>
class MEMSeeder : public ExactSeeder<NodeType> {
  public:
    typedef typename ISeeder<NodeType>::Seed Seed;

    template <typename... Args>
    MEMSeeder(Args&&... args) : ExactSeeder<NodeType>(std::forward<Args>(args)...) {}

    virtual ~MEMSeeder() {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;

    virtual const bitmap& get_mem_terminator() const = 0;
};

template <typename NodeType = uint64_t>
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

    virtual ~UniMEMSeeder() {}

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
            dbg_succ_(get_base_dbg_succ(this->graph_)) {
        assert(this->config_.min_seed_length < this->graph_.get_k());
    }

    virtual ~SuffixSeeder() {}

    void call_seeds(std::function<void(Seed&&)> callback) const override;

    BaseSeeder& get_base_seeder() { return dynamic_cast<BaseSeeder&>(*this); }

  private:
    const DBGSuccinct &dbg_succ_;

    static const DBGSuccinct& get_base_dbg_succ(const DeBruijnGraph &graph);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_SEEDER_METHODS_HPP__
