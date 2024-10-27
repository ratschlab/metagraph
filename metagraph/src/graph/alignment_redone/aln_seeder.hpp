#pragma once

#include <tsl/hopscotch_set.h>

#include "aln_match.hpp"
#include "aln_query.hpp"
#include "aligner_config.hpp"

namespace mtg::graph::align_redone {

class Seeder {
  public:
    Seeder(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    virtual ~Seeder() {}

    virtual std::vector<Anchor> get_anchors() const = 0;
    virtual std::vector<Alignment> get_inexact_anchors() const = 0;

    const Query& get_query() const { return query_; }

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;
};

class ExactSeeder : public Seeder {
  public:
    template <typename... Args>
    ExactSeeder(Args&&... args) : Seeder(std::forward<Args>(args)...) {}

    std::vector<Anchor> get_anchors() const override;
    std::vector<Alignment> get_inexact_anchors() const override;

  private:
    tsl::hopscotch_map<DeBruijnGraph::node_index, tsl::hopscotch_map<DeBruijnGraph::node_index, std::vector<std::tuple<size_t, size_t, DBGAlignerConfig::score_t>>>>
    get_node_dists(std::vector<Anchor> &anchors) const;
};

class Extender {
  public:
    Extender(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    void extend(const Alignment &aln, const std::function<void(Alignment&&)> &callback) const;

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;
};

} // namespace mtg::graph::align