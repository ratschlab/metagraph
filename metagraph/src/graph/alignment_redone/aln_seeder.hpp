#pragma once

#include "aln_match.hpp"
#include "aln_query.hpp"
#include "aligner_config.hpp"

namespace mtg::graph::align {

class Seeder {
  public:
    Seeder(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    virtual ~Seeder() {}

    virtual std::vector<Anchor> get_anchors() const = 0;
    virtual std::vector<Alignment> get_inexact_anchors() const;

    const Query& get_query() const { return query_; }

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;
};

class ExactSeeder : public Seeder {
  public:
    std::vector<Anchor> get_anchors() const override;
    std::vector<Alignment> get_inexact_anchors() const override;
};

} // namespace mtg::graph::align