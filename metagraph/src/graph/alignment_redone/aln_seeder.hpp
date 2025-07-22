#pragma once

#include <tsl/hopscotch_map.h>

#include "aln_match.hpp"
#include "aln_query.hpp"
#include "aligner_config.hpp"
#include "annotation_buffer.hpp"

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
};

class LabeledSeeder : public ExactSeeder {
  public:
    template <typename... Args>
    LabeledSeeder(AnnotationBuffer &anno_buffer, Args&&... args)
        : ExactSeeder(std::forward<Args>(args)...),
          anno_buffer_(anno_buffer) {}

    std::vector<Anchor> get_anchors() const override;

  private:
    AnnotationBuffer &anno_buffer_;
};

class Extender {
  public:
    Extender(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    void extend(const Alignment &aln,
                const std::function<void(Alignment&&)> &callback,
                bool no_bwd = false,
                bool no_fwd = false) const;

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;
};

} // namespace mtg::graph::align