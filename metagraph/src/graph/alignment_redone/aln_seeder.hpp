#pragma once

#include <tsl/hopscotch_map.h>

#include "aln_match.hpp"
#include "aln_query.hpp"
#include "aligner_config.hpp"
#include "annotation_buffer.hpp"

namespace mtg::graph::align_redone {

class AlignmentGraph;

class Seeder {
  public:
    Seeder(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    virtual ~Seeder() {}

    virtual std::vector<Anchor> get_anchors() const = 0;
    virtual std::vector<Alignment> get_inexact_anchors(bool align = true) const = 0;

    const Query& get_query() const { return query_; }

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;
};

class ExactSeeder : public Seeder {
  public:
    template <typename... Args>
    ExactSeeder(Args&&... args) : Seeder(std::forward<Args>(args)...) {}

    virtual ~ExactSeeder() {}

    std::vector<Anchor> get_anchors() const override;
    std::vector<Alignment> get_inexact_anchors(bool align = true) const override;

  protected:
    virtual AlignmentGraph make_aln_graph(Anchor::label_class_t) const;
};

class LabeledSeeder : public ExactSeeder {
  public:
    template <typename... Args>
    LabeledSeeder(AnnotationBuffer &anno_buffer, Args&&... args)
        : ExactSeeder(std::forward<Args>(args)...),
          anno_buffer_(anno_buffer) {}

    std::vector<Anchor> get_anchors() const override;
    AnnotationBuffer& get_buffer() const { return anno_buffer_; }

  protected:
    virtual AlignmentGraph make_aln_graph(Anchor::label_class_t target) const override;

  private:
    AnnotationBuffer &anno_buffer_;
};

class Extender {
  public:
    Extender(const Query &query, const DBGAlignerConfig &config)
          : query_(query), config_(config) {}

    virtual ~Extender() {}

    void extend(const Alignment &aln,
                const std::function<void(Alignment&&)> &callback,
                bool no_bwd = false,
                bool no_fwd = false) const;

  protected:
    const Query &query_;
    const DBGAlignerConfig &config_;

    virtual AlignmentGraph make_aln_graph(Anchor::label_class_t) const;
};

class LabeledExtender : public Extender {
  public:
    template <typename... Args>
    LabeledExtender(AnnotationBuffer &anno_buffer, Args&&... args)
        : Extender(std::forward<Args>(args)...),
          anno_buffer_(anno_buffer) {}

  private:
    AnnotationBuffer &anno_buffer_;

    virtual AlignmentGraph make_aln_graph(Anchor::label_class_t target) const override;
};

void align_query(const Query &query,
                 const Seeder &seeder,
                 const Extender &extender,
                 const std::function<void(Alignment&&)> &callback,
                 bool connect_anchors_in_chain = true);

void align_query(const Query &query,
                 const Seeder &seeder,
                 const std::function<void(Alignment&&)> &callback,
                 bool connect_anchors_in_chain = true);

} // namespace mtg::graph::align