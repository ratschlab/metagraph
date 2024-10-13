#pragma once

#include <string_view>

#include "aln_cigar.hpp"
#include "aligner_config.hpp"

#include "graph/annotated_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"

namespace mtg::graph::align {

std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path);

class Match {
  public:
    using node_index = SequenceGraph::node_index;
    using score_t = DBGAlignerConfig::score_t;

    static constexpr score_t ninf = DBGAlignerConfig::ninf;

    virtual ~Match() {}

    std::string_view get_query() const { return query_; }
    std::string_view get_seed() const { return seed_; }

    const std::vector<node_index>& get_path() const { return path_; }

    virtual std::string_view get_spelling() const = 0;
    virtual std::string get_path_spelling() const = 0;
    virtual std::string_view get_trim_spelling() const = 0;

    score_t get_score() const { return score_; }
    bool get_orientation() const { return orientation_; }

    virtual size_t get_end_trim() const = 0;

    bool empty() const { return path_.empty(); }

    size_t get_clipping() const { return seed_.data() - query_.data(); }
    size_t get_end_clipping() const {
        return query_.data() + query_.size() - (seed_.data() + seed_.size());
    }

    bool is_spelling_valid(const DeBruijnGraph &graph) const {
        return get_path_spelling() == spell_path(graph, path_);
    }

  protected:
    Match(std::string_view query,
          size_t begin,
          size_t end,
          bool orientation,
          std::vector<node_index>&& path,
          score_t score)
          : query_(query), orientation_(orientation), seed_(query_.substr(begin, end - begin)),
            path_(std::move(path)), score_(score) {
        assert(is_spelling_valid());
    }

    std::string_view query_;
    bool orientation_;
    std::string_view seed_;
    std::vector<node_index> path_;
    score_t score_;
};

class Anchor : public Match {
  public:
    using label_class_t = uint64_t;
    static constexpr label_class_t nlabel = std::numeric_limits<uint64_t>::max();
    static constexpr int64_t ncoord = std::numeric_limits<int64_t>::max();

    Anchor(Anchor&&) = default;
    Anchor(const Anchor&) = default;
    Anchor& operator=(const Anchor&) = default;
    Anchor& operator=(Anchor&&) = default;

    Anchor(std::string_view query,
           size_t begin,
           size_t end,
           bool orientation,
           std::vector<node_index>&& path,
           const DBGAlignerConfig &config,
           std::string&& suffix = "",
           label_class_t label_class = nlabel,
           int64_t coord = ncoord)
          : Match(query, begin, end, orientation,
                  std::move(path), config.match_score(query)),
            label_class_(label_class), coord_(coord), suffix_(std::move(suffix)) {}

    std::string_view get_spelling() const override final { return seed_; }
    std::string get_path_spelling() const override final { return std::string(seed_) + suffix_; }
    std::string_view get_trim_spelling() const override final { return suffix_; }

    size_t get_end_trim() const override final { return suffix_.size(); }

    void append(const Anchor &other, const DeBruijnGraph *graph = nullptr);

    label_class_t get_label_class() const { return label_class_; }

  private:
    label_class_t label_class_;
    int64_t coord_;
    std::string suffix_;
};

class Alignment : public Match {
  public:
    Alignment(const Alignment&) = default;
    Alignment(Alignment&&) = default;
    Alignment& operator=(const Alignment&) = default;
    Alignment& operator=(Alignment&&) = default;

    Alignment(const Anchor &anchor)
          : Match(anchor.get_query(), anchor.get_seed().data() - anchor.get_query().data(),
                  anchor.get_query().data() + anchor.get_query().size() - (anchor.get_seed().data() + anchor.get_seed().size()),
                  anchor.get_orientation(),
                  std::vector<node_index>(anchor.get_path()),
                  anchor.get_score()),
            end_trim_(anchor.get_end_trim()),
            path_spelling_(anchor.get_path_spelling()),
            cigar_(std::string(get_clipping(), 'S')
                      + std::string(seed_.size(), '=')
                      + std::string(get_end_clipping(), 'S')),
            label_classes_(path_.size(), anchor.get_label_class()) {
        assert(path_spelling_.size() == query_.size() + end_trim_);
    }

    const Cigar& get_cigar() const { return cigar_; }

    size_t get_end_trim() const override final { return end_trim_; }

    std::string_view get_spelling() const override final {
        return std::string_view(path_spelling_.c_str(), path_spelling_.size() - end_trim_);
    }

    std::string get_path_spelling() const override final { return path_spelling_; }
    std::string_view get_trim_spelling() const override final {
        return std::string_view(path_spelling_.c_str() + path_spelling_.size() - end_trim_,
                                end_trim_);
    }

  private:
    size_t end_trim_;
    std::string path_spelling_;
    Cigar cigar_;
    std::vector<size_t> label_classes_;
};

} // namespace mtg::graph::align