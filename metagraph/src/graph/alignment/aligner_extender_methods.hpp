#ifndef __DBG_EXTENDER_METHODS_HPP__
#define __DBG_EXTENDER_METHODS_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "aligner_alignment.hpp"
#include "common/aligned_vector.hpp"


namespace mtg {
namespace graph {
namespace align {

class IExtender {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    virtual ~IExtender() {}

    std::vector<Alignment>
    get_extensions(const Alignment &seed,
                   score_t min_path_score = std::numeric_limits<score_t>::min()) {
        return set_seed(seed) ? extend(min_path_score) : std::vector<Alignment>{};
    }

    virtual void set_graph(const DeBruijnGraph &graph) = 0;

    virtual size_t num_explored_nodes() const = 0;

  protected:
    virtual const Alignment& get_seed() const = 0;
    virtual bool set_seed(const Alignment &seed) = 0;

    virtual std::vector<Alignment> extend(score_t min_path_score) = 0;
};

class SeedFilteringExtender : public IExtender {
  public:
    SeedFilteringExtender(std::string_view query) : query_size_(query.size()) {}

    virtual ~SeedFilteringExtender() {}

    virtual void set_graph(const DeBruijnGraph &) override { conv_checker_.clear(); }

    virtual size_t num_explored_nodes() const override { return conv_checker_.size(); }

  protected:
    const Alignment *seed_ = nullptr;
    size_t query_size_;

    typedef std::pair<size_t, AlignedVector<score_t>> ScoreVec;
    tsl::hopscotch_map<node_index, ScoreVec> conv_checker_;

    virtual const Alignment& get_seed() const override final { return *seed_; }
    virtual bool set_seed(const Alignment &seed) override;

    virtual bool update_seed_filter(node_index node,
                                    size_t query_start,
                                    const score_t *s_begin,
                                    const score_t *s_end);

    virtual bool filter_nodes(node_index node, size_t query_start, size_t query_end);
};

class DefaultColumnExtender : public SeedFilteringExtender {
  public:
    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual void set_graph(const DeBruijnGraph &graph) override {
        SeedFilteringExtender::set_graph(graph);
        graph_ = &graph;
    }

  protected:
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;
    std::string_view query_;

    // During extension, a tree is constructed from the graph starting at the
    // seed, then the query is aligned against this tree.
    // Each Column object represents the alignment of a substring of the query
    // against a node in the tree.
    // The horizontal concatenation (hstack) of all of the columns along a path
    // in this tree is analogous to a Needleman-Wunsch dynamic programming score matrix.
    using Column = std::tuple<AlignedVector<score_t> /* S (best score) */,
                              AlignedVector<score_t> /* E (best score after insert) */,
                              AlignedVector<score_t> /* F (best score after delete) */,
                              node_index /* node */,
                              size_t /* parent index in table */,
                              char /* last char of node */,
                              ssize_t /* offset (distance from start of the first node) */,
                              ssize_t /* absolute index of maximal value*/,
                              ssize_t /* trim (starting absolute index of array) */>;
    // e.g., the maximal value is located at S[std::get<7>(col) - std::get<8>(col)]
    std::vector<Column> table;

    tsl::hopscotch_set<size_t> prev_starts;

    virtual std::vector<Alignment> extend(score_t min_path_score) override;

    // backtracking helpers
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &extensions) const {
        return extensions.size() >= config_.num_alternative_paths;
    }

    virtual bool terminate_backtrack() const { return false; }

    virtual bool skip_backtrack_start(size_t table_i) {
        return !prev_starts.emplace(table_i).second;
    }

    virtual void call_alignments(score_t cur_cell_score,
                                 score_t end_score,
                                 score_t min_path_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> & /* trace */,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 const std::function<void(Alignment&&)> &callback) {
        assert(path.size());
        assert(ops.size());

        if (cur_cell_score == 0 && ops.back().first == Cigar::MATCH
                && end_score >= min_path_score) {
            callback(construct_alignment(ops, clipping, window, path, match,
                                         end_score, offset));
        }
    }

    virtual void init_backtrack() {}

    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char)> &callback);

    Alignment construct_alignment(Cigar cigar,
                                  size_t clipping,
                                  std::string_view window,
                                  std::vector<node_index> final_path,
                                  std::string match,
                                  score_t score,
                                  size_t offset) const;

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<score_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;

    // backtrack through the DP table to reconstruct alignments
    std::vector<Alignment> backtrack(score_t min_path_score,
                                        std::string_view window);
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
