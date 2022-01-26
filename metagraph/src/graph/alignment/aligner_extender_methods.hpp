#ifndef __DBG_EXTENDER_METHODS_HPP__
#define __DBG_EXTENDER_METHODS_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "aligner_alignment.hpp"
#include "common/aligned_vector.hpp"
#include "common/hashers/hash.hpp"


namespace mtg {
namespace graph {
namespace align {

class IExtender {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    virtual ~IExtender() {}

    // If force_fixed_seed is true, then all alignments must have the seed as a
    // prefix. Otherwise, only the first node of the seed is used as a starting node.
    std::vector<Alignment>
    get_extensions(const Alignment &seed, score_t min_path_score, bool force_fixed_seed) {
        if (set_seed(seed))
            return extend(min_path_score, force_fixed_seed);

        return {};
    }

    virtual void set_graph(const DeBruijnGraph &graph) = 0;

    // report alignment extension statistics
    virtual size_t num_explored_nodes() const = 0;
    virtual size_t num_extensions() const = 0;

    // return true if the nodes in this seed have been traversed previously with
    // better or equal scores
    virtual bool check_seed(const Alignment &seed) const = 0;

  protected:
    // set the seed from which to extend and return true if extension can proceed
    virtual bool set_seed(const Alignment &seed) = 0;

    // Helper function for running the extension.
    // If force_fixed_seed is true, then all alignments must have the seed as a
    // prefix. Otherwise, only the first node of the seed is used as a starting node.
    virtual std::vector<Alignment> extend(score_t min_path_score, bool force_fixed_seed) = 0;

    // returns whether the seed must be a prefix of an extension
    virtual bool fixed_seed() const { return true; }
};

class GreedyExtender : public IExtender {
  public:
    GreedyExtender(const DeBruijnGraph &graph,
                   const DBGAlignerConfig &config,
                   std::string_view query,
                   const std::vector<Alignment> &seeds = {});

    virtual ~GreedyExtender() {}

    virtual void set_graph(const DeBruijnGraph &graph) override final { graph_ = &graph; }

    virtual size_t num_explored_nodes() const override final { return num_explored_nodes_; }
    virtual size_t num_extensions() const override final { return num_extensions_; }

    virtual bool check_seed(const Alignment &seed) const override final;

    const tsl::hopscotch_set<node_index>& get_explored_nodes() const { return nodes_; }

  private:
    const DeBruijnGraph *graph_;
    DBGAlignerConfig config_;
    std::string_view query_;
    Alignment seed_;
    size_t num_explored_nodes_ = 0;
    size_t num_extensions_ = 0;
    tsl::hopscotch_set<node_index> nodes_;
    tsl::hopscotch_map<std::pair<node_index, size_t>, size_t, utils::Hash<std::pair<node_index, size_t>>> seeds_;

    virtual bool set_seed(const Alignment &seed) override final;

    virtual std::vector<Alignment> extend(score_t min_path_score, bool force_fixed_seed) override final;

    virtual bool terminate(uint64_t query_i, uint64_t ref_i) const;
};

class SeedFilteringExtender : public IExtender {
  public:
    SeedFilteringExtender(const DBGAlignerConfig &config, std::string_view query)
          : config_(config), query_size_(query.size()) {
        assert(config_.check_config_scores());
    }

    virtual ~SeedFilteringExtender() {}

    virtual size_t num_explored_nodes() const override {
        return explored_nodes_previous_ + conv_checker_.size();
    }

    virtual bool check_seed(const Alignment &seed) const override;

    virtual bool filter_nodes(node_index node, size_t query_start, size_t query_end);

    void clear_conv_checker() {
        explored_nodes_previous_ += conv_checker_.size();
        conv_checker_.clear();
    }

  protected:
    const DBGAlignerConfig &config_;
    const Alignment *seed_ = nullptr;
    size_t query_size_;

    typedef std::pair<size_t, AlignedVector<score_t>> ScoreVec;
    tsl::hopscotch_map<node_index, ScoreVec> conv_checker_;

    size_t explored_nodes_previous_ = 0;

    virtual bool set_seed(const Alignment &seed) override;

    virtual score_t update_seed_filter(node_index node,
                                       size_t query_start,
                                       const score_t *s_begin,
                                       const score_t *s_end);

    virtual const DeBruijnGraph& get_graph() const = 0;
};

class DefaultColumnExtender : public SeedFilteringExtender {
  public:
    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    virtual ~DefaultColumnExtender() {}

    virtual void set_graph(const DeBruijnGraph &graph) override {
        graph_ = &graph;
    }

    virtual size_t num_extensions() const override final { return num_extensions_; }

  protected:
    const DeBruijnGraph *graph_;
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
                              ssize_t /* trim (starting absolute index of array) */,
                              size_t /* xdrop_cutoffs_ index */,
                              size_t /* last fork index in table */,
                              score_t /* traversal penalty */>;
    // e.g., the maximal value is located at S[std::get<7>(col) - std::get<8>(col)]
    std::vector<Column> table;
    size_t table_size_bytes_;

    tsl::hopscotch_set<size_t> prev_starts;
    std::vector<std::pair<size_t, score_t>> xdrop_cutoffs_;
    size_t num_extensions_ = 0;

    virtual const DeBruijnGraph& get_graph() const override final {
        return *graph_;
    }

    // If force_fixed_seed is true, then all alignments must have the seed as a
    // prefix. Otherwise, only the first node of the seed is used as a starting node.
    virtual std::vector<Alignment> extend(score_t min_path_score, bool force_fixed_seed) override;

    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char, score_t)> &callback,
                               size_t table_i,
                               bool force_fixed_seed = false);

    // Remove an element from the DP table. This method can break the table if
    // there exist j,j' s.t. i<=j<j' and j is the parent of j'
    virtual void pop(size_t i) { table.erase(table.begin() + i); }

    // backtrack through the DP table to reconstruct alignments
    virtual std::vector<Alignment> backtrack(score_t min_path_score, std::string_view window);

    /**
     * Backtracking helpers
     */
    // halt the currently-running backtracking
    virtual bool terminate_backtrack() const { return false; }

    // stop considering new points from which to start backtracking
    virtual bool terminate_backtrack_start(const std::vector<Alignment> &extensions) const {
        return extensions.size() >= config_.num_alternative_paths;
    }

    // skip a backtracking start point
    virtual bool skip_backtrack_start(size_t table_i) {
        return !prev_starts.emplace(table_i).second;
    }

    // This method calls at most one alignment, but can be overridden by a child
    // class to call multiple alignments.
    virtual void call_alignments(score_t cur_cell_score,
                                 score_t end_score,
                                 score_t min_path_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> & /* trace */,
                                 size_t table_i,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 score_t extra_penalty,
                                 const std::function<void(Alignment&&)> &callback) {
        assert(path.size());
        assert(ops.size());

        if (end_score < min_path_score)
            return;

        if (clipping && cur_cell_score != 0)
            return;

        if (!clipping && cur_cell_score != std::get<0>(table[0])[0])
            return;

        if (!config_.allow_left_trim) {
            if (!table_i) {
                callback(construct_alignment(ops, clipping, window, path, match,
                                             end_score, offset, extra_penalty));
            }

            return;
        }

        callback(construct_alignment(ops, clipping, window, path, match,
                                     end_score, offset, extra_penalty));
    }

    Alignment construct_alignment(Cigar cigar,
                                  size_t clipping,
                                  std::string_view window,
                                  std::vector<node_index> final_path,
                                  std::string match,
                                  score_t score,
                                  size_t offset,
                                  score_t extra_penalty) const;

  private:
    // compute perfect match scores for all suffixes
    // used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    tsl::hopscotch_map<char, AlignedVector<score_t>> profile_score_;
    tsl::hopscotch_map<char, AlignedVector<Cigar::Operator>> profile_op_;

    std::vector<score_t> scores_reached_;
    score_t min_cell_score_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
