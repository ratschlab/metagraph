#ifndef __DBG_EXTENDER_METHODS_HPP__
#define __DBG_EXTENDER_METHODS_HPP__

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "alignment.hpp"
#include "common/aligned_vector.hpp"


namespace mtg {
namespace graph {
namespace align {

class IDBGAligner;

class SeedFilteringExtender {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef Alignment::score_t score_t;

    SeedFilteringExtender(const DeBruijnGraph *graph,
                          const DBGAlignerConfig &config,
                          std::string_view query)
          : graph_(graph), config_(config), query_size_(query.size()) {
        assert(config_.check_config_scores());
    }

    virtual ~SeedFilteringExtender() {}

    /**
     * If force_fixed_seed is true, then all alignments must have the seed as a
     * prefix. Otherwise, only the first node of the seed is used as a starting node.
     * If target_length and target_node are set, then terminate the extension
     * once target_length nucleotides have been aligned to and only backtrack
     * from target_node
     */
    std::vector<Alignment> get_extensions(const Alignment &seed,
                                          score_t min_path_score,
                                          bool force_fixed_seed,
                                          size_t target_length = 0,
                                          node_index target_node = DeBruijnGraph::npos,
                                          bool trim_offset_after_extend = true,
                                          size_t trim_query_suffix = 0,
                                          score_t added_xdrop = 0) {
        if (!set_seed(seed))
            return {};

        return extend(min_path_score, force_fixed_seed, target_length, target_node,
                      trim_offset_after_extend, trim_query_suffix, added_xdrop);
    }

    virtual std::string_view get_query() const = 0;

    // report alignment extension statistics
    virtual size_t num_explored_nodes() const = 0;

    // return true if the nodes in this seed have been traversed previously with
    // better or equal scores
    virtual bool check_seed(const Alignment &seed) const;

    virtual bool filter_nodes(node_index node, size_t query_start, size_t query_end);

    virtual void clear_conv_checker() = 0;

    void set_graph(const DeBruijnGraph &graph) { graph_ = &graph; }

    void extend_seed_end(const Alignment &seed,
                         const std::function<void(Alignment&&)> &callback,
                         bool force_fixed_seed,
                         score_t min_path_score = DBGAlignerConfig::ninf);

    void rc_extend_rc(const Alignment &seed,
                      const std::function<void(Alignment&&)> &callback,
                      bool force_fixed_seed,
                      score_t min_path_score = DBGAlignerConfig::ninf);

  protected:
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig &config_;
    const Alignment *seed_ = nullptr;
    size_t query_size_;

    typedef std::pair<size_t, AlignedVector<score_t>> ScoreVec;
    tsl::hopscotch_map<node_index, ScoreVec> conv_checker_;

    size_t explored_nodes_previous_ = 0;

    /**
     * Helper function for running the extension.
     * If force_fixed_seed is true, then all alignments must have the seed as a
     * prefix. Otherwise, only the first node of the seed is used as a starting node.
     * If target_length and target_node are set, then terminate the extension
     * once target_length nucleotides have been aligned to and only backtrack
     * from target_node
     */
    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed,
                                          size_t target_length = 0,
                                          node_index target_node = DeBruijnGraph::npos,
                                          bool trim_offset_after_extend = true,
                                          size_t trim_query_suffix = 0,
                                          score_t added_xdrop = 0) = 0;

    virtual bool set_seed(const Alignment &seed);

    virtual score_t update_seed_filter(node_index node,
                                       size_t query_start,
                                       const score_t *s_begin,
                                       const score_t *s_end);
};

class DefaultColumnExtender : public SeedFilteringExtender {
  public:
    // to ensure that SIMD operations on arrays don't read out of bounds
    static const size_t kPadding = 5;

    DefaultColumnExtender(const DeBruijnGraph &graph,
                          const DBGAlignerConfig &config,
                          std::string_view query);

    DefaultColumnExtender(const IDBGAligner &aligner, std::string_view query);

    virtual ~DefaultColumnExtender() {}

    size_t num_extensions() const { return num_extensions_; }
    std::string_view get_query() const { return query_; }

    /**
     * During extension, a tree is constructed from the graph starting at the
     * seed, then the query is aligned against this tree.
     * Each DPTColumn object represents the alignment of a substring of the query
     * against a node in the tree.
     * The horizontal concatenation (hstack) of all of the columns along a path
     * in this tree is analogous to a Needleman-Wunsch dynamic programming score matrix.
     */
    struct DPTColumn {
        AlignedVector<score_t> S; // best score
        AlignedVector<score_t> E; // best score after insert
        AlignedVector<score_t> F; // best score after delete
        node_index node; // graph node represented by the column
        size_t parent_i; // index of the parent column
        char c; // the last character of the node's k-mer
        ssize_t offset; // distance from the starting character of the first node
        ssize_t max_pos; // absolute index of the maximal value
        ssize_t trim; // first index allocated by the vectors
                      // i.e., the maximal value is located at S[max_pos - trim]
        size_t xdrop_cutoff_i; // corresponding index in the xdrop_cutoff vector
        score_t score; // added score when traversing to this node (typically a negative penalty)

        // allocate and initialize with padding to ensure that SIMD operations don't
        // read/write out of bounds
        template <typename... RestArgs>
        static DPTColumn create(size_t size, RestArgs&&... args);
    };

    void clear_conv_checker() {
        explored_nodes_previous_ += table.size();
        table.clear();
        conv_checker_.clear();
    }

    size_t num_explored_nodes() const {
        return explored_nodes_previous_ + table.size();
    }

  protected:
    std::string_view query_;
    std::vector<DPTColumn> table;
    size_t table_size_bytes_;

    tsl::hopscotch_set<size_t> prev_starts;
    std::vector<std::pair<size_t, score_t>> xdrop_cutoffs_;
    size_t num_extensions_ = 0;

    // If force_fixed_seed is true, then all alignments must have the seed as a
    // prefix. Otherwise, only the first node of the seed is used as a starting node.
    // If target_length and target_node are set, then terminate the extension
    // once target_length nucleotides have been aligned to and only backtrack
    // from target_node
    virtual std::vector<Alignment> extend(score_t min_path_score,
                                          bool force_fixed_seed,
                                          size_t target_length = 0,
                                          node_index target_node = DeBruijnGraph::npos,
                                          bool trim_offset_after_extend = true,
                                          size_t trim_query_suffix = 0,
                                          score_t added_xdrop = 0) override;

    virtual void call_outgoing(node_index node,
                               size_t max_prefetch_distance,
                               const std::function<void(node_index, char, score_t)> &callback,
                               size_t table_i,
                               bool force_fixed_seed = false);

    // Remove an element from the DP table. This method can break the table if
    // there exist j,j' s.t. i<=j<j' and j is the parent of j'
    virtual void pop(size_t i) { table.erase(table.begin() + i); }

    // Backtrack through the DP table to reconstruct alignments. If target_node
    // is set, then only backtrack from target_node.
    virtual std::vector<Alignment> backtrack(score_t min_path_score,
                                             std::string_view window,
                                             score_t right_end_bonus,
                                             const std::vector<size_t> &tips,
                                             node_index target_node = DeBruijnGraph::npos);

    /**
     * Backtracking helpers
     */

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
    virtual void call_alignments(score_t end_score,
                                 const std::vector<node_index> &path,
                                 const std::vector<size_t> & /* trace */,
                                 const std::vector<score_t> &score_trace,
                                 const Cigar &ops,
                                 size_t clipping,
                                 size_t offset,
                                 std::string_view window,
                                 const std::string &match,
                                 score_t extra_score,
                                 const std::function<void(Alignment&&)> &callback) {
        callback(construct_alignment(ops, clipping, window, path, match, end_score,
                                     offset, score_trace, extra_score));
    }

    Alignment construct_alignment(Cigar cigar,
                                  size_t clipping,
                                  std::string_view window,
                                  std::vector<node_index> final_path,
                                  std::string match,
                                  score_t score,
                                  size_t offset,
                                  const std::vector<score_t> &score_trace,
                                  score_t extra_score) const;

  private:
    // compute perfect match scores for all suffixes used for branch and bound checks
    std::vector<score_t> partial_sums_;

    // a quick lookup table of char pair match/mismatch scores for the current query
    std::vector<AlignedVector<score_t>> profile_score_;
    std::vector<AlignedVector<Cigar::Operator>> profile_op_;

    std::vector<score_t> scores_reached_;
    score_t min_cell_score_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __DBG_EXTENDER_METHODS_HPP__
