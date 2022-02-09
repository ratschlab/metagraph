#ifndef __ALIGNER_ALIGNMENT_HPP__
#define __ALIGNER_ALIGNMENT_HPP__

#include <cassert>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <json/json.h>

#include "aligner_cigar.hpp"
#include "aligner_config.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "common/vector.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace graph {
namespace align {

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Alignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef DBGAlignerConfig::score_t score_t;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef annot::matrix::MultiIntMatrix::Tuple Tuple;
    typedef Vector<Column> LabelSet;
    typedef Vector<Tuple> CoordinateSet;

    Alignment() : score_(0), orientation_(false), offset_(0) {}

    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              Cigar&& cigar,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
            score_(score), cigar_(Cigar::CLIPPED, clipping), orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    // Used for constructing gapless Alignments
    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0);

    // Append other to the end of the current alignment. In this process, alignment
    // labels are intersected. If coordinates are present, then the append is only
    // successful if at least one coordinate of other immediately proceeds the
    // one of the coordinates in this. If this operation is unsuccessful, then
    // *this == {} afterwards.
    // Returns true if the label or coordinate set of this changed.
    bool append(Alignment&& other);

    size_t size() const { return nodes_.size(); }
    bool empty() const { return nodes_.empty(); }
    const std::vector<node_index>& get_nodes() const { return nodes_; }

    score_t get_score() const { return score_; }

    std::string_view get_query() const { return query_; }

    void extend_query_begin(const char *begin) {
        const char *full_query_begin = query_.data() - get_clipping();
        assert(full_query_begin >= begin);
        if (full_query_begin > begin)
            cigar_.extend_clipping(full_query_begin - begin);
    }

    void extend_query_end(const char *end) {
        const char *full_query_end = query_.data() + query_.size() + get_end_clipping();
        assert(full_query_end <= end);
        if (full_query_end < end)
            cigar_.append(Cigar::CLIPPED, end - full_query_end);
    }

    inline size_t trim_clipping() { return cigar_.trim_clipping(); }
    inline size_t trim_end_clipping() { return cigar_.trim_end_clipping(); }

    size_t trim_offset();

    size_t trim_query_prefix(size_t n,
                             const DeBruijnGraph &graph,
                             const DBGAlignerConfig &config,
                             bool trim_excess_deletions = true);

    size_t trim_query_suffix(size_t n,
                             const DeBruijnGraph &graph,
                             const DBGAlignerConfig &config,
                             bool trim_excess_deletions = true);

    size_t trim_reference_prefix(size_t n,
                                 const DeBruijnGraph &graph,
                                 const DBGAlignerConfig &config,
                                 bool trim_excess_insertions = true);

    size_t trim_reference_suffix(size_t n,
                                 const DeBruijnGraph &graph,
                                 const DBGAlignerConfig &config,
                                 bool trim_excess_insertions = true);

    // When chaining together two alignments, use this method to adapt the prefix
    // of this alignment so it can be appended to the first one.
    // a negative gap indicates an overlap
    void insert_gap_prefix(ssize_t gap_length, const DeBruijnGraph &graph, const DBGAlignerConfig &config);

    void reverse_complement(const DeBruijnGraph &graph, std::string_view query_rev_comp);

    const std::string& get_sequence() const { return sequence_; }
    const Cigar& get_cigar() const { return cigar_; }
    Cigar& get_cigar() { return cigar_; }
    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }
    Cigar::LengthType get_end_clipping() const { return cigar_.get_end_clipping(); }

    bool operator==(const Alignment &other) const {
        return orientation_ == other.orientation_
            && offset_ == other.offset_
            && score_ == other.score_
            && query_ == other.query_
            && sequence_ == other.sequence_
            && cigar_ == other.cigar_
            && nodes_ == other.nodes_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    Json::Value to_json(std::string_view query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        std::string_view name = {},
                        std::string_view label = {}) const;

    // returns a shared_ptr of the query string which is referenced in this object
    std::shared_ptr<const std::string> load_from_json(const Json::Value &alignment,
                                                      const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

    static bool coordinates_less(const Alignment &a, const Alignment &b);

    LabelSet label_columns;
    score_t extra_penalty = 0;

    // for each column in |label_columns|, store a vector of coordinate ranges
    CoordinateSet label_coordinates;

    const annot::LabelEncoder<> *label_encoder = nullptr;

    std::string format_coords() const;

  private:
    Json::Value path_json(size_t node_size, std::string_view label = {}) const;

    // TODO: rename to query_view_
    std::string_view query_;
    std::vector<node_index> nodes_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
};

std::ostream& operator<<(std::ostream& out, const Alignment &alignment);

struct LocalAlignmentLess {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is less, or
        // 2) more of the query is covered, or
        // 3) if it is in the reverse orientation, or
        // 4) if the starting point is later in the query
        return std::make_tuple(b.get_score(), a.get_query().size(),
                               a.get_orientation(), a.get_clipping())
            > std::make_tuple(a.get_score(), b.get_query().size(),
                              b.get_orientation(), b.get_clipping());
    }
};

struct LocalAlignmentGreater {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is higher, or
        // 2) less of the query is covered, or
        // 3) if it is in the forward orientation, or
        // 4) if the starting point is earlier in the query
        return std::make_tuple(a.get_score(), b.get_query().size(),
                               b.get_orientation(), b.get_clipping())
            > std::make_tuple(b.get_score(), a.get_query().size(),
                              a.get_orientation(), a.get_clipping());
    }
};

// A container holding many alignments to a shared query sequence.
// Each alignment only holds a string_view to the query, so this class
// ensures that the query sequence is always accessible.
class AlignmentResults {
  public:
    explicit AlignmentResults(std::string_view query, bool is_reverse_complement = false);

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        assert(alignments_.back().get_query().data()
            >= get_query(alignments_.back().get_orientation()).c_str());
        assert(alignments_.back().get_query().data() + alignments_.back().get_query().size()
            <= get_query(alignments_.back().get_orientation()).c_str()
                + get_query(alignments_.back().get_orientation()).size());
    }

    const std::string& get_query(bool reverse_complement = false) const {
        return !reverse_complement ? query_ : query_rc_;
    }

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }
    const Alignment& operator[](size_t i) const { return alignments_[i]; }

    auto begin() const { return alignments_.begin(); }
    auto end() const { return alignments_.end(); }

  private:
    std::string query_;
    std::string query_rc_;
    std::vector<Alignment> alignments_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif  // __ALIGNER_ALIGNMENT_HPP__
