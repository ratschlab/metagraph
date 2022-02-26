#ifndef __ALIGNER_ALIGNMENT_HPP__
#define __ALIGNER_ALIGNMENT_HPP__

#include <cassert>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <json/json.h>
#include <spdlog/fmt/fmt.h>

#include "aligner_cigar.hpp"
#include "aligner_config.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "common/vector.hpp"


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
    typedef SmallVector<int64_t> Tuple;
    typedef Vector<Column> Columns;
    typedef Vector<Tuple> CoordinateSet;
    static const score_t ninf = DBGAlignerConfig::ninf;

    Alignment() : score_(0), orientation_(false), offset_(0) {}

    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              Cigar&& cigar,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_view_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
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

    // Append |next| to the end of the current alignment. In this process, alignment
    // labels are intersected. If coordinates are present, then the append is only
    // successful if at least one coordinate of |next| immediately proceeds the
    // one of the coordinates in this. If this operation is unsuccessful, then
    // *this == {} afterwards.
    // Returns true if the label or coordinate set of this changed.
    bool append(Alignment&& next);

    bool splice(Alignment&& other) {
        if (empty()) {
            std::swap(*this, other);
            return label_columns.size();
        }

        trim_end_clipping();
        other.trim_clipping();
        return append(std::move(other));
    }

    size_t size() const { return nodes_.size(); }
    bool empty() const { return nodes_.empty(); }
    const std::vector<node_index>& get_nodes() const { return nodes_; }

    score_t get_score() const { return score_; }

    std::string_view get_query_view() const { return query_view_; }

    void extend_query_begin(const char *begin) {
        const char *full_query_begin = query_view_.data() - get_clipping();
        assert(full_query_begin >= begin);
        if (full_query_begin > begin)
            cigar_.extend_clipping(full_query_begin - begin);
    }

    void extend_query_end(const char *end) {
        const char *full_query_end
            = query_view_.data() + query_view_.size() + get_end_clipping();
        assert(full_query_end <= end);
        if (full_query_end < end)
            cigar_.append(Cigar::CLIPPED, end - full_query_end);
    }

    inline size_t trim_clipping() { return cigar_.trim_clipping(); }
    inline size_t trim_end_clipping() { return cigar_.trim_end_clipping(); }

    size_t trim_offset();

    size_t trim_query_prefix(size_t n,
                             size_t node_overlap,
                             const DBGAlignerConfig &config,
                             bool trim_excess_deletions = true);

    size_t trim_query_suffix(size_t n,
                             const DBGAlignerConfig &config,
                             bool trim_excess_deletions = true);

    size_t trim_reference_prefix(size_t n,
                                 size_t node_overlap,
                                 const DBGAlignerConfig &config,
                                 bool trim_excess_insertions = true);

    size_t trim_reference_suffix(size_t n,
                                 const DBGAlignerConfig &config,
                                 bool trim_excess_insertions = true);

    // When chaining together two alignments, use this method to adapt the prefix
    // of this alignment so it can be appended to the first one.
    // a negative gap indicates an overlap
    void insert_gap_prefix(ssize_t gap_length, size_t node_overlap, const DBGAlignerConfig &config);

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
            && query_view_ == other.query_view_
            && sequence_ == other.sequence_
            && cigar_ == other.cigar_
            && nodes_ == other.nodes_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    Json::Value to_json(size_t node_size,
                        bool is_secondary = false,
                        const std::string &name = {},
                        const std::string &label = {}) const;

    // writes to |query_str| the string which will be referenced in this object
    void load_from_json(const Json::Value &alignment,
                        const DeBruijnGraph &graph,
                        std::string *query_str);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

    static bool coordinates_less(const Alignment &a, const Alignment &b);

    Columns label_columns;
    score_t extra_penalty = 0;

    // for each column in |label_columns|, store a vector of coordinate ranges
    CoordinateSet label_coordinates;

    const annot::LabelEncoder<> *label_encoder = nullptr;

    std::string format_coords() const;

  private:
    std::string_view query_view_;
    std::vector<node_index> nodes_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
};

inline std::ostream& operator<<(std::ostream &out, const Alignment &a) {
    return out << fmt::format("{}", a);
}

struct LocalAlignmentLess {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is less, or
        // 2) more of the query is covered, or
        // 3) if it is in the reverse orientation, or
        // 4) if the starting point is later in the query
        return std::make_tuple(b.get_score(), a.get_query_view().size(),
                               a.get_orientation(), a.get_clipping())
            > std::make_tuple(a.get_score(), b.get_query_view().size(),
                              b.get_orientation(), b.get_clipping());
    }
};

struct LocalAlignmentGreater {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is higher, or
        // 2) less of the query is covered, or
        // 3) if it is in the forward orientation, or
        // 4) if the starting point is earlier in the query
        return std::make_tuple(a.get_score(), b.get_query_view().size(),
                               b.get_orientation(), b.get_clipping())
            > std::make_tuple(b.get_score(), a.get_query_view().size(),
                              a.get_orientation(), a.get_clipping());
    }
};

// A container holding many alignments to a shared query sequence.
// Each alignment only holds a string_view to the query, so this class
// ensures that the query sequence is always accessible.
class AlignmentResults {
  public:
    AlignmentResults(std::string_view query = {}, bool is_reverse_complement = false);

    // Copy constructors are disabled to ensure that the string_view pointers
    // in the Alignment objects stay valid
    AlignmentResults(const AlignmentResults&) = delete;
    AlignmentResults& operator=(const AlignmentResults&) = delete;
    AlignmentResults(AlignmentResults&&) = default;
    AlignmentResults& operator=(AlignmentResults&&) = default;

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        assert(alignments_.back().get_query_view().data()
            >= get_query(alignments_.back().get_orientation()).c_str());
        assert(alignments_.back().get_query_view().data() + alignments_.back().get_query_view().size()
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

std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path,
                       size_t offset = 0);

struct CoordIntersection {
    CoordIntersection(int64_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + offset_ < *b_begin) {
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                *out = *a_begin;
                ++a_begin;
                ++b_begin;
                ++out;
            }
        }
    }

    int64_t offset_;
};

struct CoordDifference {
    CoordDifference(int64_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin + offset_ < *b_begin) {
                *out = *a_begin;
                ++a_begin;
                ++out;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                ++a_begin;
                ++b_begin;
            }
        }
    }

    int64_t offset_;
};

struct CoordUnion {
    CoordUnion(int64_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end || b_begin != b_end) {
            if (b_begin == b_end) {
                *out = *a_begin;
                ++a_begin;
                ++out;
            } else if (a_begin == a_end || *a_begin + offset_ > *b_begin) {
                *out = *b_begin - offset_;
                ++b_begin;
                ++out;
            } else {
                if (*a_begin + offset_ == *b_begin)
                    ++b_begin;

                *out = *a_begin;
                ++a_begin;
                ++out;
            }
        }
    }

    int64_t offset_;
};

} // namespace align
} // namespace graph
} // namespace mtg


namespace fmt {
template <> struct formatter<mtg::graph::align::Alignment> {
    // Parses format specifications of the form ['f' | 'e'].
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        // we have only one format, so nothing to parse
        return ctx.end();
    }

    template <typename FormatContext>
    auto format(const mtg::graph::align::Alignment &a, FormatContext &ctx) -> decltype(ctx.out()) {
        format_to(ctx.out(), "{}\t{}\t{}\t{}\t{}\t{}",
                  a.get_orientation() ? "-" : "+",
                  a.get_sequence(),
                  a.get_score(),
                  a.get_cigar().get_num_matches(),
                  a.get_cigar().to_string(),
                  a.get_offset());

        const auto &label_columns = a.label_columns;
        const auto &label_coordinates = a.label_coordinates;

        if (label_coordinates.size()) {
            format_to(ctx.out(), "\t{}", a.format_coords());
        } else if (label_columns.size()) {
            assert(a.label_encoder);

            std::vector<std::string> decoded_labels;
            decoded_labels.reserve(label_columns.size());

            for (size_t i = 0; i < label_columns.size(); ++i) {
                decoded_labels.emplace_back(a.label_encoder->decode(label_columns[i]));
            }

            format_to(ctx.out(), "\t{}", fmt::join(decoded_labels, ";"));
        }

        return ctx.out();
    }
};
} // namespace fmt

#endif  // __ALIGNER_ALIGNMENT_HPP__
