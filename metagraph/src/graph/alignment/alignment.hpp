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
#include "common/vector.hpp"


namespace mtg {
namespace graph {
namespace align {

class AnnotationBuffer;

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Seed {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef SmallVector<int64_t> Tuple;
    typedef size_t Columns;
    typedef Vector<Tuple> CoordinateSet;

    Seed() : orientation_(false), offset_(0), clipping_(0), end_clipping_(0) {}
    Seed(std::string_view query_view,
         std::vector<node_index>&& nodes,
         bool orientation,
         size_t offset,
         size_t clipping,
         size_t end_clipping)
          : query_view_(query_view), nodes_(std::move(nodes)), orientation_(orientation),
            offset_(offset), clipping_(clipping), end_clipping_(end_clipping) {}

    std::string_view get_query_view() const { return query_view_; }

    bool empty() const { return nodes_.empty(); }

    // The matched path in the graph
    const std::vector<node_index>& get_nodes() const { return nodes_; }

    // get the spelling of the matched path
    std::string_view get_sequence() const { return query_view_; }

    // The number of characters discarded from the prefix of the first node in the path.
    size_t get_offset() const { return offset_; }

    // The number of nodes in the matched path.
    // This is equal to get_sequence().size() - graph.get_k() + 1 + get_offset()
    size_t size() const { return nodes_.size(); }

    // The orientation of the original query sequence. Return true iff the reverse
    // complement is matched to the path.
    bool get_orientation() const { return orientation_; }

    Cigar::LengthType get_clipping() const { return clipping_; }
    Cigar::LengthType get_end_clipping() const { return end_clipping_; }

    // Add more nodes to the Seed. This assumes that the original query has |next|
    // extra characters and that the updated path still spells the query view.
    void expand(const std::vector<node_index> &next) {
        query_view_ = std::string_view(query_view_.data(), query_view_.size() + next.size());
        end_clipping_ -= next.size();
        nodes_.insert(nodes_.end(), next.begin(), next.end());
    }

    AnnotationBuffer *label_encoder = nullptr;
    bool has_annotation() const { return label_encoder; }

    Columns label_columns = 0;

    const Vector<Column>& get_columns() const;
    void set_columns(Vector<Column>&& columns);

    // for each column in |label_columns|, store a vector of coordinates for the
    // alignment's first nucleotide
    // (i.e., the first node's coordinate + the alignment offset)
    CoordinateSet label_coordinates;

  private:
    std::string_view query_view_;
    std::vector<node_index> nodes_;
    bool orientation_;
    size_t offset_;
    Cigar::LengthType clipping_;
    Cigar::LengthType end_clipping_;
};

template <class T>
struct Dereferencer {
    template <class It>
    inline T& operator()(It it) const { return *it; }
};

template <class T>
struct FirstDereferencer {
    template <class It>
    inline typename T::first_type& operator()(It it) const { return it->first; }
};

template <class It, class Transformer = Dereferencer<typename std::iterator_traits<It>::value_type>>
inline size_t get_num_char_matches_in_seeds(It begin, It end) {
    Transformer transformer;
    size_t num_matching = 0;
    size_t last_q_end = 0;
    for (auto it = begin; it != end; ++it) {
        if (transformer(it).empty())
            continue;

        size_t q_begin = transformer(it).get_clipping();
        size_t q_end = q_begin + transformer(it).get_query_view().size();
        if (q_end > last_q_end) {
            num_matching += q_end - q_begin;
            if (q_begin < last_q_end)
                num_matching -= last_q_end - q_begin;
        }

        if (size_t offset = transformer(it).get_offset()) {
            size_t clipping = transformer(it).get_clipping();
            for (++it; it != end && transformer(it).get_offset() == offset
                                    && transformer(it).get_clipping() == clipping; ++it) {}
            --it;
        }

        last_q_end = q_end;
    }
    return num_matching;
}

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Alignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef annot::binmat::BinaryMatrix::Column Column;
    typedef SmallVector<int64_t> Tuple;
    typedef size_t Columns;
    typedef Vector<Tuple> CoordinateSet;
    typedef DBGAlignerConfig::score_t score_t;
    static const score_t ninf = DBGAlignerConfig::ninf;

    Alignment(std::string_view query = {},
              std::vector<node_index>&& nodes = {},
              std::string&& sequence = "",
              score_t score = 0,
              Cigar&& cigar = {},
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_view_(query), nodes_(std::move(nodes)), orientation_(orientation),
            offset_(offset), sequence_(std::move(sequence)), score_(score),
            cigar_(Cigar::CLIPPED, clipping) { cigar_.append(std::move(cigar)); }

    Alignment(const Seed &seed, const DBGAlignerConfig &config)
          : label_encoder(seed.label_encoder), label_columns(seed.label_columns),
            label_coordinates(seed.label_coordinates), query_view_(seed.get_query_view()),
            nodes_(std::vector<node_index>(seed.get_nodes())),
            orientation_(seed.get_orientation()), offset_(seed.get_offset()),
            sequence_(query_view_),
            score_(config.match_score(query_view_) + (!seed.get_clipping() ? config.left_end_bonus : 0)
                     + (!seed.get_end_clipping() ? config.right_end_bonus : 0)),
            cigar_(Cigar::CLIPPED, seed.get_clipping()) {
        cigar_.append(Cigar::MATCH, query_view_.size());
        cigar_.append(Cigar::CLIPPED, seed.get_end_clipping());
    }

    std::string_view get_query_view() const { return query_view_; }

    bool empty() const { return nodes_.empty(); }

    // The matched path in the graph
    const std::vector<node_index>& get_nodes() const { return nodes_; }

    // get the spelling of the matched path
    std::string_view get_sequence() const { return sequence_; }

    // The number of characters discarded from the prefix of the first node in the path.
    size_t get_offset() const { return offset_; }

    // The number of nodes in the matched path.
    // This is equal to get_sequence().size() - graph.get_k() + 1 + get_offset()
    size_t size() const { return nodes_.size(); }

    // The orientation of the original query sequence. Return true iff the reverse
    // complement is matched to the path.
    bool get_orientation() const { return orientation_; }

    // Append |next| to the end of the current alignment. In this process, alignment
    // labels are intersected. If coordinates are present, then the append is only
    // successful if at least one coordinate of |next| immediately proceeds the
    // one of the coordinates in this. If this operation is unsuccessful, then
    // *this == {} afterwards.
    // Returns true if the label or coordinate set of this changed.
    bool append(Alignment&& next, score_t label_change_score = DBGAlignerConfig::ninf);

    bool splice(Alignment&& other);

    score_t get_score() const { return score_; }

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

    size_t trim_offset(size_t num_nodes = std::numeric_limits<size_t>::max());
    void extend_offset(std::vector<node_index>&& path,
                       std::vector<size_t>&& columns = {},
                       std::vector<score_t>&& scores = {});

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

    void splice_with_unknown(Alignment&& other, size_t num_unknown, size_t node_overlap, const DBGAlignerConfig &config);

    void reverse_complement(const DeBruijnGraph &graph, std::string_view query_rev_comp);

    const Cigar& get_cigar() const { return cigar_; }
    Cigar& get_cigar() { return cigar_; }
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

    std::string to_gaf(const std::string &name = {}) const;

    // writes to |query_str| the string which will be referenced in this object
    void load_from_json(const Json::Value &alignment,
                        const DeBruijnGraph &graph,
                        std::string *query_str);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

    AnnotationBuffer *label_encoder = nullptr;
    bool has_annotation() const { return label_encoder; }

    Columns label_columns = 0;

    // for each column in |label_columns|, store a vector of coordinates for the
    // alignment's first nucleotide
    // (i.e., the first node's coordinate + the alignment offset)
    CoordinateSet label_coordinates;

    static bool coordinates_less(const Alignment &a, const Alignment &b);

    std::vector<Columns> label_column_diffs;
    std::vector<score_t> extra_scores;
    score_t extra_score = 0;

    std::string format_coords() const;
    std::string format_annotations() const;

    void set_columns(Vector<Column>&& columns);
    const Vector<Column>& get_columns(size_t path_i = 0) const;
    Vector<Column> get_column_union() const;
    void merge_annotations(const Alignment &other);

    std::vector<std::string> get_decoded_labels(size_t path_i) const;

  private:
    std::string_view query_view_;
    std::vector<node_index> nodes_;
    bool orientation_;
    size_t offset_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
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
    AlignmentResults(std::string_view query = {});

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

        if (a.label_coordinates.size()) {
            format_to(ctx.out(), "\t{}", a.format_coords());
        } else if (a.has_annotation()) {
            format_to(ctx.out(), "\t{}", a.format_annotations());
        }

        return ctx.out();
    }
};
} // namespace fmt

#endif  // __ALIGNER_ALIGNMENT_HPP__
