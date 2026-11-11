#ifndef __TRAVERSAL_LABEL_ORACLE_HPP__
#define __TRAVERSAL_LABEL_ORACLE_HPP__

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "traversal_types.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/binary_matrix/row_diff/row_diff_cache.hpp"


namespace mtg {

namespace annot {
class CoordToHeader;
namespace matrix {
class IRowDiff;
} // namespace matrix
} // namespace annot

namespace graph {

class AnnotatedDBG;
class DeBruijnGraph;
class CanonicalDBG;

namespace traversal {

/**
 * Per-request view of an annotated graph for traversal: the oriented graph to walk,
 * the mapping of its nodes to annotation rows, label names, and batched row access.
 *
 * In the PRIMARY regime the oracle walks a CanonicalDBG of its own over the loaded primary
 * graph (same node ids, its own caches), so concurrent requests never share wrapper state.
 * Not thread-safe: one oracle per request.
 */
class LabelOracle {
  public:
    explicit LabelOracle(const AnnotatedDBG &anno_graph);
    ~LabelOracle();

    const AnnotatedDBG& anno_graph() const { return anno_graph_; }
    // the graph used for traversal (oriented node ids)
    const DeBruijnGraph& graph() const { return *graph_; }
    Regime regime() const { return regime_; }
    size_t get_k() const;

    bool has_coordinates() const { return tuples_; }
    // the CoordToHeader index, or nullptr when there is none or the annotation has no
    // coordinates (header labels then do not exist)
    const annot::CoordToHeader* coord_to_header() const { return coord_to_header_; }
    uint64_t num_columns() const;

    // map_to_nodes(graph(), sequence) gives the annotation keys of a sequence's k-mers in
    // every regime; AnnotationKeyMapper gives the key of a single oriented node.

    // Label resolution
    std::optional<Column> find_column(const std::string &name) const;
    const std::string& column_name(Column column) const;
    // The first (column, sequence) holding a FASTA header (CoordToHeader::find_header)
    std::optional<LabelRef> find_header(const std::string &name) const;
    const std::string& header_name(Column column, uint64_t seq_id) const;
    // Whether |column| has sequence records, and so header labels (false without a
    // CoordToHeader).
    bool has_header_records(Column column) const;
    // (seq_id, local coordinate) of a column coordinate. Requires a CoordToHeader.
    std::pair<uint64_t, uint64_t> map_coord(Column column, uint64_t coord) const;
    // Whether |name| is both a column and a header, or the header of several records
    bool is_ambiguous(const std::string &name) const;
    // A column name or a header. Throws std::invalid_argument naming the label if neither
    // exists, or if it is ambiguous.
    LabelRef resolve_label(const std::string &name) const;

    // Batched row access (rows = annotation key - 1). On a row-diff
    // annotation, |cache| keeps the decoded rows across calls: for a walk, which reads the
    // rows of neighbouring nodes call after call.
    std::vector<annot::matrix::BinaryMatrix::SetBitPositions>
    get_rows(const std::vector<Row> &rows,
             annot::matrix::RowDiffCache<annot::matrix::BinaryMatrix::SetBitPositions> *cache
                 = nullptr) const;
    std::vector<annot::matrix::MultiIntMatrix::RowTuples>
    get_row_tuples(const std::vector<Row> &rows,
                   annot::matrix::RowDiffCache<annot::matrix::MultiIntMatrix::RowTuples> *cache
                       = nullptr) const;
    // whether the annotation is row-diff transformed, and so a cache helps
    bool is_row_diff() const { return row_diff_; }

  private:
    const AnnotatedDBG &anno_graph_;
    std::shared_ptr<CanonicalDBG> local_canonical_;
    const DeBruijnGraph *graph_ = nullptr;
    Regime regime_ = Regime::BASIC;

    const annot::matrix::BinaryMatrix *matrix_ = nullptr;
    const annot::matrix::MultiIntMatrix *tuples_ = nullptr;
    const annot::matrix::IRowDiff *row_diff_ = nullptr;
    const annot::CoordToHeader *coord_to_header_ = nullptr;
};


/**
 * Batched membership of a fixed set of permitted labels at annotation keys. Column
 * labels are read from whole rows; header labels (one sequence of a column) from tuple
 * rows, a coordinate counting for a label when it falls into that label's sequence. For
 * header labels the coordinates themselves are kept too (Hits::coords), so that a walk
 * can follow a record in coordinate order instead of by membership.
 * Results are cached per key, up to |max_cache_size| keys (and 16 header coordinates per
 * key on average), then the cache is cleared. A batch with more keys than that is
 * answered without caching it. On a row-diff annotation, the decoded rows are kept too
 * (RowDiffCache), so that the rows of the next step of a walk are decoded from them.
 */
class LabelQuery {
  public:
    // the permitted labels present at a key, ascending
    using NodeHits = std::vector<LabelId>;

    // Throws std::invalid_argument for an empty or duplicate label set, or header labels
    // on an annotation without coordinates.
    LabelQuery(const LabelOracle &oracle,
               std::vector<LabelRef> labels,
               size_t max_cache_size = 1'000'000);

    // A coordinate of a HEADER label's sequence record
    struct HeaderCoord {
        LabelId label;
        uint64_t coord;
        bool operator==(const HeaderCoord &o) const {
            return label == o.label && coord == o.coord;
        }
        bool operator<(const HeaderCoord &o) const {
            return label != o.label ? label < o.label : coord < o.coord;
        }
    };
    struct Hits {
        // the permitted labels present at the key, ascending
        NodeHits labels;
        // for the HEADER labels among them, the coordinates of the key that fall into the
        // label's record, ascending by (label, coordinate)
        std::vector<HeaderCoord> coords;
    };

    // Hits for every annotation key (npos yields no hits).
    std::vector<NodeHits> fetch(const std::vector<node_index> &keys);
    // The same with the header coordinates, without copying them: the pointers are valid
    // until the next fetch() or fetch_hits().
    std::vector<const Hits*> fetch_hits(const std::vector<node_index> &keys);

    // The first and last coordinate of a HEADER label's record ({ 0, 0 } for a COLUMN label)
    const std::pair<uint64_t, uint64_t>& record_coords(LabelId label) const {
        return record_coords_[label];
    }

  private:
    struct HeaderRange {
        uint64_t first;
        uint64_t last;
        LabelId label;
    };

    void read_hits(const std::vector<node_index> &keys,
                   tsl::hopscotch_map<node_index, Hits> *out);
    void hits_from_row(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                       Hits *hits) const;
    void hits_from_tuples(const annot::matrix::MultiIntMatrix::RowTuples &row,
                          Hits *hits) const;

    const LabelOracle &oracle_;
    std::vector<LabelRef> labels_;
    bool tuples_ = false;
    size_t max_cache_size_;
    // per label, the coordinate range of its record
    std::vector<std::pair<uint64_t, uint64_t>> record_coords_;

    // column -> label id for COLUMN labels
    tsl::hopscotch_map<Column, LabelId> column_labels_;
    // per column holding HEADER labels, the coordinate ranges of the requested sequences,
    // ascending: a coordinate is looked up among them by a binary search instead of being
    // mapped to its sequence
    tsl::hopscotch_map<Column, std::vector<HeaderRange>> header_ranges_;

    tsl::hopscotch_map<node_index, Hits> cache_;
    // the header coordinates held by |cache_|: a repetitive k-mer has many per label, so
    // the cache is bounded by them as well as by its keys
    size_t cached_coords_ = 0;
    // the hits of a batch too large to cache, kept until the next fetch
    tsl::hopscotch_map<node_index, Hits> uncached_;
    const Hits empty_ {};
    annot::matrix::RowDiffCache<annot::matrix::BinaryMatrix::SetBitPositions> row_cache_;
    annot::matrix::RowDiffCache<annot::matrix::MultiIntMatrix::RowTuples> tuple_cache_;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_LABEL_ORACLE_HPP__
