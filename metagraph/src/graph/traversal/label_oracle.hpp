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
#include "graph/annotation_key_mapper.hpp"
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
    struct Counters {
        uint64_t keys_mapped = 0;          // oriented nodes mapped to annotation keys
        uint64_t rows_requested = 0;       // annotation rows requested (incl. cache hits)
        uint64_t cache_hits = 0;
        uint64_t rows_fetched = 0;         // rows reconstructed with get_rows
        uint64_t tuple_rows_fetched = 0;   // rows reconstructed with get_row_tuples
        uint64_t coords_mapped = 0;        // coordinates mapped to (sequence, local coord)
        double fetch_seconds = 0;
    };

    // |coord_to_header_override| replaces the mapping loaded with the index (tests).
    explicit LabelOracle(const AnnotatedDBG &anno_graph,
                         const annot::CoordToHeader *coord_to_header_override = nullptr);
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

    // Annotation keys of all k-mers of |sequence| (npos for k-mers not in the graph).
    std::vector<node_index> keys_of_sequence(std::string_view sequence) const;
    // Annotation key of the oriented |node| spelling |kmer|.
    node_index key_of(node_index node, std::string_view kmer) const;

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
    // The first and last column coordinates of sequence |seq_id|. Requires a CoordToHeader.
    std::pair<uint64_t, uint64_t> sequence_coords(Column column, uint64_t seq_id) const;
    // A column name or a header. Throws std::invalid_argument naming the label if neither
    // exists, or if it is ambiguous: both a column and a header, or the header of several
    // records.
    LabelRef resolve_label(const std::string &name) const;

    // Batched row access (rows = annotation key - 1). Updates the counters. On a row-diff
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

    Counters& counters() const { return counters_; }

  private:
    const AnnotatedDBG &anno_graph_;
    std::shared_ptr<CanonicalDBG> local_canonical_;
    const DeBruijnGraph *graph_ = nullptr;
    std::unique_ptr<AnnotationKeyMapper> key_mapper_;
    Regime regime_ = Regime::BASIC;

    const annot::matrix::BinaryMatrix *matrix_ = nullptr;
    const annot::matrix::MultiIntMatrix *tuples_ = nullptr;
    const annot::matrix::IRowDiff *row_diff_ = nullptr;
    const annot::CoordToHeader *coord_to_header_ = nullptr;

    mutable Counters counters_;
};


/**
 * Batched membership of a fixed set of permitted labels at annotation keys. Column
 * labels are read from whole rows; header labels (one sequence of a column) from tuple
 * rows, a coordinate counting for a label when it falls into that label's sequence.
 * Results are cached per key, up to |max_cache_size| keys, then the cache is cleared. A
 * batch with more keys than that is answered without caching it. On a row-diff
 * annotation, the decoded rows are kept too (RowDiffCache), so that the rows of the next
 * step of a walk are decoded from them.
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

    const std::vector<LabelRef>& labels() const { return labels_; }
    // which accessor is used: "rows" or "tuples"
    const char* access_path() const { return tuples_ ? "tuples" : "rows"; }

    // Hits for every annotation key (npos yields no hits).
    std::vector<NodeHits> fetch(const std::vector<node_index> &keys);

  private:
    struct HeaderRange {
        uint64_t first;
        uint64_t last;
        LabelId label;
    };

    void read_hits(const std::vector<node_index> &keys,
                   tsl::hopscotch_map<node_index, NodeHits> *out);
    void hits_from_row(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                       NodeHits *hits) const;
    void hits_from_tuples(const annot::matrix::MultiIntMatrix::RowTuples &row,
                          NodeHits *hits) const;

    const LabelOracle &oracle_;
    std::vector<LabelRef> labels_;
    bool tuples_ = false;
    size_t max_cache_size_;

    // column -> label id for COLUMN labels
    tsl::hopscotch_map<Column, LabelId> column_labels_;
    // per column holding HEADER labels, the coordinate ranges of the requested sequences,
    // ascending: a coordinate is looked up among them by a binary search instead of being
    // mapped to its sequence
    tsl::hopscotch_map<Column, std::vector<HeaderRange>> header_ranges_;

    tsl::hopscotch_map<node_index, NodeHits> cache_;
    annot::matrix::RowDiffCache<annot::matrix::BinaryMatrix::SetBitPositions> row_cache_;
    annot::matrix::RowDiffCache<annot::matrix::MultiIntMatrix::RowTuples> tuple_cache_;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_LABEL_ORACLE_HPP__
