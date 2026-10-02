#ifndef __TRAVERSAL_LABEL_ORACLE_HPP__
#define __TRAVERSAL_LABEL_ORACLE_HPP__

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <tsl/hopscotch_map.h>

#include "traversal_types.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {

namespace annot {
class CoordToHeader;
}

namespace graph {

class AnnotatedDBG;
class DeBruijnGraph;
class CanonicalDBG;
class DBGSuccinct;
class NodeFirstCache;

namespace traversal {

/**
 * Per-request access to the graph and annotation of a loaded index (the
 * "TraversalContext" of the spec). Detects the graph regime, owns the per-request
 * graph handle (a CanonicalDBG clone with a NodeFirstCache for PRIMARY indexes), maps
 * oriented nodes to annotation rows and provides batched row / tuple access.
 *
 * All reads on the shared index are const; the caches created here are private to
 * this object, so one instance must not be used by several threads concurrently.
 */
class LabelOracle {
  public:
    enum class Access { AUTO, DIRECT, ROWS };

    struct Counters {
        uint64_t keys_mapped = 0;          // oriented nodes mapped to annotation keys
        uint64_t rows_requested = 0;       // annotation rows requested (incl. cache hits)
        uint64_t cache_hits = 0;
        uint64_t rows_fetched = 0;         // rows reconstructed with get_rows
        uint64_t tuple_rows_fetched = 0;   // rows reconstructed with get_row_tuples
        uint64_t direct_reads = 0;         // single-cell reads with GetEntrySupport::get
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

    // the CanonicalDBG wrapper used for traversal in the PRIMARY regime, else nullptr
    const CanonicalDBG* canonical() const { return canonical_; }
    // DBGSuccinct base graph (for BASIC/CANONICAL/PRIMARY over DBGSuccinct), nullptr otherwise
    const DBGSuccinct* dbg_succ() const { return dbg_succ_; }
    // cache for backward traversal on DBGSuccinct, nullptr if the base is not DBGSuccinct
    const NodeFirstCache* node_first_cache() const { return node_first_cache_.get(); }

    bool has_coordinates() const { return tuples_; }
    const annot::CoordToHeader* coord_to_header() const { return coord_to_header_; }
    bool supports_direct() const { return get_entry_; }
    uint64_t num_columns() const;
    uint64_t num_rows() const;

    // Annotation keys of all k-mers of |sequence| (npos for k-mers not in the graph).
    std::vector<node_index> keys_of_sequence(std::string_view sequence) const;
    // Annotation key of the oriented |node| spelling |kmer| (|kmer| may be empty in
    // the BASIC and PRIMARY regimes, where it is not needed).
    node_index key_of(node_index node, std::string_view kmer) const;
    // Annotation keys of the oriented nodes spelled by |window| (|window| spells the
    // path of |nodes|, i.e. |window|.size() == nodes.size() + k - 1).
    std::vector<node_index> keys_of_path(const std::vector<node_index> &nodes,
                                         std::string_view window) const;

    // Label resolution
    std::optional<Column> find_column(const std::string &name) const;
    const std::string& column_name(Column column) const;
    // Look up a FASTA header in the CoordToHeader index. The reverse index is built on
    // the first call and shared by all oracles of the same loaded index.
    std::optional<LabelRef> find_header(const std::string &name) const;
    const std::string& header_name(Column column, uint64_t seq_id) const;
    // (seq_id, local coordinate) of a column coordinate. Requires a CoordToHeader.
    std::pair<uint64_t, Coord> map_coord(Column column, Coord coord) const;
    // number of k-mers in the indexed sequence. Requires a CoordToHeader.
    uint64_t num_kmers_in_sequence(Column column, uint64_t seq_id) const;
    // Resolve a label name: a column name first, then a header. Throws
    // std::invalid_argument naming the label if neither exists.
    LabelRef resolve_label(const std::string &name) const;

    // Raw batched access (rows = annotation key - 1). Updates counters.
    std::vector<annot::matrix::BinaryMatrix::SetBitPositions>
    get_rows(const std::vector<Row> &rows) const;
    std::vector<annot::matrix::MultiIntMatrix::RowTuples>
    get_row_tuples(const std::vector<Row> &rows) const;
    bool get(Row row, Column column) const;

    Counters& counters() const { return counters_; }

  private:
    const AnnotatedDBG &anno_graph_;
    const DeBruijnGraph *graph_ = nullptr;
    std::shared_ptr<CanonicalDBG> local_canonical_;
    const CanonicalDBG *canonical_ = nullptr;
    const DBGSuccinct *dbg_succ_ = nullptr;
    std::shared_ptr<NodeFirstCache> node_first_cache_;
    Regime regime_ = Regime::BASIC;
    bool sshash_ = false;

    const annot::matrix::BinaryMatrix *matrix_ = nullptr;
    const annot::matrix::GetEntrySupport *get_entry_ = nullptr;
    const annot::matrix::MultiIntMatrix *tuples_ = nullptr;
    const annot::CoordToHeader *coord_to_header_ = nullptr;

    mutable Counters counters_;
};


/**
 * Batched membership of a fixed, small set of permitted labels at annotation keys,
 * optionally with the k-mer coordinates of each label (local coordinates in the
 * indexed sequence for HEADER labels, column coordinates for COLUMN labels).
 */
class LabelQuery {
  public:
    struct Hit {
        LabelId label;
        SmallVector<Coord> coords;  // sorted, only filled if |with_coords|
        bool operator==(const Hit &other) const {
            return label == other.label && coords == other.coords;
        }
    };
    // hits sorted by label id
    using NodeHits = std::vector<Hit>;

    // Throws std::invalid_argument if the requested access path is not available
    // for this annotation (no silent fallback).
    LabelQuery(const LabelOracle &oracle,
               std::vector<LabelRef> labels,
               bool with_coords,
               LabelOracle::Access access = LabelOracle::Access::AUTO,
               size_t max_cache_size = 1'000'000);

    const std::vector<LabelRef>& labels() const { return labels_; }
    bool with_coords() const { return with_coords_; }
    // which accessor is used: "direct", "rows" or "tuples"
    const char* access_path() const;

    // Hits for every annotation key (npos yields no hits).
    std::vector<NodeHits> fetch(const std::vector<node_index> &keys);
    // Hits for one annotation key.
    const NodeHits& fetch(node_index key);
    // Fetch the misses among |keys| into the cache without materialising hits.
    // Counts the rows reconstructed but no requests or cache hits (lookahead work
    // must not change the per-request counters); evicts the cache when it would
    // overflow.
    void warm(const std::vector<node_index> &keys);

    void clear_cache() { cache_.clear(); }

  private:
    enum class Path { DIRECT, ROWS, TUPLES };

    const LabelOracle &oracle_;
    std::vector<LabelRef> labels_;
    bool with_coords_;
    Path path_;
    size_t max_cache_size_;

    // column -> label id for COLUMN labels
    tsl::hopscotch_map<Column, LabelId> column_labels_;
    // column -> (seq_id -> label id) for HEADER labels
    tsl::hopscotch_map<Column, tsl::hopscotch_map<uint64_t, LabelId>> header_labels_;
    std::vector<Column> direct_columns_;  // sorted distinct columns for DIRECT/ROWS

    tsl::hopscotch_map<node_index, NodeHits> cache_;
    NodeHits empty_;

    void fetch_uncached(const std::vector<node_index> &keys);
    void hits_from_row(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                       NodeHits *hits) const;
    void hits_from_tuples(const annot::matrix::MultiIntMatrix::RowTuples &row,
                          NodeHits *hits) const;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_LABEL_ORACLE_HPP__
