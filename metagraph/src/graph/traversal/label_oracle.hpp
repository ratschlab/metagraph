#ifndef __TRAVERSAL_LABEL_ORACLE_HPP__
#define __TRAVERSAL_LABEL_ORACLE_HPP__

#include <limits>
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

    void clear_cache() { cache_.clear(); cache_bytes_ = 0; }

    // A bound on the cache's bytes (an estimate: entries, hits and coordinates), on top
    // of the row bound, evicting wholesale like it. Unbounded unless set: a request with
    // a memory budget gives the cache a fixed allotment of that budget, which it must
    // then stay within (exceeded at most by the working set of one call).
    void set_max_cache_bytes(uint64_t bytes) { max_cache_bytes_ = bytes; }
    uint64_t cache_bytes() const { return cache_bytes_; }

  private:
    enum class Path { DIRECT, ROWS, TUPLES };

    const LabelOracle &oracle_;
    std::vector<LabelRef> labels_;
    bool with_coords_;
    Path path_;
    size_t max_cache_size_;
    uint64_t max_cache_bytes_ = std::numeric_limits<uint64_t>::max();
    uint64_t cache_bytes_ = 0;

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


/**
 * The labels PRESENT at annotation keys, with no permitted set: `labels.mode: annotate`
 * (spec §6.9). Every key costs a FULL row (or tuple row): there is nothing to narrow
 * the read to, which is why this is a verification tool for small radii and not a
 * search primitive. Labels are named on first sight, in the order the walk consumes
 * them, so the dictionary is deterministic for a deterministic walk and independent
 * of prefetching (warm() caches raw rows and assigns no ids).
 *
 * Per key the list is capped at |max_labels_per_node| — the first that many labels in
 * ascending (column, seq_id) order — and the TRUE count is returned beside it, so a cut
 * list is never mistaken for the whole set. For HEADER labels the count needs the
 * coordinates mapped to sequences, which costs one mapping per (column, sequence) the
 * k-mer occurs in, not one per coordinate (a sequence's coordinates are contiguous).
 *
 * The row cache is bounded in rows (|max_cache_size|) and in kept keys
 * (|max_cache_keys|, its memory: a row holds up to the cap), whichever trips first;
 * both evict wholesale.
 */
class LabelRecorder {
  public:
    struct NodeLabels {
        std::vector<LabelId> labels;   // ascending ids, at most max_labels_per_node
        size_t total = 0;              // distinct labels at the key; > labels.size() when cut
        bool truncated() const { return total > labels.size(); }
    };

    // Throws std::invalid_argument when |kind| is HEADER and the index has no
    // CoordToHeader, or when max_labels_per_node is 0.
    LabelRecorder(const LabelOracle &oracle,
                  LabelKind kind,
                  size_t max_labels_per_node,
                  size_t max_cache_size = 1'000'000,
                  size_t max_cache_keys = 64'000'000);

    LabelKind kind() const { return kind_; }
    size_t max_labels_per_node() const { return cap_; }
    // the labels named so far (LabelId == index)
    const std::vector<LabelRef>& labels() const { return dict_; }
    // "rows" or "tuples"
    const char* access_path() const;

    // Labels at every key (npos yields an empty list with total 0). Ids are assigned
    // here, in |keys| order.
    std::vector<NodeLabels> fetch(const std::vector<node_index> &keys);
    // Fetch the misses into the row cache without naming anything (lookahead).
    void warm(const std::vector<node_index> &keys);

    // as LabelQuery::set_max_cache_bytes
    void set_max_cache_bytes(uint64_t bytes) { max_cache_bytes_ = bytes; }
    uint64_t cache_bytes() const { return cache_bytes_; }

  private:
    // (column, seq_id); seq_id is 0 for COLUMN labels
    using Key = std::pair<Column, uint64_t>;
    struct RawRow {
        std::vector<Key> kept;         // the first |cap_| keys, ascending
        size_t total = 0;
    };

    const LabelOracle &oracle_;
    LabelKind kind_;
    size_t cap_;
    size_t max_cache_size_;
    size_t max_cache_keys_;
    size_t cached_keys_ = 0;           // sum over cached rows of max(1, |kept|)
    uint64_t max_cache_bytes_ = std::numeric_limits<uint64_t>::max();
    uint64_t cache_bytes_ = 0;
    std::vector<LabelRef> dict_;
    tsl::hopscotch_map<Column, LabelId> column_ids_;
    tsl::hopscotch_map<Column, tsl::hopscotch_map<uint64_t, LabelId>> header_ids_;
    tsl::hopscotch_map<node_index, RawRow> cache_;

    void fetch_uncached(const std::vector<node_index> &keys);
    LabelId id_of(const Key &key);
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_LABEL_ORACLE_HPP__
