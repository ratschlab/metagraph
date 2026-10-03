#ifndef __TRAVERSAL_LABEL_ORACLE_HPP__
#define __TRAVERSAL_LABEL_ORACLE_HPP__

#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "traversal_types.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {

namespace annot {
class CoordToHeader;
namespace matrix {
class IRowDiff;
}
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
    // Look up a FASTA header in the CoordToHeader index. The reverse index is owned by the
    // CoordToHeader (CoordToHeader::find_header): built on the first call, shared by all
    // oracles of the same loaded index, and gone with it.
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

    // The budget-aware reads (DESIGN-traverse-graphlet.md §14, stage 3 of §14.1): whether
    // this index has them — a row-diff annotation over BRWT or ColumnMajor, with or
    // without coordinates — and the reads themselves (IRowDiff::decode_rows /
    // decode_row_tuples: every buffer charged to |budget| before it is allocated, a read
    // that does not fit refused whole). The physical counters (timing) count the rows a read
    // returns and the time of every read, a refused one's included — also when the
    // LabelQuery / LabelRecorder fetch that issued it is refused afterwards: the decoding
    // was done. The unbudgeted reads above are untouched by these.
    bool decode_charged() const;
    // a row-diff annotation (whatever its index matrix)
    bool row_diff() const { return rd_; }
    annot::matrix::DecodeStatus
    get_rows(const std::vector<Row> &rows, annot::matrix::DecodeBudget &budget,
             std::vector<annot::matrix::BinaryMatrix::SetBitPositions> *out,
             std::vector<annot::matrix::RowCost> *costs, std::vector<uint64_t> *held) const;
    annot::matrix::DecodeStatus
    get_row_tuples(const std::vector<Row> &rows, annot::matrix::DecodeBudget &budget,
                   std::vector<annot::matrix::MultiIntMatrix::RowTuples> *out,
                   std::vector<annot::matrix::RowCost> *costs, std::vector<uint64_t> *held) const;

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
    const annot::matrix::IRowDiff *rd_ = nullptr;
    const annot::CoordToHeader *coord_to_header_ = nullptr;

    mutable Counters counters_;
};

/**
 * What a key read by the budget-aware path costs (DESIGN-traverse-graphlet.md §14, stage 3),
 * kept with it in the caches: properties of the key, the same whether a fetch, the
 * lookahead or an earlier level decoded it, so that work and memory stops do not depend on
 * annotation.batch_kmers. |dependency_units|: the work of its row-diff dependency rows (8 per
 * row and 1 per entry and coordinate they store). |demand|: an upper bound of what reading
 * it ALONE holds at its peak — the decode, what building its result holds beside the raw row,
 * and the result — which every key a budgeted fetch returns is admitted against.
 */
struct KeyCost {
    uint64_t dependency_units = 0;
    uint64_t demand = 0;
};

// The longest run of keys one budget-aware read decodes together: a run that does not fit
// is retried in halves, so this bounds the wasted decoding near the budget
constexpr size_t kMaxDecodeRun = 512;

/**
 * Why the last budget-aware fetch of a LabelQuery or LabelRecorder was refused, so that the
 * caller can state the cause truthfully (a row that does not fit is not a dictionary that
 * does not fit; review of stage 3, F2): the refused key, what was left at its position and
 * what admitting it needed.
 */
struct FetchRefusal {
    enum Cause {
        // the key's read alone — its row with its row-diff dependency rows, and building its
        // hits or list — did not fit |left|: it needs more, at least |need| bytes
        DECODE,
        // the key was read, but its standalone demand (KeyCost::demand, the upper bound every
        // returned key is admitted against) exceeds |left|
        DEMAND,
        // LabelRecorder only: its demand fits |left|, but with the |labels| dictionary labels
        // it would name first (|names_bytes|: as priced, with their provisional naming) not
        NAMES,
    };
    Cause cause = DECODE;
    size_t position = 0;
    uint64_t left = 0;          // the budget's maximum minus what was held at entry and |held|
    uint64_t held = 0;          // what the keys before it returned (their hits or lists, names)
    uint64_t demand = 0;        // DEMAND, NAMES: the key's demand
    uint64_t need = 0;          // DECODE: the least its read was seen to need (> |left|)
    uint64_t labels = 0;        // NAMES: the new labels
    uint64_t names_bytes = 0;
};

/**
 * The (column, seq_id) keys of header labels, hashed for one flat table: a map per column
 * holds a whole bucket array (~1.5 KB with tsl's neighbourhood) even for a single label,
 * which a per-label charge cannot cover (review of stage 3, F1). Column and sequence ids are
 * small consecutive integers: mixed, so that a power-of-two table does not see them raw.
 */
struct LabelKeyHash {
    size_t operator()(const std::pair<Column, uint64_t> &key) const {
        uint64_t h = (key.first + 0x632BE59BD9B4E019ULL) * 0x9E3779B97F4A7C15ULL;
        h ^= (key.second + 0x7F4A7C159E3779B9ULL) * 0xBF58476D1CE4E5B9ULL;
        return static_cast<size_t>(h ^ (h >> 31));
    }
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

    /**
     * The budget-aware fetch (LabelOracle::decode_charged()), all or nothing: the hits of
     * keys[0, n) APPENDED to |*out| and their costs to |*costs| (both reserved for them by
     * the caller, so that nothing here depends on how a level is cut into calls). Keys are
     * taken in order, each against what |budget| has left at its position — its maximum
     * minus the hits returned for the keys before it: a key is admitted when its demand
     * (KeyCost) fits, whether it is cached or decoded now, so that where a fetch stops does
     * not depend on the cache (and so on the lookahead). Misses are decoded in runs, every
     * buffer charged before it is allocated; a run that does not fit is retried in halves.
     * Returns false at the first key that does not fit, with |*refused_at| its position:
     * then nothing is appended, and the cache, the counters and |budget| are as on entry.
     * On success |budget| holds the returned hits (held_bytes()), which the caller keeps
     * charged until it frees them, and the newly decoded keys are cached with their costs
     * when they fit the byte bound (never beyond it).
     */
    bool fetch(const node_index *keys, size_t n, annot::matrix::DecodeBudget &budget,
               std::vector<NodeHits> *out, std::vector<KeyCost> *costs, size_t *refused_at);
    // The budget-aware lookahead: decodes the misses among |keys| in runs within |budget|
    // and caches them with their costs within the byte bound; a run that does not fit
    // ends the warming silently (nothing a later fetch returns depends on it)
    void warm(const std::vector<node_index> &keys, annot::matrix::DecodeBudget &budget);
    // what a copy of |hits| holds (the model of decode_budget.hpp)
    static uint64_t held_bytes(const NodeHits &hits);
    // why the last budgeted fetch was refused, and what the keys before the refused one held
    const FetchRefusal& refusal() const { return refusal_; }
    uint64_t refused_held() const { return refusal_.held; }
    // Under a byte bound (set_max_cache_bytes), what the last unbudgeted call's raw rows held
    // while their hits were built (an estimate, as the cache's): what a stage-2 read holds
    // beyond the account, observed by the walker (memory_bound_soft). 0 without a bound.
    uint64_t last_call_bytes() const { return last_call_bytes_; }

    void clear_cache() { cache_.clear(); costs_.clear(); cache_bytes_ = 0; }

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
    // (column, seq_id) -> label id for HEADER labels, and the columns they are in (one flat
    // table each, see LabelKeyHash)
    tsl::hopscotch_map<std::pair<Column, uint64_t>, LabelId, LabelKeyHash> header_labels_;
    tsl::hopscotch_set<Column> header_columns_;
    std::vector<Column> direct_columns_;  // sorted distinct columns for DIRECT/ROWS

    tsl::hopscotch_map<node_index, NodeHits> cache_;
    NodeHits empty_;

    void fetch_uncached(const std::vector<node_index> &keys);
    void hits_from_row(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                       NodeHits *hits) const;
    void hits_from_tuples(const annot::matrix::MultiIntMatrix::RowTuples &row,
                          NodeHits *hits) const;

    uint64_t last_call_bytes_ = 0;

    // ---- the budget-aware path
    // the costs of the keys cached by the budget-aware path (beside |cache_|, whose layout
    // the unbudgeted path keeps)
    tsl::hopscotch_map<node_index, KeyCost> costs_;
    FetchRefusal refusal_;
    // per-label scratch of hits_budgeted(), sized |labels_| on first use
    std::vector<uint32_t> label_count_;
    std::vector<LabelId> touched_;
    // Decode the keys[0, n) (misses, not npos) as one run and build their hits at
    // (*out)[at + i], charged; |*costs| their costs. REFUSED: nothing built, |budget| as on
    // entry. |*built| counts the keys built before a conversion did not fit (OK only when
    // all were built).
    annot::matrix::DecodeStatus decode_run(const node_index *keys, size_t n,
                                           annot::matrix::DecodeBudget &budget,
                                           std::vector<NodeHits> *out, size_t at,
                                           std::vector<KeyCost> *costs, size_t *built);
    // the hits of one raw row built with exact reservations, every buffer charged before
    // it is allocated; |*peak| what building it held beside the raw row at most
    bool hits_budgeted(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                       annot::matrix::DecodeBudget &budget, NodeHits *hits, uint64_t *peak);
    bool hits_budgeted(const annot::matrix::MultiIntMatrix::RowTuples &row,
                       annot::matrix::DecodeBudget &budget, NodeHits *hits, uint64_t *peak);
    // cache the keys of keys[0, n) that are not cached (those a successful budgeted call
    // decoded: the cache does not change during the call), within the byte bound
    void cache_budgeted(const node_index *keys, size_t n, const NodeHits *hits,
                        const KeyCost *costs);
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

    /**
     * The budget-aware fetch, as LabelQuery's: all or nothing, the labels of keys[0, n)
     * appended to |*out|, their costs to |*costs|, each key admitted against what |budget|
     * has left at its position (its maximum minus the lists returned for the keys before it
     * and the names they gave) — its demand and the names it gives first: per new label what
     * |name_bytes| prices for its name (what the caller's account will charge for the
     * dictionary label) and kNamingBytes for naming it provisionally. The names are given
     * only on success, in key order, so that a refused call names nothing. On success
     * |budget| holds the lists, the priced names and the naming charges; the last two are
     * last_names_bytes() and last_naming_bytes().
     */
    bool fetch(const node_index *keys, size_t n, annot::matrix::DecodeBudget &budget,
               std::vector<NodeLabels> *out, std::vector<KeyCost> *costs, size_t *refused_at,
               const std::function<uint64_t(std::string_view name)> &name_bytes);
    void warm(const std::vector<node_index> &keys, annot::matrix::DecodeBudget &budget);
    static uint64_t held_bytes(const NodeLabels &labels);
    // What naming one label provisionally can hold, charged per new label beside its priced
    // name: its share of the call's pending table and list, rehash and copy transients
    // included (bounds pending_bytes(m) <= m * kNamingBytes for every m), so that the charge
    // is a property of the label, not of how many the call named before it
    static constexpr uint64_t kNamingBytes = 640;
    // the most naming |m| labels provisionally holds at once (the model of decode_budget.hpp)
    static uint64_t pending_bytes(uint64_t m);
    // as LabelQuery's; and, of the last successful budgeted fetch, what its new labels' names
    // cost (by |name_bytes|; the caller's account charges them with the dictionary) and what
    // naming them provisionally was charged (freed with the call)
    const FetchRefusal& refusal() const { return refusal_; }
    uint64_t refused_held() const { return refusal_.held; }
    uint64_t last_names_bytes() const { return last_names_bytes_; }
    uint64_t last_naming_bytes() const { return last_naming_bytes_; }
    uint64_t last_call_bytes() const { return last_call_bytes_; }

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
    // one flat table for the header labels (LabelKeyHash): a map per column held a whole
    // bucket array per column, beyond what the account charges a dictionary label
    tsl::hopscotch_map<Key, LabelId, LabelKeyHash> header_ids_;
    tsl::hopscotch_map<node_index, RawRow> cache_;

    void fetch_uncached(const std::vector<node_index> &keys);
    LabelId id_of(const Key &key);

    uint64_t last_call_bytes_ = 0;

    // ---- the budget-aware path (as LabelQuery's)
    tsl::hopscotch_map<node_index, KeyCost> costs_;
    FetchRefusal refusal_;
    uint64_t last_names_bytes_ = 0;
    uint64_t last_naming_bytes_ = 0;
    // the dictionary id of |key|, if it is named
    std::optional<LabelId> named(const Key &key) const;
    // decode keys[0, n) (misses) as one run into |*rows| (the cache's form), charged
    annot::matrix::DecodeStatus decode_run(const node_index *keys, size_t n,
                                           annot::matrix::DecodeBudget &budget,
                                           std::vector<RawRow> *rows,
                                           std::vector<uint64_t> *rows_held,
                                           std::vector<KeyCost> *costs, size_t *built);
    bool raw_budgeted(const annot::matrix::BinaryMatrix::SetBitPositions &row,
                      annot::matrix::DecodeBudget &budget, RawRow *raw, uint64_t *peak);
    bool raw_budgeted(const annot::matrix::MultiIntMatrix::RowTuples &row,
                      annot::matrix::DecodeBudget &budget, RawRow *raw, uint64_t *peak);
    static uint64_t raw_bytes(const RawRow &raw);
    void cache_raw(node_index key, RawRow &&raw, const KeyCost &cost);
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __TRAVERSAL_LABEL_ORACLE_HPP__
