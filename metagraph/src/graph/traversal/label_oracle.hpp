#ifndef __TRAVERSAL_LABEL_ORACLE_HPP__
#define __TRAVERSAL_LABEL_ORACLE_HPP__

#include <algorithm>
#include <chrono>
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
#include "annotation/binary_matrix/row_diff/row_diff_cache.hpp"
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
 * Paces the annotation reads of one request under a deadline (the chunked deadlines of spec
 * §6.8). A read the deadline cannot fall into is decoded in one piece: splitting a read costs
 * what its rows share — on a row-diff annotation the rows of one call share the decoding of
 * their row-diff paths, which every chunk decodes again (one-row chunks make a 1.8 s row_diff
 * walk take 30 s) — and buys nothing when the deadline is far. A read that might reach the
 * deadline is split: a first chunk of at most |first_rows| measures this read's own rows
 * (another read's rows, another level's, may be far cheaper: a first chunk sized from them can
 * take 650 ms of a 50 ms target), each next chunk at most 4 times the previous one and sized
 * at the rate the previous one measured to take min(|target_ms|, the time left), with the
 * deadline checked before each; once the rest is predicted to end well before the deadline it
 * is one piece.
 * One chunk, at least one row, stays uninterruptible. It also measures every read, paced or
 * not: |max_read_ms| is the longest single piece of decoding (and of building what it returns)
 * the request made — what stating an uninterruptible step needs. Chunking never changes what a
 * read returns: the callers keep a read's counting and cache decisions whole and only split
 * its decoding.
 */
/**
 * One uninterruptible piece of a seed's work, as its deadline record states it (SPEC §6.8;
 * timing only): what ran between two readings of the deadline — |kind| one of "read" (an
 * annotation read in one piece), "chunk" (a chunk of a split read), "rest" (the piece that ended
 * a split read), "kmer_mapping" (the seed's k-mers mapped to nodes and keys), "coord_mapping" (a
 * derived seed's k-mer whose coordinates were mapped to headers), "derivation_step" (a derived
 * seed's k-mer otherwise), "head" (the walk between two readings of the clock, reads excluded;
 * the last one ends at the walk's stop or end), "setup" (the seed phase's own processing between
 * its other pieces: validation, resolving label names and their duplicate check, the extra
 * labels, the depth-0 state — everything up to the walk's first checkpoint that no other piece
 * covers) and "finalisation" (from the walk's end or stop to its result) — with its rows and the
 * coordinates mapped to headers in it. The pieces cover the seed's time from its start to its
 * result (a read inside a head is a piece of its own, which the head excludes).
 */
struct UninterruptiblePiece {
    double ms = 0;
    const char *kind = "";
    uint64_t rows = 0;
    uint64_t coordinates = 0;
};

struct DecodePacer {
    double target_ms = 0;          // 0: pacing off, one piece per read
    size_t first_rows = 8;         // the most rows of a split read's first chunk
    double ms_per_row = 0;         // the slowest per-row time of any piece of the request
    double max_read_ms = 0;        // the longest single piece seen (ms)
    // A read not split yet is one piece when it is predicted, at the slowest per-row time of
    // the request (small reads and chunks make that pessimistic, which costs no more than a
    // first chunk), to take less than 1/|far_factor| of the time left: rows up to that many
    // times slower than any before still end before the deadline. Once a chunk measured
    // this read's own rows, the rest is one piece when predicted at that chunk's rate to take
    // less than 1/|rest_factor| of the time left. Tests set both to infinity to split every
    // read, deadline or none, into the smallest chunks.
    double far_factor = 64;
    double rest_factor = 4;
    // The rows of the next piece of a read with |remaining| rows left: all of them when the
    // deadline cannot fall into them (above; no deadline: |ms_left| infinite), else a chunk —
    // at least 1, at most |first_rows| for the read's first (|previous| = 0) and 4 x |previous|
    // after (|previous|: the read's previous chunk's rows, |previous_ms| its time), sized to
    // take min(target_ms, |ms_left|) at the rate the previous chunk measured
    size_t next(size_t remaining, double ms_left, size_t previous, double previous_ms) const;
    // a piece of |rows| rows took |ms| (|kind| and |coordinates| as UninterruptiblePiece's;
    // the coordinate mapping of a read's rows happens inside its piece, so its rate includes
    // it)
    void record(size_t rows, double ms, const char *kind = "read", uint64_t coordinates = 0);
    // The longest uninterruptible piece of the seed being walked (the walker resets it at the
    // seed's start), reads and the other pieces the walker notes; the total time of the reads
    // (record), which the walker's head pieces exclude
    UninterruptiblePiece longest;
    double read_ms = 0;
    void note_piece(const char *kind, double ms, uint64_t rows = 0, uint64_t coordinates = 0);
    // A "head" piece (the walk between two readings of the clock for a stop, its reads
    // excluded), noted as note_piece does and kept apart as the longest of the request, as
    // |max_read_ms| keeps the reads': the attempt states the larger of the two as
    // observed_max_uninterruptible_ms. A head is as uninterruptible as a read — a cancel or a
    // walk-until that falls inside it is seen only at its end — and counting only the
    // reads would let a lookahead's 60,000-node chains run 2 s to 11 s unpolled while the
    // observation says 1 ms
    double max_head_ms = 0;
    void note_head(double ms) {
        note_piece("head", ms);
        max_head_ms = std::max(max_head_ms, ms);
    }
    // the longest read or head piece of the request (ms): what an attempt observes as its
    // longest uninterruptible piece of walking (the delivery's gaps are the server's)
    double max_uninterruptible_ms() const { return std::max(max_read_ms, max_head_ms); }
    // The seed phase's own processing. While a setup is open (the walker opens it at a seed's
    // start and closes it at the walk's first checkpoint, or at its stop or end, or on its way
    // out), every piece noted first notes the time between the end of the previous piece
    // (|setup_since_ms|, on now_ms()'s clock) and its own start as a "setup" piece: with only the
    // seed phase's reads and k-mer mappings as pieces, a request naming 2,500 headers with a
    // shared 1,024-character prefix would state a longest piece of 0.297 ms for a seed phase of
    // 126 ms, almost all of it a duplicate check between two pieces
    bool setup_open = false;
    double setup_since_ms = 0;
    void open_setup(double since_ms) {
        setup_open = true;
        setup_since_ms = since_ms;
    }
    // the last setup span, from the end of the previous piece to |until_ms|, and the setup closed
    void close_setup(double until_ms);
    // the coordinates mapped so far (LabelOracle::Counters::coords_mapped; 0 without it), for a
    // piece's coordinates
    const uint64_t *coords_mapped = nullptr;
    uint64_t coords_now() const { return coords_mapped ? *coords_mapped : 0; }
    // The clock the pieces are measured with and the walker's deadlines read (ms): the steady
    // clock, or |test_clock_ms| when a test sets one (a virtual clock its slow reads advance,
    // so that what a deadline stops does not depend on the machine's load)
    std::function<double()> test_clock_ms;
    double now_ms() const {
        if (test_clock_ms)
            return test_clock_ms();
        return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now().time_since_epoch()).count();
    }
};

// The time since construction on a pacer's clock (DecodePacer::now_ms), ms
class PacerTimer {
  public:
    explicit PacerTimer(const DecodePacer &pacer) : pacer_(pacer), start_(pacer.now_ms()) {}
    double elapsed_ms() const { return pacer_.now_ms() - start_; }

  private:
    const DecodePacer &pacer_;
    double start_;
};

/**
 * What a paced read checks between its chunks (the caller's deadlines): |ms_left| the time
 * left before the nearest one (infinity: none), |stop| whether to stop before the next chunk
 * (also before the first). A read without it (null) is one piece, measured all the same.
 * An interrupted read returns nothing and changes no counter and no result (as a refusal of the
 * budget-aware reads: nothing appended, the budget as on entry); it sets |interrupted|, and
 * |units| to the work units of the rows its finished chunks decoded (8 per key, 1 per entry
 * and coordinate, and on the budget-aware reads the rows' dependency units): decoded work the
 * caller still charges, though no row was returned.
 */
struct ReadPacing {
    std::function<double()> ms_left;
    std::function<bool()> stop;
    bool interrupted = false;
    uint64_t units = 0;
};

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
        // the row-diff path cache's physical work (timing only; all 0 without it): paths cut
        // at a cached row, stored rows read from the matrix by the reads with the cache, rows
        // kept in the cache and their bytes, and the most it held so far (not a sum)
        uint64_t path_cache_hits = 0;
        uint64_t path_cache_stored_rows = 0;
        uint64_t path_cache_rows_kept = 0;
        uint64_t path_cache_bytes_kept = 0;
        uint64_t path_cache_peak_bytes = 0;
        // seed phases (Walker): the time of the seed's validation or derivation, of its label
        // resolution (in it: a first header name builds the sequence header index), and of
        // the reads in it (part of fetch_seconds)
        double seed_phase_seconds = 0;
        double label_resolve_seconds = 0;
        double seed_fetch_seconds = 0;
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
    // The sequence of a column coordinate with its first and last column coordinates
    // (CoordToHeader::sequence_range). Counts no mapping. Requires a CoordToHeader.
    struct SeqRange {
        uint64_t seq_id = 0;
        Coord first = 1;      // empty until set
        Coord last = 0;
    };
    SeqRange sequence_range(Column column, Coord coord) const;
    // the first and last column coordinates of sequence |seq_id| of |column| (its range, as
    // sequence_range states it). Requires a CoordToHeader.
    SeqRange sequence_coords(Column column, uint64_t seq_id) const;
    /**
     * map_coord of the coordinates of one column, one after another, by the runs they form
     * (the same (seq_id, local) as map_coord). The coordinates of one
     * sequence occupy a contiguous range, so in a tuple row's sorted coordinates a sequence
     * that holds several (rRNA operons: 7 copies in a genome) costs one rank and two selects
     * instead of a rank and a select each, and a coordinate in the sequence after the
     * previous one's (a single-copy gene in consecutive genomes) one select. The ranges are
     * only used once two coordinates in a row fell into one or adjacent sequences, and only
     * for lists of 8 coordinates or more: scattered coordinates and short lists
     * are mapped as map_coord maps them, at its cost (LabelOracleCoordRuns.DISABLED_Benchmark,
     * 2,000 sequences: runs of 7 at 22 ns a coordinate against 42-46 for map_coord, single
     * copies in consecutive sequences at 32-35 ns against 48, scattered ones at 40-42 against
     * 39-42; ranges used unconditionally are 1.3-1.5 times slower than map_coord there).
     */
    class CoordRuns {
      public:
        // |count|: how many coordinates will be mapped (a list of fewer than 8 is mapped as
        // map_coord maps it: it does not repay finding its runs — the mini refseq rows' 2-6
        // a column are 10-25% slower by the runs)
        CoordRuns(const annot::CoordToHeader &cth, Column column, size_t count);
        std::pair<uint64_t, Coord> map(Coord coord);
      private:
        const annot::CoordToHeader &cth_;
        const Column column_;
        const uint64_t num_sequences_;
        const bool short_;          // a short list: no ranges
        bool ranges_ = false;       // runs were seen: use the sequence ranges
        uint32_t streak_ = 0;       // coordinates in a row in the same or the next sequence
        bool have_ = false;         // |seq_| and |first_| hold the previous coordinate's
        bool last_known_ = false;   // |last_| holds its sequence's last coordinate
        uint64_t seq_ = 0;
        Coord first_ = 0;
        Coord last_ = 0;
    };
    // map_coord of every coordinate in coords[0, n) of |column|, in order, passed to
    // f(coord, seq_id, local), by CoordRuns (counted as n mappings)
    // (no buffer: the budget-aware reads charge every buffer they hold)
    template <class F>
    void map_coords(Column column, const Coord *coords, size_t n, const F &f) const {
        CoordRuns runs = coord_runs(column, n);
        for (size_t i = 0; i < n; ++i) {
            const auto [seq_id, local] = runs.map(coords[i]);
            f(coords[i], seq_id, local);
        }
        counters_.coords_mapped += n;
    }
    // the run mapper of |count| coordinates of |column|. Requires a CoordToHeader.
    CoordRuns coord_runs(Column column, size_t count) const;
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

    // The budget-aware reads (DESIGN-traverse-graphlet.md §14): whether
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
    // counters() with the path cache's physical counters brought up to date
    void sync_path_cache_counters() const;
    // the request's read pacing and measurement (the oracle is per request)
    DecodePacer& pacer() const { return pacer_; }

    /**
     * The row-diff path cache (row_diff_cache.hpp): on a row-diff
     * annotation every read above — the default and the budget-aware ones — keeps rows it
     * reconstructs (the requested rows, the rows just after them on their row-diff paths,
     * every checkpoint row and every narrow row: RowDiffCache::keeps), so that a later read's
     * path stops at a cached row instead of decoding to its anchor again. What
     * a read returns, and every RowCost (the work and memory a budget-aware read is charged
     * and admitted by: its whole path), do not depend on it; the decoding work does. Off
     * until set_path_cache_max() (the request's bound, --traverse-path-cache-mb): the walker
     * uses it within that bound, and under a memory budget within what the label cache
     * leaves of its allotment (Walker; a seed's own then, emptied when the seed starts).
     */
    void set_path_cache_max(uint64_t bytes) const {
        path_cache_max_ = bytes;
        path_cache_.set_bound(bytes);
    }
    uint64_t path_cache_max() const { return path_cache_max_; }
    annot::matrix::RowDiffPathCache& path_cache() const { return path_cache_; }
    // whether the reads use the path cache: a row-diff annotation with it, and a bound set
    bool path_cached() const;
    // Before a cache sharing the path cache's bound (a label cache under a memory budget)
    // grows by |bytes|: the path cache is trimmed to what then remains of the bound
    void make_room(uint64_t bytes) const {
        if (!path_cache_.enabled())
            return;
        const uint64_t limit = path_cache_.limit();
        path_cache_.trim(limit > bytes ? limit - bytes : 0);
    }
    // Tests: called with the row count at every annotation read of this oracle (the default
    // and the budget-aware ones), before the read — to make the reads of a small index slow
    std::function<void(size_t rows)> test_read_hook;
    // Tests: called with the name at every resolve_label — to make resolving a seed's labels
    // slow on a virtual clock (the seed phase's setup pieces)
    std::function<void(const std::string &name)> test_resolve_hook;

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
    mutable DecodePacer pacer_;
    mutable uint64_t path_cache_max_ = 0;
    mutable annot::matrix::RowDiffPathCache path_cache_;
};

/**
 * What a key read by the budget-aware path costs (DESIGN-traverse-graphlet.md §14),
 * kept with it in the caches: properties of the key, the same whether a fetch, the
 * lookahead or an earlier level decoded it, so that work and memory stops do not depend on
 * annotation.batch_kmers. |dependency_units|: the work of its row-diff dependency rows (8 per
 * row and 1 per entry and coordinate they store). |demand|: an upper bound of what reading
 * it ALONE holds at its peak — the decode, what building its result holds beside the raw row,
 * and the result — which every key a budgeted fetch returns is admitted against. |entries|:
 * the size of the row the read decoded, before anything was taken from it — its columns, and
 * for a tuple row its columns plus their coordinates —, so that a caller can charge a read
 * of a restricted LabelQuery as the whole row it decoded (on the row-diff family asking for
 * fewer labels does not make a row cheaper), not as the hits it returned. Nothing in this
 * file charges |entries|.
 *
 * Its size is part of the memory model (the costs vectors and LabelRecorder's demand charge
 * sizeof(KeyCost) per key; the cache entries kCostEntryBytes), so |entries| takes no room:
 * both unit counts are 32 bits, saturated at 2^32 - 1 (saturate_units()). A read reaches it
 * only holding gigabytes at once: its demand counts the stored rows of its path (a row's
 * slot, and every entry and coordinate they store, at least 4 bytes) and the decoded row
 * (8 bytes per column or coordinate).
 */
struct KeyCost {
    uint32_t dependency_units = 0;
    uint32_t entries = 0;
    uint64_t demand = 0;
};
static_assert(sizeof(KeyCost) == 16, "the memory model charges sizeof(KeyCost) per key");

// |units| as KeyCost stores them: at most 2^32 - 1
inline uint32_t saturate_units(uint64_t units) {
    return static_cast<uint32_t>(std::min<uint64_t>(units, std::numeric_limits<uint32_t>::max()));
}

// The longest run of keys one budget-aware read decodes together: a run that does not fit
// is retried in halves, so this bounds the wasted decoding near the budget
constexpr size_t kMaxDecodeRun = 512;

/**
 * Why the last budget-aware fetch of a LabelQuery or LabelRecorder was refused, so that the
 * caller can state the cause truthfully (a row that does not fit is not a dictionary that
 * does not fit): the refused key, what was left at its position and
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
        // no refusal of the budget: a paced fetch's deadline came before the key's run
        // (ReadPacing), restored as a refusal is
        INTERRUPTED,
    };
    Cause cause = DECODE;
    size_t position = 0;
    uint64_t left = 0;          // the budget's maximum minus what was held at entry and |held|
    uint64_t held = 0;          // what the keys before it returned (their hits or lists, names)
    uint64_t demand = 0;        // DEMAND, NAMES: the key's demand
    uint64_t need = 0;          // DECODE: the least its read was seen to need (> |left|)
    uint64_t labels = 0;        // NAMES: the new labels
    uint64_t names_bytes = 0;
    // The work units (8 per key, 1 per entry and coordinate, and the key's dependency units,
    // as a returned key's) of the keys the refused call decoded and built — those before the
    // refused key, the refused key itself when it was built (DEMAND, NAMES), and the keys of
    // its run built after it — but not of the keys the cache held: decoding done though
    // nothing was returned, which a caller with a work budget still charges. A key whose own
    // read or build did not fit (DECODE) adds nothing: its units
    // are not known. INTERRUPTED: ReadPacing::units.
    uint64_t units = 0;
};

/**
 * The (column, seq_id) keys of header labels, hashed for one flat table: a map per column
 * holds a whole bucket array (~1.5 KB with tsl's neighbourhood) even for a single label,
 * which a per-label charge cannot cover. Column and sequence ids are
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

    // Throws std::invalid_argument if the requested access path is not available for this
    // annotation (no silent fallback). AUTO picks DIRECT for at most 16 columns when the
    // annotation has direct access, which the budget-aware fetch and warm do not read
    // (they decode whole rows); no annotation they serve has direct access, so AUTO is
    // ROWS there.
    LabelQuery(const LabelOracle &oracle,
               std::vector<LabelRef> labels,
               bool with_coords,
               LabelOracle::Access access = LabelOracle::Access::AUTO,
               size_t max_cache_size = 1'000'000);

    const std::vector<LabelRef>& labels() const { return labels_; }
    bool with_coords() const { return with_coords_; }
    // which accessor is used: "direct", "rows" or "tuples"
    const char* access_path() const;

    // Hits for every annotation key (npos yields no hits). |pacing|: the misses are decoded
    // in chunks with its stop checked before each (ReadPacing; the counters, the cache's
    // eviction and the result are those of one call); interrupted: nothing is returned, no
    // counter changes, and the chunks decoded stay cached.
    std::vector<NodeHits> fetch(const std::vector<node_index> &keys,
                                ReadPacing *pacing = nullptr);
    // Hits for one annotation key.
    const NodeHits& fetch(node_index key);
    // Cache the hits of |keys| (npos and cached keys skipped) from rows the caller decoded
    // already — whole rows (|rows|[i] the row of keys[i], ROWS or DIRECT access) or tuple
    // rows (TUPLES access) —, built as a fetch builds them, so that a later fetch of the
    // keys returns from the cache without decoding them again (/resolve discovery). Counts
    // no request, cache hit or fetched row; ignores the cache bounds (the caller's working
    // set, as one call's).
    void prime(const std::vector<node_index> &keys,
               const std::vector<annot::matrix::BinaryMatrix::SetBitPositions> &rows);
    void prime(const std::vector<node_index> &keys,
               const std::vector<annot::matrix::MultiIntMatrix::RowTuples> &rows);
    // Fetch the misses among |keys| into the cache without materialising hits.
    // Counts the rows reconstructed but no requests or cache hits (lookahead work
    // must not change the per-request counters); evicts the cache when it would
    // overflow. |pacing|: as fetch(); interrupted, the warming ends silently.
    void warm(const std::vector<node_index> &keys, ReadPacing *pacing = nullptr);

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
    // |pacing|: its runs are at most a paced chunk long, its stop checked before each; an
    // interrupted fetch is refused like a refusal (cause INTERRUPTED: nothing appended, the
    // cache, the counters and |budget| as on entry), with ReadPacing::interrupted set.
    bool fetch(const node_index *keys, size_t n, annot::matrix::DecodeBudget &budget,
               std::vector<NodeHits> *out, std::vector<KeyCost> *costs, size_t *refused_at,
               ReadPacing *pacing = nullptr);
    // The budget-aware lookahead: decodes the misses among |keys| in runs within |budget|
    // and caches them with their costs within the byte bound; a run that does not fit
    // ends the warming silently (nothing a later fetch returns depends on it); |pacing|:
    // each run decoded in paced pieces, its stop checked before each (interrupted: the run
    // is dropped and the warming ends silently). The warming also ends at the first run the
    // cache cannot keep beside the runs this warm cached (it never evicts its own runs) or
    // cannot keep at all: what was cached before the warm may be evicted for its first run
    // only
    void warm(const std::vector<node_index> &keys, annot::matrix::DecodeBudget &budget,
              ReadPacing *pacing = nullptr);
    // what a copy of |hits| holds (the model of decode_budget.hpp)
    static uint64_t held_bytes(const NodeHits &hits);
    // why the last budgeted fetch was refused, and what the keys before the refused one held
    const FetchRefusal& refusal() const { return refusal_; }
    uint64_t refused_held() const { return refusal_.held; }
    // Under a byte bound (set_max_cache_bytes), what the last unbudgeted call's raw rows held
    // while their hits were built (an estimate, as the cache's): what such an unbudgeted
    // read holds beyond the account, observed by the walker (memory_bound_soft). 0 without
    // a bound.
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
    // Permitted-range filtering: per header column, the coordinate ranges of the requested
    // sequences, ascending, with their labels — a coordinate is looked up among them (a binary
    // search over the query's own sequences of that column) instead of being mapped to its
    // sequence by a rank and a select, which every coordinate of the column would cost, those
    // of the sequences no label names included
    struct HeaderRange {
        Coord first;
        Coord last;
        LabelId label;
    };
    tsl::hopscotch_map<Column, std::vector<HeaderRange>> header_ranges_;
    // f(label, local) for each coordinate of |coords| (column |column|) in a requested
    // sequence; every coordinate counts as mapped (Counters::coords_mapped), as map_coords
    template <class F>
    void map_requested(Column column, const Coord *coords, size_t n, const F &f) const;
    std::vector<Column> direct_columns_;  // sorted distinct columns for DIRECT/ROWS

    tsl::hopscotch_map<node_index, NodeHits> cache_;
    NodeHits empty_;

    void fetch_uncached(const node_index *keys, size_t n);
    // fetch_uncached of |keys| (sorted, distinct) in paced chunks (one piece without
    // |pacing|, or when the deadline cannot fall into it), each measured, chunks taken in the
    // order of first appearance in |order| (the caller's keys); false when |pacing|'s stop came
    // first, with ReadPacing::units the work of the keys decoded
    bool fetch_uncached_paced(const std::vector<node_index> &keys,
                              const std::vector<node_index> &order, ReadPacing *pacing);
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
    // What cache_budgeted may evict to keep a call's keys that do not fit beside the cache:
    // everything (EVICT: a level's read, the walk moving forward; the cache is emptied even
    // when the keys alone do not fit, as always), everything but only when the keys then fit
    // (EVICT_IF_KEPT: a warm's first run, which drops only entries from before the warm), or
    // nothing (KEEP: a warm's later runs, which must not drop the runs before them)
    enum class Eviction { EVICT, EVICT_IF_KEPT, KEEP };
    // cache the keys of keys[0, n) that are not cached (those a successful budgeted call
    // decoded: the cache does not change during the call), within the byte bound; whether
    // they are cached now
    bool cache_budgeted(const node_index *keys, size_t n, const NodeHits *hits,
                        const KeyCost *costs, Eviction eviction = Eviction::EVICT);
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
    // here, in |keys| order. |pacing|: as LabelQuery::fetch (interrupted: nothing returned
    // and nothing named).
    std::vector<NodeLabels> fetch(const std::vector<node_index> &keys,
                                  ReadPacing *pacing = nullptr);
    // Fetch the misses into the row cache without naming anything (lookahead).
    void warm(const std::vector<node_index> &keys, ReadPacing *pacing = nullptr);

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
               const std::function<uint64_t(std::string_view name)> &name_bytes,
               ReadPacing *pacing = nullptr);
    // as LabelQuery's, the warming ending at the first run the cache cannot keep beside this
    // warm's runs or at all
    void warm(const std::vector<node_index> &keys, annot::matrix::DecodeBudget &budget,
              ReadPacing *pacing = nullptr);
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
    // Every eviction, by either path, drops the rows and their costs together, so that every
    // key with a cost is a cached key: the budget-aware path's size check (equal sizes, equal
    // keys) relies on it. Clearing the rows alone would let an ordinary fetch leave a stale
    // cost behind, and a budgeted fetch of an equal-sized cache would then find a row
    // without its cost (std::out_of_range for a valid key)
    void clear_cache() {
        cache_.clear();
        costs_.clear();
        cached_keys_ = 0;
        cache_bytes_ = 0;
    }

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

    void fetch_uncached(const node_index *keys, size_t n);
    // as LabelQuery's
    bool fetch_uncached_paced(const std::vector<node_index> &keys,
                              const std::vector<node_index> &order, ReadPacing *pacing);
    LabelId id_of(const Key &key);

    uint64_t last_call_bytes_ = 0;

    // ---- the budget-aware path (as LabelQuery's)
    tsl::hopscotch_map<node_index, KeyCost> costs_;
    FetchRefusal refusal_;
    uint64_t last_names_bytes_ = 0;
    uint64_t last_naming_bytes_ = 0;
    // the dictionary id of |key|, if it is named
    std::optional<LabelId> named(const Key &key) const;
    // decode keys[0, n) (misses) as one run into (*rows)[at + i] (the cache's form) and
    // (*rows_held)[at + i], (*costs)[at + i], charged
    annot::matrix::DecodeStatus decode_run(const node_index *keys, size_t n,
                                           annot::matrix::DecodeBudget &budget,
                                           std::vector<RawRow> *rows,
                                           std::vector<uint64_t> *rows_held,
                                           std::vector<KeyCost> *costs, size_t *built,
                                           size_t at = 0);
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
