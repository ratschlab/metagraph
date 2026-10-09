#ifndef __ROW_DIFF_CACHE_HPP__
#define __ROW_DIFF_CACHE_HPP__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include <tsl/hopscotch_map.h>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * What the budget-aware decode (IRowDiff::decode_rows) computes for a row over its WHOLE
 * row-diff path to the anchor — the aggregates of its visit (row_diff_budgeted.cpp), from
 * which a requested row's RowCost is made. Kept with a cached row, so that a path the cache
 * cuts at that row gets the costs of the path decoded to its anchor: the costs, and the work
 * and memory stops made from them, do not depend on what is cached.
 */
struct PathAggregates {
    uint64_t entries = 0;       // entries of the row's own stored row (a diff, or an anchor)
    uint32_t length = 0;        // rows on its path, itself included
    uint64_t dep_entries = 0;   // entries the other rows of its path store
    uint64_t path_stored = 0;   // bytes charged for the stored rows of its path
    uint64_t max_scratch = 0;   // the most scratch entries a row of its path took
    uint64_t recon_peak = 0;    // the most its reconstruction holds beyond the stored rows
    uint64_t full = 0;          // bytes charged for its reconstructed row
    uint64_t anchor = 0;        // the row its path ends at (BinaryMatrix::Row)
};

// What a copy of a row holds (decode_budget.hpp's model: each buffer at its exact size)
inline uint64_t row_copy_bytes(const BinaryMatrix::SetBitPositions &row) {
    return small_vector_bytes(row.size(), sizeof(BinaryMatrix::Column));
}
inline uint64_t row_copy_bytes(const MultiIntMatrix::RowTuples &row) {
    uint64_t bytes = buffer_bytes(row.size(), sizeof(row[0]));
    for (const auto &entry : row) {
        bytes += small_vector_bytes(entry.second.size(), sizeof(uint64_t));
    }
    return bytes;
}

// The entries of a row: its columns, and for a tuple row their coordinates too
inline uint64_t row_entries(const BinaryMatrix::SetBitPositions &row) { return row.size(); }
inline uint64_t row_entries(const MultiIntMatrix::RowTuples &row) {
    uint64_t entries = row.size();
    for (const auto &entry : row) {
        entries += entry.second.size();
    }
    return entries;
}

/**
 * How the cache stores a row. A binary row as it is (one buffer). A tuple row flattened into
 * three buffers — its columns, the end of each column's coordinates, and the coordinates —
 * rather than as a RowTuples copy, which allocates one buffer per column with more than two
 * coordinates and frees them one by one when the row is evicted: on a wide row of a coordinate
 * annotation (refseq33m 23S: ~11k columns, ~28k coordinates) those allocations are most of
 * what keeping a row costs. load() rebuilds the row exactly; copy_bytes is what a
 * RowTuples copy of it holds (row_copy_bytes), which the budget-aware decode charges for a hit.
 */
template <class RowT>
struct StoredRow;

template <>
struct StoredRow<BinaryMatrix::SetBitPositions> {
    BinaryMatrix::SetBitPositions row;
    void store(const BinaryMatrix::SetBitPositions &full) { row = full; }
    void load(BinaryMatrix::SetBitPositions *out) const { *out = row; }
    uint64_t bytes() const { return row_copy_bytes(row); }
    // what bytes() is once |full| is stored, computed from |full| without copying it (the
    // cache admits a row before it copies it)
    static uint64_t bytes_of(const BinaryMatrix::SetBitPositions &full) {
        return row_copy_bytes(full);
    }
};

template <>
struct StoredRow<MultiIntMatrix::RowTuples> {
    std::vector<BinaryMatrix::Column> columns;
    std::vector<uint32_t> ends;
    std::vector<uint64_t> coords;
    static uint64_t flat_bytes(uint64_t num_columns, uint64_t num_coords) {
        return buffer_bytes(num_columns, sizeof(BinaryMatrix::Column))
             + buffer_bytes(num_columns, sizeof(uint32_t))
             + buffer_bytes(num_coords, sizeof(uint64_t));
    }
    static uint64_t total_coords(const MultiIntMatrix::RowTuples &full) {
        uint64_t total = 0;
        for (const auto &entry : full) {
            total += entry.second.size();
        }
        return total;
    }
    // what bytes() is once |full| is stored (its three buffers at their exact sizes),
    // computed from |full| without copying it (the cache admits a row before it copies it)
    static uint64_t bytes_of(const MultiIntMatrix::RowTuples &full) {
        return flat_bytes(full.size(), total_coords(full));
    }
    void store(const MultiIntMatrix::RowTuples &full) {
        const uint64_t total = total_coords(full);
        assert(total <= std::numeric_limits<uint32_t>::max());
        columns.resize(full.size());
        ends.resize(full.size());
        coords.resize(total);
        uint64_t *out = coords.data();
        for (size_t i = 0; i < full.size(); ++i) {
            columns[i] = full[i].first;
            out = std::copy(full[i].second.begin(), full[i].second.end(), out);
            ends[i] = out - coords.data();
        }
    }
    void load(MultiIntMatrix::RowTuples *out) const {
        MultiIntMatrix::RowTuples row;
        row.reserve(columns.size());
        uint32_t begin = 0;
        for (size_t i = 0; i < columns.size(); ++i) {
            row.emplace_back(columns[i], MultiIntMatrix::Tuple(coords.begin() + begin,
                                                               coords.begin() + ends[i]));
            begin = ends[i];
        }
        out->swap(row);
    }
    uint64_t bytes() const {
        assert(ends.size() == columns.size());
        return flat_bytes(columns.size(), coords.size());
    }
};

/**
 * The row-diff path cache: reconstructed rows of a row-diff annotation — the rows a decode
 * call returned and the rows of their paths its retention rule selects (keeps()) — kept
 * across the decode calls of one request, so that a later call's path that meets a cached
 * row stops there instead of decoding to its anchor again. Without it, on refseq33m a fetch
 * of one warm 23S k-mer costs one whole path decode (66-86 ms) and each further row on the
 * same path 2.1 ms; a walk fetching a few keys per level pays 27-41 ms a row against 4.6 ms
 * read in one call. A row's content does not depend on how it was reached, so what a read
 * returns never depends on the cache; only the decoding work does.
 *
 * Bounded in bytes (StoredRow::bytes plus kEntryBytes per entry, an estimate of the table's
 * share), in two generations: inserts go to the current one, which becomes the older one
 * when it reaches half the bound (the older one is dropped then); a hit in the older one
 * moves the row to the current one. A walk moves forward, so the rows of the last few levels
 * are the ones its next paths meet. A generation dropped releases its table: a cleared
 * hopscotch map keeps its bucket array, sized for the most entries it ever held, so after
 * many narrow rows a few wide ones would hold that array beside their bound (74.7 MiB of
 * heap for 63.9 MiB accounted at a 64 MiB bound). Single-threaded, as the request that owns
 * it.
 */
template <class RowT>
class RowDiffCache {
  public:
    using Row = BinaryMatrix::Row;
    struct Entry {
        StoredRow<RowT> row;
        // what a copy of the row holds (row_copy_bytes): the budget-aware decode charges a
        // hit as its copy
        uint64_t copy_bytes = 0;
        uint64_t bytes = 0;
        // rows from it to its anchor, itself excluded (0: an anchor): what the retention rule
        // reads for the rows reconstructed from it (keeps())
        uint32_t depth = 0;
        // whether |path| is known: a row the budget-aware decode cached; the default decode
        // does not compute it, and the budget-aware decode treats such a row as not cached
        bool has_path = false;
        PathAggregates path;
    };
    // an entry's share of the tables besides its row's buffers: a slot of key and entry at
    // load 0.4 after a growth, with the transient of a rehash (both bucket arrays)
    static constexpr uint64_t kEntryBytes = 640;

    // The bound: |max_bytes|, or less while |room| (when set) returns less: the bytes a
    // caller sharing one allotment with another cache leaves it (re-read at every insert)
    void set_bound(uint64_t max_bytes, std::function<uint64_t()> room = {}) {
        max_bytes_ = max_bytes;
        room_ = std::move(room);
        trim(limit());
    }
    bool enabled() const { return max_bytes_ > 0; }
    uint64_t limit() const {
        return room_ ? std::min(max_bytes_, room_()) : max_bytes_;
    }
    uint64_t bytes() const { return bytes_[0] + bytes_[1]; }
    size_t size() const { return gen_[0].size() + gen_[1].size(); }
    // the tables' buckets, both generations (what the tests measure the tables by)
    size_t bucket_count() const { return gen_[0].bucket_count() + gen_[1].bucket_count(); }

    // The entry of |row| (with a known path only, if |need_path|), or null. A hit in the older
    // generation moves the row to the current one. The pointer is valid until the next call
    // of any member.
    const Entry* find(Row row, bool need_path) {
        auto it = gen_[0].find(row);
        if (it != gen_[0].end())
            return !need_path || it->second.has_path ? &it->second : nullptr;
        auto old = gen_[1].find(row);
        if (old == gen_[1].end() || (need_path && !old->second.has_path))
            return nullptr;
        // moved, not copied: the bytes stay within the bound
        Entry entry = std::move(old.value());
        gen_[1].erase(old);
        bytes_[1] -= entry.bytes;
        bytes_[0] += entry.bytes;
        return &gen_[0].emplace(row, std::move(entry)).first->second;
    }
    /**
     * Which of the rows a decode call reconstructs are kept. Keeping every row of every path
     * would copy each row of a long path into the cache: on a coordinate annotation with wide
     * rows (refseq33m 23S: ~28k coordinates in ~11k columns) a first read, whose paths no
     * later read meets, is then several times slower than without the cache, and a walk
     * reading one row per call (batch_kmers 1) churns the cache's generations so that every
     * call decodes its whole path again — 21.7 s against 2.9 s without the cache on a
     * synthetic index of 8,000 labels with paths of up to 1,000 rows. A call keeps
     *   - every row it was asked for (later paths meet them: a predecessor's path runs
     *     through its successor, and a repeated read is a hit),
     *   - the first |successors| rows after each of them on its path (a forward walk reads
     *     these next, in its following levels),
     *   - every row whose distance to its anchor is a multiple of |checkpoint|, anchors
     *     included: any later path through a row of this one stops within checkpoint - 1
     *     rows instead of running to the anchor; the distance is a property of the row (each
     *     row has one row-diff successor), so which rows these are does not depend on the
     *     calls, and
     *   - every narrow row (row_copy_bytes below |narrow_bytes|), whose copy costs little
     *     beside the descent that decoded it: on UHGG (binary rows of 10-200 bytes) keeping
     *     all rows is what made branched walks about twice as fast, and the rule keeps them
     *     all there (the same stored rows read and hits as keeping every row).
     * So a call of n rows whose paths read s stored rows keeps at most n x (successors + 3) +
     * s / checkpoint rows that are not narrow, where keeping every row would keep s of them
     * (rows_inserted and bytes_inserted count what it copies). Measured on the synthetic
     * index (wide rows; paths of up to 100 and of up to 1,000 rows), against the cache off:
     * no request slower beyond noise, isolated first reads 1.3-2.2x faster than the cache off
     * (a cut at a checkpoint replaces up to a whole path), walks with batch_kmers 1 14-15x;
     * the defaults were chosen from 16/8, 32/16, 64/32 (also 4/8 and 8/32 on UHGG), 16/8
     * being never slower than the cache off on either index.
     */
    static constexpr uint32_t kCheckpoint = 16;
    static constexpr uint32_t kSuccessors = 8;
    static constexpr uint64_t kNarrowBytes = 4096;
    // the rule's parameters (tests vary them; what is kept never changes a row or a cost;
    // narrow_bytes 0: no row is narrow)
    uint32_t checkpoint = kCheckpoint;
    uint32_t successors = kSuccessors;
    uint64_t narrow_bytes = kNarrowBytes;
    // |depth|: the row's distance to its anchor; |from_front|: its position on the path of
    // the requested row that reconstructed it (0: that row); |row|: the row
    bool keeps(uint32_t depth, size_t from_front, const RowT &row) const {
        if (from_front <= successors || depth % checkpoint == 0)
            return true;
        // the row's own buffer first, a lower bound of its copy: a wide row is told without
        // a pass over its columns (row_copy_bytes of a tuple row visits every column)
        if (row.size() * sizeof(row[0]) >= narrow_bytes)
            return false;
        return row_copy_bytes(row) < narrow_bytes;
    }

    // Cache the reconstructed row |full| of |row| at |depth| (with its path aggregates, if
    // known), within the bound; a row that does not fit even into an empty cache is not kept.
    // Admitted before it is copied: its stored size is computed from |full|, the shared bound
    // read, and the older entries evicted, and only then is the row copied — so a row the
    // cache refuses is never copied, and the copy of one it keeps is made within the bound,
    // never beside a cache that is full. Copying first would allocate a refused row whole
    // beside the bound (4.2 MB for a row of 1,048,576 columns in a cache bounded at 1 KiB)
    // with no count stating it.
    void insert(Row row, const RowT &full, uint32_t depth = 0,
                const PathAggregates *path = nullptr) {
        if (!enabled())
            return;
        auto it = gen_[0].find(row);
        if (it != gen_[0].end()) {
            if (path && !it->second.has_path) {
                it.value().has_path = true;
                it.value().path = *path;
            }
            return;
        }
        auto old = gen_[1].find(row);
        if (old != gen_[1].end()) {
            if (path && !old->second.has_path) {
                old.value().has_path = true;
                old.value().path = *path;
            }
            return;
        }
        const uint64_t b = StoredRow<RowT>::bytes_of(full) + kEntryBytes;
        const uint64_t bound = limit();
        // a shared bound may have shrunk since the last insert
        trim(bound);
        if (b > bound / 2) {
            rows_refused++;
            return;
        }
        if (bytes_[0] + b > bound / 2 || bytes() + b > bound) {
            // the current generation becomes the older one, the older one is dropped (the
            // counts move with their tables: a count left behind would let trim() and
            // make_room keep an older generation they no longer see)
            std::swap(gen_[0], gen_[1]);
            std::swap(bytes_[0], bytes_[1]);
            drop(0);
        }
        if (bytes() + b > bound)
            drop(1);
        // the copy, now that there is room for it: the cache holds bytes() + b from here on
        // (the table's growth on the emplace is kEntryBytes's share), which the peak counts
        Entry entry;
        entry.row.store(full);
        assert(entry.row.bytes() + kEntryBytes == b);
        entry.copy_bytes = row_copy_bytes(full);
        entry.bytes = b;
        entry.depth = depth;
        rows_inserted++;
        bytes_inserted += b - kEntryBytes;
        if (path) {
            entry.has_path = true;
            entry.path = *path;
        }
        gen_[0].emplace(row, std::move(entry));
        bytes_[0] += b;
        peak_bytes = std::max(peak_bytes, bytes());
    }
    // evict (the older generation first) until the cache holds at most |max| bytes; the
    // decode calls trim to limit() when they begin, so that a shared bound that shrank
    // without a trim (LabelOracle::make_room) is kept by every call that reads the cache
    void trim(uint64_t max) {
        if (bytes() <= max)
            return;
        drop(1);
        if (bytes_[0] > max)
            drop(0);
    }
    void clear() {
        drop(0);
        drop(1);
    }

    // statistics (physical, timing only): paths cut at a cached row; rows copied into the
    // cache and their bytes (row_copy_bytes)
    uint64_t hits = 0;
    // stored rows (diffs and anchors) the decode calls with the cache read from the matrix
    uint64_t stored_rows_read = 0;
    uint64_t rows_inserted = 0;
    uint64_t bytes_inserted = 0;
    // rows the rule kept that did not fit (more than half the bound): refused before any copy
    uint64_t rows_refused = 0;
    // the most the cache held at once (bytes(), at the end of an insert: rows are copied
    // only after the evictions that make room for them, so no insert holds more)
    uint64_t peak_bytes = 0;

    // A budget-aware read with the cache holds less than without it — its paths stop at
    // cached rows — so a read refused without the cache may fit with it. That makes no
    // difference to where a walk stops (a row is admitted by its whole path's demand, which a
    // read alone never exceeds), but the lookahead under a memory budget reads ahead until a
    // run does not fit, and with the cache it would read rows no level asked for (UHGG 16S,
    // annotate mode, 8 MiB: 75k rows warmed against 52-55k without the cache, 13-27% more
    // instructions). Set (the walker's, for the lookahead under a memory
    // budget), a read also charges, until its stored rows are read, what the decode to the
    // anchors would have held for the rows the cache spared it (IRowDiff::decode_budgeted)
    bool admit_as_uncached = false;

  private:
    using Map = tsl::hopscotch_map<Row, Entry>;
    // a generation dropped with its table
    void drop(size_t i) {
        Map().swap(gen_[i]);
        bytes_[i] = 0;
    }

    Map gen_[2];
    uint64_t bytes_[2] = { 0, 0 };
    uint64_t max_bytes_ = 0;
    std::function<uint64_t()> room_;
};

/**
 * The path caches of one request: one per row type (a binary row-diff annotation caches
 * SetBitPositions, a coordinate one RowTuples, from which both its rows and its tuple rows
 * are made), sharing one bound — only one of them is ever filled, by the annotation's type.
 */
struct RowDiffPathCache {
    RowDiffCache<BinaryMatrix::SetBitPositions> rows;
    RowDiffCache<MultiIntMatrix::RowTuples> tuples;

    void set_bound(uint64_t max_bytes, const std::function<uint64_t()> &room = {}) {
        rows.set_bound(max_bytes, room);
        tuples.set_bound(max_bytes, room);
    }
    bool enabled() const { return rows.enabled(); }
    uint64_t limit() const { return rows.limit(); }
    uint64_t bytes() const { return rows.bytes() + tuples.bytes(); }
    uint64_t hits() const { return rows.hits + tuples.hits; }
    uint64_t stored_rows_read() const { return rows.stored_rows_read + tuples.stored_rows_read; }
    uint64_t rows_inserted() const { return rows.rows_inserted + tuples.rows_inserted; }
    uint64_t bytes_inserted() const { return rows.bytes_inserted + tuples.bytes_inserted; }
    uint64_t peak_bytes() const { return std::max(rows.peak_bytes, tuples.peak_bytes); }
    // RowDiffCache::keeps's parameters, for both
    void set_retention(uint32_t checkpoint, uint32_t successors, uint64_t narrow_bytes) {
        rows.checkpoint = tuples.checkpoint = std::max<uint32_t>(checkpoint, 1);
        rows.successors = tuples.successors = successors;
        rows.narrow_bytes = tuples.narrow_bytes = narrow_bytes;
    }
    // RowDiffCache::admit_as_uncached
    void set_admit_as_uncached(bool on) {
        rows.admit_as_uncached = on;
        tuples.admit_as_uncached = on;
    }
    void trim(uint64_t max) {
        rows.trim(max);
        tuples.trim(max);
    }
    void clear() {
        rows.clear();
        tuples.clear();
    }
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __ROW_DIFF_CACHE_HPP__
