#ifndef __ROW_DIFF_CACHE_HPP__
#define __ROW_DIFF_CACHE_HPP__

#include <cstdint>
#include <functional>
#include <limits>
#include <utility>

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

/**
 * The row-diff path cache (the efficiency pass): reconstructed rows of a row-diff
 * annotation — the rows a decode call returned and every row on their paths, anchors
 * included — kept across the decode calls of one request, so that a later call's path that
 * meets a cached row stops there instead of decoding to its anchor again. On refseq33m a
 * fetch of one warm 23S k-mer cost one whole path decode (66-86 ms) and each further row on
 * the same path 2.1 ms; a walk fetching a few keys per level paid 27-41 ms a row against 4.6
 * ms read in one call. A row's content does not depend on how it was reached, so what a
 * read returns never depends on the cache; only the decoding work does.
 *
 * Bounded in bytes (row_copy_bytes plus kEntryBytes per entry, an estimate of the table's
 * share), in two generations: inserts go to the current one, which becomes the older one
 * when it reaches half the bound (the older one is dropped then); a hit in the older one
 * moves the row to the current one. A walk moves forward, so the rows of the last few
 * levels are the ones its next paths meet. A generation dropped releases its table: a
 * cleared hopscotch map keeps its bucket array, sized for the most entries it ever held, so
 * after many narrow rows a few wide ones held that array beside their bound (74.7 MiB of heap
 * for 63.9 MiB accounted at a 64 MiB bound; review of the efficiency pass). Single-threaded,
 * as the request that owns it.
 */
template <class RowT>
class RowDiffCache {
  public:
    using Row = BinaryMatrix::Row;
    struct Entry {
        RowT row;
        uint64_t bytes = 0;
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
    // Cache the reconstructed row |full| of |row| (with its path aggregates, if known),
    // within the bound; a row that does not fit even into an empty cache is not kept
    void insert(Row row, const RowT &full, const PathAggregates *path = nullptr) {
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
        const uint64_t b = row_copy_bytes(full) + kEntryBytes;
        const uint64_t bound = limit();
        // a shared bound may have shrunk since the last insert
        trim(bound);
        if (b > bound / 2)
            return;
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
        Entry entry;
        entry.row = full;
        entry.bytes = b;
        if (path) {
            entry.has_path = true;
            entry.path = *path;
        }
        gen_[0].emplace(row, std::move(entry));
        bytes_[0] += b;
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

    // statistics (the physical work saved): paths cut at a cached row
    uint64_t hits = 0;

    // A budget-aware read with the cache holds less than without it — its paths stop at
    // cached rows — so a read refused without the cache may fit with it. That makes no
    // difference to where a walk stops (a row is admitted by its whole path's demand, which a
    // read alone never exceeds), but the lookahead under a memory budget reads ahead until a
    // run does not fit, and with the cache it read rows no level asked for (UHGG 16S, annotate
    // mode, 8 MiB: 75k rows warmed against 52-55k without the cache, 13-27% more instructions;
    // review of the efficiency pass). Set (the walker's, for the lookahead under a memory
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
