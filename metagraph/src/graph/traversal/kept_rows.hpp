#ifndef __KEPT_ROWS_HPP__
#define __KEPT_ROWS_HPP__

#include <cstdint>
#include <functional>
#include <utility>

#include <tsl/hopscotch_map.h>

#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "graph/traversal/charged_table.hpp"
#include "graph/traversal/seeded_hash.hpp"


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The rows a discovery's pass keeps for a later occurrence of their k-mer
 * (SPEC-labeled-traversal-core.md §4.6), with what |budget| holds for them: each row's own
 * charge (the decoder's, passed in), released when the row is dropped, and the table's: the
 * model, kEntryBytes per entry charged before an entry is added (hopscotch_entry_bytes: a table
 * doubling at 80% load or above, which its seeded hash makes the case for keys the seed does
 * not favour, T3 review R2 F1: the entry, the empty buckets and the table before the growth),
 * and, through the table's allocator (TableAccount, ChargedAllocator), whatever its allocations
 * — the buckets with tsl's tail, at a growth the old and the new ones, the overflow list's
 * nodes — would hold beyond the model, charged before they are made (GPT review 7, finding 7:
 * the model alone was passed by a table of a few hundred entries, whose tail it does not
 * cover, and by keys chosen for a known seed, up to 8 times). A growth that does not fit is
 * refused before anything is allocated, and the row is then not kept. The table never shrinks
 * while it lives (erasing an entry keeps its bucket), so the entries' charge is a high-water
 * count: an entry freed by a last use is charged again only past the most entries held at
 * once, and the whole count is released only by drop(), which frees the table. |Hash| is
 * SeededHash (a test passes another to choose the buckets).
 */
template <class RowT, class Key = uint64_t, class Hash = SeededHash>
class KeptRows {
  public:
    struct Kept {
        RowT row;
        uint64_t bytes = 0;     // as row_copy_bytes counts them (the pass's kept_max)
        uint64_t charged = 0;   // the row's own charge the budget holds
    };
    using Alloc = ChargedAllocator<std::pair<Key, Kept>>;
    using Table = tsl::hopscotch_map<Key, Kept, Hash, std::equal_to<Key>, Alloc>;

    // 180 B for a row of set bits (a value of 40 B), 210 B for a row of tuples (48 B)
    static constexpr uint64_t kEntryBytes
            = hopscotch_entry_bytes(sizeof(typename Table::value_type));

    explicit KeptRows(annot::matrix::DecodeBudget &budget, const Hash &hash = Hash())
          : budget_(budget), account_(budget, /* refuse */ true), hash_(hash),
            table_(empty_table()) {}
    KeptRows(const KeptRows &) = delete;
    KeptRows& operator=(const KeptRows &) = delete;
    // what is still held is released (the pass releases it by drop() before it ends)
    ~KeptRows() { drop(); }

    bool contains(Key key) const { return table_.count(key); }
    typename Table::iterator find(Key key) { return table_.find(key); }
    typename Table::iterator end() { return table_.end(); }

    // Keep |row| (|bytes| as row_copy_bytes counts them, |charged| the decoder's charge the
    // budget holds for it, taken over): false, with nothing kept and |row| untouched, when a new
    // entry's charge or the table's growth does not fit (a new entry's charge, once made, stays
    // with the table's until drop())
    bool keep(Key key, RowT &row, uint64_t bytes, uint64_t charged) {
        if (table_.size() == entries_charged_) {
            if (!account_.cover(kEntryBytes))
                return false;
            ++entries_charged_;
        }
        try {
            // the entry first, empty (a growth refused leaves the table and the row as they
            // were), then the row moved into it (no allocation). Two statements: in one
            // assignment the row would be moved before the entry is made
            auto it = table_.try_emplace(key).first;
            it.value() = Kept { std::move(row), bytes, charged };
        } catch (const TableGrowthRefused &) {
            return false;
        }
        bytes_ += bytes;
        return true;
    }
    bool keep(Key key, RowT &&row, uint64_t bytes, uint64_t charged) {
        return keep(key, row, bytes, charged);
    }

    // The last use of a kept row: its own charge released; its entry's stays (its bucket does)
    void last_use(typename Table::iterator it) {
        bytes_ -= it->second.bytes;
        budget_.release(it->second.charged);
        table_.erase(it);
    }

    // Every row given up and the table freed: everything charged for them released (the
    // table's allocations first, then the model)
    void drop() {
        for (const auto &[key, k] : table_) {
            budget_.release(k.charged);
        }
        empty_table().swap(table_);
        account_.uncover(entries_charged_ * kEntryBytes);
        entries_charged_ = 0;
        bytes_ = 0;
    }

    // whether anything is held for the rows kept (rows, or the entries' charge)
    bool holds() const { return !table_.empty() || entries_charged_; }
    bool empty() const { return table_.empty(); }
    size_t size() const { return table_.size(); }
    uint64_t bytes() const { return bytes_; }
    size_t entries_charged() const { return entries_charged_; }
    size_t bucket_count() const { return table_.bucket_count(); }
    size_t overflow_size() const { return table_.overflow_size(); }
    // the table's account: what its allocations hold, and what the budget holds for it
    const TableAccount& account() const { return account_; }

  private:
    annot::matrix::DecodeBudget &budget_;
    // (before the table: destroyed after it, so that its buckets are freed to it)
    TableAccount account_;
    Hash hash_;
    Table table_;

    // an empty table charged to |account_| (no allocation until its first entry)
    Table empty_table() { return Table(0, hash_, std::equal_to<Key>(), Alloc(&account_)); }
    uint64_t bytes_ = 0;
    size_t entries_charged_ = 0;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __KEPT_ROWS_HPP__
