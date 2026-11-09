#ifndef __KEPT_ROWS_HPP__
#define __KEPT_ROWS_HPP__

#include <cstdint>
#include <utility>

#include <tsl/hopscotch_map.h>

#include "annotation/binary_matrix/base/decode_budget.hpp"


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The rows a discovery's pass keeps for a later occurrence of their k-mer
 * (SPEC-labeled-traversal-core.md §4.6), with what |budget| holds for them: each row's own
 * charge (the decoder's, passed in), released when the row is dropped, and kEntryBytes per
 * entry of the table, charged before the table can grow (a hash table at most 90% full that
 * doubles: the entry, the empty buckets and the table before the growth). The table never
 * shrinks while it lives (erasing an entry keeps its bucket), so the entries' charge is a
 * high-water count: an entry freed by a last use is charged again only past the most entries
 * held at once, and the whole count is released only by drop(), which frees the table.
 */
template <class RowT, class Key = uint64_t>
class KeptRows {
  public:
    static constexpr uint64_t kEntryBytes = 192;

    struct Kept {
        RowT row;
        uint64_t bytes;     // as row_copy_bytes counts them (the pass's kept_max)
        uint64_t charged;   // the row's own charge the budget holds
    };
    using Table = tsl::hopscotch_map<Key, Kept>;

    explicit KeptRows(annot::matrix::DecodeBudget &budget) : budget_(budget) {}
    KeptRows(const KeptRows &) = delete;
    KeptRows& operator=(const KeptRows &) = delete;
    // what is still held is released (the pass releases it by drop() before it ends)
    ~KeptRows() { drop(); }

    bool contains(Key key) const { return table_.count(key); }
    typename Table::iterator find(Key key) { return table_.find(key); }
    typename Table::iterator end() { return table_.end(); }

    // Keep |row| (|bytes| as row_copy_bytes counts them, |charged| the decoder's charge the
    // budget holds for it, taken over): false, with nothing kept and nothing charged, when a
    // new entry's charge does not fit
    bool keep(Key key, RowT &&row, uint64_t bytes, uint64_t charged) {
        if (table_.size() == entries_charged_) {
            if (!budget_.charge(kEntryBytes))
                return false;
            ++entries_charged_;
        }
        bytes_ += bytes;
        table_.emplace(key, Kept { std::move(row), bytes, charged });
        return true;
    }

    // The last use of a kept row: its own charge released; its entry's stays (its bucket does)
    void last_use(typename Table::iterator it) {
        bytes_ -= it->second.bytes;
        budget_.release(it->second.charged);
        table_.erase(it);
    }

    // Every row given up and the table freed: everything charged for them released
    void drop() {
        for (const auto &[key, k] : table_) {
            budget_.release(k.charged);
        }
        budget_.release(entries_charged_ * kEntryBytes);
        entries_charged_ = 0;
        bytes_ = 0;
        Table().swap(table_);
    }

    // whether anything is held for the rows kept (rows, or the entries' charge)
    bool holds() const { return !table_.empty() || entries_charged_; }
    bool empty() const { return table_.empty(); }
    size_t size() const { return table_.size(); }
    uint64_t bytes() const { return bytes_; }
    size_t entries_charged() const { return entries_charged_; }
    size_t bucket_count() const { return table_.bucket_count(); }

  private:
    annot::matrix::DecodeBudget &budget_;
    Table table_;
    uint64_t bytes_ = 0;
    size_t entries_charged_ = 0;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __KEPT_ROWS_HPP__
