#ifndef __CHARGED_TABLE_HPP__
#define __CHARGED_TABLE_HPP__

#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <type_traits>

#include "annotation/binary_matrix/base/decode_budget.hpp"


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The account of one hash table's memory against a DecodeBudget (SPEC-labeled-traversal-core.md
 * §4.6, GPT review 7, finding 7). Its owner charges the model of the table (hopscotch_entry_bytes
 * an entry: what random keys make a seeded table hold at its largest) through cover() (or
 * adopt(), for model bytes it has charged itself); the
 * table's allocator (ChargedAllocator) reports every allocation BEFORE it is made and every
 * deallocation after it, and whatever the table's allocations would hold beyond the model is
 * charged then: the budget holds max(the model, what the table holds) at every moment, its growth
 * included (tsl allocates the new buckets while the old ones are held, and keeps a tail of
 * NeighborhoodSize - 1 buckets past them), whatever the keys. A tsl::hopscotch_map doubles down to
 * 10% load when a key's neighbourhood is full, which keys chosen for a known hash seed reach:
 * up to 8 times the model at a growth (V7-B), all of it charged here.
 *
 * |refuse| true: an allocation whose excess does not fit throws TableGrowthRefused before
 * anything is allocated (tsl leaves the table as it was: it builds a grown table aside, and an
 * entry that goes to its overflow list is not inserted), for an owner that can refuse the insert.
 * false: the excess is charged past the budget's maximum (the allocation is made), for an insert
 * that cannot be refused where it happens; the budget then refuses every charge until enough is
 * released, so the owner's next check stops. With keys the seed does not favour the model covers
 * every growth (the lowest load at one, 0.82 measured, is above the model's 0.8) and nothing is
 * charged beyond it but tsl's tail where a table is small
 */
class TableGrowthRefused : public std::bad_alloc {
  public:
    const char* what() const noexcept override {
        return "a hash table's growth does not fit its memory account";
    }
};

class TableAccount {
  public:
    TableAccount(annot::matrix::DecodeBudget &budget, bool refuse)
          : budget_(budget), refuse_(refuse) {}
    TableAccount(const TableAccount &) = delete;
    TableAccount& operator=(const TableAccount &) = delete;
    // what is still charged is released (the owner frees its table first)
    ~TableAccount() {
        budget_.release(covered_ + extra_);
    }

    // Charge |bytes| of the table's model: false, with nothing charged, when they do not fit
    bool cover(uint64_t bytes) {
        if (!budget_.charge(bytes))
            return false;
        covered_ += bytes;
        rebalance();
        return true;
    }
    // Take over |bytes| of the model the owner has charged to the budget already
    void adopt(uint64_t bytes) {
        covered_ += bytes;
        rebalance();
    }
    // Charge |bytes| of the model past the budget's maximum if need be (an accumulation done)
    void cover_beyond(uint64_t bytes) {
        budget_.charge_beyond(bytes);
        covered_ += bytes;
        rebalance();
    }
    // Release |bytes| of the model (once the entries they priced are gone; the excess of the
    // table's allocations over what is left is charged)
    void uncover(uint64_t bytes) {
        covered_ -= bytes;
        budget_.release(bytes);
        rebalance();
    }

    // ChargedAllocator's: before an allocation of |bytes|, after a deallocation
    void allocating(uint64_t bytes) {
        const uint64_t need = held_ + bytes > covered_ ? held_ + bytes - covered_ : 0;
        if (need > extra_) {
            const uint64_t more = need - extra_;
            if (refuse_) {
                if (!budget_.charge(more))
                    throw TableGrowthRefused();
            } else {
                budget_.charge_beyond(more);
            }
            extra_ += more;
        }
        held_ += bytes;
    }
    void deallocated(uint64_t bytes) noexcept {
        held_ -= bytes;
        rebalance();
    }

    // what the table's allocations hold now (alloc_bytes of each)
    uint64_t held() const { return held_; }
    // the model charged by the owner
    uint64_t covered() const { return covered_; }
    // charged beyond the model: what the allocations hold above it
    uint64_t extra() const { return extra_; }
    // what the budget holds for the table: max(the model, the allocations)
    uint64_t charged() const { return covered_ + extra_; }

  private:
    annot::matrix::DecodeBudget &budget_;
    const bool refuse_;
    uint64_t covered_ = 0;
    uint64_t held_ = 0;
    uint64_t extra_ = 0;

    // extra_ = max(0, held_ - covered_): released where the model covers more, charged (past
    // the maximum if need be: the memory is held already) where it covers less
    void rebalance() noexcept {
        const uint64_t need = held_ > covered_ ? held_ - covered_ : 0;
        if (extra_ > need) {
            budget_.release(extra_ - need);
        } else if (extra_ < need) {
            budget_.charge_beyond(need - extra_);
        }
        extra_ = need;
    }
};

// The allocator of a table whose memory a TableAccount charges (none: a plain allocator).
// Propagated with its table (move, copy, swap), so that a table and the empty one it is swapped
// with are charged to the same account
template <class T>
class ChargedAllocator {
  public:
    using value_type = T;
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::false_type;

    ChargedAllocator() noexcept = default;
    explicit ChargedAllocator(TableAccount *account) noexcept : account_(account) {}
    template <class U>
    ChargedAllocator(const ChargedAllocator<U> &other) noexcept : account_(other.account()) {}

    T* allocate(std::size_t n) {
        const uint64_t bytes = annot::matrix::alloc_bytes(n * sizeof(T));
        if (account_)
            account_->allocating(bytes);
        try {
            return std::allocator<T>().allocate(n);
        } catch (...) {
            if (account_)
                account_->deallocated(bytes);
            throw;
        }
    }
    void deallocate(T *p, std::size_t n) noexcept {
        std::allocator<T>().deallocate(p, n);
        if (account_)
            account_->deallocated(annot::matrix::alloc_bytes(n * sizeof(T)));
    }

    TableAccount* account() const noexcept { return account_; }

    template <class U>
    bool operator==(const ChargedAllocator<U> &other) const noexcept {
        return account_ == other.account();
    }
    template <class U>
    bool operator!=(const ChargedAllocator<U> &other) const noexcept {
        return account_ != other.account();
    }

  private:
    TableAccount *account_ = nullptr;
};

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __CHARGED_TABLE_HPP__
