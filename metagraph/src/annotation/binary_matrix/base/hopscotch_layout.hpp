#ifndef __HOPSCOTCH_LAYOUT_HPP__
#define __HOPSCOTCH_LAYOUT_HPP__

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "annotation/binary_matrix/base/decode_budget.hpp"


namespace mtg {
namespace annot {
namespace matrix {

// The capacity of a tsl::hopscotch_map or hopscotch_set, seeded (common::SeededHash) or not,
// counted as the budget counts every buffer (alloc_bytes): the row-diff path cache's tables
// (row_diff_cache.hpp). A /resolve discovery's tables are charged through their allocator
// instead (graph/traversal/charged_table.hpp), which counts the same buffers at the same sizes

/**
 * The heap of a tsl::hopscotch_map or hopscotch_set as it is, whatever its keys and its hash
 * (GPT review 7, findings 7 and 8: a table can double at 10% load, so a per-entry model is a
 * bound only for keys spread as random keys are). A table holds its bucket array — once it has
 * buckets, bucket_count() + NeighborhoodSize - 1 of them (the last bucket's neighbourhood is
 * padding: a table of one entry holds 63 buckets) — and its overflow list's nodes (a value and
 * two pointers each). An insert doubles the bucket count (2 after none), the new array held
 * beside the old one while the entries move; an insert that doubles a table more than once
 * (keys colliding in the hash after a doubling, which the seeded hash makes improbable) holds
 * more. Each buffer is counted at its jemalloc size class (alloc_bytes, as every buffer of
 * the budget's model): what the heap holds with jemalloc 5.3 and earlier (four classes per
 * doubling), an upper bound of it under jemalloc 5.4, which gives large buffers page-granular
 * sizes on macOS (16 KiB pages: a 2,105,984 B buffer takes 2,113,536 B, counted 2,621,440)
 * and has an 8 B class (counted 16); never below the heap (R4 review, S4).
 */
template <class Table>
struct HopscotchLayout;

template <class Value, unsigned int NeighborhoodSize, bool StoreHash, class GrowthPolicy>
struct HopscotchLayoutOf {
    static_assert(std::is_same<GrowthPolicy, tsl::hh::power_of_two_growth_policy<2>>::value,
                  "the layout knows the doubling growth policy only");
    // a bucket: its value and its neighbourhood bitmap (and its stored hash, if any)
    static constexpr uint64_t kBucketBytes = sizeof(
            tsl::detail_hopscotch_hash::hopscotch_bucket<Value, NeighborhoodSize, StoreHash>);
    static constexpr uint64_t kPaddingBuckets = NeighborhoodSize - 1;
    // an overflow node of std::list: two pointers, then the value at its alignment
    static uint64_t node_bytes() {
        return annot::matrix::alloc_bytes(
                (2 * sizeof(void *) + alignof(Value) - 1) / alignof(Value) * alignof(Value)
                    + sizeof(Value));
    }

    // the bucket array of a table of |bucket_count| buckets (one buffer) at its size class:
    // what the heap holds, or an upper bound of it (jemalloc 5.4's page-granular large sizes)
    static uint64_t array_bytes(uint64_t bucket_count) {
        return bucket_count
                ? annot::matrix::alloc_bytes((bucket_count + kPaddingBuckets) * kBucketBytes)
                : 0;
    }
    // the bucket count an insert that grows a table of |bucket_count| buckets doubles it to
    static constexpr uint64_t next_bucket_count(uint64_t bucket_count) {
        return bucket_count ? 2 * bucket_count : 2;
    }
};

template <class Key, class T, class Hash, class KeyEqual, class Allocator,
          unsigned int NeighborhoodSize, bool StoreHash, class GrowthPolicy>
struct HopscotchLayout<tsl::hopscotch_map<Key, T, Hash, KeyEqual, Allocator, NeighborhoodSize,
                                          StoreHash, GrowthPolicy>>
      : HopscotchLayoutOf<std::pair<Key, T>, NeighborhoodSize, StoreHash, GrowthPolicy> {};

template <class Key, class Hash, class KeyEqual, class Allocator,
          unsigned int NeighborhoodSize, bool StoreHash, class GrowthPolicy>
struct HopscotchLayout<tsl::hopscotch_set<Key, Hash, KeyEqual, Allocator, NeighborhoodSize,
                                          StoreHash, GrowthPolicy>>
      : HopscotchLayoutOf<Key, NeighborhoodSize, StoreHash, GrowthPolicy> {};

// What |table| holds: its bucket array and its overflow nodes
template <class Table>
uint64_t hopscotch_held_bytes(const Table &table) {
    using Layout = HopscotchLayout<Table>;
    return Layout::array_bytes(table.bucket_count())
            + table.overflow_size() * Layout::node_bytes();
}

// What one insert into |table| adds at most to what it holds while it runs, when it doubles
// the table once: the new bucket array beside the old one, or one overflow node
template <class Table>
uint64_t hopscotch_growth_bytes(const Table &table) {
    using Layout = HopscotchLayout<Table>;
    return Layout::array_bytes(Layout::next_bucket_count(table.bucket_count()))
            + Layout::node_bytes();
}

/**
 * The load (entries in buckets per bucket) from which an insert is taken to be able to double
 * a seeded table. tsl doubles a table at its own threshold — an insert into a table whose
 * buckets hold size_type(bucket_count x max_load_factor) entries (0.9) rehashes first — and
 * before it, from 10% load, when no empty bucket can be brought into the key's neighbourhood.
 * Which insert finds its neighbourhood full depends on the bucket array, which tsl does not
 * expose; with the seeded hash it happens near the threshold only: over 133k growths of seeded
 * tables of 2^1 to 2^22 buckets (2,000 seeds x 4 key patterns to 2^16, 60 to 2^22; random,
 * consecutive, bit-reversed and i << 20 keys) the lowest load at a growth was 0.825, none
 * below 0.8 (T4-L4 fix round, s1/growth.cpp), and 0.822 over 280 growths of 2^20 to 2^26
 * buckets (graph/traversal/seeded_hash.hpp, hopscotch_entry_bytes). A growth below this load
 * needs keys that collide in the hash, which the seeded hash keeps a client from choosing.
 */
constexpr float kHopscotchGrowthLoad = 0.8f;

// Whether one insert into |table| can double it: at tsl's threshold, or from
// kHopscotchGrowthLoad (a full neighbourhood of a seeded table). An empty table (no buckets)
// gets its first array at its first insert: it can
template <class Table>
bool hopscotch_insert_can_grow(const Table &table) {
    const size_t buckets = table.bucket_count();
    const size_t in_buckets = table.size() - table.overflow_size();
    return in_buckets >= size_t(float(buckets) * table.max_load_factor())
            || in_buckets >= size_t(float(buckets) * kHopscotchGrowthLoad);
}

// What one insert into |table| adds at most to what it holds while it runs: its doubling
// (hopscotch_growth_bytes) when it can grow the table (hopscotch_insert_can_grow), else one
// overflow node (an entry put into a bucket adds nothing). Ids that double a table below
// kHopscotchGrowthLoad add more: the caller reads what the table holds after the insert
template <class Table>
uint64_t hopscotch_insert_bytes(const Table &table) {
    return hopscotch_insert_can_grow(table) ? hopscotch_growth_bytes(table)
                                            : HopscotchLayout<Table>::node_bytes();
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __HOPSCOTCH_LAYOUT_HPP__
