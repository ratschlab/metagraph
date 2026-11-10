#ifndef __SEEDED_HASH_HPP__
#define __SEEDED_HASH_HPP__

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The hashes of a /resolve discovery's hash tables whose keys a request chooses: node indices
 * (the pass's setup table, its batch table and its kept rows, the selection's distinct keys)
 * and records (the record index). A tsl::hopscotch_map takes a key's bucket from the hash's low
 * bits and doubles, down to 10% load, when no empty bucket can be brought into the key's
 * neighbourhood. std::hash of an integer is the identity, so a query whose k-mers have chosen
 * node indices could double a table at each 63 keys sharing low bits, to about 20 times its
 * entries (T3 review R2, F1: 1,344 MiB for 1 M kept rows charged 183 MiB). Seeded once per
 * process, the hash gives any key set the buckets of random keys, at which the budget's
 * per-entry charges are measured (hopscotch_entry_bytes). Nothing an answer states depends on
 * the seed: these tables are looked up, never iterated for an answer.
 */

// The process's seed, drawn once (std::random_device) when it is first used
inline uint64_t process_hash_seed() {
    static const uint64_t seed = []() {
        std::random_device device;
        return (uint64_t(device()) << 32) ^ uint64_t(device());
    }();
    return seed;
}

// MurmurHash3's 64-bit finalizer: a bijection whose every output bit depends on every input bit.
// Cheaper mixes spread keys worse: over random and consecutive keys a 128-bit multiply fold let a
// table double at 42% load, a one-multiply xorshift at 63%, this one at 84% at the lowest
// (H2 t3fix, grow_mix.cpp); its cost is about a nanosecond a lookup
inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xFF51AFD7ED558CCDULL;
    x ^= x >> 33;
    x *= 0xC4CEB9FE1A85EC53ULL;
    x ^= x >> 33;
    return x;
}

// An integer key (a node index), seeded
struct SeededHash {
    uint64_t seed = process_hash_seed();

    size_t operator()(uint64_t key) const { return static_cast<size_t>(mix64(key ^ seed)); }
};

// A pair of integer keys (a record: its column and its sequence id), seeded: the first rotated
// into the high half beside the second (one to one while both are below 2^32), then one mix
struct SeededPairHash {
    uint64_t seed = process_hash_seed();

    template <class A, class B>
    size_t operator()(const std::pair<A, B> &key) const {
        const uint64_t first = static_cast<uint64_t>(key.first);
        const uint64_t second = static_cast<uint64_t>(key.second);
        return static_cast<size_t>(mix64(((first << 32) | (first >> 32)) ^ second ^ seed));
    }
};

// The budget's charge per entry of a tsl::hopscotch_map (seeded) whose value_type has
// |value_bytes| bytes: a bucket is the value and its neighbourhood bitmap (8 B); at a growth the
// table holds its buckets and twice as many new ones, 3 buckets per bucket, and a table is taken
// to double at 80% load or above, so an entry is charged 3 buckets / 0.8. Measured with random
// keys (the seeded hash's case; a full neighbourhood doubles a table before its 90% load
// factor): the lowest load at a growth was 0.833 over 3,600 growths of tables of 2^11 to 2^22
// buckets (median 0.90), 0.822 over 280 growths of 2^20 to 2^26 buckets (median 0.87)
constexpr uint64_t hopscotch_entry_bytes(size_t value_bytes) {
    return (3 * (value_bytes + sizeof(uint64_t)) * 5 + 3) / 4;
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __SEEDED_HASH_HPP__
