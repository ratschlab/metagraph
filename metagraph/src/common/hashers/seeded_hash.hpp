#ifndef __COMMON_SEEDED_HASH_HPP__
#define __COMMON_SEEDED_HASH_HPP__

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>


namespace mtg {
namespace common {

/**
 * Hashes for the hash tables whose keys a request chooses, seeded once per process: a
 * /resolve discovery's tables (graph/traversal/seeded_hash.hpp states why and what is charged
 * for them) and the row-diff path cache (row_diff_cache.hpp). A tsl::hopscotch_map takes a
 * key's bucket from the hash's low bits; std::hash of an integer is the identity, so chosen
 * keys could double a table far below its load factor. Seeded, a key set chosen without the
 * seed gets the buckets of random keys. Nothing an answer states depends on the seed.
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

} // namespace common
} // namespace mtg

#endif // __COMMON_SEEDED_HASH_HPP__
