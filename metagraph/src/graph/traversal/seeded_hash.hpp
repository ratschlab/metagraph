#ifndef __SEEDED_HASH_HPP__
#define __SEEDED_HASH_HPP__

#include <cstddef>
#include <cstdint>

#include "common/hashers/seeded_hash.hpp"


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
 * process, the hash gives a key set chosen without the seed the buckets of random keys, at which
 * the budget's per-entry model is measured (hopscotch_entry_bytes). Keys chosen for a known seed
 * can still double a table early, so the pass's tables are charged as allocated beyond the model
 * (charged_table.hpp, GPT review 7, finding 7). Nothing an answer states depends on the seed:
 * these tables are looked up, never iterated for an answer; a memory stop's k-mer can, where a
 * table's allocations pass the model (small tables, chosen keys).
 */

// The hashes themselves (common/hashers/seeded_hash.hpp: the row-diff path cache uses them too)
using common::process_hash_seed;
using common::mix64;
using common::SeededHash;
using common::SeededPairHash;

// The budget's model per entry of a tsl::hopscotch_map (seeded) whose value_type has
// |value_bytes| bytes: a bucket is the value and its neighbourhood bitmap (8 B); at a growth the
// table holds its buckets and twice as many new ones, 3 buckets per bucket, and a table is taken
// to double at 80% load or above, so an entry is charged 3 buckets / 0.8. Measured with random
// keys (the seeded hash's case; a full neighbourhood doubles a table before its 90% load
// factor): the lowest load at a growth was 0.833 over 3,600 growths of tables of 2^11 to 2^22
// buckets (median 0.90), 0.822 over 280 growths of 2^20 to 2^26 buckets (median 0.87). A model,
// not a bound whatever the keys: tsl's tail of 61 buckets passes it below a few hundred entries,
// and keys chosen for a known seed double a table down to 10% load (up to 8 times the model at a
// growth, V7-B); a table charged through a TableAccount is charged that excess as it is allocated
constexpr uint64_t hopscotch_entry_bytes(size_t value_bytes) {
    return (3 * (value_bytes + sizeof(uint64_t)) * 5 + 3) / 4;
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __SEEDED_HASH_HPP__
