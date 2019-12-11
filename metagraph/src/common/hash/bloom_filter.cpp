#include "bloom_filter.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "utils/serialization.hpp"

constexpr uint32_t BLOCK_MASK = 0b111111111;
constexpr uint32_t SHIFT = 9;

BloomFilter::BloomFilter(size_t filter_size, uint32_t num_hash_functions)
      : filter_(filter_size ? (((filter_size + BLOCK_MASK) >> SHIFT) << SHIFT) : BLOCK_MASK + 1),
        num_hash_functions_(num_hash_functions) {
    assert(filter_.size() >= filter_size);
    assert(filter_.size() > BLOCK_MASK);
    assert(!(filter_.size() & BLOCK_MASK));
}

BloomFilter::BloomFilter(size_t filter_size,
                         size_t expected_num_elements,
                         uint32_t max_num_hash_functions)
      : BloomFilter(filter_size,
                    std::min(max_num_hash_functions,
                             optim_h(filter_size, expected_num_elements))) {}

// TODO: rewrite some of these with SIMD instructions

void BloomFilter::insert(uint64_t hash) {
    const auto size = filter_.size();

    // use the 64-bit hash to select a 512-bit block
    size_t offset = ((hash % size) >> SHIFT) << SHIFT;

    // split 64-bit hash into two 32-bit hashes
    uint32_t h1 = hash & 0xFFFFFFFF;
    uint32_t h2 = hash >> 32;

    /*
     * Use two 32-bit hashes to generate num_hash_functions hashes
     * Kirsch, A., & Mitzenmacher, M. (2006, September).
     * Less hashing, same performance: building a better bloom filter.
     * In European Symposium on Algorithms (pp. 456-467). Springer, Berlin, Heidelberg.
     */
    for (uint32_t i = 0; i < num_hash_functions_; ++i) {
        // set bits within block
        filter_[offset + ((h1 + i * h2) & BLOCK_MASK)] = true;
    }

    assert(check(hash));
}

bool BloomFilter::check(uint64_t hash) const {
    const auto size = filter_.size();

    // use the 64-bit hash to select a 512-bit block
    size_t offset = ((hash % size) >> SHIFT) << SHIFT;

    // split 64-bit hash into two 32-bit hashes
    uint32_t h1 = hash & 0xFFFFFFFF;
    uint32_t h2 = hash >> 32;

    bool found = true;
    for (uint32_t i = 0; i < num_hash_functions_; ++i) {
        // check bits within block
        found &= filter_[offset + ((h1 + i * h2) & BLOCK_MASK)];
    }

    return found;
}


// constants for batch insert and check
#ifdef __AVX2__
const __m128i zeros = _mm_set1_epi32(0);

// used to restrict indices to the size of a block
const __m128i blockmask = _mm_set1_epi32(BLOCK_MASK);
#endif

void BloomFilter::batch_insert(const uint64_t hashes[], size_t len) {
    const uint64_t *hs = hashes;
    const uint64_t *end = hashes + len;

#ifdef __AVX2__
    // compute Bloom filter hashes in batches of 4
    const auto size = filter_.size();

    uint64_t indices[4] __attribute__ ((aligned (32)));
    uint32_t *hh;

    __m256i offsets;
    __m128i lhashes, mult;

    for (; hs + 4 <= end; hs += 4) {
        // compute offsets
        // offset = ((hash % size) >> SHIFT) << SHIFT;
        offsets = _mm256_setr_epi64x(hs[0] % size,
                                     hs[1] % size,
                                     hs[2] % size,
                                     hs[3] % size);
        offsets = _mm256_srli_epi64(offsets, SHIFT);
        offsets = _mm256_slli_epi64(offsets, SHIFT);

        _mm256_store_si256((__m256i*)indices, offsets);

        // clean up after AVX2 instructions
        _mm256_zeroupper();

        for (size_t j = 0; j < 4; ++j) {
            uint32_t k = 0;

            // hash functions in pairs
            for (; k + 2 <= num_hash_functions_; k += 2) {
                // load hash
                lhashes = _mm_set1_epi64x(hs[j]);

                // compute hashes
                // hash = (h1 + k * h2) & BLOCK_MASK
                mult = _mm_setr_epi32(1, k, 1, k + 1);
                lhashes = _mm_mullo_epi32(lhashes, mult);
                lhashes = _mm_hadd_epi32(lhashes, zeros);
                lhashes = _mm_and_si128(lhashes, blockmask);

                // add to filter
                hh = (uint32_t*)&lhashes;
                filter_[indices[j] + hh[0]] = true;
                filter_[indices[j] + hh[1]] = true;
            }

            // if num_hash_functions is odd, add the last one
            if (num_hash_functions_ & 1) {
                hh = (uint32_t*)&hs[j];
                filter_[indices[j]
                            + ((hh[0] + (num_hash_functions_ - 1) * hh[1])
                                & BLOCK_MASK)] = true;
            }
        }

        assert(std::all_of(hs, hs + 4, [&](uint64_t hash) { return check(hash); }));
    }

#endif

    for (; hs < end; ++hs) {
        insert(*hs);
    }
}

sdsl::bit_vector BloomFilter
::batch_check(const std::vector<std::pair<uint64_t, size_t>> &hash_index,
             size_t length) const {
    const auto num_elements = hash_index.size();
    sdsl::bit_vector presence(length, false);
    size_t i = 0;

#ifdef __AVX2__
    // compute Bloom filter hashes in batches of 4

    const auto size = filter_.size();

    uint64_t hs[4] __attribute__ ((aligned (32)));
    uint64_t indices[4] __attribute__ ((aligned (32)));
    uint32_t *hh;

    __m256i offsets;
    __m128i hashes, mult;

    for (; i + 4 <= num_elements; i += 4) {
        // copy hashes
        hs[0] = hash_index[i].first;
        hs[1] = hash_index[i + 1].first;
        hs[2] = hash_index[i + 2].first;
        hs[3] = hash_index[i + 3].first;

        // compute offsets
        // offset = ((hash % size) >> SHIFT) << SHIFT;
        offsets = _mm256_setr_epi64x(hs[0] % size,
                                     hs[1] % size,
                                     hs[2] % size,
                                     hs[3] % size);
        offsets = _mm256_srli_epi64(offsets, SHIFT);
        offsets = _mm256_slli_epi64(offsets, SHIFT);

        _mm256_store_si256((__m256i*)indices, offsets);

        // clean up after AVX2 instructions
        _mm256_zeroupper();

        for (size_t j = 0; j < 4; ++j) {
            bool found = true;
            uint32_t k = 0;

            // hash functions in pairs
            for (; found && k + 2 <= num_hash_functions_; k += 2) {
                // load hash
                hashes = _mm_set1_epi64x(hs[j]);

                // compute hashes
                // hash = (h1 + k * h2) & BLOCK_MASK
                mult = _mm_setr_epi32(1, k, 1, k + 1);
                hashes = _mm_mullo_epi32(hashes, mult);
                hashes = _mm_hadd_epi32(hashes, zeros);
                hashes = _mm_and_si128(hashes, blockmask);

                // check hashes
                hh = (uint32_t*)&hashes;
                found &= filter_[indices[j] + hh[0]] & filter_[indices[j] + hh[1]];
            }

            // if num_hash_functions is odd, check the last one
            if (found && (num_hash_functions_ & 1)) {
                hh = (uint32_t*)&hs[j];
                found &= filter_[indices[j]
                                    + ((hh[0] + (num_hash_functions_ - 1) * hh[1])
                                        & BLOCK_MASK)];
            }

            if (found) {
                presence[hash_index[i + j].second] = true;
                assert(check(hs[j]));
            }
        }
    }

#endif

    for (; i < num_elements; ++i) {
        presence[hash_index[i].second] = check(hash_index[i].first);
    }

    return presence;
}

void BloomFilter::serialize(std::ostream &out) const {
    filter_.serialize(out);
    serialize_number(out, num_hash_functions_);
}

bool BloomFilter::load(std::istream &in) {
    try {
        filter_.load(in);
        num_hash_functions_ = load_number(in);
        return true;
    } catch (...) {
        return false;
    }
}
