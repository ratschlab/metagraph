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
#ifdef __AVX__
const __m256i zeros = _mm256_set1_epi32(0);

// used to restrict indices to the size of a block
const __m256i blockmask = _mm256_set1_epi32(BLOCK_MASK);

// permute 32-bit hashes into correct order after _mm256_hadd_epi32
const __m256i permute = _mm256_setr_epi32(0, 2, 1, 3, 4, 6, 5, 7);
#endif

void BloomFilter::batch_insert(const uint64_t hashes[], size_t len) {
    size_t j = 0;

#ifdef __AVX2__
    // compute Bloom filter hashes in batches of 4

    const auto size = filter_.size();

    uint64_t hs[4] __attribute__ ((aligned (16)));
    uint64_t indices[4] __attribute__ ((aligned (16)));

    __m256i offsets, inds, mult;

    for (; j + 4 <= len; j += 4) {
        // copy hashes
        memcpy(hs, hashes + j, sizeof(uint64_t) * 4);

        // compute offsets
        // offset = ((hash % size) >> SHIFT) << SHIFT;
        offsets = _mm256_setr_epi64x(hs[0] % size,
                                     hs[1] % size,
                                     hs[2] % size,
                                     hs[3] % size);
        offsets = _mm256_srli_epi64(offsets, SHIFT);
        offsets = _mm256_slli_epi64(offsets, SHIFT);

        for (uint32_t k = 0; k < num_hash_functions_; ++k) {
            // load hashes to ymm register
            inds = _mm256_load_si256((__m256i*)hs);

            // compute the kth hashes
            // hash = (h1 + k * h2) & BLOCK_MASK
            mult = _mm256_setr_epi32(1, k, 1, k, 1, k, 1, k);
            inds = _mm256_mullo_epi32(inds, mult);
            inds = _mm256_hadd_epi32(inds, zeros);
            inds = _mm256_and_si256(inds, blockmask);
            inds = _mm256_permutevar8x32_epi32(inds, permute);

            // add offsets to hashes
            // id += offset
            inds = _mm256_add_epi64(offsets, inds);

            // store indices
            _mm256_store_si256((__m256i*)indices, inds);

            // add to filter
            for (size_t i = 0; i < 4; ++i) {
                filter_[indices[i]] = true;
            }
        }
    }

#endif

    for (; j < len; ++j) {
        insert(hashes[j]);
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

    uint64_t indices[4] __attribute__ ((aligned (16)));
    uint64_t hs[4] __attribute__ ((aligned (16)));

    bool found[4];

    __m256i offsets, inds, mult;

    size_t found_count = 0;

    for (; i + 4 <= num_elements; i += 4) {
        // reset found array
        found[0] = true;
        found[1] = true;
        found[2] = true;
        found[3] = true;

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

        for (uint32_t k = 0; k < num_hash_functions_; ++k) {
            // load hashes to ymm register
            inds = _mm256_load_si256((__m256i*)hs);

            // compute the kth hashes
            // hash = (h1 + k * h2) & BLOCK_MASK
            mult = _mm256_setr_epi32(1, k, 1, k, 1, k, 1, k);
            inds = _mm256_mullo_epi32(inds, mult);
            inds = _mm256_hadd_epi32(inds, zeros);
            inds = _mm256_and_si256(inds, blockmask);
            inds = _mm256_permutevar8x32_epi32(inds, permute);

            // add offsets to hashes
            // id += offset
            inds = _mm256_add_epi64(offsets, inds);

            // store indices
            _mm256_store_si256((__m256i*)indices, inds);

            // check indices
            found_count = 0;
            for (size_t j = 0; j < 4; ++j) {
                if (found[j]) {
                    found_count += (found[j] &= filter_[indices[j]]);
                }
            }

            if (!found_count)
                break;
        }

        for (size_t j = 0; j < 4; ++j) {
            if (found[j])
                presence[hash_index[i + j].second] = true;
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
