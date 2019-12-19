#include "bloom_filter.hpp"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "utils/serialization.hpp"

constexpr uint32_t BLOCK_MASK = 0b111111111;
constexpr uint32_t SHIFT = 9;

BloomFilter::BloomFilter(size_t filter_size, uint32_t num_hash_functions)
      : filter_(filter_size ? (((filter_size + BLOCK_MASK) >> SHIFT) << SHIFT)
                            : BLOCK_MASK + 1),
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

static_assert(sizeof(long long unsigned int) == sizeof(uint64_t));

// fast map of h uniformly to the region [0, size)
// (h * size) >> 64
inline uint64_t restrict_to(long long unsigned int h, size_t size) {
#ifdef __BMI2__
    _mulx_u64(h, size, &h);

    assert(h < size);
    return h;
#else
    assert(((static_cast<__uint128_t>(h) * size) >> 64) < size);

    return (static_cast<__uint128_t>(h) * size) >> 64;
#endif
}

void BloomFilter::insert(uint64_t hash) {
    // use the 64-bit hash to select a 512-bit block
    const size_t offset = (restrict_to(hash, filter_.size()) >> SHIFT) << SHIFT;

    // split 64-bit hash into two 32-bit hashes
    const uint32_t h1 = hash & 0xFFFFFFFF;
    const uint32_t h2 = hash >> 32;

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
    // use the 64-bit hash to select a 512-bit block
    const size_t offset = (restrict_to(hash, filter_.size()) >> SHIFT) << SHIFT;

    // split 64-bit hash into two 32-bit hashes
    const uint32_t h1 = hash & 0xFFFFFFFF;
    const uint32_t h2 = hash >> 32;

    for (uint32_t i = 0; i < num_hash_functions_; ++i) {
        // check bits within block
        if (!filter_[offset + ((h1 + i * h2) & BLOCK_MASK)])
            return false;
    }

    return true;
}


#ifdef __AVX2__

// compute Bloom filter hashes in batches of 4
const uint64_t* batch_insert_avx2(const BloomFilter &bloom,
                                  sdsl::bit_vector &filter_,
                                  const uint32_t num_hash_functions_,
                                  const uint64_t *hs,
                                  const uint64_t *end) {
    // only used for assert
    std::ignore = bloom;

    // used to restrict indices to the size of a block
    const __m128i blockmask = _mm_set1_epi32(BLOCK_MASK);

    // used to set offset from block start
    const __m128i shift = _mm_setr_epi32(6, 6, 0, 0);

    // used to select bit
    const __m128i andmask = _mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0x3F, 0x3F);

    const __m128i ones = _mm_set1_epi64x(1);

    const size_t size = filter_.size();

    uint64_t indices[4] __attribute__ ((aligned (32)));
    uint32_t ht[2] __attribute__ ((aligned (8)));
    uint64_t updated_block[2] __attribute__ ((aligned (16)));
    uint32_t *hh;
    uint64_t *block;
    uint32_t offset;

    __m256i block_indices;
    __m128i hashes_init, hashes, mult;
    __m64 *harray = (__m64*)&hashes;

    for (; hs + 4 <= end; hs += 4) {
        // check next batch to see if all hashes are absent
        // compute offsets
        block_indices = _mm256_setr_epi64x(
            restrict_to(hs[0], size),
            restrict_to(hs[1], size),
            restrict_to(hs[2], size),
            restrict_to(hs[3], size)
        );
        block_indices = _mm256_srli_epi64(block_indices, SHIFT);
        block_indices = _mm256_slli_epi64(block_indices, SHIFT - 6);

        _mm256_store_si256((__m256i*)indices, block_indices);

        // clean up after AVX2 instructions
        _mm256_zeroupper();

        for (size_t j = 0; j < 4; ++j) {
            uint32_t k = 0;

            block = filter_.data() + indices[j];

            // ensure the start of the 512-bit block is prefetched, prepare for write
            __builtin_prefetch(block, 1);

            hashes_init = _mm_set1_epi64x(hs[j]);

            // hash functions in pairs
            for (; k + 2 <= num_hash_functions_; k += 2) {
                // load hash
                hashes = hashes_init;

                // compute hashes
                // hash = (h1 + k * h2) & BLOCK_MASK
                mult = _mm_setr_epi32(1, k, 1, k + 1);
                hashes = _mm_mullo_epi32(hashes, mult);
                hashes = _mm_hadd_epi32(hashes, hashes);
                hashes = _mm_and_si128(hashes, blockmask);

                // add to filter
                // block[hash / 64] |= 1llu << (hash % 64)
                hashes = _mm_srlv_epi32(hashes, shift);
                hashes = _mm_and_si128(hashes, andmask);

                _mm_store_si128(
                    (__m128i*)updated_block,
                    _mm_sllv_epi64(ones, _mm_cvtepi32_epi64(_mm_set1_epi64(harray[1])))
                );
                _mm_store_si128((__m128i*)ht, hashes);

                assert(ht[0] < 8);
                assert(ht[1] < 8);

                block[ht[0]] |= updated_block[0];
                block[ht[1]] |= updated_block[1];
            }

            // if num_hash_functions is odd, add the last one
            if (num_hash_functions_ & 1) {
                hh = (uint32_t*)&hs[j];
                offset = (hh[0] + (num_hash_functions_ - 1) * hh[1]) & BLOCK_MASK;
                block[offset >> 6] |= 1llu << (offset & 0x3F);
            }
        }

        assert(std::all_of(hs, hs + 4, [&](uint64_t hash) { return bloom.check(hash); }));
    }

    return hs;
}

#endif

void BloomFilter::insert(const uint64_t *hashes_begin, const uint64_t *hashes_end) {
#ifdef __AVX2__
    hashes_begin = batch_insert_avx2(*this,
                                     filter_,
                                     num_hash_functions_,
                                     hashes_begin,
                                     hashes_end);
#endif

    // insert residual
    for (; hashes_begin < hashes_end; ++hashes_begin) {
        insert(*hashes_begin);
    }
}

#ifdef __AVX2__

// compute Bloom filter hashes in batches of 4
uint64_t batch_check_avx2(const BloomFilter &bloom,
                          const uint64_t *hashes_begin,
                          const uint64_t *hashes_end,
                          sdsl::bit_vector &presence,
                          const sdsl::bit_vector &filter_,
                          const uint32_t num_hash_functions_) {
    // only used for assert
    std::ignore = bloom;

    // used to restrict indices to the size of a block
    const __m128i blockmask = _mm_set1_epi32(BLOCK_MASK);

    // used to set offset from block start
    const __m128i shift = _mm_setr_epi32(6, 6, 0, 0);

    // used to select bit
    const __m128i andmask = _mm_setr_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0x3F, 0x3F);

    const __m128i ones = _mm_set1_epi64x(1);

    const size_t size = filter_.size();

    uint64_t indices[4] __attribute__ ((aligned (32)));
    uint32_t *hh;
    const uint64_t *block;
    uint32_t offset;

    __m256i block_indices;
    __m128i hashes_init, hashes, mult;
    __m64 *harray = (__m64*)&hashes;

    size_t i = 0;
    const size_t num_elements = hashes_end - hashes_begin;
    for (; i + 4 <= num_elements; i += 4) {
        // copy hashes
        block_indices = _mm256_setr_epi64x(
            restrict_to(hashes_begin[0], size),
            restrict_to(hashes_begin[1], size),
            restrict_to(hashes_begin[2], size),
            restrict_to(hashes_begin[3], size)
        );
        block_indices = _mm256_srli_epi64(block_indices, SHIFT);
        block_indices = _mm256_slli_epi64(block_indices, SHIFT - 6);

        _mm256_store_si256((__m256i*)indices, block_indices);

        // clean up after AVX2 instructions
        _mm256_zeroupper();

        for (size_t j = 0; j < 4; ++j) {
            bool found = true;
            uint32_t k = 0;

            block = filter_.data() + indices[j];

            // ensure the start of the 512-bit block is prefetched, prepare for read
            __builtin_prefetch(block, 0);

            hashes_init = _mm_set1_epi64x(hashes_begin[j]);

            // hash functions in pairs
            for (; found && k + 2 <= num_hash_functions_; k += 2) {
                // load hash
                hashes = hashes_init;

                // compute hashes
                // hash = (h1 + k * h2) & BLOCK_MASK
                mult = _mm_setr_epi32(1, k, 1, k + 1);
                hashes = _mm_mullo_epi32(hashes, mult);
                hashes = _mm_hadd_epi32(hashes, hashes);
                hashes = _mm_and_si128(hashes, blockmask);

                // check hashes
                // found &= bool(block[hash / 64] & (1llu << (hash % 64)))
                hashes = _mm_srlv_epi32(hashes, shift);
                hashes = _mm_and_si128(hashes, andmask);
                found &= _mm_testc_si128(
                    _mm_srlv_epi64(
                        _mm_i32gather_epi64((long long int*)block, hashes, 8),
                        _mm_cvtepi32_epi64(_mm_set1_epi64(harray[1]))
                    ),
                    ones
                );
            }

            // if num_hash_functions is odd, check the last one
            if (found && (num_hash_functions_ & 1)) {
                hh = (uint32_t*)&hashes_begin[j];
                offset = (hh[0] + (num_hash_functions_ - 1) * hh[1]) & BLOCK_MASK;
                found = block[offset >> 6] & (1llu << (offset & 0x3F));
            }

            if (found) {
                presence[i + j] = true;
                assert(bloom.check(hashes_begin[j]));
            }
        }

        hashes_begin += 4;
    }

    return i;
}

#endif

sdsl::bit_vector BloomFilter
::check(const uint64_t *hashes_begin, const uint64_t *hashes_end) const {
    assert(hashes_end >= hashes_begin);

    sdsl::bit_vector presence(hashes_end - hashes_begin, false);
    size_t i = 0;
    const size_t num_elements = hashes_end - hashes_begin;

#ifdef __AVX2__
    i = batch_check_avx2(*this,
                         hashes_begin,
                         hashes_end,
                         presence,
                         filter_,
                         num_hash_functions_);
#endif

    // check residual
    for (; i < num_elements; ++i) {
        if (check(hashes_begin[i]))
            presence[i] = true;
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
