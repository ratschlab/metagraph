#include "bloom_filter.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/serialization.hpp"

// used to implement size % 512
constexpr uint64_t BLOCK_MASK = 0x1FF;
constexpr uint64_t BLOCK_MASK_OUT = ~BLOCK_MASK;


BloomFilter::BloomFilter(size_t filter_size, uint32_t num_hash_functions)
      : filter_(filter_size ? (filter_size + BLOCK_MASK) & BLOCK_MASK_OUT
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

void BloomFilter::insert(uint64_t hash) {
    // use the 64-bit hash to select a 512-bit block
    const size_t offset = restrict_to(hash, filter_.size()) & BLOCK_MASK_OUT;

    // split 64-bit hash into two 32-bit hashes
    const uint64_t h_hi = hash >> 32;
    const uint64_t h_lo = hash & 0xFFFFFFFF;

    uint64_t id;

    __builtin_prefetch(filter_.data() + (offset >> 6));

    /*
     * Use two 32-bit hashes to generate num_hash_functions hashes
     * Kirsch, A., & Mitzenmacher, M. (2006, September).
     * Less hashing, same performance: building a better bloom filter.
     * In European Symposium on Algorithms (pp. 456-467). Springer, Berlin, Heidelberg.
     */
    for (uint64_t i = 0; i < num_hash_functions_; ++i) {
        // set bits within block
        id = offset + ((h_lo * i + h_hi) & BLOCK_MASK);
        filter_.data()[id >> 6] |= 1llu << (id & 0x3F);
    }

    assert(check(hash));
}

bool BloomFilter::check(uint64_t hash) const {
    // use the 64-bit hash to select a 512-bit block
    const size_t offset = restrict_to(hash, filter_.size()) & BLOCK_MASK_OUT;

    // split 64-bit hash into two 32-bit hashes
    const uint64_t h_hi = hash >> 32;
    const uint64_t h_lo = hash & 0xFFFFFFFF;

    uint64_t id;

    __builtin_prefetch(filter_.data() + (offset >> 6));

    for (uint64_t i = 0; i < num_hash_functions_; ++i) {
        // check bits within block
        id = offset + ((h_lo * i + h_hi) & BLOCK_MASK);
        if (!((filter_.data()[id >> 6] >> (id & 0x3F)) & 1))
            return false;
    }

    return true;
}


#ifdef __AVX2__

/**
 * If num_hash_functions_ == 1, then four elements can he inserted/checked at
 * once. If num_hash_functions_ > 1, then each element is inserted/checked
 * individually, but four hash functions at a time.
 */

// compute Bloom filter hashes in batches of 4
inline const uint64_t*
batch_insert_avx2(BloomFilter &bloom,
                  const uint32_t num_hash_functions_,
                  const uint64_t *hashes_begin,
                  const uint64_t *hashes_end) {
    assert(num_hash_functions_);

    auto &filter_ = bloom.data();

    static_assert(sizeof(long long int) == sizeof(uint64_t));
    long long int *filter_cast = reinterpret_cast<long long int*>(filter_.data());
    long long int *offset;

    const size_t size = filter_.size();

    const __m256i block_mask = _mm256_set1_epi64x(BLOCK_MASK);
    const __m256i block_mask_out = _mm256_set1_epi64x(BLOCK_MASK_OUT);
    const __m256i mod_mask = _mm256_set1_epi64x(0x3F);
    const __m256i ones = _mm256_set1_epi64x(1);
    const __m256i add = _mm256_setr_epi64x(0, 1, 2, 3);
    const __m256i numhash = _mm256_set1_epi64x(num_hash_functions_);

    uint64_t ids_store[4] __attribute__ ((aligned (32)));
    uint64_t updates_store[4] __attribute__ ((aligned (32)));
    __m256i *ids_cast = reinterpret_cast<__m256i*>(&ids_store);
    __m256i *updates_cast = reinterpret_cast<__m256i*>(&updates_store);

    // check four input elements (represented by hashes) at a time
    for (; hashes_begin + 4 <= hashes_end; hashes_begin += 4) {
        __m256i block_indices = restrict_to_mask_epi64(hashes_begin, size, block_mask_out);
        const uint64_t *block_index_array = reinterpret_cast<const uint64_t*>(&block_indices);

        if (num_hash_functions_ > 1) {
            for (size_t j = 0; j < 4; ++j) {
                __m256i h = _mm256_set1_epi64x(hashes_begin[j]);
                __m256i h_hi = _mm256_srli_epi64(h, 32);

                offset = filter_cast + (block_index_array[j] >> 6);

                // ensure that the entire 512-bit block is cached
                __builtin_prefetch(offset, 0);

                // compute hash functions four at a time
                for (uint32_t k = 0; k < num_hash_functions_; k += 4) {
                    // hash = (h_hi + k * h_lo) & BLOCK_MASK
                    __m256i mult = _mm256_add_epi64(_mm256_set1_epi64x(k), add);
                    __m256i hash = _mm256_and_si256(
                        _mm256_add_epi64(_mm256_mul_epu32(h, mult), h_hi),
                        block_mask
                    );

                    _mm256_store_si256(ids_cast, _mm256_srli_epi64(hash, 6));
                    _mm256_store_si256(
                        updates_cast,
                        _mm256_and_si256(
                            _mm256_sllv_epi64(ones, _mm256_and_si256(hash, mod_mask)),
                            _mm256_cmpgt_epi64(numhash, mult)
                        )
                    );

                    // don't do these with vector operations since values in
                    // ids_store may match
                    offset[ids_store[0]] |= updates_store[0];
                    offset[ids_store[1]] |= updates_store[1];
                    offset[ids_store[2]] |= updates_store[2];
                    offset[ids_store[3]] |= updates_store[3];
                }
            }

        } else {
            // insert all four elements in parallel
            // hash = h_hi & BLOCK_MASK
            __m256i hash = _mm256_add_epi64(
                _mm256_and_si256(
                    _mm256_srli_epi64(_mm256_loadu_si256((__m256i*)hashes_begin), 32),
                    block_mask
                ),
                block_indices
            );

            _mm256_store_si256(ids_cast, _mm256_srli_epi64(hash, 6));
            _mm256_store_si256(updates_cast,
                               _mm256_sllv_epi64(ones, _mm256_and_si256(hash, mod_mask)));

            // don't do these with vector operations since values in
            // ids_store may match
            filter_cast[ids_store[0]] |= updates_store[0];
            filter_cast[ids_store[1]] |= updates_store[1];
            filter_cast[ids_store[2]] |= updates_store[2];
            filter_cast[ids_store[3]] |= updates_store[3];
        }

        assert(std::all_of(hashes_begin, hashes_begin + 4,
                           [&](uint64_t h) { return bloom.check(h); }));
    }

    return hashes_begin;
}


// compute Bloom filter hashes in batches of 4
inline uint64_t
batch_check_avx2(const BloomFilter &bloom,
                 const uint64_t *hashes_begin,
                 const uint64_t *hashes_end,
                 const uint64_t num_hash_functions_,
                 const std::function<void(size_t)> &present_index_callback) {
    assert(num_hash_functions_);

    static_assert(sizeof(long long int) == sizeof(uint64_t));
    const auto *filter_cast = reinterpret_cast<const long long int*>(bloom.data().data());
    const long long int *offset;

    const size_t size = bloom.size();

    const __m256i block_mask = _mm256_set1_epi64x(BLOCK_MASK);
    const __m256i block_mask_out = _mm256_set1_epi64x(BLOCK_MASK_OUT);
    const __m256i mod_mask = _mm256_set1_epi64x(0x3F);
    const __m256i ones = _mm256_set1_epi64x(1);
    const __m256i all = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
    const __m256i add = _mm256_setr_epi64x(0, 1, 2, 3);
    const __m256i numhash = _mm256_set1_epi64x(num_hash_functions_);

    __m256i block_indices;
    const uint64_t *block_index_array = reinterpret_cast<const uint64_t*>(&block_indices);

    // check four input elements (represented by hashes) at a time
    size_t i = 0;
    for (; hashes_begin + 4 <= hashes_end; hashes_begin += 4) {
        block_indices = restrict_to_mask_epi64(hashes_begin, size, block_mask_out);

        if (num_hash_functions_ > 1) {
            for (size_t j = 0; j < 4; ++j) {
                __m256i h = _mm256_set1_epi64x(hashes_begin[j]);
                __m256i h_hi = _mm256_srli_epi64(h, 32);

                offset = filter_cast + (block_index_array[j] >> 6);

                // ensure that the entire 512-bit block is cached
                __builtin_prefetch(offset, 0);

                // compute hash functions four at a time
                bool found = true;
                uint32_t k = 0;
                for (; found && k < num_hash_functions_; k += 4) {
                    // hash = (h_hi + k * h_lo) & BLOCK_MASK
                    __m256i mult = _mm256_add_epi64(_mm256_set1_epi64x(k), add);
                    __m256i hash = _mm256_and_si256(
                        _mm256_add_epi64(_mm256_mul_epu32(h, mult), h_hi),
                        block_mask
                    );

                    // test all four hashes
                    // word = offset[uint32_t(hash / 64)] >> (hash % 64)
                    // For hash indices greater than num_hash_functions_ - 1,
                    // fill the word with -1 so that they test true
                    found &= _mm256_testc_si256(
                        _mm256_srlv_epi64(
                            _mm256_mask_i32gather_epi64(
                                all,
                                offset,
                                cvtepi64_epi32(_mm256_srli_epi64(hash, 6)),
                                _mm256_cmpgt_epi64(numhash, mult),
                                8
                            ),
                            _mm256_and_si256(hash, mod_mask)
                        ),
                        ones
                    );
                }

                assert(found == bloom.check(hashes_begin[j]));

                if (found)
                    present_index_callback(i + j);
            }

        } else {
            // check all four elements in parallel
            __m256i hash = _mm256_add_epi64(
                _mm256_and_si256(
                    _mm256_srli_epi64(_mm256_loadu_si256((__m256i*)hashes_begin), 32),
                    block_mask
                ),
                block_indices
            );

            // the ith byte is set to FF iff the ith input hash is present in the filter
            int32_t test = _mm256_movemask_epi8(_mm256_cmpeq_epi64(
                _mm256_and_si256(_mm256_srlv_epi64(
                    _mm256_i64gather_epi64(filter_cast, _mm256_srli_epi64(hash, 6), 8),
                    _mm256_and_si256(hash, mod_mask)
                ), ones),
                ones
            ));

            assert(bool(test & 0x1) == bloom.check(hashes_begin[0]));
            assert(bool(test & 0x100) == bloom.check(hashes_begin[1]));
            assert(bool(test & 0x10000) == bloom.check(hashes_begin[2]));
            assert(bool(test & 0x1000000) == bloom.check(hashes_begin[3]));

            if (test & 0x1)
                present_index_callback(i);

            if (test & 0x100)
                present_index_callback(i + 1);

            if (test & 0x10000)
                present_index_callback(i + 2);

            if (test & 0x1000000)
                present_index_callback(i + 3);

        }

        i += 4;
    }

    return i;
}

#endif


void BloomFilter::insert(const uint64_t *hashes_begin, const uint64_t *hashes_end) {
    if (!num_hash_functions_) {
        assert(std::all_of(hashes_begin, hashes_end,
                           [&](uint64_t hash) { return check(hash); }));

        return;
    }

    const uint64_t *it = hashes_begin;

#ifdef __AVX2__
    it = batch_insert_avx2(*this, num_hash_functions_, it, hashes_end);
#endif

    // insert residual
    for (; it < hashes_end; ++it) {
        insert(*it);
        assert(check(*it));
    }

    assert(std::all_of(hashes_begin, hashes_end,
                       [&](uint64_t hash) { return check(hash); }));
}

void BloomFilter::check(const uint64_t *hashes_begin,
                        const uint64_t *hashes_end,
                        const std::function<void(size_t)> &present_index_callback) const {
    assert(hashes_end >= hashes_begin);

    const size_t num_elements = hashes_end - hashes_begin;

    if (!num_hash_functions_) {
        for (size_t i = 0; i < num_elements; ++i) {
            present_index_callback(i);
        }

        return;
    }

    size_t i = 0;

#ifdef __AVX2__
    i = batch_check_avx2(*this,
                         hashes_begin,
                         hashes_end,
                         num_hash_functions_,
                         present_index_callback);
#endif

    // check residual
    for (; i < num_elements; ++i) {
        if (check(hashes_begin[i]))
            present_index_callback(i);
    }
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
