#include "bloom_filter.hpp"

#include <vector>
#include <immintrin.h>

#include "utils/serialization.hpp"

constexpr size_t BLOCK_MASK = 0b111111111;
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
    uint32_t h1 = hash;
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

void BloomFilter::batch_insert(const uint64_t hashes[], size_t len) {
    const auto size = filter_.size();

    std::vector<std::pair<size_t, uint64_t>> sorted_hashes(len);
    std::transform(
        hashes,
        hashes + len,
        sorted_hashes.begin(),
        [&](auto hash) {
            return std::make_pair(((hash % size) >> SHIFT) << SHIFT, hash);
        }
    );
    std::sort(sorted_hashes.begin(), sorted_hashes.end());

    // use the 64-bit hash to select a 512-bit block
    for (const auto &[offset, hash] : sorted_hashes) {
        // split 64-bit hash into two 32-bit hashes
        uint32_t h1 = hash;
        uint32_t h2 = hash >> 32;

        for (uint32_t i = 0; i < num_hash_functions_; ++i) {
            filter_[offset + ((h1 + i * h2) & BLOCK_MASK)] = true;
        }

        assert(check(hash));
    }
}

bool BloomFilter::check(uint64_t hash) const {
    const auto size = filter_.size();

    // use the 64-bit hash to select a 512-bit block
    size_t offset = ((hash % size) >> SHIFT) << SHIFT;

    // split 64-bit hash into two 32-bit hashes
    uint32_t h1 = hash;
    uint32_t h2 = hash >> 32;

    bool found = true;
    for (uint32_t i = 0; i < num_hash_functions_; ++i) {
        // check bits within block
        found &= filter_[offset + ((h1 + i * h2) & BLOCK_MASK)];
    }

    return found;
}

sdsl::bit_vector BloomFilter::batch_check(const uint64_t hashes[], size_t len) const {
    sdsl::bit_vector presence(len, false);
    size_t i = 0;

    const auto size = filter_.size();
    uint64_t offsets[4];
    uint32_t h[8];

    for (; i + 4 <= len; i += 4) {
        // offset = ((hash % size) >> SHIFT) << SHIFT;
        offsets[0] = hashes[i] % size;
        offsets[1] = hashes[i + 1] % size;
        offsets[2] = hashes[i + 2] % size;
        offsets[3] = hashes[i + 3] % size;

        // __m256i ymm = _mm256_castpd_si256(_mm256_loadu_pd((const double*)offsets));
        // ymm = _mm256_srli_epi64(ymm, SHIFT);
        // ymm = _mm256_slli_epi64(ymm, SHIFT);

        // // dump hashes
        // _mm256_storeu_ps((float*)h, _mm256_castsi256_ps(ymm));

        offsets[0] = (offsets[0] >> SHIFT) << SHIFT;
        offsets[1] = (offsets[1] >> SHIFT) << SHIFT;
        offsets[2] = (offsets[2] >> SHIFT) << SHIFT;
        offsets[3] = (offsets[3] >> SHIFT) << SHIFT;
        memcpy(h, (const uint32_t*)offsets, sizeof(uint32_t) * 8);

        for (size_t j = 0; j < 4; ++j) {
            bool found = true;
            for (size_t k = 0; found && k < num_hash_functions_; ++k) {
                found &= filter_[offsets[j] + ((h[j * 2] + k * h[j * 2 + 1]) & BLOCK_MASK)];
            }

            if (found)
                presence[i + j] = true;
        }
    }

    for (; i < len; ++i) {
        presence[i] = check(hashes[i]);
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
