#include "bloom_filter.hpp"

#include <vector>

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
