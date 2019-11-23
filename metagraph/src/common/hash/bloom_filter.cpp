#include "bloom_filter.hpp"

#include <vector>

#include "utils/serialization.hpp"

constexpr size_t BLOCK_MASK = 0b111111111;
constexpr uint32_t SHIFT = 9;

BloomFilter::BloomFilter(size_t filter_size, size_t num_hash_functions)
      : filter_(((filter_size + BLOCK_MASK) >> SHIFT) << SHIFT),
        num_hash_functions_(num_hash_functions) {
    assert(filter_.size() >= filter_size);
    assert(!filter_.size() || filter_.size() > BLOCK_MASK);
    assert(!(filter_.size() & BLOCK_MASK));
}

BloomFilter::BloomFilter(size_t filter_size,
                         size_t expected_num_elements,
                         size_t max_num_hash_functions)
      : BloomFilter(filter_size,
                    std::min(max_num_hash_functions,
                             optim_h(filter_size, expected_num_elements))) {}

void BloomFilter::insert(uint64_t hash) {
    const auto size = filter_.size();
    if (!size)
        return;

    auto offset = ((hash % size) >> SHIFT) << SHIFT;
    uint32_t base = hash, jump = hash >> 32;

    /*
     * Kirsch, A., & Mitzenmacher, M. (2006, September).
     * Less hashing, same performance: building a better bloom filter.
     * In European Symposium on Algorithms (pp. 456-467). Springer, Berlin, Heidelberg.
     */
    for (size_t i = 0; i < num_hash_functions_; ++i) {
        filter_[offset + ((base + i * jump) & BLOCK_MASK)] = true;
    }

    assert(check(hash));
}

void BloomFilter::batch_insert(const uint64_t hashes[], size_t len) {
    const auto size = filter_.size();
    if (!size)
        return;

    std::vector<std::pair<uint64_t, uint64_t>> sorted_hashes(len);
    std::transform(
        hashes,
        hashes + len,
        sorted_hashes.begin(),
        [&](auto hash) {
            return std::make_pair(((hash % size) >> SHIFT) << SHIFT, hash);
        }
    );
    std::sort(sorted_hashes.begin(), sorted_hashes.end());

    for (auto [offset, hash] : sorted_hashes) {
        uint32_t base = hash, jump = hash >> 32;

        for (size_t i = 0; i < num_hash_functions_; ++i) {
            filter_[offset + ((base + i * jump) & BLOCK_MASK)] = true;
        }

        assert(check(hash));
    }
}

bool BloomFilter::check(uint64_t hash) const {
    const auto size = filter_.size();
    if (!size)
        return true;

    bool found = true;

    auto offset = ((hash % size) >> SHIFT) << SHIFT;
    uint32_t base = hash, jump = hash >> 32;

    for (size_t i = 0; i < num_hash_functions_; ++i) {
        found &= filter_[offset + ((base + i * jump) & BLOCK_MASK)];
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
