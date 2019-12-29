#include "bloom_filter.hpp"

#include "common/serialization.hpp"


BloomFilter::BloomFilter(size_t filter_size, size_t num_hash_functions)
      : filter_(filter_size),
        num_hash_functions_(num_hash_functions) {}

BloomFilter::BloomFilter(size_t filter_size,
                         size_t expected_num_elements,
                         size_t max_num_hash_functions)
      : BloomFilter(filter_size,
                    std::min(max_num_hash_functions,
                             optim_h(filter_size, expected_num_elements))) {}

void BloomFilter::insert(uint64_t hash1, uint64_t hash2) {
    const auto size = filter_.size();
    if (!size)
        return;

    // Kirsch, A., & Mitzenmacher, M. (2006, September).
    // Less hashing, same performance: building a better bloom filter.
    // In European Symposium on Algorithms (pp. 456-467). Springer, Berlin, Heidelberg.
    for (size_t i = 0; i < num_hash_functions_; ++i) {
        // TODO: do some locking here to make this multithreaded
        const auto hash = hash1 + i * hash2;
        filter_[hash - hash / size * size] = true;
    }

    assert(check(hash1, hash2));
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
