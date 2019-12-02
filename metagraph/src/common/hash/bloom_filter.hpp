#ifndef __BLOOM_FILTER_HPP__
#define __BLOOM_FILTER_HPP__

#include <iostream>
#include <cmath>

#include <sdsl/int_vector.hpp>


class BloomFilter {
  public:
    BloomFilter(size_t filter_size = 0, size_t num_hash_functions = 0);

    BloomFilter(size_t filter_size,
                size_t expected_num_elements,
                size_t max_num_hash_functions);

    // TODO: pass all hashes `insert(uint64_t hashes[])`
    void insert(uint64_t hash1, uint64_t hash2);

    // TODO: pass all hashes `insert(uint64_t hashes[])`
    inline bool check(uint64_t hash1, uint64_t hash2) const {
        const auto size = filter_.size();
        if (!size)
            return true;

        for (size_t i = 0; i < num_hash_functions_; ++i) {
            const auto hash = hash1 + i * hash2;
            if (!filter_[hash - hash / size * size])
                return false;
        }

        return true;
    }

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return num_hash_functions_; }

    constexpr static size_t optim_size(double false_positive_prob,
                                       size_t expected_num_elements) {
        if (false_positive_prob <= 0.0)
            throw std::runtime_error("False positive probability must be > 0.0");

        return -std::log2(false_positive_prob) * expected_num_elements / M_LN2 / M_LN2;
    }

    constexpr static size_t optim_h(double false_positive_prob) {
        if (false_positive_prob <= 0.0)
            throw std::runtime_error("False positive probability must be > 0.0");

        return std::ceil(-std::log2(false_positive_prob));
    }

    constexpr static size_t optim_h(size_t filter_size, size_t expected_num_elements) {
        return expected_num_elements
            ? std::ceil(M_LN2 * filter_size / expected_num_elements)
            : static_cast<size_t>(-1);
    }

  private:
    sdsl::bit_vector filter_;
    size_t num_hash_functions_;
};


#endif // __BLOOM_FILTER_HPP__
