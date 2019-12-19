#ifndef __BLOOM_FILTER_HPP__
#define __BLOOM_FILTER_HPP__

#include <iostream>
#include <vector>
#include <cmath>

#include <sdsl/int_vector.hpp>


/**
 * BloomFilter is a probabilistic data structure for supporting approximate
 * membership queries.
 *
 * This structure is implemented as a more cache-efficient blocked Bloom filter,
 * in which all filter indices computed from an input single hash value are
 * in the same 512-bit block.
 */
class BloomFilter {
  public:
    /**
     * Constructs a Bloom filter from the given parameters.
     * @param filter_size set the filter size to be the next greatest multiple of
     * 512 and filter_size
     * @param num_hash_functions the number of hash functions to compute per element
     */
    BloomFilter(size_t filter_size = 0, uint32_t num_hash_functions = 0);

    /**
     * Constructs a Bloom filter with optimal parameters.
     * @param filter_size set the filter size to be the next greatest multiple of
     * 512 and filter_size
     * @param expected_num_elements the number of elements expected to be inserted
     * @param max_num_hash_functions upper bound on the number of hash functions
     * to use per element
     */
    BloomFilter(size_t filter_size,
                size_t expected_num_elements,
                uint32_t max_num_hash_functions);

    /**
     * Insert an element into the Bloom filter.
     * @param hash the hash of the inserted element
     */
    void insert(uint64_t hash);

    /**
     * Insert a batch of elements into the Bloom filter by hash value
     */
    void insert(const uint64_t *hashes_begin, const uint64_t *hashes_end);

    /**
     * Check an element for presence in the Bloom filter. Returns false if
     * the element is not in the set and true if the element may be in the set.
     * @param hash the hash of the inserted element
     */
    bool check(uint64_t hash) const;

    /**
     * Check a batch of elements in the Bloom filter and report their
     * presence/absence.
     * @param hash_index a vector of pairs defining the hash value and corresponding
     * index in the returned bit_vector for each element. Each index must be
     * less than length.
     * @param length the length of the returned bit_vector
     */
    sdsl::bit_vector check(const uint64_t *hashes_begin, const uint64_t *hashes_end) const;

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    size_t size() const { return filter_.size(); }
    uint32_t num_hash_functions() const { return num_hash_functions_; }

    const sdsl::bit_vector& data() const { return filter_; }

    /**
     * Compute the optimal size of a Bloom filter given the parameters.
     */
    constexpr static size_t optim_size(double false_positive_prob,
                                       size_t expected_num_elements) {
        if (false_positive_prob <= 0.0)
            throw std::runtime_error("False positive probability must be > 0.0");

        return -std::log2(false_positive_prob) * expected_num_elements / M_LN2;
    }

    /**
     * Compute the optimal number of hash functions given the input parameters
     */
    constexpr static uint32_t optim_h(double false_positive_prob) {
        if (false_positive_prob <= 0.0)
            throw std::runtime_error("False positive probability must be > 0.0");

        return std::ceil(-std::log2(false_positive_prob));
    }

    /**
     * Compute the optimal number of hash functions given the input parameters
     */
    constexpr static uint32_t optim_h(size_t filter_size, size_t expected_num_elements) {
        return expected_num_elements
            ? std::ceil(M_LN2 * filter_size / expected_num_elements)
            : static_cast<size_t>(-1);
    }

  private:
    sdsl::bit_vector filter_;
    uint32_t num_hash_functions_;
};


#endif // __BLOOM_FILTER_HPP__
