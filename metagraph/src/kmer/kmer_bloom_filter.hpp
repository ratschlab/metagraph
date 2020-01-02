#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include "common/hash/rolling_hasher.hpp"
#include "common/hash/bloom_filter.hpp"


// Bloom filter for approximate membership queries on k-mers
template <class KmerHasher = RollingMultiHash<2>>
class KmerBloomFilter {
  public:
    typedef KmerHasher KmerHasherType;

    template <typename... Args>
    KmerBloomFilter(size_t k, bool canonical_mode, Args&&... args)
          : filter_(std::forward<Args>(args)...),
            canonical_mode_(canonical_mode),
            k_(k),
            hasher_(k_) {}

    // Add the k-mers of the sequence to the Bloom filter
    void add_sequence(std::string_view sequence);

    // Checks for k-mer presence in the Bloom filter
    sdsl::bit_vector check_kmer_presence(std::string_view sequence) const;

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }
    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); }

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    KmerHasherType get_hasher() const { return hasher_; }

  private:
    BloomFilter filter_;
    bool canonical_mode_;
    size_t k_;
    const KmerHasherType hasher_;
};


/**
 * Construct a callback, where the i-th call returns `true` if the i-th k-mer
 * in the sequence is not rejected by the Bloom filter and `false` otherwise.
 * If `bloom_filter` is NULL, assume that none of the k-mers are missing.
 */
std::function<bool()> get_missing_kmer_skipper(const KmerBloomFilter<> *bloom_filter,
                                               std::string_view sequence);

#endif // __KMER_BLOOM_FILTER_HPP__
