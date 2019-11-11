#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include <iostream>

#include "rolling_hasher.hpp"
#include "bloom_filter.hpp"


// Bloom filter for approximate membership queries on k-mers
template <class KmerHasher = RollingMultiHasher<2>>
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
    void add_sequence(const char *begin, const char *end);
    void add_sequence(const std::string &sequence) {
        add_sequence(sequence.data(), sequence.data() + sequence.size());
    }

    // Checks for k-mer presence in the Bloom filter
    sdsl::bit_vector check_kmer_presence(const char *begin,
                                         const char *end) const;

    sdsl::bit_vector check_kmer_presence(const std::string &sequence) const {
        return check_kmer_presence(sequence.data(), sequence.data() + sequence.size());
    }

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }
    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); }

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    void print_stats() const {
        std::cout << "Bloom filter parameters" << std::endl
                  << "Size:\t\t\t" << size() << " bits" << std::endl
                  << "Num hash functions:\t" << num_hash_functions() << std::endl;
    }

  private:
    void call_kmers(const char *begin, const char *end,
                    const std::function<void(size_t /* position */,
                                             uint64_t /* hash1 */,
                                             uint64_t /* hash2 */)> &callback) const;

    BloomFilter filter_;
    bool canonical_mode_;
    size_t k_;
    KmerHasherType hasher_;
};


#endif // __KMER_BLOOM_FILTER_HPP__
