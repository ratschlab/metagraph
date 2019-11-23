#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include "common/hash/rolling_hasher.hpp"
#include "common/hash/bloom_filter.hpp"


// Bloom filter for approximate membership queries on k-mers
template <class KmerHasher = RollingHash<>>
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

    typedef std::function<void(const std::string&)> CallString;
    void add_sequences(const std::function<void(const CallString&)> &generate_sequences);

    bool find(const char *begin,
              const char *end,
              double discovery_fraction = 1.0) const;

    bool find(const std::string &sequence,
              double discovery_fraction = 1.0) const {
        return find(sequence.c_str(), sequence.c_str() + sequence.length(),
                    discovery_fraction);
    }

    // Checks for k-mer presence in the Bloom filter
    sdsl::bit_vector check_kmer_presence(const char *begin, const char *end) const;

    sdsl::bit_vector check_kmer_presence(const std::string &sequence) const {
        return check_kmer_presence(sequence.data(), sequence.data() + sequence.size());
    }

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }
    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); }

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    void print_stats() const;

  private:
    void call_kmers(const char *begin, const char *end,
                    const std::function<void(size_t /* position */,
                                             uint64_t /* hash */)> &callback,
                    const std::function<bool()> &terminate
                        = []() { return false; }) const;

    void call_kmer_presence(const char *begin,
                            const char *end,
                            const std::function<void(bool)> &callback,
                            const std::function<bool()> &terminate
                                = []() { return false; }) const;

    BloomFilter filter_;
    bool canonical_mode_;
    size_t k_;
    const KmerHasherType hasher_;
};


#endif // __KMER_BLOOM_FILTER_HPP__
