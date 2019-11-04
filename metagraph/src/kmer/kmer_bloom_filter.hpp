#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include <fstream>
#include <cmath>
#include <memory>

#include <sdsl/int_vector.hpp>
#include <cyclichash.h>


template <typename TAlphabet = uint8_t>
class RollingKmerHasher {
  public:
    static constexpr TAlphabet MAXVAL = std::numeric_limits<TAlphabet>::max();

    // Note: this constructor is expensive. Try to construct it once and
    // make copies of the object.
    explicit RollingKmerHasher(size_t k, uint32_t seed1 = 0, uint32_t seed2 = 1)
          : hash_(k, seed1, seed2, 64) {
        // initialize
        for (int i = 0; i < hash_.n; ++i) {
            hash_.eat(MAXVAL);
        }

        // store initialized value for faster reset
        reset_value_ = hash_.hashvalue;
    }

    // After calling reset, the first k values of prev (next) in shift left (right)
    // should be MAXVAL
    void reset() { hash_.hashvalue = reset_value_; }

    void shift_left(TAlphabet next, TAlphabet prev) { hash_.update(prev, next); }

    void shift_right(TAlphabet prev, TAlphabet next) { hash_.reverse_update(prev, next); }

    bool operator<(const RollingKmerHasher &other) const {
        return hash_.hashvalue < other.hash_.hashvalue;
    }

    bool operator>(const RollingKmerHasher &other) const {
        return hash_.hashvalue > other.hash_.hashvalue;
    }

    bool operator==(const RollingKmerHasher &other) const {
        return hash_.hashvalue == other.hash_.hashvalue;
    }

    size_t get_k() const { return hash_.n; }

    uint64_t get_hash() const { return hash_.hashvalue; }

  private:
    CyclicHash<uint64_t, TAlphabet> hash_;
    decltype(hash_.hashvalue) reset_value_;
};


template <int h,
          typename TAlphabet = uint8_t,
          class RollingKmerHasher = ::RollingKmerHasher<TAlphabet>>
class RollingKmerMultiHasher {
  public:
    static constexpr TAlphabet MAXVAL = RollingKmerHasher::MAXVAL;

    explicit RollingKmerMultiHasher(size_t k) {
        static_assert(h);

        hashers_.reserve(h);
        for (size_t i = 0; i < h; ++i) {
            hashers_.emplace_back(k, i * 2, i * 2 + 1);
        }

        assert(hashers_.size() == h);
    }

    // After calling reset, the first k values of prev (next) in shift left (right)
    // should be MAXVAL
    void reset() {
        assert(hashers_.size() == h);
        for (auto &hasher : hashers_) {
            hasher.reset();
        }
    }

    void shift_left(TAlphabet next, TAlphabet prev) {
        assert(hashers_.size() == h);
        for (auto &hasher : hashers_) {
            hasher.shift_left(next, prev);
        }
    }

    void shift_right(TAlphabet prev, TAlphabet next) {
        assert(hashers_.size() == h);
        for (auto &hasher : hashers_) {
            hasher.shift_right(prev, next);
        }
    }

    bool operator<(const RollingKmerMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ < other.hashers_;
    }

    bool operator>(const RollingKmerMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ > other.hashers_;
    }

    bool operator==(const RollingKmerMultiHasher &other) const {
        assert(hashers_.size() == h);
        return hashers_ == other.hashers_;
    }

    size_t get_k() const {
        assert(hashers_.size() == h);
        return hashers_.at(0).get_k();
    }

    template <int j>
    uint64_t get_hash() const {
        static_assert(j < h);
        assert(hashers_.size() == h);
        return hashers_.at(j).get_hash();
    }

  private:
    std::vector<RollingKmerHasher> hashers_;
};

class BloomFilter {
  public:
    BloomFilter(size_t filter_size = 0, size_t num_hash_functions = 0);

    BloomFilter(size_t filter_size,
                size_t expected_num_elements,
                size_t max_num_hash_functions);

    void insert(uint64_t hash1, uint64_t hash2);
    bool check(uint64_t hash1, uint64_t hash2) const;

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

// Bloom filter for approximate membership queries on k-mers
// Sizes are round up to the nearest power of 2 to allow for
// faster querying and insertion
template <class KmerHasher = RollingKmerMultiHasher<2>>
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
        add_sequence(&*sequence.begin(), &*sequence.end());
    }

    // Checks for k-mer presence in the Bloom filter
    sdsl::bit_vector check_kmer_presence(const char *begin,
                                         const char *end) const;
    sdsl::bit_vector check_kmer_presence(const std::string &sequence) const {
        return check_kmer_presence(&*sequence.begin(), &*sequence.end());
    }

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }
    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); }

    void serialize(const std::string &filename_base) const;
    void serialize(std::ostream &out) const;
    bool load(const std::string &filename_base);
    bool load(std::istream &in);

    void print_stats() const {
        std::cout << "Bloom filter parameters" << std::endl
                  << "Size:\t\t\t" << size() << " bits" << std::endl
                  << "Num hash functions:\t" << num_hash_functions() << std::endl;
    }

    static std::string file_extension() { return kExtension; }

  private:
    void call_kmers(const char *begin, const char *end,
                    const std::function<void(size_t /* position */,
                                             uint64_t /* hash1 */,
                                             uint64_t /* hash2 */)> &callback) const;

    BloomFilter filter_;
    bool canonical_mode_;
    size_t k_;
    KmerHasherType hasher_;
    static constexpr auto kExtension = ".bloom";
};


#endif // __KMER_BLOOM_FILTER_HPP__
