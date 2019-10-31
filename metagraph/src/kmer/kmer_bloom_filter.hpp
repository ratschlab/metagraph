#ifndef __KMER_BLOOM_FILTER_HPP__
#define __KMER_BLOOM_FILTER_HPP__

#include <memory>
#include <sdsl/int_vector.hpp>


template <typename WordType, typename CharType>
class CyclicHash;

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
        for (size_t i = 0; i < h; ++i) {
            if (hashers_.at(i) < other.hashers_.at(i))
                return true;

            if (hashers_.at(i) > other.hashers_.at(i))
                return false;
        }

        return false;
    }

    bool operator==(const RollingKmerMultiHasher &other) const {
        assert(hashers_.size() == h);
        return std::equal(hashers_.begin(), hashers_.end(),
                          other.hashers_.begin(), other.hashers_.end());
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


// Bloom filter for approximate membership queries on k-mers
// Sizes are round up to the nearest power of 2 to allow for
// faster querying and insertion
class IKmerBloomFilter {
  public:
    virtual ~IKmerBloomFilter() {}

    static std::unique_ptr<IKmerBloomFilter>
    initialize_from_fpr(size_t k,
                        double false_positive_prob,
                        size_t expected_num_elements,
                        size_t max_num_hash_functions = -1,
                        bool canonical_mode = false);

    static std::unique_ptr<IKmerBloomFilter>
    initialize(size_t k,
               size_t filter_size,
               size_t num_hash_functions,
               bool canonical_mode = false);

    static std::unique_ptr<IKmerBloomFilter>
    initialize(size_t k,
               size_t filter_size = 0,
               size_t expected_num_elements = 0,
               size_t max_num_hash_functions = -1,
               bool canonical_mode = false);

    // Add the k-mers of the sequence to the Bloom filter
    virtual void add_sequence(const char *begin, const char *end) = 0;
    void add_sequence(const std::string &sequence) {
        add_sequence(&*sequence.begin(), &*sequence.end());
    }

    // Checks for k-mer presence in the Bloom filter
    virtual sdsl::bit_vector check_kmer_presence(const char *begin,
                                                 const char *end) const = 0;
    sdsl::bit_vector check_kmer_presence(const std::string &sequence) {
        return check_kmer_presence(&*sequence.begin(), &*sequence.end());
    }

    virtual bool is_canonical_mode() const = 0;

    virtual size_t get_k() const = 0;
    virtual size_t size() const = 0;
    virtual size_t num_hash_functions() const = 0;

    virtual void serialize(const std::string &filename_base) const = 0;
    virtual void serialize(std::ostream &out) const = 0;
    virtual bool load(const std::string &filename_base) = 0;
    virtual bool load(std::istream &in) = 0;

    virtual void print_stats() const {
        std::cout << "Bloom filter parameters" << std::endl
                  << "Size:\t\t\t" << size() << " bits" << std::endl
                  << "Num hash functions:\t" << num_hash_functions() << std::endl;
    }

    static std::string file_extension() { return kExtension; }

  private:
    static constexpr auto kExtension = ".bloom";
};


#endif // __KMER_BLOOM_FILTER_HPP__
