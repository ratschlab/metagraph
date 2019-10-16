#include "kmer_bloom_filter.hpp"

#include <fstream>
#include <cmath>

#include <cyclichash.h>

#include "kmer_extractor.hpp"
#include "serialization.hpp"
#include "utils.hpp"


typedef KmerExtractorBOSS KmerDef;
typedef KmerDef::TAlphabet TAlphabet;

std::vector<TAlphabet> reverse_complement(const std::vector<TAlphabet> &seq) {
    return KmerDef::reverse_complement(seq);
}

TAlphabet code(char c) { return KmerDef::encode(c); }

template <class StringIt>
std::vector<TAlphabet> encode(StringIt begin, StringIt end) {
    assert(end >= begin);

    std::vector<TAlphabet> encoded(end - begin, 0);
    std::transform(begin, end, encoded.begin(), code);

    return encoded;
}

std::string decode(const std::vector<TAlphabet> &encoded) {
    return KmerDef::decode(encoded);
}


template <class T>
class BloomFilterWrapper {
  public:
    virtual ~BloomFilterWrapper() {}

    virtual void insert(const T&) = 0;
    virtual bool check(const T&) const = 0;

    virtual void serialize(ostream &out) const = 0;
    virtual bool load(istream &in) = 0;

    virtual size_t size() const = 0;
    virtual size_t num_hash_functions() const = 0;
};

template <class T>
class BloomFilter : public BloomFilterWrapper<T> {
  public:
    BloomFilter(double false_positive_prob,
               size_t expected_num_elements,
               size_t max_num_hash_functions)
          : BloomFilter(optim_size(false_positive_prob, expected_num_elements),
                        expected_num_elements,
                        max_num_hash_functions) {}

    BloomFilter(size_t filter_size, size_t num_hash_functions)
          : filter_(size_t(1) << utils::code_length(filter_size)),
            num_hash_functions_(num_hash_functions),
            filter_mask_(filter_.size() - 1) {}

    BloomFilter(size_t filter_size,
                size_t expected_num_elements,
                size_t max_num_hash_functions)
          : filter_(size_t(1) << utils::code_length(filter_size)),
            num_hash_functions_(std::min(max_num_hash_functions,
                                         optim_h(filter_.size(), expected_num_elements))),
            filter_mask_(filter_.size() - 1) {}

    void insert(const T &obj) {
        assert(filter_.size());
        size_t d = obj.template get_hash<1>();
        // Kirsch, A., & Mitzenmacher, M. (2006, September).
        // Less hashing, same performance: building a better bloom filter.
        // In European Symposium on Algorithms (pp. 456-467). Springer, Berlin, Heidelberg.
        for (size_t i = 0, h = obj.template get_hash<0>(); i < num_hash_functions_; ++i, h += d) {
            // This only works if the filter size is a power of 2
            // TODO: do some locking here to make this multithreaded
            filter_[h & filter_mask_] = 1;
        }

        assert(check(obj));
    }

    bool check(const T &obj) const {
        assert(filter_.size());
        size_t d = obj.template get_hash<1>();
        for (size_t i = 0, h = obj.template get_hash<0>(); i < num_hash_functions_; ++i, h += d) {
            // This only works if the filter size is a power of 2
            if (!filter_[h & filter_mask_])
                return false;
        }

        return true;
    }

    void serialize(ostream &out) const {
        filter_.serialize(out);
        serialize_number(out, num_hash_functions_);
        serialize_number(out, filter_mask_);
    }

    bool load(istream &in) {
        try {
            filter_.load(in);
            num_hash_functions_ = load_number(in);
            filter_mask_ = load_number(in);
            return true;
        } catch (...) {
            return false;
        }
    }

    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return num_hash_functions_; }

  private:
    constexpr static size_t optim_size(double false_positive_prob, size_t expected_num_elements) {
        return -std::log2(false_positive_prob) * expected_num_elements / M_LN2 / M_LN2;
    }

    constexpr static size_t optim_h(double false_positive_prob) {
        return std::ceil(-std::log2(false_positive_prob));
    }

    constexpr static size_t optim_h(size_t filter_size, size_t expected_num_elements) {
        return std::ceil(M_LN2 * filter_size / expected_num_elements);
    }

    sdsl::bit_vector filter_;
    size_t num_hash_functions_;
    size_t filter_mask_;
};

template <class KmerHasher = RollingKmerMultiHasher<2, KmerDef::TAlphabet>,
          class BloomFilter = BloomFilter<KmerHasher>>
class KmerBloomFilter : public IKmerBloomFilter {
  public:
    typedef KmerHasher KmerHasherType;

    // Initialize a k-mer Bloom filter
    explicit KmerBloomFilter(size_t k,
                             double false_positive_prob,
                             size_t expected_num_elements,
                             size_t max_num_hash_functions = -1,
                             bool canonical_mode = false,
                             uint64_t seed = 0x100000000)
          : filter_(false_positive_prob, expected_num_elements, max_num_hash_functions),
            canonical_mode_(canonical_mode),
            seed_(seed),
            k_(k),
            hasher_(k_, seed_) {}

    explicit KmerBloomFilter(size_t k,
                             size_t filter_size,
                             size_t num_hash_functions,
                             bool canonical_mode = false,
                             uint64_t seed = 0x100000000)
          : filter_(filter_size, num_hash_functions),
            canonical_mode_(canonical_mode),
            seed_(seed),
            k_(k),
            hasher_(k_, seed_) {}

    explicit KmerBloomFilter(size_t k,
                             size_t filter_size,
                             size_t expected_num_elements,
                             size_t max_num_hash_functions = -1,
                             bool canonical_mode = false,
                             uint64_t seed = 0x100000000)
          : filter_(filter_size, expected_num_elements, max_num_hash_functions),
            canonical_mode_(canonical_mode),
            seed_(seed),
            k_(k),
            hasher_(k_, seed_) {}

    void add_sequence(const char *begin, const char *end) {
        assert(begin + k_ > begin && begin + k_ <= end);

        #ifndef NDEBUG
        uint64_t counter = 0;
        #endif

        call_kmers(begin, end, [&](auto, const auto &kmer) {
            filter_.insert(kmer);
            assert(filter_.check(kmer));
            #ifndef NDEBUG
            counter++;
            #endif
        });

        assert(sdsl::util::cnt_one_bits(check_kmer_presence(begin, end)) == counter);
        assert(!canonical_mode_
                || sdsl::util::cnt_one_bits(IKmerBloomFilter::check_kmer_presence(
                       decode(reverse_complement(encode(begin, end)))
                   )) == counter);
    }

    sdsl::bit_vector check_kmer_presence(const char *begin, const char *end) const {
        if (begin + k_ > end)
            return sdsl::bit_vector();

        sdsl::bit_vector check(end - begin - k_ + 1);
        call_kmers(begin, end, [&](auto i, const auto &kmer) {
            assert(i < check.size());

            if (filter_.check(kmer))
                check[i] = true;
        });

        return check;
    }

    bool is_canonical_mode() const { return canonical_mode_; }

    size_t get_k() const { return k_; }

    size_t size() const { return filter_.size(); }
    size_t num_hash_functions() const { return filter_.num_hash_functions(); };

    void serialize(std::ostream &out) const {
        if (!out.good())
            throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

        serialize_number(out, k_);
        serialize_number(out, seed_);
        serialize_number(out, canonical_mode_);
        filter_.serialize(out);
    }

    bool load(std::istream &in) {
        if (!in.good())
            return false;

        try {
            k_ = load_number(in);
            seed_ = load_number(in);
            canonical_mode_ = load_number(in);
            hasher_ = KmerHasherType(k_, seed_);

            return filter_.load(in);
        } catch (...) {
            return false;
        }
    }

    void serialize(const std::string &filename_base) const {
        std::ofstream out(utils::remove_suffix(filename_base, file_extension())
                              + file_extension(),
                          std::ios::binary);
        serialize(out);
    }

    bool load(const std::string &filename_base) {
        std::ifstream in(utils::remove_suffix(filename_base, file_extension())
                             + file_extension(),
                         std::ios::binary);
        return load(in);
    }

  private:
    void call_kmers(const char *begin, const char *end,
                    const std::function<void(size_t i /* position */,
                                             const KmerHasherType&)> &callback) const {
        auto fwd = hasher_;

        if (is_canonical_mode()) {
            auto rev = hasher_;

            if (begin + k_ > end || begin + k_ < begin)
                return;

            std::vector<TAlphabet> coded(end - begin + k_, KmerHasherType::MAXVAL);
            std::transform(begin, end, coded.begin() + k_, code);
            auto rc_coded = reverse_complement(coded);

            size_t chars = 0;
            for (TAlphabet *i = coded.data(),
                           *j = rc_coded.data() + rc_coded.size() - 1;
                            i < &*coded.end();
                            ++i, --j) {
                assert(j >= rc_coded.data());
                if (*i >= KmerDef::alphabet.size()) {
                    chars = 0;
                    continue;
                }
                assert(*j < KmerDef::alphabet.size());

                fwd.shift_left(i, i - k_);
                rev.shift_right(j, j + k_);

                if (++chars >= k_)
                    callback(i - coded.data() + 1 - (k_ << 1), std::min(fwd, rev));
            }

        } else {
            if (begin + k_ > end || begin + k_ < begin)
                return;

            std::vector<TAlphabet> coded(end - begin + k_, KmerHasherType::MAXVAL);
            std::transform(begin, end, coded.begin() + k_, code);

            size_t chars = 0;
            for (TAlphabet *i = coded.data(); i < &*coded.end(); ++i) {
                if (*i >= KmerDef::alphabet.size()) {
                    chars = 0;
                    continue;
                }

                fwd.shift_left(i, i - k_);

                if (++chars >= k_)
                    callback(i - coded.data() + 1 - (k_ << 1), fwd);
            }
        }
    }

    BloomFilter filter_;
    bool canonical_mode_;
    size_t seed_;
    size_t k_;
    KmerHasherType hasher_;
};

template <class KmerHasher>
std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter
::initialize(size_t k,
             double false_positive_prob,
             size_t expected_num_elements,
             size_t max_num_hash_functions,
             bool canonical_mode,
             uint64_t seed) {
    auto bloom_filter = std::make_unique<KmerBloomFilter<KmerHasher>>(
        k,
        false_positive_prob,
        expected_num_elements,
        max_num_hash_functions,
        canonical_mode,
        seed
    );

    if (utils::get_verbose())
        bloom_filter->print_stats();

    return bloom_filter;
}

template std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter
::initialize<RollingKmerMultiHasher<2>>(size_t, double, size_t, size_t, bool, uint64_t);

template <class KmerHasher>
std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter::
initialize(size_t k,
           size_t filter_size,
           size_t num_hash_functions,
           bool canonical_mode,
           uint64_t seed) {
    auto bloom_filter = std::make_unique<KmerBloomFilter<KmerHasher>>(
        k,
        filter_size,
        num_hash_functions,
        canonical_mode,
        seed
    );

    if (utils::get_verbose() && filter_size)
        bloom_filter->print_stats();

    return bloom_filter;
}

template std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter
::initialize<RollingKmerMultiHasher<2>>(size_t, size_t, size_t, bool, uint64_t);

template <class KmerHasher>
std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter::
initialize(size_t k,
           size_t filter_size,
           size_t expected_num_elements,
           size_t max_num_hash_functions,
           bool canonical_mode,
           uint64_t seed) {
    auto bloom_filter = std::make_unique<KmerBloomFilter<KmerHasher>>(
        k,
        filter_size,
        expected_num_elements,
        max_num_hash_functions,
        canonical_mode,
        seed
    );

    if (utils::get_verbose() && filter_size)
        bloom_filter->print_stats();

    return bloom_filter;
}

template std::unique_ptr<IKmerBloomFilter> IKmerBloomFilter
::initialize<RollingKmerMultiHasher<2>>(size_t, size_t, size_t, size_t, bool, uint64_t);


void IKmerBloomFilter::print_stats() const {
    std::cout << "Bloom filter parameters" << std::endl
              << "Size:\t\t\t" << size() << " bits" << std::endl
              << "Num hash functions:\t" << num_hash_functions() << std::endl;
}
