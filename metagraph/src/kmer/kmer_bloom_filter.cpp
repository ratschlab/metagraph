#include "kmer_bloom_filter.hpp"


#include "serialization.hpp"
#include "utils.hpp"


template <class StringIt>
std::vector<TAlphabet>
encode_with_sentinel_prefix(StringIt begin, StringIt end,
                            size_t k,
                            TAlphabet default_value = 0) {
    assert(end >= begin);

    std::vector<TAlphabet> encoded(end - begin + k, default_value);
    std::transform(begin, end, encoded.begin() + k,
                   [](char c) { return KmerDef::encode(c); });

    return encoded;
}

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
        // This only works if the filter size is a power of 2
        // TODO: do some locking here to make this multithreaded
        const auto hash = hash1 + i * hash2;
        filter_[hash - hash / size * size] = true;
    }

    assert(check(hash1, hash2));
}

bool BloomFilter::check(uint64_t hash1, uint64_t hash2) const {
    const auto size = filter_.size();
    if (!size)
        return true;

    for (size_t i = 0; i < num_hash_functions_; ++i) {
        // This only works if the filter size is a power of 2
        const auto hash = hash1 + i * hash2;
        if (!filter_[hash - hash / size * size])
            return false;
    }

    return true;
}

void BloomFilter::serialize(ostream &out) const {
    filter_.serialize(out);
    serialize_number(out, num_hash_functions_);
}

bool BloomFilter::load(istream &in) {
    try {
        filter_.load(in);
        num_hash_functions_ = load_number(in);
        return true;
    } catch (...) {
        return false;
    }
}


template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::add_sequence(const char *begin, const char *end) {
    assert(begin + k_ > begin && begin + k_ <= end);

    #ifndef NDEBUG
    uint64_t counter = 0;
    #endif

    call_kmers(begin, end, [&](auto, auto hash1, auto hash2) {
        filter_.insert(hash1, hash2);
        assert(filter_.check(hash1, hash2));
        #ifndef NDEBUG
        counter++;
        #endif
    });

    assert(sdsl::util::cnt_one_bits(check_kmer_presence(begin, end)) == counter);
    assert(!canonical_mode_
            || sdsl::util::cnt_one_bits(check_kmer_presence(
                   KmerDef::decode(KmerDef::reverse_complement(KmerDef::encode(std::string(begin, end))))
               )) == counter);
}

template <class KmerHasher>
sdsl::bit_vector KmerBloomFilter<KmerHasher>
::check_kmer_presence(const char *begin, const char *end) const {
    if (begin + k_ > end)
        return sdsl::bit_vector();

    sdsl::bit_vector check_vec(end - begin - k_ + 1);
    call_kmers(begin, end, [&](auto i, auto hash1, auto hash2) {
        assert(i < check_vec.size());

        if (filter_.check(hash1, hash2))
            check_vec[i] = true;
    });

    return check_vec;
}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Error: trying to dump graph to a bad stream");

    serialize_number(out, k_);
    serialize_number(out, canonical_mode_);
    filter_.serialize(out);
}

template <class KmerHasher>
bool KmerBloomFilter<KmerHasher>
::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        k_ = load_number(in);
        canonical_mode_ = load_number(in);
        hasher_ = KmerHasherType(k_);

        return filter_.load(in);
    } catch (...) {
        return false;
    }
}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::call_kmers(const char *begin, const char *end,
             const std::function<void(size_t /* position */,
                                      uint64_t /* hash1 */,
                                      uint64_t /* hash2 */)> &callback) const {
    auto fwd = hasher_;

    if (is_canonical_mode()) {
        auto rev = hasher_;

        if (begin + k_ > end || begin + k_ < begin)
            return;

        auto coded = encode_with_sentinel_prefix(begin, end, k_, fwd.get_default_val());
        auto rc_coded = KmerDef::reverse_complement(coded);

        size_t chars = 0;
        for (size_t i = 0, j = coded.size() - 1; i < coded.size(); ++i, --j) {
            if (coded.at(i) >= KmerDef::alphabet.size()) {
                chars = 0;
                continue;
            }

            assert(rc_coded.at(j) < KmerDef::alphabet.size());

            assert(i >= k_);
            fwd.shift_left(coded.at(i), coded.at(i - k_));

            assert(j + k_ < coded.size());
            rev.shift_right(rc_coded.at(j), rc_coded.at(j + k_));

            if (++chars >= k_) {
                const auto &canonical = std::min(fwd, rev);
                callback(i + 1 - (k_ << 1),
                         canonical.template get_hash<0>(),
                         canonical.template get_hash<1>());
            }
        }

    } else {
        if (begin + k_ > end || begin + k_ < begin)
            return;

        auto coded = encode_with_sentinel_prefix(begin, end, k_, fwd.get_default_val());

        size_t chars = 0;
        for (size_t i = 0; i < coded.size(); ++i) {
            if (coded.at(i) >= KmerDef::alphabet.size()) {
                chars = 0;
                continue;
            }

            assert(i >= k_);
            fwd.shift_left(coded.at(i), coded.at(i - k_));

            if (++chars >= k_)
                callback(i + 1 - (k_ << 1),
                         fwd.template get_hash<0>(),
                         fwd.template get_hash<1>());
        }
    }
}

template class KmerBloomFilter<RollingKmerMultiHasher<2>>;
