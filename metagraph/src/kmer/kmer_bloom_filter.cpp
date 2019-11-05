#include "kmer_bloom_filter.hpp"

#include "serialization.hpp"
#include "utils.hpp"

typedef KmerDef::TAlphabet TAlphabet;


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
    assert(end >= begin && static_cast<size_t>(end - begin) >= k_);

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
    if (begin >= end || static_cast<size_t>(end - begin) < k_)
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
    if (begin >= end || static_cast<size_t>(end - begin) < k_)
        return;

    std::vector<TAlphabet> coded(end - begin);
    std::transform(begin, end,
                   coded.begin(),
                   [](char c) { return KmerDef::encode(c); });

    // determine the number of invalid characters in the first k-mer
    size_t chars = std::find_if(coded.rend() - k_, coded.rend(),
                                [](TAlphabet c) { return c >= KmerDef::alphabet.size(); })
                       - (coded.rend() - k_);
    assert(chars <= k_);

    auto fwd = hasher_;
    fwd.reset(coded.data());

    if (is_canonical_mode()) {
        std::vector<TAlphabet> rc_coded(end - begin);
        std::transform(coded.begin(), coded.end(),
                       rc_coded.begin(),
                       [](TAlphabet c) { return KmerDef::complement(c); });

        auto rev = hasher_;
        rev.reverse_reset(rc_coded.data());

        if (chars == k_) {
            const auto &canonical = std::min(fwd, rev);
            callback(0,
                     canonical.template get_hash<0>(),
                     canonical.template get_hash<1>());
        }

        for (size_t i = k_, j = k_; i < coded.size(); ++i, ++j) {
            if (coded.at(i) >= KmerDef::alphabet.size()) {
                chars = 0;
                continue;
            }

            assert(rc_coded.at(j) < KmerDef::alphabet.size());

            assert(i >= k_);
            fwd.next(coded.at(i));
            rev.prev(rc_coded.at(j));

            if (++chars >= k_) {
                const auto &canonical = std::min(fwd, rev);
                callback(i - k_ + 1,
                         canonical.template get_hash<0>(),
                         canonical.template get_hash<1>());
            }
        }

    } else {
        if (chars == k_)
            callback(0, fwd.template get_hash<0>(), fwd.template get_hash<1>());

        for (size_t i = k_; i < coded.size(); ++i) {
            if (coded.at(i) >= KmerDef::alphabet.size()) {
                chars = 0;
                continue;
            }

            fwd.next(coded.at(i));

            if (++chars >= k_)
                callback(i - k_ + 1,
                         fwd.template get_hash<0>(),
                         fwd.template get_hash<1>());
        }
    }
}

template class KmerBloomFilter<RollingKmerMultiHasher<2>>;
