#include "kmer_bloom_filter.hpp"

#ifndef NDEBUG
#include "common/seq_tools/reverse_complement.hpp"
#endif
#include "utils/serialization.hpp"
#include "utils/algorithms.hpp"
#include "kmer/kmer_extractor.hpp"

// TODO: switch to KmerExtractor once it supports all alphabets
typedef KmerExtractorBOSS KmerDef;
typedef KmerDef::TAlphabet TAlphabet;


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
#ifndef NDEBUG
    std::string rev_comp(begin, end);
    reverse_complement(rev_comp.begin(), rev_comp.end());
    assert(!canonical_mode_
            || sdsl::util::cnt_one_bits(check_kmer_presence(rev_comp)) == counter);
#endif
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

    const auto max_encoded_val = KmerDef::alphabet.size();

    std::vector<TAlphabet> coded(end - begin);
    std::transform(begin, end,
                   coded.begin(),
                   [](char c) { return KmerDef::encode(c); });
    auto invalid = utils::drag_and_mark_segments(
        coded, max_encoded_val, k_
    );

    auto fwd = hasher_;
    fwd.reset(coded.data());

    if (is_canonical_mode()) {
        std::vector<TAlphabet> rc_coded(end - begin);
        std::transform(coded.begin(), coded.end(),
                       rc_coded.rbegin(),
                       [](TAlphabet c) { return KmerDef::complement(c); });

        auto rev = hasher_;
        rev.reset(rc_coded.data() + rc_coded.size() - k_);

        if (LIKELY(!invalid[k_ - 1])) {
            const auto &canonical = std::min(fwd, rev);
            callback(0,
                     canonical.template get_hash<0>(),
                     canonical.template get_hash<1>());
        }

        for (size_t i = k_, j = coded.size() - k_ - 1; i < coded.size(); ++i, --j) {
            if (UNLIKELY(coded.at(i) >= max_encoded_val))
                continue;

            assert(rc_coded.at(j) < max_encoded_val);

            fwd.next(coded.at(i));
            rev.prev(rc_coded.at(j));

            if (LIKELY(!invalid[i])) {
                assert(i + 1 >= k_);
                const auto &canonical = std::min(fwd, rev);
                callback(i + 1 - k_,
                         canonical.template get_hash<0>(),
                         canonical.template get_hash<1>());
            }
        }

    } else {
        if (LIKELY(!invalid[k_ - 1]))
            callback(0, fwd.template get_hash<0>(), fwd.template get_hash<1>());

        for (size_t i = k_; i < coded.size(); ++i) {
            if (UNLIKELY(coded.at(i) >= max_encoded_val))
                continue;

            fwd.next(coded.at(i));

            if (LIKELY(!invalid[i])) {
                assert(i + 1 >= k_);
                callback(i + 1 - k_,
                         fwd.template get_hash<0>(),
                         fwd.template get_hash<1>());
            }
        }
    }
}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::print_stats() const {
    std::cout << "Bloom filter parameters" << std::endl
              << "Size:\t\t\t" << size() << " bits" << std::endl
              << "Num hash functions:\t" << num_hash_functions() << std::endl;
}

template class KmerBloomFilter<>;


std::function<bool()> get_missing_kmer_skipper(const KmerBloomFilter<> *bloom_filter,
                                               const char *begin,
                                               const char *end) {
    if (!bloom_filter)
        return []() { return false; };

    if (begin + bloom_filter->get_k() > end)
        return []() { return true; };

    // use shared_ptr to prevent copying this vector and keep it alive for the
    // returned callback
    auto bloom_check = std::make_shared<sdsl::bit_vector>(
        bloom_filter->check_kmer_presence(begin, end)
    );

    assert(begin + bloom_check->size() == end - bloom_filter->get_k() + 1);

    auto it = bloom_check->begin();

    // these need to be specified explicitly to ensure that they're copied
    return [it,bloom_check]() mutable {
        assert(it < bloom_check->end());
        bool in_bloom = *it;
        ++it;
        return !in_bloom;
    };
}
