#include "kmer_bloom_filter.hpp"

#include <Eigen/StdVector>

#ifndef NDEBUG
#include "common/seq_tools/reverse_complement.hpp"
#endif
#include "common/threading.hpp"
#include "kmer/kmer_extractor.hpp"
#include "utils/serialization.hpp"
#include "utils/algorithms.hpp"

// TODO: switch to KmerExtractor once it supports all alphabets
typedef KmerExtractorBOSS KmerDef;
typedef KmerDef::TAlphabet TAlphabet;

template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;


template <class KmerBF>
inline void call_kmers(const KmerBF &kmer_bloom,
                       const char *begin,
                       const char *end,
                       const std::function<void(uint64_t)> &callback) {
    const auto k = kmer_bloom.get_k();
    if (begin >= end || static_cast<size_t>(end - begin) < k)
        return;

    const auto max_encoded_val = KmerDef::alphabet.size();
    const auto &filter = kmer_bloom.get_filter();

    std::vector<TAlphabet> coded(end - begin);
    std::transform(begin, end,
                   coded.begin(),
                   [](char c) { return KmerDef::encode(c); });

    auto invalid = utils::drag_and_mark_segments(coded, max_encoded_val, k);

    auto fwd = kmer_bloom.get_hasher();
    fwd.reset(coded.data());

    if (kmer_bloom.is_canonical_mode()) {
        std::vector<TAlphabet> rc_coded(end - begin);
        std::transform(coded.begin(), coded.end(),
                       rc_coded.rbegin(),
                       [](TAlphabet c) { return KmerDef::complement(c); });

        auto rev = kmer_bloom.get_hasher();
        rev.reset(rc_coded.data() + rc_coded.size() - k);

        callback(!invalid[k - 1]
            ? filter.set_present(std::min(fwd, rev))
            : filter.ABSENCE_CHECK);

        for (size_t i = k, j = coded.size() - k - 1; i < coded.size(); ++i, --j) {
            if (coded.at(i) < max_encoded_val) {
                assert(rc_coded.at(j) < max_encoded_val);

                fwd.next(coded.at(i));
                rev.prev(rc_coded.at(j));

                assert(invalid[i] || i + 1 >= k);
                callback(!invalid[i]
                    ? filter.set_present(std::min(fwd, rev))
                    : filter.ABSENCE_CHECK);
            } else {
                callback(filter.ABSENCE_CHECK);
            }
        }

    } else {
        callback(!invalid[k - 1] ? filter.set_present(fwd) : filter.ABSENCE_CHECK);

        for (size_t i = k; i < coded.size(); ++i) {
            if (coded.at(i) < max_encoded_val) {
                fwd.next(coded.at(i));

                assert(invalid[i] || i + 1 >= k);
                callback(!invalid[i] ? filter.set_present(fwd) : filter.ABSENCE_CHECK);
            } else {
                callback(filter.ABSENCE_CHECK);
            }
        }
    }
}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::add_sequence(const char *begin, const char *end) {
    assert(end >= begin && static_cast<size_t>(end - begin) >= k_);

#ifndef NDEBUG
    size_t counter = 0;
#endif

    size_t i = 0;
    AlignedVector<uint64_t> hashes(end - begin - k_ + 1);
    call_kmers(*this, begin, end, [&](auto hash) {
        hashes[i++] = hash;

#ifndef NDEBUG
        counter += !BloomFilter::is_absent(hash);
#endif

    });

    assert(i == end - begin - k_ + 1);

    filter_.insert(hashes.data(), hashes.data() + hashes.size());

    assert(sdsl::util::cnt_one_bits(check_kmer_presence(begin, end)) == counter);

#ifndef NDEBUG
    std::string rev_comp(begin, end);
    reverse_complement(rev_comp.begin(), rev_comp.end());
    assert(!canonical_mode_
            || sdsl::util::cnt_one_bits(check_kmer_presence(rev_comp)) == counter);
#endif

}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::add_sequences(const std::function<void(const CallString&)> &generate_sequences) {
    constexpr size_t buffer_size = 1000000;
    AlignedVector<uint64_t> buffer;
    buffer.reserve(buffer_size);

    generate_sequences([&](const std::string &sequence) {
        if (sequence.size() < k_)
            return;

        if (buffer.capacity() < buffer.size() + sequence.size() - k_ + 1) {
            if (buffer.capacity() < sequence.size() - k_ + 1) {
                buffer.reserve(buffer.size() + sequence.size() - k_ + 1);
            } else {
                filter_.insert(buffer.data(), buffer.data() + buffer.size());
                buffer.clear();
            }
        }

        call_kmers(*this, sequence.c_str(), sequence.c_str() + sequence.size(),
                   [&](auto hash) { buffer.emplace_back(hash); });
    });

    filter_.insert(buffer.data(), buffer.data() + buffer.size());
}

template <class KmerHasher>
sdsl::bit_vector KmerBloomFilter<KmerHasher>
::check_kmer_presence(const char *begin, const char *end) const {
    if (begin >= end || static_cast<size_t>(end - begin) < k_)
        return sdsl::bit_vector();

    // aggregate hashes, then batch check
    size_t i = 0;
    AlignedVector<uint64_t> hashes(end - begin - k_ + 1);
    call_kmers(*this, begin, end, [&](auto hash) { hashes[i++] = hash; });
    assert(i == end - begin - k_ + 1);

    return filter_.check(hashes.data(), hashes.data() + hashes.size());
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
        const_cast<KmerHasher&>(hasher_) = KmerHasher(k_);

        return filter_.load(in);
    } catch (...) {
        return false;
    }
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
