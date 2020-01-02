#include "kmer_bloom_filter.hpp"

#ifndef NDEBUG
#include "common/seq_tools/reverse_complement.hpp"
#endif
#include "common/serialization.hpp"
#include "common/algorithms.hpp"
#include "kmer/kmer_extractor.hpp"

// TODO: switch to KmerExtractor once it supports all alphabets
typedef KmerExtractorBOSS KmerDef;
typedef KmerDef::TAlphabet TAlphabet;


template <class KmerBF, typename Callback>
inline void call_kmers(const KmerBF &kmer_bloom,
                       std::string_view sequence,
                       Callback callback) {
    const auto k = kmer_bloom.get_k();
    if (sequence.size() < k)
        return;

    const auto max_encoded_val = KmerDef::alphabet.size();

    std::vector<TAlphabet> coded(sequence.size());
    std::transform(sequence.begin(), sequence.end(),
                   coded.begin(),
                   [](char c) { return KmerDef::encode(c); });
    auto invalid = utils::drag_and_mark_segments(
        coded, max_encoded_val, k
    );

    auto fwd = kmer_bloom.get_hasher();
    fwd.reset(coded.data());

    if (kmer_bloom.is_canonical_mode()) {
        std::vector<TAlphabet> rc_coded(sequence.size());
        std::transform(coded.begin(), coded.end(),
                       rc_coded.rbegin(),
                       [](TAlphabet c) { return KmerDef::complement(c); });

        auto rev = kmer_bloom.get_hasher();
        rev.reset(rc_coded.data() + rc_coded.size() - k);

        if (!invalid[k - 1]) {
            const auto &canonical = std::min(fwd, rev);
            callback(0,
                     canonical.template get_hash<0>(),
                     canonical.template get_hash<1>());
        }

        for (size_t i = k, j = coded.size() - k - 1; i < coded.size(); ++i, --j) {
            if (coded.at(i) >= max_encoded_val)
                continue;

            assert(rc_coded.at(j) < max_encoded_val);

            fwd.next(coded.at(i));
            rev.prev(rc_coded.at(j));

            if (!invalid[i]) {
                assert(i + 1 >= k);
                const auto &canonical = std::min(fwd, rev);
                callback(i + 1 - k,
                         canonical.template get_hash<0>(),
                         canonical.template get_hash<1>());
            }
        }

    } else {
        if (!invalid[k - 1])
            callback(0, fwd.template get_hash<0>(), fwd.template get_hash<1>());

        for (size_t i = k; i < coded.size(); ++i) {
            if (coded.at(i) >= max_encoded_val)
                continue;

            fwd.next(coded.at(i));

            if (!invalid[i]) {
                assert(i + 1 >= k);
                callback(i + 1 - k,
                         fwd.template get_hash<0>(),
                         fwd.template get_hash<1>());
            }
        }
    }
}

template <class KmerHasher>
void KmerBloomFilter<KmerHasher>
::add_sequence(std::string_view sequence) {
    assert(sequence.size() >= k_);

#ifndef NDEBUG
    uint64_t counter = 0;
#endif

    // TODO: insert all hashes in arbitrary order in a single call
    call_kmers(*this, sequence, [&](size_t, auto hash1, auto hash2) {
        filter_.insert(hash1, hash2);
        assert(filter_.check(hash1, hash2));
#ifndef NDEBUG
        counter++;
#endif
    });

    assert(sdsl::util::cnt_one_bits(check_kmer_presence(sequence)) == counter);
#ifndef NDEBUG
    std::string rev_comp(sequence);
    reverse_complement(rev_comp.begin(), rev_comp.end());
    assert(!canonical_mode_
            || sdsl::util::cnt_one_bits(check_kmer_presence(rev_comp)) == counter);
#endif
}

template <class KmerHasher>
sdsl::bit_vector KmerBloomFilter<KmerHasher>
::check_kmer_presence(std::string_view sequence) const {
    if (sequence.size() < k_)
        return sdsl::bit_vector();

    sdsl::bit_vector check_vec(sequence.size() - k_ + 1);
    call_kmers(*this, sequence, [&](size_t i, auto hash1, auto hash2) {
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
        const_cast<KmerHasherType&>(hasher_) = KmerHasherType(k_);

        return filter_.load(in);
    } catch (...) {
        return false;
    }
}

template class KmerBloomFilter<>;


std::function<bool()> get_missing_kmer_skipper(const KmerBloomFilter<> *bloom_filter,
                                               std::string_view sequence) {
    if (!bloom_filter)
        return []() { return false; };

    if (sequence.size() < bloom_filter->get_k())
        return []() { return true; };

    // use shared_ptr to prevent copying this vector and keep it alive for the
    // returned callback
    auto bloom_check = std::make_shared<sdsl::bit_vector>(
        bloom_filter->check_kmer_presence(sequence)
    );

    assert(bloom_check->size() + bloom_filter->get_k() - 1 == sequence.size());

    auto it = bloom_check->begin();

    // these need to be specified explicitly to ensure that they're copied
    return [it,bloom_check]() mutable {
        assert(it < bloom_check->end());
        bool in_bloom = *it;
        ++it;
        return !in_bloom;
    };
}
