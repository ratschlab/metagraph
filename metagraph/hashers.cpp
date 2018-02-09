#include "hashers.hpp"
#include <cyclichash.h>


namespace annotate {

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                               const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "ORing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] | b[i];
    }
    return merged;
}

std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                const std::vector<uint64_t> &b) {

    assert(a.size() == b.size() && "ANDing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] & b[i];
    }
    return merged;
}

uint64_t popcount(const std::vector<uint64_t> &a) {
    uint64_t popcount = 0;
    for (auto value : a) {
        popcount += __builtin_popcountl(value);
    }
    return popcount;
}

bool equal(const std::vector<uint64_t> &a,
           const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "Checking different sizes");
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

bool test_bit(const std::vector<uint64_t> &a, size_t col) {
    return a[col >> 6] & (1llu << (col % 64));
}

void set_bit(std::vector<uint64_t> &a, size_t col) {
    a[col >> 6] |= 1llu << (col % 64);
}

void print(const std::vector<uint64_t> &a) {
    for (auto it = a.begin(); it != a.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

//HASH ITERATORS

HashIterator::HashIterator(const std::string &sequence,
                           size_t num_hash, size_t k)
      : seq_begin(sequence.c_str()),
        seq_cur(seq_begin),
        seq_end(sequence.c_str() + sequence.length()),
        hashes_(num_hash),
        k_(k) {
    assert(sequence.length() >= k_);
}

HashIterator::HashIterator(size_t num_hash, size_t k)
      : seq_begin(NULL),
        seq_cur(NULL),
        seq_end(NULL),
        hashes_(num_hash),
        k_(k) {}

MultiHash HashIterator::get_hash() {
    assert(hashes_.size());
    return MultiHash(operator*(), hashes_.size());
}

std::vector<MultiHash> HashIterator::generate_hashes() {
    std::vector<MultiHash> hashes;
    for (; *this != end(); operator++()) {
        hashes.emplace_back(operator*(), hashes_.size());
    }
    return hashes;
}

//MurmurHash
void MurmurHashIterator::compute_hashes() {
    for (size_t i = 0; i < hashes_.size(); ++i) {
        //use index as seed
        hashes_[i] = compute_murmur_hash(seq_cur, k_, i);
    }
}

MurmurHashIterator& MurmurHashIterator::operator++() {
    if (seq_cur + k_ <= seq_end)
        compute_hashes();
    seq_cur++;
    return *this;
}

MurmurHashIterator::MurmurHashIterator(const std::string &kmer, size_t num_hash)
  : HashIterator(num_hash, kmer.length()),
    cache_(kmer.begin(), kmer.end()),
    back_(kmer.length() - 1) {
    for (size_t i = 0; i < hashes_.size(); ++i) {
        hashes_[i] = compute_murmur_hash(kmer.c_str(), kmer.length(), i);
    }
    assert(kmer.length() == k_ || *this != end());
}

MurmurHashIterator& MurmurHashIterator::update(char next) {
    back_ = (back_ + 1) % cache_.size();
    cache_[back_] = next;
    std::string kmer =
        std::string(cache_.begin() + ((back_ + 1) % cache_.size()), cache_.end())
        + std::string(cache_.begin(), cache_.begin() + ((back_ + 1) % cache_.size()));
    for (size_t i = 0; i < hashes_.size(); ++i) {
        hashes_[i] = compute_murmur_hash(kmer.c_str(), kmer.length(), i);
    }
    return *this;
}

MurmurHashIterator& MurmurHashIterator::reverse_update(char prev) {
    cache_[back_] = prev;
    back_ = (back_ + cache_.size() - 1) % cache_.size();
    std::string kmer =
        std::string(cache_.begin() + ((back_ + 1) % cache_.size()), cache_.end())
        + std::string(cache_.begin(), cache_.begin() + ((back_ + 1) % cache_.size()));
    for (size_t i = 0; i < hashes_.size(); ++i) {
        hashes_[i] = compute_murmur_hash(kmer.c_str(), kmer.length(), i);
    }
    return *this;
}

//CyclicHashIterator
void CyclicHashIterator::init(const std::string &sequence) {
    chashers_.reserve(hashes_.size());
    for (uint32_t j = 0; j < hashes_.size(); ++j) {
        chashers_.push_back(new CyclicHash<uint64_t>(k_, j, j + 1, 64lu));
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_.back());
        for (size_t i = 0; i < k_; ++i) {
            cyclic_hash->eat(sequence[i]);
        }
        hashes_[j] = cyclic_hash->hashvalue;
    }
}

void CyclicHashIterator::compute_hashes() {
    assert(seq_cur + k_ <= seq_end);
    assert(seq_cur > seq_begin);
    for (size_t i = 0; i < hashes_.size(); ++i) {
        hashes_[i] = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->hashvalue;
    }
}

CyclicHashIterator& CyclicHashIterator::operator++() {
    if (seq_cur + k_ <= seq_end) {
        for (size_t i = 0; i < hashes_.size(); ++i) {
            reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i])->update(
                *(seq_cur - 1),
                *(seq_cur - 1 + k_)
            );
        }
        compute_hashes();
    }
    seq_cur++;
    return *this;
}

CyclicHashIterator::CyclicHashIterator(const std::string &kmer,
                                       size_t num_hash)
      : HashIterator(num_hash, kmer.length()),
    cache_(kmer.begin(), kmer.end()),
    back_(kmer.length() - 1) {
    init(kmer);
    assert(kmer.length() == k_ || *this != end());
}

CyclicHashIterator::CyclicHashIterator(const std::string &sequence,
                                       size_t num_hash, size_t k)
      : HashIterator(sequence, num_hash, k) {
    init(sequence);
    seq_cur++;
    assert(sequence.length() == k_ || *this != end());
}

CyclicHashIterator::~CyclicHashIterator() {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        delete reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
    }
}

CyclicHashIterator& CyclicHashIterator::update(char next) {
    back_ = (back_ + 1) % cache_.size();
    for (size_t i = 0; i < chashers_.size(); ++i) {
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
        cyclic_hash->update(cache_[back_], next);
        hashes_[i] = cyclic_hash->hashvalue;
    }
    cache_[back_] = next;
    return *this;
}

CyclicHashIterator& CyclicHashIterator::reverse_update(char prev) {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
        cyclic_hash->reverse_update(prev, cache_[back_]);
        hashes_[i] = cyclic_hash->hashvalue;
    }
    cache_[back_] = prev;
    back_ = (back_ + cache_.size() - 1) % cache_.size();
    return *this;
}

} // namespace annotate
