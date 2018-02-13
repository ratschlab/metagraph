#include "hashers.hpp"

#include "cyclichash.h"


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

bool test_bit(const std::vector<uint64_t> &a, size_t i) {
    return a[i >> 6] & (1llu << (i % 64));
}

void set_bit(std::vector<uint64_t> &a, size_t i) {
    a[i >> 6] |= 1llu << (i % 64);
}

void print(const std::vector<uint64_t> &a) {
    for (auto it = a.begin(); it != a.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}


//CyclicHashIterator
CyclicHashIterator::CyclicHashIterator(const char *begin, const char *end,
                                       size_t num_hash, size_t k)
      : hashes_(num_hash),
        k_(k),
        seq_cur(begin), seq_end(end),
        cache_(begin, end),
        back_(k - 1) {

    assert(k_);
    assert(end - begin >= static_cast<int>(k_));

    init(begin);
    seq_cur++;

    assert(end - begin == static_cast<int>(k_) || !is_end());
}

CyclicHashIterator::~CyclicHashIterator() {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        delete reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
    }
}

void CyclicHashIterator::init(const char *data) {
    chashers_.reserve(hashes_.size());
    for (uint32_t j = 0; j < hashes_.size(); ++j) {
        chashers_.push_back(new CyclicHash<uint64_t>(k_, j, j + 1, 64lu));
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_.back());
        for (size_t i = 0; i < k_; ++i) {
            cyclic_hash->eat(data[i]);
        }
        hashes_[j] = cyclic_hash->hashvalue;
    }
}

void CyclicHashIterator::compute_hashes() {
    assert(seq_cur + k_ <= seq_end);
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
