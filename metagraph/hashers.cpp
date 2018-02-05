#include "hashers.hpp"

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

bool equal(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "Checking different sizes");
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

bool test_bit(const std::vector<uint64_t> &a, const size_t col) {
    return (a[col >> 6] | (1llu << (col % 64))) == a[col >> 6];
}

void set_bit(std::vector<uint64_t> &a, const size_t col) {
    a[col >> 6] |= 1llu << (col % 64);
}

template <class HashIt>
std::vector<MultiHash> hash(HashIt&& hash_it, const size_t num_hash) {
    std::vector<MultiHash> hashes;
    for (;hash_it != hash_it.end(); ++hash_it) {
        hashes.emplace_back(*hash_it, num_hash);
    }
    return hashes;
}

std::vector<MultiHash> hash_murmur(
        const std::string &sequence,
        const size_t num_hash,
        const size_t k) {
    auto hashes = hash(HashIterator(sequence, num_hash, k), num_hash);
    assert(hashes.size() == sequence.length() - k + 1);
    return hashes;
}

std::vector<size_t> annotate(const MultiHash &multihash, const size_t max_size) {
    std::vector<size_t> annotation((max_size >> 6) + 1);
    for (auto it = multihash.begin(); it != multihash.end(); ++it) {
        set_bit(annotation, *it % max_size);
    }
    return annotation;
}

template std::vector<MultiHash> hash(HashIterator&&, const size_t);
template std::vector<MultiHash> hash(HashIterator&, const size_t);
template std::vector<MultiHash> hash(ntHashIterator&&, const size_t);
template std::vector<MultiHash> hash(ntHashIterator&, const size_t);

} // namespace annotate
