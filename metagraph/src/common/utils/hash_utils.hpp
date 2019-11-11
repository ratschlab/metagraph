#ifndef __HASH_UTILS_HPP__
#define __HASH_UTILS_HPP__

#include <utility>
#include <bitset>


namespace utils {

template <typename T>
struct Hash {
    size_t operator()(const T &x) const {
        return hasher(reinterpret_cast<const std::bitset<sizeof(T) * 8>&>(x));
    }

    std::hash<std::bitset<sizeof(T) * 8>> hasher;
};

struct VectorHash {
    template <class Vector>
    std::size_t operator()(const Vector &vector) const {
        uint64_t hash = 0;
        for (uint64_t value : vector) {
            hash ^= value + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return static_cast<std::size_t>(hash);
    }
};

} // namespace utils

#endif // __HASH_UTILS_HPP__
