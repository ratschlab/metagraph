#ifndef __VECTOR_SET_HPP__
#define __VECTOR_SET_HPP__

#include <cstdint>

#include <tsl/ordered_set.h>

template <typename T,
          class Hash = std::hash<T>,
          typename IndexType = uint64_t,
          class EqualTo = std::equal_to<T>,
          class Allocator = std::allocator<T>,
          class Container = std::vector<T, Allocator>>
using VectorSet = tsl::ordered_set<T, Hash, EqualTo, Allocator, Container, IndexType>;

template <typename T,
          class Hash = std::hash<T>,
          typename IndexType = uint64_t,
          class EqualTo = std::equal_to<T>,
          class Allocator = std::allocator<T>,
          class Container = std::vector<T, Allocator>>
Container to_vector(VectorSet<T, Hash, IndexType, EqualTo, Allocator, Container> &&set) {
    Container vector;
    vector.swap(const_cast<Container&>(set.values_container()));
    return vector;
}

#endif // __VECTOR_SET_HPP__
