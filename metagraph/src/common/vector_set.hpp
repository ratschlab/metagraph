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

// Use when the set is no longer needed and you want to extract the underlying container
template <typename T, class Hash, typename IndexType, class EqualTo, class Allocator, class Container>
Container to_vector(VectorSet<T, Hash, IndexType, EqualTo, Allocator, Container> &&set) {
    Container container;
    container.swap(const_cast<Container &>(set.values_container()));
    return container;
}

#endif // __VECTOR_SET_HPP__
