#ifndef __VECTOR_SET_HPP__
#define __VECTOR_SET_HPP__

#include <cstdint>

#include <tsl/ordered_set.h>


template <typename Key, class Hash = std::hash<Key>, typename IndexType = uint64_t,
          class EqualTo = std::equal_to<Key>, class Allocator = std::allocator<Key>,
          class Container = std::vector<Key, Allocator>>
using VectorSet = tsl::ordered_set<Key, Hash, EqualTo, Allocator, Container, IndexType>;

#endif // __VECTOR_SET_HPP__
