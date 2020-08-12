#ifndef __VECTOR_MAP_HPP__
#define __VECTOR_MAP_HPP__

#include <cstdint>

#include <tsl/ordered_map.h>

template <typename Key, typename T, typename IndexType = uint64_t>
using VectorMap = tsl::ordered_map<Key, T,
                                   std::hash<Key>, std::equal_to<Key>,
                                   std::allocator<std::pair<Key, T>>,
                                   std::vector<std::pair<Key, T>>,
                                   IndexType>;

#endif // __VECTOR_MAP_HPP__
