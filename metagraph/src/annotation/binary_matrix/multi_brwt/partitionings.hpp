#ifndef __PARTITIONINGS_HPP__
#define __PARTITIONINGS_HPP__

#include <vector>

#include "common/vectors/bit_vector.hpp"

// Clustering columns for Multi-BRWT

// input: columns
// output: partition -- a set of column pairs greedily matched
std::vector<std::vector<uint64_t>>
greedy_matching(const std::vector<const bit_vector *> &columns,
                size_t num_threads = 1);

#endif // __PARTITIONINGS_HPP__
