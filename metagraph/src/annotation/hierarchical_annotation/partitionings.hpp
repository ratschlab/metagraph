#ifndef __PARTITIONINGS_HPP__
#define __PARTITIONINGS_HPP__

#include "BRWT_builders.hpp"


// Partitionings for BRWT

// input: columns
// output: partition, for instance -- a set of column pairs
BRWTBottomUpBuilder::Partition
parallel_binary_grouping_greedy(const BRWTBottomUpBuilder::VectorsPtr &columns,
                                size_t num_threads);

BRWTBottomUpBuilder::Partition
binary_grouping_greedy(const BRWTBottomUpBuilder::VectorsPtr &columns);

std::function<BRWTBottomUpBuilder::Partition(const BRWTBottomUpBuilder::VectorsPtr &)>
get_parallel_binary_grouping_greedy(size_t num_threads);


#endif // __PARTITIONINGS_HPP__
