#ifndef __ALIGNED_VECTOR_HPP__
#define __ALIGNED_VECTOR_HPP__

#include <vector>

#include <Eigen/StdVector>


template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;


#endif // __ALIGNED_VECTOR_HPP__
