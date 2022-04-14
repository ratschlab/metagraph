#ifndef __ALIGNED_VECTOR_HPP__
#define __ALIGNED_VECTOR_HPP__

#include <vector>

#include <boost/align/aligned_allocator.hpp>


template <typename T>
using AlignedVector = std::vector<T, boost::alignment::aligned_allocator<T>>;


#endif // __ALIGNED_VECTOR_HPP__
