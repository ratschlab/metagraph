#ifndef __CLUSTERING_HPP__
#define __CLUSTERING_HPP__

#include <vector>

#include <Eigen/Dense>

#include "common/vectors/bit_vector.hpp"


namespace mtg {
namespace annot {
namespace binmat {

typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, 4> LinkageMatrix;

// Clustering of columns for Multi-BRWT

// input: columns
// output: partition -- a set of column pairs greedily matched
std::vector<std::vector<uint64_t>>
greedy_matching(const std::vector<sdsl::bit_vector> &columns,
                size_t num_threads = 1);

// Format resembling the Z matrix from scipy.cluster.hierarchy.linkage
// result: (x - 1) by 4 matrix
// Points result[i, 0] and result[i, 1] are merged into result[i, 3]
// result[i, 2] = dist(result[i, 0], result[i, 1])
LinkageMatrix
agglomerative_greedy_linkage(std::vector<sdsl::bit_vector>&& columns,
                             size_t num_threads = 1);

// Merges points in their original order
LinkageMatrix agglomerative_linkage_trivial(size_t num_columns);

std::vector<uint64_t>
sample_row_indexes(uint64_t num_rows, uint64_t size, int seed = 1);

std::vector<sdsl::bit_vector>
random_submatrix(const std::vector<const bit_vector *> &columns,
                 uint64_t num_rows_sampled,
                 size_t num_threads = 1, int seed = 1);

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __CLUSTERING_HPP__
