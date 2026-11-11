#ifndef __ANNOTATION_REMAP_HPP__
#define __ANNOTATION_REMAP_HPP__

#include <filesystem>
#include <memory>
#include <string>

#include <sdsl/int_vector.hpp>


namespace mtg {

namespace graph {
class DeBruijnGraph;
} // namespace graph

namespace annot {

template <typename Label>
class MultiLabelAnnotation;

/**
 * Match the annotation rows of `source` to the rows of `target` by k-mer.
 * Both graphs must be in BASIC mode and contain the same k-mers. Source nodes
 * that do not represent a k-mer, such as the dummy nodes of a succinct graph,
 * are dropped.
 *
 * Returns a vector with source.max_index() entries. Entry r holds the target row
 * + 1 for source row r, or 0 if the row is dropped.
 */
sdsl::int_vector<> map_rows_to_graph(const graph::DeBruijnGraph &source,
                                     const graph::DeBruijnGraph &target,
                                     size_t num_threads = 1);

/**
 * Same, but writes the vector to `filename` in the sdsl::int_vector<> format
 * without holding it in memory. Uses temporary files in `tmp_dir` of about
 * twice the size of the result.
 */
void map_rows_to_graph(const graph::DeBruijnGraph &source,
                       const graph::DeBruijnGraph &target,
                       const std::string &filename,
                       const std::filesystem::path &tmp_dir,
                       size_t num_threads = 1);

/**
 * Check a saved row map against both graphs. Returns false if the row counts
 * differ, or if a kept row maps to another target row than its k-mer. The
 * k-mer check covers num_samples evenly spaced rows and skips dropped ones.
 */
bool check_row_map(const graph::DeBruijnGraph &source,
                   const graph::DeBruijnGraph &target,
                   const sdsl::int_vector<> &row_map,
                   uint64_t num_samples = 10000);

/**
 * Same, but reads the row map from a file without holding it in memory.
 */
bool check_row_map(const graph::DeBruijnGraph &source,
                   const graph::DeBruijnGraph &target,
                   const std::string &row_map_filename,
                   uint64_t num_samples = 10000);

/**
 * Move the rows of a BRWT or RowDiff<BRWT> annotation built for `source` to a graph
 * with target_num_rows rows, using a row map from map_rows_to_graph. The BRWT
 * topology and column order are kept. For RowDiff, every kept row keeps its
 * successor in `source`, which is stored as a successor bitmap for the target graph.
 */
std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const graph::DeBruijnGraph &source,
                      const sdsl::int_vector<> &row_map,
                      uint64_t target_num_rows,
                      size_t num_threads = 1);

/**
 * Same, but reads the row map from a file without holding it in memory.
 */
std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const graph::DeBruijnGraph &source,
                      const std::string &row_map_filename,
                      uint64_t target_num_rows,
                      size_t num_threads = 1);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_REMAP_HPP__
