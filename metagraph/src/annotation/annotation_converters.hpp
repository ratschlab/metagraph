#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <vector>
#include <filesystem>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg {
namespace annot {

template <typename Label>
class RowCompressed;

template <typename LabelType>
class MultiLabelEncoded;


template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert(RowCompressed<Label>&& annotation);


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
std::unique_ptr<StaticAnnotation> convert(ColumnCompressed<Label>&& annotation);

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation> convert(const std::string &filename);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_simple_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);

template <class StaticAnnotation, typename Label>
typename std::unique_ptr<StaticAnnotation>
convert_to_greedy_BRWT(ColumnCompressed<Label>&& annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);

template <class StaticAnnotation>
typename std::unique_ptr<StaticAnnotation>
convert_to_BRWT(const std::vector<std::string> &annotation_files,
                const std::string &linkage_matrix_file,
                size_t num_parallel_nodes = 1,
                size_t num_threads = 1,
                const std::filesystem::path &tmp_dir = "");

template <class StaticAnnotation>
void relax_BRWT(StaticAnnotation *annotation,
                size_t relax_max_arity,
                size_t num_threads = 1);

template <class StaticAnnotation>
typename std::unique_ptr<StaticAnnotation>
convert_to_RbBRWT(const std::vector<std::string> &annotation_files,
                  size_t max_brwt_arity);

template <class ToAnnotation, typename Label>
void merge(std::vector<std::unique_ptr<MultiLabelEncoded<Label>>>&& annotators,
           const std::vector<std::string> &filenames,
           const std::string &outfile);

template <typename Label>
void convert_to_row_annotator(const ColumnCompressed<Label> &annotator,
                              const std::string &outfbase,
                              size_t num_threads = 1);

template <typename Label>
void convert_to_row_annotator(const ColumnCompressed<Label> &annotator,
                              RowCompressed<Label> *target,
                              size_t num_threads = 1);

template <typename Label>
std::unique_ptr<ColumnDiffAnnotator>
convert_to_column_diff(const graph::DBGSuccinct &graph,
                       const ColumnCompressed<Label> &annotation,
                       const std::string &outfbase,
                       uint32_t max_depth);

/**
 * Converts a row-compressed representation to a row-diff representation based on the given graph
 * @param graph the graph used for selecting rows to diff
 * @param annotation input annotation
 * @param num_threads number of threads to use when traversing the graph
 * @param max_depth maximum allowed path length before forcing a terminal node
 * @return newly constructed row-diff annotation
 */
std::unique_ptr<RowDiffAnnotator> convert_to_row_diff(const graph::DBGSuccinct &graph,
                                                      RowCompressed<std::string> &&annotation,
                                                      uint32_t num_threads = 1,
                                                      uint32_t max_depth = -1);


} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_CONVERTERS_HPP__
