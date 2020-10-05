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

std::unique_ptr<MultiBRWTAnnotator>
convert_to_simple_BRWT(ColumnCompressed<std::string>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);

std::unique_ptr<BRWTRowDiffAnnotator>
convert_to_simple_BRWT(RowDiffAnnotator &&annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);

std::unique_ptr<MultiBRWTAnnotator>
convert_to_greedy_BRWT(ColumnCompressed<std::string>&& annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);

std::unique_ptr<BRWTRowDiffAnnotator>
convert_to_greedy_BRWT(RowDiffAnnotator &&annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation>
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
std::unique_ptr<StaticAnnotation>
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

/**
 * Sparsifies annotations in #ColumnCompressed format by storing diffs between sucessive
 * nodes rather than the actual annotation.
 * Since adjacent nodes in a graph have a high probability of sharing annotations, storing
 * diffs rather than the full annotation results in a sparser matrix.
 * @param graph the graph used for selecting adjacent nodes
 * @param sources list of annotations in #ColumnCompressed format. Typically, each source
 * contains a single column, but this is not a requirement
 * @param outfbase directory where the pred/succ/term vectors that determine the diff
 * succession will be dumped
 * @param max_depth maximum distance between terminal nodes, i.e. nodes whose annotation
 * is fully stored
 * @return the sparsified columns, wrapped in a static annotator that contains a #RowDiff
 * that wraps the column in #ColumnMajor format
 */
template <typename Label>
std::vector<std::unique_ptr<RowDiffAnnotator>>
convert_to_row_diff(const graph::DBGSuccinct &graph,
                       const std::vector<std::unique_ptr<ColumnCompressed<Label>>> &sources,
                       const std::string &outfbase,
                       uint32_t max_depth);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_CONVERTERS_HPP__
