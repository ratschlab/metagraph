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

std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_simple_BRWT(RowDiffAnnotator &&annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);

std::unique_ptr<MultiBRWTAnnotator>
convert_to_greedy_BRWT(ColumnCompressed<std::string>&& annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);

std::unique_ptr<RowDiffBRWTAnnotator>
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

void relax_BRWT(binmat::BRWT *annotation, size_t relax_max_arity, size_t num_threads = 1);

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
 * @param files annotations in #ColumnCompressed format. Typically, each source
 * contains a single column, but this is not a requirement
 * @param graph the graph used for selecting adjacent nodes
 * @param graph_path directory where the pred/succ/term vectors that determine the diff
 * succession will be dumped
 * @param mem_bytes available memory for loading columns; the function will load as many
 * columns as can fit in memory in one go and transform the entire batch, then repeat with
 * the next batch until all input files are exhausted
 * @param max_path_length maximum distance between terminal nodes, i.e. nodes whose annotation
 * is fully stored
 * @param out_dir directory where the transformed columns will be dumped. Filenames are
 * kept, extension is changed from 'column.annodbg' to 'row_diff.annodbg'
 */
void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string& graph_fname,
                         size_t mem_bytes,
                         uint32_t max_path_length,
                         const std::filesystem::path &dest_dir);

void convert_row_diff_to_col_compressed(const std::vector<std::string> &files,
                                        const std::string &outfbase);

/**
 * Converts a RowDiff annotation into RowSparse.
 */
std::unique_ptr<RowSparseAnnotator> convert(const RowDiffAnnotator &annotator);

template <class BinaryMatrix>
void wrap_in_row_diff(const std::string &anno_file,
                      const std::string &graph_file,
                      const std::string &out_dir);


} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_CONVERTERS_HPP__
