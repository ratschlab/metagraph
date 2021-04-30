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

// TODO: remove
std::unique_ptr<MultiBRWTAnnotator>
convert_to_simple_BRWT(ColumnCompressed<std::string>&& annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);
// TODO: remove
std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_simple_BRWT(RowDiffColumnAnnotator &&annotation,
                       size_t grouping_arity = 2,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1);
// TODO: remove
std::unique_ptr<MultiBRWTAnnotator>
convert_to_greedy_BRWT(ColumnCompressed<std::string>&& annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);
// TODO: remove
std::unique_ptr<RowDiffBRWTAnnotator>
convert_to_greedy_BRWT(RowDiffColumnAnnotator &&annotation,
                       size_t num_parallel_nodes = 1,
                       size_t num_threads = 1,
                       uint64_t num_rows_subsampled = 1'000'000);

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation>
convert_to_BRWT(const std::vector<std::string> &annotation_files,
                const std::vector<std::vector<uint64_t>> &linkage_matrix,
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
 * @param swap_dir directory for temporary files
 */
enum class RowDiffStage { COUNT_LABELS = 0, COMPUTE_REDUCTION, CONVERT };
void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string &graph_fname,
                         size_t mem_bytes,
                         uint32_t max_path_length,
                         std::filesystem::path out_dir,
                         std::filesystem::path swap_dir,
                         RowDiffStage construction_stage,
                         std::filesystem::path count_vector_fname = "");

void convert_row_diff_to_col_compressed(const std::vector<std::string> &files,
                                        const std::string &outfbase);

/**
 * Converts a RowDiff annotation into RowDiff<RowSparse>.
 */
std::unique_ptr<RowDiffRowSparseAnnotator>
convert_row_diff_to_RowDiffSparse(const std::vector<std::string> &filenames);

/**
 * Wraps an existing annotation (e.g. BRWT) into a RowDiff annotation. Typically this
 * happens when transforming RowDiff columns back to column compress, manipulate the
 * column compressed into some other format, and then wrapping the result back into a
 * RowDiff.
 */
void wrap_in_row_diff(MultiLabelEncoded<std::string> &&anno,
                      const std::string &graph_file,
                      const std::string &out_file);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_CONVERTERS_HPP__
