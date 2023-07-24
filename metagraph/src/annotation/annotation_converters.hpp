#ifndef __ANNOTATION_CONVERTERS_HPP__
#define __ANNOTATION_CONVERTERS_HPP__

#include <memory>
#include <vector>
#include <filesystem>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

namespace mtg {
namespace annot {

template <typename LabelType>
class MultiLabelAnnotation;


template <typename Label>
class ColumnCompressed;

template <class StaticAnnotation, typename Label>
std::unique_ptr<StaticAnnotation> convert(ColumnCompressed<Label>&& annotation);

template <class StaticAnnotation, typename Label>
void convert(ColumnCompressed<Label>&& annotation, const std::string &outfbase) {
    convert<StaticAnnotation, Label>(std::move(annotation))->serialize(outfbase);
}

template <>
void convert<RowFlatAnnotator, std::string>(ColumnCompressed<std::string>&& annotator,
                                            const std::string &outfbase);

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

void relax_BRWT(matrix::BRWT *annotation, size_t relax_max_arity, size_t num_threads = 1);

template <class StaticAnnotation>
std::unique_ptr<StaticAnnotation>
convert_to_RbBRWT(const std::vector<std::string> &annotation_files,
                  size_t max_brwt_arity);

// For RowDiffDisk/RowDiffDiskCoord/RowDiffRowFlat/RowDiffRowSparse -Annotator
template <class RowDiffAnnotator>
void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string &anchors_file_fbase,
                         const std::string &outfbase,
                         size_t num_threads,
                         size_t mem_bytes);

void merge_row_compressed(const std::vector<std::string> &filenames,
                          const std::string &outfile);

void merge_brwt(const std::vector<std::string> &filenames,
                const std::string &outfile);

// transform to RowCompressed<Label>
template <typename Label>
void convert_to_row_annotator(const ColumnCompressed<Label> &annotator,
                              const std::string &outfbase);

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
 * @param with_values row-diff transform with k-mer counts/attributes
 * @param with_coordinates row-diff transform with k-mer coordinates/attributes
 * @param num_coords_per_seq assume a constant length of each sequence (0: off)
 */
enum class RowDiffStage { COUNT_LABELS = 0, COMPUTE_REDUCTION, CONVERT };
void convert_to_row_diff(const std::vector<std::string> &files,
                         const std::string &graph_fname,
                         size_t mem_bytes,
                         uint32_t max_path_length,
                         std::filesystem::path out_dir,
                         std::filesystem::path swap_dir,
                         RowDiffStage construction_stage,
                         std::filesystem::path count_vector_fname = "",
                         bool with_values = false,
                         bool with_coordinates = false,
                         size_t num_coords_per_seq = 0);

void convert_row_diff_to_col_compressed(const std::vector<std::string> &files,
                                        const std::string &outfbase);

std::pair<std::string, std::string> get_anchors_and_fork_fnames(const std::string &fbase);

template <class Annotator>
StaticBinRelAnnotator<matrix::TupleCSCMatrix<typename Annotator::binary_matrix_type>, std::string>
load_coords(Annotator&& anno, const std::vector<std::string> &files);

} // namespace annot
} // namespace mtg

#endif // __ANNOTATION_CONVERTERS_HPP__
