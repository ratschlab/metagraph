#pragma once

#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {
namespace annot {

void count_labels_per_row(const std::vector<std::string> &source_files,
                          const std::string &row_count_fname);

void build_pred_succ(const std::string &graph_filename,
                     const std::string &outfbase,
                     const std::string &count_vectors_dir,
                     const std::string &row_count_extension,
                     uint32_t num_threads);

void assign_anchors(const std::string &graph_filename,
                    const std::string &outfbase,
                    const std::filesystem::path &dest_dir,
                    uint32_t max_length,
                    const std::string &row_reduction_extension,
                    uint32_t num_threads);

void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &col_out_dir,
                               const std::filesystem::path &swap_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size_bytes,
                               bool compute_row_reduction = true);

} // namespace annot
} // namespace mtg
