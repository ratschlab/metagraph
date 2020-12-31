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

void build_successor(const std::string &graph_filename,
                     const std::string &outfbase,
                     uint32_t max_length,
                     uint32_t num_threads);

void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::string &anchors_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size,
                               bool compute_row_reduction = true);

void optimize_anchors_in_row_diff(const std::string &graph_fname,
                                  const std::filesystem::path &dest_dir,
                                  const std::string &row_reduction_extension);

} // namespace annot
} // namespace mtg
