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

void convert_batch_to_row_diff(const std::string &graph_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir,
                               uint32_t max_depth);

} // namespace annot
} // namespace mtg
