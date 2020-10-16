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

using node_index = graph::boss::BOSS::node_index;
using CallOnes = std::function<void(const bit_vector &source_col,
                                    node_index row_idx,
                                    node_index start_idx,
                                    uint64_t source_idx,
                                    uint64_t col_idx,
                                    const std::vector<uint64_t> &succ_chunk,
                                    const std::vector<uint64_t> &pred_chunk,
                                    const std::vector<uint64_t> &pred_chunk_idx)>;

void traverse_anno_chunked(
        const std::string &name,
        const graph::DBGSuccinct &graph,
        const std::string &outfbase,
        const std::vector<std::unique_ptr<annot::ColumnCompressed<>>> &col_annotations,
        const std::function<void()> &before_chunk,
        const CallOnes &call_ones,
        const std::function<void(node_index start, uint64_t size)> &after_chunk);

void convert_batch_to_row_diff(const graph::DBGSuccinct &graph,
                               const std::string &graph_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir,
                               uint32_t max_depth);

}
}
