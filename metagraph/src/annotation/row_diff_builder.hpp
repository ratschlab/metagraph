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

/** Marker type to indicate a value represent a node index in a BOSS graph */
using node_index = graph::boss::BOSS::node_index;

/**
 * Callback invoked by #traverse_anno_chunked for each set bit in the annotation matrix.
 * @param source_col the column for which the callback was invoked, in bit_vector format
 * @param row_idx the row in which the bit is set
 * @param row_idx_chunk relative index of the row in the current chunk
 * @param source_idx index of the source file for the current column
 * @param coll_idx index of the column in the current source file (typically, the source
 *        files contain a single column each, but that's not a requirement)
 * @param succ_chunk current chunk of successor values (indexed by #row_idx_chunk)
 * @param pred_chunk current chunk of predecessor values (indexed by #row_idx_chunk)
 * @param pred_chunk_idx indexes pred_chunk. The predecessors of #row_idx are located
 *        in #pred_chunk between pred_chunk_idx[row_idx_chunk] and
 *        pred_chunk_idx[row_idx_chunk + 1]
 */
using CallOnes = std::function<void(const bit_vector &source_col,
                                    node_index row_idx,
                                    node_index row_idx_chunk,
                                    uint64_t source_idx,
                                    uint64_t col_idx,
                                    const std::vector<uint64_t> &succ_chunk,
                                    const std::vector<uint64_t> &pred_chunk,
                                    const std::vector<uint64_t> &pred_chunk_idx)>;

/**
 * Traverses a group of column compressed annotations (loaded in memory) in chunks of
 * 1'000'000 rows at a time and invokes #call_ones for each set bit.
 * @param log_header label to be displayed in the progress bar
 * @param num_rows number of rows in the annotation
 * @param pred_succ_fprefix prefix for the pred/succ files containg the predecessors and
 * the successor for each node
 * @param col_annotations the annotations to be traversed
 * @param before_chunk callback to invoke before a chunk is traversed
 * @param call_ones callback to invoke on a set bit
 * @param after_chunk callback to invoke after a chunk is traversed
 */
void traverse_anno_chunked(
        const std::string &log_header,
        uint64_t num_rows,
        const std::string &pred_succ_fprefix,
        const std::vector<std::unique_ptr<annot::ColumnCompressed<>>> &col_annotations,
        const std::function<void()> &before_chunk,
        const CallOnes &call_ones,
        const std::function<void(node_index start, uint64_t size)> &after_chunk);

void convert_batch_to_row_diff(const std::string &graph_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir);

} // namespace annot
} // namespace mtg
