#include "row_diff.hpp"

#include <filesystem>
#include <fstream>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"

namespace mtg {
namespace annot {
namespace binmat {

void IRowDiff::load_anchor(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
        common::logger->error("Can't read anchor file: {}", filename);
        std::exit(1);
    }
    std::ifstream f(filename, ios::binary);
    if (!f.good()) {
        common::logger->error("Could not open anchor file {}", filename);
        std::exit(1);
    }
    anchor_.load(f);
}

void IRowDiff::load_fork_succ(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
        common::logger->error("Can't read fork successor file: {}", filename);
        std::exit(1);
    }
    std::ifstream f(filename, ios::binary);
    if (!f.good()) {
        common::logger->error("Could not open fork successor file {}", filename);
        std::exit(1);
    }
    fork_succ_.load(f);
}

std::pair<std::vector<BinaryMatrix::Row>, std::vector<std::vector<size_t>>>
IRowDiff::get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids) const {
    assert(graph_ && "graph must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());

    using Row = BinaryMatrix::Row;

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    rd_ids.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());

    const graph::boss::BOSS &boss = graph_->get_boss();
    const bit_vector &rd_succ = fork_succ_.size() ? fork_succ_ : boss.get_last();

    for (size_t i = 0; i < row_ids.size(); ++i) {
        Row row = row_ids[i];

        graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(
                graph::AnnotatedSequenceGraph::anno_to_graph_index(row));

        while (true) {
            row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                    graph_->boss_to_kmer_index(boss_edge));

            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            rd_paths_trunc[i].push_back(it.value());

            // If a node had been reached before, we interrupt the diff path.
            // The annotation for that node will have been reconstructed earlier
            // than for other nodes in this path as well. Thus, we will start
            // reconstruction from that node and don't need its successors.
            if (!is_new)
                break;

            rd_ids.push_back(row);

            if (anchor_[row])
                break;

            boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
        }
    }

    return std::make_pair(std::move(rd_ids), std::move(rd_paths_trunc));
}

} // namespace binmat
} // namespace annot
} // namespace mtg
