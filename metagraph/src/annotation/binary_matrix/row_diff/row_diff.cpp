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

std::pair<std::vector<BinaryMatrix::Row>, std::vector<std::vector<std::pair<size_t, size_t>>>>
IRowDiff::get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids) const {
    assert(graph_ && "graph must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->num_nodes() + 1);

    using Row = BinaryMatrix::Row;

    const size_t RD_PATH_RESERVE_SIZE = 2;

    // diff rows annotating nodes along the row-diff paths
    std::vector<Row> rd_ids;
    rd_ids.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<std::pair<size_t, size_t>>> rd_paths_trunc(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        std::vector<std::pair<size_t, size_t>> &rd_path = rd_paths_trunc[i];

        std::vector<size_t> path;
        Vector<std::pair<size_t, Row>> queue;
        queue.emplace_back(0, row_ids[i]);

        while (queue.size()) {
            size_t depth = queue.back().first;
            Row row = queue.back().second;
            queue.pop_back();
            while (depth < path.size()) {
                assert(path.size() > 1);
                rd_path.emplace_back(*(path.rbegin() + 1), *path.rbegin());
                path.pop_back();
            }
            auto [it, is_new] = node_to_rd.try_emplace(row, rd_ids.size());
            path.push_back(it.value());
            // If a node had been reached before, we interrupt the diff path.
            // The annotation for that node will have been reconstructed earlier
            // than for other nodes in this path as well. Thus, we will start
            // reconstruction from that node and don't need its successors.
            if (!is_new)
                continue;

            rd_ids.push_back(row);

            if (anchor_[row])
                continue;

            auto node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);
            graph_->call_row_diff_successors(node, fork_succ_, [&](auto succ) {
                queue.emplace_back(depth + 1, graph::AnnotatedSequenceGraph::graph_to_anno_index(succ));
            });
        }

        while (path.size() > 1) {
            rd_path.emplace_back(*(path.rbegin() + 1), *path.rbegin());
            path.pop_back();
        }
        assert(path.size());
        rd_path.emplace_back(-1, path[0]);
    }

    return std::make_pair(std::move(rd_ids), std::move(rd_paths_trunc));
}

} // namespace binmat
} // namespace annot
} // namespace mtg
