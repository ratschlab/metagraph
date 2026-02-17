#include "row_diff.hpp"

#include <filesystem>
#include <fstream>

#include <tsl/hopscotch_set.h>
#include <ips4o.hpp>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "common/utils/file_utils.hpp"

namespace mtg {
namespace annot {

using node_index = graph::DeBruijnGraph::node_index;

node_index row_diff_successor(const graph::DeBruijnGraph &graph,
                              node_index node,
                              const bit_vector &rd_succ) {
    if (auto* dbg_succ = dynamic_cast<graph::DBGSuccinct const*>(&graph)) {
        return dbg_succ->get_boss().row_diff_successor(
            node,
            rd_succ.size() ? rd_succ : dbg_succ->get_boss().get_last()
        );
    } else {
        assert(rd_succ.size());
        node_index succ = graph::DeBruijnGraph::npos;
        graph.adjacent_outgoing_nodes(node, [&](node_index adjacent_node) {
            if (rd_succ[adjacent_node]) {
                succ = adjacent_node;
            }
        });
        assert(graph.in_graph(succ) && "a row diff successor must exist");
        return succ;
    }
}


namespace matrix {

void IRowDiff::load_anchor(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
        common::logger->error("Can't read anchor file: {}", filename);
        std::exit(1);
    }
    std::unique_ptr<std::ifstream> f = utils::open_ifstream(filename);
    if (!f->good()) {
        common::logger->error("Could not open anchor file {}", filename);
        std::exit(1);
    }
    anchor_.load(*f);
}

void IRowDiff::load_fork_succ(const std::string &filename) {
    if (!std::filesystem::exists(filename)) {
        common::logger->error("Can't read fork successor file: {}", filename);
        std::exit(1);
    }
    std::unique_ptr<std::ifstream> f = utils::open_ifstream(filename);
    if (!f->good()) {
        common::logger->error("Could not open fork successor file {}", filename);
        std::exit(1);
    }
    fork_succ_.load(*f);
}

std::tuple<std::vector<BinaryMatrix::Row>, std::vector<std::vector<size_t>>, std::vector<size_t>>
IRowDiff::get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids, size_t num_threads) const {
    assert(graph_ && "graph must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    using Row = BinaryMatrix::Row;

    // map row index to its index in |rd_rows|
    VectorMap<Row, size_t> node_to_rd;
    node_to_rd.reserve(row_ids.size() * RD_PATH_RESERVE_SIZE);

    // Truncated row-diff paths, indexes to |rd_rows|.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());

    const size_t kMaxBlockSize = 1000;
    const size_t block_size = num_threads > 1
            ? std::min<size_t>(kMaxBlockSize, (row_ids.size() + num_threads - 1) / num_threads)
            : row_ids.size();
    tsl::hopscotch_set<Row> visited_rows;
    #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic) private(visited_rows)
    for (size_t begin = 0; begin < row_ids.size(); begin += block_size) {
        size_t visited_before = visited_rows.size();
        // Phase 1: Parallel path tracing.
        // Thread-private node_to_rd provides intra-thread dedup (early path
        // termination when a node was already visited by the same thread).
        // Cross-thread dedup is handled in Phase 2.
        const size_t end = std::min<size_t>(begin + block_size, row_ids.size());
        for (size_t i = begin; i < end; ++i) {
            Row row = row_ids[i];

            node_index node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);

            while (true) {
                assert(graph_->in_graph(node));
                row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);

                bool is_new = visited_rows.insert(row).second;
                rd_paths_trunc[i].push_back(row);

                // If a node had been reached before, we interrupt the diff path.
                // The annotation for that node will have been reconstructed earlier
                // than for other nodes in this path as well. Thus, we will start
                // reconstruction from that node and don't need its successors.
                if (!is_new)
                    break;

                if (anchor_[row])
                    break;

                node = row_diff_successor(*graph_, node, fork_succ_);
            }
        }

        // Phase 2: Serial index assignment with cross-thread deduplication.
        // Re-processes stored row IDs through a global map to assign canonical
        // indices and truncate paths that converge across thread boundaries.
        #pragma omp ordered
        {
            // insert new rows visited by other threads during this thread's execution
            for (auto it = node_to_rd.begin() + visited_before; it != node_to_rd.end(); ++it) {
                visited_rows.insert(it->first);
            }
            for (size_t i = begin; i < end; ++i) {
                for (size_t j = 0; j < rd_paths_trunc[i].size(); ++j) {
                    auto [it, is_new] = node_to_rd.try_emplace(rd_paths_trunc[i][j], node_to_rd.size());
                    rd_paths_trunc[i][j] = it.value();
                    if (!is_new) {
                        rd_paths_trunc[i].resize(j + 1);
                        break;
                    }
                }
            }
            assert(visited_rows.size() == node_to_rd.size());  // this thread's cache is synced now
        }
    }

    auto m = to_vector(std::move(node_to_rd));
    // sort by indexes of rd rows
    // the second value points to the index in batch
    ips4o::parallel::sort(m.begin(), m.end(), utils::LessFirst(), num_threads);
    // collect an array of rd rows
    std::vector<Row> rd_ids(m.size());
    for (size_t i = 0; i < m.size(); ++i) {
        rd_ids[i] = m[i].first;
    }
    // make m[].first map indexes in batch to new indexes (where rows are sorted)
    for (size_t i = 0; i < m.size(); ++i) {
        m[m[i].second].first = i;
    }

    // keeps how many times rows in |rd_rows| will be queried
    std::vector<size_t> times_traversed(rd_ids.size(), 0);

    for (size_t i = 0; i < row_ids.size(); ++i) {
        for (auto &j : rd_paths_trunc[i]) {
            j = m[j].first;
            times_traversed[j]++;
        }
    }

    return { rd_ids, rd_paths_trunc, times_traversed };
}

} // namespace matrix
} // namespace annot
} // namespace mtg
