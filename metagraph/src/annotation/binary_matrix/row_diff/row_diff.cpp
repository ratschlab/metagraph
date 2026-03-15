#include "row_diff.hpp"

#include <filesystem>
#include <fstream>

#include <ips4o.hpp>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vector_set.hpp"

namespace mtg {
namespace annot {

using node_index = graph::DeBruijnGraph::node_index;
using common::logger;

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

std::tuple<std::vector<BinaryMatrix::Row>,
           std::vector<std::vector<size_t>>,
           std::vector<size_t>,
           std::vector<std::vector<size_t>>>
IRowDiff::get_rd_ids(const std::vector<BinaryMatrix::Row> &row_ids, size_t num_threads) const {
    assert(graph_ && "graph must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);
    num_threads = std::max<size_t>(1, num_threads);

    using Row = BinaryMatrix::Row;

    // Thread-local sets indexing rows.
    std::vector<VectorSet<Row>> rows_visited_local(num_threads);

    // Truncated row-diff paths, store thread-local indexes of deduped rows.
    // The last index in each path points to an anchor or to a row which had
    // been reached before, and thus, will be reconstructed before this one.
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());
    // Records which paths were traced by each thread (for offset remapping in Phase 2).
    std::vector<std::vector<size_t>> groups(num_threads);

    const size_t block_size = get_chunk_size(row_ids.size(), 10000 /* max_chunk_size */, num_threads);

    // Phase 1: Parallel path tracing.
    // |rows_visited_local| provides intra-thread dedup (early path termination
    // when a node was already visited by the same thread).
    // FYI: It's essential that the omp loop is ORDERED so that each thread gets rows in order.
    #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic, block_size)
    for (size_t i = 0; i < row_ids.size(); ++i) {
        const int t = omp_get_thread_num();
        auto &rows_visited = rows_visited_local[t];

        groups[t].push_back(i);

        Row row = row_ids[i];

        node_index node = graph::AnnotatedSequenceGraph::anno_to_graph_index(row);

        while (true) {
            assert(graph_->in_graph(node));
            row = graph::AnnotatedSequenceGraph::graph_to_anno_index(node);

            auto [it, is_new] = rows_visited.emplace(row);
            rd_paths_trunc[i].push_back(it - rows_visited.begin());

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

    // Phase 2: Merge thread-local maps and remap path indices.
    //
    // Note: We skip the cross-thread dedup because it it has to be sequential
    //       and that adds too much overhead. Besides, the intra-thread dedup
    //       already significantly reduces redundant paths.
    std::vector<size_t> offsets(groups.size() + 1);
    for (size_t g = 0; g < groups.size(); ++g) {
        offsets[g + 1] = offsets[g] + rows_visited_local[g].size();
    }
    std::vector<std::pair<Row, size_t>> m(offsets.back());
    std::vector<Row> rd_ids(offsets.back());
    std::vector<size_t> new_idx(offsets.back());
    std::vector<size_t> times_traversed(offsets.back());
    #pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (size_t g = 0; g < groups.size(); ++g) {
        const size_t offset = offsets[g];
        const auto &vec = rows_visited_local[g].values_container();
        for (size_t i = 0; i < vec.size(); ++i) {
            m[offset + i] = { vec[i], offset + i };
        }
    }

    // sort by indexes of rd rows
    // the second value points to the index in batch
    // (row_4, 4), (row_4, 0), (row_7, 2), ...
    ips4o::parallel::sort(m.begin(), m.end(), utils::LessFirst(), num_threads);

    // map indexes in batch (used in `rd_paths_trunc`) to new indexes
    // (where rows are sorted).
    #pragma omp parallel for num_threads(num_threads) schedule(static)
    for (size_t i = 0; i < m.size(); ++i) {
        rd_ids[i] = m[i].first;
        new_idx[m[i].second] = i;
    }

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t g = 0; g < groups.size(); ++g) {
        for (size_t i : groups[g]) {
            for (auto &j : rd_paths_trunc[i]) {
                j = new_idx[j + offsets[g]];
                times_traversed[j]++;
            }
        }
    }

    return { std::move(rd_ids), std::move(rd_paths_trunc),
             std::move(times_traversed), std::move(groups) };
}

} // namespace matrix
} // namespace annot
} // namespace mtg
