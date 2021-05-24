#include "row_diff.hpp"

#include <filesystem>
#include <fstream>

#include <tsl/ordered_set.h>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
void RowDiff<BaseMatrix>::load_anchor(const std::string &filename) {
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
    f.close();
}

template <typename Key, class Hash = std::hash<Key>, class EqualTo = std::equal_to<Key>,
          class Allocator = std::allocator<Key>, class Container = std::vector<Key, Allocator>,
          typename Size = uint64_t>
using VectorSet = tsl::ordered_set<Key, Hash, EqualTo, Allocator, Container, Size>;

template <class BaseMatrix>
std::vector<sdsl::bit_vector>
RowDiff<BaseMatrix>::has_column(const std::vector<Row> &row_ids,
                                const SetBitPositions &columns) const {
    std::vector<sdsl::bit_vector> results;
    results.reserve(columns.size());
    for (size_t i = 0; i < columns.size(); ++i) {
        results.emplace_back(row_ids.size(), false);
    }

    if (row_ids.empty())
        return results;

    const graph::boss::BOSS &boss = graph_->get_boss();
    const bit_vector &rd_succ = fork_succ_.size() ? fork_succ_ : boss.get_last();

    VectorSet<Row> row_set;
    std::vector<std::vector<size_t>> rd_paths_trunc(row_ids.size());
    for (size_t i = 0; i < row_ids.size(); ++i) {
        Row row = row_ids[i];

        graph::boss::BOSS::edge_index boss_edge = graph_->kmer_to_boss_index(
                graph::AnnotatedSequenceGraph::anno_to_graph_index(row));

        while (true) {
            row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                    graph_->boss_to_kmer_index(boss_edge));

            auto [it, is_new] = row_set.emplace(row);
            rd_paths_trunc[i].push_back(it - row_set.begin());
            assert(row_set.values_container()[rd_paths_trunc[i].back()] == row);

            // If a node had been reached before, we interrupt the diff path.
            // The annotation for that node will have been reconstructed earlier
            // than for other nodes in this path as well. Thus, we will start
            // reconstruction from that node and don't need its successors.
            if (!is_new || anchor()[row])
                break;

            boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
        }
    }

    std::vector<sdsl::bit_vector> row_set_masks = diffs_.has_column(
        row_set.values_container(), columns
    );

    for (size_t j = 0; j < columns.size(); ++j) {
        sdsl::bit_vector &result = results[j];
        sdsl::bit_vector &row_set_mask = row_set_masks[j];
        for (size_t i = 0; i < row_ids.size(); ++i) {
            bool found = row_set_mask[rd_paths_trunc[i].back()];
            for (auto it = rd_paths_trunc[i].rbegin() + 1;
                    it != rd_paths_trunc[i].rend(); ++it) {
                row_set_mask[*it] = (found ^= row_set_mask[*it]);
            }

            result[i] = found;
        }
    }

    return results;
}

template <class BaseMatrix>
void RowDiff<BaseMatrix>::load_fork_succ(const std::string &filename) {
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
    f.close();
}

template
class RowDiff<ColumnMajor>;
template
class RowDiff<BRWT>;
template
class RowDiff<RowSparse>;

} // namespace binmat
} // namespace annot
} // namespace mtg
