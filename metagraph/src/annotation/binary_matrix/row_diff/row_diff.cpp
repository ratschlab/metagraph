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
void RowDiff<BaseMatrix>::serialize(const std::string &filename) const {
    std::ofstream f(filename, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(const std::string &filename) {
    std::ifstream f(filename, ios::binary);
    bool result = load(f);
    return result;
}

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

template <class BaseMatrix>
Vector<size_t>
RowDiff<BaseMatrix>::extend_maximal(const std::vector<Row> &rows,
                                    const SetBitPositions &columns) const {
    return BinaryMatrix::extend_maximal(rows, columns);

    if (rows.empty() || columns.empty())
        return {};

    size_t num_anchors = 0;
    sdsl::bit_vector is_anchor(rows.size(), false);
    for (size_t i = 0; i < rows.size(); ++i) {
        is_anchor[i] = anchor_[rows[i]];
        num_anchors += is_anchor[i];
    }

    const graph::boss::BOSS &boss = graph_->get_boss();

    std::vector<graph::boss::BOSS::edge_index> edges(rows.size(), 0);
    std::vector<graph::boss::BOSS::TAlphabet> edge_chars(rows.size(), 0);
    for (size_t i = 0; i < rows.size(); ++i) {
        edges[i] = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(rows[i]));
        edge_chars[i] = boss.get_W(edges[i]) % boss.alph_size;
    }

    bool forward = true;
    for (size_t i = 1; i < edges.size() && forward; ++i) {
        forward = (boss.fwd(edges[i - 1], edge_chars[i - 1]) == edges[i]);
    }

    if (forward) {
        if (num_anchors && !is_anchor[rows.size() - 1])
            return BinaryMatrix::extend_maximal(rows, columns);

        Vector<size_t> result(columns.size());
        for (size_t i = 0; i < columns.size(); ++i) {
            sdsl::bit_vector mask = diffs_.has_column(rows, columns[i]);
            if (!num_anchors) {
                result[i] = std::min(rows.size(), sdsl::util::next_bit(mask, 0) + 1);
            } else {
                bool val = mask[mask.size() - 1];
                size_t res = rows.size();
                for (size_t j = rows.size() - 1; j > 0; --j) {
                    if (!is_anchor[j - 1]) {
                        val ^= mask[j - 1];
                    } else {
                        val = mask[j - 1];
                    }

                    if (!val)
                        res = std::min(res, j);
                }
                assert(val);
                result[i] = res;
            }
        }

        return result;
    }

    bool reverse = true;
    for (size_t i = edges.size() - 1; i > 0 && reverse; --i) {
        reverse = (boss.fwd(edges[i], edge_chars[i]) == edges[i - 1]);
    }

    if (reverse) {
        Vector<size_t> result(columns.size());
        for (size_t i = 0; i < columns.size(); ++i) {
            sdsl::bit_vector mask = diffs_.has_column(rows, columns[i]);
            if (!num_anchors) {
                result[i] = sdsl::util::next_bit(mask, 1);
            } else {
                bool val = true;
                for (size_t j = 1; j < rows.size(); ++j) {
                    if (!is_anchor[j]) {
                        val ^= mask[j];
                    } else {
                        val = mask[j];
                    }

                    if (!val) {
                        result[i] = j;
                        break;
                    }
                }
            }
        }

        return result;
    }

    return BinaryMatrix::extend_maximal(rows, columns);
}

template <typename Key, class Hash = std::hash<Key>, class EqualTo = std::equal_to<Key>,
          class Allocator = std::allocator<Key>, class Container = std::vector<Key, Allocator>,
          typename Size = uint64_t>
using VectorSet = tsl::ordered_set<Key, Hash, EqualTo, Allocator, Container, Size>;

template <class BaseMatrix>
sdsl::bit_vector
RowDiff<BaseMatrix>::has_column(const std::vector<Row> &row_ids, Column column) const {
    sdsl::bit_vector result(row_ids.size(), false);

    if (row_ids.empty())
        return result;

    const graph::boss::BOSS &boss = graph_->get_boss();

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

            graph::boss::BOSS::TAlphabet w = boss.get_W(boss_edge);
            assert(boss_edge > 1 && w != 0);
            // fwd always selects the last outgoing edge for a given node
            boss_edge = boss.fwd(boss_edge, w % boss.alph_size);
        }
    }

    sdsl::bit_vector row_set_mask = diffs_.has_column(row_set.values_container(), column);

    for (size_t i = 0; i < row_ids.size(); ++i) {
        bool found = row_set_mask[rd_paths_trunc[i].back()];
        for (auto it = rd_paths_trunc[i].rbegin() + 1; it != rd_paths_trunc[i].rend(); ++it) {
            row_set_mask[*it] = (found ^= row_set_mask[*it]);
        }

        result[i] = found;
    }

    return result;
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
