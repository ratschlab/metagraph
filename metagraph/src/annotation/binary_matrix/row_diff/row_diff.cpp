#include "row_diff.hpp"

#include <filesystem>
#include <fstream>

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


/**
 * Returns the given column.
 */
template <class BaseMatrix>
std::vector<BinaryMatrix::Row> RowDiff<BaseMatrix>::get_column(Column column) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    const graph::boss::BOSS &boss = graph_->get_boss();
    assert(!fork_succ_.size() || fork_succ_.size() == boss.get_last().size());

    std::vector<Row> result;

    // // TODO: implement a more efficient algorithm
    // for (Row row = 0; row < num_rows(); ++row) {
    //     auto edge = graph_->kmer_to_boss_index(
    //         graph::AnnotatedSequenceGraph::anno_to_graph_index(row)
    //     );

    //     if (boss.get_W(edge) && get(row, column))
    //         result.push_back(row);
    // }
    // return result;

    diffs_.call_columns({ column }, [&](size_t, const bitmap &starts) {
        starts.call_ones([&](uint64_t row) {
            std::vector<uint64_t> rows;
            uint64_t boss_edge = graph_->kmer_to_boss_index(
                    graph::AnnotatedSequenceGraph::anno_to_graph_index(row));
            const bit_vector &rd_succ = fork_succ_.size() ? fork_succ_ : boss.get_last();

            while (!anchor_[row]) {
                rows.emplace_back(row);
                boss_edge = boss.row_diff_successor(boss_edge, rd_succ);
                row = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                        graph_->boss_to_kmer_index(boss_edge));
            }

            rows.emplace_back(row);
            std::vector<bool> vals(rows.size(), diffs_.get(row, column));
            auto it = rows.rbegin() + 1;
            auto jt = vals.rbegin() + 1;
            for ( ; it != rows.rend(); ++it, ++jt) {
                *jt = *(jt - 1) ^ diffs_.get(*it, column);
            }

            for (size_t i = 0; i < rows.size(); ++i) {
                if (vals[i])
                    result.emplace_back(rows[i]);
            }

            auto populate = [&](uint64_t start_row) {
                std::vector<uint64_t> stack {
                    graph::AnnotatedSequenceGraph::anno_to_graph_index(start_row)
                };
                while (stack.size()) {
                    uint64_t node = stack.back();
                    stack.pop_back();

                    graph_->call_incoming_kmers(node, [&](auto prev_node, char c) {
                        uint64_t row = graph::AnnotatedSequenceGraph::graph_to_anno_index(prev_node);
                        if (c != boss.kSentinel && row) {
                            if (!anchor_[row]
                                    && rd_succ[graph_->kmer_to_boss_index(prev_node)]
                                    && !diffs_.get(row, column)) {
                                result.emplace_back(row);
                                stack.emplace_back(prev_node);
                            }
                        }
                    });
                }
            };

            if (vals.size() > 1) {
                auto kt = rows.rbegin() + 1;
                for (auto it = vals.rbegin() + 1; it != vals.rend(); ++it, ++kt) {
                    if (*it && !*(it - 1)) {
                        auto jt = it;
                        auto lt = kt;
                        while (jt != vals.rend() && *jt
                                && graph_->has_single_incoming(graph::AnnotatedSequenceGraph::anno_to_graph_index(*lt))) {
                            ++jt;
                            ++lt;
                        }

                        if (jt == vals.rend() || !*jt) {
                            populate(*(lt - 1));
                            it = jt - 1;
                            kt = lt - 1;
                        } else {
                            populate(*lt);
                            it = jt;
                            kt = lt;
                        }
                    }
                }
            } else if (rows[0]) {
                populate(rows[0]);
            }
        });
    });

    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
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
