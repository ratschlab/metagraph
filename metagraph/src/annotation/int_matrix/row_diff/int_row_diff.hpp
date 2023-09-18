#ifndef __INT_ROW_DIFF_HPP__
#define __INT_ROW_DIFF_HPP__

#include <algorithm>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vector_map.hpp"
#include "common/vector.hpp"
#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"


namespace mtg {
namespace annot {
namespace matrix {

/**
 * Convert deltas to positive integer for enable compression:
 *      0  -> X (not allowed, zero diffs must be skipped)
 *      1  -> 0
 *      -1 -> 1
 *      2  -> 2
 *      -2 -> 3
 *      ...
 */
inline uint64_t encode_diff(int64_t x) {
    assert(x);
    return (std::abs(x) - 1) * 2 + (x < 0);
}

inline int64_t decode_diff(uint64_t c) {
    return !(c & 1) ? c / 2 + 1 : -((c + 1) / 2);
}

template <class BaseMatrix>
class IntRowDiff : public IRowDiff, public BinaryMatrix, public IntMatrix {
  public:
    static_assert(std::is_convertible<BaseMatrix*, IntMatrix*>::value);

    template <typename... Args>
    IntRowDiff(const graph::DBGSuccinct *graph = nullptr, Args&&... args)
        : diffs_(std::forward<Args>(args)...) { graph_ = graph; }

    std::vector<Row> get_column(Column j) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    // query integer values
    std::vector<RowValues> get_row_values(const std::vector<Row> &rows) const override;

    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    const BinaryMatrix& get_binary_matrix() const override { return *this; }

  private:
    static void decode_diffs(RowValues *diffs);
    static void add_diff(const RowValues &diff, RowValues *row);

    BaseMatrix diffs_;
};


template <class BaseMatrix>
std::vector<BinaryMatrix::Row> IntRowDiff<BaseMatrix>::get_column(Column j) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    const graph::boss::BOSS &boss = graph_->get_boss();
    assert(!fork_succ_.size() || fork_succ_.size() == boss.get_last().size());

    // TODO: implement a more efficient algorithm
    std::vector<Row> result;
    for (Row i = 0; i < num_rows(); ++i) {
        auto edge = graph_->kmer_to_boss_index(
            graph::AnnotatedSequenceGraph::anno_to_graph_index(i)
        );

        if (!boss.get_W(edge))
            continue;

        SetBitPositions set_bits = get_rows({ i })[0];
        if (std::binary_search(set_bits.begin(), set_bits.end(), j))
            result.push_back(i);
    }
    return result;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
IntRowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows;
    rows.reserve(row_ids.size());
    for (const auto &row : get_row_values(row_ids)) {
        rows.emplace_back();
        rows.back().reserve(row.size());
        for (const auto &[j, _] : row) {
            rows.back().push_back(j);
        }
    }
    return rows;
}

template <class BaseMatrix>
std::vector<IntMatrix::RowValues>
IntRowDiff<BaseMatrix>::get_row_values(const std::vector<Row> &row_ids) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");
    assert(!fork_succ_.size() || fork_succ_.size() == graph_->get_boss().get_last().size());

    // get row-diff paths
    auto [rd_ids, rd_paths_trunc, times_traversed] = get_rd_ids(row_ids);

    std::vector<RowValues> rd_rows = diffs_.get_row_values(rd_ids);
    for (auto &row : rd_rows) {
        decode_diffs(&row);
    }

    rd_ids = std::vector<Row>();

    // reconstruct annotation rows from row-diff
    std::vector<RowValues> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        RowValues &result = rows[i];
        // propagate back and reconstruct full annotations for predecessors
        for (auto it = rd_paths_trunc[i].rbegin(); it != rd_paths_trunc[i].rend(); ++it) {
            std::sort(rd_rows[*it].begin(), rd_rows[*it].end());
            add_diff(rd_rows[*it], &result);
            // replace diff row with full reconstructed annotation
            if (--times_traversed[*it]) {
                rd_rows[*it] = result;
            } else {
                // free memory
                rd_rows[*it] = {};
            }
        }
        assert(std::all_of(result.begin(), result.end(),
                           [](auto &p) { return p.second; }));
        assert(std::all_of(result.begin(), result.end(),
                           [](auto &p) { return (int64_t)p.second > 0; }));
    }

    return rows;
}

template <class BaseMatrix>
bool IntRowDiff<BaseMatrix>::load(std::istream &in) {
    std::string version(4, '\0');
    in.read(version.data(), 4);
    return anchor_.load(in) && fork_succ_.load(in) && diffs_.load(in);
}

template <class BaseMatrix>
void IntRowDiff<BaseMatrix>::serialize(std::ostream &out) const {
    out.write("v2.0", 4);
    anchor_.serialize(out);
    fork_succ_.serialize(out);
    diffs_.serialize(out);
}

template <class BaseMatrix>
void IntRowDiff<BaseMatrix>::decode_diffs(RowValues *diffs) {
    for (auto &[j, value] : *diffs) {
        value = decode_diff(value);
    }
}

template <class BaseMatrix>
void IntRowDiff<BaseMatrix>::add_diff(const RowValues &diff, RowValues *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.empty())
        return;

    RowValues result;
    result.reserve(row->size() + diff.size());

    auto it = row->begin();
    auto it2 = diff.begin();
    while (it != row->end() && it2 != diff.end()) {
        if (it->first < it2->first) {
            result.push_back(*it);
            ++it;
        } else if (it->first > it2->first) {
            result.push_back(*it2);
            ++it2;
        } else {
            if (uint64_t sum = it->second + it2->second)
                result.emplace_back(it->first, sum);
            ++it;
            ++it2;
        }
    }
    std::copy(it, row->end(), std::back_inserter(result));
    std::copy(it2, diff.end(), std::back_inserter(result));

    row->swap(result);
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_ROW_DIFF_HPP__
