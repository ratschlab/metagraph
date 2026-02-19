#ifndef __TUPLE_ROW_DIFF_HPP__
#define __TUPLE_ROW_DIFF_HPP__

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

template <class BaseMatrix>
class TupleRowDiff : public IRowDiff, public BinaryMatrix, public MultiIntMatrix {
  public:
    static_assert(std::is_convertible<BaseMatrix*, MultiIntMatrix*>::value);
    static const int SHIFT = 1; // coordinates increase by 1 at each edge

    template <typename... Args>
    TupleRowDiff(const graph::DeBruijnGraph *graph = nullptr, Args&&... args)
        : diffs_(std::forward<Args>(args)...) { graph_ = graph; }

    std::vector<Row> get_column(Column j) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<SetBitPositions>
    get_rows_dict(std::vector<Row> *rows, size_t num_threads) const override;
    std::vector<RowTuples> get_row_tuples(const std::vector<Row> &rows,
                                          size_t num_threads = 1) const override;

    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_attributes() const override { return diffs_.num_attributes(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    const BinaryMatrix& get_binary_matrix() const override { return *this; }

  private:
    static void decode_diffs(RowTuples *diffs);
    static void add_diff(const RowTuples &diff, RowTuples *row);

    BaseMatrix diffs_;
};


template <class BaseMatrix>
std::vector<BinaryMatrix::Row> TupleRowDiff<BaseMatrix>::get_column(Column j) const {
    assert(graph_ && "graph must be loaded");
    assert(anchor_.size() == diffs_.num_rows() && "anchors must be loaded");

    assert(!fork_succ_.size() || fork_succ_.size() == graph_->max_index() + 1);

    // TODO: implement a more efficient algorithm
    std::vector<Row> result;
    graph_->call_nodes([&](auto node) {
        auto i = graph::AnnotatedDBG::graph_to_anno_index(node);
        SetBitPositions set_bits = get_rows({ i })[0];
        if (std::binary_search(set_bits.begin(), set_bits.end(), j))
            result.push_back(i);
    });
    return result;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
TupleRowDiff<BaseMatrix>::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows;
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](const RowTuples &row) { rows.push_back(utils::get_firsts<SetBitPositions>(row)); },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
TupleRowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    VectorSet<SetBitPositions, utils::VectorHash> unique_rows;
    size_t i = 0;
    call_rows(*rows,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](const RowTuples &row) {
            auto it = unique_rows.emplace(utils::get_firsts<SetBitPositions>(row)).first;
            (*rows)[i++] = it - unique_rows.begin();
        },
        num_threads
    );
    return to_vector(std::move(unique_rows));
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::RowTuples>
TupleRowDiff<BaseMatrix>::get_row_tuples(const std::vector<Row> &row_ids, size_t num_threads) const {
    std::vector<RowTuples> rows;
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](const RowTuples &row) { rows.push_back(row); },
        num_threads
    );
    return rows;
}

template <class BaseMatrix>
bool TupleRowDiff<BaseMatrix>::load(std::istream &in) {
    std::string version(4, '\0');
    in.read(version.data(), 4);
    return anchor_.load(in) && fork_succ_.load(in) && diffs_.load(in);
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::serialize(std::ostream &out) const {
    out.write("v2.0", 4);
    anchor_.serialize(out);
    fork_succ_.serialize(out);
    diffs_.serialize(out);
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::decode_diffs(RowTuples *diffs) {
    std::ignore = diffs;
    // no encoding
}

template <class BaseMatrix>
void TupleRowDiff<BaseMatrix>::add_diff(const RowTuples &diff, RowTuples *row) {
    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::is_sorted(diff.begin(), diff.end()));

    if (diff.size()) {
        RowTuples result;
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
                if (it2->second.size()) {
                    result.emplace_back(it->first, Tuple{});
                    std::set_symmetric_difference(it->second.begin(), it->second.end(),
                                                  it2->second.begin(), it2->second.end(),
                                                  std::back_inserter(result.back().second));
                }
                ++it;
                ++it2;
            }
        }
        std::copy(it, row->end(), std::back_inserter(result));
        std::copy(it2, diff.end(), std::back_inserter(result));

        row->swap(result);
    }

    assert(std::is_sorted(row->begin(), row->end()));
    assert(std::all_of(row->begin(), row->end(),
                       [](auto &p) { return p.second.size(); }));
    for (auto &[j, tuple] : *row) {
        assert(std::is_sorted(tuple.begin(), tuple.end()));
        for (uint64_t &c : tuple) {
            c -= SHIFT;
        }
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_ROW_DIFF_HPP__
