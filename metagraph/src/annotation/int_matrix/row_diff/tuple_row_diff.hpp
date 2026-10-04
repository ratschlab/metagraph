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
    // no deduplication: see class comment on RowDiff::get_rows_dict (speed vs limited size win)
    std::vector<SetBitPositions>
    get_rows_dict(std::vector<Row> *rows, size_t num_threads) const override;
    std::vector<RowValues> get_row_values(const std::vector<Row> &rows,
                                          size_t num_threads = 1) const override;
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

    // the budget-aware decode path (IRowDiff), over TupleCSCMatrix<BRWT | ColumnMajor>;
    // decode_rows() decodes the tuples and keeps their columns, as get_rows() does
    bool supports_budgeted_decode() const override { return HasRowTuples<BaseMatrix>::value; }
    DecodeStatus decode_rows(const std::vector<Row> &rows, DecodeBudget &budget,
                             std::vector<SetBitPositions> *out, std::vector<RowCost> *costs,
                             std::vector<uint64_t> *held,
                             RowDiffPathCache *cache = nullptr) const override;
    DecodeStatus decode_row_tuples(const std::vector<Row> &rows, DecodeBudget &budget,
                                   std::vector<RowTuples> *out, std::vector<RowCost> *costs,
                                   std::vector<uint64_t> *held,
                                   RowDiffPathCache *cache = nullptr) const override;

    // the default decode with the path cache (IRowDiff): the cache holds full tuple rows,
    // from which get_rows_cached() keeps the columns, as get_rows() does
    bool supports_path_cache() const override { return true; }
    std::vector<SetBitPositions> get_rows_cached(const std::vector<Row> &rows,
                                                 RowDiffPathCache &cache) const override;
    std::vector<RowTuples> get_row_tuples_cached(const std::vector<Row> &rows,
                                                 RowDiffPathCache &cache) const override;

  private:
    static void decode_diffs(RowTuples *diffs);
    static void add_diff(const RowTuples &diff, RowTuples *row);

    BaseMatrix diffs_;
};

// The stored rows of a coordinate row-diff matrix for the budget-aware decode path: the
// tuple matrix's single-row read (TupleCSCMatrix::row_tuples), which charges each tuple
// before its coordinates are copied
template <class BaseMatrix>
class BudgetedTupleFetcher : public RowFetcher<MultiIntMatrix::RowTuples> {
  public:
    using ColRank = std::pair<BinaryMatrix::Column, uint64_t>;

    explicit BudgetedTupleFetcher(const BaseMatrix &matrix) : matrix_(matrix) {}

    bool fetch(BinaryMatrix::Row r, MultiIntMatrix::RowTuples *row, uint64_t *bytes,
               DecodeBudget &budget) override {
        return matrix_.row_tuples(r, &scratch_, row, bytes, budget);
    }
    uint64_t scratch_charged() const override { return scratch_.charged(); }
    uint64_t scratch_peak(uint64_t entries) const override {
        return ChargedBuffer<ColRank>::peak_for(entries);
    }
    uint64_t scratch_entries(const MultiIntMatrix::RowTuples &row) const override {
        return row.size();
    }

  private:
    const BaseMatrix &matrix_;
    ChargedBuffer<ColRank> scratch_;
};

template <class BaseMatrix>
DecodeStatus TupleRowDiff<BaseMatrix>::decode_row_tuples(const std::vector<Row> &rows,
                                                         DecodeBudget &budget,
                                                         std::vector<RowTuples> *out,
                                                         std::vector<RowCost> *costs,
                                                         std::vector<uint64_t> *held,
                                                         RowDiffPathCache *cache) const {
    if constexpr(HasRowTuples<BaseMatrix>::value) {
        BudgetedTupleFetcher<BaseMatrix> fetcher(diffs_);
        return decode_budgeted<RowTuples>(rows, budget, fetcher, out, costs, held,
                                          cache ? &cache->tuples : nullptr);
    } else {
        return DecodeStatus::UNSUPPORTED;
    }
}

template <class BaseMatrix>
DecodeStatus TupleRowDiff<BaseMatrix>::decode_rows(const std::vector<Row> &rows,
                                                   DecodeBudget &budget,
                                                   std::vector<SetBitPositions> *out,
                                                   std::vector<RowCost> *costs,
                                                   std::vector<uint64_t> *held,
                                                   RowDiffPathCache *cache) const {
    if constexpr(HasRowTuples<BaseMatrix>::value) {
        const uint64_t at_entry = budget.held();
        std::vector<RowTuples> tuples;
        std::vector<RowCost> row_costs;
        std::vector<uint64_t> tuples_held;
        const DecodeStatus status = decode_row_tuples(rows, budget, &tuples, &row_costs,
                                                      &tuples_held, cache);
        if (status != DecodeStatus::OK)
            return status;
        // the columns replace the tuples row by row, as get_rows() converts them; a row's
        // demand adds what converting it alone holds beside its tuples
        const size_t n = rows.size();
        if (!budget.charge(output_bytes<SetBitPositions>(n))) {
            budget.restore(at_entry);
            return DecodeStatus::REFUSED;
        }
        std::vector<SetBitPositions> result(n);
        std::vector<uint64_t> result_held(n);
        for (size_t i = 0; i < n; ++i) {
            const uint64_t bytes = small_vector_bytes(tuples[i].size(), sizeof(Column));
            if (!budget.charge(bytes)) {
                budget.restore(at_entry);
                return DecodeStatus::REFUSED;
            }
            result[i] = utils::get_firsts<SetBitPositions>(tuples[i]);
            result_held[i] = bytes;
            budget.release(tuples_held[i]);
            RowTuples().swap(tuples[i]);
            row_costs[i].demand += bytes + output_bytes<SetBitPositions>(1);
        }
        budget.release(output_bytes<RowTuples>(n));
        out->swap(result);
        costs->swap(row_costs);
        held->swap(result_held);
        return DecodeStatus::OK;
    } else {
        return DecodeStatus::UNSUPPORTED;
    }
}


template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
TupleRowDiff<BaseMatrix>::get_rows_cached(const std::vector<Row> &row_ids,
                                          RowDiffPathCache &cache) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows_cached<RowTuples>(row_ids, cache.tuples,
        [this](const std::vector<Row> &rd_ids) { return diffs_.get_row_tuples(rd_ids, 1); },
        add_diff, decode_diffs,
        [&](size_t i, RowTuples &&row) { rows[i] = utils::get_firsts<SetBitPositions>(row); }
    );
    return rows;
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::RowTuples>
TupleRowDiff<BaseMatrix>::get_row_tuples_cached(const std::vector<Row> &row_ids,
                                                RowDiffPathCache &cache) const {
    std::vector<RowTuples> rows(row_ids.size());
    call_rows_cached<RowTuples>(row_ids, cache.tuples,
        [this](const std::vector<Row> &rd_ids) { return diffs_.get_row_tuples(rd_ids, 1); },
        add_diff, decode_diffs,
        [&](size_t i, RowTuples &&row) { rows[i] = std::move(row); }
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::Row> TupleRowDiff<BaseMatrix>::get_column(Column j) const {
    assert(graph_ && "graph must be loaded");
    assert(diffs_.num_rows() == graph_->max_index());
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
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowTuples &row) { rows[i] = utils::get_firsts<SetBitPositions>(row); },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
TupleRowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    std::vector<SetBitPositions> rows_dict(rows->size());
    call_rows(*rows,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowTuples &row) {
            rows_dict[i] = utils::get_firsts<SetBitPositions>(row);
            (*rows)[i] = i;
        },
        num_threads
    );
    return rows_dict;
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::RowValues>
TupleRowDiff<BaseMatrix>::get_row_values(const std::vector<Row> &row_ids, size_t num_threads) const {
    std::vector<RowValues> rows(row_ids.size());
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowTuples &row) {
            RowValues &row_values = rows[i];
            row_values.reserve(row.size());
            for (const auto &[j, tuple] : row) {
                row_values.emplace_back(j, tuple.size());
            }
        },
        num_threads
    );
    return rows;
}

template <class BaseMatrix>
std::vector<MultiIntMatrix::RowTuples>
TupleRowDiff<BaseMatrix>::get_row_tuples(const std::vector<Row> &row_ids, size_t num_threads) const {
    std::vector<RowTuples> rows(row_ids.size());
    call_rows(row_ids,
        [this](const std::vector<Row> &rd_ids, size_t num_threads) {
            return diffs_.get_row_tuples(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowTuples &row) { rows[i] = row; },
        num_threads
    );
    return rows;
}

template <class BaseMatrix>
bool TupleRowDiff<BaseMatrix>::load(std::istream &in) {
    std::string version(4, '\0');
    in.read(version.data(), 4);
    auto anchor_start = in.tellg();
    if (!anchor_.load(in) || !fork_succ_.load(in))
        return false;
    // anchor_ / fork_succ_ are accessed randomly, hint the kernel.
    utils::madvise_random_range(in, anchor_start, in.tellg() - anchor_start);
    return diffs_.load(in);
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
                    // just for safety, normally rows without coordinates shouldn't be annotated
                    if (result.back().second.empty())
                        result.pop_back();
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
        for (uint64_t &c : tuple) {
            c -= SHIFT;
        }
        assert(std::is_sorted(tuple.begin(), tuple.end()));
    }
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __TUPLE_ROW_DIFF_HPP__
