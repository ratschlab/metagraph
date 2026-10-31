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

// Count matrices with the single-row read of values the budget-aware decode path needs
// (CSCMatrix over BRWT or ColumnMajor: row_values())
template <class M, class = void>
struct HasRowValues : std::false_type {};
template <class M>
struct HasRowValues<M, std::void_t<decltype(std::declval<const M&>().get_binary_matrix().row_column_ranks(
        BinaryMatrix::Row(),
        static_cast<ChargedBuffer<std::pair<BinaryMatrix::Column, uint64_t>>*>(nullptr),
        std::declval<DecodeBudget&>())), decltype(&M::row_values)>> : std::true_type {};

template <class BaseMatrix>
class IntRowDiff : public IRowDiff, public BinaryMatrix, public IntMatrix {
  public:
    static_assert(std::is_convertible<BaseMatrix*, IntMatrix*>::value);

    template <typename... Args>
    IntRowDiff(const graph::DeBruijnGraph *graph = nullptr, Args&&... args)
        : diffs_(std::forward<Args>(args)...) { graph_ = graph; }

    std::vector<Row> get_column(Column column) const override {
        return get_column_from_rows(*this, column);
    }
    std::vector<SetBitPositions> get_rows(std::span<const Row> rows) const override;
    // no deduplication: see class comment on RowDiff::get_rows_dict (speed vs limited size win)
    std::vector<SetBitPositions>
    get_rows_dict(std::vector<Row> *rows, size_t num_threads) const override;
    // query integer values
    std::vector<RowValues> get_row_values(std::span<const Row> rows,
                                          size_t num_threads = 1) const override;

    uint64_t num_columns() const override { return diffs_.num_columns(); }
    uint64_t num_relations() const override { return diffs_.num_relations(); }
    uint64_t num_rows() const override { return diffs_.num_rows(); }

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    const BaseMatrix& diffs() const { return diffs_; }
    BaseMatrix& diffs() { return diffs_; }

    const BinaryMatrix& get_binary_matrix() const override { return *this; }

    // The budget-aware decode path (IRowDiff), over CSCMatrix<BRWT | ColumnMajor, counts>
    // (row_diff_int_brwt). Presence only, as get_rows() returns: whether a column is in a
    // row is whether its count there is non-zero, the sum of the diffs on the row's path, so
    // the stored rows are read with their values (the index matrix's single-row descent and
    // the counts of the entries it finds), the counts reconstructed (the sum rule of
    // add_diff), and their columns returned. No path cache (|cache| is not used): the
    // default decode of counts has none either.
    bool supports_budgeted_decode() const override { return HasRowValues<BaseMatrix>::value; }
    DecodeStatus decode_rows(const std::vector<Row> &rows, DecodeBudget &budget,
                             std::vector<SetBitPositions> *out, std::vector<RowCost> *costs,
                             std::vector<uint64_t> *held,
                             RowDiffPathCache *cache = nullptr) const override;

  private:
    static void decode_diffs(RowValues *diffs);
    static void add_diff(const RowValues &diff, RowValues *row);

    BaseMatrix diffs_;
};

// The stored rows of a count row-diff matrix for the budget-aware decode path: the count
// matrix's single-row read (CSCMatrix::row_values), which charges the row before its values
// are copied, then the stored diffs decoded in place (IntRowDiff::decode_diffs, which the
// default decode applies to every stored row, anchors included)
template <class BaseMatrix>
class BudgetedValuesFetcher : public RowFetcher<IntMatrix::RowValues> {
  public:
    using ColRank = std::pair<BinaryMatrix::Column, uint64_t>;

    explicit BudgetedValuesFetcher(const BaseMatrix &matrix) : matrix_(matrix) {}

    bool fetch(BinaryMatrix::Row r, IntMatrix::RowValues *row, uint64_t *bytes,
               DecodeBudget &budget) override {
        if (!matrix_.row_values(r, &scratch_, row, bytes, budget))
            return false;
        for (auto &[j, value] : *row) {
            value = decode_diff(value);
        }
        return true;
    }
    uint64_t scratch_charged() const override { return scratch_.charged(); }
    uint64_t scratch_peak(uint64_t entries) const override {
        return ChargedBuffer<ColRank>::peak_for(entries);
    }
    uint64_t scratch_entries(const IntMatrix::RowValues &row) const override {
        return row.size();
    }

  private:
    const BaseMatrix &matrix_;
    ChargedBuffer<ColRank> scratch_;
};

template <class BaseMatrix>
DecodeStatus IntRowDiff<BaseMatrix>::decode_rows(const std::vector<Row> &rows,
                                                 DecodeBudget &budget,
                                                 std::vector<SetBitPositions> *out,
                                                 std::vector<RowCost> *costs,
                                                 std::vector<uint64_t> *held,
                                                 RowDiffPathCache *) const {
    if constexpr(HasRowValues<BaseMatrix>::value) {
        const uint64_t at_entry = budget.held();
        BudgetedValuesFetcher<BaseMatrix> fetcher(diffs_);
        std::vector<RowValues> values;
        std::vector<RowCost> row_costs;
        std::vector<uint64_t> values_held;
        const DecodeStatus status = decode_budgeted<RowValues>(rows, budget, fetcher, &values,
                                                               &row_costs, &values_held);
        if (status != DecodeStatus::OK)
            return status;
        // the columns replace the counts row by row, as get_rows() converts them; a row's
        // demand adds what converting it alone holds beside its counts (as
        // TupleRowDiff::decode_rows does for tuples)
        const size_t n = rows.size();
        if (!budget.charge(output_bytes<SetBitPositions>(n))) {
            budget.restore(at_entry);
            return DecodeStatus::REFUSED;
        }
        std::vector<SetBitPositions> result(n);
        std::vector<uint64_t> result_held(n);
        for (size_t i = 0; i < n; ++i) {
            const uint64_t bytes = small_vector_bytes(values[i].size(), sizeof(Column));
            if (!budget.charge(bytes)) {
                budget.restore(at_entry);
                return DecodeStatus::REFUSED;
            }
            result[i] = utils::get_firsts<SetBitPositions>(values[i]);
            result_held[i] = bytes;
            budget.release(values_held[i]);
            RowValues().swap(values[i]);
            row_costs[i].demand += bytes + output_bytes<SetBitPositions>(1);
        }
        budget.release(output_bytes<RowValues>(n));
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
IntRowDiff<BaseMatrix>::get_rows(std::span<const Row> row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    call_rows(row_ids,
        [this](std::span<const Row> rd_ids, size_t num_threads) {
            return diffs_.get_row_values(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowValues &row) { rows[i] = utils::get_firsts<SetBitPositions>(row); },
        1
    );
    return rows;
}

template <class BaseMatrix>
std::vector<BinaryMatrix::SetBitPositions>
IntRowDiff<BaseMatrix>::get_rows_dict(std::vector<Row> *rows, size_t num_threads) const {
    std::vector<SetBitPositions> rows_dict(rows->size());
    call_rows(*rows,
        [this](std::span<const Row> rd_ids, size_t num_threads) {
            return diffs_.get_row_values(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowValues &row) {
            rows_dict[i] = utils::get_firsts<SetBitPositions>(row);
            (*rows)[i] = i;
        },
        num_threads
    );
    return rows_dict;
}

template <class BaseMatrix>
std::vector<IntMatrix::RowValues>
IntRowDiff<BaseMatrix>::get_row_values(std::span<const Row> row_ids, size_t num_threads) const {
    std::vector<RowValues> rows(row_ids.size());
    call_rows(row_ids,
        [this](std::span<const Row> rd_ids, size_t num_threads) {
            return diffs_.get_row_values(rd_ids, num_threads);
        },
        add_diff, decode_diffs,
        [&](size_t i, const RowValues &row) { rows[i] = row; },
        num_threads
    );
    return rows;
}

template <class BaseMatrix>
bool IntRowDiff<BaseMatrix>::load(std::istream &in) {
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
    assert(std::all_of(row->begin(), row->end(),
                       [](auto &p) { return (int64_t)p.second > 0; }));

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

    assert(std::all_of(row->begin(), row->end(),
                       [](auto &p) { return (int64_t)p.second > 0; }));
}

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __INT_ROW_DIFF_HPP__
