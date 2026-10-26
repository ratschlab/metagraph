#include "column_major.hpp"

#include <numeric>

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/serialization.hpp"


namespace mtg {
namespace annot {
namespace matrix {

uint64_t ColumnMajor::num_rows() const {
    if (!columns_.size()) {
        return explicit_num_rows_if_empty_;
    } else {
        assert(columns_[0]);
        return columns_[0]->size();
    }
}

bool ColumnMajor::get(Row row, Column column) const {
    assert(column < columns_.size());
    assert(columns_[column]);
    assert(row < columns_[column]->size());
    return (*columns_[column])[row];
}

std::vector<ColumnMajor::SetBitPositions>
ColumnMajor::get_rows(std::span<const Row> row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());
    if (columns_.empty())
        return rows;
#ifndef NDEBUG
    for (auto row : row_ids)
        assert(row < num_rows());
#endif

    // Group nearby requests, retaining their original positions in the result.
    std::vector<size_t> order(row_ids.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return row_ids[a] > row_ids[b];
    });

    std::vector<size_t> group_ends;
    for (size_t begin = 0; begin < order.size();) {
        const auto base = row_ids[order[begin]] & ~uint64_t(63);
        size_t end = begin + 1;
        while (end < order.size() && row_ids[order[end]] >= base)
            ++end;
        group_ends.push_back(end);
        begin = end;
    }

    for (size_t j = 0; j < columns_.size(); ++j) {
        assert(columns_[j]);
        const bit_vector *column = columns_[j].get();
        if (const auto *adaptive = dynamic_cast<const bit_vector_adaptive*>(column))
            column = &adaptive->data();
        const bit_vector &col = *column;

        size_t begin = 0;
        for (size_t end : group_ends) {
            if (end == begin + 1) {
                const auto i = order[begin];
                if (col[row_ids[i]])
                    rows[i].push_back(j);
            } else {
                const auto first_row = row_ids[order[end - 1]];
                const auto width = row_ids[order[begin]] - first_row + 1;
                const auto bits = col.get_int(first_row, width);
                for (size_t pos = begin; pos < end; ++pos) {
                    const auto i = order[pos];
                    if ((bits >> (row_ids[i] - first_row)) & 1)
                        rows[i].push_back(j);
                }
            }
            begin = end;
        }
    }

    return rows;
}

std::vector<Vector<std::pair<ColumnMajor::Column, uint64_t>>>
ColumnMajor::get_column_ranks(std::span<const Row> rows, size_t num_threads) const {
    using T = Vector<std::pair<Column, uint64_t>>;

    return get_row_data_parallel<T>(rows, num_threads, [&](std::span<const Row> row_ids) {
        std::vector<T> result(row_ids.size());

        for (size_t j = 0; j < columns_.size(); ++j) {
            assert(columns_[j]);
            const bit_vector &col = *columns_[j];

            for (size_t i = 0; i < row_ids.size(); ++i) {
                assert(row_ids[i] < num_rows());

                if (uint64_t r = col.conditional_rank1(row_ids[i]))
                    result[i].emplace_back(j, r);
            }
        }

        return result;
    });
}

void ColumnMajor::call_columns(const std::vector<Column> &columns,
                               const std::function<void(size_t, const bitmap&)> &callback,
                               size_t num_threads) const {
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t j = 0; j < columns.size(); ++j) {
        callback(j, *columns_[columns[j]]);
    }
}

std::vector<ColumnMajor::Row> ColumnMajor::get_column(Column column) const {
    assert(column < columns_.size());
    assert(columns_[column]);

    std::vector<Row> result;
    columns_[column]->call_ones([&result](auto i) { result.push_back(i); });
    return result;
}

bool ColumnMajor::load(std::istream &in) {
    if (!in.good())
        return false;

    columns_.clear();
    explicit_num_rows_if_empty_ = 0;

    try {
        columns_.resize(load_number(in));

        for (auto &c : columns_) {
            assert(!c);

            // TODO: switch to bit_vector_smart?
            c = std::make_unique<bit_vector_sd>();
            if (!c->load(in))
                return false;
        }
        // optional trailing row count for empty matrices
        if (columns_.empty() && in.good() && in.peek() != std::istream::traits_type::eof()) {
            explicit_num_rows_if_empty_ = load_number(in);
        }
        return true;
    } catch (...) {
        return false;
    }
}

void ColumnMajor::serialize(std::ostream &out) const {
    serialize_number(out, columns_.size());

    for (const auto &c : columns_) {
        assert(c);
        // TODO: better serialize in the original format and add a loader that is able
        // to infer the format from the binary data
        if (const auto *sd_ptr = dynamic_cast<const bit_vector_sd *>(c.get())) {
            sd_ptr->serialize(out);
        } else {
            c->copy_to<bit_vector_sd>().serialize(out);
        }
    }
    if (columns_.empty() && explicit_num_rows_if_empty_ > 0) {
        serialize_number(out, explicit_num_rows_if_empty_);
    }
}

// number of ones in the matrix
uint64_t ColumnMajor::num_relations() const {
    uint64_t num_set_bits = 0;

    for (const auto &c : columns_) {
        assert(c);
        num_set_bits += c->num_set_bits();
    }
    return num_set_bits;
}

std::vector<std::pair<ColumnMajor::Column, size_t /* count */>>
ColumnMajor::sum_rows(const std::vector<std::pair<Row, size_t>> &index_counts,
                      size_t min_count) const {
    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<std::pair<Column, size_t>> result;
    result.reserve(num_columns());

    for (size_t j = 0; j < num_columns(); ++j) {
        size_t total_checked = 0;
        size_t total_matched = 0;

        const bit_vector &column = *columns_[j];

        for (auto [i, count] : index_counts) {
            total_checked += count;
            total_matched += count * column[i];

            if (total_matched + (total_sum_count - total_checked) < min_count)
                break;
        }

        if (total_matched >= min_count)
            result.emplace_back(j, total_matched);
    }

    return result;
}

} // namespace matrix
} // namespace annot
} // namespace mtg
