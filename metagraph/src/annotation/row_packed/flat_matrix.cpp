#include "flat_matrix.hpp"

#include <algorithm>
#include <functional>

#include "utils.hpp"
#include "serialization.hpp"


template <typename BitVector>
RowConcatenated<BitVector>::
RowConcatenated(const std::function<void(RowCallback)> &call_rows,
                uint64_t num_columns,
                uint64_t num_rows,
                uint64_t num_set_bits)
      : num_columns_(num_columns) {
    compressed_rows_.reset(new BitVector(
        [&](const std::function<void(uint64_t)> &callback) {
            uint64_t pos = 0;
            call_rows([&](const auto &column_indices) {
                for (const auto &a : column_indices) {
                    callback(pos + a);
                }
                pos += num_columns_;
            });
        }, num_columns_ * num_rows, num_set_bits)
    );
}

template <typename BitVector>
bool RowConcatenated<BitVector>::get(Row row, Column column) const {
    return (*compressed_rows_)[row * num_columns_ + column];
}

template <typename BitVector>
std::vector<typename RowConcatenated<BitVector>::Column>
RowConcatenated<BitVector>::get_row(Row row) const {
    assert(row * num_columns_ < compressed_rows_->size());
    std::vector<Column> columns;

    uint64_t offset = row * num_columns_;
    uint64_t first = compressed_rows_->rank1(offset) + !(*compressed_rows_)[offset];
    uint64_t last = compressed_rows_->rank1(offset + num_columns_ - 1);
    for (uint64_t i = first; i <= last; ++i) {
        columns.emplace_back(compressed_rows_->select1(i) - offset);
    }
    return columns;
}

template <typename BitVector>
std::vector<typename RowConcatenated<BitVector>::Row>
RowConcatenated<BitVector>::get_column(Column column) const {
    assert(column < num_columns_);
    std::vector<Row> rows;

    for (uint64_t i = column; i < compressed_rows_->size(); i += num_columns_) {
        if ((*compressed_rows_)[i])
            rows.emplace_back(i / num_columns_);
    }
    return rows;
}

template <typename BitVector>
bool RowConcatenated<BitVector>::load(std::istream &in) {
    compressed_rows_.reset(new BitVector());
    try {
        num_columns_ = load_number(in);
        compressed_rows_->load(in);
    } catch (...) {
        return false;
    }
    return true;
}

template <typename BitVector>
void RowConcatenated<BitVector>::serialize(std::ostream &out) const {
    serialize_number(out, num_columns_);
    compressed_rows_->serialize(out);
}

template class RowConcatenated<bit_vector_sd>;
template class RowConcatenated<bit_vector_small>;
