#include "flat_matrix.hpp"

#include <algorithm>
#include <functional>

#include "common/serialization.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"


namespace mtg {
namespace annot {
namespace binmat {

template <typename BitVector>
RowFlat<BitVector>::
RowFlat(const std::function<void(const RowCallback&)> &call_rows,
        uint64_t num_columns,
        uint64_t num_rows,
        uint64_t num_set_bits)
          : num_columns_(num_columns),
            num_rows_(num_columns ? num_rows : 0) {
    compressed_rows_.reset(new BitVector(
        [&](const std::function<void(uint64_t)> &callback) {
            uint64_t pos = 0;
            call_rows([&](const auto &column_indices) {
                for (const auto &a : column_indices) {
                    callback(pos + a);
                }
                pos += num_columns_;
            });
        }, num_columns_ * num_rows_, num_set_bits)
    );
    num_relations_ = num_set_bits;
}

template <typename BitVector>
bool RowFlat<BitVector>::get(Row row, Column column) const {
    assert(compressed_rows_.get());
    return (*compressed_rows_)[row * num_columns_ + column];
}

template <typename BitVector>
typename RowFlat<BitVector>::SetBitPositions
RowFlat<BitVector>::get_row(Row row) const {
    assert(compressed_rows_.get());
    assert(row * num_columns_ < compressed_rows_->size());
    SetBitPositions columns;
    compressed_rows_->call_ones_in_range(row * num_columns_, (row + 1) * num_columns_,
                                         [&](uint64_t j) { columns.emplace_back(j); });
    return columns;
}

template <typename BitVector>
std::vector<typename RowFlat<BitVector>::Row>
RowFlat<BitVector>::get_column(Column column) const {
    assert(compressed_rows_.get());
    assert(column < num_columns_);
    std::vector<Row> rows;

    for (uint64_t i = column; i < compressed_rows_->size(); i += num_columns_) {
        if ((*compressed_rows_)[i])
            rows.emplace_back(i / num_columns_);
    }
    return rows;
}

template <typename BitVector>
bool RowFlat<BitVector>::load(std::istream &in) {
    compressed_rows_.reset(new BitVector());
    try {
        num_columns_ = load_number(in);
        if (!compressed_rows_->load(in))
            return false;

        num_relations_ = compressed_rows_->num_set_bits();

        // backwards compatibility
        try {
            num_rows_ = load_number(in);
        } catch (...) {
            num_rows_ = compressed_rows_->size() / num_columns_;
        }
    } catch (...) {
        return false;
    }
    return true;
}

template <typename BitVector>
void RowFlat<BitVector>::serialize(std::ostream &out) const {
    assert(compressed_rows_.get());
    serialize_number(out, num_columns_);
    compressed_rows_->serialize(out);
    serialize_number(out, num_rows_);
}

template <typename BitVector>
void RowFlat<BitVector>::serialize(const std::function<void(const RowCallback&)> &call_rows,
                                   uint64_t num_columns,
                                   uint64_t num_rows,
                                   uint64_t num_set_bits,
                                   const std::string &filename, bool append_file) {
    std::ofstream out(filename, std::ios::binary | (append_file ? std::ios::app : 0));
    serialize_number(out, num_columns);
    out.close();
    auto call_bits = [&](auto callback) {
        uint64_t offset = 0;
        call_rows([&](const auto &row) {
            for (auto j : row) {
                callback(offset + j);
            }
            offset += num_columns;
        });
    };
    BitVector::serialize(call_bits, num_columns * num_rows, num_set_bits, filename, true);
    out.open(filename, std::ios::binary | std::ios::app);
    serialize_number(out, num_rows);
}

template class RowFlat<bit_vector_sd>;
template class RowFlat<bit_vector_small>;

} // namespace binmat
} // namespace annot
} // namespace mtg
