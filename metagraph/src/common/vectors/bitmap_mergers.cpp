#include "bitmap_mergers.hpp"

#include <cassert>

#include "common/vector.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/vectors/bit_vector_dyn.hpp"


namespace utils {

RowsFromColumnsTransformer::VectorBitStream
::VectorBitStream(const bit_vector &vector,
                  uint64_t begin,
                  uint64_t end)
      : vector_(vector),
        begin_(begin),
        current_rank_(begin ? vector_.rank1(begin - 1) : 0) {
    assert(begin <= end);
    assert(current_rank_ <= begin);
    max_rank_ = (end == static_cast<uint64_t>(-1)
                    ? vector_.num_set_bits()
                    : (end ? vector_.rank1(end - 1) : 0));
}

uint64_t RowsFromColumnsTransformer::VectorBitStream::next_value() {
    assert(current_rank_ < max_rank_);
    return vector_.select1(++current_rank_) - begin_;
}

RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<const bit_vector*> &columns) {
    initialize(columns);
}

RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::unique_ptr<bit_vector>> &columns) {
    std::vector<const bit_vector*> cols;
    for (const auto &col_ptr : columns) {
        cols.push_back(col_ptr.get());
    }
    initialize(cols);
}

RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const std::vector<std::shared_ptr<bit_vector>> &columns) {
    std::vector<const bit_vector*> cols;
    for (const auto &col_ptr : columns) {
        cols.push_back(col_ptr.get());
    }
    initialize(cols);
}

RowsFromColumnsTransformer
::RowsFromColumnsTransformer(const bit_vector &columns_concatenated,
                             uint64_t column_size) {
    assert(columns_concatenated.size() % column_size == 0);

    if (!columns_concatenated.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = column_size;
    uint64_t num_columns = columns_concatenated.size() / column_size;
    streams_.reserve(num_columns);
    for (uint64_t i = 0; i < num_columns; ++i) {
        streams_.emplace_back(new VectorBitStream(columns_concatenated,
                                                  column_size * i,
                                                  column_size * (i + 1)));
    }

    init_heap();
}

void RowsFromColumnsTransformer
::initialize(const std::vector<const bit_vector*> &columns) {
    if (!columns.size()) {
        num_rows_ = 0;
        return;
    }

    num_rows_ = columns.at(0)->size();
    for (const auto &column : columns) {
        if (column->size() != num_rows_)
            throw std::runtime_error("Error: columns are not of the same size");
    }

    streams_.reserve(columns.size());
    for (const auto &col_ptr : columns) {
        streams_.emplace_back(new VectorBitStream(*col_ptr));
    }

    init_heap();
}

void RowsFromColumnsTransformer::init_heap() {
    for (uint64_t column = 0; column < streams_.size(); ++column) {
        assert(streams_[column].get());

        num_set_bits_left_ += streams_[column]->values_left();

        if (streams_[column]->values_left()) {
            index_heap_.emplace(kmer_label_pair {
                streams_[column]->next_value(),
                column
            });
        }
    }
}

void RowsFromColumnsTransformer::call_next(ValueCallback callback) {
    assert(values_left() && index_heap_.size());

    auto index = std::move(index_heap_.top());

    index_heap_.pop();
    num_set_bits_left_--;

    // TODO: Important!
    // Load multiple values from the stream to improve cache efficiency
    if (streams_[index.col_id]->values_left()) {
        index_heap_.emplace(kmer_label_pair {
            streams_[index.col_id]->next_value(),
            index.col_id
        });
    }

    assert(index_heap_.size() || !values_left());
    callback(index.row_id, index.col_id);
}

std::tuple<uint64_t, uint64_t> RowsFromColumnsIterator::next_set_bit() {
    uint64_t row;
    uint64_t column;
    transformer_->call_next([&](uint64_t row_, uint64_t column_) {
        row = row_;
        column = column_;
    });
    return std::make_tuple(row, column);
}

template <typename Vector>
Vector RowsFromColumnsIterator::next_row() {
    Vector indices;

    if (i_ > 0 && (row_ == i_)) {
        indices.push_back(column_);
    }

    if (!values_left() || row_ > i_) {
        i_++;
        return indices;
    }

    while (values_left()) {
        std::tie(row_, column_) = next_set_bit();
        if (row_ != i_)
            break;
        indices.push_back(column_);
    }
    i_++;
    return indices;
}

template Vector<uint64_t> RowsFromColumnsIterator::next_row<Vector<uint64_t>>();


template <typename Vector>
void RowsFromColumnsTransformer
::call_rows(const std::function<void(const Vector &)> &callback) {
    uint64_t cur_row = 0;
    Vector indices;

    while (values_left()) {
        call_next([&](uint64_t row, uint64_t column) {
            while (cur_row < row) {
                callback(indices);
                indices.clear();
                cur_row++;
            }
            indices.push_back(column);
        });
    }

    while (cur_row++ < rows()) {
        callback(indices);
        indices.clear();
    }

    assert(!values_left());
}

template void RowsFromColumnsTransformer::call_rows<Vector<uint64_t>>(
        const std::function<void(const Vector<uint64_t> &)> &callback);

template void RowsFromColumnsTransformer::call_rows<SmallVector<uint32_t>>(
        const std::function<void(const SmallVector<uint32_t> &)> &callback);


template <class BitVectorType>
std::vector<std::unique_ptr<bit_vector>>
transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix) {
    std::vector<std::unique_ptr<bit_vector>> transposed;
    if (!matrix.size())
        return transposed;

    uint64_t num_rows = matrix.size();
    uint64_t num_columns = matrix[0]->size();
    transposed.reserve(num_columns);

    RowsFromColumnsTransformer transformer(matrix);
    //TODO: use call_next directly, without creating Vector<uint64_t>
    transformer.call_rows<Vector<uint64_t>>(
        [&](const Vector<uint64_t> &column_indices) {
            sdsl::bit_vector bv(num_rows, false);
            for (const auto &row_id : column_indices) {
                bv[row_id] = true;
            }
            transposed.emplace_back(new BitVectorType(std::move(bv)));
        }
    );

    return transposed;
}

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_stat>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_small>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

template
std::vector<std::unique_ptr<bit_vector>>
transpose<bit_vector_dyn>(const std::vector<std::unique_ptr<bit_vector>> &matrix);

}
