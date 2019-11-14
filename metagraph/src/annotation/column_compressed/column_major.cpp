#include "column_major.hpp"

#include "serialization.hpp"


ColMajorCompressed
::ColMajorCompressed(const std::vector<std::unique_ptr<bit_vector_sd>> &columns) {
    columns_.reserve(columns.size());
    for (auto &column : columns) {
        columns_.emplace_back(new bit_vector_sd(*column));
    }
}

ColMajorCompressed
::ColMajorCompressed(std::vector<std::unique_ptr<bit_vector_sd>>&& columns)
      : columns_(std::move(columns)) {}

uint64_t ColMajorCompressed::num_rows() const {
    if (!columns_.size()) {
        return 0;
    } else {
        assert(columns_[0].get());
        return columns_[0]->size();
    }
}

bool ColMajorCompressed::get(Row row, Column column) const {
    assert(column < columns_.size());
    assert(columns_[column].get());
    assert(row < columns_[column]->size());
    return (*columns_[column])[row];
}

ColMajorCompressed::SetBitPositions
ColMajorCompressed::get_row(Row row) const {
    assert(row < num_rows());

    SetBitPositions result;
    for (size_t i = 0; i < columns_.size(); ++i) {
        assert(columns_[i].get());

        if ((*columns_[i])[row])
            result.push_back(i);
    }
    return result;
}

std::vector<ColMajorCompressed::SetBitPositions>
ColMajorCompressed::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t j = 0; j < columns_.size(); ++j) {
        assert(columns_[j].get());
        const auto &col = (*columns_[j]);

        for (size_t i = 0; i < row_ids.size(); ++i) {
            assert(row_ids[i] < num_rows());

            if (col[row_ids[i]])
                rows[i].push_back(j);
        }
    }

    return rows;
}

std::vector<ColMajorCompressed::Row>
ColMajorCompressed::get_column(Column column) const {
    assert(column < columns_.size());
    assert(columns_[column].get());

    std::vector<Row> result;
    columns_[column]->call_ones([&result](auto i) { result.push_back(i); });
    return result;
}

bool ColMajorCompressed::load(std::istream &in) {
    if (!in.good())
        return false;

    columns_.clear();

    try {
        columns_.resize(load_number(in));

        for (auto &column : columns_) {
            assert(!column.get());

            auto next = std::make_unique<bit_vector_sd>();
            if (!next->load(in))
                return false;
            column = std::move(next);
        }
        return true;
    } catch (...) {
        return false;
    }
}

void ColMajorCompressed::serialize(std::ostream &out) const {
    serialize_number(out, columns_.size());

    for (const auto &column : columns_) {
        assert(column.get());
        column->serialize(out);
    }
}

// number of ones in the matrix
uint64_t ColMajorCompressed::num_relations() const {
    uint64_t num_set_bits = 0;

    for (const auto &column : columns_) {
        assert(column.get());
        num_set_bits += column->num_set_bits();
    }
    return num_set_bits;
}
