#include "column_major.hpp"

#include "common/vectors/bit_vector_sd.hpp"
#include "common/serialization.hpp"


namespace mtg {
namespace annot {
namespace binmat {

uint64_t ColumnMajor::num_rows() const {
    if (!columns_.size()) {
        return 0;
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

ColumnMajor::SetBitPositions ColumnMajor::get_row(Row row) const {
    assert(row < num_rows() || !columns_.size());

    SetBitPositions result;
    for (size_t i = 0; i < columns_.size(); ++i) {
        assert(columns_[i]);

        if ((*columns_[i])[row])
            result.push_back(i);
    }
    return result;
}

std::vector<ColumnMajor::SetBitPositions>
ColumnMajor::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t j = 0; j < columns_.size(); ++j) {
        assert(columns_[j]);
        const bit_vector &col = *columns_[j];

        for (size_t i = 0; i < row_ids.size(); ++i) {
            assert(row_ids[i] < num_rows());

            if (col[row_ids[i]])
                rows[i].push_back(j);
        }
    }

    return rows;
}

Vector<std::pair<ColumnMajor::Column, uint64_t>>
ColumnMajor::get_column_ranks(Row row) const {
    assert(row < num_rows() || !columns_.size());

    Vector<std::pair<Column, uint64_t>> result;
    for (size_t i = 0; i < columns_.size(); ++i) {
        assert(columns_[i]);

        if (uint64_t r = columns_[i]->conditional_rank1(row))
            result.emplace_back(i, r);
    }
    return result;
}

std::vector<Vector<std::pair<ColumnMajor::Column, uint64_t>>>
ColumnMajor::get_column_ranks(const std::vector<Row> &row_ids) const {
    std::vector<Vector<std::pair<Column, uint64_t>>> result(row_ids.size());

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
}

std::vector<ColumnMajor::Column>
ColumnMajor::slice_rows(const std::vector<Row> &row_ids) const {
    std::vector<Column> slice;

    for (const auto &row : get_rows(row_ids)) {
        for (uint64_t j : row) {
            slice.push_back(j);
        }
        slice.push_back(std::numeric_limits<Column>::max());
    }

    return slice;
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

    try {
        columns_.resize(load_number(in));

        for (auto &c : columns_) {
            assert(!c);

            // TODO: switch to bit_vector_smart?
            c = std::make_unique<bit_vector_sd>();
            if (!c->load(in))
                return false;
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
                       size_t min_count,
                       size_t count_cap) const {
    assert(count_cap >= min_count);

    if (!count_cap)
        return {};

    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<std::pair<uint64_t, size_t>> result;
    result.reserve(num_columns());

    for (size_t j = 0; j < num_columns(); ++j) {
        size_t total_checked = 0;
        size_t total_matched = 0;

        const bit_vector &column = *columns_[j];

        for (auto [i, count] : index_counts) {
            total_checked += count;
            total_matched += count * column[i];

            if (total_matched >= count_cap
                    || total_matched + (total_sum_count - total_checked) < min_count)
                break;
        }

        if (total_matched >= min_count)
            result.emplace_back(j, std::min(total_matched, count_cap));
    }

    return result;
}

} // namespace binmat
} // namespace annot
} // namespace mtg
