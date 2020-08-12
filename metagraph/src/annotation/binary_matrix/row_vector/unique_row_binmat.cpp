#include "unique_row_binmat.hpp"

#include <sdsl/int_vector.hpp>
#include <tsl/ordered_set.h>

#include "common/hashers/hash.hpp"
#include "common/algorithms.hpp"
#include "common/serialization.hpp"


namespace mtg {
namespace annot {
namespace binmat {

UniqueRowBinmat::UniqueRowBinmat(uint64_t num_rows)
      : unique_rows_(1), row_rank_(num_rows, 0) {}

UniqueRowBinmat::UniqueRowBinmat(std::vector<SetBitPositions>&& unique_rows,
                                 std::vector<uint32_t>&& row_rank,
                                 uint32_t num_columns)
      : num_columns_(num_columns),
        unique_rows_(std::move(unique_rows)),
        row_rank_(std::move(row_rank)) {
    // make sure there are no columns with indexes greater than num_labels
    assert(std::all_of(unique_rows_.begin(), unique_rows_.end(), [&](const auto &row) {
        return std::all_of(row.begin(), row.end(),
                           [num_columns](uint32_t j) { return j < num_columns; });
    }));

    for (uint32_t r : row_rank_) {
        assert(r < unique_rows_.size());
        num_relations_ += unique_rows_[r].size();
    }
}

UniqueRowBinmat
::UniqueRowBinmat(const std::function<void(const RowCallback &)> &call_rows,
                  uint32_t num_columns) {
    num_columns_ = num_columns;

    using RowSet = tsl::ordered_set<SetBitPositions,
                                    utils::VectorHash,
                                    std::equal_to<SetBitPositions>,
                                    std::allocator<SetBitPositions>,
                                    std::vector<SetBitPositions>,
                                    uint32_t>;
    RowSet unique_rows;

    call_rows([&](const SetBitPositions &row) {
        num_relations_ += row.size();
        auto it = unique_rows.emplace(row).first;
        row_rank_.push_back(it - unique_rows.begin());
        if (unique_rows.size() == std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("There must be less than 2^32 unique rows");
    });

    unique_rows_ = const_cast<std::vector<SetBitPositions>&&>(
        unique_rows.values_container()
    );
}

bool UniqueRowBinmat::get(Row i, Column j) const {
    assert(i < row_rank_.size());
    assert(row_rank_[i] < unique_rows_.size());
    const auto &row = unique_rows_[row_rank_[i]];
    return std::find(row.begin(), row.end(), j) != row.end();
}

UniqueRowBinmat::SetBitPositions UniqueRowBinmat::get_row(Row i) const {
    assert(i < row_rank_.size());
    assert(row_rank_[i] < unique_rows_.size());
    return unique_rows_[row_rank_[i]];
}

std::vector<BinaryMatrix::SetBitPositions>
UniqueRowBinmat::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows;
    std::vector<size_t> ranks;
    rows.reserve(row_ids.size());
    ranks.reserve(row_ids.size());
    for (Row i : row_ids) {
        ranks.push_back(row_rank_[i]);
    }
    for (size_t r : ranks) {
        rows.push_back(unique_rows_[r]);
    }
    return rows;
}

std::vector<UniqueRowBinmat::Row> UniqueRowBinmat::get_column(Column j) const {
    // first, find all unique rows with `1` in the j-th column
    tsl::ordered_set<uint32_t> row_ranks;
    for (uint32_t r = 0; r < unique_rows_.size(); ++r) {
        const auto &row = unique_rows_[r];
        if (std::find(row.begin(), row.end(), j) != row.end())
            row_ranks.insert(r);
    }
    // fetch all rows that point to the unique rows collected above
    std::vector<Row> result;
    for (uint64_t i = 0; i < row_rank_.size(); ++i) {
        if (row_ranks.count(row_rank_[i]))
            result.push_back(i);
    }
    return result;
}

bool UniqueRowBinmat::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        size_t num_unique_rows = load_number(instream);
        num_columns_ = load_number(instream);
        num_relations_ = load_number(instream);
        unique_rows_.clear();
        unique_rows_.shrink_to_fit();
        unique_rows_.resize(num_unique_rows);

        sdsl::int_vector<> full_vector;
        full_vector.load(instream);

        for (size_t k = 0, i = 0; k < full_vector.size(); ++k) {
            if (full_vector[k]) {
                unique_rows_[i].push_back(full_vector[k] - 1);
                if (unique_rows_[i].back() >= num_columns_)
                    return false;
            } else {
                i++;
            }
        }

        return load_number_vector(instream, &row_rank_);
    } catch (...) {
        return false;
    }
}

void UniqueRowBinmat::serialize(std::ostream &outstream) const {
    serialize_number(outstream, unique_rows_.size());
    serialize_number(outstream, num_columns_);
    serialize_number(outstream, num_relations_);

    uint64_t unique_rows_size = unique_rows_.size();
    for (const auto &row : unique_rows_) {
        unique_rows_size += row.size();
    }

    sdsl::int_vector<> full_vector(unique_rows_size, 0,
                                   sdsl::bits::hi(num_columns()) + 1);

    for (uint64_t i = 0, p = 0; i < unique_rows_.size(); ++i) {
        for (uint64_t value : unique_rows_[i]) {
            full_vector[p++] = value + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);

    serialize_number_vector(outstream, row_rank_);
}

// matrix density
double UniqueRowBinmat::density() const {
    return static_cast<double>(num_relations()) / num_columns() / num_rows();
}

} // namespace binmat
} // namespace annot
} // namespace mtg
