#include "vector_row_binmat.hpp"

#include <cstdio>

#include <sdsl/int_vector_buffer.hpp>

#include "common/serialization.hpp"
#include "common/threads/threading.hpp"
#include "common/algorithms.hpp"


namespace mtg {
namespace anno {
namespace binmat {

template <typename RowType>
VectorRowBinMat<RowType>::VectorRowBinMat(std::vector<RowType>&& rows,
                                          uint64_t num_columns)
      : num_columns_(num_columns), vector_(std::move(rows)) {
    // make sure there are no columns with indexes greater than num_labels
    assert(std::all_of(vector_.begin(), vector_.end(), [&](const auto &row) {
        return std::all_of(row.begin(), row.end(),
                           [num_columns](uint64_t col_id) { return col_id < num_columns; });
    }));
}

template <typename RowType>
bool VectorRowBinMat<RowType>::get(Row row, Column column) const {
    assert(row < vector_.size());
    return std::find(vector_[row].begin(), vector_[row].end(), column)
                != vector_[row].end();
}

template <typename RowType>
void VectorRowBinMat<RowType>::set(Row row, Column column) {
    assert(row < vector_.size());

    if (!get(row, column))
        vector_[row].push_back(column);

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

template <typename RowType>
void VectorRowBinMat<RowType>::force_set(Row row, Column column) {
    assert(row < vector_.size());

    vector_[row].push_back(column);

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

template <typename RowType>
void VectorRowBinMat<RowType>::standardize_rows() {
    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < vector_.size(); ++i) {
        std::sort(vector_[i].begin(), vector_[i].end());
        vector_[i].erase(std::unique(vector_[i].begin(), vector_[i].end()),
                         vector_[i].end());
    }
}

template <typename RowType>
typename VectorRowBinMat<RowType>::SetBitPositions
VectorRowBinMat<RowType>::get_row(Row row) const {
    assert(row < vector_.size());
    if constexpr(std::is_same_v<RowType, SetBitPositions>) {
        return vector_[row];
    } else {
        const auto &v = vector_[row];
        return SetBitPositions(v.begin(), v.end());
    }
}

template <typename RowType>
void VectorRowBinMat<RowType>::clear_row(Row row) {
    assert(row < vector_.size());
    vector_[row].clear();
}

template <typename RowType>
std::vector<typename VectorRowBinMat<RowType>::Row>
VectorRowBinMat<RowType>::get_column(Column column) const {
    std::vector<Row> result;
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        if (get(i, column))
            result.push_back(i);
    }
    return result;
}

template <typename RowType>
void VectorRowBinMat<RowType>::insert_rows(const std::vector<Row> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));
    utils::insert(&vector_, rows, {});
}

template <typename RowType>
bool VectorRowBinMat<RowType>::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        size_t num_rows = load_number(instream);
        num_columns_ = load_number(instream);
        vector_.clear();
        vector_.shrink_to_fit();
        vector_.resize(num_rows);

        sdsl::int_vector<> full_vector;
        full_vector.load(instream);

        for (size_t k = 0, i = 0; k < full_vector.size(); ++k) {
            if (full_vector[k]) {
                vector_[i].push_back(full_vector[k] - 1);
                if (vector_[i].back() >= num_columns_)
                    return false;
            } else {
                i++;
            }
        }

        return true;
    } catch (...) {
        return false;
    }
}

template <typename RowType>
void VectorRowBinMat<RowType>::serialize(std::ostream &outstream) const {
    serialize_number(outstream, num_rows());
    serialize_number(outstream, num_columns());

    sdsl::int_vector<> full_vector(num_relations() + num_rows(),
                                   0,
                                   sdsl::bits::hi(num_columns()) + 1);

    for (uint64_t i = 0, p = 0; i < vector_.size(); ++i) {
        for (uint64_t value : vector_[i]) {
            full_vector[p++] = value + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

// number of ones in the matrix
template <typename RowType>
uint64_t VectorRowBinMat<RowType>::num_relations() const {
    return std::accumulate(
        vector_.begin(), vector_.end(), uint64_t(0),
        [](uint64_t sum, const auto &v) { return sum + v.size(); }
    );
}

// matrix density
template <typename RowType>
double VectorRowBinMat<RowType>::density() const {
    return static_cast<double>(num_relations()) / num_columns() / num_rows();
}

template class VectorRowBinMat<SmallVector<uint32_t>>;
template class VectorRowBinMat<Vector<uint64_t>>;

} // namespace binmat
} // namespace anno
} // namespace mtg
