#include "eigen_spmat.hpp"

#include "common/serialization.hpp"


namespace mtg {
namespace anno {
namespace binmat {

EigenSpMat::EigenSpMat(uint64_t num_rows, uint64_t max_num_cols) {
    initialize(num_rows, max_num_cols);
}

void EigenSpMat::initialize(size_t num_rows, uint64_t max_num_cols) {
    mat_ = Eigen::SparseMatrix<bool, Eigen::RowMajor>(num_rows, max_num_cols);
    mat_.reserve(num_rows);
}

bool EigenSpMat::get(Row row, Column column) const {
    assert(row < num_rows());
    return mat_.coeff(row, column);
}

void EigenSpMat::set(Row row, Column column) {
    assert(row < num_rows());
    if (!get(row, column))
        mat_.coeffRef(row, column) = 1;

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

EigenSpMat::SetBitPositions
EigenSpMat::get_row(Row row) const {
    assert(row < num_rows());

    SetBitPositions result;
    for (decltype(mat_)::InnerIterator it(mat_, row); it; ++it) {
        result.push_back(it.index());
    }
    return result;
}

std::vector<EigenSpMat::Row>
EigenSpMat::get_column(Column column) const {
    std::vector<Row> result;

    for (uint64_t i = 0, n_rows = num_rows(); i < n_rows; ++i) {
        if (get(i, column))
            result.push_back(i);
    }
    return result;
}

void EigenSpMat::clear_row(Row row) {
    assert(row < num_rows());

    for (decltype(mat_)::InnerIterator it(mat_, row); it; ++it) {
        it.valueRef() = 0;
    }
}

void EigenSpMat::insert_rows(const std::vector<Row> &) {
    throw std::runtime_error("Error: Not implemented");
}

bool EigenSpMat::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        size_t num_rows = load_number(instream);
        num_columns_ = load_number(instream);

        initialize(num_rows, std::max(num_columns_,
                                      static_cast<uint64_t>(mat_.cols())));

        sdsl::int_vector<> full_vector;
        full_vector.load(instream);

        for (size_t k = 0, i = 0; k < full_vector.size(); ++k) {
            if (full_vector[k]) {
                if (full_vector[k] - 1 >= num_columns_)
                    return false;
                mat_.coeffRef(i, full_vector[k] - 1) = true;
            } else {
                i++;
            }
        }

        return true;
    } catch (...) {
        return false;
    }
}

void EigenSpMat::serialize(std::ostream &outstream) const {
    serialize_number(outstream, num_rows());
    serialize_number(outstream, num_columns());

    sdsl::int_vector<> full_vector(num_relations() + num_rows(),
                                   0,
                                   sdsl::bits::hi(num_columns()) + 1);

    for (uint64_t i = 0, p = 0, n_rows = num_rows(); i < n_rows; ++i) {
        for (decltype(mat_)::InnerIterator it(mat_, i); it; ++it) {
            full_vector[p++] = it.index() + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

// number of ones in the matrix
uint64_t EigenSpMat::num_relations() const {
    uint64_t num_set_bits = 0;
    for (uint64_t i = 0, n_rows = num_rows(); i < n_rows; ++i) {
        num_set_bits += mat_.innerVector(i).nonZeros();
    }
    return num_set_bits;
}

// matrix density
double EigenSpMat::density() const {
    return static_cast<double>(num_relations()) / num_columns() / num_rows();
}

} // namespace binmat
} // namespace anno
} // namespace mtg
