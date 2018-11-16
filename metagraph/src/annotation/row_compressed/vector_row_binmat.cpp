#include "vector_row_binmat.hpp"


bool VectorRowBinMat::get(Row row, Column column) const {
    assert(row < vector_.size());
    return std::find(vector_[row].begin(), vector_[row].end(), column)
                != vector_[row].end();
}

void VectorRowBinMat::set(Row row, Column column) {
    assert(row < vector_.size());

    if (!get(row, column))
        vector_[row].push_back(column);

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

std::vector<VectorRowBinMat::Column>
VectorRowBinMat::get_row(Row row) const {
    assert(row < vector_.size());
    const auto &v = vector_[row];
    return std::vector<Column>(v.begin(), v.end());
}

void VectorRowBinMat::clear_row(Row row) {
    assert(row < vector_.size());
    vector_[row].clear();
}

std::vector<VectorRowBinMat::Row>
VectorRowBinMat::get_column(Column column) const {
    std::vector<Row> result;
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        if (get(i, column))
            result.push_back(i);
    }
    return result;
}

void VectorRowBinMat::insert_rows(const std::vector<Row> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));
    utils::insert_default_values(rows, &vector_);
}

bool VectorRowBinMat::load(std::istream &instream) {
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

void VectorRowBinMat::serialize(std::ostream &outstream) const {
    serialize_number(outstream, num_rows());
    serialize_number(outstream, num_columns());

    sdsl::int_vector<> full_vector(num_relations() + num_rows(),
                                   0,
                                   utils::code_length(num_columns()));

    for (uint64_t i = 0, p = 0; i < vector_.size(); ++i) {
        for (uint64_t value : vector_[i]) {
            full_vector[p++] = value + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

// number of ones in the matrix
uint64_t VectorRowBinMat::num_relations() const {
    return std::accumulate(
        vector_.begin(), vector_.end(), uint64_t(0),
        [](uint64_t sum, const auto &v) { return sum + v.size(); }
    );
}

// matrix density
double VectorRowBinMat::density() const {
    return static_cast<double>(num_relations()) / num_columns() / num_rows();
}
