#ifndef __HYBRID_MATRIX__
#define __HYBRID_MATRIX__

namespace mtg {
namespace annot {
namespace binmat {

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/logger.hpp"

template<class BaseMatrix1, class BaseMatrix2>
class HybridMatrix : public BinaryMatrix {
public:
    typedef bit_vector_smallrank stored_in_matrix1_bv_type;

    HybridMatrix() = default;

    uint64_t num_columns() const override {
        return matrix1_.num_columns();
    }

    uint64_t num_rows() const override {
        return matrix1_.num_rows() + matrix2_.num_rows();
    }

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<Row> get_column(Column column) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override {
        return matrix1_.num_relations() + matrix2_.num_relations();
    }

private:
    BaseMatrix1 matrix1_; //dense
    BaseMatrix2 matrix2_; //sparse
    stored_in_matrix1_bv_type in_matrix1_;
};

template<class BaseMatrix1, class BaseMatrix2>
bool HybridMatrix<BaseMatrix1, BaseMatrix2>::get(Row row, Column column) const {
    assert(column < num_columns());
    assert(row < num_rows());
    if (in_matrix1_[row]) {
        Row row_in_matrix1 = in_matrix1_.rank1(row) - 1;
        return matrix1_.get(row_in_matrix1, column);
    } else {
        Row row_in_matrix2 = in_matrix1_.rank0(row) - 1;
        return matrix2_.get(row_in_matrix2, column);
    }
}

template<class BaseMatrix1, class BaseMatrix2>
BinaryMatrix::SetBitPositions
HybridMatrix<BaseMatrix1, BaseMatrix2>::get_row(Row row) const {
    if (in_matrix1_[row]) {
        Row row_in_matrix1 = in_matrix1_.rank1(row) - 1;
        return matrix1_.get_row(row_in_matrix1);
    } else {
        Row row_in_matrix2 = in_matrix1_.rank0(row) - 1;
        return matrix2_.get_row(row_in_matrix2);
    }
}

template<class BaseMatrix1, class BaseMatrix2>
std::vector<BinaryMatrix::SetBitPositions>
HybridMatrix<BaseMatrix1, BaseMatrix2>::
get_rows(const std::vector<Row> &row_ids) const {
    std::vector<Row> rows_matrix1;
    std::vector<Row> rows_matrix2;

    for (size_t i = 0 ; i < row_ids.size(); ++i) {
        if (in_matrix1_[row_ids[i]])
            rows_matrix1.push_back(in_matrix1_.rank1(row_ids[i]) - 1);
        else
            rows_matrix2.push_back(in_matrix1_.rank0(row_ids[i]) - 1);
    }

    auto res_matrix1 = matrix1_.get_rows(rows_matrix1);
    auto res_matrix2 = matrix2_.get_rows(rows_matrix2);

    std::vector<BinaryMatrix::SetBitPositions> res(row_ids.size());
    size_t pos_matrix1{};
    size_t pos_matrix2{};
    for (size_t i = 0 ; i < row_ids.size(); ++i) {
        if (in_matrix1_[row_ids[i]])
            res[i] = std::move(res_matrix1[pos_matrix1++]);
        else
            res[i] = std::move(res_matrix2[pos_matrix2++]);
    }

    return res;
}

template<class BaseMatrix1, class BaseMatrix2>
std::vector<BinaryMatrix::Row>
HybridMatrix<BaseMatrix1, BaseMatrix2>::
get_column(Column column) const {
    assert(column < num_columns());

    std::vector<Row> result;

    auto col1 = matrix1_.get_column(column);
    auto col2 = matrix2_.get_column(column);

    assert(std::is_sorted(col1.begin(), col1.end()));
    assert(std::is_sorted(col2.begin(), col2.end()));

    uint64_t pos_col1 = 0;
    uint64_t pos_col2 = 0;

    uint64_t val_col1 = 0;
    uint64_t val_col2 = 0;

    if (col1.size() && col2.size()) {
        val_col1 = in_matrix1_.select1(col1[pos_col1++] + 1);
        val_col2 = in_matrix1_.select0(col2[pos_col2++] + 1);

        while (true) {
            assert(val_col1 !=val_col2);
            if (val_col1 < val_col2) {
                result.push_back(val_col1);
                if (pos_col1 == col1.size()) {
                    result.push_back(val_col2);
                    break;
                }
                val_col1 = in_matrix1_.select1(col1[pos_col1++] + 1);
            } else {
                result.push_back(val_col2);
                if (pos_col2 == col2.size()) {
                    result.push_back(val_col1);
                    break;
                }
                val_col2 = in_matrix1_.select0(col2[pos_col2++] + 1);
            }
        }
    }

    while (pos_col1 < col1.size())
        result.push_back(in_matrix1_.select1(col1[pos_col1++] + 1));

    while (pos_col2 < col2.size())
        result.push_back(in_matrix1_.select0(col2[pos_col2++] + 1));

    return result;
}

template<class BaseMatrix1, class BaseMatrix2>
bool HybridMatrix<BaseMatrix1, BaseMatrix2>::load(std::istream &in) {
    return in_matrix1_.load(in) &&
           matrix1_.load(in) &&
           matrix2_.load(in);
}

template<class BaseMatrix1, class BaseMatrix2>
void HybridMatrix<BaseMatrix1, BaseMatrix2>::serialize(std::ostream &out) const {
    in_matrix1_.serialize(out);
    matrix1_.serialize(out);
    matrix2_.serialize(out);
}


} // namespace binmat
} // namespace annot
} // namespace mtg


#endif // __HYBRID_MATRIX__
