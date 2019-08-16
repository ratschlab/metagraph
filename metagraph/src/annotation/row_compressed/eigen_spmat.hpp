#ifndef __EIGEN_SPMAT_HPP__
#define __EIGEN_SPMAT_HPP__

#include <vector>

#include <Eigen/Sparse>

#include "binary_matrix.hpp"


class EigenSpMat : public BinaryMatrixRowDynamic {
  public:
    EigenSpMat(uint64_t num_rows = 0, uint64_t max_num_cols = 10'000'000);

    EigenSpMat(const EigenSpMat &other) = default;
    EigenSpMat& operator=(const EigenSpMat &other) = default;

    EigenSpMat(EigenSpMat&& other) = default;
    EigenSpMat& operator=(EigenSpMat&& other) = default;

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return mat_.rows(); }

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    void set(Row row, Column column);
    void clear_row(Row row);

    void insert_rows(const std::vector<Row> &rows);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;
    // matrix density
    double density() const;

  private:
    void initialize(size_t num_rows, uint64_t max_num_cols);

    uint64_t num_columns_ = 0;
    Eigen::SparseMatrix<bool, Eigen::RowMajor> mat_;
};

#endif // __EIGEN_SPMAT_HPP__