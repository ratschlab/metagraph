#ifndef __UNIQUE_ROW_BINMAT_HPP__
#define __UNIQUE_ROW_BINMAT_HPP__

#include <vector>

#include "annotation/binary_matrix/base/binary_matrix.hpp"


class UniqueRowBinmat : public BinaryMatrix {
  public:
    typedef SetBitPositions row_type;

    explicit UniqueRowBinmat(uint64_t num_rows = 0);

    UniqueRowBinmat(std::vector<row_type>&& unique_rows,
                    std::vector<uint64_t>&& row_rank,
                    uint64_t num_columns);

    UniqueRowBinmat(const std::function<void(const RowCallback &)> &call_rows,
                    uint64_t num_columns);

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return row_rank_.size(); }

    bool get(Row row, Column column) const;
    SetBitPositions get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;
    // matrix density
    double density() const;

  private:
    uint64_t num_columns_ = 0;
    std::vector<row_type> unique_rows_;
    std::vector<uint64_t> row_rank_;
};

#endif // __UNIQUE_ROW_BINMAT_HPP__
