#ifndef __VECTOR_ROW_BINMAT_HPP__
#define __VECTOR_ROW_BINMAT_HPP__

#include <vector>

#include "binary_matrix.hpp"
#include "utils.hpp"


class VectorRowBinMat : public BinaryMatrixRowDynamic {
  public:
    VectorRowBinMat(uint64_t num_rows = 0) : vector_(num_rows) {}

    VectorRowBinMat(const VectorRowBinMat &other) = default;
    VectorRowBinMat& operator=(const VectorRowBinMat &other) = default;

    VectorRowBinMat(VectorRowBinMat&& other) = default;
    VectorRowBinMat& operator=(VectorRowBinMat&& other) = default;

    uint64_t num_columns() const { return num_columns_; }
    uint64_t num_rows() const { return vector_.size(); }

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    void set(Row row, Column column);
    // do not check for existance before adding the label
    void force_set(Row row, Column column);
    void clear_row(Row row);

    // sort all indexes in rows and leave only distinct ones
    void standardize_rows();

    void insert_rows(const std::vector<Row> &rows);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    class StreamRows {
      public:
        StreamRows(std::ifstream &instream, const std::string filename);
        std::unique_ptr<std::vector<Row> > next_row();
      private:
        sdsl::int_vector_buffer<> *inbuf_;
        uint64_t i_ = 0;
    };
    static void write_rows(std::ofstream &outstream,
                    const std::string &filename,
                    const std::function<void (const std::function<void (const std::vector<uint64_t> &)>&)> &callback,
                    uint64_t num_rows,
                    uint64_t num_cols);

    // number of ones in the matrix
    uint64_t num_relations() const;
    // matrix density
    double density() const;

  private:
    uint64_t num_columns_ = 0;
    std::vector<SmallVector> vector_;
};

#endif // __VECTOR_ROW_BINMAT_HPP__
