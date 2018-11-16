#ifndef __BIN_REL_WT_SDSL_HPP__
#define __BIN_REL_WT_SDSL_HPP__

#include <sdsl/wt_int.hpp>

#include "binary_matrix.hpp"


class BinRelWT_sdsl : public BinaryMatrix {
  public:
    BinRelWT_sdsl() : delimiters_(1, 1) {};

    typedef std::vector<Column> RowSetBits;
    typedef std::function<void(const RowSetBits &)> RowCallback;

    BinRelWT_sdsl(const std::function<void(const RowCallback &)> &generate_rows,
                  uint64_t num_set_bits,
                  uint64_t num_columns);

    explicit BinRelWT_sdsl(const BinRelWT_sdsl &other) = default;
    BinRelWT_sdsl& operator=(const BinRelWT_sdsl &other) = default;

    BinRelWT_sdsl(BinRelWT_sdsl&& other) = default;
    BinRelWT_sdsl& operator=(BinRelWT_sdsl&& other) = default;

    uint64_t num_columns() const;
    uint64_t num_rows() const;

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

  private:
    sdsl::wt_int<sdsl::rrr_vector<63>> wt_;
    bit_vector_rrr<> delimiters_;
    uint64_t num_columns_ = 0;
};

#endif // __BIN_REL_WT_SDSL_HPP__
