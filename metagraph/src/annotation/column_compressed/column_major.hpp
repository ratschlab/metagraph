#ifndef __COLUMN_MAJOR_HPP__
#define __COLUMN_MAJOR_HPP__

#include <map>

#include "binary_matrix.hpp"
#include "bit_vector.hpp"


class ColMajorCompressed : public BinaryMatrix {
  public:
    ColMajorCompressed() {}
    ColMajorCompressed(const std::vector<std::unique_ptr<bit_vector_sd>> &columns);
    ColMajorCompressed(std::vector<std::unique_ptr<bit_vector_sd>>&& columns);

    ColMajorCompressed(const ColMajorCompressed &other) = default;
    ColMajorCompressed& operator=(const ColMajorCompressed &other) = default;

    ColMajorCompressed(ColMajorCompressed&& other) = default;
    ColMajorCompressed& operator=(ColMajorCompressed&& other) = default;

    uint64_t num_columns() const { return columns_.size(); }
    uint64_t num_rows() const;

    bool get(Row row, Column column) const;
    std::vector<Column> get_row(Row row) const;
    std::vector<Row> get_column(Column column) const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    // number of ones in the matrix
    uint64_t num_relations() const;

  private:
    std::vector<std::unique_ptr<bit_vector_sd>> columns_;
};

#endif // __COLUMN_MAJOR_HPP__
