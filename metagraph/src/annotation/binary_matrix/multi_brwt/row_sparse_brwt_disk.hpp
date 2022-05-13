#ifndef __ROW_SPARSE_BRWT_DISK_HPP__
#define __ROW_SPARSE_BRWT_DISK_HPP__


#include <vector>
#include <unordered_map>
#include <memory>

#include "common/logger.hpp"

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/range_partition.hpp"
#include "common/ifstream_with_name_and_offset.hpp"

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_disk.hpp"
#include "annotation/binary_matrix/row_sparse_disk/row_sparse_disk.hpp"

using mtg::common::logger;

namespace mtg {
namespace annot {
namespace binmat {

using NodeDepth = size_t;

class RowSparseBRWT_Disk : public BinaryMatrix {
  BRWT_Disk brwt_disk;
  RowSparseDisk row_sparse_disk;
  mutable uint64_t n_row_sparse_queried_rows{}, n_brwt_queried_rows{}; // mkokot, TODO: debugowe usunac
  public:

    //typedef bit_vector_small stored_as_row_sparse_bv_type;   // mkokot, TODO: consider using bit_vector_smallrank
    typedef bit_vector_smallrank stored_as_row_sparse_bv_type; // mkokot, TODO: verify if it is not less performant because it seems that for 24750 BIGSI subset when queyring I got worse times:
    // slice_rows time: 4.733126812 <- time for bit_vector_smallrank
    // slice_rows time: 3.458182989 <- time for bit_vector_small
    
    RowSparseBRWT_Disk() = default;

    uint64_t num_columns() const override { return row_sparse_disk.num_columns(); }

    uint64_t num_rows() const override { return row_sparse_disk.num_rows() + brwt_disk.num_rows(); }        

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<Row> get_column(Column column) const override;
    
    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    void set_brwt_max_anno_mem(size_t _brwt_max_anno_mem) { brwt_disk.set_brwt_max_anno_mem(_brwt_max_anno_mem); }

    ~RowSparseBRWT_Disk()
    {
      std::cerr << "num row sparse queried rows: " <<  n_row_sparse_queried_rows << "\n"; // mkokot, TODO: remove
      std::cerr << "num brwt queried rows: " <<  n_brwt_queried_rows << "\n"; // mkokot, TODO: remove      
    }
  private:    
    stored_as_row_sparse_bv_type stored_as_row_sparse;
  
};

} // namespace binmat
} // namespace annot
} // namespace mtg


#endif //__ROW_SPARSE_BRWT_DISK_HPP__